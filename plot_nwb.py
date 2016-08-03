#!/usr/local/python-2.7.8/bin/python
#
# plot_nwb.py
#
# Developed by Gennady Denisov on 2016-08-03.

import sys
import os

import h5py
import matplotlib.pyplot as plt
from pylab import *
import re

import numpy as np
import optparse

# ------------------------------------------------------------------------------
success = False
data_path = ""
time_path = ""

def plot_nwb_command_line_parser(parser):
    parser.add_option("-D","--debug",   action="store_true",  dest="debug", help="output debugging info", default=False)
    parser.add_option("-e","--no_error_handling", action="store_false", dest="handle_errors", help="handle_errors", default=True)
    parser.add_option("-n","--ts_name", dest="ts_name",help="plot the data for the specified timeseries", default="")
    parser.add_option("-t", "--tint",   dest="time_interval", help="time intervals to be used", metavar="time_interval", default="")
    parser.add_option("-T", "--trial",  dest="trial_id", help="trial for which data will be visualized", metavar="trial", default="")
    parser.add_option("-o","--outfolder",dest="output_folder", help="output folder (default=same as input folder)",metavar="output_folder",default="")
    parser.add_option("-r","--replace",  action="store_true", dest="replace", help="if the output file already exists, replace it", default=False)
    parser.add_option("-s","--summary",action="store_true",dest="summary",help="output summary table", default=False)
    parser.add_option("-v","--verbose",  action="store_true", dest="verbose", help="increase the verbosity level of output", default=False)

    return parser

# ------------------------------------------------------------------------------
# Retuen a list of all keyNames and corresponding timeseries_names

def update_summary(summary_table, trials, nwb_object, item_path, options):
    item_pointer = nwb_object[item_path]
    item_name   = item_path.split("/")[len(item_path.split("/"))-1]
    if options.debug:        
        print "item_path=", item_path
    try:
        keys = item_pointer.keys()
        # item is a group
        if options.verbose:
            print "item is a group"
        if not re.search("Trial", item_name):
            for k in keys:       
                if item_path == "/":
                   item_path1 = item_path + k
                else: 
                   item_path1 = item_path + '/' + k
                summary_table, trials = update_summary(summary_table, trials, nwb_object, item_path1, options)
        else:
            if not item_name in trials:
                trials.append(str(item_name))
    except: 
        # item is dataset
        if options.verbose:
            print "item is a dataset"
        series_name = item_path.split("/")[len(item_path.split("/"))-2]

        if item_name == "data" and not re.search("unit", series_name):
            if "keyName" in nwb_object[item_path].attrs.keys():
                keyName = nwb_object[item_path].attrs["keyName"]
            else:
                keyName = "          "

            data, timestamps = get_data_by_ts_name(series_name, nwb_object, options)
            time_range = str(np.min(timestamps)) + "-" + str(np.max(timestamps))

            if len(series_name) < 8:
                series_name = series_name + "    "
            if len(keyName) < 8:
                keyName = keyName + "    "
            summary_table.append([series_name, keyName, time_range, item_path])
    return (summary_table, trials)

# ------------------------------------------------------------------------------

def parse_item(target_name, nwb_object, item_path, options):
    global success
    global data_path
    global time_path 
    item_pointer = nwb_object[item_path]
    try:
        # item is group
        keys = item_pointer.keys()
        if options.debug:        
            print "item_path=", item_path, " keys=", keys
        if len(keys) > 0:
            for k in keys:
                if success:
                    break
                if item_path == "/":
                    item_path1 = "/" + k
                else:
                    item_path1 = item_path + '/' + k
                parse_item(target_name, nwb_object, item_path1, options)
    except:
        # item is dataset
        item_name = item_path.split("/")[len(item_path.split("/"))-1]
        if item_name == "data":
            series_name = item_path.split("/")[len(item_path.split("/"))-2]
            if options.debug:       
                print "item_path=", item_path, " item_name=", item_name, " series_name=", series_name, " target_name=", target_name
            if series_name == target_name:
                data_path = item_path[0:(len(item_path)-len(item_name))] + "data"
                time_path = item_path[0:(len(item_path)-len(item_name))] + "timestamps"
                success = True

# ------------------------------------------------------------------------------

def get_time_range(ts_name, nwb_object, options):
    global time_path

    timestamps = []

    parse_item(ts_name, nwb_object, "/", options)

    if len(data_path) > 0 and len(time_path) > 0:
        timestamps = np.array(nwb_object[time_path])

    time_min   = min(timestamps[~np.isnan(timestamps)])
    time_max   = max(timestamps[~np.isnan(timestamps)])
    time_range = str(time_min) + "-" + str(time_max)

    return time_range           

# ------------------------------------------------------------------------------

def get_data_by_ts_name(ts_name, nwb_object, options):
    global data_path
    global time_path

    data = []
    timestamps = []

    parse_item(ts_name, nwb_object, "/", options)

    if len(data_path) > 0 and len(time_path) > 0:
        data       = np.array(nwb_object[data_path])
        timestamps = np.array(nwb_object[time_path])
        
        # Eliminate nan values
        good_inds  = ~np.isnan(data      ) & ~np.isnan(timestamps)
        data       = data[      good_inds]
        timestamps = timestamps[good_inds]
    return (data, timestamps)

# ------------------------------------------------------------------------------

def output_summary(nwb_object, options):
    summary_table = []
    trials = []
    summary_table, trials = update_summary(summary_table, trials, nwb_object, "/", options)
    print "\nSummary table:\n"
    sys.stdout.write("%s\t%s\t\t%s\t\t%s\n" % ("time_series", "keyName", "time_range", "path"))
    sys.stdout.write("%s" % ("-----------------------------------------------------------------------\n"))
    for i in range(len(summary_table)):
        sys.stdout.write("%s\t%s\t%s\t%s\n" % (summary_table[i][0], summary_table[i][1], summary_table[i][2], summary_table[i][3]))
    sys.stdout.write("\n" % ())
    print "Trial ids:\n\n", trials, "\n"

# ------------------------------------------------------------------------------

def plot_nwb(nwb_object, options):

    data, timestamps = get_data_by_ts_name(options.ts_name, nwb_object, options)

    if len(data) > 0 and len(data) == len(timestamps):       
        ibeg = 0
        iend = len(timestamps)-1
        if len(options.time_interval) > 0:
            tmin = np.float32(options.time_interval.split("-")[0])
            tmax = np.float32(options.time_interval.split("-")[1])
            ibeg = np.min(np.where(timestamps >= tmin))
            iend = np.max(np.where(timestamps <= tmax))
        elif len(options.trial_id) > 0:
            # Extract time interval for a given trial
            try:
                start_path = "/epochs/" + str(options.trial_id) + "/start_time"
                stop_path  = "/epochs/" + str(options.trial_id) + "/stop_time"
                tmin = np.array(nwb_object[start_path])*1000
                tmax = np.array(nwb_object[stop_path])*1000
                if options.verbose:
                    print "tmin=", tmin, " tmax=", tmax
                ibeg = np.min(np.where(timestamps >= tmin))
                iend = np.max(np.where(timestamps <= tmax))
            except:
                try:
                    time_path = "/epochs/" + str(options.trial_id) + "/" + options.ts_name + "/timeseries/timestamps"
                    print "time_path=", time_path
                    timestamps1 = np.array(nwb_object[time_path])
                    tmin = min(timestamps1[~isnan(timestamps1)])
                    tmax = max(timestamps1[~isnan(timestamps1)])
                    ibeg = np.min(np.where(timestamps >= tmin))
                    iend = np.max(np.where(timestamps <= tmax))
                except:
                    print "Data for Trial_"  + str(options.trial_id) + " not found"
        if options.verbose:
            print "tmin=", tmin, " tmax=", tmax, " min(timestamps)=", min(timestamps), " max(timestamps)=", max(timestamps), \
                  " ibeg=", ibeg, " iend=", iend, " num_points=", iend-ibeg+1
        if min(data[~np.isnan(data)]) < max(data[~np.isnan(data)]):
            if options.verbose:
                print "Plotting data vs timestamps (len=", iend-ibeg+1, "):\n"
            plt.plot(timestamps[ibeg:iend], data[ibeg:iend], 'b')
            plt.xlabel('time ')
            plt.ylabel(options.ts_name)
            plt.title('Plot of data vs timestamps')
            plt.grid(True)
            plt.savefig("data.svg")
            plt.show()
        else:
            if options.verbose:
                print "Plotting timestamps (len=", iend-ibeg+1, "):\n"
            figure()
            for i in range(ibeg, iend+1):
                x = [timestamps[i],timestamps[i]]
                y = [0, 1]
                plot(x, y, 'b')
            axis([timestamps[ibeg], timestamps[iend], 0, 2])
            show()
    else:
        "No data to plot\n"

# ------------------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog <data_file>.nwb -t ts_name [other options (-h to list)]"

    parser = optparse.OptionParser(usage=usage)
    parser = plot_nwb_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.verbose:
        print "len(args)=", len(args)

    if len(args) == 1:
        data_path     = args[0]
        if not re.search(".nwb", data_path):
            sys.exit("\nScript plot_mwb.py accepts as input a NWB data file only")
        
        data_basename = os.path.basename(data_path)
        if len(options.output_folder) == 0:
            options.output_folder = os.path.dirname(data_path)

        nwb_object = h5py.File(data_path, "r")

        if options.summary:            
            output_summary(nwb_object, options)
        elif len(options.ts_name) > 0:
            plot_nwb(nwb_object, options)

        nwb_object.close()
    else:
        parser.print_usage()
        sys.exit(2)

