#!/usr/bin/python
# DM      Sigma      Time (s)     Sample    Downfact


import numpy as np
import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from pulsar_tools import disp_delay
import math
import sys
import pandas as pd

from bokeh.io import output_file, show, vform, push_notebook
from bokeh.models import CustomJS, ColumnDataSource, Slider, HoverTool
from bokeh.plotting import *
from ipywidgets import interact

output_notebook()



def load_file(filename):
    if filename == None:
        print "No filename supplied to read..."

    elif filename.endswith('.singlepulse'):
        print "Generating singlepulse file Data Frame..."
        data = pd.read_csv(filename)
        return data

    elif filename.endswith('.flag'):
        print "Generating flags list..."

        flags = np.genfromtxt(filename ,comments="#", autostrip=True)

        if len(flags) == 0:
            print "No flags/bad times provided. No times in final output will be masked."

        return flags

    else:
        print "File name suplied is not recognised. Must be either .singlepulse, .bad or .flag"
        




"""
def load_file(filename):
    if filename==None:
        print "No filename supplied to read..."

    elif filename.endswith('.singlepulse'):
        DM       = np.genfromtxt(filename, comments="#", autostrip=True, usecols=0, skip_header=1)
        Sigma    = np.genfromtxt(filename, comments="#", autostrip=True, usecols=1, skip_header=1)
        Time     = np.genfromtxt(filename, comments="#", autostrip=True, usecols=2, skip_header=1)
        Sample   = np.genfromtxt(filename, comments="#", autostrip=True, usecols=3, skip_header=1)
        Downfact = np.genfromtxt(filename, comments="#", autostrip=True, usecols=4, skip_header=1)
        inf_file = np.genfromtxt(filename, comments="#", autostrip=True, usecols=5, dtype=str, skip_header=1)

        sp = [SinglePulse(dm, sig, time, samp, dfact, inf) for dm, sig, time, samp, dfact, inf \
                  in zip(DM, Sigma, Time, Sample, Downfact, inf_file)]

        return SPList(sp)
        #return sp

    elif filename.endswith('.flag'):
        flags = np.genfromtxt(filename ,comments="#", autostrip=True)

        if len(flags) == 0:
            print "No flags/bad times provided. Not times in final output will be masked."

        return flags

    else:
        print "File name suplied is not recognised. Must be either .singlepulse, .bad or .flag"
"""

#def load_flags(filename):
#    if filename==None:
#        print "No filename supplied to read into flags..."
#
#    flags = np.genfromtxt(filename ,comments="#", autostrip=True)
#    if len(flags)==0:
#        print "No flags provided. Not times in final output will be hidden."
#
#    return flags


def obs_stats(time, flags):
    # Not doing total time correctly, depends on last single pulse detection instead of observation time
    flag_time = 0

    # BWM: if there is only 1 masked region, flags is a list, 
    # if there are 2+ masked regions, flags is a list of lists.
    if any(isinstance(l, np.ndarray) for l in flags):
        for flag in flags:
            flag_time += (float(flag[1]) - float(flag[0]))
    else:
        flag_time = float(flags[1]) - float(flags[0])
    print "{0:.2f} seconds flagged from {1:.2f} seconds of data ({2:.2f} percent)".format(flag_time, time, flag_time * 100 / time)




def flagfile(basename, max_DM=2097.2, freq_l=0.169615, freq_h=0.200335, padding=3):
    """
    This function takes in a text file of bad 0 DM times and
    writes out one flagged over the correct de-dispersive smearing
    times, looking for overlaps along the way. There must be a text file named
    basename.bad with rows indicating bad times for this to work. 
    """
    from subprocess import check_call

    # BWM: originally planned to move this to the load_file function, 
    #  but left it incase we JUST want to call flagfile
    bads = np.genfromtxt(basename+'.bad', comments='#', autostrip=True)

    # BWM: again because how np.genfromtxt works, if there is only 1 bad line, we get a list, 
    # if there are 2+ bad lines we get a list of lists. So have to check for np.ndarray 
    # instances and change method accordingly.
    

    i = 0 # initialize counter for new list
    flags = []

    if any(isinstance(b, np.ndarray) for b in bads):
        for bad in bads:
            start = bad[0] - (padding + disp_delay(freq1=freq_l, freq2=freq_h, DM=max_DM)/1000.0)
            if start < 0:
                start = 0
            stop = bad[1] + padding
            if len(flags) > 0:
                if start <= flags[-1][1]:
                    flags[-1][1] = stop
                else:
                    flags.append([start, stop])
            else:
                flags.append([start, stop])
    
    else:
        # if there is a no bad regions (defaulted) then don't put any padding in
        if bads[0] == bads[1]:
            padding = 0
            max_DM = 0

        start = bads[0] - (padding + disp_delay(freq1=freq_l, freq2=freq_h, DM=max_DM)/1000.0)
        if start < 0:
            start = 0
        stop = bads[1] + padding
        if len(flags) > 0:
            if start <= flags[-1][1]:
                flags[-1][1] = stop
            else:
                flags.append([start, stop])
        else:
            flags.append([start, stop])
    # save new file  as basename.flag
    np.savetxt(basename+'.flag', flags, fmt='%d')
    # call flag.sh script to creat masked singlepulse file
    check_call(['flag.sh', basename])
    #Popen(['flag.sh', basename]).communicate()[0]


def singlepulse_plot(basename=None, DMvTime=1, StatPlots=False, raw=False, threshold=5.0, movie=False):
    """
    Plots up the flagged data, should switch to using genfromtxt when I have the time.
       BWM: switched to using load_file to load singlepulse and flags. Uses pandas data frame/genfromtxt.
    """

    print "Make sure you have run sort_singlepulse.py to gather the single pulse events into the one file {0}.singelpulse".format(basename)

    if raw:
        data = load_file(basename + '.singlepulse')
        #flag_times = False
    else:
        #flag_times = load_file(basename+'.bad')
        try:
            flagfile(basename) # BWM: should we be providing appropriate freqs and DM for this?
        except:
            print "No {}.bad file given. Creating one with entry [0 0]".format(basename)
            f=open('{}.bad'.format(basename),'w')
            f.write('0 0')
            f.close()
            print "Saved {}.bad".format(basename)
            print "Retrying..."
            flagfile(basename)

        data = load_file(basename + '_flagged.singlepulse')
        flags = load_file(basename + '.flag')

    # ensure only using data points with sigma > threshold
    data = data[data['Sigma']>threshold]

    cm = plt.cm.get_cmap('gist_rainbow')

    Size = 5 * data['Sigma'].values/4.0
    percentile = np.percentile(Size, 99.5)
    Size[np.where(Size > 5 * percentile)] = 5 * percentile/4.0

    TOOLS = "resize,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select"
    color = 'navy'
    color_array= [color]*len(Size)
    source = ColumnDataSource(data=dict(time=data['Time'].values,\
                                        dm=data['DM'].values,\
                                        sigma=data['Sigma'].values,\
                                        size=Size,color=color_array,\
                                        downfact=data['Downfact'].values.astype(float)))
    

    source2 = ColumnDataSource(data=dict(time=[],dm=[],sigma=[]))

    timeseries = figure(plot_width=900, plot_height=400, webgl=True,tools=TOOLS)
    timeseries.scatter('time','dm',source=source, marker='o',size='size',color='color')
    timeseries.xaxis.axis_label = 'Time'
    timeseries.yaxis.axis_label = 'DM'

    if StatPlots:
        top_left = figure(plot_width=300, plot_height=300 ,webgl=True, tools=TOOLS)
        hist,edges = np.histogram(data['Sigma'].values, bins=50)
        #source.add(hist,name='sigma_hist')
        top_left.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],line_color=color,fill_color=color)
        top_left.xaxis.axis_label = 'S/N'
        top_left.yaxis.axis_label = 'Number of pulses'

        top_mid = figure(plot_width=300, plot_height=300, webgl=True,tools=TOOLS)
        hist,edges = np.histogram(data['DM'].values, bins=50)
        #source.add(hist,name='dm_hist')
        top_mid.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],line_color=color,fill_color=color)
        top_mid.xaxis.axis_label = 'DM'
        top_mid.yaxis.axis_label = 'Number of pulses'

        top_right = figure(plot_width=300, plot_height=300, webgl=True,tools=TOOLS,\
                           x_range=top_mid.x_range)
        top_right.scatter('dm','sigma',source=source, marker='o',color='color')
        top_right.xaxis.axis_label = 'DM'
        top_right.yaxis.axis_label = 'S/N'

        

    if StatPlots:
        plots = [gridplot([[top_left,top_mid,top_right]]),timeseries]
        show(vplot(*plots))
    else:
        show(timeseries)



            

if __name__ == '__main__':

    modes = ['interactive','movie']
    from optparse import OptionParser, OptionGroup
    parser = OptionParser(description="A python tool to plot, flag, and do otherwise with singlepulse search data from PRESTO")
    parser.add_option("-m", "--mode", type="choice", choices=['interactive','movie'], help="Mode you want to run. {0}".format(modes))
    parser.add_option("--dm_range", action="store", type="string", nargs=2, default=(0,2000), help="(Not yet implemented) The lowest and highest DM to plot. [default=%default]")
    parser.add_option("--obsid", action="store", type="string", help="Observation ID or other basename for files. [No default]")
    parser.add_option("--threshold", action="store", type="float", default=5.0, help="S/N threshold. [default=%default]")
	
    (opts, args) = parser.parse_args()
    if opts.mode == 'movie':
    	singlepulse_plot(basename=opts.obsid, DMvTime=1, StatPlots=True, raw = False, threshold=opts.threshold, movie=True)
    elif opts.mode == 'interactive':
        singlepulse_plot(basename=opts.obsid, DMvTime=1, StatPlots=True, raw=False, threshold=opts.threshold, movie=False)
    else:
        print "Somehow your non-standard mode snuck through. Try again with one of {0}".format(modes)
        quit()


