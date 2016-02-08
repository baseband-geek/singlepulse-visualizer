#!/usr/bin/python
# DM      Sigma      Time (s)     Sample    Downfact


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pulsar_tools import disp_delay
import math




def loadfile(filename):
    if filename==None:
        print "No filename supplied to read..."

    with open(filename) as f:
        data = f.read()

    data = data.split('\n')

    while not data[-1]:
        data = data[:-1]
    return data

def obs_stats(time, flags):
    # Not doing total time correctly, depends on last single pulse detection instead of observation time
    flag_time = 0

    for flag in flags:
        flag_time+=(float(flag.split()[1])-float(flag.split()[0]))
    print "%.2f seconds flagged from %.2f seconds of data (%.2f percent)" % ( flag_time, time[-1], flag_time/time[-1]*100)


def flagfile(basename, max_DM=2097.2, freq_l=0.169615, freq_h=0.200335, padding=3):
    """This function takes in a text file of bad 0 DM times and
    writes out one flagged over the correct de-dispersive smearing
    times, looking for overlaps along the way. There must be a text file named
    basename.bad with rows indicating bad times for this to work. """
    from subprocess import Popen


    bads = np.genfromtxt(basename+'.bad')
    i = 0 # initialize counter for new list
    flags = []
    for bad in bads:
        start = bad[0] - (padding + disp_delay(freq1=freq_l, freq2=freq_h, DM=max_DM)/1000)
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
    np.savetxt(basename+'.flag', flags, fmt='%d')
    Popen(['flag.sh', basename]).communicate()[0]



def singlepulse_plot(basename=None, DMvTime=1, StatPlots=False, raw = False, threshold=5.0	):
    """Plots up the flagged data, should switch to using genfromtxt when I have the time"""
    if raw:
        data=loadfile(basename+'.singlepulse')[1:]
        flags = False
    else:
        data = loadfile(basename+'_flagged.singlepulse')[1:]
        flags = loadfile(basename+'.flag')


    DM = [float(row.split()[0]) for row in data if float(row.split()[1]) >= threshold]
    Sigma = [float(row.split()[1]) for row in data if float(row.split()[1]) >= threshold]
    Time = [float(row.split()[2]) for row in data if float(row.split()[1]) >= threshold]
    Sample = [int(row.split()[3]) for row in data if float(row.split()[1]) >= threshold]
    Downfact = [int(row.split()[4]) for row in data if float(row.split()[1]) >= threshold]


    Sigma_float = [float(value) for value in Sigma]
    Size = [value**1.7 for value in Sigma_float]
    Downfact_float = [float(value) for value in Downfact]

    fig = plt.figure()

    if StatPlots:
        ax0 = fig.add_subplot(231)
        plt.hist(Sigma, histtype='step', bins=60)
        ax0.set_xlabel('Signal-to-Noise')
        ax0.set_ylabel('Number of Pulses')

        ax1 = fig.add_subplot(232)
        plt.hist(DM, histtype='step', bins=int(0.5*len(set(DM))))
        ax1.set_xlabel('DM (pc cm^-3)')
        ax1.set_ylabel('Number of Pulses')

        ax2 = fig.add_subplot(233)
        plt.scatter(DM, Sigma, c=Downfact_float, alpha=0.9)
        ax2.set_ylabel('Signal-to-Noise')
        ax2.set_xlabel('DM (pc cm^-3)')
        ax2.set_xlim([min(DM), max(DM)])
        ax2.set_ylim([min(Sigma), max(Sigma)])

        ax3 = fig.add_subplot(212)
    else:
        ax3 = fig.add_subplot(111)

#	ax3.set_title("Single Pulse Sigma")    
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('DM (pc cm^-3)')
    ax3.set_ylim([min(DM), max(DM)])
    ax3.set_xlim([min(Time), max(Time)])
    cm = plt.cm.get_cmap('gist_rainbow')
#    sc=ax3.scatter(Time,DM, c=Sigma_float, vmin=min(Sigma_float), vmax=max(Sigma_float), cmap=cm, picker=1)
    sc=ax3.scatter(Time,DM, s=Size, alpha=0.33, c=Downfact_float, cmap=cm, vmin=min(Size), vmax=max(Size), picker=1)
#	leg = ax1.legend()
    plt.colorbar(sc, label="Downfact")
    if not raw:
        for flag in flags:
            flag_area = matplotlib.patches.Rectangle((float(flag.split()[0]), min(DM)), float(flag.split()[1])-float(flag.split()[0]), max(DM)-min(DM),  edgecolor='0', facecolor='0.66')
            ax3.add_patch(flag_area)

    def onpick(event):
        points = event.artist
        ind = event.ind
        mouseevent = event.mouseevent
        print '\n'
        print "Information for data points around click event %.4f, %.4f:" % (mouseevent.xdata,  mouseevent.ydata)
        for i in ind:
            if ( DM[i] < 150):
                boxcar = Downfact[i]
            elif ( 150<= DM[i] < 823.2 ):
                boxcar = Downfact[i] * 2
            elif ( 823.2 <= DM[i] < 1486.2):
                boxcar = Downfact[i] * 5
            elif ( 1426.2 <= DM[i] < 2100):
                boxcar = Downfact[i] * 10
            print "%.2f seconds, %.2f Sigma event detected at a DM of %.2f with a boxcar of: %d ms" % (Time[i], Sigma[i], DM[i], boxcar)

    fig.canvas.mpl_connect('pick_event', onpick)

    '''
    ax2 = fig.add_subplot(122)
    ax2.set_title("Single Pulse Boxcar")
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('DM (pc cm^-3)')
    cm = plt.cm.get_cmap('RdYlBu')
    sc2=ax2.scatter(Time,DM, c=Downfact_float, vmin=min(Downfact_float), vmax=max(Downfact_float), cmap=cm)
#	leg = ax1.legend()
    plt.colorbar(sc2)
    if not raw:
        for flag in flags:
            flag_area = matplotlib.patches.Rectangle((float(flag.split()[0]), min(DM)), float(flag.split()[1])-float(flag.split()[0]), max(DM)-min(DM),  edgecolor='0', facecolor='0.66')
            ax2.add_patch(flag_area)
    '''
    fig.suptitle('Single Pulse Search results for '+basename)
    plt.show()

    obs_stats(Time, flags)




def slice(infile, dm=None, timerange=None, sigma=None, downfact=None):
    # Not properly implemented yet

    data = read_singlepulse(infile)

    slices = [None]*5

    slice_map = {'dm':0, 'sigma':1, 'timerange':2, 'sample':3, 'downfact':4}




    DM = [row.split()[0] for row in data]
    Sigma = [row.split()[1] for row in data]
    Time = [row.split()[2] for row in data]
    Sample = [row.split()[3] for row in data]
    Downfact = [row.split()[4] for row in data]

    if dm:
        if type(dm) == type(0) or type(0.0):
            data = [row   for row in data if dm <= row.split()[0]]
        elif type(dm) == type([]):
            data = [row  for row in data if dm[0] <= row.split()[0] <= dm[1]]
    if sigma:
        if type(sigma) == type(0) or type(0.0):
            data = [row for row in data if sigma <= row.split()[1]  ]
        elif type(sigma) == type([]):
            data = [row for row in data if sigma[0] <= row.split()[1] <= sigma[1]]




