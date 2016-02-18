#!/usr/bin/python
# DM      Sigma      Time (s)     Sample    Downfact


import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
#from matplotlib.widgets import Slider, Button, RadioButtons
from pulsar_tools import disp_delay
import math




class SinglePulse:
    """
    A class to contain all the relevant information for each single pulse detection 
    (i.e. S/N, box-car window size, DM, etc.). This is for ease of access during 
    plotting/other interactive stuff.
    """
    def __init__(self, DM, sig, t, samp, dfact):
        self.dm = DM          # float
        self.sigma = sig      # float
        self.time = t         # float
        self.sample = samp    # integer
        self.downfact = dfact # integer

    def print_params(self):
        print "DM:",self.dm
        print "Sigma:",self.sigma
        print "Time:",self.time
        print "Sample:",self.sample
        print "Downfactor:",self.downfact


class SPList:
    """
    A class to contain a number of SinglePulse objects in a numpy.array and grant easy acces to paramter 
    lists of those objects. Contains the original list of SinglePulse objects, and a list 
    of each object's: 
         DM, sigma, time, sample and downfactor.
    """
    def __init__(self, sp_list):
        self.list = np.array(sp_list)
        self.dm_list = np.array([sp.dm for sp in sp_list])
        self.sigma_list = np.array([sp.sigma for sp in sp_list])
        self.time_list = np.array([sp.time for sp in sp_list])
        self.sample_list = np.array([sp.sample for sp in sp_list])
        self.downfact_list = np.array([sp.downfact for sp in sp_list])
        



#def make_data_points(data):
#    """
#    Function to create all the SinglePulse objects for each single pulse in the data file.
#    Also creates an SPList object to containall the SinglePulse events.
#    Returns an SPList object.
#    """
#    # TODO: Need to re-work and see if there's a smarter way to populate all of the SinglePulse objects
#    DM = [float(row.split()[0]) for row in data]
#    Sigma = [float(row.split()[1]) for row in data]
#    Time = [float(row.split()[2]) for row in data]
#    Sample = [int(row.split()[3]) for row in data]
#    Downfact = [int(row.split()[4]) for row in data]
#
#    sp = [SinglePulse(dm, sig, time, samp, dfact) for dm, sig, time, samp, dfact \
#                       in zip(DM, Sigma, Time, Sample, Downfact)]
#    sp_list = SPList(sp)
#
#    return sp_list



def load_file(filename):
    if filename==None:
        print "No filename supplied to read..."

    elif filename.endswith('.singlepulse'):
        DM = np.genfromtxt(filename, comments="#", autostrip=True, usecols=0)
        Sigma = np.genfromtxt(filename, comments="#", autostrip=True, usecols=1)
        Time = np.genfromtxt(filename, comments="#", autostrip=True, usecols=2)
        Sample = np.genfromtxt(filename, comments="#", autostrip=True, usecols=3)
        Downfact = np.genfromtxt(filename, comments="#", autostrip=True, usecols=4)
    
        sp = [SinglePulse(dm, sig, time, samp, dfact) for dm, sig, time, samp, dfact \
                  in zip(DM, Sigma, Time, Sample, Downfact)]

        return SPList(sp)

    elif filename.endswith('.flag'):
        flags = np.genfromtxt(filename ,comments="#", autostrip=True)

        if len(flags)==0:
            print "No flags/bad times provided. Not times in final output will be masked."

        return flags

    else:
        print "File name suplied is not recognised. Must be either .singlepulse, .bad or .flag"


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
    print "%.2f seconds flagged from %.2f seconds of data (%.2f percent)" % ( flag_time, time[-1], flag_time/time[-1]*100)


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
    print bads

    if any(isinstance(b, np.ndarray) for b in bads):
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
    
    else:
        start = bads[0] - (padding + disp_delay(freq1=freq_l, freq2=freq_h, DM=max_DM)/1000)
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

    np.savetxt(basename+'.flag', flags, fmt='%d')
    check_call(['flag.sh', basename])
    #Popen(['flag.sh', basename]).communicate()[0]


def singlepulse_plot(basename=None, DMvTime=1, StatPlots=False, raw = False, threshold=5.0, movie=False):
    """
    Plots up the flagged data, should switch to using genfromtxt when I have the time.
       BWM: switched to using load_file to load singlepulse and flags. Uses genfromtext.
    """
    if raw:
        data = load_file(basename+'.singlepulse')
        flag_times = False
    else:
        #flag_times = load_file(basename+'.bad')
        flagfile(basename)
        data = load_file(basename+'_flagged.singlepulse')
        flags = load_file(basename+'.flag')


    data = SPList(data.list[np.where(data.sigma_list >= threshold)])

    #DM = [float(row.split()[0]) for row in data if float(row.split()[1]) >= threshold]
    #Sigma = [float(row.split()[1]) for row in data if float(row.split()[1]) >= threshold]
    #Time = [float(row.split()[2]) for row in data if float(row.split()[1]) >= threshold]
    #Sample = [int(row.split()[3]) for row in data if float(row.split()[1]) >= threshold]
    #Downfact = [int(row.split()[4]) for row in data if float(row.split()[1]) >= threshold]
    
    DM = data.dm_list
    Sigma = data.sigma_list
    Time = data.time_list
    Downfact = data.downfact_list

    Downfact_float = Downfact.astype(float)

    fig = plt.figure()
    cm = plt.cm.get_cmap('gist_rainbow')

    if StatPlots:
        ax0 = fig.add_subplot(231)
        plt.hist(Sigma, histtype='step', bins=60)
        ax0.set_xlabel('Signal-to-Noise',fontsize=18)
        ax0.set_ylabel('Number of Pulses',fontsize=18)

        ax1 = fig.add_subplot(232)
        plt.hist(DM, histtype='step', bins=int(0.5*len(set(DM))))
        ax1.set_xlabel('DM ($\mathrm{pc\, cm^{-3}}$)',fontsize=18)
        ax1.set_ylabel('Number of Pulses',fontsize=18)

        ax2 = fig.add_subplot(233, sharex=ax1) 
        # BWM: now shares x-axis with ax1, so changing DM on one will change range on the other
        plt.scatter(DM, Sigma, c=Downfact_float, cmap=cm, alpha=0.9)
        ax2.set_ylabel('Signal-to-Noise',fontsize=18)
        ax2.set_xlabel('DM ($\mathrm{p\, cm^{-3}}$)',fontsize=18)
        ax2.set_xlim([min(DM), max(DM)])
        ax2.set_ylim([min(Sigma), 1.1*max(Sigma)])

        ax3 = fig.add_subplot(212)

    else:
        ax3 = fig.add_subplot(111)

    # TODO: need to figure out how (if at all) we can make the axis sharing work 
    #       for x-axis to y-axis share

#	ax3.set_title("Single Pulse Sigma")    
    ax3.set_xlabel('Time (s)',fontsize=18)
    ax3.set_ylabel('DM ($\mathrm{pc\, cm^{-3}}$)',fontsize=18)
    ax3.set_ylim([min(DM), max(DM)])
    ax3.set_xlim([min(Time), max(Time)])
    #cm = plt.cm.get_cmap('gist_rainbow')

    # BWM: grab axis3 size to allocate marker sizes
    # TODO: is this really the best way to contrain the sizes?
    bbox_pix = ax3.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox_pix.width, bbox_pix.height
    area = width * height # axes area in inches^2 (apparently)

    Size = (Sigma/max(Sigma))**2
    Size = area*fig.dpi*(Size/max(Size))
    
    obs_stats(Time, flags)
#    sc=ax3.scatter(Time,DM, c=Sigma_float, vmin=min(Sigma_float), vmax=max(Sigma_float), cmap=cm, picker=1)
    sc = ax3.scatter(Time, DM, s=Size, c=Downfact_float, cmap=cm, vmin=min(Downfact_float),\
                       vmax=max(Downfact_float), picker=1, facecolor='none')
#	leg = ax1.legend()
         
    plt.colorbar(sc, label="Downfact", pad=0.01)
    # BWM: can't seem to get the bottom plot to extend the entire width when the color bar is active.
    fig.subplots_adjust(hspace=0.2, wspace=0.5)

    if not raw:
        if any(isinstance(l, np.ndarray) for l in flags):
            for flag in flags:
                flag_area = patches.Rectangle((float(flag[0]), min(DM)), \
                                                         float(flag[1])-float(flag[0]), \
                                                         max(DM)-min(DM), edgecolor='0', facecolor='0.66')
                ax3.add_patch(flag_area)

        else:
            flag_area = patches.Rectangle((float(flags[0]), min(DM)), \
                                                         float(flags[1])-float(flags[0]), \
                                                         max(DM)-min(DM), edgecolor='0', facecolor='0.66')
            ax3.add_patch(flag_area)




    def onpick(event):
        points = event.artist
        ind = event.ind
        mouseevent = event.mouseevent
        print '\n'
        print "Information for data points around click event %.4f, %.4f:" % (mouseevent.xdata,  mouseevent.ydata)
        for i in ind:		# These are fudge factors to turn samples into ms. 
            if ( DM[i] < 150):
                boxcar = Downfact[i]
            elif ( 150<= DM[i] < 823.2 ):
                boxcar = Downfact[i] * 2
            elif ( 823.2 <= DM[i] < 1486.2):
                boxcar = Downfact[i] * 2
            elif ( 1426.2 <= DM[i] < 2100):
                boxcar = Downfact[i] * 2
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
    #plt.tight_layout(w_pad=0.1, h_pad=0.1)
    plt.show()

    #obs_stats(Time, flags)




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
            

if __name__ == '__main__':

	modes=['interactive','movie']
    from optparse import OptionParser, OptionGroup
    parser=OptionParser(description="A python tool to plot, flag, and do otherwise with singlepulse search data from PRESTO")
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


