#!/usr/bin/env python
# DM      Sigma      Time (s)     Sample    Downfact


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from joblib import Parallel, delayed

#from matplotlib.widgets import Slider, Button, RadioButtons
from pulsar_tools import disp_delay
import math
from subprocess import call, check_call, Popen, PIPE




class SinglePulse:
    """
    A class to contain all the relevant information for each single pulse detection 
    (i.e. S/N, box-car window size, DM, etc.). This is for ease of access during 
    plotting/other interactive stuff.
    """
    def __init__(self, DM, sig, t, samp, dfact):
        self.dm       = DM          # float
        self.sigma    = sig      # float
        self.time     = t         # float
        self.sample   = samp    # integer
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
        self.list          = np.array(sp_list)
        self.dm_list       = np.array([sp.dm for sp in sp_list])
        self.sigma_list    = np.array([sp.sigma for sp in sp_list])
        self.time_list     = np.array([sp.time for sp in sp_list])
        self.sample_list   = np.array([sp.sample for sp in sp_list])
        self.downfact_list = np.array([sp.downfact for sp in sp_list])
    def __iter__(self):
    	return self
    def __next__(self):
    	try:
    		result = self.text[self.index].upper()
    	except IndexError:
    		raise StopIteration
    	self.index += 1
    	return result





def load_file(filename):
    if filename==None:
        print "No filename supplied to read..."

    elif filename.endswith('.singlepulse'):
        DM       = np.genfromtxt(filename, comments="#", autostrip=True, usecols=0)
        Sigma    = np.genfromtxt(filename, comments="#", autostrip=True, usecols=1)
        Time     = np.genfromtxt(filename, comments="#", autostrip=True, usecols=2)
        Sample   = np.genfromtxt(filename, comments="#", autostrip=True, usecols=3)
        Downfact = np.genfromtxt(filename, comments="#", autostrip=True, usecols=4)
    
        sp = [SinglePulse(dm, sig, time, samp, dfact) for dm, sig, time, samp, dfact \
                  in zip(DM, Sigma, Time, Sample, Downfact)]

        return SPList(sp)

    elif filename.endswith('.flag'):
        flags = np.genfromtxt(filename ,comments="#", autostrip=True)

        if len(flags) == 0:
            print "No flags/bad times provided. Not times in final output will be masked."

        return flags

    else:
        print "File name suplied is not recognised. Must be either .singlepulse, .bad or .flag"



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
    # call flag.sh script to create masked singlepulse file
    check_call(['flag.sh', basename])
    #Popen(['flag.sh', basename]).communicate()[0]


def singlepulse_plot(basename=None, DMvTime=1, StatPlots=True, raw = False, threshold=5.0, obsid=None):
    """
    Plots up the flagged data, should switch to using genfromtxt when I have the time.
       BWM: switched to using load_file to load singlepulse and flags. Uses genfromtext.
    """
    if raw:
        data = load_file(basename + '.singlepulse')
        flag_times = False
    else:
        #flag_times = load_file(basename+'.bad')
        flagfile(basename)
        data = load_file(basename + '_flagged.singlepulse')
        flags = load_file(basename + '.flag')


    data = SPList(data.list[np.where(data.sigma_list >= threshold)])

    Downfact_float = data.downfact_list.astype(float)

    
    cm = plt.cm.get_cmap('gist_rainbow')
    
    t_0 = 0 
    time_step = 2
    num_steps = 4
    time_span =  200
    last_time = data.time_list.max()
    
    def make_frame(t_frame):
#    for t in range(t_0,t_0+num_steps):
#    for t in range(t_0, int(last_time/time_step)):
		print "Making frome {0} out of {1}".format((t_frame+1), int(last_time/time_step))
		time= time_step*t_frame
		fig = plt.figure()
#		temp_data = SPList(data.list[np.where( data.time_list <= (time + time_span))])
		temp_data = SPList(data.list[ (time <= data.time_list) & ( data.time_list <= (time + time_span))])
		#    temp_data=[pulse  for pulse in data if (timerange[0] <= pulse.list.time <= timerange[1])] 
		if StatPlots:
			ax0 = fig.add_subplot(231)
			try:
				plt.hist(temp_data.sigma_list, histtype='step', bins=60)
			except ValueError:
				print "Plot contains no points"
			ax0.set_xlabel('Signal-to-Noise')
			ax0.set_xlim(data.sigma_list.min()-0.5)
			ax0.set_ylabel('Number of Pulses')	

			ax1 = fig.add_subplot(232)
			try:
				plt.hist(temp_data.dm_list, histtype='step', bins=int(0.5 * len(set(temp_data.dm_list))))
			except ValueError:
				print "Plot contains no points"
			ax1.set_xlabel('DM ($\mathrm{pc\, cm^{-3}}$)')
			ax1.set_ylabel('Number of Pulses')
			ax2 = fig.add_subplot(233, sharex=ax1)	 
			# BWM: now shares x-axis with ax1, so changing DM on one will change range on the other
			try:
				plt.scatter(temp_data.dm_list, temp_data.sigma_list, alpha=0.9)
			except ValueError:
				print "Plot contains no points"
			#plt.scatter(temp_data.dm_list, temp_data.sigma_list, c=Downfact_float, cmap=cm, alpha=0.9)
			ax2.set_ylabel('Signal-to-Noise')
			ax2.set_xlabel('DM ($\mathrm{p\, cm^{-3}}$)')
			try:
				ax2.set_xlim([0, temp_data.dm_list.max()])
			except ValueError:
				ax2.set_xlim([0, data.dm_list.max()])
			try:
				ax2.set_ylim([data.sigma_list.min(), 1.1 * temp_data.sigma_list.max()])	
			except ValueError:
				ax2.set_ylim([data.sigma_list.min(), 1.1 * data.sigma_list.max()])

			ax3 = fig.add_subplot(212)	

		else:
			ax3 = fig.add_subplot(111)
		if not raw:
			if any(isinstance(l, np.ndarray) for l in flags):
				for flag in flags:
					flag_area = patches.Rectangle((0, float(flag[0]) ), \
													  (2000),\
														  (float(flag[1]) - float(flag[0])), \
														  edgecolor='0', facecolor='0.66')
					ax3.add_patch(flag_area)

			else:
				flag_area = patches.Rectangle((float(flags[0]), 0), \
												  (float(flags[1]) - float(flags[0])), \
												  (2000), \
												  edgecolor='0', facecolor='0.66')
				ax3.add_patch(flag_area)

		# TODO: need to figure out how (if at all) we can make the axis sharing work 
		#       for x-axis to y-axis share	

		#	ax3.set_title("Single Pulse Sigma")    
		ax3.set_ylabel('Time (s)')
		ax3.set_xlabel('DM ($\mathrm{pc\, cm^{-3}}$)')
		ax3.set_xlim([0, 2000])
		ax3.set_ylim([time, time + time_span])
		#cm = plt.cm.get_cmap('gist_rainbow')

		# grab axis3 size to allocate marker sizes
		bbox_pix = ax3.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		width, height = bbox_pix.width, bbox_pix.height
		area = width * height # axes area in inches^2 (apparently)

		#TODO: need to try and use something like percentiles to make sure that just one
		# big pulse doesn't swamp the sizes or colorbars.
		
		Size = (3. * area / 2.) * (data.sigma_list**2 / np.percentile(data.sigma_list, 99.5))
		
		
		#    sc=ax3.scatter(Time,DM, s=Size, c=Sigma, vmin=min(Sigma), vmax=max(Sigma),\
		#                       cmap=cm, picker=1)
		sc = ax3.scatter(data.dm_list, data.time_list, s=Size, c=Downfact_float, cmap=cm, \
							 vmin=Downfact_float.min(), vmax=Downfact_float.max(), facecolor='none')
		#	leg = ax1.legend()
	 
		#    plt.colorbar(sc, label="Sigma", pad=0.01) 
		#    plt.colorbar(sc, label="Downfact", pad=0.01)
		# BWM: can't seem to get the bottom plot to extend the entire width when the color bar is active. SET: removing bar for now for just this reason
		fig.subplots_adjust(hspace=0.2, wspace=0.5)

		fig.suptitle('Single Pulse Search results for ' + basename) 
	#		plt.tight_layout(w_pad=0.1, h_pad=0.1)
	#		writer.grab_frame()
		f_name="_tmp{0:0>10}.png".format(t_frame)
		plt.savefig(f_name)
	#		plt.savefig(f_name, bbox_inches='tight')
		plt.close(fig)
    
    [make_frame(t) for t in range(t_0, int(last_time/time_step))]
#    Parallel(n_jobs=4)(delayed(make_frame)(t) for t in range(t_0,t_0+num_steps))
    make_movie="avconv -f image2 -i _tmp%10d.png {0}_sp.mp4".format(obsid)
    call(make_movie,shell=True)
    obs_stats(data.time_list, flags)




if __name__ == '__main__':

    from optparse import OptionParser, OptionGroup
    parser = OptionParser(description="A python tool to plot, flag, and do otherwise with singlepulse search data from PRESTO")
    parser.add_option("--dm_range", action="store", type="string", nargs=2, default=(0,2000), help="(Not yet implemented) The lowest and highest DM to plot. [default=%default]")
    parser.add_option("--obsid", action="store", type="string", help="Observation ID or other basename for files. [No default]")
    parser.add_option("--threshold", action="store", type="float", default=5.0, help="S/N threshold. [default=%default]")
    parser.add_option("--raw", action="store_true", default=False, help="Plots the data without any flagging [default=%default]")
    
    (opts, args) = parser.parse_args()
	
    singlepulse_plot(basename=opts.obsid, DMvTime=1, StatPlots=True, raw = opts.raw, threshold=opts.threshold, obsid=opts.obsid)

