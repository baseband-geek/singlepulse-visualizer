#!/usr/bin/env python
# DM      Sigma      Time (s)     Sample    Downfact


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as patches
import matplotlib.pyplot as plt

#from matplotlib.widgets import Slider, Button, RadioButtons
from pulsar_tools import disp_delay
import math
from subprocess import call, check_call, Popen, PIPE
import time
import multiprocessing as mp

class AsyncPlotter():
    # Written by Thomas Robitaille and pulled from https://gist.github.com/astrofrog/1453933
    
    def __init__(self, processes=mp.cpu_count()):

        self.manager = mp.Manager()
        self.nc = self.manager.Value('i', 0)
        self.pids = []
        self.processes = processes

    def async_plotter(self, nc, fig, filename, processes):
        while nc.value >= processes:
            time.sleep(0.1)
        nc.value += 1
        print "Plotting " + filename
        fig.savefig(filename, format='png')
        plt.close(fig)
        nc.value -= 1

    def save(self, fig, filename):
        p = mp.Process(target=self.async_plotter,
                       args=(self.nc, fig, filename, self.processes))
        p.start()
        self.pids.append(p)

    def join(self):
        for p in self.pids:
            p.join()


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


def singlepulse_plot(basename=None, DMvTime=1, StatPlots=True, raw = False, threshold=5.0, obsid=None, dm_pairs=[[0,2000]], directory="./"):
    """
    Plots up the flagged data, should switch to using genfromtxt when I have the time.
       BWM: switched to using load_file to load singlepulse and flags. Uses genfromtext.
    """
    if raw:
        data = load_file(directory + basename + '.singlepulse')
        flag_times = False
    else:
        #flag_times = load_file(basename+'.bad')
        flagfile(directory + basename)
        data = load_file(directory + basename + '_flagged.singlepulse')
        flags = load_file(directory + basename + '.flag')


    data = SPList(data.list[data.sigma_list >= threshold])

    Downfact_float = data.downfact_list.astype(float)

    
    #cm = plt.cm.get_cmap('gist_rainbow')
    
    t_0 = 0 
    time_step = 2
    num_steps = 4
    time_span =  200
    last_time = data.time_list.max()
    
    asplot = AsyncPlotter()
    
    def make_frame(t_frame):
        print "Making fromes for timestep {0} out of {1}".format((t_frame+1), int(last_time/time_step))
        time= time_step*t_frame
        timerange_data = SPList(data.list[ (time <= data.time_list) & ( data.time_list <= (time + time_span))])
        for dm_pair in dm_pairs:
            fig = plt.figure(figsize=(11.69,8.27))
            temp_data = SPList(timerange_data.list[(dm_pair[0] <= timerange_data.dm_list) & ( timerange_data.dm_list <= dm_pair[1])])
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
                ax1.set_xlim(dm_pair[0], dm_pair[1])
                ax2 = fig.add_subplot(233, sharex=ax1) # BWM: now shares x-axis with ax1, so changing DM on one will change range on the other
                try:
                    plt.scatter(temp_data.dm_list, temp_data.sigma_list, alpha=0.9)
                except ValueError:
                    print "Plot contains no points"
                ax2.set_ylabel('Signal-to-Noise')
                ax2.set_xlabel('DM ($\mathrm{p\, cm^{-3}}$)')
                ax2.set_xlim([dm_pair[0], dm_pair[1]])
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
                        flag_area = patches.Rectangle((int(dm_pair[0]), float(flag[0]) ), \
														  (int(dm_pair[1])),\
															  (float(flag[1]) - float(flag[0])), \
															  edgecolor='0', facecolor='0.66')
                        ax3.add_patch(flag_area)

                else:
                    flag_area = patches.Rectangle((float(flags[0]), dm_pair[0]), \
													  (float(flags[1]) - float(flags[0])), \
													  (dm_pair[1]), \
													  edgecolor='0', facecolor='0.66')
                    ax3.add_patch(flag_area)

			# TODO: need to figure out how (if at all) we can make the axis sharing work 
			#       for x-axis to y-axis share	

			#	ax3.set_title("Single Pulse Sigma")    
            ax3.set_ylabel('Time (s)')
            ax3.set_xlabel('DM ($\mathrm{pc\, cm^{-3}}$)')
            ax3.set_xlim([dm_pair[0], dm_pair[1]])
            ax3.set_ylim([time, time + time_span])
            cm = plt.cm.get_cmap('gist_rainbow')
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
            #fig.subplots_adjust(hspace=0.2, wspace=0.5)

            fig.suptitle('Single Pulse Search results for ' + basename) 
            #		plt.tight_layout(w_pad=0.1, h_pad=0.1)
            #		writer.grab_frame()
            f_name="_dm{0:0>4}_{1:0>4}_tmp{2:0>10}.png".format(dm_pair[0],dm_pair[1],t_frame)
            #		plt.savefig(f_name)
            #plt.savefig(f_name, bbox_inches='tight', dpi=400)
            asplot.save(fig, f_name)
            #plt.savefig(f_name, format='png')
            plt.close(fig)
    
    [make_frame(t) for t in range(t_0, int((last_time-time_span)/time_step))]
    for dm_pair in dm_pairs:
        make_movie="avconv -f image2 -i _dm{0:0>4}_{1:0>4}_tmp%10d.png {2}_dm_{0}_{1}_sp.mp4".format(dm_pair[0], dm_pair[1], obsid)
        call(make_movie,shell=True)
    obs_stats(data.time_list, flags)


def galaxy_batchjobs(obs_id=None, DMvTime=1, StatPlots=True, raw = False, threshold=5.0, obsid=None, dm_pairs=[[0,2000]]):
    workdir = "/scratch2/mwaops/stremblay/data/observations/{0}/".format(obs_id)
    for dm_pair in dm_pairs:
        batch_name = "{0}{1}_sp_movie_{2}_{3}.batch".format(workdir+"batch_files/",obs_id, dm_pair[0], dm_pair[1])
        batch_out = batch_name[:-5]+"out"
        with open(batch_name,'w') as batch_file:
            batch_line = "#!/bin/bash -l\n#SBATCH --time=12:00:00\n#SBATCH \n#SBATCH --output={0}\n#SBATCH --export=NONE\n#SBATCH -p workq\n".format(batch_out)
            batch_file.write(batch_line)
            batch_line = "aprun make_sp_movie.py --obsid {0} --dm_ranges {1},{2}\n".format(obs_id, dm_pair[0], dm_pair[1])
            batch_file.write(batch_line)
            
        submit_line = "sbatch --partition=workq --workdir={0} {1}\n".format(workdir+"singlepulse/",batch_name)
        submit_cmd = Popen(submit_line,shell=True,stdout=PIPE)
        jobid=""
        for line in submit_cmd.stdout:
            if "Submitted" in line:
                (word1,word2,word3,jobid) = line.split()

if __name__ == '__main__':

    from optparse import OptionParser, OptionGroup
    parser = OptionParser(description="A python tool to plot, flag, and do otherwise with singlepulse search data from PRESTO")
    parser.add_option("--dm_ranges", action="store", type="string", default="0,500,1000,1500,2000", help="A list of the DM ranges you want plotted. Write as a comma separated integers. [default=%default]")
    parser.add_option("--obsid", action="store", type="string", help="Observation ID or other basename for files. [No default]")
    parser.add_option("--threshold", action="store", type="float", default=5.0, help="S/N threshold. [default=%default]")
    parser.add_option("--raw", action="store_true", default=False, help="Plots the data without any flagging [default=%default]")
    parser.add_option("--galaxy", action="store_true", default=False, help="Boolean trigger controlling splitting each DM pair into separate batch jobs [default=%defalut]")
    parser.add_option("--local", action="store_true", default=False, help="Use this flag to bypass 'standard' directory structure [default=%default]")
    
    (opts, args) = parser.parse_args()	# Parse string into a list, then build a list of pairs for dm ranges to be plotted
    opts.dm_ranges = [dm for dm in opts.dm_ranges.split(",")]
    dm_pairs = [[int(i),int(j)] for i,j in zip(opts.dm_ranges[:-1], opts.dm_ranges[1:])]
    if len(dm_pairs) > 2:
        dm_pairs.insert(0,[int(opts.dm_ranges[0]),int(opts.dm_ranges[-1])])
    if opts.galaxy:
        galaxy_batchjobs(obs_id=opts.obsid, DMvTime=1, StatPlots=True, raw = opts.raw, threshold=opts.threshold, obsid=opts.obsid, dm_pairs=dm_pairs)
        quit()
        
    if opts.local:
        wdir = "./"
    else:
        wdir = "/scratch2/mwaops/stremblay/data/observations/{0}/singlepulse/".format(opts.obsid)
    start=time.time()
    singlepulse_plot(basename=opts.obsid, DMvTime=1, StatPlots=True, raw = opts.raw, threshold=opts.threshold, obsid=opts.obsid, dm_pairs=dm_pairs, directory=wdir)
    print "Took {0} seconds".format(time.time()-start)


