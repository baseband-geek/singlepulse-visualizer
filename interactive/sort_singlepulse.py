#!/usr/bin/env python

import os
import sys
import argparse
from operator import itemgetter
import numpy as np

# The following class is taken directly from PRESTO, written by Scott Ransom. 
# It allows easy access to the information in the .inf files produced by PRESTO's subroutines. 
class infodata:
    def __init__(self, filenm):
        self.breaks = 0
        for line in open(filenm):
            if line.startswith(" Data file name"):
                self.basenm = line.split("=")[-1].strip()
                continue
            if line.startswith(" Telescope"):
                self.telescope = line.split("=")[-1].strip()
                continue
            if line.startswith(" Instrument"):
                self.instrument = line.split("=")[-1].strip()
                continue
            if line.startswith(" Object being observed"):
                self.object = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Right Ascension"):
                self.RA = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Declination"):
                self.DEC = line.split("=")[-1].strip()
                continue
            if line.startswith(" Data observed by"):
                self.observer = line.split("=")[-1].strip()
                continue
            if line.startswith(" Epoch"):
                self.epoch = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Barycentered?"):
                self.bary = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of bins"):
                self.N = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Width of each time series bin"):
                self.dt = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Any breaks in the data?"):
                self.breaks = int(line.split("=")[-1].strip())
                if self.breaks:
                    self.onoff = []
                continue
            if line.startswith(" On/Off bin pair"):
                vals = line.split("=")[-1].strip().split(",")
                self.onoff.append((int(vals[0]), int(vals[1])))
                continue
            if line.startswith(" Type of observation"):
                self.waveband = line.split("=")[-1].strip()
                continue
            if line.startswith(" Beam diameter"):
                self.beam_diam = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Dispersion measure"):
                self.DM = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Central freq of low channel"):
                self.lofreq = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Total bandwidth"):
                self.BW = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of channels"):
                self.numchan = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Channel bandwidth"):
                self.chan_width = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Data analyzed by"):
                self.analyzer = line.split("=")[-1].strip()
                continue



def sort_singlepulse(basename, directory, verbose):
    """
    Accepts the base name (usually Observation ID) and the directory where the relevant files are located. 
    If no directory argument is given, assumes all files are in current working directory (via os.getcwd())

    Creates a total singlepulse file from all singlpulse files with the given base name. 
    Ensures unique entries only, and the file output is sorted in time (and therefore sample).
    """

    # grab all files with relevant basename in current directory
    base_files = sorted([f for f in os.listdir(directory) if f.startswith(basename)])
    sp_files = [s for s in base_files if s.endswith('.singlepulse') and '_DM' in s]
    
    # create a list of single pulse events from the .singlepulse file
    sp_events = []
    for sp in sp_files:
        empty = False
        if verbose: print "loading data from: {0}".format(sp)
        
        # load data from files, unles file is empty in which case do nothing and move on
        try:
            data = np.genfromtxt(sp, comments='#', skip_header=1)
        except:
            if verbose: print "empty file. not appending to events list."
            empty = True
        
        if empty is False:
            inf_file = sp.replace('.singlepulse', '.inf')

            # create info object with each paramter in the .inf file as an attribute
            info = infodata(inf_file)
            
            # create list of info to append to singlepulse data
            inf_list = [inf_file, info.telescope, info.RA, info.DEC, info.epoch,\
                        info.N, info.dt, info.lofreq, info.BW,\
                        info.numchan, info.chan_width]

            if any(isinstance(d, np.ndarray) for d in data):

                # contains 2 or more lines
                for d in data:
                    sp_events.append(np.append(d, inf_list).tolist())
    
            else:
                # is only a single line file
                sp_events.append(np.append(d, inf_list).tolist())



    # annoying but seemingly necessary type conversions for sorting
    for s in sp_events:
        s[0]=float(s[0])
        s[1]=float(s[1])
        s[2]=float(s[2])
        s[3]=int(float(s[3]))
        s[4]=int(float(s[4]))

    # create a list of tuples, keeping only unique pulse events. Output is a list of tuples
    ordered = sorted(set(map(tuple, sp_events)), key=itemgetter(2))

    if verbose: print "writing {0}.singlepulse".format(basename)
    
    with open(basename+'.singlepulse','wb') as f:
        f.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15}\n'.format('DM','Sigma','Time','Sample',\
                                                                                                 'Downfact','inf_filename','telescope',\
                                                                                                 'point_RA','point_DEC','epoch',\
                                                                                                 'total_samples','dt','lochan_centre',\
                                                                                                 'BW','num_chan','chan_width'))
        for event in ordered:
            f.write(''.join('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15}\n'.format(*event)))




parser = argparse.ArgumentParser(description=\
"Sort PRESTO singlepulse files into one comma-separated-values (csv) file.",\
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("basename", type=str, help="basename of singlepulse files, e.g. basename_DM100.00.singlepulse")

parser.add_argument('-wdir', type=str, action='store', help="Working directory where singlepulse files are located.",\
 default=os.getcwd())

parser.add_argument('-v', action='store_true', help="Use verbose mode.", default=False)



args = parser.parse_args()

sort_singlepulse(args.basename, args.wdir, args.v)
