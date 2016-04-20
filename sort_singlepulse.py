#!/usr/bin/env python

import os
import sys
import argparse
from operator import itemgetter
import numpy as np

def sort_singlepulse(basename, directory=os.getcwd(), verbose=False):
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
            if any(isinstance(d, np.ndarray) for d in data):
                # contains 2 or more lines
                for d in data:
                    sp_events.append(np.append(d, sp.replace('.singlepulse', '.inf')).tolist())
    
            else:
                # is only a single line file
                sp_events.append(np.append(d, sp.replace('.singlepulse', '.inf')).tolist())

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
        f.write('{0},{1},{2},{3},{4},{5}\n'.format('DM','Sigma','Time(s)','Sample',\
                                                                'Downfact','inf_file'))
        for event in ordered:
            f.write(''.join('{0},{1},{2},{3},{4},{5}\n'.format(*event)))




parser = argparse.ArgumentParser(description=\
"Sort PRESTO singlepulse files into one comma-separated-values (csv) file.",\
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("basename", type=str, help="basename of singlepulse files, e.g. basename_DM100.00.singlepulse")

parser.add_argument('-wdir', type=str, action='store', help="Working directory where singlepulse files are located.",\
 default="current working directory, {0}".format(os.getcwd()))

parser.add_argument('-v', action='store_true', help="Use verbose mode.", default=False)



args = parser.parse_args()

sort_singlepulse(args.basename, args.work_dir, args.verbose)
