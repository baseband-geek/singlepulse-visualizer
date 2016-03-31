#!/usr/bin/env python

import os
import sys

def sort_singlepulse(basename, directory='.',verbose=False):
    """
    Accepts the base name (usually Observation ID) and the directory where the relevant files are located. 
    If no directory argument is given, assumes all files are in current working directory (.)

    Creates a total singlepulse file from all singlpulse files with the given base name. 
    Ensures unique entries only, and the file output is sorted in time (and therefore sample).
    """
    from operator import itemgetter
    from os import listdir
    import numpy as np

    if verbose not in ['True','False']:
        print "Running sort_singlepulse with verbose mode off..."
        verbose = False
    # grab all files with relevant basename in current directory
    base_files = sorted([f for f in listdir(directory) if f.startswith(basename)])
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

    #sys.exit(0)
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
        f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format('DM','Sigma','Time(s)','Sample',\
                                                                'Downfact','inf_file'))
        for event in ordered:
            f.write(''.join('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(*event)))


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print "If using as stand alone script, supply file basename to sort as first argument, the directory in which they are located as the second and whether you want verbose output as the third."
        print "basename here is defined as: basename_DM*.fits"
        print "Verbose mode is set to False by default. Send third argument as True if output is desired."
        print "e.g.               python sort_singlepulse.py filebasename /path/to/singlepulse/files True"

        sys.exit(0)
    elif len(sys.argv) == 2:
        print "Assuming singlepulse files with basename {0} are in the current working directory: {1}".format(sys.argv[1],os.getcwd())
        sort_singlepulse(sys.argv[1])
    elif len(sys.argv) == 3:
        print "Sorting singlepulse events from basename {0} in {1}".format(sys.argv[1],sys.argv[2])
        sort_singlepulse(sys.argv[1],sys.argv[2])
    elif len(sys.argv) == 4:
        print "Sorting singlepulse events from basename {0} in {1}, with verbose output set to {2}".format(sys.argv[1],sys.argv[2],sys.argv[3])
        sort_singlepulse(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print "Too many argument supplied. Just need: basename and (optionally) the directory where the singlepulse files are located."
        sys.exit(0)
 





    
