#!/usr/bin/env python3

"""
Copyright (C) 2019 by
The Salk Institute for Biological Studies and
Pittsburgh Supercomputing Center, Carnegie Mellon University

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

For the complete terms of the GNU General Public License, please see this URL:
http://www.gnu.org/licenses/gpl-2.0.html
"""


import sys
import utils
import random
import multiprocessing
from collections import Counter

MAIN_MDL_FILE = 'Scene.main.mdl' 

# TODO: find automatically, allow to override
MCELL_EXECUTABLE = '/nadata/cnl/home/ahusar/src4/mcell_tools/work/build_mcell/mcell'


def run_mcell(seed):
    cmd = [MCELL_EXECUTABLE, '-seed', str(seed), MAIN_MDL_FILE] 
    ec = utils.run(cmd, fout_name='mcell.out')
    if ec != 0:
        print("Fatal error, run with seed " + str(seed) + "failed.")
        sys.exit(1)
        

def get_molecule_counts_from_ascii_file(filename):
    counts = {}
    with open(filename, 'r') as fin:
        for line in fin:
            # the first item is the ID
            id = line.split(' ')[0]
            if id in counts:
                counts[id] = counts[id] + 1
            else:
                counts[id] = 1
                
    return counts


# single threaded execution for now
def get_molecule_counts_for_multiple_runs(seeds):
    counts = {}
    current_run = 1
    #for s in seeds:
    #    print("RUN " + str(current_run) + ", SEED " + str(s))
    #    current_run += 1
        
    #    run_mcell(s)
    
    # Set up the parallel task pool to use all available processors
    count = 12 # multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    
    # Run the jobs
    pool.map(run_mcell, seeds)

    for s in seeds:
        # TODO: find the last .ascii file automatically
        curr_counts = get_molecule_counts_from_ascii_file('viz_data/seed_' + str(s).zfill(5) + '/Scene.ascii.50.dat')
        
        # add values with common key
        counts = Counter(counts) + Counter(curr_counts) 

    return counts


def generate_seeds(count):
    res  = []
    
    random.seed(a=200)
    for i in range(0, count):
        res.append(random.randint(1, 65535))
        
    return res

# just one variant for now, but this could grow int oa more versatile tool
def main():
    
    nr_runs = 1
    if len(sys.argv) == 2:
        nr_runs = int(sys.argv[1]) 
    
    seeds = generate_seeds(nr_runs)
    res = get_molecule_counts_for_multiple_runs(seeds)
    print(str(res))


if __name__ == "__main__":
    main()
    