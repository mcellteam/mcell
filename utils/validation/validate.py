#!/usr/bin/env python3

"""
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to [http://unlicense.org] 
"""

import sys
import utils
import random
import multiprocessing
from collections import Counter

MAIN_MDL_FILE = 'Scene.main.mdl' 

# TODO: find automatically, allow to override
MCELL_EXECUTABLE = '/home/ahusar/src4_display/mcell/build/release/mcell'

g_extra_arg = ''

def run_mcell(seed):
    cmd = [MCELL_EXECUTABLE, '-seed', str(seed), MAIN_MDL_FILE] 
    if g_extra_arg:
        cmd.append(g_extra_arg)
    ec = utils.run(cmd, fout_name='mcell_' + str(seed) + '.out')
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
        curr_counts = get_molecule_counts_from_ascii_file('viz_data/seed_' + str(s).zfill(5) + '/Scene.ascii.1000.dat')
        
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
    if len(sys.argv) >= 2:
        nr_runs = int(sys.argv[1]) 

    if len(sys.argv) >= 3:
        global g_extra_arg 
        g_extra_arg = sys.argv[2] 
    
    seeds = generate_seeds(nr_runs)
    res = get_molecule_counts_for_multiple_runs(seeds)
    print(str(res))


if __name__ == "__main__":
    main()
    