#!/usr/bin/env python

"""
Copyright (C) 2020 by
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

import os
import sys
import shutil
import glob
import pathlib
import argparse
import itertools
import multiprocessing
import re
import subprocess

from mcell4_runner import file_markers_start, file_markers_finish, generate_seeds

TEST_BNGL = 'test.bngl'

class Options:
    def __init__(self):
        self.seeds_str = None
        self.bng2pl_path = None
        self.max_cores = None
        self.main_model_file = None
    
    
def create_argparse():
    parser = argparse.ArgumentParser(description='MCell4 Runner')
    parser.add_argument(
        '-s', '--seeds', type=str, 
        help='seeds in the form first:last:step, e.g. 1:100:2 will use seeds 1 through 100 in steps of 2, '
             'model must accept "-seed N" argument, '
             'the output directories are different for different seeds, ODE model ignores this argument')
    parser.add_argument('-j', '--max-cores', type=int, 
        help='sets maximum number of cores for running, default is all if -j is not used')
    parser.add_argument('-b', '--bng2pl', type=str, 
        help='sets path to BNG2.pl')
    parser.add_argument('main_model_file',  
        help='sets path to the BNGL model')
    return parser


def process_opts():
    parser = create_argparse()
    args = parser.parse_args()
    opts = Options()
    
    if args.seeds:
        opts.seeds_str = args.seeds
    else:
        print("Error: argument -s/--seeds must be specified.")
        sys.exit(1)

    if args.bng2pl:
        if os.path.exists(args.bng2pl): 
            opts.bng2pl_path = args.bng2pl
        else:
            print("Error: file " + args.bng2pl + " does not exist.")
            sys.exit(1)
    else:
        print("Error: argument -b/--bng2pl must be specified.")
        sys.exit(1)

    if args.max_cores:
        opts.max_cores = args.max_cores
        
    if args.main_model_file:
        if os.path.exists(args.main_model_file): 
            opts.main_model_file = args.main_model_file
        else:
            print("Error: file " + args.main_model_file + " does not exist.")
            sys.exit(1)
    else:
        print("Error: main model file must be specified as a positional argument.")
        sys.exit(1) 
        
    return opts
        

def run_bng(abs_dir, opts):
    
    os.chdir(abs_dir)
    
    cmd_str = 'perl ' + opts.bng2pl_path + ' ' + TEST_BNGL
    print("Running " + cmd_str + " in " + abs_dir)
        
    log_name = opts.main_model_file + '_' + os.path.basename(abs_dir) + '.bng2pl.log'
    
    exit_code = 1
    with open(log_name, "w") as f:
        proc = subprocess.Popen(cmd_str, shell=True, cwd=os.getcwd(), stdout=f, stderr=subprocess.STDOUT)
        proc.communicate()
        exit_code = proc.returncode

    with open(log_name, "a") as f:
        f.write("DIR:" + os.getcwd() + "\n")
        f.write("CMD:" + cmd_str + "\n")

    if exit_code != 0:
        print("BNG2.pl failed, see '" + os.path.join(os.getcwd(), log_name) + "'.")
        return exit_code
    else:
        return 0
    
    
def find_in_file(fname, search_for):
    lines = []
    with open(fname, "r") as infile:
        for line in infile:
            if search_for in line:
                return line
    return ''   
 
 
def replace_in_file(fname, search_for, replace_with):
    lines = []
    with open(fname, "r") as infile:
        for line in infile:
            line = line.replace(search_for, replace_with)
            lines.append(line)
    with open(fname, "w") as outfile:
        for line in lines:
            outfile.write(line) 
 
        
def run_bng_parallel(opts, seeds):
    
    # nfsim or ode?
    # does not handle comments
    line = find_in_file(opts.main_model_file, 'method=>"nf"')
    if not line:
        # ODE
        dir = 'ode'
        os.mkdir(dir)
        shutil.copy(TEST_BNGL, dir)
                 
        self.run_bng(dir, opts)
        
    else:
        # NFSim - multiple runs are needed
        cwd = os.getcwd()
        dirs = []
        for s in seeds:
            dir = 'nf_' + str(s).zfill(5)
            if not os.path.exists(dir):
                os.mkdir(dir)
            shutil.copy(opts.main_model_file, os.path.join(dir, TEST_BNGL))
            
            # copy also all other .bngl files from the main file's directory
            files = glob.iglob(os.path.join(os.path.dirname(opts.main_model_file), "*.bngl"))
            for file in files:
                if file != opts.main_model_file and os.path.isfile(file):
                    shutil.copy2(file, dir)
        
            # update seed value
            replace_in_file(os.path.join(dir, TEST_BNGL), 'seed=>1', 'seed=>' + str(s))
            dirs.append(os.path.join(cwd, dir))
        
        if opts.max_cores:
            # maximum number of processes specified
            cpu_count = int(opts.max_cores)
        else:
            cpu_count = multiprocessing.cpu_count()
        
        pool = multiprocessing.Pool(processes=cpu_count)

        # run in parallel        
        res_codes = pool.starmap(run_bng, zip(dirs, itertools.repeat(opts)), 1)
        
        num_total = 0
        num_failed = 0
        for i in range(len(res_codes)):
            c = res_codes[i]
            if c != 0:
                print("BNG run with seed '" + str(seeds[i]) + "' failed with exit code " + str(c) + ".")
                num_failed += 1
    
        if num_failed == 0:
            print("Finished, all runs passed.")
            return 0
        else:
            print("Finished with errors, " + str(num_failed) + "/" + str(num_total) + " runs failed.")
            return 1   

                
if __name__ == '__main__':
    file_markers_start()
    
    opts = process_opts()
    
    seeds = generate_seeds(opts.seeds_str)

    exit_code = run_bng_parallel(opts, seeds)
    
    file_markers_finish()
    
    sys.exit(exit_code)
    