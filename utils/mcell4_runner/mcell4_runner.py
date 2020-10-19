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

MCELL_DIR = ''
RUNNING_MARKER = 'running.marker'
FINISHED_MARKER = 'finished.marker'

class Options:
    def __init__(self):
        self.seeds_str = None
        self.args_file = None
        self.max_cores = None
        self.main_model_file = None
    
    
def create_argparse():
    parser = argparse.ArgumentParser(description='MCell4 Runner')
    parser.add_argument(
        '-s', '--seeds', type=str, 
        help='seeds in the form first:last:step, e.g. 1:100:2 will use seeds 1 through 100 in steps of 2, '
             'model must accept "-seed N" argument, '
             'model must make sure that the output directories are different for different seeds')
    parser.add_argument(
        '-a', '--args-file', type=str, 
        help='arguments file where each line contains arguments passed to the main model file, '
             'model must make sure that the output directories are different for different arguments')
    parser.add_argument('-j', '--max-cores', type=int, 
        help='sets maximum number of cores for running, default is all if -s is not used')
    parser.add_argument('main_model_file',  
        help='sets path to the MCell4 model')
    return parser


def process_opts():
    parser = create_argparse()
    args = parser.parse_args()
    opts = Options()
    
    if args.seeds and args.args_file:
        print("Error: only one argument -s/--seeds or -a/--args-file must be specified, not both.")
        sys.exit(1)

    if not args.seeds and not args.args_file:
        print("Error: one of arguments -s/--seeds or -a/--args-file must be specified.")
        sys.exit(1)
         
    if args.seeds:
        opts.seeds_str = args.seeds
    
    if args.args_file:
        if os.path.exists(args.args_file): 
            opts.args_file = args.args_file
        else:
            print("Error: file " + args.args_file + " does not exist.")
            sys.exit(1)
         
    if args.max_cores:
        seed.max_cores = args.max_cores
        
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
        

def check_prerequisites():
    global MCELL_DIR
    MCELL_DIR = os.environ.get('MCELL_DIR', '')
    if MCELL_DIR:
        sys.path.append(os.path.join(MCELL_DIR, 'lib'))
    else:
        print("Error: system variable MCELL_DIR that is used to find the mcell library was not set.")
        sys.exit(1)
    
    mcell_so_path = os.path.join(MCELL_DIR, 'lib', 'mcell.so')
    if not os.path.exists(mcell_so_path):
        print("Could not find library '" + mcell_so_path + ".")
        sys.exit(1)
    
            
def generate_seeds(seeds_str):
    seeds_info = seeds_str.split(':')
    if len(seeds_info) != 3 or \
        not seeds_info[0].isdigit() or \
        not seeds_info[1].isdigit() or \
        not seeds_info[2].isdigit():
        print("Error: invalid seed string, must be in for form min:max:step, given '" + seeds_str + "'.")
        sys.exit(1)
    
    res  = []
    for i in range(int(seeds_info[0]), int(seeds_info[1]) + 1, int(seeds_info[2])):
        res.append(i)
    
    return res


def prepare_args(opts):
    if opts.seeds_str:
        seeds = generate_seeds(opts.seeds_str)
        return ['-seed ' + str(i) for i in seeds ]
    elif opts.args_file:
        res = []
        with open(opts.args_file, 'r') as f:
            for line in f:
                if not line.isspace():
                    if line.endswith('\n'):
                        res.append(line[:-1])
                    else:
                        res.append(line[:-1])
        return res
        
    
def run_mcell4(args_str_w_opts):
    args_str = args_str_w_opts[0]
    opts = args_str_w_opts[1]
    
    cmd_str = sys.executable + ' ' + opts.main_model_file + ' '
    cmd_str += args_str
    
    print("Running " + cmd_str)
    
    args_log = args_str.replace(' ', '_').replace('-', '_')
    log_name = os.path.splitext(opts.main_model_file)[0] + '_' + args_log + '.mcell4.log'
    
    exit_code = 1
    with open(log_name, "w") as f:
        proc = subprocess.Popen(cmd_str, shell=True, cwd=os.getcwd(), stdout=f, stderr=subprocess.STDOUT)
        proc.communicate()
        exit_code = proc.returncode

    with open(log_name, "a") as f:
        f.write("DIR:" + os.getcwd() + "\n")
        f.write("CMD:" + cmd_str + "\n")
    
    if exit_code != 0:
        print("MCell4 failed, see '" + os.path.join(os.getcwd(), log_name) + "'.")
        return exit_code
    else:
        return 0

    
def run_mcell4_parallel(opts, args):
    
    # set up the parallel task pool to use all available processors ot the 
    # maximum specified
    if opts.max_cores:
        cpu_count = int(opts.max_cores)
    else:
        cpu_count = multiprocessing.cpu_count()
     
    pool = multiprocessing.Pool(processes=cpu_count)
    
    # run the jobs
    args_str_w_opts = zip(args, itertools.repeat(opts))
    res_codes = pool.map(run_mcell4, args_str_w_opts)

    num_total = 0
    num_failed = 0
    for i in range(len(res_codes)):
        c = res_codes[i]
        if c != 0:
            print("MCell run with args '" + args[i] + "' failed with exit code " + str(c) + ".")
            num_failed += 1
        
        num_total += 1

    if num_failed == 0:
        print("Finished, all runs passed.")
        return 0
    else:
        print("Finished with errors, " + str(num_failed) + "/" + str(num_total) + " runs failed.")
        return 1    
    
            
def file_markers_start():    
    if (os.path.exists(FINISHED_MARKER)):
        os.remove(FINISHED_MARKER)
    pathlib.Path(RUNNING_MARKER).touch()


def file_markers_finish():    
    if (os.path.exists(RUNNING_MARKER)):
        os.remove(RUNNING_MARKER)
    pathlib.Path(FINISHED_MARKER).touch()


if __name__ == '__main__':
    file_markers_start()
    
    check_prerequisites()
    
    opts = process_opts()
        
    args = prepare_args(opts)
    
    exit_code = run_mcell4_parallel(opts, args)
    
    file_markers_finish()
    
    sys.exit(exit_code)