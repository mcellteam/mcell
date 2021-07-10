#!/usr/bin/env python

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

import os
import sys
import shutil
import glob
import argparse
import itertools
import multiprocessing
import re
import subprocess

MCELL_PATH = ''
RUNNING_MARKER = 'running.marker'
FINISHED_MARKER = 'finished.marker'
LOGS_DIR = 'logs'

class Options:
    def __init__(self):
        self.seeds_str = None
        self.extra_arg = ''
        self.args_file = None
        self.max_cores = None
        self.main_model_file = None
    
    
def create_argparse(name):
    parser = argparse.ArgumentParser(description=name + ' Runner')
    parser.add_argument(
        '-s', '--seeds', type=str, 
        help='seeds in the form first:last:step, e.g. 1:100:2 will use seeds 1 through 100 in steps of 2, '
             'model must accept "-seed N" argument, '
             'model must make sure that the output directories are different for different seeds, '
             'can be used with argument -x')
    parser.add_argument(
        '-x', '--extra-arg', type=str, 
        help='optional extra argument to be passed when -s is used')
    parser.add_argument(
        '-a', '--args-file', type=str, 
        help='arguments file where each line contains arguments passed to the main model file, '
             'model must make sure that the output directories are different for different arguments')
    parser.add_argument('-j', '--max-cores', type=int, 
        help='sets maximum number of cores for running, default is all if -j is not used')
    parser.add_argument('main_model_file',  
        help='sets path to the MCell4 model')
    return parser


def process_opts(name = 'MCell4'):
    parser = create_argparse(name)
    args = parser.parse_args()
    opts = Options()
    
    if args.seeds and args.args_file:
        sys.exit("Error: only one argument -s/--seeds or -a/--args-file must be specified, not both.")

    if not args.seeds and not args.args_file:
        sys.exit("Error: one of arguments -s/--seeds or -a/--args-file must be specified.")
         
    if args.seeds:
        opts.seeds_str = args.seeds
    
    if args.extra_arg:
        opts.extra_arg = args.extra_arg
    
    if args.args_file:
        if os.path.exists(args.args_file): 
            opts.args_file = args.args_file
        else:
            sys.exit("Error: file " + args.args_file + " does not exist.")
            
        if args.extra_arg:
            sys.exit("Error: argument -x/-extra-args cannot be used with -a/--args-file.")
         
    if args.max_cores:
        opts.max_cores = args.max_cores
        
    if args.main_model_file:
        if os.path.exists(args.main_model_file): 
            opts.main_model_file = args.main_model_file
        else:
            sys.exit("Error: file " + args.main_model_file + " does not exist.")
    else:
        sys.exit("Error: main model file must be specified as a positional argument.")
        
    return opts
        

def check_prerequisites():
    global MCELL_PATH
    MCELL_PATH = os.environ.get('MCELL_PATH', '')
    if MCELL_PATH:
        sys.path.append(os.path.join(MCELL_PATH, 'lib'))
    else:
        sys.exit("Error: system variable MCELL_PATH that is used to find the mcell library was not set.")
    
    if os.name == 'nt':
        ext = '.pyd'
    else:
        ext = '.so'
    
    mcell_so_path = os.path.join(MCELL_PATH, 'lib', 'mcell' + ext)
    
    
    if not os.path.exists(mcell_so_path):
        sys.exit("Could not find library '" + mcell_so_path + ".")
    
            
def generate_seeds(seeds_str):
    seeds_info = seeds_str.split(':')
    if len(seeds_info) != 3 or \
        not seeds_info[0].isdigit() or \
        not seeds_info[1].isdigit() or \
        not seeds_info[2].isdigit():
        sys.exit("Error: invalid seed string, must be in for form min:max:step, given '" + seeds_str + "'.")
    
    res  = []
    for i in range(int(seeds_info[0]), int(seeds_info[1]) + 1, int(seeds_info[2])):
        res.append(i)
    
    return res


def prepare_args(opts):
    if opts.seeds_str:
        seeds = generate_seeds(opts.seeds_str)
        return ['-seed ' + str(i) + ' ' + opts.extra_arg for i in seeds ]
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
    
    args_log = args_str.replace(' ', '_').replace('-', '_').replace('\\', '_').replace('/', '_')
    log_name = os.path.join(
        LOGS_DIR, 
        os.path.splitext(os.path.basename(opts.main_model_file))[0] + '_' + args_log + '.mcell4.log')

    with open(log_name, "w") as f:
        f.write("DIR:" + os.getcwd() + "\n")
        f.write("CMD:" + cmd_str + "\n")
    
    exit_code = 1
    with open(log_name, "a") as f:
        proc = subprocess.Popen(cmd_str, shell=True, cwd=os.getcwd(), stdout=f, stderr=subprocess.STDOUT)
        proc.communicate()
        exit_code = proc.returncode

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
     
    # create logs directory
    if not os.path.exists(LOGS_DIR):
        os.mkdir(LOGS_DIR)
    
    # run the jobs
    pool = multiprocessing.Pool(processes=cpu_count)
    args_str_w_opts = zip(args, itertools.repeat(opts))
    res_codes = pool.map(run_mcell4, args_str_w_opts, 1)

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

    
def my_touch(fname):
    # emulates 'touch', pathlib.Path(FINISHED_MARKER).touch() may not be available
    try:
        if os.path.exists(fname):
            os.utime(fname, None)
        else:
            open(fname, 'a').close()
    except:
        print("Warning: could not 'touch' file " + fname + ".")
    
            
def file_markers_start():    
    if (os.path.exists(FINISHED_MARKER)):
        os.remove(FINISHED_MARKER)
    my_touch(RUNNING_MARKER)


def file_markers_finish():    
    if (os.path.exists(RUNNING_MARKER)):
        os.remove(RUNNING_MARKER)
    my_touch(FINISHED_MARKER)


if __name__ == '__main__':
    file_markers_start()
    
    check_prerequisites()
    
    opts = process_opts()
        
    args = prepare_args(opts)
    
    exit_code = run_mcell4_parallel(opts, args)
    
    file_markers_finish()
    
    sys.exit(exit_code)