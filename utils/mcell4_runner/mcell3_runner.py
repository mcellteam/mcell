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
import pathlib
import argparse
import itertools
import multiprocessing
import re
import subprocess

from mcell4_runner import *
    
def run_mcell3(args_str_w_opts):
    args_str = args_str_w_opts[0]
    opts = args_str_w_opts[1]
    
    if os.name == 'nt':
        ext = '.exe'
    else:
        ext = ''
    
    mcell_path = os.environ.get('MCELL_PATH', '')
    cmd_str = os.path.join(mcell_path, 'mcell' + ext) + ' ' + opts.main_model_file + ' '
    cmd_str += args_str
    
    print("Running " + cmd_str)
    
    args_log = args_str.replace(' ', '_').replace('-', '_')
    log_name = os.path.splitext(opts.main_model_file)[0] + '_' + args_log + '.mcell3.log'
    
    exit_code = 1
    with open(log_name, "w") as f:
        proc = subprocess.Popen(cmd_str, shell=True, cwd=os.getcwd(), stdout=f, stderr=subprocess.STDOUT)
        proc.communicate()
        exit_code = proc.returncode

    with open(log_name, "a") as f:
        f.write("DIR:" + os.getcwd() + "\n")
        f.write("CMD:" + cmd_str + "\n")
    
    if exit_code != 0:
        print("MCell3 failed, see '" + os.path.join(os.getcwd(), log_name) + "'.")
        return exit_code
    else:
        return 0

    
def run_mcell3_parallel(opts, args):
    
    # set up the parallel task pool to use all available processors ot the 
    # maximum specified
    if opts.max_cores:
        cpu_count = int(opts.max_cores)
    else:
        cpu_count = multiprocessing.cpu_count()
     
    pool = multiprocessing.Pool(processes=cpu_count)
    
    # run the jobs
    args_str_w_opts = zip(args, itertools.repeat(opts))
    res_codes = pool.map(run_mcell3, args_str_w_opts, 1)

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


if __name__ == '__main__':
    file_markers_start()
    
    check_prerequisites()
    
    opts = process_opts('MCell3')
        
    args = prepare_args(opts)
    
    exit_code = run_mcell3_parallel(opts, args)
    
    file_markers_finish()
    
    sys.exit(exit_code)