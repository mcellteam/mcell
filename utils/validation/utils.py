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
import subprocess
import shutil
from threading import Timer
from subprocess import Popen, PIPE
          

def get_cwd_no_link():
    # get current working directory even though we are in a symlinked directory
    # the shell argument must be set to True in this case
    cwd = Popen(['pwd'], stdout=PIPE, shell=True).communicate()[0].strip()
    return cwd.decode('ascii')
          

def kill_proc(proc, f, timeout_is_fatal):
    proc.kill()
    f.write("Terminated after timeout")
    if timeout_is_fatal:
        sys.exit(1)

    
def print_file(file_name):
    with open(file_name, "r") as fin:
        for line in fin:
            print(line)


def run_with_ascii_output(cmd, cwd):
    # does not return exit code, neither provides timeout
    return Popen(cmd, cwd=cwd, stdout=PIPE).communicate()[0].strip().decode('ascii')
    


def execute(cmd, cwd, timeout_sec, timeout_is_fatal, outfile, shell=False):
    if shell:
        # for shell=True, the command must be a single string
        cmd = str.join(" ", cmd)
    
    proc = Popen(cmd, shell=shell, cwd=cwd, stdout=outfile, stderr=subprocess.STDOUT)
    timer = Timer(timeout_sec, kill_proc, [proc, outfile, timeout_is_fatal])
    try:
        timer.start()
        exit_code = proc.wait()
    finally:
        timer.cancel()
        
    return exit_code


# can be simplified by using subprocess.run from Python 3.5
def run(
        cmd, 
        cwd=os.getcwd(),
        fout_name="", 
        append_path_to_output=False, 
        print_redirected_output=False, 
        timeout_sec=60,
        timeout_is_fatal = True, 
        verbose=True,
        shell=False 
        ):
    if verbose:
        log("    Executing: '" + str.join(" ", cmd) + "' " + str(cmd) + " in '" + cwd + "'")

    if fout_name:
        if append_path_to_output:
            full_fout_path = os.path.join(cwd, fout_name)
        else:
            full_fout_path = fout_name
    
        with open(full_fout_path, "w") as f:
            f.write("cwd: " + cwd + "\n")
            f.write(str.join(" ", cmd) + "\n")  # first item is the command being executed
            
            # run the actual command
            exit_code = execute(cmd, cwd, timeout_sec, timeout_is_fatal, f, shell=shell)

        if (print_redirected_output):
            print_file(full_fout_path)
            
    else:
        exit_code = execute(cmd, cwd, timeout_sec, timeout_is_fatal, sys.stdout, shell=shell)

    if verbose:
        log("Exit code: " + str(exit_code))
    return exit_code


def log(msg):
    print("* " + msg)
    sys.stdout.flush()


def warning(msg):
    print("* Warning: " + msg)
    sys.stdout.flush()


def fatal_error(msg):
    print("* Error: " + msg)
    sys.stdout.flush()
    sys.exit(1)    
   

def check_ec(ec, cmd):
    if ec != 0:
        cmd_str = str.join(" ", cmd)
        fatal_error("Error: command '" + cmd_str + "' failed, terminating.")
        