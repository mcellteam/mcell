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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import glob

def load_dat_files(dat_dir):
    counts = {}
    
    res = pd.DataFrame
    
    # read all .dat files
    dat_files = sorted(glob.glob(os.path.join(dat_dir, "*.dat"))) 
    for file in dat_files:
        observable = os.path.splitext(os.path.basename(file))[0]
        df = pd.read_csv(file, sep=' ', index_col='time', names=['time', observable])

        # use the first data frame as basis and the join with new observables 
        # to create a single data frame        
        if res.empty:
            res = df
        else:
            res = res.join(df)

    return res


def main():
    if len(sys.argv) != 2:
        sys.exit("Expecting exactly one argument that is the path to directory with .dat files.")
    
    # load all .dat files in directory passed as the first argument
    dat_dir = sys.argv[1]
    if not os.path.exists(dat_dir):
        sys.exit("Directory " + dat_dir + " does not exist.")
    
    
    df = load_dat_files(dat_dir)
    
    df.plot(kind='line')
    plt.show()


if __name__ == '__main__':
    main()
    
    