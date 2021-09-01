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
from io import StringIO

def main():
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        sys.exit("Expecting exactly one argument that is the .gdat " 
                 "file and optionally comma-separated column indices.")
    
    
    
    # load all .dat files in directory passed as the first argument
    file = sys.argv[1]
    if not os.path.exists(file):
        sys.exit("File " + file + " does not exist.")
    
    wholefile = open(file).read()
    
    df = pd.read_csv(StringIO(wholefile[1:]), index_col=0, skipinitialspace=True, delim_whitespace=True)
    
    
    #df = df.set_index('time')
    print(df)
    
    col_names = pd.DataFrame
    if len(sys.argv) == 3:
        sel = [int(i) for i in sys.argv[2].split(',') ]
        col_names = df.columns[sel]
    print(col_names)
    
    if not col_names.empty:
        df_sel = df[col_names]
    else:
        df_sel = df
    
    df_sel.plot(kind='line')
    plt.show()


if __name__ == '__main__':
    main()
    
    