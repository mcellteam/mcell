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
import pandas as pd
import numpy as np

# set this to change the precision of the values printed, it must be a string
DIGITS = '6'

def write_bngl_function(input_fname, param):
    # load the input file, make sure that the values are floats
    df = pd.read_csv(
        input_fname, delim_whitespace=True, comment='#', 
        names=['time', 'value'], dtype={'time':np.float64, 'value':np.float64})
                                        
    #print(df['time'])
    #print(df['value'])
    
    if df['time'][0] != 0:
        print("First time must be 0")
        sys.exit(1)
    
    ind = "    "
    for i in range(1, len(df)):
        print(
            ind +
            "(if((" + param + "<" + str.format('{0:.' + DIGITS + '}', df['time'][i]) + "), " + 
            "(" + str.format('{0:.' + DIGITS + '}', df['value'][i - 1]) + "), \\") 

    # last value
    print(ind + "  (" + str.format('{0:.' + DIGITS + '}', df['value'].iloc[-1]) + ") \\") 

    print(ind, end='') 
    for i in range(1, len(df)):
        print("))", end = '')
    print()

def main():
    if len(sys.argv) != 3:
        print("Expecting input file as 1st argument and input variable name (usually time) as 2nd argument")
        sys.exit(1)
    
    write_bngl_function(sys.argv[1], sys.argv[2])
    

if __name__ == "__main__":
    main()
    