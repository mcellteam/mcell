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
import os
import pandas as pd

MCELL_PATH = os.environ.get('MCELL_PATH', '')
if MCELL_PATH:
    sys.path.append(os.path.join(MCELL_PATH, 'lib'))
else:
    print("Error: variable MCELL_PATH that is used to find the mcell library was not set.")
    sys.exit(1)

import mcell as m

# utility module to load ASCII viz_output .dat file from 
import viz_reader


# arbitrary values are used here
WEIGHTS = {
    'ampar_tarp': 1.5,
    'psd95': 2.5,
    'syngap': 3.5,
}
    
    
def analyze_dat_file(file_name):
    # read the .dat file and parse complexes to the internal MCell
    # representation 
    complex_counts = viz_reader.read_dat_file(file_name)
    
    #print(complex_counts[0][0])
    #print(complex_counts[0][1])
    
    # process the read data
    for (complex, count) in complex_counts:
        
        weight = 0.0
        
        # iterate over elementary molecules from which the 
        # complex is composed
        for mi in complex.elementary_molecules:
            name = mi.elementary_molecule_type.name
            if name not in WEIGHTS:
                print("Error: unknown molecular weight of " + name)
                sys.exit(1)
                
            weight += WEIGHTS[name]
            
        print("----------------------------------")
        print(complex.to_bngl_str())
        print("")
        print("weight: " + str(weight) + ", copy nr.: " + str(count))
    
    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Expected .dat file name as argument")
        sys.exit(1)
        
    analyze_dat_file(sys.argv[1])
    