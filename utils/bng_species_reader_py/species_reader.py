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


def read_species_file(file_name):
    # load the .species file as a list of pairs (mcell.Complex, int), the second item is count
    res = []
    with open(file_name, 'r') as f:
        for line in f:
            if not line.strip():
                continue 
            
            items = line.split()
            assert len(items) == 2, "Invalid input file contents " + line
            
            # constructor m.Complex parses the BNGL representaion into 
            # a MCell4 API representation 
            # (see https://cnl.salk.edu/~ahusar/mcell4_documentation/generated/subsystem.html#complex
            cplx = m.Complex(items[0])
            res.append((cplx, str(items[1])))
    
    return res