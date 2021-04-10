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


def load_counts_from_dat_file(file_name):
    # using pandas read_cvs because it is usually faster 
    # we could also read the file line by line and count the number of occurences 
    df = pd.read_csv(
        file_name, 
        sep=' ', 
        names=['name', 'id', 'x', 'y', 'z', 'nx', 'ny', 'nz'])

    # create a dataframe with counts per unique species
    counts = df['name'].value_counts().rename_axis('species').reset_index(name='count')
    
    # need to sort the data to get consistent results
    counts = counts.sort_values(['count','species'], ascending=(False, True))
    return counts


def parse_bngl_strings_to_complex_representations(counts_df):
    # returns a list of pairs (mcell.Complex, int), the second item is count 
    res = []
    for index, row in counts_df.iterrows():
        # constructor m.Complex parses the BNGL representaion into 
        # a MCell4 API representation 
        # (see https://cnl.salk.edu/~ahusar/mcell4_documentation/generated/subsystem.html#complex
        cplx = m.Complex(row['species'])
        res.append((cplx, row['count']))
    
    return res


def read_dat_file(file_name):
    # returns a list of pairs (mcell.Complex, int), the second item is count 
    
    # load the .dat file as a pandas dataframe
    counts_df = load_counts_from_dat_file(file_name)
    
    # parse the BNGL representations
    complex_counts = parse_bngl_strings_to_complex_representations(counts_df)
    
    return complex_counts
