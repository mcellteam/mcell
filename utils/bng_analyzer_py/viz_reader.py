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
    
    # the MCell4 API does not have currently a function to parse and convert 
    # BNGL string to the internal Complex representation, 
    # so we will create a BNGL file and load it instead
    with open('tmp.bngl', 'w') as f:
        f.write('begin seed species\n')
        for index, row in counts_df.iterrows():
            f.write('  ' + row['species'] + ' ' + str(row['count']) + '\n')  
        f.write('\nend seed species\n')      

    model = m.Model()
    obj = m.geometry_utils.create_box('box', 1)
    model.load_bngl('tmp.bngl', '', obj)
    
    # the information was generated as seed species so we just read it from 
    # the release sites
    res = []
    for rel_site in model.release_sites:
        res.append((rel_site.complex, rel_site.number_to_release))
    
    return res


def read_dat_file(file_name):
    # returns a list of pairs (mcell.Complex, int), the second item is count 
    
    # load the .dat file as a pandas dataframe
    counts_df = load_counts_from_dat_file(file_name)
    
    # parse the BNGL representations
    complex_counts = parse_bngl_strings_to_complex_representations(counts_df)
    
    return complex_counts
