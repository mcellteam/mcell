import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from itertools import count

# for now a simple scripts that prints averages of last n values of react_data files
# usage:
# python3 reac_data_analyzer.py N
#     - N is the number of the last observable counts from react_data outputs that will be averaged
# 

# function collectreturns 
def get_last_observable_counts(num_last_samples_to_avg=1):
    counts = {}
    num_seeds = 0
    seed_dirs = os.listdir()
        
    # go through all seed_* directories
    for seed_dir in seed_dirs:
        if not seed_dir.startswith('seed_'):
            continue
        
        # we need the total number of directories to compute an average later
        num_seeds += 1
        
        # go through all *.dat files in the seed directory
        file_list = os.listdir(seed_dir)
        for file in file_list:
            file_path = os.path.join(seed_dir, file)
            
            # we care only about .dat files
            if os.path.isfile(file_path) and file.endswith('.dat'):
                print("Processing " + file_path)
                observable = os.path.splitext(file)[0]
                if observable.endswith('_MDLString'):
                    observable = observable[:-len('_MDLString')]
                
                # read the .dat file into a pandas dataframe
                df = pd.read_csv(file_path, sep=' ', names=['time', 'count'])
                
                #print(df.tail(num_last_samples_to_avg)['count'].mean())
                #print("From " + observable)
                #sys.exit()
                
                # get average of the last N items
                avg_cnt = df.tail(num_last_samples_to_avg)['count'].mean()
                
                # and accumulate the observable count
                if observable in counts:
                    counts[observable] += avg_cnt
                else:
                    counts[observable] = avg_cnt
    
    # compute average of the sums of averages we computed above
    res = {}
    for v,c in sorted(counts.items()):
        res[v] = c / num_seeds
    
    return res
    

# process argument
num_last_samples_to_avg = 1    
if len(sys.argv) == 2:
    num_last_samples_to_avg = int(sys.argv[1])
    
# read all *.dat files in the current directory and return a  
# dictionary observable -> average count
avg_counts = get_last_observable_counts(num_last_samples_to_avg)


for v,c in sorted(avg_counts.items()):
    print(v + ": " + str(c))

