import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from itertools import count

# for now a simple scripts that prints acerages of last values of react_data files

def get_last_observable_counts():
    counts = {}
    num_seeds = 0
    seed_dirs = os.listdir()
        
    for seed_dir in seed_dirs:
        if not seed_dir.startswith('seed_'):
            continue
        
        num_seeds += 1
        file_list = os.listdir(seed_dir)
        for file in file_list:
            file_path = os.path.join(seed_dir, file)
            if os.path.isfile(file_path) and file.endswith('.dat'):
                print("Processing " + file_path)
                observable = os.path.splitext(file)[0]
                if observable.endswith('_MDLString'):
                    observable = observable[:-len('_MDLString')]
                    
                df = pd.read_csv(file_path, sep=' ', names=['time', 'count'])
                #print(df.tail(1).iloc[0]['count'])
                #print("From " + observable)
                #sys.exit()
                cnt = df.tail(1).iloc[0]['count']
                if observable in counts:
                    counts[observable] += cnt
                else:
                    counts[observable] = cnt
    
    return counts, num_seeds            
    
counts, num_seeds = get_last_observable_counts()
print(counts)

for v,c in sorted(counts.items()):
    print(v + ":" + str(c / num_seeds))

