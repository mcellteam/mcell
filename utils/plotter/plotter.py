"""
Copyright (C) 2020 by
The Salk Institute for Biological Studies and
Pittsburgh Supercomputing Center, Carnegie Mellon University

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

For the complete terms of the GNU General Public License, please see this URL:
http://www.gnu.org/licenses/gpl-2.0.html
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import argparse


class Options:
    def __init__(self):
        self.mcell3_dir = None
        self.mcell4_dir = None
        self.bng_dir = None
        self.single_bng_run = False

        
def create_argparse():
    parser = argparse.ArgumentParser(description='MCell4 Runner')
    parser.add_argument('-m4', '--mcell4', type=str, help='mcell4 react_data directory')
    parser.add_argument('-m3', '--mcell3', type=str, help='mcell3 react_data directory')
    parser.add_argument('-b', '--bng', type=str, help='bionetgen directory')
    parser.add_argument('-s', '--single_bng_run', action='store_true', help='the bionetgen directory contains only a single .gdat file')
    return parser


def process_opts():
    parser = create_argparse()
    args = parser.parse_args()
    opts = Options()
    
    if args.mcell4:
        opts.mcell4_dir = args.mcell4 

    if args.mcell3:
        opts.mcell3_dir = args.mcell3 

    if args.bng:
        opts.bng_dir = args.bng 

    if args.single_bng_run:
        opts.single_bng_run = args.single_bng_run 

    return opts


def get_mcell_observables_counts(dir):
    counts = {}
    seed_dirs = os.listdir(dir)
    
    for seed_dir in seed_dirs:
        if not seed_dir.startswith('seed_'):
            continue
        
        file_list = os.listdir(os.path.join(dir, seed_dir))
        for file in file_list:
            file_path = os.path.join(dir, seed_dir, file)
            if os.path.isfile(file_path) and file.endswith('.dat'):
                observable = os.path.splitext(file)[0]
                if observable.endswith('_MDLString'):
                    observable = observable[:-len('_MDLString')]
                
                if observable not in counts:
                    index = 0
                else:
                    index = counts[observable].shape[1] - 1 
                
                col_name = 'count' + str(index)
                df = pd.read_csv(file_path, sep=' ', names=['time', col_name])

                if observable not in counts:
                    counts[observable] = df
                else:
                    # add new column
                    counts[observable][col_name] = df[col_name]
                
    return counts 


def get_bng_observables_counts(file, counts):
    if not os.path.exists(file):
        print("Expected file " + file + " not found, skipping it")
        return
    
    with open(file, 'r') as f:
        first_line = f.readline()
        header = first_line.split()[1:]
    df = pd.read_csv(file, delim_whitespace=True, comment='#', names=header)
    return df
    

def process_nsfim_gdat_file(full_dir, counts):
    df = get_bng_observables_counts(os.path.join(full_dir, 'test.gdat'), counts)
    
    # transform into separate dataframes based on observable 
    for i in range(1, df.shape[1]):
        observable = df.columns[i] 
        if observable not in counts:
            col_name = 'count0'  
            # select time and the current observable
            counts[observable] = pd.DataFrame()           
            counts[observable]['time'] = df.iloc[:, 0]
            counts[observable][col_name] = df.iloc[:, i]
        else:
            col_name = 'count' + str(counts[observable].shape[1] - 1)            
            counts[observable][col_name] = df.iloc[:, i]             
                
                
def get_nfsim_observables_counts(opts):
    single_bng_run = opts.single_bng_run
    dir = opts.bng_dir
    counts = {}
    
    if not single_bng_run:
        nf_dirs = os.listdir(dir)
        
        for nf_dir in nf_dirs:
            full_dir = os.path.join(dir, nf_dir)
            if not nf_dir.startswith('nf_') or not os.path.isdir(full_dir):
                continue
            process_nsfim_gdat_file(full_dir, counts)
    else:
        process_nsfim_gdat_file(dir, counts)

    return counts 


def main():
    
    opts = process_opts()
    
    counts = []
    if opts.mcell4_dir:
        if os.path.exists(opts.mcell4_dir):
            print("Reading MCell data from " + opts.mcell4_dir)
            counts.append(get_mcell_observables_counts(opts.mcell4_dir))
        else:
            print("Directory " + opts.mcell4_dir + " with MCell4 data not found, ignored")
            sys.exit(1)
    else:
        counts.append({})

    if opts.mcell3_dir:
        if os.path.exists(opts.mcell3_dir):
            print("Reading MCell data from " + opts.mcell3_dir)
            counts.append(get_mcell_observables_counts(opts.mcell3_dir))
        else:
            print("Error: directory " + opts.mcell3_dir + " with MCell3 data not found, ignored")
            sys.exit(1)
    else:
        counts.append({})
    
    # get_nfsim_observables_counts may return an empty dict
    if opts.bng_dir:
        if os.path.exists(opts.bng_dir):
            print("Reading BNG data from " + opts.bng_dir)
            counts.append(get_nfsim_observables_counts(opts))
        else:
            print("Error: directory " + opts.bng_dir + " with BNG data not found, ignored")
            sys.exit(1)
    else:
        counts.append({})
            
    names = ['MCell4', 'MCell3R', 'BNG']
    clrs = ['b', 'g', 'r'] 

    all_observables = set(counts[0].keys())
    all_observables = all_observables.union(set(counts[1].keys()))
    all_observables = all_observables.union(set(counts[2].keys()))
    
    for obs in sorted(all_observables): 
        print("Processing observable " + obs)
        
        fig,ax = plt.subplots()
        ax.set_title(obs)
        
        legend_names = []
        for i in range(len(counts)):
            if obs not in counts[i]:
                continue
            
            data = counts[i][obs]
            
            df = pd.DataFrame()           
            df['time'] = data.iloc[:, 0]
            df['means'] = data.iloc[:, 1:].mean(axis=1)
            df['mean_minus_std'] = df['means'] - data.iloc[:, 1:].std(axis=1)
            df['mean_plus_std'] = df['means'] + data.iloc[:, 1:].std(axis=1)
    
            # free collected data to decrease memory consumption        
            del data
                    
            ax.plot(df['time'], df['means'], label=obs + "1", c=clrs[i])
            ax.fill_between(
                df['time'], 
                df['mean_minus_std'], df['mean_plus_std'],
                alpha=0.2, facecolor=clrs[i])
    
            legend_names.append(names[i])
    
        plt.legend(legend_names)
        
        plt.savefig(obs + '.png', dpi=600)
        plt.close(fig)


if __name__ == '__main__':
    main()
    
    