#!/usr/bin/env python3

"""
Copyright (C) 2019 by
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

import os
import sys

def get_last_line_count(file):
    last_line = ''
    with open(file, 'r') as f:
        for line in f:
            if line:
                last_line = line
                
    return float(last_line.split()[1])


def add_last_line_observables_counts(dir, counts):
    seed_dir = os.path.join(os.getcwd(), dir)
    file_list = os.listdir(seed_dir)
    
    for file in file_list:
        file_path = os.path.join(seed_dir, file)
        if os.path.isfile(file_path) and file.endswith('.dat'):
            cnt = get_last_line_count(file_path)
            observable = os.path.splitext(file)[0]
            if observable.endswith('_MDLString'):
                observable = observable[:-len('_MDLString')] 
            if observable in counts: 
                counts[observable] += cnt
            else:
                counts[observable] = cnt


# just one variant for now, but this could grow int oa more versatile tool
def main():
    
    # go through all seed_xxx directories in the current directory
    dir_list = os.listdir(os.getcwd())
    num_dirs = 0
    counts = {}
    for dir in dir_list:
        if os.path.isdir(dir) and dir.startswith('seed_'):
            add_last_line_observables_counts(dir, counts)
            num_dirs += 1
    
    for name, count in sorted(counts.items()):
        print(name.ljust(20) + ': ' + format(float(count)/num_dirs, '.4f'))


if __name__ == "__main__":
    main()
    