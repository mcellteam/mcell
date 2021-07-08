#!/usr/bin/env python3

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
    