#!/usr/bin/env python3

"""
Auxiliary script to generate constant definitiosn from a JSON FILE
"""

import json
import sys


def process(json_object, keys, values):
    if type(json_object) == dict:
        for key, value in json_object.items():
            keys.append(key)
            process(value, keys, values)
    elif type(json_object) == list:
        for v in json_object:
            process(v, keys, values)
    else:
        values.append(str(json_object))
        
def print_as_const(lst, prefix):
    
    for s in lst:
        if not s.isidentifier():
            continue
        
        print("const char* const " + prefix + s.upper() + " = \"" + s + "\";")

if __name__ == '__main__':
    with open(sys.argv[1], "r") as f:
        json_object = json.load(f)
    
    # using list so that we get related keys and values close to each other
    keys = []
    values = []
    process(json_object, keys, values)

    print_as_const(keys, "KEY_")
    print_as_const(values, "VALUE_")