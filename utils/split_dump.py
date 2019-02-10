import os
import sys
import re


with open(sys.argv[1]) as fin:
    line = fin.readline()
    index = 0
    while line:
        if "************ (START)" in line:
            # starting a new file
            type = fin.readline()
            with open("dump_" + string(index).zfill(3) + "_" + type.strip() + ".txt", "w") as fout:
                fout.write(type) 
                line = fin.readline()
                while line and not "************ (END)" in line:
                    four.write(line)
