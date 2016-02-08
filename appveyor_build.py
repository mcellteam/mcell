#!/usr/bin/env python
# Build the project on AppVeyor.

import os
import shutil
from subprocess import check_call

mcell_src = './src'
build_command = ['gcc.exe', '-mconsole', '-std=c99', '-O3', '-fno-schedule-insns2', '-o', 'mcell.exe', '*.c']
files = ["config.h", "version.h", "mdllex.c", "mdlparse.h", "mdlparse.c"]

for f in files:
    shutil.move("./appveyor_windows/%s" % f, mcell_src)
os.chdir(mcell_src)
check_call(build_command)
