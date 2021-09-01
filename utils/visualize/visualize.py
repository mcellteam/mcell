# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
# 
# For more information, please refer to [http://unlicense.org] 

import sys
import os
import platform
import subprocess

MCELL_PATH = os.environ.get('MCELL_PATH', '')
if MCELL_PATH:
    lib_path = os.path.join(MCELL_PATH, 'lib')
    if os.path.exists(os.path.join(lib_path, 'mcell.so')) or \
        os.path.exists(os.path.join(lib_path, 'mcell.pyd')):
        sys.path.append(lib_path)
    else:
        print("Error: Python module mcell.so or mcell.pyd was not found in "
              "directory '" + lib_path + "' constructed from system variable "
              "MCELL_PATH.")
        sys.exit(1)
else:
    print("Error: system variable MCELL_PATH that is used to find the mcell "
          "library was not set.")
    sys.exit(1)
    
    
if len(sys.argv) != 2:
    sys.exit("Expecing one argument that is the path to the viz output directory, e.g. viz_data/seed_00001/")

    
REL_BLENDER_PATH = None


if platform.system() == 'Darwin':
    REL_BLENDER_PATH = os.path.join(MCELL_PATH, '..', '..', '..', '..', '..', '..', '..', '..', '..')
else:
    REL_BLENDER_PATH = os.path.join(MCELL_PATH, '..', '..', '..', '..', '..', '..')
    
ABS_BLENDER_PATH = os.path.abspath(REL_BLENDER_PATH)


if platform.system() == 'Darwin':
    VIZ_MCELL_SCRIPT = \
        os.path.join(ABS_BLENDER_PATH, 'blender.app', 'Contents', 'Resources', '2.93', 'scripts', 'addons', 'cellblender', 'developer_utilities', 'mol_viz_scripts', 'viz_mcell_run.py')
else:
    VIZ_MCELL_SCRIPT = \
        os.path.join(ABS_BLENDER_PATH, '2.93', 'scripts', 'addons', 'cellblender', 'developer_utilities', 'mol_viz_scripts', 'viz_mcell_run.py')


if 'Windows' in platform.system():
    CMD = os.path.join(ABS_BLENDER_PATH, 'blender.exe')    
else:
    CMD = 'bash ' + os.path.join(ABS_BLENDER_PATH, 'my_blender')

CMD += ' -P ' + VIZ_MCELL_SCRIPT + ' -- ' + sys.argv[1]

subprocess.run(CMD, shell=True)
   