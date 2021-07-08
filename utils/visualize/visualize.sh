#!/bin/bash


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

if [ "$MCELL_PATH" == "" ]; then
  echo "System variable MCELL_PATH is not set, don't know where to find CellBlender"
  exit 1
fi

if [ "$1" == "" ]; then
  echo "Expecing one argument that is the path to the viz output directory, e.g. viz_data/seed_00001/"
  exit 1
fi

if [ "$(uname)" == "Darwin" ]; then
  # path on MacOS
  # /Applications/Blender-2.93-CellBlender/blender.app/Contents/Resources/2.93/scripts/addons/cellblender/extensions/mcell/
  REL_BLENDER_PATH=$MCELL_PATH/../../../../../../../../../
else
  REL_BLENDER_PATH=$MCELL_PATH/../../../../../../
fi

# get absolute path to the REL_BLENDER_PATH to deal with links
ABS_BLENDER_PATH=`cd $REL_BLENDER_PATH; pwd`

if [ "$(uname)" == "Darwin" ]; then
  VIZ_MCELL_SCRIPT=$ABS_BLENDER_PATH/blender.app/Contents/Resources/2.93/scripts/addons/cellblender/developer_utilities/mol_viz_scripts/viz_mcell_run.py
else
  VIZ_MCELL_SCRIPT=$ABS_BLENDER_PATH/2.93/scripts/addons/cellblender/developer_utilities/mol_viz_scripts/viz_mcell_run.py
fi


if [ "$(uname)" == "Darwin" -o "$(uname)" == "Linux" ]; then 
  $ABS_BLENDER_PATH/my_blender -P $VIZ_MCELL_SCRIPT -- $1
else
  # Windows
  $ABS_BLENDER_PATH/blender.exe -P $VIZ_MCELL_SCRIPT -- $1
fi
