#!/bin/bash

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
  # /Applications/Blender-2.79-CellBlender/blender.app/Contents/Resources/2.79/scripts/addons/cellblender/extensions/mcell/
  REL_BLENDER_PATH=$MCELL_PATH/../../../../../../../../../
else
  REL_BLENDER_PATH=$MCELL_PATH/../../../../../../
fi

# get absolute path to the REL_BLENDER_PATH to deal with links
ABS_BLENDER_PATH=`cd $REL_BLENDER_PATH; pwd`

if [ "$(uname)" == "Darwin" ]; then
  VIZ_MCELL_SCRIPT=$ABS_BLENDER_PATH/blender.app/Contents/Resources/2.79/scripts/addons/cellblender/developer_utilities/mol_viz_scripts/viz_mcell_run.py
else
  VIZ_MCELL_SCRIPT=$ABS_BLENDER_PATH/2.79/scripts/addons/cellblender//developer_utilities/mol_viz_scripts/viz_mcell_run.py
fi


if [ "$(uname)" == "Darwin" -o "$(uname)" == "Linux" ]; then 
  $ABS_BLENDER_PATH/my_blender -P $VIZ_MCELL_SCRIPT -- $1
else
  # Windows
  $ABS_BLENDER_PATH/blender.exe -P $VIZ_MCELL_SCRIPT -- $1
fi
