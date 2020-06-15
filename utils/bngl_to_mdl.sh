#!/bin/bash

#this link is required in cellblender
#ln -s '/home/ahusar/src4_ref/mcell_tools/work/build_mcell' '/home/ahusar/src4_ref/cellblender/extensions/mcell'

python3 '/home/ahusar/src4_ref/cellblender/bng/bngl_to_data_model.py' $1 data_model.json || exit 1

python3 '/home/ahusar/src4_ref/cellblender/mdl/data_model_to_mdl.py' data_model.json Scene.main.mdl -fail-on-error || exit 1

sed -i 's/CELLBLENDER/ASCII/g' Scene.viz_output.mdl

echo
echo "*** Conversion passed ***"
echo