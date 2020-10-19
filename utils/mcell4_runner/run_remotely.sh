#!/bin/bash

# Example: 
# first 'cd' to the directory with your model then run:
# run_remotely.sh $USER@<machine>:/home/$USER/tmp/my_model <path to mcell4_runner>/mcell4_runner.py . "model.py -s 1:10:1"


USAGE=\
"$1 - remote machine with user and working directory e.g. user@a.server.com:/home/user/tmp\n"\
"$2 - local path to mcell4_runner.py script (including the script name)\n"\
"$3 - directory with model\n"\
"$4 - arguments for mcell4_runner.py script (enclosed in quotes)\n"  

if [ "$#" != "4" ]; then
    echo "Error: expected 4 arguments:\n $USAGE"
    exit 1
fi

DST=$1
RUNNER_PY=$2
SRC=$3
ARGS=$4

# copy model directory
rsync --delete -arv $SRC $DST || exit 1

# copy script 
rsync --delete -arv $2 $DST || exit 1

URL=`echo $DST | cut -d: -f1`
DIR=`echo $DST | cut -d: -f2`
CMD="conda activate py35; cd $DIR; nohup python mcell4_runner.py $4"
ssh $URL $CMD
