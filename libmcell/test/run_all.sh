#!/bin/bash

for F in *.py; do
    if [[ "$F" == *"error"* ]]; then
        NEGATIVE=1
    else
        NEGATIVE=0
    fi 
    
    echo "*** Running $F (negative=$NEGATIVE) ***"
    ./$F
    EC=$?
    if [ "$NEGATIVE" == "0" -a "$EC" != "0" ]; then
        echo "ERROR (test should pass)"
        exit 1
    fi 

    if [ "$NEGATIVE" == "1" -a "$EC" == "0" ]; then
        echo "ERROR (test should not pass)"
        exit 1
    fi 
done

echo "PASSED"