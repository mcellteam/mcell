#!/bin/bash

for F in *.py; do
    echo "*** Running $F ***"
    ./$F
    EC=$?
    if [ "$EC" != "0" ]; then
        echo "ERROR"
        exit 1
    fi 
done

echo "PASSED"