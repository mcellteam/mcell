#!/bin/bash
# this is an auxiliary script used e.g. when API changes
find ./ -type f -exec sed -i "s/double/double/g" {} \;