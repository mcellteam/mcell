#!/bin/bash 

g++ \
	-isystem /nadata/cnl/home/ahusar/src4/mcell/libs/ \
	-I`pwd`/../../src4 -I`pwd`/../../src -I`pwd`/../../include -I`pwd`/../../../mcell_tools/work/build_mcell/deps/ \
	-DNOSWIG=1 \
	-O0 -g3 -Wall -shared -std=c++17 -fPIC `python3 -m pybind11 --includes` -o mcell.so \
	../api/bindings.cpp ../generated/gen_species.cpp ../generated/gen_release_site.cpp
	

