#!/bin/bash

mkdir -p build/build_gperftools
mkdir -p build/install_gperftools

cd build/build_gperftools

../../libs/gperftools/configure --prefix=`pwd`/../install_gperftools --disable-dependency-tracking
make -j 8
make install
