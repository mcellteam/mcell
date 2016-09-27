# MCell

MCell (Monte Carlo Cell) development is supported by the NIGMS-funded
(P41GM103712) National Center for Multiscale Modeling of Biological Systems
(MMBioS).

MCell is a program that uses spatially realistic 3D cellular models and
specialized Monte Carlo algorithms to simulate the movements and reactions of
molecules within and between cells—cellular microphysiology. 

[![Build Status](https://travis-ci.org/mcellteam/mcell.svg?branch=master)](https://travis-ci.org/mcellteam/mcell)
[![Build status](https://ci.appveyor.com/api/projects/status/github/mcellteam/mcell?branch=master&svg=true)](https://ci.appveyor.com/project/jczech/mcell/branch/master)
<a href="https://scan.coverity.com/projects/mcellteam-mcell">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/8521/badge.svg"/>
</a>

## Build Requirements:

### Ubuntu 14.04:

Run the following command:

    sudo apt-get install cmake build-essential bison flex

## Building from Source:

### CMake

To build MCell for Macs or Linux, run the following commands from the main
mcell directory:

    mkdir build
    cd build
    cmake ..
    make

### Autoconf and Automake (Deprecated)

The old build system is still available and can be used by issuing the 
following commands:

    cd ./src
    ./bootstrap
    cd ..
    mkdir build
    cd build
    ../src/configure CC=gcc CFLAGS='-O2 -Wall' 
    make

You only need to bootstrap (first three steps) when starting from a fresh
branch or checkout. Depending on your needs, you may have to change the
build options slightly.

See the [Windows
Development](https://github.com/mcellteam/mcell/wiki/Windows-Development) page
on the github wiki for information about building MCell on Windows.

## How To Test:

[nutmeg](https://github.com/mcellteam/nutmeg) is a regression test
framework for MCell. Installation and usage instructions are listed on the
nutmeg project page.
