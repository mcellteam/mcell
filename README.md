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

### Ubuntu 16.04:

Run the following command:

    sudo apt-get install cmake build-essential bison flex python3-dev

## Building MCell executable from Source:

### CMake

If this is your first time cloning the repo, you'll want to do this first:

    git submodule init
    git submodule udpate

To build MCell for Macs or Linux, run the following commands from the main
mcell directory:

    mkdir build
    cd build
    cmake ..
    make

See the [Windows
Development](https://github.com/mcellteam/mcell/wiki/Windows-Development) page
on the github wiki for information about building MCell on Windows.

## Building pyMCell (MCell Python library):

You will need swig and some version of Python 3 (preferably 3.5) Run the
following command:

  python setupy.py build

## How To Test:

[nutmeg](https://github.com/mcellteam/nutmeg) is a regression test
framework for MCell. Installation and usage instructions are listed on the
nutmeg project page.
