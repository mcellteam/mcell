# MCell

## Overview

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

## Download Latest Test Builds

These builds are the from the head of this branch and are not guaranteed to be
stable. Use at your own risk.

* [Linux](https://bintray.com/jczech/mcell/download_file?file_path=mcell-linux-gcc.tgz)
* [OSX](https://bintray.com/jczech/mcell/download_file?file_path=mcell-osx-gcc.tgz)

## Build Requirements:

### Ubuntu 16.04:

Run the following commands:

    sudo apt-get update
    sudo apt-get install cmake build-essential bison flex python3-dev swig libboost-all-dev

### Windows

Install the package manager chocolatey from https://chocolatey.org/install.
Follow the instructions listed on that page, which will have you copy a command
into PowerShell as an Administrator.

Once chocolatey is installed, open up a new PowerShell terminal as an
Administrator, and run the following commands:

    choco install -y git
    choco install -y cmake
    choco install -y ninja
    choco install -y msys2
    choco install -y winflexbison
    choco install -y swig
    
Open a new PowerShell Administrator terminal and enter the following:

    msys2

This will open a different MSYS2 terminal (not a PowerShell terminal). In this
MSYS2 terminal, enter the following commands:

    pacman -Syuu
    pacman -S --needed base-devel mingw-w64-i686-toolchain mingw-w64-x86_64-toolchain \
                    git subversion mercurial \
                    mingw-w64-i686-cmake mingw-w64-x86_64-cmake
    pacman -Sy mingw-w64-i686-boost mingw-w64-x86_64-boost

You may have to explicitly add some of the executables to your path.

## Building MCell Executable from Source (OSX, Linux, Windows)

Open a terminal (PowerShell for Windows users) and clone the repo and checkout
the appropriate branch:

    git clone https://github.com/mcellteam/mcell
    cd mcell
    git checkout nfsim_dynamic_meshes_pymcell


### CMake

If this is your first time cloning the repo, you'll want to do this first:

    git submodule init
    git submodule update

To build MCell and pyMCell for Mac or Linux, run the following commands from
the main mcell directory:

    mkdir build
    cd build
    cmake ..
    make

If you're building on Windows with Ninja, change the last two steps to this:

    cmake -G Ninja ..
    ninja

Note: pyMCell does not currently build on Windows.

## Alternative (non-CMake) Method to Build pyMCell:

PyMCell is an experimental MCell-Python library. You can build it using the
traditional CMake method above or this distutils based method, which requires
swig and a newer version of Python 3 (preferably 3.5 or greater). Run the
following command:

    python3 setupy.py build

## How to Test:

### Testing with nutmeg

[nutmeg](https://github.com/mcellteam/nutmeg) is a regression test
framework for MCell. Installation and usage instructions are listed on the
nutmeg project page.

### Testing MCellR

MCellR testing hasn't been incorporated into nutmeg yet, but you can test
MCellR functionality directly after building MCell. Simply run the following
commands (starting at the top level of the "mcell" project directory):

    python3 -m venv mcell_venv
    source mcell_venv/bin/activate
    pip install -r requirements.txt
    cd build
    python mdlr2mdl.py -ni ./fceri_files/fceri.mdlr -o ./fceri_files/fceri -r
