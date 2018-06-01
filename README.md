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

It's easiest to get all the dependencies using [chocolatey](https://chocolatey.org/). Once that's installed, open up a powershell terminal, and run the following commands:

    choco install cmake
    choco install ninja
    choco install msys2
    choco install winflexbison
    
Once you have msys2 installed, run the msys2 command (which will open a different non-Powershell terminal) and execute the following commands:

    pacman -Syuu
    pacman -S --needed base-devel mingw-w64-i686-toolchain mingw-w64-x86_64-toolchain \
                    git subversion mercurial \
                    mingw-w64-i686-cmake mingw-w64-x86_64-cmake
    pacman -Sy mingw-w64-i686-boost mingw-w64-x86_64-boost

You may have to explicitly add some of the executables to your path.

## Building from Source:

To build MCell, run the following commands from the main
mcell directory:

    mkdir build
    cd build
    cmake ..
    make
    
If you're building on Windows with Ninja, change the last two steps to this:

    cmake -G Ninja ..
    ninja


See the [Windows
Development](https://github.com/mcellteam/mcell/wiki/Windows-Development) page
on the github wiki for information about building MCell on Windows.

## How To Test:

[nutmeg](https://github.com/mcellteam/nutmeg) is a regression test
framework for MCell. Installation and usage instructions are listed on the
nutmeg project page.
