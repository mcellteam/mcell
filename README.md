# MCell

## Overview

MCell (Monte Carlo Cell) development is supported by the NIGMS-funded
(P41GM103712) National Center for Multiscale Modeling of Biological Systems
(MMBioS).

MCell is a program that uses spatially realistic 3D cellular models and
specialized Monte Carlo algorithms to simulate the movements and reactions of
molecules within and between cells—cellular microphysiology. 

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
    choco install -y python3
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

Add CMake and MinGW64 to your [Windows PATH Environment
Variable](https://helpdeskgeek.com/windows-10/add-windows-path-environment-variable/).
If you have installed these dependencies with Chocolatey and your top level
drive is named `C:`, you can append (or prepend) the following paths into your
Windows PATH (Search "Control Panel"; View by small icons; System; Advanced System Settings;
              Environment Variables; User Variables; Path; Edit):

    C:\tools\msys64\mingw64\bin
    C:\Program Files\CMake\bin

## Building MCell Executable from Source (OSX, Linux, Windows)

Open a terminal (non-Administrator PowerShell for Windows users), clone the
mcell_tools repo and run:

  git clone https://github.com/mcellteam/mcell\_tools.git
  cd mcell\_tools
  python run.py

This will clone all the required repositories and run build of all the components.
Running 'python run.py --help' shows other options.

To build just MCell, run these commands after you ran the 'python run.py'.

  cd mcell\_tools/work/build\_mcell
  cmake ../../../mcell -DCMAKE_BUILD_TYPE=Release
  make 

### Testing MCell

  cd mcell\_tests
  python run_tests.py

Running 'python run_tests.py --help' shows other options. 
