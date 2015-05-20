MCell
=====

MCell: Monte Carlo Simulator of Cellular Microphysiology

[![Build Status](https://travis-ci.org/mcellteam/mcell.svg?branch=master)](https://travis-ci.org/mcellteam/mcell)


MCell3 Build Requirements:
--------------------------

flex newer than 2.5.6, due to the 'reentrant' option. Extensive testing has
been done using Flex 2.5.33.


How To Bootstrap:
-----------------

When starting from a fresh branch or checkout, cd into ./src and do the
following: 

        ./bootstrap


How To Build:
-------------

To build MCell for Macs or Linux, run the following commands from the main
mcell directory:

        mkdir build
        cd build
        ../src/configure CC=gcc CFLAGS='-O2 -Wall' 
        make

Depending on your needs, you may have to change the build options slightly.

See the [Windows
Development](https://github.com/mcellteam/mcell/wiki/Windows-Development) page
on the github wiki for information about building MCell on Windows.

How To Test:
------------

[nutmeg](https://github.com/haskelladdict/nutmeg) is a regression test
framework for MCell. Installation and usage instructions are listed on the
nutmeg project page.
