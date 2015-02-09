MCell
=====

MCell: Monte Carlo Simulator of Cellular Microphysiology

[![Build Status](https://travis-ci.org/mcellteam/mcell.svg?branch=master)](https://travis-ci.org/mcellteam/mcell)


MCell3 Build Requirements:
--------------------------

1. flex newer than 2.5.6, due to the 'reentrant' option.  Extensive testing
     has been done using Flex 2.5.33.


How To Bootstrap:
-----------------

1. When starting from a fresh branch or checkout, cd into ./src and do the following: 

        ./bootstrap


How To Build:
-------------

1. To build for a given platform do something like this:

        mkdir my_build_dir
        cd my_build_dir
        /path/to/configure CC=gcc CFLAGS='-O3 -march=core2 -Wall -g' 
        make


How To Test:
------------

[nutmeg](https://github.com/haskelladdict/nutmeg) is a regression test
framework for MCell. Installation and usage instructions are listed on the
nutmeg project page.
