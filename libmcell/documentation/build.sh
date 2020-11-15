#!/bin/bash
sphinx-build -M html . build/

if [ "$1" == "-u" ]; then
    rsync -r build/html/* huxley.snl.salk.edu:/home/ahusar/public_html/mcell4_documentation
fi
