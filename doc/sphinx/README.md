The README for the MCell Quick Reference Guide

The MCell Quick Reference Guide located at http://mcell.org/documentation/qrg
is generated with Sphinx (http://sphinx-doc.org/index.html) from
reStructuredText files. If Sphinx is installed, the html files can be generated
by using the command "make html" in this directory and the PDF can be generated
using "make latexpdf". Here are some instructions for new developers who want
to contribute to the QRG but are unfamiliar with ReST and Sphinx:

* Read this page to develop a basic understanding of ReST syntax:
  http://sphinx-doc.org/rest.html
    * The official "quick start" guide:
      http://docutils.sourceforge.net/docs/user/rst/quickstart.html
    * The official ReST reference guide:
      http://docutils.sourceforge.net/docs/user/rst/quickref.html
    * Sphinx specific syntax: http://sphinx-doc.org/markup/index.html
* Install sphinx.
    * For Debian/Ubuntu machines use this command:
        * apt-get install python-sphinx
* Follow this tutorial, which will briefly explain how Sphinx projects work:
  http://sphinx-doc.org/tutorial.html
* The current implementation of the QRG uses a lot of grid tables which can be
  difficult to edit out of the box with some text editors. Here are some
  suggestions.
    * Vim
        * The vim-rst-tables plugin can be very useful, although it sometimes
          behaves oddly.
        * Consider setting one of the following options
            * set ve=block
            * set ve=all
