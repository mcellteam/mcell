#!/usr/bin/env python

# commands to run setup of swig files 
	# swig -python mcell.i  
	# vim mcell_swig_setup.py 

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


mcell_module = Extension('_mcellSwig',
                           sources=['mcellSwig_wrap.c', 'mcell_misc.c', 'mcell_objects.c', 'mcell_react_out.c', 'mcell_reactions.c', 'mcell_release.c', 'mcell_species.c', 'mcell_viz.c', 'mcell_surfclass.c', 'react_output.c', 'strfunc.c']
                           )

setup (name = 'mcellSwig',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [mcell_module],
       py_modules = ["mcellSwig"],
       )
