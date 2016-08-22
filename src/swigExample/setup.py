#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


example_module = Extension('_example',
                           sources=['example_wrap.c', 'mcell_misc.c', 'mcell_objects.c', 'mcell_react_out.c', 'mcell_reactions.c', 'mcell_release.c', 'mcell_species.c', 'mcell_viz.c', 'mcell_surfclass.c', 'strfunc.c']
                           )

setup (name = 'example',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [example_module],
       py_modules = ["example"],
       )
