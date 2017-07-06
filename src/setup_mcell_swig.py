#!/usr/bin/env python

# commands to run setup of swig files 
# swig -python pymcell.i  

"""
setup.py file for pyMCell
"""

from distutils.core import setup, Extension


mcell_module = Extension(
    '_pymcell',
    sources=[
        'pymcell_wrap.c',
        'mcell_react_out.c',
        'mcell_misc.c',
        'util.c',
        'argparse.c',
        'mcell_objects.c',
        'mcell_reactions.c',
        'mcell_release.c',
        'mcell_species.c',
        'mcell_viz.c',
        'mcell_surfclass.c',
        'strfunc.c',
        'mcell_init.c',
        'chkpt.c',
        'init.c',
        'react_output.c',
        'version_info.c',
        'sched_util.c',
        'wall_util.c',
        'mem_util.c',
        'react_trig.c',
        'react_cond.c',
        'count_util.c',
        'vol_util.c',
        'diffuse.c',
        'vector.c',
        'diffuse_util.c',
        'sym_table.c',
        'triangle_overlap.c',
        'logging.c',
        'grid_util.c',
        'react_outc.c',
        'react_util.c',
        'react_outc_trimol.c',
        'diffuse_trimol.c',
        'isaac64.c',
        'rng.c',
        'viz_output.c',
        'mcell_run.c',
        'volume_output.c'])                           

setup (name = 'pymcell',
       version = '0.1',
       author      = "The MCell team",
       description = """Python bindings to libmcell""",
       author_email = "mcell-devel@salk.edu",
       ext_modules = [mcell_module],
       license = 'GPL v2',
       py_modules = ["pymcell"],
       )
