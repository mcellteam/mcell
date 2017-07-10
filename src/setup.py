#!/usr/bin/env python

# commands to run setup of swig files 
# swig -python -py3 -ltypemaps pymcell.i  
# python setup.py build_ext --inplace

"""
setup.py file for pyMCell
"""

from distutils.core import setup, Extension
from distutils.command.build import build


# class CustomBuild(build):
#     sub_commands = [
#         ('build_ext', build.has_ext_modules),
#         ('build_py', build.has_pure_modules),
#         ('build_clib', build.has_c_libraries),
#         ('build_scripts', build.has_scripts),
#     ]



mcell_module = Extension(
    '_pymcell',
    sources=[
        'pymcell_wrap.c',
        'argparse.c',
        'chkpt.c',
        'count_util.c',
        'diffuse.c',
        'diffuse_trimol.c',
        'isaac64.c',
        'diffuse_util.c',
        'dyngeom.c',
        'dyngeom_lex.c',
        'dyngeom_parse_extras.c',
        'dyngeom_yacc.c',
        'grid_util.c',
        'react_outc.c',
        'init.c',
        'logging.c',
        'mcell_init.c',
        'mcell_misc.c',
        'mcell_objects.c',
        'mcell_react_out.c',
        'mcell_reactions.c',
        'mcell_release.c',
        'mcell_run.c',
        'mcell_species.c',
        'mcell_surfclass.c',
        'mcell_viz.c',
        'mcell_dyngeom.c',
        'mem_util.c',
        # 'pymcell.i',
        'react_cond.c',
        'react_outc_trimol.c',
        'react_output.c',
        'react_trig.c',
        'react_util.c',
        'rng.c',
        'sched_util.c',
        'strfunc.c',
        'sym_table.c',
        'triangle_overlap.c',
        'util.c',
        'vector.c',
        'version_info.c',
        'viz_output.c',
        'vol_util.c',
        'volume_output.c',
        'wall_util.c',
        ],
    swig_opts=['-ltypemaps', '-py3'],
    extra_compile_args=['-O2'])                           

setup (name = 'pymcell',
       version = '0.1',
       author = "The MCell team",
       description = """Python bindings to libmcell""",
       author_email = "mcell-devel@salk.edu",
       ext_modules = [mcell_module],
       license = 'GPL v2',
       py_modules = ["pymcell2"],
       # cmdclass={'build': CustomBuild},
       )
