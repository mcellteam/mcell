#!/usr/bin/env python

# command to build pymcell
# python setup.py build
# command to create pymcell tarball
# python setup.py sdist

"""
setup.py file for pyMCell
"""

from distutils.core import setup, Extension
from distutils.command.build import build
from distutils.command.sdist import sdist
import shutil
import sys


def disallow_python2():
    if sys.version_info[0] == 2:
        sys.exit("Sorry, Python 2 is not supported.")


class CustomBuild(build):
    def run(self):
        disallow_python2()
        shutil.copy("./appveyor_windows/config.h", "./src")
        shutil.copy("./appveyor_windows/version.h", "./src")
        shutil.copy("./src/pymcell_helpers.py", ".")
        shutil.copy("./src/data_model_import.py", ".")
        self.run_command('build_ext')
        shutil.copy("./src/pymcell.py", ".")
        build.run(self)


class CustomSDist(sdist):
    def run(self):
        disallow_python2()
        sdist.run(self)


mcell_module = Extension(
    '_pymcell',
    sources=[
        './src/argparse.c',
        './src/chkpt.c',
        './src/count_util.c',
        './src/diffuse.c',
        './src/diffuse_trimol.c',
        './src/isaac64.c',
        './src/diffuse_util.c',
        './src/dyngeom.c',
        './src/dyngeom_lex.c',
        './src/dyngeom_parse_extras.c',
        './src/dyngeom_yacc.c',
        './src/grid_util.c',
        './src/react_outc.c',
        './src/init.c',
        './src/logging.c',
        './src/mcell_init.c',
        './src/mcell_misc.c',
        './src/mcell_objects.c',
        './src/mcell_react_out.c',
        './src/mcell_reactions.c',
        './src/mcell_release.c',
        './src/mcell_run.c',
        './src/mcell_species.c',
        './src/mcell_surfclass.c',
        './src/mcell_viz.c',
        './src/mcell_dyngeom.c',
        './src/mem_util.c',
        './src/pymcell.i',
        './src/react_cond.c',
        './src/react_outc_trimol.c',
        './src/react_output.c',
        './src/react_trig.c',
        './src/react_util.c',
        './src/rng.c',
        './src/sched_util.c',
        './src/strfunc.c',
        './src/sym_table.c',
        './src/triangle_overlap.c',
        './src/util.c',
        './src/vector.c',
        './src/version_info.c',
        './src/viz_output.c',
        './src/vol_util.c',
        './src/volume_output.c',
        './src/wall_util.c',
        ],
    swig_opts=['-ltypemaps', '-py3'],
    extra_compile_args=['-O2'])

setup(name='pymcell',
      version='0.2',
      author="The MCell team",
      description="""Python bindings to libmcell""",
      author_email="mcell-devel@salk.edu",
      url="https://github.com/mcellteam/mcell",
      download_url="https://github.com/mcellteam/mcell/archive/pymcell_0.2.tar.gz",
      ext_modules=[mcell_module],
      license='GPL v2',
      py_modules=["pymcell"],
      cmdclass={'build': CustomBuild, 'sdist': CustomSDist},
      )
