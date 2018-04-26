#!/usr/bin/env python

# command to build pymcell
# python setup.py build
# command to create pymcell tarball
# python setup.py sdist

"""
setup.py file for pyMCell
"""

import shutil
import os
import re
import sys
import platform
import subprocess
from setuptools import setup, Extension, find_packages
from distutils.command.build import build
from distutils.command.sdist import sdist
from distutils.version import LooseVersion
from setuptools.command.build_ext import build_ext

# def disallow_python2():
#     if sys.version_info[0] == 2:
#         sys.exit("Sorry, Python 2 is not supported.")

# Code for cmake extensions from (c45488d  on Jun 10, 2016), with minor
# modifications:
# https://github.com/pybind/cmake_example/blob/master/setup.py

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir='', extra_cxx_flags=[], *args, **kw):
        Extension.__init__(self, name, sources=[], *args, **kw)
        self.sourcedir = os.path.abspath(sourcedir)
        self.extra_cxx_flags = extra_cxx_flags


class CMakeBuild(build_ext):
    def __init__(self, dist, *args, **kw):
        build_ext.__init__(self, dist, *args, **kw)

    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        subprocess.call(['ln', '-s', self.build_lib, 'lib'])

        os.chdir('nfsimCInterface')
        subprocess.call(['ln', '-s', os.path.join('..', self.build_lib), 'lib'])
        subprocess.call(['ln', '-s', os.path.join('..', 'nfsim', 'include')])
        os.chdir('..')

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        print("building", ext.name)
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        self.build_temp_lib = self.build_temp + '/' + ext.name
        if not os.path.exists(self.build_temp_lib):
            os.makedirs(self.build_temp_lib)

        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp_lib, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp_lib)

class CustomBuild(build):
    def run(self):
        # TODO: some of these lines can be trimmed down eventually
        if not os.path.exists("build"):
            os.makedirs("build")
        shutil.copy("./appveyor_windows/config.h", "./src")
        shutil.copy("./appveyor_windows/version.h", "./src")
        shutil.copy("./src/pymcell_helpers.py", ".")
        shutil.copy("./src/data_model_import.py", ".")

        self.run_command('build_ext')
        shutil.copy("./src/pymcell.py", ".")
        build.run(self)
        print("ran successfully!")


class CustomSDist(sdist):
    def run(self):
        # disallow_python2()
        sdist.run(self)

# mcell_module = Extension(
#     '_pymcell',
#     include_dirs=['./include'],
#     libraries=['nfsim_c', 'NFsim'],
#     library_dirs=['./build/lib'],
#     sources=[
#         './src/argparse.c',
#         './src/chkpt.c',
#         './src/count_util.c',
#         './src/diffuse.c',
#         './src/diffuse_trimol.c',
#         './src/diffuse_util.c',
#         './src/dyngeom.c',
#         './src/dyngeom_lex.c',
#         './src/dyngeom_parse_extras.c',
#         './src/dyngeom_yacc.c',
#         './src/grid_util.c',
#         './src/hashmap.c',
#         './src/init.c',
#         './src/isaac64.c',
#         './src/logging.c',
#         './src/mcell_dyngeom.c',
#         './src/mcell_init.c',
#         './src/mcell_misc.c',
#         './src/mcell_objects.c',
#         './src/mcell_react_out.c',
#         './src/mcell_reactions.c',
#         './src/mcell_release.c',
#         './src/mcell_run.c',
#         './src/mcell_species.c',
#         './src/mcell_surfclass.c',
#         './src/mcell_viz.c',
#         './src/mem_util.c',
#         './src/nfsim_func.c',
#         './src/pymcell.i',
#         './src/react_cond.c',
#         './src/react_outc.c',
#         './src/react_outc_nfsim.c',
#         './src/react_outc_trimol.c',
#         './src/react_output.c',
#         './src/react_trig.c',
#         './src/react_trig_nfsim.c',
#         './src/react_util.c',
#         './src/react_util_nfsim.c',
#         './src/rng.c',
#         './src/sched_util.c',
#         './src/strfunc.c',
#         './src/sym_table.c',
#         './src/triangle_overlap.c',
#         './src/util.c',
#         './src/vector.c',
#         './src/version_info.c',
#         './src/viz_output.c',
#         './src/vol_util.c',
#         './src/volume_output.c',
#         './src/wall_util.c',
#         ],
#     swig_opts=['-ltypemaps', '-py3'],
#     extra_compile_args=['-O2'])

ext_modules = [CMakeExtension('nfsim', sourcedir='nfsim'),
               CMakeExtension('nfsimCInterface', sourcedir='nfsimCInterface'),
               CMakeExtension('_pymcell')]

setup(name='pymcell',
      py_modules=['pymcell'],
      version='0.1',
      author="The MCell team",
      python_requires='>=3',
      description="""Python bindings to libmcell""",
      author_email="mcell-devel@salk.edu",
      url="https://github.com/mcellteam/mcell",
      download_url="https://github.com/mcellteam/mcell/archive/pymcell_0.1.tar.gz",
      ext_modules=ext_modules,
      license='GPL v2',
      packages=find_packages(),
      cmdclass={'build_ext': CMakeBuild, 'build': CustomBuild, 'sdist': CustomSDist},
      )
