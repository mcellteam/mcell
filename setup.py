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

def disallow_python2():
    if sys.version_info[0] == 2:
        sys.exit("Sorry, Python 2 is not supported.")

# Code for cmake extensions from (c45488d  on Jun 10, 2016), with minor
# modifications:
# https://github.com/pybind/cmake_example/blob/master/setup.py

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir='', sources=[], *args, **kw):
        Extension.__init__(self, name, sources=sources, *args, **kw)
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def __init__(self, dist, *args, **kw):
        build_ext.__init__(self, dist, *args, **kw)

    def run(self):
        if subprocess.call(["git", "branch"], stderr=subprocess.STDOUT, stdout=open(os.devnull, 'w')) == 0:
            subprocess.call(['git', 'submodule', 'update', '--init'])
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        if not os.path.exists(self.build_lib):
            os.makedirs(self.build_lib)

        ext = self.extensions[0]
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        os.chdir('nfsimCInterface')
        subprocess.call(['ln', '-s', os.path.join('..', 'nfsim', 'include')])
        subprocess.call(['ln', '-s', os.path.join('..', extdir), 'lib'])
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
        disallow_python2()
        self.run_command('build_ext')
        build.run(self)

class CustomSDist(sdist):
    def run(self):
        disallow_python2()
        # clone submodules if this is a git repo:
        if subprocess.call(["git", "branch"], stderr=subprocess.STDOUT, stdout=open(os.devnull, 'w')) == 0:
            subprocess.call(['git', 'submodule', 'update', '--init'])
        sdist.run(self)


# NOTE: Might be better to rename _pymcell to just mcell, but I do not want
# to break backwards compatibility.
ext_modules = [CMakeExtension('nfsim', sourcedir='nfsim'),
               CMakeExtension('nfsimCInterface', sourcedir='nfsimCInterface'),
               CMakeExtension('_pymcell', sourcedir='mcell')]

setup(name='pymcell',
      packages=['pymcell'],
      version='0.1',
      author="The MCell team",
      description="""Python bindings to libmcell""",
      author_email="mcell-devel@salk.edu",
      url="https://github.com/mcellteam/mcell",
      download_url="https://github.com/mcellteam/mcell/archive/pymcell_0.1.tar.gz",
      ext_modules=ext_modules,
      license='GPL v2',
      cmdclass={'build_ext': CMakeBuild, 'build':CustomBuild, 'sdist': CustomSDist},
      provides=['pymcell']
      )
