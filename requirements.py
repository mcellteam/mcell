import subprocess
import os
import shutil
import sys


def nfsim():
    os.chdir('nfsim')

    try:
        os.mkdir('lib')
    except OSError:
        pass
    os.chdir('lib')
    subprocess.call(['cmake', '..'])
    subprocess.call(['make'])
    os.chdir('..')
    os.chdir('..')


def nfsim_lib():
    os.chdir('nfsimCInterface')

    try:
        os.mkdir('lib')
    except OSError:
        pass

    os.chdir('lib')
    extension = "so"
    if (sys.platform == 'darwin'):
        extension = "dylib"
    subprocess.call(['ln', '-s', os.path.join('..', '..', 'nfsim', 'lib', 'libNFsim.{}'.format(extension))])
    os.chdir('..')

    subprocess.call(['ln', '-s', os.path.join('..', 'nfsim', 'include')])

    try:
        os.mkdir('build')
    except OSError:
        pass
    os.chdir('build')
    subprocess.call(['cmake', '..'])
    subprocess.call(['make'])
    os.chdir('..')
    os.chdir('..')


def copy_library_files():
    try:
        os.mkdir('build/lib')
    except OSError:
        pass
    os.chdir('build/lib')
    extension = "so"
    if (sys.platform == 'darwin'):
        extension = "dylib"
    shutil.copy(os.path.join('..', '..', 'nfsim', 'lib', 'libNFsim.{}'.format(extension)), ".")
    shutil.copy(os.path.join('..', '..', 'nfsimCInterface', 'build', 'libnfsim_c.{}'.format(extension)), ".")


if __name__ == "__main__":
    subprocess.call(['git', 'submodule', 'init'])
    subprocess.call(['git', 'submodule', 'update'])
    os.chdir("..")
    nfsim()
    nfsim_lib()
    copy_library_files()
