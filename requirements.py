import subprocess
import os
import shutil
import sys

def compile_it():
    if (sys.platform == 'win32'):
        subprocess.call(['cmake', '-G', 'Ninja', '..'])
        subprocess.call(['ninja'])
    else:
        subprocess.call(['cmake', '..'])
        subprocess.call(['make'])


def get_extension():
    if (sys.platform == 'linux') or (sys.platform == 'linux2'):
        extension = "so"
    elif (sys.platform == 'darwin'):
        extension = "dylib"
    elif (sys.platform == 'win32'):
        extension = "dll"
    else:
        raise Exception("Unexpected platform: {0}".format(sys.platform))
    return extension


def nfsim():
    os.chdir('nfsim')

    try:
        os.mkdir('lib')
    except OSError:
        pass
    os.chdir('lib')
    compile_it()
    os.chdir('..')
    os.chdir('..')


def nfsim_lib():
    os.chdir('nfsimCInterface')

    try:
        os.mkdir('lib')
    except OSError:
        pass

    os.chdir('lib')
    extension = get_extension()

    if (sys.platform == 'win32'):
        subprocess.call(['cmd', '/c', 'mklink', '.\libNFsim.{0}'.format(extension), os.path.join('..', '..', 'nfsim', 'lib', 'libNFsim.{0}'.format(extension))])
        os.chdir('..')
        subprocess.call(['cmd', '/c', 'mklink', '/D', '.\include', os.path.join('..', 'nfsim', 'include')])
    else:
        subprocess.call(['ln', '-s', os.path.join('..', '..', 'nfsim', 'lib', 'libNFsim.{0}'.format(extension))])
        os.chdir('..')
        subprocess.call(['ln', '-s', os.path.join('..', 'nfsim', 'include')])


    try:
        os.mkdir('build')
    except OSError:
        pass
    os.chdir('build')
    compile_it()
    os.chdir('..')
    os.chdir('..')


def copy_library_files():
    try:
        os.mkdir('build/lib')
    except OSError:
        pass
    os.chdir('build')

    extension = get_extension()

    # XXX: For some reason, Windows will build MCell with the files in
    # ./build/lib but will fail at runtime if they aren't in the same directory
    # as the executable itself (i.e. the libs must be in ./build at runtime)
    if (sys.platform == 'win32'):
        shutil.copy(os.path.join('..', 'nfsim', 'lib', 'libNFsim.{0}'.format(extension)), ".")
        shutil.copy(os.path.join('..', 'nfsimCInterface', 'build', 'libnfsim_c.{0}'.format(extension)), ".")

    shutil.copy(os.path.join('..', 'nfsim', 'lib', 'libNFsim.{0}'.format(extension)), "lib")
    shutil.copy(os.path.join('..', 'nfsimCInterface', 'build', 'libnfsim_c.{0}'.format(extension)), "lib")


if __name__ == "__main__":
    subprocess.call(['git', 'submodule', 'init'])
    subprocess.call(['git', 'submodule', 'update'])
    os.chdir("..")
    nfsim()
    nfsim_lib()
    copy_library_files()
