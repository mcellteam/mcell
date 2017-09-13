import subprocess
import os

def nfsim():
    subprocess.call(['git', 'clone', 'https://github.com/RuleWorld/nfsim.git'])
    os.chdir('nfsim')
    subprocess.call(['git', 'checkout', 'nfsim_lib_shared_ptr'])
    subprocess.call(['git', 'pull'])

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
    subprocess.call(['git', 'clone', 'https://github.com/jjtapia/nfsimCInterface'])
    os.chdir('nfsimCInterface')
    subprocess.call(['git', 'pull'])

    try:
        os.mkdir('lib')
    except OSError:
        pass

    os.chdir('lib')
    subprocess.call(['ln','-s',os.path.join('..','..','nfsim','lib','libNFsim.so')])
    os.chdir('..')

    subprocess.call(['ln','-s',os.path.join('..','nfsim','include')])


    try:
        os.mkdir('build')
    except OSError:
        pass
    os.chdir('build')
    subprocess.call(['cmake', '..'])
    subprocess.call(['make'])
    os.chdir('..')
    os.chdir('..')

def create_links():
    os.chdir('..')
    try:
        os.mkdir('lib')
    except OSError:
        pass
    os.chdir('lib')
    subprocess.call(['ln','-s',os.path.join('..','build','nfsim','lib','libNFsim.so')])
    subprocess.call(['ln','-s',os.path.join('..','build','nfsimCInterface','build','libnfsim_c.so')])
    os.chdir('..')
    os.chdir('build')

if __name__ == "__main__":
    nfsim()
    nfsim_lib()
    create_links()
