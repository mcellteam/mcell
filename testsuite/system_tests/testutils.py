###################################################################################
#                                                                                 #
# Copyright (C) 2006-2013 by                                                      #
# The Salk Institute for Biological Studies and                                   #
# Pittsburgh Supercomputing Center, Carnegie Mellon University                    #
#                                                                                 #
# This program is free software; you can redistribute it and/or                   #
# modify it under the terms of the GNU General Public License                     #
# as published by the Free Software Foundation; either version 2                  #
# of the License, or (at your option) any later version.                          #
#                                                                                 #
# This program is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                  #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   #
# GNU General Public License for more details.                                    #
#                                                                                 #
# You should have received a copy of the GNU General Public License               #
# along with this program; if not, write to the Free Software                     #
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. #
#                                                                                 #
###################################################################################

"""Basic MCell testing utilities.

Miscellaneous useful utilities for testing are placed in this module.  In
particular, this module contains classes and functions related to invoking and
testing generic MCell runs.  Most in-depth output-testing utilities have been
placed in other modules (such as viz_output or reaction_output).

Most new tests will use McellTest (or a subclass of MCellTest) to handle
starting the MCell run in a clean directory, running it, checking a handful of
small invariant criteria that should hold for all MCell runs (for instance,
that it did not segfault.)

"""

import unittest
import stat
import os
import sys
import ConfigParser
import string
#import types
import random
import re
import shutil
import subprocess
import windows_util


def check_assertions_enabled():
    """Check if assertions are enabled."""

    assertions_enabled = 0

    try:
        assert 1 == 2
    except AssertionError:
        assertions_enabled = 1

    return assertions_enabled


def require_assertions_enabled():
    """Complain and exit if assertions are not enabled."""
    if check_assertions_enabled() == 0:
        print "Please turn off Python optimization (remove the -O flag from the command line)."
        print "The optimization feature turns off assertions, which are vital to the functioning"
        print "of the test suite."
        raise Exception("Assertions are disabled, testsuite cannot run")

require_assertions_enabled()


def get_caller_filename(lvl=1):
    """Get the filename of the source file containing the function which calls
    this function."""
    frame = sys._getframe()
    for i in range(0, lvl):
        frame = frame.f_back
    return frame.f_code.co_filename


def safe_concat(l1, l2):
    """safe_concat(l1, l2) -> combined list.

    Utility to safely concatenate two iterables, even if one or the other may
    be None, or the two iterables have different types.  If they have different
    types, the return type will be a tuple.

    """
    if l1 is None:
        return l2
    elif l2 is None:
        return l1
    elif type(l1) is not type(l2):
        return tuple(l1) + tuple(l2)
    else:
        return l1 + l2


def get_output_dir():
    """Get the preferred output directory."""
    g = globals()
    if not "test_output_dir" in g:
        global test_output_dir
        test_output_dir = "./test_results"
    return g["test_output_dir"]


def crange(l, h, s=1):
    """ Utility to generate a closed range.

    crange(start, stop[, step]) -> list of integers.

    Works like 'range', but is closed on both ends -- i.e. includes both
    endpoints.

    """
    if s > 0:
        return range(l, h + 1, s)
    else:
        return range(l, h - 1, s)


def cleandir(directory):
    """Ensure that the speficied directory exists, and that it is empty.

    Make sure a directory exists, and contains no garbage. Used to ensure a
    clean working directory for each run of the testsuite.

    """
    if os.path.exists(directory):
        shutil.rmtree(directory)
    #try:
    os.mkdir(directory)
    #except:
    #pass
    #for root, dirs, files in os.walk(directory, topdown=False):
    #    map(lambda f: os.unlink(os.path.join(root, f)), files)
    #    map(lambda d: os.rmdir(os.path.join(root, d)), dirs)


def assertFileExists(fname):
    """Raise an error if the specified file does not exist. """
    try:
        os.stat(fname)
    except:
        assert False, "Expected file '%s' was not created" % fname


class RequireFileExists:
    def __init__(self, name):
        self.name = name

    def check(self):
        assertFileExists(self.name)


def assertFileNotExists(fname):
    """Raise an error if the specified file exists."""
    try:
        os.stat(fname)
    except:
        return
    assert False, "Specifically unexpected file '%s' was created" % fname


class RequireFileNotExists:
    def __init__(self, name):
        self.name = name

    def check(self):
        assertFileNotExists(self.name)


def assertFileEmpty(fname):
    """Give an error if a file doesn't exist or isn't empty."""
    try:
        sb = os.stat(fname)
    except:
        assert False, "Expected empty file '%s' was not created" % fname
    assert sb.st_size == 0, "Expected file '%s' should be empty but has length %d" % (fname, sb.st_size)


class RequireFileEmpty:
    def __init__(self, name):
        self.name = name

    def check(self):
        assertFileEmpty(self.name)


def assertFileNonempty(fname, expect=None):
    """Give an error if a file doesn't exist or is empty."""
    try:
        sb = os.stat(fname)
    except:
        assert False, "Expected file '%s' was not created" % fname
    assert sb.st_size != 0, "Expected file '%s' shouldn't be empty but is" % fname
    if expect is not None:
        assert sb.st_size == expect, "Expected file '%s' is the wrong size (%d bytes instead of %d bytes)" % (fname, sb.st_size, expect)


class RequireFileNonempty:
    def __init__(self, name, expect=None):
        self.name = name
        self.expect = expect

    def check(self):
        assertFileNonempty(self.name, self.expect)


def assertFileEquals(fname, contents):
    """ Give an error if a file doesn't exist or if its contents differ from
    the contents passed in.

    """
    try:
        got_contents = open(fname).read()
    except:
        assert False, "Expected file '%s' was not created" % fname
    assert got_contents == contents, "Expected file '%s' had incorrect contents" % fname


class RequireFileEquals:
    def __init__(self, name, contents):
        self.name = name
        self.contents = contents

    def check(self):
        assertFileEquals(self.name, self.contents)


def assertFileMatches(fname, regex, expectMinMatches=1, expectMaxMatches=sys.maxint):
    """Give an error if a file doesn't exist or if its contents match (or do
    not match) the regular expression passed in.

    """
    try:
        got_contents = open(fname).read()
    except:
        assert False, "Expected file '%s' was not created" % fname
    matches = re.findall(regex, got_contents)
    if len(matches) == 1:
        plural = ""
    else:
        plural = "es"
    assert len(matches) >= expectMinMatches, "Expected file '%s' had incorrect contents (found %d match%s for regex '%s', expected at least %d)" % (fname, len(matches), plural, regex, expectMinMatches)
    assert len(matches) <= expectMaxMatches, "Expected file '%s' had incorrect contents (found %d match%s for regex '%s', expected at most %d)" % (fname, len(matches), plural, regex, expectMaxMatches)


class RequireFileMatches:
    def __init__(self, name, regex, expectMinMatches=1, expectMaxMatches=sys.maxint):
        self.name = name
        self.regex = regex
        self.expectMinMatches = expectMinMatches
        self.expectMaxMatches = expectMaxMatches

    def check(self):
        assertFileMatches(self.name, self.regex, self.expectMinMatches, self.expectMaxMatches)


def assertFileSymlink(fname, target=None):
    """ Give an error if a file isn't a symlink, or (optionally) if the
    destination of the symlink isn't as specified.

    """
    assert os.path.exists(fname), "Expected symlink '%s' was not created" % fname
    assert os.path.xislink(fname), "Expected symlink '%s' is not a symlink" % fname
    if target is not None:
        got_target = os.xreadlink(fname)
        assert got_target == target, "Expected symlink '%s' should point to '%s', but instead points to '%s'" % (fname, target, got_target)


class RequireFileSymlink:
    def __init__(self, name, target=None):
        self.name = name
        self.target = target

    def check(self):
        assertFileSymlink(self.name, self.target)


def assertFileDir(fname):
    """Give an error if a file isn't a directory."""
    try:
        sb = os.lstat(fname)
    except:
        assert False, "Expected directory '%s' was not created" % fname
    assert stat.S_ISDIR(sb.st_mode), "Expected directory '%s' is not a directory" % fname


class RequireFileDir:
    def __init__(self, name, target=None):
        self.name = name

    def check(self):
        assertFileDir(self.name)


class test_config(object):
    """Small wrapper class to handle loading test suite configuration"""
    def __init__(self, filepath="./test.cfg"):
        dict = {
            "mcellpath": "mcell"
        }
        self.config = ConfigParser.ConfigParser(dict)
        self.filepath = filepath

        try:
            self.config.readfp(open(self.filepath))
        except:
            print "ERROR: invalid file path '%s' to the configuration file" % self.filepath
            sys.exit(0)

    def get(self, sect, val):
        if self.config.has_section(sect):
            return self.config.get(sect, val)
        else:
            return self.config.get(ConfigParser.DEFAULTSECT, val)


class test_run_context(object):
    """Class to handle running command-line executables as PyUnit tests."""
    testidx = 1

    def __init__(self, cmd, args):
        self.command = cmd
        self.args = args
        self.testidx = test_run_context.testidx
        self.check_stdout = 0
        self.check_stderr = 0
        self.check_stdin = 0
        self.expect_exitcode = 0
        self.cleaned = False
        test_run_context.testidx += 1

    def set_check_std_handles(self, i, o, e):
        self.check_stdin = i
        self.check_stdout = o
        self.check_stderr = e

    def set_expected_exit_code(self, ec):
        self.expect_exitcode = ec

    def invoke(self, sandboxdir):
        curdir = os.getcwd()
        testpath = '%s/test-%04d' % (sandboxdir, self.testidx)
        if not self.cleaned:
            cleandir(testpath)
            self.cleaned = True
        os.chdir(testpath)
        try:
            try:
                self.__run()
                self.__check_results()
            except AssertionError, e:
                e.args = ("%s: %s" % (testpath, e.args[0])),
                raise
        finally:
            os.chdir(curdir)

    def check_stdout_valid(self, fullpath):
        assertFileEmpty(fullpath)

    def check_stderr_valid(self, fullpath):
        assertFileEmpty(fullpath)

    def check_output_files(self):
        pass

    def __exit_code(self):
        if hasattr(os, 'WIFEXITED') and callable(os.WIFEXITED) and os.WIFEXITED(self.got_exitcode):
            return os.WEXITSTATUS(self.got_exitcode)
        else:
            return self.got_exitcode

    def __check_results(self):
        assert self.got_exitcode >= 0, "Process died due to signal %d" % -self.got_exitcode
        assert self.__exit_code() == self.expect_exitcode, "Expected exit code %d, got exit code %d" % (self.expect_exitcode, self.__exit_code())
        if self.check_stdout:
            self.check_stdout_valid(os.path.join(os.getcwd(), "stdout"))
        if self.check_stderr:
            self.check_stderr_valid(os.path.join(os.getcwd(), "stderr"))
        self.check_output_files()

    def __run(self):
        try:
            f = open("cmdline.txt", "w", 0644)
            f.write("executable: ")
            f.write(self.command)
            f.write('\n')
            f.write("full cmdline: ")
            f.write(string.join(self.args))
            f.write('\n')
        finally:
            f.close()
        new_stdout = os.open("./stdout", os.O_CREAT | os.O_WRONLY | os.O_EXCL, 0644)
        new_stderr = os.open("./stderr", os.O_CREAT | os.O_WRONLY | os.O_EXCL, 0644)
        self.got_exitcode = subprocess.call(self.args, executable=self.command, stdin=None, stdout=new_stdout, stderr=new_stderr)
        os.close(new_stdout)
        os.close(new_stderr)


class McellTest(test_run_context):
    """Utility base class for MCell tests.

    Specialized class for running MCell jobs as PyUnit tests

    This class will build up the command-line, choosing a random seed,
    looking up the MCell executable in a configuration file, redirect output
    to log/err files using -logfile and -errfile, and do a handful of
    global validations (criteria which must be true for any MCell run).

    """
    rand = random.Random()

    def __init__(self, cat, f, args=[]):
        """Create a new MCell test runner.

        'cat' determines the section of the config file to search for the
            MCell executable path.
        'file' is the name of the MDL file to use.
        'args' is a list of arguments other than the MDL file to pass.

        """
        path = os.path.dirname(os.path.realpath(get_caller_filename(2)))
        try:
            os.stat(os.path.join(path, f))
        except:
            assert False, "Didn't find MDL file '%s' in the expected location (%s)." % (f, path)
        mcell = McellTest.config.get(cat, "mcellpath")
        try:
            os.stat(mcell)
        except:
            print "ERROR: path '%s' to mcell executable in configuration file is invalid" % mcell
            sys.exit(0)
        real_args = [mcell]
        real_args.extend(["-seed", str(McellTest.rand.randint(0, 50000))])
        real_args.extend(["-logfile", "realout"])
        real_args.extend(["-errfile", "realerr"])
        real_args.extend(args)
        real_args.append(os.path.join(path, f))
        test_run_context.__init__(self, mcell, real_args)
        self.set_check_std_handles(1, 1, 1)
        self.set_expected_exit_code(0)
        self.extra_checks = []

    ## Add an extra check to be performed upon completion of the run
    def add_extra_check(self, chk):
        """Add an extra post-run check to this test. The check must be an
        object that has a 'check' method on it, which contains appropriate
        'assert' statements to perform the check.

        Typically, this is used as follows:
                o = McellTest(...)
                o.add_extra_check(RequireFileEquals(filename, contents))
                o.invoke()

        This is especially useful for types of checks that cannot be added via
        the utility methods given below.

        """
        self.extra_checks.append(chk)

    ## Add one or more expected output files (either a string or an iterable)
    def add_exist_file(self, e):
        """Utility to add a file or list of files whose existence is a
        requisite post-condition for a successful run of the given MDL file.

        """
        #if type(e) == types.StringType:
        if isinstance(e, str):
            self.extra_checks.append(RequireFileExists(e))
        else:
            self.extra_checks.extend([RequireFileExists(x) for x in e])
        return self

    ## Add one or more expected empty output files (either a string or an
    ## iterable)
    def add_empty_file(self, e):
        """Utility to add a file or list of files whose existence (and
        emptiness) is a requisite post-condition for a successful run of the
        given MDL file.

        """
        #if type(e) == types.StringType:
        if isinstance(e, str):
            self.extra_checks.append(RequireFileEmpty(e))
        else:
            self.extra_checks.extend([RequireFileEmpty(x) for x in e])
        return self

    ## Add one or more expected non-empty output files (either a string or an
    ## iterable)
    def add_nonempty_file(self, e, expected_size=None):
        """Utility to add a file or list of files whose existence (and
        non-emptiness) is a requisite post-condition for a successful run of
        the given MDL file.

        """
        #if type(e) == types.StringType:
        if isinstance(e, str):
            self.extra_checks.append(RequireFileNonempty(e, expected_size))
        else:
            self.extra_checks.extend([RequireFileNonempty(x, expected_size) for x in e])
        return self

    ## Add one or more expected constant output files (either fname and cnt must
    ## be strings, or they must be iterables with the same number of items)
    def add_constant_file(self, fname, cnt):
        """Utility to add a file or list of files whose existence is a
        requisite post-condition for a successful run of the given MDL file,
        and whose contents must exactly match a specified string.

        """
        #if type(fname) == types.StringType:
        if isinstance(fname, str):
            self.extra_checks.append(RequireFileEquals(fname, cnt))
        else:
            self.extra_checks.extend([RequireFileEquals(x, cnt) for x in fname])
        return self

    ## Add one or more expected symlinks, with optional target. fname may be a
    ## string or an iterable returning a string. target may be a string, or if
    ## fname is an iterable, an iterable returning a string.
    def add_symlink(self, fname, target=None):
        """Utility to add a file or list of filenames whose existence is a
        requisite post-condition for a successful run of the given MDL file,
        and which must be symlinks, optionally checking the target of the
        symlink against a pre-computed value.

        """
        #if type(fname) == types.StringType:
        if isinstance(fname, str):
            self.extra_checks.append(RequireFileSymlink(fname, target))
        else:
            if target is not None:
                #if type(target) == types.StringType:
                if isinstance(target, str):
                    self.extra_checks.extend([RequireFileSymlink(x, target) for x in fname])
                else:
                    assert len(fname) == len(target)
                    self.extra_checks.extend([RequireFileSymlink(*a) for a in zip(fname, target)])
        return self

    def check_output_files(self):
        """Callback from test_run_context which is responsible for checking the
        state of all expected output files from the executable.

        """
        assertFileExists("realout")
        assertFileExists("realerr")
        for chk in self.extra_checks:
            chk.check()
