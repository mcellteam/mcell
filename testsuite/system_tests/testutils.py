"""Basic MCell testing utilities.

Miscellaneous useful utilities for testing are placed in this module.  In
particular, this module contains classes and functions related to invoking and
testing generic MCell runs.  Most in-depth output-testing utilities have been
placed in other modules (such as viz_output or reaction_output).

Most new tests will use MCellTest (or a subclass of MCellTest) to handle
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
import types
import random
import re

###################################################################
# Check if assertions are enabled
###################################################################
def check_assertions_enabled():
  assertions_enabled = 0
  
  try:
    assert 1 == 2
  except AssertionError:
    assertions_enabled = 1

  return assertions_enabled

###################################################################
# Complain and exit if assertions are not enabled
###################################################################
def require_assertions_enabled():
  if check_assertions_enabled() == 0:
    print "Please turn off Python optimization (remove the -O flag from the command line)."
    print "The optimization feature turns off assertions, which are vital to the functioning"
    print "of the test suite."
    raise Exception("Assertions are disabled, testsuite cannot run")

require_assertions_enabled()

###################################################################
# Get the filename of the source file containing the function which calls this
# function.
###################################################################
def get_caller_filename(lvl=1):
  frame = sys._getframe()
  for i in range(0, lvl):
    frame = frame.f_back
  return frame.f_code.co_filename

###################################################################
# Utility to safely concatenate two iterables, even if one or the other may be
# None, or the two iterables have different types.  If they have different
# types, the return type will be a tuple.
###################################################################
def safe_concat(l1, l2):
  """safe_concat(l1, l2) -> combined list.

  Concatenate two sequences safely, even if one or both may be None, or the
  types may differ.
  """
  if l1 is None:
    return l2
  elif l2 is None:
    return l1
  elif type(l1) is not type(l2):
    return tuple(l1) + tuple(l2)
  else:
    return l1 + l2

###################################################################
# Get the preferred output directory.
###################################################################
def get_output_dir():
  g = globals()
  if not g.has_key("test_output_dir"):
    global test_output_dir
    test_output_dir = "./test_results"
  return g["test_output_dir"]

###################################################################
# Utility to generate a closed range
###################################################################
def crange(l, h, s=1):
  """crange(start, stop[, step]) -> list of integers.

  Works like 'range', but is closed on both ends -- i.e. includes both
  endpoints.
  """
  if s > 0:
    return range(l, h+1, s)
  else:
    return range(l, h-1, s)

###################################################################
# Make sure a directory exists, and contains no garbage.  Used to ensure a
# clean working directory for each run of the testsuite.
###################################################################
def cleandir(directory):
  """Ensure that the speficied directory exists, and that it is empty.
  """
  try:
    os.mkdir(directory)
  except:
    pass
  for root, dirs, files in os.walk(directory, topdown=False):
    map(lambda f: os.unlink(os.path.join(root, f)), files)
    map(lambda d: os.rmdir(os.path.join(root, d)), dirs)

###################################################################
# Give an error if a file doesn't exist
###################################################################
def assertFileExists(fname):
  """Raise an error if the specified file does not exist.
  """
  try:
    os.stat(fname)
  except:
    assert False, "Expected file '%s' was not created" % fname

class RequireFileExists:
  def __init__(self, name):
    self.name = name

  def check(self):
    assertFileExists(self.name)

###################################################################
# Give an error if a file does exist
###################################################################
def assertFileNotExists(fname):
  """Raise an error if the specified file exists.
  """
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

###################################################################
# Give an error if a file doesn't exist or isn't empty
###################################################################
def assertFileEmpty(fname):
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

###################################################################
# Give an error if a file doesn't exist or is empty
###################################################################
def assertFileNonempty(fname, expect=None):
  try:
    sb = os.stat(fname)
  except:
    assert False, "Expected file '%s' was not created" % fname
  assert sb.st_size != 0, "Expected file '%s' shouldn't be empty but is" % fname
  if expect != None:
    assert sb.st_size == expect, "Expected file '%s' is the wrong size (%d bytes instead of %d bytes)" % (fname, sb.st_size, expect)

class RequireFileNonempty:
  def __init__(self, name, expect=None):
    self.name = name
    self.expect = expect

  def check(self):
    assertFileNonempty(self.name, self.expect)

###################################################################
# Give an error if a file doesn't exist or if its contents differ from the
# contents passed in.
###################################################################
def assertFileEquals(fname, contents):
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

###################################################################
# Give an error if a file doesn't exist or if its contents match (or do not
# match) the regular expression passed in.
###################################################################
def assertFileMatches(fname, regex, expectMinMatches=1, expectMaxMatches=sys.maxint):
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
  assert len(matches) <= expectMaxMatches, "Expected file '%s' had incorrect contents (found %d match%s for regex '%s', expected at most %d)"  % (fname, len(matches), plural, regex, expectMaxMatches)

class RequireFileMatches:
  def __init__(self, name, regex, expectMinMatches=1, expectMaxMatches=sys.maxint):
    self.name = name
    self.regex = regex
    self.expectMinMatches = expectMinMatches
    self.expectMaxMatches = expectMaxMatches

  def check(self):
    assertFileMatches(self.name, self.regex, self.expectMinMatches, self.expectMaxMatches)

###################################################################
# Give an error if a file isn't a symlink, or (optionally) if the destination
# of the symlink isn't as specified.
###################################################################
def assertFileSymlink(fname, target=None):
  try:
    sb = os.lstat(fname)
  except:
    assert False, "Expected symlink '%s' was not created" % fname
  assert stat.S_ISLNK(sb.st_mode), "Expected symlink '%s' is not a symlink" % fname
  if target != None:
    got_target = os.readlink(fname)
    assert got_target == target, "Expected symlink '%s' should point to '%s', but instead points to '%s'" % (fname, target, got_target)

class RequireFileSymlink:
  def __init__(self, name, target=None):
    self.name = name
    self.target = target

  def check(self):
    assertFileSymlink(self.name, self.target)

###################################################################
# Give an error if a file isn't a directory.
###################################################################
def assertFileDir(fname):
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

###################################################################
# Small wrapper class to handle loading test suite configuration
###################################################################
class test_config(object):
  def __init__(self, filepath="./test.cfg"):
    dict = {
        "mcellpath": "mcell"
    }
    self.config = ConfigParser.ConfigParser(dict)
    try:
      self.config.read(filepath)
    except:
      pass

  def get(self, sect, val):
    if self.config.has_section(sect):
      return self.config.get(sect, val)
    else:
      return self.config.get(ConfigParser.DEFAULTSECT, val)

###################################################################
# Class to handle running command-line executables as PyUnit tests.
###################################################################
class test_run_context(object):
  testidx = 1

  def __init__(self, cmd, args):
    self.command = cmd
    self.args    = args
    self.testidx = test_run_context.testidx
    self.check_stdout = 0
    self.check_stderr = 0
    self.check_stdin  = 0
    self.expect_exitcode = 0
    self.cleaned = False
    test_run_context.testidx += 1

  def set_check_std_handles(self, i, o, e):
    self.check_stdin  = i
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

  def __check_results(self):
    assert not os.WIFSIGNALED(self.got_exitcode), "Process died due to signal %d" % os.WTERMSIG(self.got_exitcode)
    assert os.WEXITSTATUS(self.got_exitcode) == self.expect_exitcode, "Expected exit code %d, got exit code %d" % (self.expect_exitcode, os.WEXITSTATUS(self.got_exitcode))
    if self.check_stdout:
      self.check_stdout_valid(os.path.join(os.getcwd(), "stdout"))
    if self.check_stderr:
      self.check_stderr_valid(os.path.join(os.getcwd(), "stderr"))
    self.check_output_files()

  def __run(self):
    pid = os.fork()
    if pid != 0:
      pid, self.got_exitcode = os.waitpid(pid, 0)
    else:
      newout = os.dup(1)

      # Close stdin
      if self.check_stdin:
        os.close(0)

      # Save stdout
      new_stdout = os.open("./stdout", os.O_CREAT | os.O_WRONLY | os.O_EXCL, 0644)
      os.dup2(new_stdout, 1)

      # Save stderr
      new_stderr = os.open("./stderr", os.O_CREAT | os.O_WRONLY | os.O_EXCL, 0644)
      os.dup2(new_stderr, 2)

      self.__runchild(newout)

  def __runchild(self, newout):
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

    try:
      os.execvp(self.command, self.args)
    except:
      os.dup2(newout, 1)
      sys.__excepthook__(* sys.exc_info())
      sys.exit(127)

###################################################################
# Specialized class for running MCell jobs as PyUnit tests
###################################################################
class McellTest(test_run_context):
  """Utility base class for MCell tests.
  
  This class will build up the command-line, choosing a random seed,
  looking up the MCell executable in a configuration file, redirect output
  to log/err files using -logfile and -errfile, and do a handful of
  global validations (criteria which must be true for any MCell run).
  """

  rand = random.Random()
  config = test_config("./test.cfg")

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
    """Add an extra post-run check to this test.  The check must be
    an object that has a 'check' method on it, which contains appropriate
    'assert' statements to perform the check.

    Typically, this is used as follows:
        o = MCellTest(...)
        o.add_extra_check(RequireFileEquals(filename, contents))
        o.invoke()

    This is especially useful for types of checks that cannot be added via the
    utility methods given below.
    """
    self.extra_checks.append(chk)

  ## Add one or more expected output files (either a string or an iterable)
  def add_exist_file(self, e):
    """Utility to add a file or list of files whose existence is a requisite
    post-condition for a successful run of the given MDL file.
    """
    if type(e) == types.StringType:
      self.extra_checks.append(RequireFileExists(e))
    else:
      self.extra_checks.extend([RequireFileExists(x) for x in e])
    return self

  ## Add one or more expected empty output files (either a string or an
  ## iterable)
  def add_empty_file(self, e):
    """Utility to add a file or list of files whose existence (and emptiness)
    is a requisite post-condition for a successful run of the given MDL file.
    """
    if type(e) == types.StringType:
      self.extra_checks.append(RequireFileEmpty(e))
    else:
      self.extra_checks.extend([RequireFileEmpty(x) for x in e])
    return self

  ## Add one or more expected non-empty output files (either a string or an
  ## iterable)
  def add_nonempty_file(self, e, expected_size=None):
    """Utility to add a file or list of files whose existence (and
    non-emptiness) is a requisite post-condition for a successful run of the
    given MDL file.
    """
    if type(e) == types.StringType:
      self.extra_checks.append(RequireFileNonempty(e, expected_size))
    else:
      self.extra_checks.extend([RequireFileNonempty(x, expected_size) for x in e])
    return self

  ## Add one or more expected constant output files (either fname and cnt must
  ## be strings, or they must be iterables with the same number of items)
  def add_constant_file(self, fname, cnt):
    """Utility to add a file or list of files whose existence is a requisite
    post-condition for a successful run of the given MDL file, and whose
    contents must exactly match a specified string.
    """
    if type(fname) == types.StringType:
      self.extra_checks.append(RequireFileEquals(fname, cnt))
    else:
      self.extra_checks.extend([RequireFileEquals(x, cnt) for x in fname])
    return self

  ## Add one or more expected symlinks, with optional target.  fname may be a
  ## string or an iterable returning a string.  target may be a string, or if
  ## fname is an iterable, an iterable returning a string.
  def add_symlink(self, fname, target=None):
    """Utility to add a file or list of filenames whose existence is a
    requisite post-condition for a successful run of the given MDL file, and
    which must be symlinks, optionally checking the target of the symlink
    against a pre-computed value.
    """
    if type(fname) == types.StringType:
      self.extra_checks.append(RequireFileSymlink(fname, target))
    else:
      if target != None:
        if type(target) == types.StringType:
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
