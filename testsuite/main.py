#! /usr/bin/env python

# If Psyco is installed, fire it up so that the test suite will run more
# quickly
try:
  import psyco
  psyco.full()
except ImportError:
  pass

import sys
import os
import unittest
from optparse import OptionParser

# Add the system_tests directory to our python path
mypath = os.path.dirname(sys.argv[0])
sys.path.append(os.path.join(mypath, "system_tests"))
import testutils

def all(iterable):
  """
  Return True if all elements of the iterable are true 
  (or if the iterable is empty).
  
  iterable - any object with an iterable context.
  """

  for i in iterable:
    if not i:
      return False
  return True

def exist(d, iterable):
  """
  Return generator object consisting of Boolean values.

  d - dictionary
  iterable - any object with an iterable context.
             At present a list or single item.
  """
  for i in iterable:
    yield d.has_key(i)

def uniq(iterable):
  """
  Return generator object consisting
  from the uniq items from the 'iterable'.
  'iterable' should be presorted.

  iterable - any object with iterable context
  """
  last = None
  for i in iterable:
    if i != last:
      last = i
      yield i

class TestSet(object):
  """
  Describe an hierarchy of aggregate test names/test collections

  """ 
  def __init__(self, key, desc):
    self.key = key
    self.desc = desc
    self.children = None
    self.suite = None

def clean_pyc_files(dir):
  """
  Remove files ended on 'pyc' or 'pyo' from directory 'dir'

  dir - directory name
  """
  def path_visitor(a, d, f):
    for i in f:
      fpath = os.path.join(d, i)
      if fpath.lower().endswith(".pyc")  or  fpath.lower().endswith(".pyo"):
        os.unlink(fpath)
    del f[:]
  os.path.walk(dir, path_visitor, None)

def find_all_tests(path, name='', create_all=False):
  """
  Return a dictionary of TestSet objects.
  Searches recursively the supplied directory tree
  for the names/descriptions of the tests on the directory level.

  path - path to the testsuite location
  name - a string in extracting components of the file path name
  create_all - flag
  """
  locals = dict()
  sys.path.append(path)
  clean_pyc_files(path)
  execfile(os.path.join(path, "test_info.py"), globals(), locals)
  basename = os.path.basename(path)
  tests = {}
  individual = []
  if locals.has_key("tests"):
    for k, v in locals["tests"].items():
      if locals.has_key(k):
        t = TestSet(k, v)
        t.suite = locals[k]()
        tests[k] = t
        individual.append(t)
      else:
        print 'Ignoring test "%s", as no corresponding test object was found.' % k
  if locals.has_key("collections"):
    for k, v in locals["collections"].items():
      subtests = v[1]
      desc = v[0]
      if all(exist(tests, subtests)):
        t = TestSet(k, desc)
        t.children = {}
        for i in subtests:
          t.children[i] = None
        tests[k] = t
      else:
        print 'Ignoring test "%s", as one or more subtests were missing.' % k
  if locals.has_key("subdirs"):
    for k, v in locals["subdirs"].items():
      subname = os.path.join(name, k)
      subtests = find_all_tests(os.path.join(path, k), subname, True)
      if subtests.has_key('all'):
        t = subtests['all']
        del subtests['all']
        t.key = subname
        t.desc = v
        t.children = subtests
        tests[subname] = t
        individual.append(t)
  if create_all  and  len(individual) != 0:
    t = TestSet("all", "All tests")
    tests["all"] = t
  return tests

def print_tests(t, indent='  '):
  """
  Print out the dictionary of tests names/descriptions
  on the directory level (aggregated).
  Used when running with option '-list or -l'

  t - dictionary of test name: TestSet object pairs
  indent - string used for formatted output
  """
  for n, v in t.items():
    if v is None:
      print indent + "- (" + n + ")"
    else:
      print indent + "- " + n + " : " + v.desc
      if v.children is not None:
        print_tests(v.children, indent + '  ')

def include_test(all_tests, rt, inc):
  """
  Added test names specified in 'inc'
  to the list of tests names collection 'rt'
  
  all_tests - dictionary of test name: TestSet object pairs
  rt - list of aggregate test names
  inc - name of the aggregate test directory or test collection
  """
  path = inc.split('/')
  if not all_tests.has_key(path[0]):
    return
  rt.append(inc)

def append_expanded(t, rt, prefix):
  """
  Recursively adds names of the test collections held in the TestSet.children
  to the list of test collections names 'rt'

  it - TestSet object
  rt - list of tests names (they may be directory names)
  prefix - part of the path name used for extracting test
           collection names on the different level of granularity
  """
  new_prefix = os.path.join(prefix, t.key)
  if t.suite is not None:
    rt.append(new_prefix)

  if t.children is not None:
    for i, child in t.children.items():
      if child is None:
        rt.append(os.path.join(prefix, i))
      else:
        append_expanded(child, rt, new_prefix)

def expand_test(all_tests, rt, exp, prefix=''):
  """
  Recursively expands the list 'rt' to include
  test collections on the levels below the names/paths
  as specified in 'rt'

  all_tests - dictionary of test name: TestSet object pairs
  rt - list of aggregate test names
  exp - name of the aggregate test
  prefix - string used to extract parts of the file path name
  """
  full = os.path.join(prefix, exp)
  try:
    while True:
      rt.remove(full)
  except:
    pass
  s = exp.split('/')
  if len(s) > 1 and all_tests.has_key(s[0]):
    t = all_tests[s[0]].children
    expand_test(t, rt, str.join('/', s[1:]), os.path.join(prefix, s[0]))
  elif all_tests.has_key(exp):
    append_expanded(all_tests[exp], rt, prefix)

def exclude_test(all_tests, rt, exc):
  """
  Remove the aggregate test name specified in 'exc'
  from the list of aggregate test names 'rt'
  
  all_tests - dictionary of test name: TestSet object pairs
  rt - list of aggregate test names
  exc - list of aggregate test names to remove from 'rt'
  """
  path = exc.split('/')
  if not all_tests.has_key(path[0]):
    return

  names = [exc]
  expand_test(all_tests, names, exc)
  for i in names:
    try:
      while True:
        rt.remove(i)
    except:
      pass

def generate_run_tests(all_tests, rt, inc, exc):
  """
  Create list of aggregate test names/test collections

  all_tests - dictionary of test name: TestSet object pairs
  rt - list of aggregate test names
  inc - list of aggregate test names to include into 'rt'
  exc - list of aggregate test names to remove from 'rt'
  """
  if inc is not None:
    for i in inc:
      include_test(all_tests, rt, i)
  ort = list(rt)
  for r in ort:
    expand_test(all_tests, rt, r)
  if exc is not None:
    for e in exc:
      exclude_test(all_tests, rt, e)

def add_to_test_suite(suite, all_tests, r):
  """
  Recursively add to the test suite aggregate
  tests described in the path name 'r'

  suite - object of unittest.TestSuite()
  all_tests - dictionary of test name: TestSet object pairs
  r - name of the aggregate test to add to the suite
  """
  path = r.split('/')
  while len(path) > 1:
    all_tests = all_tests[path[0]].children
    del path[0]
  suite.addTest(all_tests[path[0]].suite)

def build_test_suite(all_tests, run_tests):
  """
  Create a TestSuite() suite

  all_tests - dictionary of test name: TestSet object pairs
  run_tests - list of aggregate test names
  """
  suite = unittest.TestSuite()
  for r in run_tests:
    add_to_test_suite(suite, all_tests, r)
  return suite

op = OptionParser()
op.add_option("-c", "--config",   dest="config",    default="./test.cfg",        help="load configuration from CONFIG")
op.add_option("-T", "--testpath", dest="testpath",  default="../mdl/testsuite",  help="run tests in MDL files located under TESTPATH")
op.add_option("-l", "--list",     dest="list",      action="store_true",         help="display a list of all tests found under the test directory")
op.add_option("-i", "--include",  dest="include",   action="append",             help="comma-separated list of tests to include (default: all tests)")
op.add_option("-e", "--exclude",  dest="exclude",   action="append",             help="comma-separated list of tests to exclude (default: none)")
op.add_option("-v", "--verbose",  dest="verbosity", action="count",              help="increase the verbosity of the test suite")
op.add_option("-r", "--results",  dest="results",   default="./test_results",    help="run all MCell tests under the directory RESULTS (WARNING: directory will be wiped clean first!)")
(options, args) = op.parse_args()

all_tests = find_all_tests(options.testpath)

# List all tests, if requested
if options.list:
  print "Found tests:"
  print_tests(all_tests)
  sys.exit(0)

# Load our configuration.
test_conf = testutils.test_config(options.config)

# Check our configuration for which tests to run
try:
  run_tests = [x.strip() for x in test_conf.get("main", "run_tests").split(',')]
except:
  if options.include is None:
    run_tests = all_tests.keys()
  else:
    run_tests = []

# Parse include options
include = None
if options.include and len(options.include) != 0:
  include = []
  [include.extend(x) for x in [y.split(',') for y in options.include]]

# Parse exclude options
exclude = None
if options.exclude and len(options.exclude) != 0:
  exclude = []
  [exclude.extend(x) for x in [y.split(',') for y in options.exclude]]

# Generate the final list of tests to run
generate_run_tests(all_tests, run_tests, include, exclude)
run_tests.sort()
run_tests = list(uniq(run_tests))

# Check that run_tests is not empty...
if len(run_tests) == 0:
  print "No tests to run."
  sys.exit(0)

# Report on which tests are being run
print "Running tests:"
for i in run_tests:
  print "  - " + i

# Build a top-level suite containing all relevant tests
suite = build_test_suite(all_tests, run_tests)
testutils.cleandir(options.results)
unittest.TextTestRunner(verbosity=options.verbosity).run(suite)
