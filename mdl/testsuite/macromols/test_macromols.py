#! /usr/bin/env python

from testutils import McellTest
from testutils import get_output_dir
from testutils import cleandir
from testutils import crange
from reaction_output import RequireCountConstraints
from reaction_output import RequireCountEquilibrium
#from viz_output import RequireVizAscii, RequireVizRK, RequireVizDX
#from viz_output import RequireVizDreammV3
#from viz_output import RequireVizDreammV3MolsBin
#from viz_output import RequireVizDreammV3MeshBin
#from viz_output import RequireVizDreammV3MolsAscii
#from viz_output import RequireVizDreammV3MeshAscii
#from viz_output import RequireVizDreammV3Grouped
from macromol_test_types import InvalidMacromolParserTest
import unittest

###################################################################
# Test cases for macromol parser error handling
###################################################################
class TestMacromolParseInvalid(unittest.TestCase):
  pass

def make_invalid_test(i):
  methname = "test%02d" % i
  filename = "invalid-%02d.mdl" % i
  func = lambda self: InvalidMacromolParserTest(filename).invoke(get_output_dir())
  setattr(TestMacromolParseInvalid, methname, func)

## Bulk generate invalid test cases 1...40
for i in crange(1, 40):
  make_invalid_test(i)

###################################################################
# Test cases for numeric validation of macromolecules
###################################################################
class TestMacromolNumeric(unittest.TestCase):

  def test_volume_mixed(self):
    t = McellTest("macromols", "01-macro.mdl", ["-quiet"])
    t.set_check_std_handles(1, 1, 1)
    t.add_extra_check(RequireCountConstraints("dat/01-macro/counts.dat",
                            [(1, 1,  1, 1,  0,  0, 0),          # 180
                             (0, 1, -1, 0,  0,  0, 0),          # 0
                             (1, 1,  0, 0, -1,  0, 0),          # 0
                             (0, 0,  1, 1,  0, -1, 0),          # 0
                             (0, 1,  0, 1,  0,  0, 1)],         # 800
                            [180, 0, 0, 0, 800],
                            header=True))
    t.add_extra_check(RequireCountEquilibrium("dat/01-macro/counts.dat",
                            [62.2, 30.0, 30.0, 57.8, 92.1, 87.9, 712.1],
                            [ 3.0,  1.5,  1.5,  3.0,  4.5,  4.5,  10.0],
                            min_time=0.001,
                            header=True))
    t.invoke(get_output_dir())

  def test_volume_nonmixed(self):
    t = McellTest("macromols", "01-macro_nonmixed.mdl", ["-quiet"])
    t.set_check_std_handles(1, 1, 1)
    t.add_extra_check(RequireCountConstraints("dat/01-macro_nonmixed/counts.dat",
                            [(1, 1,  1, 1,  0,  0, 0),          # 180
                             (0, 1, -1, 0,  0,  0, 0),          # 0
                             (1, 1,  0, 0, -1,  0, 0),          # 0
                             (0, 0,  1, 1,  0, -1, 0),          # 0
                             (0, 1,  0, 1,  0,  0, 1)],         # 800
                            [180, 0, 0, 0, 800],
                            header=True))
    t.add_extra_check(RequireCountEquilibrium("dat/01-macro_nonmixed/counts.dat",
                            [62.2, 30.0, 30.0, 57.8, 92.1, 87.9, 712.1],
                            [ 3.0,  1.5,  1.5,  3.0,  4.5,  4.5,  10.0],
                            min_time=0.01,
                            header=True))
    t.invoke(get_output_dir())

###################################################################
# Generate a test suite for all invalid mdl files
###################################################################
def errorsuite():
  return unittest.makeSuite(TestMacromolParseInvalid, "test")

###################################################################
# Generate a test suite for all numeric tests
###################################################################
def numericsuite():
  return unittest.makeSuite(TestMacromolNumeric, "test")

###################################################################
# Default use of this file will invoke all tests
###################################################################
if __name__ == "__main__":
  cleandir(get_output_dir())
  unittest.main()
