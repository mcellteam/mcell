#! /usr/bin/env python

from testutils import McellTest
from testutils import get_output_dir
from testutils import cleandir
from reaction_output import RequireCountConstraints
from reaction_output import RequireCountEquilibrium
from reaction_output import RequireCountRxnRate
import unittest
import math

class TestReactionsNumeric(unittest.TestCase):

  def test_volume(self):
    t = McellTest("reactions", "01-volume_highconc.mdl", ["-quiet"])
    t.set_check_std_handles(1, 1, 1)
    t.add_extra_check(RequireCountConstraints("dat/01-volume_highconc/V_out.dat",
                            [(0, 1, -1, 0,  0,  0, 0, 0, 0),          # 0
                             (0, 0,  0, 1, -1,  0, 0, 0, 0),          # 0
                             (0, 0,  0, 1,  0, -1, 0, 0, 0),          # 0
                             (1, 0,  0, 0,  0,  0, 1, 0, 0),          # 1000
                             (0, 1,  0, 0,  0,  0, 0, 1, 0),          # 1000
                             (0, 0,  0, 1,  0,  0, 0, 0, 1)],         # 1000
                            [0, 0, 0, 1000, 1000, 1000],
                            header=True))
    t.add_extra_check(RequireCountEquilibrium("dat/01-volume_highconc/V_out.dat",
                            [500, 500, 500, 500, 500, 500, 500, 500, 500],
                            [ 25,  25,  25,  25,  25,  25,  25,  25,  25],
                            header=True))
    t.add_extra_check(RequireCountRxnRate("dat/01-volume_highconc/rxn_out.dat",
                            values=    [1e5,   1e5,   1e5,   1e5,   1e5,   1e5],
                            tolerances=[1.5e4, 1.5e4, 1.5e4, 1.5e4, 1.5e4, 1.5e4],
                            min_time=5e-3,
                            base_time=0.0,
                            header=True))
    t.invoke(get_output_dir())

###################################################################
# Generate a test suite for all numeric tests
###################################################################
def numericsuite():
  return unittest.makeSuite(TestReactionsNumeric, "test")

###################################################################
# Default use of this file will invoke all tests
###################################################################
if __name__ == "__main__":
  cleandir(get_output_dir())
  unittest.main()
