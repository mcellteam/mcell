#! /usr/bin/env python

from testutils import McellTest
from testutils import get_output_dir
from testutils import cleandir
from testutils import crange
from testutils import RequireFileMatches
from reaction_output import RequireCountConstraints
from reaction_output import RequireCountEquilibrium

import unittest

class TestRegressions(unittest.TestCase):

  def test_001(self):
    mt = McellTest("regression", "01-remove_per_species_list_from_ht.mdl", ["-quiet"])
    mt.set_check_std_handles(1, 1, 1)
    mt.invoke(get_output_dir())

  def test_002(self):
    mt = McellTest("regression", "02-orientflipflip_rxn.mdl", ["-quiet"])
    mt.set_check_std_handles(1, 1, 1)
    mt.add_extra_check(RequireCountConstraints("counts.txt",
                                               [(1, 0,  0,  0, -1,  0),    # 0
                                                (0, 1,  0, -1,  0,  0),    # 0
                                                (1, 1, -1,  0,  0,  0),    # 0
                                                (0, 0,  0,  1,  1, -1),    # 0
                                                (0, 0,  1,  0,  0,  0),    # 1000
                                                (0, 0,  0,  0,  0,  1)],   # 1000
                                               [0, 0, 0, 0, 1000, 1000],
                                               header=True))
    mt.invoke(get_output_dir())

  def test_003(self):
    mt = McellTest("regression", "03-coincident_surfaces.mdl", ["-quiet"])
    mt.set_check_std_handles(1, 1, 1)
    mt.add_extra_check(RequireCountConstraints("cannonballs.txt",
                                               [(1, 1, -1,  0,  0,  0),    # 0
                                                (0, 0,  0,  1,  1, -1),    # 0
                                                (0, 0,  1,  0,  0,  0),    # 500
                                                (0, 0,  0,  0,  0,  1)],   # 500
                                               [0, 0, 500, 500],
                                               header=True))
    mt.invoke(get_output_dir())

  def test_004(self):
    mt = McellTest("regression", "04-rx_flipflip.mdl", ["-quiet"])
    mt.set_check_std_handles(1, 1, 1)
    mt.add_extra_check(RequireCountConstraints("counts.txt",
                                               [(1, 1,  0,  0,  0),    # 300
                                                (0, 0,  1,  1,  0)],   # 300
                                               [300, 300],
                                               header=True))
    mt.add_extra_check(RequireCountConstraints("counts.txt",
                                               [],
                                               [],
                                               minimums=[ 40,  40,  40,  40, 10000],
                                               maximums=[260, 260, 260, 260, 1e300],
                                               min_time=5e-5,
                                               header=True))
    mt.invoke(get_output_dir())

  def test_005(self):
    mt = McellTest("regression", "05-rx_dissociate_inwards.mdl", ["-quiet"])
    mt.set_check_std_handles(1, 1, 1)
    mt.add_extra_check(RequireCountConstraints("molecules.txt",
                                               [(1, 1,  0),         # 1000
                                                (0, 1, -1)],        # 0
                                               [1000, 0],
                                               header=True))
    mt.invoke(get_output_dir())

  def test_006(self):
    mt = McellTest("regression", "06-misreporting_rxn_products.mdl", ["-quiet"])
    mt.set_check_std_handles(1, 1, 1)
    mt.add_extra_check(RequireFileMatches("realout", '\s*Probability.*set for a\{0\} \+ b\{0\} -> c\{0\}', expectMaxMatches=1))
    mt.add_extra_check(RequireFileMatches("realout", '\s*Probability.*set for a\{0\} \+ b\{0\} -> d\{0\}', expectMaxMatches=1))
    mt.add_extra_check(RequireFileMatches("realout", '\s*Probability.*set for a\{0\} \+ b\{0\} -> e\{0\}', expectMaxMatches=1))
    mt.invoke(get_output_dir())

def suite():
  return unittest.makeSuite(TestRegressions, "test")
