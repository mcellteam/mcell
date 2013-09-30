#! /usr/bin/env python

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

from testutils import McellTest
from testutils import get_output_dir
from testutils import cleandir
from testutils import crange
from testutils import RequireFileMatches
from reaction_output import RequireCountConstraints
from reaction_output import RequireCountEquilibrium
from reaction_output import RequireCounts
from reaction_output import RequireCountsPositive
from reaction_output import RequireCountsLessThan

import os
import unittest

class TestRegressions(unittest.TestCase):

  def test_001(self):
    mt = McellTest("regression", "01-remove_per_species_list_from_ht.mdl", ["-quiet"])
    mt.invoke(get_output_dir())

  def test_002(self):
    mt = McellTest("regression", "02-orientflipflip_rxn.mdl", ["-quiet"])
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
    mt.add_extra_check(RequireCountConstraints("molecules.txt",
                                               [(1, 1,  0),         # 1000
                                                (0, 1, -1)],        # 0
                                               [1000, 0],
                                               header=True))
    mt.invoke(get_output_dir())

  def test_006(self):
    mt = McellTest("regression", "06-misreporting_rxn_products.mdl", ["-quiet"])
    mt.add_extra_check(RequireFileMatches("realout", '\s*Probability.*set for a\{0\} \+ b\{0\} -> c\{0\}', expectMaxMatches=1))
    mt.add_extra_check(RequireFileMatches("realout", '\s*Probability.*set for a\{0\} \+ b\{0\} -> d\{0\}', expectMaxMatches=1))
    mt.add_extra_check(RequireFileMatches("realout", '\s*Probability.*set for a\{0\} \+ b\{0\} -> e\{0\}', expectMaxMatches=1))
    mt.invoke(get_output_dir())

  def test_007(self):
    mt = McellTest("regression", "07-volvol_crash.mdl", ["-quiet"])
    mt.invoke(get_output_dir())

  def test_008(self):
    mt = McellTest("regression", "08-find_corresponding_region.mdl", ["-quiet"])
    mt.invoke(get_output_dir())

  def __rename(self, p1, p2):
    os.rename(p1, p2)

  def test_009(self):
    mt = McellTest("regression", "09-incorrect_times_in_chkpt_A.mdl", ["-quiet"])
    testpath = '%s/test-%04d' % (get_output_dir(), mt.testidx)
    mt.invoke(get_output_dir())
    self.__rename(os.path.join(testpath, 'cmdline.txt'), os.path.join(testpath, 'cmdline_0.txt'))
    self.__rename(os.path.join(testpath, 'realout'), os.path.join(testpath, 'realout.0'))
    self.__rename(os.path.join(testpath, 'realerr'), os.path.join(testpath, 'realerr.0'))
    self.__rename(os.path.join(testpath, 'stdout'),  os.path.join(testpath, 'stdout.0'))
    self.__rename(os.path.join(testpath, 'stderr'),  os.path.join(testpath, 'stderr.0'))

    mt.args[-1] = os.path.join(os.path.dirname(mt.args[-1]), "09-incorrect_times_in_chkpt_B.mdl")
    mt.invoke(get_output_dir())
    self.__rename(os.path.join(testpath, 'cmdline.txt'), os.path.join(testpath, 'cmdline_1.txt'))
    self.__rename(os.path.join(testpath, 'realout'), os.path.join(testpath, 'realout.1'))
    self.__rename(os.path.join(testpath, 'realerr'), os.path.join(testpath, 'realerr.1'))
    self.__rename(os.path.join(testpath, 'stdout'),  os.path.join(testpath, 'stdout.1'))
    self.__rename(os.path.join(testpath, 'stderr'),  os.path.join(testpath, 'stderr.1'))

    times = [i*1e-7 for i in range(0, 6)] + [(2*i+1)*1e-7 for i in range(3, 13)]
    values = zip(times, [0]*len(times))
    mt.add_extra_check(RequireCounts('A_World.dat', values))
    mt.invoke(get_output_dir())
    self.__rename(os.path.join(testpath, 'cmdline.txt'), os.path.join(testpath, 'cmdline_2.txt'))
    self.__rename(os.path.join(testpath, 'realout'), os.path.join(testpath, 'realout.2'))
    self.__rename(os.path.join(testpath, 'realerr'), os.path.join(testpath, 'realerr.2'))
    self.__rename(os.path.join(testpath, 'stdout'),  os.path.join(testpath, 'stdout.2'))
    self.__rename(os.path.join(testpath, 'stderr'),  os.path.join(testpath, 'stderr.2'))

  def test_010(self):
    mt = McellTest("regression", "10-counting_crashes_on_coincident_wall.mdl", ["-quiet"])
    mt.invoke(get_output_dir())

  def test_011(self):
    mt = McellTest("regression", "11-quoted_tickmark_counts_parse_error.mdl", ["-quiet"])
    mt.invoke(get_output_dir())

  def test_012(self):
    mt = McellTest("regression", "12-no_waypoints_counting_fail.mdl", ["-quiet"])
    mt.add_extra_check(RequireCountConstraints("counts.txt",
                                                # A+ ----------------------   B+ ----------------------   C+ ----------------------   A- ----------   B- ----------   C- ----------
                                               [( 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 1,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 500
                                                ( 0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0),  # 0
                                                ( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1)], # 0
                                               [0]*6 + [500] + [0]*25))
    mt.add_extra_check(RequireCountConstraints("counts.txt",
                                               constraints=[],
                                               totals=[],
                                               min_time=9e-5,
                                               minimums=[1]*21 + [0]*12))
    mt.invoke(get_output_dir())

  def test_013(self):
    mt = McellTest("regression", "13-find_corresponding_region_error_crash_1.mdl", ["-quiet"])
    mt.add_extra_check(RequireFileMatches('realerr', r"Can't find new region corresponding to foo\.boxB,ALL for world\.bar\.rs \(copy of foo\.bar\.rs\)"))
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())

  def test_014(self):
    mt = McellTest("regression", "14-find_corresponding_region_error_crash_2.mdl", ["-quiet"])
    mt.add_extra_check(RequireFileMatches('realerr', r"Can't find new region corresponding to foo\.boxB,ALL for world\.bar\.rs \(copy of foo\.bar\.rs\)"))
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())

  def test_015(self):
    mt = McellTest("regression", "15-mol_grid_grid_crash.mdl", ["-quiet"])
    mt.invoke(get_output_dir())

  def test_016(self):
    mt = McellTest("regression", "16-mol_surf_crash.mdl", ["-quiet"])
    mt.invoke(get_output_dir())

  def test_017(self):
    mt = McellTest("regression", "17-uninstantiated_reference_crash.mdl", ["-quiet"])
    mt.add_extra_check(RequireFileMatches('realerr', r"Region neither instanced nor grouped with release site."))
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())

  def test_018(self):
    mt = McellTest("regression", "18-uninstantiated_reference_crash_2.mdl", ["-quiet"])
    mt.add_extra_check(RequireFileMatches('realerr', r"Cannot produce visualization for the uninstantiated object 'uninstantiated'"))
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())

  def test_019(self):
    mt = McellTest("regression", "19-enclosed_surfmol_miscount.mdl", ["-quiet"])
    mt.add_extra_check(RequireCountConstraints("A.dat",
                                               [(1, -1, -1)],
                                               [0],
                                               minimums=[50, 0, 0],
                                               maximums=[50, 50, 50]))
    mt.add_extra_check(RequireCountEquilibrium("A.dat",
                                               values=[50.0, 25.0, 25.0],
                                               tolerances=[0.0, 5.0, 5.0]))
    mt.invoke(get_output_dir())
  
  def test_020(self):
    mt = McellTest("regression", "20-reaction_null_products_crash.mdl", ["-quiet"])
    mt.invoke(get_output_dir())

  def test_021(self):
    mt = McellTest("regression", "21-enclosed_meshes_with_different_properties.mdl", ["-quiet"])
    mt.add_extra_check(RequireCounts("A.dat",[(f*1e-6, 0) for f in range (0,101)]))
    mt.invoke(get_output_dir())

  def test_022(self):
    mt = McellTest("regression", "22-rx_reflective_surface_bug.mdl", ["-quiet"])
    mt.add_extra_check(RequireCountsPositive("dat/22-rx_reflective_surface_bug/refl.dat", "# Iteration_# sm_L sm_M"))
    mt.invoke(get_output_dir())

  def test_023(self):
    mt = McellTest("regression", "23-walls_overlap_two_objects.mdl", ["-quiet"])
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())

  def test_024(self):
    mt = McellTest("regression", "24-walls_overlap_partial.mdl", ["-quiet"])
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())

  def test_025(self):
    mt = McellTest("regression", "25-walls_overlap_complete.mdl", ["-quiet"])
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())
  
  def test_026(self):
    mt = McellTest("regression", "26-walls_overlap_share_edge.mdl", ["-quiet"])
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())

  def test_027(self):
    mt = McellTest("regression", "27-walls_coincident.mdl", ["-quiet"])
    mt.set_expected_exit_code(1)
    mt.invoke(get_output_dir())

  def test_028(self):
    mt = McellTest("regression", "28-unimolecular_reaction_after_chkpt.mdl", ["-quiet"])
    testpath = '%s/test-%04d' % (get_output_dir(), mt.testidx)
    mt.invoke(get_output_dir())
    self.__rename(os.path.join(testpath, 'cmdline.txt'), os.path.join(testpath, 'cmdline_0.txt'))
    self.__rename(os.path.join(testpath, 'realout'), os.path.join(testpath, 'realout.0'))
    self.__rename(os.path.join(testpath, 'realerr'), os.path.join(testpath, 'realerr.0'))
    self.__rename(os.path.join(testpath, 'stdout'),  os.path.join(testpath, 'stdout.0'))
    self.__rename(os.path.join(testpath, 'stderr'),  os.path.join(testpath, 'stderr.0'))

    mt.invoke(get_output_dir())
    self.__rename(os.path.join(testpath, 'cmdline.txt'), os.path.join(testpath, 'cmdline_1.txt'))
    self.__rename(os.path.join(testpath, 'realout'), os.path.join(testpath, 'realout.1'))
    self.__rename(os.path.join(testpath, 'realerr'), os.path.join(testpath, 'realerr.1'))
    self.__rename(os.path.join(testpath, 'stdout'),  os.path.join(testpath, 'stdout.1'))
    self.__rename(os.path.join(testpath, 'stderr'),  os.path.join(testpath, 'stderr.1'))
    mt.add_extra_check(RequireCountsLessThan("B_World.dat", 10))
   
  def test_029(self):
    mt = McellTest("regression", "29-unimolecular_FOREVER.mdl", ["-quiet"])
    mt.add_extra_check(RequireCountEquilibrium("S1_down.dat", [10] *1, 
                               [3] *1, 1e-3, 1.5e-3, header=True))
    mt.invoke(get_output_dir())

 
def suite():
  return unittest.makeSuite(TestRegressions, "test")
