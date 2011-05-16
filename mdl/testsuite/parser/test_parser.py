#! /usr/bin/env python

from testutils import get_output_dir
from testutils import cleandir
from testutils import crange
from testutils import RequireFileEquals
from viz_output import RequireVizAscii, RequireVizRK, RequireVizDX
from viz_output import RequireVizDreammV3
from viz_output import RequireVizDreammV3MolsBin
from viz_output import RequireVizDreammV3MeshBin
from viz_output import RequireVizDreammV3MolsAscii
from viz_output import RequireVizDreammV3MeshAscii
from viz_output import RequireVizDreammV3Grouped
from parser_test_types import InvalidParserTest
from parser_test_types import ValidParserTest
from parser_test_types import CheckpointTest
from parser_test_types import KitchenSinkParserTest
from reaction_output import RequireCounts
import unittest

###################################################################
# Test cases for valid "kitchen sink" parses, no viz output
###################################################################
class TestParseValid(unittest.TestCase):

  def test_vanilla1(self):
    KitchenSinkParserTest("01-kitchen_sink").invoke(get_output_dir())

  def test_vanilla2(self):
    KitchenSinkParserTest("01-kitchen_sink_viz_NONE").invoke(get_output_dir())

  def test_fullyrandom(self):
    KitchenSinkParserTest("01-kitchen_sink_fully_random").invoke(get_output_dir())

  def test_microrev_surface(self):
    KitchenSinkParserTest("01-kitchen_sink_mr_surf").invoke(get_output_dir())

  def test_microrev_all(self):
    KitchenSinkParserTest("01-kitchen_sink_mr_true").invoke(get_output_dir())

  def test_microrev_volume(self):
    KitchenSinkParserTest("01-kitchen_sink_mr_vol").invoke(get_output_dir())

  def test_space_step(self):
    KitchenSinkParserTest("01-kitchen_sink_space_step").invoke(get_output_dir())

  def test_silent(self):
    KitchenSinkParserTest("01-kitchen_sink_silent", silent=True).invoke(get_output_dir())
    
  def test_utility_cmds(self):
    t = ValidParserTest("01-kitchen_sink_utility_cmds_grammar.mdl")
    t.add_nonempty_file("my_file.dat")
    t.add_nonempty_file("exp.dat")
    t.add_extra_check(RequireCounts("exp.dat", [(f*1e-6,10) for f in range (0,2)]))
    t.invoke(get_output_dir())

###################################################################
# Test cases for valid "kitchen sink" parses, ASCII/RK viz modes
###################################################################
class TestParseVizAscii(unittest.TestCase):

  def test_viz_ascii1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_ascii_iterframedata")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.invoke(get_output_dir())

  def test_viz_ascii2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_ascii_iter")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.invoke(get_output_dir())

  def test_viz_ascii3(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_ascii_time")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.invoke(get_output_dir())

  def test_viz_rk1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_rk1")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.invoke(get_output_dir())

  def test_viz_rk2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_rk2")
    t.add_extra_check(RequireVizRK("viz_dat/%s/molecules.rk.dat" % t.basename, 10))
    t.invoke(get_output_dir())

###################################################################
# Test cases for valid "kitchen sink" parses, old DX viz modes
###################################################################
class TestParseVizDx(unittest.TestCase):

  def test_viz_dx1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dx_iterframedata")
    t.add_extra_check(RequireVizDX("viz_dat/%s" % t.basename,
                                   molfile="molecules",
                                   objprefixes=["objects_a", "objects_b"], 
                                   mpositers=[1] + crange(3, 100, 10),
                                   mstateiters=[11, 25],
                                   epositers=crange(1, 100, 10),
                                   estateiters=[1] + crange(2, 100, 10),
                                   opositers=[1, 2],
                                   ostateiters=[1, 3]))
    t.invoke(get_output_dir())

  def test_viz_dx2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dx_iter")
    t.add_extra_check(RequireVizDX("viz_dat/%s" % t.basename,
                                   molfile="molecules",
                                   objprefixes=["objects_a", "objects_b"], 
                                   alliters=crange(1, 100, 10)))
    t.invoke(get_output_dir())

  def test_viz_dx3(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dx_time")
    t.add_extra_check(RequireVizDX("viz_dat/%s" % t.basename,
                                   molfile="molecules",
                                   objprefixes=["objects_a", "objects_b"], 
                                   alliters=crange(1, 100, 10)))
    t.invoke(get_output_dir())

###################################################################
# Test cases for valid "kitchen sink" parses, DREAMM viz modes
###################################################################
class TestParseVizDreamm(unittest.TestCase):
  ALL_MOLECULES = [
      "%s_%d" % (m, i)
        for (m, lim) in [("c_s", 0),
                         ("c_v", 1),
                         ("m_g", 11),
                         ("m_v", 11),
                         ("s_g", 11),
                         ("s_v", 11)]
        for i in crange(0, lim)
  ]

  ALL_WORLDS = [ "world", "world2" ]
  ALL_MESHES = [ "box1", "box2", "box3", "newbox", "poly1", "poly2" ]

  ALL_OBJECTS = [
      "%s.big_object.%s" % (w, o)
        for w in ALL_WORLDS
        for o in ALL_MESHES
  ]

  ALL_OBJECTS_WITH_REGIONS = [
      "%s.big_object.%s" % (w, o)
        for w in ALL_WORLDS
        for o in ["box1", "box2", "box3", "poly2"]
  ]

  def test_viz_default1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_default_iter")
    outdir = "viz_dat/%s/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outdir, "world", 101, 101))
    t.add_extra_check(RequireVizDreammV3MolsBin(outdir, crange(0, 100)))
    t.add_extra_check(RequireVizDreammV3MeshBin(outdir, crange(0, 100),
                                                positers=[1, 2, 3],
                                                regioniters=[1, 3, 4]))
    t.invoke(get_output_dir())

  def test_viz_default2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_default_time")
    outdir = "viz_dat/%s/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outdir, "world", 101, 101))
    t.add_extra_check(RequireVizDreammV3MolsAscii(outdir, crange(0, 100), TestParseVizDreamm.ALL_MOLECULES))
    t.add_extra_check(RequireVizDreammV3MeshBin(outdir, crange(0, 100),
                                                positers=[1, 2, 3],
                                                regioniters=[1, 3, 4]))
    t.invoke(get_output_dir())

  def test_viz_dreamm1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dreamm_time")
    outdir = "viz_dat/%s/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outdir, "world", 101, 101))
    t.add_extra_check(RequireVizDreammV3MolsBin(outdir, crange(0, 100)))
    t.add_extra_check(RequireVizDreammV3MeshAscii(outdir, crange(0, 100), TestParseVizDreamm.ALL_OBJECTS,
                                                  objswithregions=TestParseVizDreamm.ALL_OBJECTS_WITH_REGIONS,
                                                  positers=[1, 3, 4],
                                                  regioniters=[1, 2, 4]))
    t.invoke(get_output_dir())

  def test_viz_dreammgrouped1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dreammgrouped_time")
    outdir = "viz_dat/%s" % t.basename
    t.add_extra_check(RequireVizDreammV3Grouped(outdir, "world", n_iters=101))
    t.invoke(get_output_dir())

  def test_viz_dreammgrouped2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dreammgrouped_wildcard")
    outdir = "viz_dat/%s" % t.basename
    t.add_extra_check(RequireVizDreammV3Grouped(outdir, "world", n_iters=101, molstate=True, molsnonempty=False))
    t.add_extra_check(RequireVizDreammV3Grouped(outdir, "world2", n_iters=101))
    t.invoke(get_output_dir())

  def test_viz_everything(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_everything")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/01 - ASCII old/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.add_extra_check(RequireVizAscii("viz_dat/%s/02 - ASCII old/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.add_extra_check(RequireVizAscii("viz_dat/%s/03 - ASCII old/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.add_extra_check(RequireVizAscii("viz_dat/%s/04 - ASCII new/world" % t.basename, [0, 2, 3] + crange(10, 100, 10), [2147483647], [2147483647]))
    t.add_extra_check(RequireVizAscii("viz_dat/%s/05 - ASCII new/world" % t.basename, [0, 2, 3] + crange(10, 100, 10), [2147483647], [2147483647]))
    t.add_extra_check(RequireVizAscii("viz_dat/%s/11 - RK old ASCII/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.add_extra_check(RequireVizRK("viz_dat/%s/12 - RK old/molecules.rk.dat" % t.basename, 10))
    t.add_extra_check(RequireVizAscii("viz_dat/%s/13 - RK new ASCII/world" % t.basename, [0, 2, 3] + crange(10, 100, 10), [2147483647], [2147483647]))
    t.add_extra_check(RequireVizRK("viz_dat/%s/14 - RK new/world.rk.dat" % t.basename, 11))
    t.add_extra_check(RequireVizDX("viz_dat/%s/21 - DX old" % t.basename,
                                   molfile="molecules",
                                   objprefixes=["objects_a", "objects_b", "objects_c"], 
                                   mpositers=[1] + crange(3, 100, 10),
                                   mstateiters=[11, 25],
                                   epositers=crange(1, 100, 10),
                                   estateiters=[1] + crange(2, 100, 10),
                                   opositers=[1, 2],
                                   ostateiters=[1, 3]))
    t.add_extra_check(RequireVizDX("viz_dat/%s/22 - DX old" % t.basename,
                                   molfile="molecules",
                                   objprefixes=["objects_a", "objects_b"], 
                                   alliters=crange(1, 100, 10)))
    t.add_extra_check(RequireVizDX("viz_dat/%s/23 - DX old" % t.basename,
                                   molfile="molecules",
                                   objprefixes=["objects_a", "objects_b"], 
                                   alliters=crange(1, 100, 10)))
    t.add_extra_check(RequireVizDX("viz_dat/%s/24 - DX new" % t.basename,
                                   molfile="world",
                                   objprefixes=["objects_a", "objects_b"],
                                   mpositers=[0, 2, 3] + crange(10, 100, 10),
                                   mstateiters=crange(0, 100, 10),
                                   epositers=[0, 2, 3] + crange(10, 100, 10),
                                   estateiters=crange(0, 100, 10),
                                   opositers=[1, 2, 3],
                                   ostateiters=[1]))
    t.add_extra_check(RequireVizDX("viz_dat/%s/25 - DX new" % t.basename,
                                   molfile="world",
                                   objprefixes=["objects_a", "objects_b"],
                                   mpositers=[0, 2, 3] + crange(10, 100, 10),
                                   mstateiters=crange(0, 100, 10),
                                   epositers=[0, 2, 3] + crange(10, 100, 10),
                                   estateiters=crange(0, 100, 10),
                                   opositers=[1, 2, 3],
                                   ostateiters=[1]))
    outsubdir = "viz_dat/%s/34 - DREAMM new/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outsubdir, "world", 14, 15))
    t.add_extra_check(RequireVizDreammV3MolsBin(outsubdir, [0, 2, 3, 4] + crange(10, 100, 10),
                                                surfpositers=[0, 2, 3] + crange(10, 100, 10),
                                                volpositers=[0, 2, 3] + crange(10, 100, 10),
                                                surforientiters=[2, 3, 4],
                                                volorientiters=[2, 3, 4]))
    t.add_extra_check(RequireVizDreammV3MeshBin(outsubdir, crange(1, 4),
                                                positers=[1, 2, 3],
                                                regioniters=[1, 3, 4]))
    outsubdir = "viz_dat/%s/35 - DREAMM new/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outsubdir, "world", 14, 15))
    t.add_extra_check(RequireVizDreammV3MolsAscii(outsubdir, [0, 2, 3, 4] + crange(10, 100, 10),
                      ["c_s_0", "c_v_0", "c_v_1"] + 
                      ["%c_%c_%d" % (k, j, i) for i in crange(0, 11) for j in ['v', 'g'] for k in ['m', 's']],
                      positers=[0, 2, 3] + crange(10, 100, 10),
                      orientiters=[2, 3, 4]))
    t.add_extra_check(RequireVizDreammV3MeshBin(outsubdir, crange(1, 4),
                                                positers=[1, 3, 4],
                                                regioniters=[1, 2, 4]))
    outsubdir = "viz_dat/%s/36 - DREAMM new/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outsubdir, "world", 14, 15))
    t.add_extra_check(RequireVizDreammV3MolsBin(outsubdir, [0, 2, 3, 4] + crange(10, 100, 10),
                                                surfpositers=[0, 2, 3] + crange(10, 100, 10),
                                                volpositers=[0, 2, 3] + crange(10, 100, 10),
                                                surforientiters=[2, 3, 4],
                                                volorientiters=[2, 3, 4]))
    t.add_extra_check(RequireVizDreammV3MeshAscii(outsubdir, crange(1, 4),
                        ["%s.big_object.%s" % (j, i) for i in ["box1", "box2", "box3", "newbox", "poly1", "poly2"] for j in ["world", "world2"]],
                        objswithregions=["%s.big_object.%s" % (j, i) for i in ["box1", "box2", "box3", "poly2"] for j in ["world", "world2"]],
                        positers=[1, 2, 3],
                        regioniters=[1, 3, 4]))
    outsubdir = "viz_dat/%s/37 - DREAMM new/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outsubdir, "world", 14, 15))
    t.add_extra_check(RequireVizDreammV3MolsBin(outsubdir, [0, 2, 3, 4] + crange(10, 100, 10),
                                                surfpositers=[0, 2, 3] + crange(10, 100, 10),
                                                volpositers=[0, 2, 3] + crange(10, 100, 10),
                                                surforientiters=[2, 3, 4],
                                                volorientiters=[2, 3, 4]))
    t.add_extra_check(RequireVizDreammV3MeshBin(outsubdir, crange(1, 4),
                                                positers=[1, 2, 3],
                                                regioniters=[1, 3, 4]))
    t.add_extra_check(RequireVizDreammV3Grouped("viz_dat/%s/44 - DREAMM grouped new" % t.basename, "world", n_iters=14, n_times=15))
    t.add_extra_check(RequireVizDreammV3Grouped("viz_dat/%s/45 - DREAMM grouped new" % t.basename, "world2", n_iters=14, n_times=15))
    t.add_extra_check(RequireVizDreammV3Grouped("viz_dat/%s/46 - DREAMM grouped new" % t.basename, "world", n_iters=14, n_times=15, molstate=True, molsnonempty=False))
    t.invoke(get_output_dir())

###################################################################
# Test cases for RT checkpoint parsing
###################################################################
class TestParseRtCheckpoint(unittest.TestCase):

  def test01(self):
    CheckpointTest("02-rtcheckpoint_seconds").invoke(get_output_dir())

  def test02(self):
    CheckpointTest("02-rtcheckpoint_minutes").invoke(get_output_dir())

  def test03(self):
    CheckpointTest("02-rtcheckpoint_minutes2").invoke(get_output_dir())

  def test04(self):
    CheckpointTest("02-rtcheckpoint_hours").invoke(get_output_dir())

  def test05(self):
    CheckpointTest("02-rtcheckpoint_hours2").invoke(get_output_dir())

  def test06(self):
    CheckpointTest("02-rtcheckpoint_days").invoke(get_output_dir())

  def test07(self):
    CheckpointTest("02-rtcheckpoint_days2").invoke(get_output_dir())

  def test08(self):
    CheckpointTest("02-rtcheckpoint_noexitspec").invoke(get_output_dir())

###################################################################
# Test cases for parser error handling
###################################################################
class TestParseInvalid(unittest.TestCase):

  def test025(self):
    InvalidParserTest("invalid-025.mdl").add_empty_file("invalid-025.tmp").invoke(get_output_dir())

def make_invalid_test(i):
  methname = "test%03d" % i
  filename = "invalid-%03d.mdl" % i
  func = lambda self: InvalidParserTest(filename).invoke(get_output_dir())
  setattr(TestParseInvalid, methname, func)

## Bulk generate invalid test cases 1...23, 26...27, 29...85
## 25 is specially specified, and 24 and 28 do not presently exist.
for i in crange(1, 23) + crange(26, 85):
  make_invalid_test(i)

###################################################################
# Generate a test suite for non-DREAMM viz output modes
###################################################################
def oldvizsuite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestParseVizAscii, "test"))
  suite.addTest(unittest.makeSuite(TestParseVizDx,    "test"))
  return suite

###################################################################
# Generate a test suite for DREAMM viz output modes
###################################################################
def vizsuite():
  return unittest.makeSuite(TestParseVizDreamm, "test")

###################################################################
# Generate a test suite for all invalid mdl files
###################################################################
def errorsuite():
  return unittest.makeSuite(TestParseInvalid, "test")

###################################################################
# Generate a test suite for all valid "kitchen sink" (non-viz) files
###################################################################
def kitchensinksuite():
  return unittest.makeSuite(TestParseValid, "test")

###################################################################
# Generate a test suite for all realtime checkpoint files
###################################################################
def rtcheckpointsuite():
  return unittest.makeSuite(TestParseRtCheckpoint, "test")

###################################################################
# Generate a test suite for all fast-running tests except for the error tests)
###################################################################
def quicksuite():
  suite = unittest.TestSuite()
  suite.addTest(TestParseVizDreamm("test_viz_dreamm1"))
  suite.addTest(TestParseRtCheckpoint("test01"))
  return suite

###################################################################
# Generate a test suite with all tests included.
###################################################################
def fullsuite():
  suite = unittest.TestSuite()
  suite.addTest(allvizsuite())
  suite.addTest(errorsuite())
  suite.addTest(kitchensinksuite())
  suite.addTest(rtcheckpointsuite())
  return suite

###################################################################
# Default use of this file will invoke all tests
###################################################################
if __name__ == "__main__":
  cleandir(get_output_dir())
  unittest.main()
