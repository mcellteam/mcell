#! /usr/bin/env python

from testutils import cleandir
from testutils import crange
from viz_output import RequireVizAscii, RequireVizRK, RequireVizDX
from viz_output import RequireVizDreammV3
from viz_output import RequireVizDreammV3MolsBin
from viz_output import RequireVizDreammV3MeshBin
from viz_output import RequireVizDreammV3MolsAscii
from viz_output import RequireVizDreammV3MeshAscii
from viz_output import RequireVizDreammV3Grouped
from parser_test_types import InvalidParserTest
from parser_test_types import CheckpointTest
from parser_test_types import KitchenSinkParserTest
import unittest

###################################################################
# Test cases for valid "kitchen sink" parses, no viz output
###################################################################
class TestParseValid(unittest.TestCase):

  def test_vanilla1(self):
    KitchenSinkParserTest("01-kitchen_sink").invoke("./test_results")

  def test_vanilla2(self):
    KitchenSinkParserTest("01-kitchen_sink_viz_NONE").invoke("./test_results")

  def test_fullyrandom(self):
    KitchenSinkParserTest("01-kitchen_sink_fully_random").invoke("./test_results")

  def test_microrev_surface(self):
    KitchenSinkParserTest("01-kitchen_sink_mr_surf").invoke("./test_results")

  def test_microrev_all(self):
    KitchenSinkParserTest("01-kitchen_sink_mr_true").invoke("./test_results")

  def test_microrev_volume(self):
    KitchenSinkParserTest("01-kitchen_sink_mr_vol").invoke("./test_results")

  def test_space_step(self):
    KitchenSinkParserTest("01-kitchen_sink_space_step").invoke("./test_results")

  def test_silent(self):
    KitchenSinkParserTest("01-kitchen_sink_silent", silent=True).invoke("./test_results")

###################################################################
# Test cases for valid "kitchen sink" parses, ASCII/RK viz modes
###################################################################
class TestParseVizAscii(unittest.TestCase):

  def test_viz_ascii1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_ascii_iterframedata")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.invoke("./test_results")

  def test_viz_ascii2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_ascii_iter")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.invoke("./test_results")

  def test_viz_ascii3(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_ascii_time")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.invoke("./test_results")

  def test_viz_rk1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_rk1")
    t.add_extra_check(RequireVizAscii("viz_dat/%s/molecules" % t.basename, crange(1, 100, 10), [4], [3]))
    t.invoke("./test_results")

  def test_viz_rk2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_rk2")
    t.add_extra_check(RequireVizRK("viz_dat/%s/molecules.rk.dat" % t.basename, 10))
    t.invoke("./test_results")

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
    t.invoke("./test_results")

  def test_viz_dx2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dx_iter")
    t.add_extra_check(RequireVizDX("viz_dat/%s" % t.basename,
                                   molfile="molecules",
                                   objprefixes=["objects_a", "objects_b"], 
                                   alliters=crange(1, 100, 10)))
    t.invoke("./test_results")

  def test_viz_dx3(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dx_time")
    t.add_extra_check(RequireVizDX("viz_dat/%s" % t.basename,
                                   molfile="molecules",
                                   objprefixes=["objects_a", "objects_b"], 
                                   alliters=crange(1, 100, 10)))
    t.invoke("./test_results")

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
    t.add_extra_check(RequireVizDreammV3(outdir, "world", 101))
    t.add_extra_check(RequireVizDreammV3MolsBin(outdir, crange(0, 100)))
    t.add_extra_check(RequireVizDreammV3MeshBin(outdir, crange(0, 100),
                                                positers=[1, 2, 3],
                                                regioniters=[1, 3, 4]))
    t.invoke("./test_results")

  def test_viz_default2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_default_time")
    outdir = "viz_dat/%s/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outdir, "world", 101))
    t.add_extra_check(RequireVizDreammV3MolsAscii(outdir, crange(0, 100), TestParseVizDreamm.ALL_MOLECULES))
    t.add_extra_check(RequireVizDreammV3MeshBin(outdir, crange(0, 100),
                                                positers=[1, 2, 3],
                                                regioniters=[1, 3, 4]))
    t.invoke("./test_results")

  def test_viz_dreamm1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dreamm_time")
    outdir = "viz_dat/%s/world_viz_data" % t.basename
    t.add_extra_check(RequireVizDreammV3(outdir, "world", 101))
    t.add_extra_check(RequireVizDreammV3MolsBin(outdir, crange(0, 100)))
    t.add_extra_check(RequireVizDreammV3MeshAscii(outdir, crange(0, 100), TestParseVizDreamm.ALL_OBJECTS,
                                                  objswithregions=TestParseVizDreamm.ALL_OBJECTS_WITH_REGIONS,
                                                  positers=[1, 3, 4],
                                                  regioniters=[1, 2, 4]))
    t.invoke("./test_results")

  def test_viz_dreammgrouped1(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dreammgrouped_time")
    outdir = "viz_dat/%s" % t.basename
    t.add_extra_check(RequireVizDreammV3Grouped(outdir, "world", n_iters=101))
    t.invoke("./test_results")

  def test_viz_dreammgrouped2(self):
    t = KitchenSinkParserTest("01-kitchen_sink_viz_dreammgrouped_wildcard")
    outdir = "viz_dat/%s" % t.basename
    t.add_extra_check(RequireVizDreammV3Grouped(outdir, "world", n_iters=101, molstate=True))
    t.invoke("./test_results")

###################################################################
# Test cases for RT checkpoint parsing
###################################################################
class TestParseRtCheckpoint(unittest.TestCase):

  def test01(self):
    CheckpointTest("02-rtcheckpoint_seconds").invoke("./test_results")

  def test02(self):
    CheckpointTest("02-rtcheckpoint_minutes").invoke("./test_results")

  def test03(self):
    CheckpointTest("02-rtcheckpoint_minutes2").invoke("./test_results")

  def test04(self):
    CheckpointTest("02-rtcheckpoint_hours").invoke("./test_results")

  def test05(self):
    CheckpointTest("02-rtcheckpoint_hours2").invoke("./test_results")

  def test06(self):
    CheckpointTest("02-rtcheckpoint_days").invoke("./test_results")

  def test07(self):
    CheckpointTest("02-rtcheckpoint_days2").invoke("./test_results")

  def test08(self):
    CheckpointTest("02-rtcheckpoint_noexitspec").invoke("./test_results")

###################################################################
# Test cases for parser error handling
###################################################################
class TestParseInvalid(unittest.TestCase):

  def test025(self):
    InvalidParserTest("invalid-025.mdl").add_empty_file("invalid-025.tmp").invoke("./test_results")

def make_invalid_test(i):
  methname = "test%03d" % i
  filename = "invalid-%03d.mdl" % i
  func = lambda self: InvalidParserTest(filename).invoke("./test_results")
  setattr(TestParseInvalid, methname, func)

## Bulk generate invalid test cases 1...23, 26...27, 29...85
## 25 is specially specified, and 24 and 28 do not presently exist.
for i in crange(1, 23) + crange(26, 27) + crange(29, 85):
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
# Generate a test suite for all viz output modes
###################################################################
def allvizsuite():
  suite = unittest.TestSuite()
  suite.addTest(oldvizsuite())
  suite.addTest(vizsuite())
  return suite

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
# Generate a test suite for all fast-running tests (currently, just the parser
# error tests)
###################################################################
def quicksuite():
  suite = errorsuite()
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
  cleandir("./test_results")
  unittest.main()
