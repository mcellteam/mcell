from testutils import assertFileExists, McellTest
from testutils import assertFileEquals
from reaction_output import assertCounts, assertValidTriggerOutput
import time
import re
import os

###################################################################
# Base class for parser tests.  Sets common options.
###################################################################
class ParserTest(McellTest):
  def __init__(self, file, iters=1):
    McellTest.__init__(self, "parser", file, ["-quiet", "-iterations", str(iters)])
    self.set_check_std_handles(1, 1, 1)

###################################################################
# Class for parser error tests.
###################################################################
class InvalidParserTest(ParserTest):
  def __init__(self, file):
    ParserTest.__init__(self, file)
    self.set_expected_exit_code(1)
    self.add_empty_file("realout")
    self.add_nonempty_file("realerr")

###################################################################
# Class for parser checkpoint tests.
###################################################################
class CheckpointTest(ParserTest):
  def __init__(self, file):
    McellTest.__init__(self, "parser", file + ".mdl", ["-quiet", "-iterations", "2000000000"])
    self.set_check_std_handles(1, 1, 1)
    self.add_empty_file("realout")
    self.add_empty_file("realerr")
    self.basename = file
    self.add_nonempty_file(file + ".stamp")
    self.add_nonempty_file(file + ".cp")

  def check_output_files(self):
    ParserTest.check_output_files(self)
    stampfile = self.basename + ".stamp"
    cpfile    = self.basename + ".cp"

    stamptime = os.stat(stampfile).st_mtime
    cptime    = os.stat(cpfile).st_mtime

    assert (cptime - stamptime) < 35, "Realtime checkpoint scheduled for 30 seconds, but time between timestamp and checkpoint exceeds 35 seconds (%d seconds)" % (cptime - stamptime)
    assert (cptime - stamptime) > 25, "Realtime checkpoint scheduled for 30 seconds, but time between timestamp and checkpoint is less than seconds (%d seconds)" % (cptime - stamptime)

###################################################################
# Class for parser "kitchen sink" tests.
###################################################################
class KitchenSinkParserTest(ParserTest):

  DAT_FILE_TEMPLATE = """43110, world
Message should have been "43110, world"
Current day of the week is %s
"""
  def __init__(self, file, silent=False):
    ParserTest.__init__(self, file + ".mdl", 100)
    self.basename = file
    self.preday = time.strftime('%A', time.localtime())
    self.add_nonempty_file(file + "-expressions.dat")
    self.add_nonempty_file(file + ".dat")
    self.add_nonempty_file(["dat/%s/counting-%d.txt" % (file, i) for i in range(1, 18)])
    self.add_nonempty_file(["vol_dat/%s/vo.%d.dat" % (file, i) for i in [0, 100]])
    self.add_nonempty_file(["vol_dat/%s/ovo.%d.dat" % (file, i) for i in [1, 10, 100]])
    self.add_nonempty_file("vol_dat/%s/nvo.100.dat" % file)
    self.add_nonempty_file("%s.cp" % file)
    if not silent:
      self.add_nonempty_file("realout")
    else:
      self.add_empty_file("realout")
    if not silent:
      self.add_constant_file("realerr", "Warning: could not release 6 of s_g_1 (surface full)\n")
    else:
      self.add_empty_file("realerr")

  def check_output_files(self):
    ParserTest.check_output_files(self)

    # Check contents of <basename>-expressions.dat
    self.__check_expressions()

    # Check contents of <basename>.dat
    self.__check_formatted_io()

    # Check contents of counting-*.txt?
    self.__check_reaction_output()

    # Check contents of vo.*.dat/ovo.*.dat/nvo.*.dat
    for vo in ["vol_dat/%s/vo.%d.dat" % (self.basename, i) for i in [0, 100]]:
      self.__check_volume_output(vo, 25, 25, 25)
    for vo in ["vol_dat/%s/ovo.%d.dat" % (self.basename, i) for i in [1, 10, 100]]:
      self.__check_volume_output(vo, 1, 2, 4)
    self.__check_volume_output("vol_dat/%s/nvo.100.dat" % self.basename, 1, 1, 1)

  def __check_formatted_io(self):
    postday = time.strftime('%A', time.localtime())
    if self.preday == postday:
      assertFileEquals(self.basename + ".dat", KitchenSinkParserTest.DAT_FILE_TEMPLATE % postday)
    else:
      try:
        assertFileEquals(self.basename + ".dat", KitchenSinkParserTest.DAT_FILE_TEMPLATE % self.preday)
      except AssertionError:
        assertFileEquals(self.basename + ".dat", KitchenSinkParserTest.DAT_FILE_TEMPLATE % postday)

  def __check_expressions(self):
    lines = open(self.basename + "-expressions.dat").read().split('\n')
    exact    = [l[5:] for l in lines if l.find(" == ") != -1]
    gaussian = [l[5:] for l in lines if l.find(" ~= ") != -1]
    for e in exact:
      assert eval(e), "Expression evaluated incorrectly: %s" % e
    gauss_exp = re.compile('''([0-9e.-]+) *~= *([0-9e.-]+) */ *([0-9e.-]+)''')
    for g in gaussian:
      m = gauss_exp.search(g)
      assert m != None, 'Expression formatted incorrectly: %s' % g
      val  = float(m.group(1))
      mean = float(m.group(2))
      std  = float(m.group(3))
      if val < (mean - 2.0*std)  or  val > (mean + 2.0*std):
        print 'Warning: Gaussian value outside of 95% confidence interval.'
        print '         Expected to occur occasionally.'
        print '         Value: %.15g  (mean=%.15g, stdev=%.15g)' % (val, mean, std)

  ###
  ### Check reaction output.  This should probably do a lot more than it does.
  ###
  def __check_reaction_output(self):
    # 1...  Constant count of 5, times spaced every 1e-6 from 0 to 1e-4
    file = "dat/%s/counting-1.txt" % self.basename
    assertCounts(file, [(f*1e-6, 5) for f in range(0, 101)])

    # 2...  Constant count of 6, times spaced every 1e-6 from 0 to 1e-4
    # Exact time turned off -- should have no effect on counts
    file = "dat/%s/counting-2.txt" % self.basename
    assertCounts(file, [(f*1e-6, 6) for f in range(0, 101)])

    # 3...  Constant count of 7, times spaced every 1e-6 from 0 to 1e-4
    # Header explicitly turned off -- should have
    file = "dat/%s/counting-3.txt" % self.basename
    assertCounts(file, [(f*1e-6, 7) for f in range(0, 101)])

    # 4...  Constant count of 8, times spaced every 1e-6 from 0 to 1e-4
    # Exact time turned on -- should have no effect on counts
    file = "dat/%s/counting-4.txt" % self.basename
    assertCounts(file, [(f*1e-6, 8) for f in range(0, 101)])

    # 5...  Constant count of 8, times spaced every 1e-6 from 0 to 1e-4
    # Header turned off (YES/NO syntax, rather than NONE)
    file = "dat/%s/counting-5.txt" % self.basename
    assertCounts(file, [(f*1e-6, 9) for f in range(0, 101)])

    # 6...  Constant count of 10, times spaced every 1e-6 from 0 to 1e-4
    # Header turned on, but no custom header given
    file = "dat/%s/counting-6.txt" % self.basename
    assertCounts(file, [(f*1e-6, 10) for f in range(0, 101)], "Seconds constant #6")

    # 7...  Constant counts of 11, 12, 13, times spaced every 1e-6 from 0 to 1e-4
    # Custom header "Duck, Duck, Goose" given
    file = "dat/%s/counting-7.txt" % self.basename
    assertCounts(file, [(f*1e-6, 11, 12, 13) for f in range(0, 101)], "Duck, Duck, GooseSeconds constant #6 constant #7 constant #8")

    # Testsuite needs some work before the next several will be validatable...
    # 8...  Counts of all molecules in the world, times every 1e-6 from 0 to 1e-4
    # 9...  Counts of s_g_1 and s_g_2 molecules in the world as two columns,
    #       every 2nd iteration.  s_g_1 and s_g_2 specified as wildcards, but
    #       with no special characters (i.e. quoted)
    # 10... Total of all molecules in the world, every 2nd iteration
    # 11... Multicolumn, exercising count expression functionality, every 2nd
    #       iteration
    # 12... Count of s_g_1' on a region, every 2nd iteration
    # 13... Count of named reaction pathway rxn1, every 2nd iteration
    # 14... Count of s_g_1{0} on a region, every 2nd iteration
    # 15... Count of s_g_1{2} on a region, every 2nd iteration
    # 16... Count of s_g_0 subunits on a c_s_0 complex in the world, every 2nd
    #       iteration
    # 17... Count of s_g_0 subunits on a c_s_0 complex, including two related
    #       subunit restrictions, in the world, every 2nd iteration
    # 18... Count of c_s_0 in the world, time specified as a TIME_LIST, every
    #       7e-6 seconds from 1e-6 to 20e-6.

    # 19... Trigger for count of s_v_4 in world.big_object.newbox.  Check that
    #       it exists.  Could try to validate against counting-8.txt...  Exact
    #       time is left at default (ON).
    file = "dat/%s/counting-19.txt" % self.basename
    assertFileExists(file)
    assertValidTriggerOutput(file, 2, True, xrange=(-20, 20), yrange=(-20, 20), zrange=(-20, 20))

    # 20... Trigger for count of s_v_4 in world.big_object.newbox.  Check that
    #       it exists.  Could try to validate against counting-8.txt...  Exact
    #       time is explicitly turned OFF.
    file = "dat/%s/counting-20.txt" % self.basename
    assertFileExists(file)
    assertValidTriggerOutput(file, 2, False, xrange=(-20, 20), yrange=(-20, 20), zrange=(-20, 20))

    # 21... Trigger for count of s_v_0 in world.big_object.newbox.  Check that
    #       it exists.  Could try to validate against counting-8.txt...  Exact
    #       time is explicitly turned ON.
    file = "dat/%s/counting-21.txt" % self.basename
    assertFileExists(file)
    assertValidTriggerOutput(file, 2, True, xrange=(-20, 20), yrange=(-20, 20), zrange=(-20, 20))

  def __check_volume_output(self, vo, nx, ny, nz):
    lines = open(vo).read().split('\n')
    header_exp = re.compile('''# *nx=([0-9]+) *ny=([0-9]+) *nz=([0-9]+)''')
    m = header_exp.search(lines[0])
    assert m != None, "Volume output file '%s' has no header line" % vo
    assert nx == int(m.group(1)), "Volume output file '%s' has incorrect x dimensions (%s instead of %d)" % (vo, m.group(1), nx)
    assert ny == int(m.group(2)), "Volume output file '%s' has incorrect y dimensions (%s instead of %d)" % (vo, m.group(2), ny)
    assert nz == int(m.group(3)), "Volume output file '%s' has incorrect z dimensions (%s instead of %d)" % (vo, m.group(3), nz)
    expect_numlines = (ny+1)*nz + 1 + 1     # Extra +1 because split creates '' after final '\n'
    assert len(lines) == expect_numlines, "Volume output file '%s' has incorrect number of lines (%d instead of %d)" % (vo, len(lines), expect_numlines)

