from testutils import assertFileExists, McellTest
from testutils import assertFileEquals
from reaction_output import assertCounts, assertValidTriggerOutput
import time
import re
import os

###################################################################
# Base class for macromol parser tests.  Sets common options.
###################################################################
class MacromolParserTest(McellTest):
  def __init__(self, file, iters=1):
    McellTest.__init__(self, "mmparser", file, ["-quiet", "-iterations", str(iters)])
    self.set_check_std_handles(1, 1, 1)

###################################################################
# Class for parser error tests.
###################################################################
class InvalidMacromolParserTest(MacromolParserTest):
  def __init__(self, file):
    MacromolParserTest.__init__(self, file)
    self.set_expected_exit_code(1)
    self.add_empty_file("realout")
    self.add_nonempty_file("realerr")

