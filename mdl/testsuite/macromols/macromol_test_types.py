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

