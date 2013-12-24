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

from test_parser import oldvizsuite, vizsuite, errorsuite, quicksuite, kitchensinksuite, rtcheckpointsuite

tests = {
    "oldvizsuite"       : "VIZ output tests for ASCII/DX modes",
    "vizsuite"          : "VIZ output tests for DREAMM V3 modes",
    "errorsuite"        : "Test error handling for invalid MDL files",
    "quicksuite"        : "A few quick running tests which cover most valid MDL options",
    "kitchensinksuite"  : "Kitchen Sink Test: (very nearly) every parser option",
    "rtcheckpointsuite" : "Basic test of timed checkpoint functionality"
}
collections = {
    "allvizsuite"  : ("VIZ output tests for all modes (old+new)", ["oldvizsuite", "vizsuite"]),
    "fasttests"    : ("All quick running tests (valid+invalid MDL)", ["errorsuite", "quicksuite"]),
}
