/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#ifndef API_GEN_CONSTANTS
#define API_GEN_CONSTANTS

#include <string>

namespace MCell {
namespace API {

const std::string STATE_UNSET = "state_unset";
const int STATE_UNSET_INT = -1;
const int BOND_UNBOUND = 0;
const int BOND_BOUND = -1;
const int ORIENTATION_DOWN = -1;
const int ORIENTATION_NONE = 0;
const int ORIENTATION_UP = 1;
const int ORIENTATION_NOT_SET = 2;

void define_pybinding_constants(py::module& m);

} // namespace API
} // namespace MCell

#endif // API_GEN_CONSTANTS

