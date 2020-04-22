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

#include "../api/common.h"

namespace MCell {
namespace API {

void define_pybinding_constants(py::module& m) {
  m.attr("STATE_UNSET") = py::str("state_unset");
  m.attr("STATE_UNSET_INT") = py::int_(-1);
  m.attr("BOND_UNBOUND") = py::int_(0);
  m.attr("BOND_BOUND") = py::int_(-1);
}

} // namespace API
} // namespace MCell

