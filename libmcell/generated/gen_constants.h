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

enum class Orientation {
  Down = -1,
  None = 0,
  Up = 1,
  NotSet = 2
};


static inline  std::ostream& operator << (std::ostream& out, const Orientation v) {
  switch (v) {
    case Orientation::Down: out << "Orientation.Down (-1)"; break;
    case Orientation::None: out << "Orientation.None (0)"; break;
    case Orientation::Up: out << "Orientation.Up (1)"; break;
    case Orientation::NotSet: out << "Orientation.NotSet (2)"; break;
  }
  return out;
};

enum class Notification {
  None = 0,
  Brief = 1,
  Full = 2
};


static inline  std::ostream& operator << (std::ostream& out, const Notification v) {
  switch (v) {
    case Notification::None: out << "Notification.None (0)"; break;
    case Notification::Brief: out << "Notification.Brief (1)"; break;
    case Notification::Full: out << "Notification.Full (2)"; break;
  }
  return out;
};

void define_pybinding_constants(py::module& m);

} // namespace API
} // namespace MCell

#endif // API_GEN_CONSTANTS

