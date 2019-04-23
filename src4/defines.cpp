/******************************************************************************
 *
 * Copyright (C) 2019 by
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

#include "defines.h"

#include <sstream>

using namespace std;

namespace mcell {

std::ostream & operator<<(std::ostream &out, const vec3_t &a) {
  out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return out;
}

string vec3_t::to_string() const {
  stringstream ss;
  ss << *this;
  return ss.str();
}

void vec3_t::dump(const std::string extra_comment, const std::string ind) const {
  cout << ind << extra_comment << *this << "\n";
}

void world_constants_t::dump() {
  cout << "time_unit: \t\t" << time_unit << " [float_t] \t\t\n";
  cout << "length_unit: \t\t" << length_unit << " [float_t] \t\t\n";
  cout << "rx_radius_3d: \t\t" << rx_radius_3d << " [float_t] \t\t\n";
  cout << "partition_edge_length: \t\t" << partition_edge_length << " [float_t] \t\t\n";
  cout << "subpartitions_per_partition_dimension: \t\t" << subpartitions_per_partition_dimension << " [uint32_t] \t\t\n";
  cout << "subpartition_edge_length: \t\t" << subpartition_edge_length << " [float_t] \t\t\n";
}

} // namespace mcell
