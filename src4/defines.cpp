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

#include <iostream>
#include <sstream>

#ifdef DWITHGPERFTOOLS
// using longer path to avoid collisions
#include "install_gperftools/include/profiler.h"
#endif

using namespace std;

namespace MCell {

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


std::ostream & operator<<(std::ostream &out, const vec2_t &a) {
  out << "(" << a.u << ", " << a.v << ")";
  return out;
}


string vec2_t::to_string() const {
  stringstream ss;
  ss << *this;
  return ss.str();
}


void vec2_t::dump(const std::string extra_comment, const std::string ind) const {
  cout << ind << extra_comment << *this << "\n";
}


void SimulationStats::dump() {
  cout << "Total number of ray-subvolume intersection tests (number of ray_trace calls): " << ray_voxel_tests << "\n";
  cout << "Total number of ray-polygon intersection tests: " << ray_polygon_tests << "\n";
  cout << "Total number of ray-polygon intersections: " << ray_polygon_colls << "\n";
  cout << "Total number of molecule moves between walls: " << mol_moves_between_walls << "\n";
}


void SimulationConfig::dump() {
  cout << "time_unit: \t\t" << time_unit << " [float_t] \t\t\n";
  cout << "length_unit: \t\t" << length_unit << " [float_t] \t\t\n";
  cout << "rx_radius_3d: \t\t" << rx_radius_3d << " [float_t] \t\t\n";
  cout << "partition_edge_length: \t\t" << partition_edge_length << " [float_t] \t\t\n";
  cout << "subpartitions_per_partition_dimension: \t\t" << subpartitions_per_partition_dimension << " [uint] \t\t\n";
  cout << "subpartition_edge_length: \t\t" << subpartition_edge_length << " [float_t] \t\t\n";
}

} // namespace mcell
