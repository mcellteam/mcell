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


void uint_set_t::dump() {
  cout << "Indices contained in a subpartition: ";
  int cnt = 0;
  for (uint idx: *this) {
    cout << idx << ", ";

    if (cnt %20 == 0 && cnt != 0) {
      cout << "\n";
    }
    cnt++;
  }
  cout << "\n";
}


void simulation_stats_t::dump() {
  cout << "Total number of ray-subvolume intersection tests (number of ray_trace calls): " << ray_voxel_tests << "\n";
  cout << "Total number of ray-polygon intersection tests: " << ray_polygon_tests << "\n";
  cout << "Total number of ray-polygon intersections: " << ray_polygon_colls << "\n";
}

} // namespace mcell
