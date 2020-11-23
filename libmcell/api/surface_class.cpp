/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#include "api/surface_class.h"
#include "api/surface_property.h"

#include <set>

using namespace std;

namespace MCell {
namespace API {

bool SurfaceClass::__eq__(const SurfaceClass& other) const {
  if (!eq_nonarray_attributes(other)) {
    return false;
  }

  // are surface properties the same?
  std::set<SurfaceProperty> s1;
  for (auto& c: properties) {
    s1.insert(*c);
  }
  std::set<SurfaceProperty> s2;
  for (auto& c: other.properties) {
    s2.insert(*c);
  }
  return s1 == s2;
}


} // namespace API
} // namespace MCell
