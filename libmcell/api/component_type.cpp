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

#include "api/component_type.h"

#include <set>

using namespace std;

namespace MCell {
namespace API {

bool ComponentType::__eq__(const ComponentType& other) const {
  return
      eq_nonarray_attributes(other) &&
      std::set<string>(states.begin(), states.end()) ==
      std::set<string>(other.states.begin(), other.states.end());
}


std::string ComponentType::to_bngl_str() const {
  std::string res;

  res = name;
  for (const string& s: states) {
    res += "~" + s;
  }

  return res;
}


// useful when we need to put component types to a set
bool ComponentType::operator < (const ComponentType& other) const {
  if (name == other.name) {
    if (!__eq__(other)) {
      throw RuntimeError(
          "Cannot define ordering (less) between " + to_bngl_str() + " and " + other.to_bngl_str() + ", "
          "they have the same name but different states. " +
          "Error might have occurred due to a call to " + NAME_CLASS_ELEMENTARY_MOLECULE_TYPE + ".__eq__().");
    }
    return false;
  }
  else {
    return name < other.name;
  }
}


std::string ComponentType::get_canonical_name() const {
  std::string res;

  res = name;

  set<string> sorted(states.begin(), states.end());
  for (const string& s: states) {
    res += "~" + s;
  }

  return res;
}

} // namespace API
} // namespace MCell
