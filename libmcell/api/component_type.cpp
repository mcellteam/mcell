/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
