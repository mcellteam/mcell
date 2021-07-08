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

#include "api/component.h"
#include "api/component_type.h"


using namespace std;

namespace MCell {
namespace API {

bool Component::operator < (const Component& other) const {
  if (name != other.name) {
    return name < other.name;
  }

  if (state != other.state) {
    return state < other.state;
  }

  return bond < other.bond;
}


std::string Component::to_bngl_str() const {
  std::string res;

  res = component_type->name;

  if (state != STATE_UNSET) {
    res += "~" + state;
  }

  if (bond != BOND_UNBOUND) {
    if (bond == BOND_BOUND) {
      res += "!+";
    }
    else if (bond == BOND_ANY) {
      res += "!?";
    }
    else {
      res += "!" + std::to_string(bond);
    }
  }

  return res;
}

} // namespace API
} // namespace MCell
