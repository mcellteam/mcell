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

#include "api/elementary_molecule_instance.h"
#include "api/elementary_molecule_type.h"
#include "api/component_type.h"
#include "api/component_instance.h"
#include "api/complex.h"

using namespace std;

namespace MCell {
namespace API {


bool ElementaryMoleculeInstance::__eq__(const ElementaryMoleculeInstance& other) const {

  // do we have the same mol type?
  if (*elementary_molecule_type != *other.elementary_molecule_type) {
    return false;
  }

  // we must sort the components,
  // canonicalization of BNGL strings cannot be used because it maintains
  // the component ordering and it is not known when processing a single string
  std::set<ComponentInstance> s1;
  for (auto& c: components) {
    s1.insert(*c);
  }
  std::set<ComponentInstance> s2;
  for (auto& c: other.components) {
    s2.insert(*c);
  }
  return s1 == s2;
}


std::string ElementaryMoleculeInstance::to_bngl_str() const {
  std::string res;

  res = elementary_molecule_type->name;

  if (!components.empty()) {
    res += "(";
    for (size_t i = 0; i < components.size(); i++) {
      res += components[i]->to_bngl_str();
      if (i + 1 != components.size()) {
        res += ",";
      }
    }
    res += ")";
  }

  return res;
}

} // namespace API
} // namespace MCell
