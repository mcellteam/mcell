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

#include "api/elementary_molecule.h"
#include "api/elementary_molecule_type.h"
#include "api/component_type.h"
#include "api/component.h"
#include "api/complex.h"

using namespace std;

namespace MCell {
namespace API {


bool ElementaryMolecule::__eq__(const ElementaryMolecule& other) const {

  // do we have the same mol type?
  if (!eq_nonarray_attributes(other)) {
    return false;
  }

  // are components the same (order does not matter)
  std::set<Component> s1;
  for (auto& c: components) {
    s1.insert(*c);
  }
  std::set<Component> s2;
  for (auto& c: other.components) {
    s2.insert(*c);
  }
  return s1 == s2;
}


std::string ElementaryMolecule::to_bngl_str(const bool with_compartment) const {
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

  if (with_compartment && is_set(compartment_name)) {
    res += "@" + compartment_name;
  }

  return res;
}


bool ElementaryMolecule::is_surf() const {
  assert(is_set(elementary_molecule_type));
  return is_set(elementary_molecule_type->diffusion_constant_2d);
}

} // namespace API
} // namespace MCell
