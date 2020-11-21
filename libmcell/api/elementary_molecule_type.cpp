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

#include "api/elementary_molecule_type.h"
#include "api/component_type.h"

using namespace std;

namespace MCell {
namespace API {


static std::string get_components_str(
    const std::vector<std::shared_ptr<ComponentType>>& components,
    const bool canonical = false
) {
  string res;
  if (!components.empty()) {
    res += "(";
    for (size_t i = 0; i < components.size(); i++) {
      if (!canonical) {
        res += components[i]->to_bngl_str();
      }
      else {
        res += components[i]->get_canonical_name();
      }
      if (i + 1 != components.size()) {
        res += ",";
      }
    }
    res += ")";
  }
  return res;
}


std::string ElementaryMoleculeType::get_canonical_name() const {
  std::vector<std::shared_ptr<ComponentType>> sorted;
  sorted = components;
  std::sort(sorted.begin(), sorted.end(),
      [](const std::shared_ptr<ComponentType>& a, const std::shared_ptr<ComponentType>& b) -> bool {
          return *a < *b;
      });
  return name + get_components_str(sorted);
}


bool ElementaryMoleculeType::__eq__(const ElementaryMoleculeType& other) const {

  if (!eq_nonarray_attributes(other)) {
    return false;
  }

  return get_canonical_name() == other.get_canonical_name();
}


std::string ElementaryMoleculeType::to_bngl_str() const {
  return name + get_components_str(components);
}


} // namespace API
} // namespace MCell
