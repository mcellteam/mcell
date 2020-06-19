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

using namespace std;

namespace MCell {
namespace API {

std::string ElementaryMoleculeInstance::to_bngl_str() {
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
