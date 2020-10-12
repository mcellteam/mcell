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

#include <api/complex.h>
#include "api/species.h"
#include "api/elementary_molecule_instance.h"
#include "api/elementary_molecule_type.h"

using namespace std;

namespace MCell {
namespace API {


bool Complex::is_surf() const {
  for (auto em: elementary_molecule_instances) {
    if (is_set(em->elementary_molecule_type->diffusion_constant_2d)) {
      return true;
    }
  }
  return false;
}


std::string Complex::to_bngl_str() {
  if (is_set(name)) {
    return name;
  }
  else {
    std::string res;
    for (size_t i = 0; i < elementary_molecule_instances.size(); i++) {
      res += elementary_molecule_instances[i]->to_bngl_str();
      if (i + 1 != elementary_molecule_instances.size()) {
        res += ".";
      }
    }

    if (orientation == Orientation::UP) {
      res += "'";
    }
    else if (orientation == Orientation::DOWN) {
      res += ",";
    }

    return res;
  }
}


bool Complex::is_species_object() const {
  return dynamic_cast<const Species*>(this) != nullptr;
}


std::shared_ptr<Species> Complex::as_species() {
  return std::make_shared<Species>(*this);
}


} // namespace API
} // namespace MCell
