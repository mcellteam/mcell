/******************************************************************************
 *
 * Copyright (C) 2020 by
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


#include "species_info.h"
#include "datamodel_defines.h"

using namespace std;

namespace MCell {

void SpeciesInfo::to_data_model(Json::Value& mcell) const {
  Json::Value& define_molecules = mcell[KEY_DEFINE_MOLECULES];
  json_add_version(define_molecules, JSON_DM_VERSION_1638);
  Json::Value& molecule_list = define_molecules[KEY_MOLECULE_LIST];

  for (const Species &s: species) {
    if (s.name == ALL_MOLECULES || s.name == ALL_VOLUME_MOLECULES || s.name == ALL_SURFACE_MOLECULES) {
      continue;
    }
    Json::Value species_value;
    s.to_data_model(species_value);
    molecule_list.append(species_value);
  }
}

} // namespace mcell
