/******************************************************************************
 *
 * Copyright (C) 2019 by
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

#include <iostream>
#include <iomanip> // needed for std::setprecision()
#include "mdlparse_util.h"

#include "species.h"
#include "datamodel_defines.h"

using namespace std;

namespace MCell {

void Species::dump(const string ind) const {
  cout << ind <<"species_id: \t\t" << species_id << " [uint16_t] \t\t/* Unique ID for this species */\n";
  cout << ind <<"mcell_species_id: \t\t" << mcell3_species_id << " [uint] \t\t/* Unique ID for this species from mcell3 representation*/\n";
  cout << ind <<"name: *\t\t" << name << " [string] \t\t/* Symbol table entry (name) */\n";
  cout << ind <<"D: \t\t" << D << " [float_t] \t\t/* Diffusion constant */\n";
  cout << ind <<"space_step: \t\t" << space_step << " [float_t] \t\t/* Characteristic step length */\n";
  cout << ind <<"time_step: \t\t" << time_step << " [float_t] \t\t/* Minimum (maximum?) sensible timestep */\n";
}

void Species::dump_array(const vector<Species>& vec) {
  cout << "Species array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump("  ");
  }
}

void Species::to_data_model(Json::Value& species) const{

  json_add_version(species, JSON_DM_VERSION_1632);

  Json::Value& display = species[KEY_DISPLAY];
  display[KEY_EMIT] = 0.0;
  Json::Value& color_value = display[KEY_COLOR];
  color_value.append(color.r);
  color_value.append(color.g);
  color_value.append(color.b);
  display[KEY_GLYPH] = VALUE_GLYPH_SPHERE_1;
  display[KEY_SCALE] = scale;


  species[KEY_BNGL_COMPONENT_LIST] = Json::Value(Json::arrayValue); // empty array;
  species[KEY_MOL_BNGL_LABEL] = "";
  species[KEY_DESCRIPTION] = "";
  species[KEY_MOL_TYPE] = is_vol() ? VALUE_MOL_TYPE_3D : VALUE_MOL_TYPE_2D;
  species[KEY_EXPORT_VIZ] = false; // true causes error in cellblender
  species[KEY_CUSTOM_SPACE_STEP] = "";
  species[KEY_MAXIMUM_STEP_LENGTH] = "";
  species[KEY_TARGET_ONLY] = false;
  species[KEY_DIFFUSION_CONSTANT] = to_string(D);
  species[KEY_SPATIAL_STRUCTURE] = "None";
  species[KEY_MOL_NAME] = name;
  species[KEY_CUSTOM_TIME_STEP] = "";
}

// sets mol display color
void Species::set_color(float_t r, float_t g, float_t b) {
  assert((r <= 1 && r >= 0) && (g <= 1 && g >= 0) && (b <= 1 && b >= 0));
  color.r = r;
  color.g = g;
  color.b = b;
}

// sets mol display size
void Species::set_scale(float_t s) {
  assert(s > 0);
  scale = s;
}

} // namespace mcell
