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

#include "bng_data_to_datamodel_converter.h"

#include "datamodel_defines.h"
#include "bng/bng.h"

using namespace std;
using namespace DMUtil;
using Json::Value;

namespace MCell {

BngDataToDatamodelConverter::BngDataToDatamodelConverter() {
  reset();

  // list from https://www.rapidtables.com/web/color/RGB_Color.html
  colors.push_back(Vec3(1, 0, 0));
  colors.push_back(Vec3(0, 1, 0));
  colors.push_back(Vec3(0, 0, 1));
  colors.push_back(Vec3(1, 1, 0));
  colors.push_back(Vec3(0, 1, 1));
  colors.push_back(Vec3(1, 0, 1));
  // ...
}


void BngDataToDatamodelConverter::reset() {
  bng_engine = nullptr;
  next_color_index = 0;
}


void BngDataToDatamodelConverter::to_data_model(
    const BNG::BNGEngine& bng_engine_, Value& mcell_node) {

  reset();

  bng_engine = &bng_engine_;

  convert_species(mcell_node);

  //convert_reactions(mcell_node);

  return;
}


Vec3 BngDataToDatamodelConverter::get_next_color() {
  Vec3 res = colors[next_color_index];
  next_color_index++;
  if (next_color_index == colors.size()) {
    next_color_index = 0;
  }
  return res;
}


void BngDataToDatamodelConverter::convert_species(Value& mcell_node) {
  //
  Value& define_molecules = mcell_node[KEY_DEFINE_MOLECULES];
  json_add_version(define_molecules, JSON_DM_VERSION_1638);
  Value& molecule_list = define_molecules[KEY_MOLECULE_LIST];

  for (const BNG::Species& s: bng_engine->get_all_species().get_species_vector()) {
    if (s.name == ALL_MOLECULES || s.name == ALL_VOLUME_MOLECULES || s.name == ALL_SURFACE_MOLECULES) {
      continue;
    }
    Value species_node;
    convert_single_species(s, species_node);
    molecule_list.append(species_node);
  }
}

void BngDataToDatamodelConverter::convert_single_species(const BNG::Species& s, Value& species_node) {

  json_add_version(species_node, JSON_DM_VERSION_1632);

  Value& display = species_node[KEY_DISPLAY];
  display[KEY_EMIT] = 0.0;
  Value& color_value = display[KEY_COLOR];
  Vec3 c = get_next_color();
  color_value.append(c.r);
  color_value.append(c.g);
  color_value.append(c.b);
  display[KEY_GLYPH] = VALUE_GLYPH_SPHERE_1;
  display[KEY_SCALE] = 1;

  species_node[KEY_BNGL_COMPONENT_LIST] = Value(Json::arrayValue); // empty array;
  species_node[KEY_MOL_BNGL_LABEL] = "";
  species_node[KEY_DESCRIPTION] = "";
  species_node[KEY_MOL_TYPE] = s.is_vol() ? VALUE_MOL_TYPE_3D : VALUE_MOL_TYPE_2D;
  species_node[KEY_EXPORT_VIZ] = false; // true causes error in cellblender
  species_node[KEY_CUSTOM_SPACE_STEP] = "";
  species_node[KEY_MAXIMUM_STEP_LENGTH] = "";
  species_node[KEY_TARGET_ONLY] = false;
  species_node[KEY_DIFFUSION_CONSTANT] = to_string(s.D);
  species_node[KEY_SPATIAL_STRUCTURE] = "None";
  species_node[KEY_MOL_NAME] = s.name;
  species_node[KEY_CUSTOM_TIME_STEP] = "";
}


} // namespace MCell

