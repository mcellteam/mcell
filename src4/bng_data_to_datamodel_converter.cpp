#if 0

temporary storage of a file removed during merge


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

#ifndef SRC4_SPECIES_H_
#define SRC4_SPECIES_H_

#include <string>
#include <vector>
#include "defines.h"

#include "mcell_structs.h"

namespace Json {
class Value;
}

namespace MCell {

// same as in mcell_structs but renamed to make sure it is used correctly
enum species_flag_t {

  SPECIES_FLAG_ON_GRID = ON_GRID, // 0x01
  SPECIES_FLAG_IS_SURFACE = IS_SURFACE, // 0x02 - not stored
  SPECIES_FLAG_CAN_VOLVOL = CAN_VOLVOL, // 0x10  - can vol vol react?
  SPECIES_FLAG_CAN_VOLSURF = CAN_VOLSURF, // 0x20
  SPECIES_FLAG_CAN_SURFSURF = CAN_SURFSURF, // 0x80
  SPECIES_FLAG_CANT_INITIATE = CANT_INITIATE, // 0x400 - not sure what to do with this yet
  SPECIES_FLAG_CAN_SURFSURFSURF = CAN_SURFSURFSURF, // 0x20000 - not supported
  SPECIES_FLAG_SET_MAX_STEP_LENGTH = SET_MAX_STEP_LENGTH, // 0x80000
  SPECIES_FLAG_CAN_REGION_BORDER = CAN_REGION_BORDER, // 0x100000
  SPECIES_FLAG_EXTERNAL_SPECIES = EXTERNAL_SPECIES // 0x400000 - not supported
};

/**
 * Holds information on one species type.
 */
class Species {
public:
  Species()
  : species_id(SPECIES_ID_INVALID), mcell3_species_id(0),
    D(FLT_INVALID), space_step(FLT_INVALID), time_step(FLT_INVALID),
    flags(0), color(1, 0, 0), scale(1)
    {
  }

  species_id_t species_id;

  uint mcell3_species_id;
  float_t D; // diffusion constant
  std::string name;

  // TODO: make private
  float_t space_step;
  float_t time_step; // in standard time

  uint flags;

  // DMFIXME: how to represent colors? what is the color format in datamodel?
  Vec3 color;  // mol color default is red
  float_t scale; // scale = 1 by default

  bool has_flag(species_flag_t flag) const {
    assert(__builtin_popcount(flag) == 1);
    return (flags & flag) != 0;
  }

  bool is_surf() const {
    return has_flag(SPECIES_FLAG_ON_GRID);
  }

  bool is_vol() const {
    return !has_flag(SPECIES_FLAG_ON_GRID);
  }

  bool is_reactive_surface() const {
    return has_flag(SPECIES_FLAG_IS_SURFACE);
  }

  // true if can interact with edge of an border
  bool can_interact_with_border() const {
    return has_flag(SPECIES_FLAG_CAN_REGION_BORDER);
  }

  bool can_diffuse() const {
    return D != DIFFUSION_CONSTANT_ZER0;
  }

  float_t get_time_step() const {
    return time_step;
  }

  float_t get_space_step() const {
    return space_step;
  }

  void dump(const std::string ind) const;
  static void dump_array(const std::vector<Species>& vec);

  void to_data_model(Json::Value& species) const;

  void set_color(float_t r, float_t g, float_t b);
  void set_scale(float_t s);
};

} // namespace mcell

#endif // SRC4_SPECIES_H_


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


#endif
