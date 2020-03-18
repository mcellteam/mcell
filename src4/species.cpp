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

#define ONE_DECIMAL 1 // DMFIXME: this is wrong, we will cut off too many decimal places

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

void Species::to_data_model(std::ostream& out) const{
  // stores default ostream output precision
  int old_precision = out.precision();

  // DMFIXME: remove setprecision

  out <<
      "{" <<
      "\n\"display\": {" <<
      "\n\"emit\": " <<
      "1.0," <<
      "\n\"color\": [";
  // configure ostream for 1 decimanl place output of floating points
  out << setprecision(ONE_DECIMAL) << fixed;
  out <<
      "\n" << color.r << "," <<
      "\n" << color.g << "," <<
      "\n" << color.b;
  // return ostream output to default precision settings
  out << setprecision(old_precision) << defaultfloat;
  out <<
      "\n]," <<
      "\n\"glyph\": " <<
      "\"" <<
      "Sphere_1" << // Sphere_1 by default for now.
      "\",";
      /*
       * This should only come up if the glyph is called "Letter"
      "\n\"letter\": " <<
      "\"\",";
      */
  out << setprecision(ONE_DECIMAL) << fixed;
  out <<
      "\n\"scale\": " <<
      scale;
  out << setprecision(old_precision) << defaultfloat;
  out <<
      "\n}," <<
      "\n\"bngl_component_list\": " <<
      "\[]," <<  // No bngl_component_list for now
      "\n\"mol_bngl_label\": " <<
      "\"\"," <<
      "\n\"data_model_version\": " <<
      JSON_DM_VERSION_1632 <<
      "," <<
      "\n\"description\": " <<
      "\"\"," <<
      "\n\"mol_type\": ";
  if (is_vol()) {
    out << "\"3D\",";
  }
  else {
    out << "\"2D\",";
  }
  out <<
      "\n\"export_viz\": " <<
      "false," << // or else errors in cellblender
      "\n\"custom_space_step\": " <<
      "\"\"," <<
      "\n\"maximum_step_length\": " <<
      "\"\"," <<
      "\n\"target_only\": " <<
      "false," << // or else errors in cellblender
      "\n\"diffusion_constant\": " <<
      "\"" << D << "\"," <<
      "\n\"spatial_structure\": " <<
      // this will always be "None" until mcell4 updates this feature
      "\"None\"," <<
      "\n\"mol_name\": " <<
      "\"" << name << "\"," <<
      "\n\"custom_time_step\": " <<
      "\"\"" <<
      "\n}";
}

// sets mol display color
void Species::set_color(float_t r, float_t g, float_t b) {
  assert((r <= 1 && r >= 0) && (g <= 1 && g >= 0) && (b <= 0 && b >= 0));
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
