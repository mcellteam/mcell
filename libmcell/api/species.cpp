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

#include "api/complex.h"
#include "api/species.h"
#include "api/elementary_molecule.h"
#include "api/elementary_molecule_type.h"
#include "api/python_export_utils.h"

using namespace std;

namespace MCell {
namespace API {

// TODO: how to make this consistent with API definition?
bool Species::__eq__(const Species& other) const {
  return
    name == other.name &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    custom_time_step == other.custom_time_step &&
    custom_space_step == other.custom_space_step &&
    target_only == other.target_only &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name &&

    Complex::__eq__(other) // make canonical comparison
    ;
}


// species name is sufficient
std::string Species::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }

  if (is_species_superclass()) {
    if (name == BNG::ALL_MOLECULES) {
      return S(MDOT) + NAME_CV_AllMolecules;
    }
    else if (name == BNG::ALL_VOLUME_MOLECULES) {
      return S(MDOT) + NAME_CV_AllVolumeMolecules;
    }
    else if (name == BNG::ALL_SURFACE_MOLECULES) {
      return S(MDOT) + NAME_CV_AllSurfaceMolecules;
    }
    else {
      assert(false);
    }
  }

  std::string exported_name = "species_";
  if (name.size() > MAX_SPECIES_NAME_LENGTH) {
    exported_name += to_string(ctx.postinc_counter("species"));
  }
  else {
    exported_name += fix_id(name);
  }
  ctx.add_exported(this, exported_name);

  std::stringstream ss;
  ss << exported_name << " = m.Species(\n";
  assert(is_set(name));
  ss << "  name = " << "'" << name << "'" << "\n";

  ss << ")\n\n";
  out << ss.str();
  return exported_name;
}

} // namespace API
} // namespace MCell
