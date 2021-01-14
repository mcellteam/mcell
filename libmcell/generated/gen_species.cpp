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

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_species.h"
#include "api\species.h"
#include "api\complex.h"
#include "api\elementary_molecule.h"
#include "api\species.h"

namespace MCell {
namespace API {

void GenSpecies::check_semantics() const {
}

void GenSpecies::set_initialized() {
  vec_set_initialized(elementary_molecules);
  initialized = true;
}

void GenSpecies::set_all_attributes_as_default_or_unset() {
  class_name = "Species";
  name = STR_UNSET;
  diffusion_constant_2d = FLT_UNSET;
  diffusion_constant_3d = FLT_UNSET;
  custom_time_step = FLT_UNSET;
  custom_space_step = FLT_UNSET;
  target_only = false;
  name = STR_UNSET;
  elementary_molecules = std::vector<std::shared_ptr<ElementaryMolecule>>();
  orientation = Orientation::DEFAULT;
  compartment_name = STR_UNSET;
}

bool GenSpecies::__eq__(const Species& other) const {
  return
    name == other.name &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    custom_time_step == other.custom_time_step &&
    custom_space_step == other.custom_space_step &&
    target_only == other.target_only &&
    name == other.name &&
    vec_ptr_eq(elementary_molecules, other.elementary_molecules) &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name;
}

bool GenSpecies::eq_nonarray_attributes(const Species& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    custom_time_step == other.custom_time_step &&
    custom_space_step == other.custom_space_step &&
    target_only == other.target_only &&
    (ignore_name || name == other.name) &&
    true /*elementary_molecules*/ &&
    orientation == other.orientation &&
    compartment_name == other.compartment_name;
}

std::string GenSpecies::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "diffusion_constant_2d=" << diffusion_constant_2d << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d << ", " <<
      "custom_time_step=" << custom_time_step << ", " <<
      "custom_space_step=" << custom_space_step << ", " <<
      "target_only=" << target_only << ", " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "elementary_molecules=" << vec_ptr_to_str(elementary_molecules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation << ", " <<
      "compartment_name=" << compartment_name;
  return ss.str();
}

py::class_<Species> define_pybinding_Species(py::module& m) {
  return py::class_<Species, Complex, std::shared_ptr<Species>>(m, "Species")
      .def(
          py::init<
            const std::string&,
            const float_t,
            const float_t,
            const float_t,
            const float_t,
            const bool,
            const std::vector<std::shared_ptr<ElementaryMolecule>>,
            const Orientation,
            const std::string&
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("diffusion_constant_2d") = FLT_UNSET,
          py::arg("diffusion_constant_3d") = FLT_UNSET,
          py::arg("custom_time_step") = FLT_UNSET,
          py::arg("custom_space_step") = FLT_UNSET,
          py::arg("target_only") = false,
          py::arg("elementary_molecules") = std::vector<std::shared_ptr<ElementaryMolecule>>(),
          py::arg("orientation") = Orientation::DEFAULT,
          py::arg("compartment_name") = STR_UNSET
      )
      .def("check_semantics", &Species::check_semantics)
      .def("__str__", &Species::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Species::__eq__, py::arg("other"))
      .def("inst", &Species::inst, py::arg("orientation") = Orientation::DEFAULT, py::arg("compartment_name") = STR_UNSET)
      .def("dump", &Species::dump)
      .def_property("name", &Species::get_name, &Species::set_name)
      .def_property("diffusion_constant_2d", &Species::get_diffusion_constant_2d, &Species::set_diffusion_constant_2d)
      .def_property("diffusion_constant_3d", &Species::get_diffusion_constant_3d, &Species::set_diffusion_constant_3d)
      .def_property("custom_time_step", &Species::get_custom_time_step, &Species::set_custom_time_step)
      .def_property("custom_space_step", &Species::get_custom_space_step, &Species::set_custom_space_step)
      .def_property("target_only", &Species::get_target_only, &Species::set_target_only)
    ;
}

std::string GenSpecies::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("species") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("species")));
  if (!export_even_if_already_exported()) {
    ctx.add_exported(this, exported_name);
  }

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.Species(" << nl;
  if (name != STR_UNSET) {
    ss << ind << "name = " << "'" << name << "'" << "," << nl;
  }
  if (elementary_molecules != std::vector<std::shared_ptr<ElementaryMolecule>>() && !skip_vectors_export()) {
    ss << ind << "elementary_molecules = " << export_vec_elementary_molecules(out, ctx, exported_name) << "," << nl;
  }
  if (orientation != Orientation::DEFAULT) {
    ss << ind << "orientation = " << orientation << "," << nl;
  }
  if (compartment_name != STR_UNSET) {
    ss << ind << "compartment_name = " << "'" << compartment_name << "'" << "," << nl;
  }
  if (diffusion_constant_2d != FLT_UNSET) {
    ss << ind << "diffusion_constant_2d = " << f_to_str(diffusion_constant_2d) << "," << nl;
  }
  if (diffusion_constant_3d != FLT_UNSET) {
    ss << ind << "diffusion_constant_3d = " << f_to_str(diffusion_constant_3d) << "," << nl;
  }
  if (custom_time_step != FLT_UNSET) {
    ss << ind << "custom_time_step = " << f_to_str(custom_time_step) << "," << nl;
  }
  if (custom_space_step != FLT_UNSET) {
    ss << ind << "custom_space_step = " << f_to_str(custom_space_step) << "," << nl;
  }
  if (target_only != false) {
    ss << ind << "target_only = " << target_only << "," << nl;
  }
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

std::string GenSpecies::export_vec_elementary_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < elementary_molecules.size(); i++) {
    const auto& item = elementary_molecules[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

