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
#include "libs/pybind11/include/pybind11/stl.h"
#include "api/python_export_utils.h"
#include "gen_elementary_molecule_type.h"
#include "api/elementary_molecule_type.h"
#include "api/component.h"
#include "api/component_type.h"
#include "api/elementary_molecule.h"

namespace MCell {
namespace API {

void GenElementaryMoleculeType::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

void GenElementaryMoleculeType::set_initialized() {
  vec_set_initialized(components);
  initialized = true;
}

void GenElementaryMoleculeType::set_all_attributes_as_default_or_unset() {
  class_name = "ElementaryMoleculeType";
  name = STR_UNSET;
  components = std::vector<std::shared_ptr<ComponentType>>();
  diffusion_constant_2d = FLT_UNSET;
  diffusion_constant_3d = FLT_UNSET;
  custom_time_step = FLT_UNSET;
  custom_space_step = FLT_UNSET;
  target_only = false;
}

bool GenElementaryMoleculeType::__eq__(const ElementaryMoleculeType& other) const {
  return
    name == other.name &&
    vec_ptr_eq(components, other.components) &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    custom_time_step == other.custom_time_step &&
    custom_space_step == other.custom_space_step &&
    target_only == other.target_only;
}

bool GenElementaryMoleculeType::eq_nonarray_attributes(const ElementaryMoleculeType& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    true /*components*/ &&
    diffusion_constant_2d == other.diffusion_constant_2d &&
    diffusion_constant_3d == other.diffusion_constant_3d &&
    custom_time_step == other.custom_time_step &&
    custom_space_step == other.custom_space_step &&
    target_only == other.target_only;
}

std::string GenElementaryMoleculeType::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "components=" << vec_ptr_to_str(components, ind + "  ") << ", " << "\n" << ind + "  " <<
      "diffusion_constant_2d=" << diffusion_constant_2d << ", " <<
      "diffusion_constant_3d=" << diffusion_constant_3d << ", " <<
      "custom_time_step=" << custom_time_step << ", " <<
      "custom_space_step=" << custom_space_step << ", " <<
      "target_only=" << target_only;
  return ss.str();
}

py::class_<ElementaryMoleculeType> define_pybinding_ElementaryMoleculeType(py::module& m) {
  return py::class_<ElementaryMoleculeType, std::shared_ptr<ElementaryMoleculeType>>(m, "ElementaryMoleculeType")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<ComponentType>>,
            const float_t,
            const float_t,
            const float_t,
            const float_t,
            const bool
          >(),
          py::arg("name"),
          py::arg("components") = std::vector<std::shared_ptr<ComponentType>>(),
          py::arg("diffusion_constant_2d") = FLT_UNSET,
          py::arg("diffusion_constant_3d") = FLT_UNSET,
          py::arg("custom_time_step") = FLT_UNSET,
          py::arg("custom_space_step") = FLT_UNSET,
          py::arg("target_only") = false
      )
      .def("check_semantics", &ElementaryMoleculeType::check_semantics)
      .def("__str__", &ElementaryMoleculeType::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ElementaryMoleculeType::__eq__, py::arg("other"))
      .def("inst", &ElementaryMoleculeType::inst, py::arg("components") = std::vector<std::shared_ptr<Component>>())
      .def("to_bngl_str", &ElementaryMoleculeType::to_bngl_str)
      .def("dump", &ElementaryMoleculeType::dump)
      .def_property("name", &ElementaryMoleculeType::get_name, &ElementaryMoleculeType::set_name)
      .def_property("components", &ElementaryMoleculeType::get_components, &ElementaryMoleculeType::set_components)
      .def_property("diffusion_constant_2d", &ElementaryMoleculeType::get_diffusion_constant_2d, &ElementaryMoleculeType::set_diffusion_constant_2d)
      .def_property("diffusion_constant_3d", &ElementaryMoleculeType::get_diffusion_constant_3d, &ElementaryMoleculeType::set_diffusion_constant_3d)
      .def_property("custom_time_step", &ElementaryMoleculeType::get_custom_time_step, &ElementaryMoleculeType::set_custom_time_step)
      .def_property("custom_space_step", &ElementaryMoleculeType::get_custom_space_step, &ElementaryMoleculeType::set_custom_space_step)
      .def_property("target_only", &ElementaryMoleculeType::get_target_only, &ElementaryMoleculeType::set_target_only)
    ;
}

std::string GenElementaryMoleculeType::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("elementary_molecule_type") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("elementary_molecule_type")));
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
  ss << "m.ElementaryMoleculeType(" << nl;
  ss << ind << "name = " << "'" << name << "'" << "," << nl;
  if (components != std::vector<std::shared_ptr<ComponentType>>() && !skip_vectors_export()) {
    ss << ind << "components = " << export_vec_components(out, ctx, exported_name) << "," << nl;
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

std::string GenElementaryMoleculeType::export_vec_components(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < components.size(); i++) {
    const auto& item = components[i];
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

