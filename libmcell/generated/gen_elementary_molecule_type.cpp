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

std::shared_ptr<ElementaryMoleculeType> GenElementaryMoleculeType::copy_elementary_molecule_type() const {
  if (initialized) {
    throw RuntimeError("Object of class ElementaryMoleculeType cannot be cloned with 'copy' after this object was used in model initialization.");
  }

  std::shared_ptr<ElementaryMoleculeType> res = std::make_shared<ElementaryMoleculeType>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->components = components;
  res->diffusion_constant_2d = diffusion_constant_2d;
  res->diffusion_constant_3d = diffusion_constant_3d;
  res->custom_time_step = custom_time_step;
  res->custom_space_step = custom_space_step;
  res->target_only = target_only;

  return res;
}

std::shared_ptr<ElementaryMoleculeType> GenElementaryMoleculeType::deepcopy_elementary_molecule_type(py::dict) const {
  if (initialized) {
    throw RuntimeError("Object of class ElementaryMoleculeType cannot be cloned with 'deepcopy' after this object was used in model initialization.");
  }

  std::shared_ptr<ElementaryMoleculeType> res = std::make_shared<ElementaryMoleculeType>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  for (const auto& item: components) {
    res->components.push_back((is_set(item)) ? item->deepcopy_component_type() : nullptr);
  }
  res->diffusion_constant_2d = diffusion_constant_2d;
  res->diffusion_constant_3d = diffusion_constant_3d;
  res->custom_time_step = custom_time_step;
  res->custom_space_step = custom_space_step;
  res->target_only = target_only;

  return res;
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
  return py::class_<ElementaryMoleculeType, std::shared_ptr<ElementaryMoleculeType>>(m, "ElementaryMoleculeType", "An elementary molecule type is a base indivisible entity. It is the same as  \na molecule type in BNGL entered in section molecule types. \nThe 'elementary' prefix was added to distinguish it clearly from molecules in \nsimulation.\n")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<ComponentType>>,
            const double,
            const double,
            const double,
            const double,
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
      .def("__copy__", &ElementaryMoleculeType::copy_elementary_molecule_type)
      .def("__deepcopy__", &ElementaryMoleculeType::deepcopy_elementary_molecule_type, py::arg("memo"))
      .def("__str__", &ElementaryMoleculeType::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ElementaryMoleculeType::__eq__, py::arg("other"))
      .def("inst", &ElementaryMoleculeType::inst, py::arg("components") = std::vector<std::shared_ptr<Component>>(), py::arg("compartment_name") = STR_UNSET, "Create an elementary molecule based on this elementary molecule type.\n- components: Instances of components for the the created elementary molecule.\nNot all components need to be specified in case when the elementary \nmolecule is used in a pattern.\n \n\n\n- compartment_name: Optional specification of compartment name for the created elementary molecule. \n\n\n")
      .def("to_bngl_str", &ElementaryMoleculeType::to_bngl_str, "Creates a string that corresponds to its BNGL representation.")
      .def("dump", &ElementaryMoleculeType::dump)
      .def_property("name", &ElementaryMoleculeType::get_name, &ElementaryMoleculeType::set_name, "Name of this elementary molecule type.")
      .def_property("components", &ElementaryMoleculeType::get_components, &ElementaryMoleculeType::set_components, py::return_value_policy::reference, "List of components used by this elementary molecule type.")
      .def_property("diffusion_constant_2d", &ElementaryMoleculeType::get_diffusion_constant_2d, &ElementaryMoleculeType::set_diffusion_constant_2d, "Elementary molecule based on this type is constrained to a surface\nand diffuses with the specified diffusion constant.\nD can be zero, in which case the molecule doesn’t move. \nThe units of D are cm^2 /s.\n")
      .def_property("diffusion_constant_3d", &ElementaryMoleculeType::get_diffusion_constant_3d, &ElementaryMoleculeType::set_diffusion_constant_3d, "Elementary molecule based on this type diffuses in space with the \nspecified diffusion constant D. \nD can be zero, in which case the molecule doesn’t move. \nThe units of D are cm^2 /s.\n")
      .def_property("custom_time_step", &ElementaryMoleculeType::get_custom_time_step, &ElementaryMoleculeType::set_custom_time_step, "This molecule should take timesteps of length custom_time_step (in seconds). \nUse either this or custom_time_step, not both.\n")
      .def_property("custom_space_step", &ElementaryMoleculeType::get_custom_space_step, &ElementaryMoleculeType::set_custom_space_step, "This molecule should take steps of average length given by the custom_space_step value (in microns). \nUse either this or custom_time_step, not both.\n")
      .def_property("target_only", &ElementaryMoleculeType::get_target_only, &ElementaryMoleculeType::set_target_only, "This molecule will not initiate reactions when it runs into other molecules. This\nsetting can speed up simulations when applied to a molecule at high concentrations \nthat reacts with a molecule at low concentrations (it is more efficient for\nthe low-concentration molecule to trigger the reactions). This directive does\nnot affect unimolecular reactions.      \n")
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

