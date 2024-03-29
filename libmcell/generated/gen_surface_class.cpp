/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_surface_class.h"
#include "api/surface_class.h"
#include "api/complex.h"
#include "api/surface_property.h"

namespace MCell {
namespace API {

void GenSurfaceClass::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

void GenSurfaceClass::set_initialized() {
  vec_set_initialized(properties);
  if (is_set(affected_complex_pattern)) {
    affected_complex_pattern->set_initialized();
  }
  initialized = true;
}

void GenSurfaceClass::set_all_attributes_as_default_or_unset() {
  class_name = "SurfaceClass";
  name = STR_UNSET;
  properties = std::vector<std::shared_ptr<SurfaceProperty>>();
  type = SurfacePropertyType::UNSET;
  affected_complex_pattern = nullptr;
  concentration = FLT_UNSET;
}

std::shared_ptr<SurfaceClass> GenSurfaceClass::copy_surface_class() const {
  std::shared_ptr<SurfaceClass> res = std::make_shared<SurfaceClass>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->properties = properties;
  res->type = type;
  res->affected_complex_pattern = affected_complex_pattern;
  res->concentration = concentration;

  return res;
}

std::shared_ptr<SurfaceClass> GenSurfaceClass::deepcopy_surface_class(py::dict) const {
  std::shared_ptr<SurfaceClass> res = std::make_shared<SurfaceClass>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  for (const auto& item: properties) {
    res->properties.push_back((is_set(item)) ? item->deepcopy_surface_property() : nullptr);
  }
  res->type = type;
  res->affected_complex_pattern = is_set(affected_complex_pattern) ? affected_complex_pattern->deepcopy_complex() : nullptr;
  res->concentration = concentration;

  return res;
}

bool GenSurfaceClass::__eq__(const SurfaceClass& other) const {
  return
    name == other.name &&
    vec_ptr_eq(properties, other.properties) &&
    type == other.type &&
    (
      (is_set(affected_complex_pattern)) ?
        (is_set(other.affected_complex_pattern) ?
          (affected_complex_pattern->__eq__(*other.affected_complex_pattern)) : 
          false
        ) :
        (is_set(other.affected_complex_pattern) ?
          false :
          true
        )
     )  &&
    concentration == other.concentration;
}

bool GenSurfaceClass::eq_nonarray_attributes(const SurfaceClass& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    true /*properties*/ &&
    type == other.type &&
    (
      (is_set(affected_complex_pattern)) ?
        (is_set(other.affected_complex_pattern) ?
          (affected_complex_pattern->__eq__(*other.affected_complex_pattern)) : 
          false
        ) :
        (is_set(other.affected_complex_pattern) ?
          false :
          true
        )
     )  &&
    concentration == other.concentration;
}

std::string GenSurfaceClass::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "properties=" << vec_ptr_to_str(properties, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "type=" << type << ", " <<
      "\n" << ind + "  " << "affected_complex_pattern=" << "(" << ((affected_complex_pattern != nullptr) ? affected_complex_pattern->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "concentration=" << concentration;
  return ss.str();
}

py::class_<SurfaceClass> define_pybinding_SurfaceClass(py::module& m) {
  return py::class_<SurfaceClass, SurfaceProperty, std::shared_ptr<SurfaceClass>>(m, "SurfaceClass", "Defining a surface class allows surfaces to behave like species. For instance, one may wish \nto specify that a surface does not block the diffusion of molecules. Each type of surface is defined\nby name, and each surface name must be unique in the simulation and should not match any molecule names.\nTo define a reaction with a surface class, use constructor call m.Complex(name) as one of the reactants.\n")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<SurfaceProperty>>,
            const SurfacePropertyType,
            std::shared_ptr<Complex>,
            const double
          >(),
          py::arg("name"),
          py::arg("properties") = std::vector<std::shared_ptr<SurfaceProperty>>(),
          py::arg("type") = SurfacePropertyType::UNSET,
          py::arg("affected_complex_pattern") = nullptr,
          py::arg("concentration") = FLT_UNSET
      )
      .def("check_semantics", &SurfaceClass::check_semantics)
      .def("__copy__", &SurfaceClass::copy_surface_class)
      .def("__deepcopy__", &SurfaceClass::deepcopy_surface_class, py::arg("memo"))
      .def("__str__", &SurfaceClass::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &SurfaceClass::__eq__, py::arg("other"))
      .def("dump", &SurfaceClass::dump)
      .def_property("name", &SurfaceClass::get_name, &SurfaceClass::set_name, "Name of the surface class.")
      .def_property("properties", &SurfaceClass::get_properties, &SurfaceClass::set_properties, py::return_value_policy::reference, "A surface class can either have a list of properties or just one property.\nIn the usual case of having one property, one can set the attributes \ntype, affected_species, etc. inherited from SurfaceProperty directly.\n")
    ;
}

std::string GenSurfaceClass::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("surface_class") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("surface_class")));
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
  ss << "m.SurfaceClass(" << nl;
  if (type != SurfacePropertyType::UNSET) {
    ss << ind << "type = " << type << "," << nl;
  }
  if (is_set(affected_complex_pattern)) {
    ss << ind << "affected_complex_pattern = " << affected_complex_pattern->export_to_python(out, ctx) << "," << nl;
  }
  if (concentration != FLT_UNSET) {
    ss << ind << "concentration = " << f_to_str(concentration) << "," << nl;
  }
  ss << ind << "name = " << "'" << name << "'" << "," << nl;
  if (properties != std::vector<std::shared_ptr<SurfaceProperty>>() && !skip_vectors_export()) {
    ss << ind << "properties = " << export_vec_properties(out, ctx, exported_name) << "," << nl;
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

std::string GenSurfaceClass::export_vec_properties(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < properties.size(); i++) {
    const auto& item = properties[i];
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

