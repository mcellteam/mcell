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
#include "gen_surface_property.h"
#include "api/surface_property.h"
#include "api/complex.h"

namespace MCell {
namespace API {

void GenSurfaceProperty::check_semantics() const {
}

void GenSurfaceProperty::set_initialized() {
  if (is_set(affected_complex_pattern)) {
    affected_complex_pattern->set_initialized();
  }
  initialized = true;
}

void GenSurfaceProperty::set_all_attributes_as_default_or_unset() {
  class_name = "SurfaceProperty";
  type = SurfacePropertyType::UNSET;
  affected_complex_pattern = nullptr;
  concentration = FLT_UNSET;
}

std::shared_ptr<SurfaceProperty> GenSurfaceProperty::copy_surface_property() const {
  std::shared_ptr<SurfaceProperty> res = std::make_shared<SurfaceProperty>(DefaultCtorArgType());
  res->class_name = class_name;
  res->type = type;
  res->affected_complex_pattern = affected_complex_pattern;
  res->concentration = concentration;

  return res;
}

std::shared_ptr<SurfaceProperty> GenSurfaceProperty::deepcopy_surface_property(py::dict) const {
  std::shared_ptr<SurfaceProperty> res = std::make_shared<SurfaceProperty>(DefaultCtorArgType());
  res->class_name = class_name;
  res->type = type;
  res->affected_complex_pattern = is_set(affected_complex_pattern) ? affected_complex_pattern->deepcopy_complex() : nullptr;
  res->concentration = concentration;

  return res;
}

bool GenSurfaceProperty::__eq__(const SurfaceProperty& other) const {
  return
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

bool GenSurfaceProperty::eq_nonarray_attributes(const SurfaceProperty& other, const bool ignore_name) const {
  return
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

std::string GenSurfaceProperty::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "type=" << type << ", " <<
      "\n" << ind + "  " << "affected_complex_pattern=" << "(" << ((affected_complex_pattern != nullptr) ? affected_complex_pattern->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "concentration=" << concentration;
  return ss.str();
}

py::class_<SurfaceProperty> define_pybinding_SurfaceProperty(py::module& m) {
  return py::class_<SurfaceProperty, std::shared_ptr<SurfaceProperty>>(m, "SurfaceProperty", "Single property for a SurfaceClass.")
      .def(
          py::init<
            const SurfacePropertyType,
            std::shared_ptr<Complex>,
            const double
          >(),
          py::arg("type") = SurfacePropertyType::UNSET,
          py::arg("affected_complex_pattern") = nullptr,
          py::arg("concentration") = FLT_UNSET
      )
      .def("check_semantics", &SurfaceProperty::check_semantics)
      .def("__copy__", &SurfaceProperty::copy_surface_property)
      .def("__deepcopy__", &SurfaceProperty::deepcopy_surface_property, py::arg("memo"))
      .def("__str__", &SurfaceProperty::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &SurfaceProperty::__eq__, py::arg("other"))
      .def("dump", &SurfaceProperty::dump)
      .def_property("type", &SurfaceProperty::get_type, &SurfaceProperty::set_type, "Must be set. See SurfacePropertyType for options.\n")
      .def_property("affected_complex_pattern", &SurfaceProperty::get_affected_complex_pattern, &SurfaceProperty::set_affected_complex_pattern, "A complex pattern with optional orientation must be set.\nDefault orientation means that the pattern matches any orientation.\nFor concentration or flux clamp the orientation specifies on which side  \nwill be the concentration held (UP is front or outside, DOWN is back or \ninside, and DEFAULT, ANY or NONE is on both sides).\nThe complex pattern must not use compartments.\n")
      .def_property("concentration", &SurfaceProperty::get_concentration, &SurfaceProperty::set_concentration, "Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or \nSurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side \nof the wall that has this concentration or flux clamped.\n")
    ;
}

std::string GenSurfaceProperty::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "surface_property_" + std::to_string(ctx.postinc_counter("surface_property"));
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
  ss << "m.SurfaceProperty(" << nl;
  if (type != SurfacePropertyType::UNSET) {
    ss << ind << "type = " << type << "," << nl;
  }
  if (is_set(affected_complex_pattern)) {
    ss << ind << "affected_complex_pattern = " << affected_complex_pattern->export_to_python(out, ctx) << "," << nl;
  }
  if (concentration != FLT_UNSET) {
    ss << ind << "concentration = " << f_to_str(concentration) << "," << nl;
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

} // namespace API
} // namespace MCell

