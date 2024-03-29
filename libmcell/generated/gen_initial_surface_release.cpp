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
#include "gen_initial_surface_release.h"
#include "api/initial_surface_release.h"
#include "api/complex.h"

namespace MCell {
namespace API {

void GenInitialSurfaceRelease::check_semantics() const {
  if (!is_set(complex)) {
    throw ValueError("Parameter 'complex' must be set.");
  }
}

void GenInitialSurfaceRelease::set_initialized() {
  if (is_set(complex)) {
    complex->set_initialized();
  }
  initialized = true;
}

void GenInitialSurfaceRelease::set_all_attributes_as_default_or_unset() {
  class_name = "InitialSurfaceRelease";
  complex = nullptr;
  number_to_release = INT_UNSET;
  density = FLT_UNSET;
}

std::shared_ptr<InitialSurfaceRelease> GenInitialSurfaceRelease::copy_initial_surface_release() const {
  std::shared_ptr<InitialSurfaceRelease> res = std::make_shared<InitialSurfaceRelease>(DefaultCtorArgType());
  res->class_name = class_name;
  res->complex = complex;
  res->number_to_release = number_to_release;
  res->density = density;

  return res;
}

std::shared_ptr<InitialSurfaceRelease> GenInitialSurfaceRelease::deepcopy_initial_surface_release(py::dict) const {
  std::shared_ptr<InitialSurfaceRelease> res = std::make_shared<InitialSurfaceRelease>(DefaultCtorArgType());
  res->class_name = class_name;
  res->complex = is_set(complex) ? complex->deepcopy_complex() : nullptr;
  res->number_to_release = number_to_release;
  res->density = density;

  return res;
}

bool GenInitialSurfaceRelease::__eq__(const InitialSurfaceRelease& other) const {
  return
    (
      (is_set(complex)) ?
        (is_set(other.complex) ?
          (complex->__eq__(*other.complex)) : 
          false
        ) :
        (is_set(other.complex) ?
          false :
          true
        )
     )  &&
    number_to_release == other.number_to_release &&
    density == other.density;
}

bool GenInitialSurfaceRelease::eq_nonarray_attributes(const InitialSurfaceRelease& other, const bool ignore_name) const {
  return
    (
      (is_set(complex)) ?
        (is_set(other.complex) ?
          (complex->__eq__(*other.complex)) : 
          false
        ) :
        (is_set(other.complex) ?
          false :
          true
        )
     )  &&
    number_to_release == other.number_to_release &&
    density == other.density;
}

std::string GenInitialSurfaceRelease::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "complex=" << "(" << ((complex != nullptr) ? complex->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "number_to_release=" << number_to_release << ", " <<
      "density=" << density;
  return ss.str();
}

py::class_<InitialSurfaceRelease> define_pybinding_InitialSurfaceRelease(py::module& m) {
  return py::class_<InitialSurfaceRelease, std::shared_ptr<InitialSurfaceRelease>>(m, "InitialSurfaceRelease", "Defines molecules to be released onto a SurfaceRegion right when simulation starts")
      .def(
          py::init<
            std::shared_ptr<Complex>,
            const int,
            const double
          >(),
          py::arg("complex"),
          py::arg("number_to_release") = INT_UNSET,
          py::arg("density") = FLT_UNSET
      )
      .def("check_semantics", &InitialSurfaceRelease::check_semantics)
      .def("__copy__", &InitialSurfaceRelease::copy_initial_surface_release)
      .def("__deepcopy__", &InitialSurfaceRelease::deepcopy_initial_surface_release, py::arg("memo"))
      .def("__str__", &InitialSurfaceRelease::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &InitialSurfaceRelease::__eq__, py::arg("other"))
      .def("dump", &InitialSurfaceRelease::dump)
      .def_property("complex", &InitialSurfaceRelease::get_complex, &InitialSurfaceRelease::set_complex, "Defines the species of the molecule that will be released.\n")
      .def_property("number_to_release", &InitialSurfaceRelease::get_number_to_release, &InitialSurfaceRelease::set_number_to_release, "Number of molecules to be released onto a region,\nonly one of number_to_release and density can be set.\n")
      .def_property("density", &InitialSurfaceRelease::get_density, &InitialSurfaceRelease::set_density, "Density of molecules to be released onto a region,\nonly one of number_to_release and density can be set.\n")
    ;
}

std::string GenInitialSurfaceRelease::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "initial_surface_release_" + std::to_string(ctx.postinc_counter("initial_surface_release"));
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
  ss << "m.InitialSurfaceRelease(" << nl;
  ss << ind << "complex = " << complex->export_to_python(out, ctx) << "," << nl;
  if (number_to_release != INT_UNSET) {
    ss << ind << "number_to_release = " << number_to_release << "," << nl;
  }
  if (density != FLT_UNSET) {
    ss << ind << "density = " << f_to_str(density) << "," << nl;
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

