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
#include "gen_molecule_release_info.h"
#include "api/molecule_release_info.h"
#include "api/complex.h"

namespace MCell {
namespace API {

void GenMoleculeReleaseInfo::check_semantics() const {
  if (!is_set(complex)) {
    throw ValueError("Parameter 'complex' must be set.");
  }
  if (!is_set(location)) {
    throw ValueError("Parameter 'location' must be set and the value must not be an empty list.");
  }
}

void GenMoleculeReleaseInfo::set_initialized() {
  if (is_set(complex)) {
    complex->set_initialized();
  }
  initialized = true;
}

void GenMoleculeReleaseInfo::set_all_attributes_as_default_or_unset() {
  class_name = "MoleculeReleaseInfo";
  complex = nullptr;
  location = std::vector<double>();
}

std::shared_ptr<MoleculeReleaseInfo> GenMoleculeReleaseInfo::copy_molecule_release_info() const {
  std::shared_ptr<MoleculeReleaseInfo> res = std::make_shared<MoleculeReleaseInfo>(DefaultCtorArgType());
  res->class_name = class_name;
  res->complex = complex;
  res->location = location;

  return res;
}

std::shared_ptr<MoleculeReleaseInfo> GenMoleculeReleaseInfo::deepcopy_molecule_release_info(py::dict) const {
  std::shared_ptr<MoleculeReleaseInfo> res = std::make_shared<MoleculeReleaseInfo>(DefaultCtorArgType());
  res->class_name = class_name;
  res->complex = is_set(complex) ? complex->deepcopy_complex() : nullptr;
  res->location = location;

  return res;
}

bool GenMoleculeReleaseInfo::__eq__(const MoleculeReleaseInfo& other) const {
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
    location == other.location;
}

bool GenMoleculeReleaseInfo::eq_nonarray_attributes(const MoleculeReleaseInfo& other, const bool ignore_name) const {
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
    true /*location*/;
}

std::string GenMoleculeReleaseInfo::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "complex=" << "(" << ((complex != nullptr) ? complex->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "location=" << vec_nonptr_to_str(location, all_details, ind + "  ");
  return ss.str();
}

py::class_<MoleculeReleaseInfo> define_pybinding_MoleculeReleaseInfo(py::module& m) {
  return py::class_<MoleculeReleaseInfo, std::shared_ptr<MoleculeReleaseInfo>>(m, "MoleculeReleaseInfo", "Defines a pair (molecule, location). Used in ReleaseSite when its shape is Shape.LIST.\n")
      .def(
          py::init<
            std::shared_ptr<Complex>,
            const std::vector<double>
          >(),
          py::arg("complex"),
          py::arg("location")
      )
      .def("check_semantics", &MoleculeReleaseInfo::check_semantics)
      .def("__copy__", &MoleculeReleaseInfo::copy_molecule_release_info)
      .def("__deepcopy__", &MoleculeReleaseInfo::deepcopy_molecule_release_info, py::arg("memo"))
      .def("__str__", &MoleculeReleaseInfo::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &MoleculeReleaseInfo::__eq__, py::arg("other"))
      .def("dump", &MoleculeReleaseInfo::dump)
      .def_property("complex", &MoleculeReleaseInfo::get_complex, &MoleculeReleaseInfo::set_complex, "Complex instance defining the molecule that will be released.\nOrientation of the complex instance is used to define orientation of the released molecule,\nwhen Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and\nsurface molecules are released with Orientation.UP.\nCompartment must not be set because this specific release definition states the location.  \n")
      .def_property("location", &MoleculeReleaseInfo::get_location, &MoleculeReleaseInfo::set_location, py::return_value_policy::reference, "3D position where the molecule will be released. \nIf a molecule has a 2D diffusion constant, it will be\nplaced on the surface closest to the coordinate given. \nArgument must have exactly three floating point values [x, y, z].\n  \n")
    ;
}

std::string GenMoleculeReleaseInfo::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "molecule_release_info_" + std::to_string(ctx.postinc_counter("molecule_release_info"));
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
  ss << "m.MoleculeReleaseInfo(" << nl;
  ss << ind << "complex = " << complex->export_to_python(out, ctx) << "," << nl;
  ss << ind << "location = " << export_vec_location(out, ctx, exported_name) << "," << nl;
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
}

std::string GenMoleculeReleaseInfo::export_vec_location(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < location.size(); i++) {
    const auto& item = location[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << f_to_str(item) << ", ";
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

