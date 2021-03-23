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
#include "gen_release_site.h"
#include "api/release_site.h"
#include "api/complex.h"
#include "api/molecule_release_info.h"
#include "api/region.h"
#include "api/release_pattern.h"

namespace MCell {
namespace API {

void GenReleaseSite::check_semantics() const {
  if (!is_set(name)) {
    throw ValueError("Parameter 'name' must be set.");
  }
}

void GenReleaseSite::set_initialized() {
  if (is_set(complex)) {
    complex->set_initialized();
  }
  vec_set_initialized(molecule_list);
  if (is_set(release_pattern)) {
    release_pattern->set_initialized();
  }
  if (is_set(region)) {
    region->set_initialized();
  }
  initialized = true;
}

void GenReleaseSite::set_all_attributes_as_default_or_unset() {
  class_name = "ReleaseSite";
  name = STR_UNSET;
  complex = nullptr;
  molecule_list = std::vector<std::shared_ptr<MoleculeReleaseInfo>>();
  release_time = 0;
  release_pattern = nullptr;
  shape = Shape::UNSET;
  region = nullptr;
  location = std::vector<float_t>();
  site_diameter = 0;
  site_radius = FLT_UNSET;
  number_to_release = FLT_UNSET;
  density = FLT_UNSET;
  concentration = FLT_UNSET;
}

bool GenReleaseSite::__eq__(const ReleaseSite& other) const {
  return
    name == other.name &&
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
    vec_ptr_eq(molecule_list, other.molecule_list) &&
    release_time == other.release_time &&
    (
      (is_set(release_pattern)) ?
        (is_set(other.release_pattern) ?
          (release_pattern->__eq__(*other.release_pattern)) : 
          false
        ) :
        (is_set(other.release_pattern) ?
          false :
          true
        )
     )  &&
    shape == other.shape &&
    (
      (is_set(region)) ?
        (is_set(other.region) ?
          (region->__eq__(*other.region)) : 
          false
        ) :
        (is_set(other.region) ?
          false :
          true
        )
     )  &&
    location == other.location &&
    site_diameter == other.site_diameter &&
    site_radius == other.site_radius &&
    number_to_release == other.number_to_release &&
    density == other.density &&
    concentration == other.concentration;
}

bool GenReleaseSite::eq_nonarray_attributes(const ReleaseSite& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
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
    true /*molecule_list*/ &&
    release_time == other.release_time &&
    (
      (is_set(release_pattern)) ?
        (is_set(other.release_pattern) ?
          (release_pattern->__eq__(*other.release_pattern)) : 
          false
        ) :
        (is_set(other.release_pattern) ?
          false :
          true
        )
     )  &&
    shape == other.shape &&
    (
      (is_set(region)) ?
        (is_set(other.region) ?
          (region->__eq__(*other.region)) : 
          false
        ) :
        (is_set(other.region) ?
          false :
          true
        )
     )  &&
    true /*location*/ &&
    site_diameter == other.site_diameter &&
    site_radius == other.site_radius &&
    number_to_release == other.number_to_release &&
    density == other.density &&
    concentration == other.concentration;
}

std::string GenReleaseSite::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "complex=" << "(" << ((complex != nullptr) ? complex->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "molecule_list=" << vec_ptr_to_str(molecule_list, ind + "  ") << ", " << "\n" << ind + "  " <<
      "release_time=" << release_time << ", " <<
      "\n" << ind + "  " << "release_pattern=" << "(" << ((release_pattern != nullptr) ? release_pattern->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "shape=" << shape << ", " <<
      "\n" << ind + "  " << "region=" << "(" << ((region != nullptr) ? region->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "location=" << vec_nonptr_to_str(location, ind + "  ") << ", " <<
      "site_diameter=" << site_diameter << ", " <<
      "site_radius=" << site_radius << ", " <<
      "number_to_release=" << number_to_release << ", " <<
      "density=" << density << ", " <<
      "concentration=" << concentration;
  return ss.str();
}

py::class_<ReleaseSite> define_pybinding_ReleaseSite(py::module& m) {
  return py::class_<ReleaseSite, std::shared_ptr<ReleaseSite>>(m, "ReleaseSite", "Defines a release site that specifies where, when and how should molecules be released. \n")
      .def(
          py::init<
            const std::string&,
            std::shared_ptr<Complex>,
            const std::vector<std::shared_ptr<MoleculeReleaseInfo>>,
            const float_t,
            std::shared_ptr<ReleasePattern>,
            const Shape,
            std::shared_ptr<Region>,
            const std::vector<float_t>,
            const float_t,
            const float_t,
            const float_t,
            const float_t,
            const float_t
          >(),
          py::arg("name"),
          py::arg("complex") = nullptr,
          py::arg("molecule_list") = std::vector<std::shared_ptr<MoleculeReleaseInfo>>(),
          py::arg("release_time") = 0,
          py::arg("release_pattern") = nullptr,
          py::arg("shape") = Shape::UNSET,
          py::arg("region") = nullptr,
          py::arg("location") = std::vector<float_t>(),
          py::arg("site_diameter") = 0,
          py::arg("site_radius") = FLT_UNSET,
          py::arg("number_to_release") = FLT_UNSET,
          py::arg("density") = FLT_UNSET,
          py::arg("concentration") = FLT_UNSET
      )
      .def("check_semantics", &ReleaseSite::check_semantics)
      .def("__str__", &ReleaseSite::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ReleaseSite::__eq__, py::arg("other"))
      .def("dump", &ReleaseSite::dump)
      .def_property("name", &ReleaseSite::get_name, &ReleaseSite::set_name, "Name of the release site")
      .def_property("complex", &ReleaseSite::get_complex, &ReleaseSite::set_complex, "Defines the species of the molecule that will be released. Not used for the LIST shape. \nMust be set when molecule_list is empty and unset when molecule_list is not empty.\nOrientation of the complex instance is used to define orientation of the released molecule,\nwhen Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and\nsurface molecules are released with Orientation.UP.\nWhen compartment is specified, this sets shape to Shape.COMPARTMENT and the molecules are released \ninto the compartment.\n")
      .def_property("molecule_list", &ReleaseSite::get_molecule_list, &ReleaseSite::set_molecule_list, "Used for LIST shape release mode. \nOnly one of number_to_release, density, concentration or molecule_list can be set.\n")
      .def_property("release_time", &ReleaseSite::get_release_time, &ReleaseSite::set_release_time, "Specifies time in seconds when the release event is executed.\nIn case when a release pattern is used, this is the time of the first release.      \nEquivalent to MDL's RELEASE_PATTERN command DELAY.\n")
      .def_property("release_pattern", &ReleaseSite::get_release_pattern, &ReleaseSite::set_release_pattern, "Use the release pattern to define schedule of releases. \nThe default is to release the specified number of molecules at the set release_time. \n")
      .def_property("shape", &ReleaseSite::get_shape, &ReleaseSite::set_shape, "Defines how the molecules shoudl be released. \nSet automatically for these cases to the following values: \nregion is set - Shape.REGION_EXPR,\nregion is not set and complex uses a compartment - Shape.COMPARTMENT,\nmolecule_list is set - Shape.LIST,\nlocation is set - Shape.SPHERICAL.\n")
      .def_property("region", &ReleaseSite::get_region, &ReleaseSite::set_region, "Defines a volume or surface region where to release molecules. \nSetting it sets shape to Shape.REGION_EXPR. \n")
      .def_property("location", &ReleaseSite::get_location, &ReleaseSite::set_location, "Defines center of a sphere where to release molecules. \nSetting it sets shape to Shape.SPHERICAL.\n")
      .def_property("site_diameter", &ReleaseSite::get_site_diameter, &ReleaseSite::set_site_diameter, "For a geometrical release site, this releases molecules uniformly within\na radius r computed as site_diameter/2. \nUsed only when shape is Shape.SPHERICAL.\nMaximum one of site_diameter or site_radius may be set.\n")
      .def_property("site_radius", &ReleaseSite::get_site_radius, &ReleaseSite::set_site_radius, "For a geometrical release site, this releases molecules uniformly within\na radius site_radius.\nUsed only when shape is Shape.SPHERICAL.\nMaximum one of site_diameter or site_radius may be set.\n")
      .def_property("number_to_release", &ReleaseSite::get_number_to_release, &ReleaseSite::set_number_to_release, "Sets number of molecules to release. Cannot be set when shape is Shape.LIST. \nOnly one of number_to_release, density, concentration or molecule_list can be set.\nValue is truncated (floored) to an integer.\n")
      .def_property("density", &ReleaseSite::get_density, &ReleaseSite::set_density, "Unit is molecules per square micron (for surfaces). \nOnly one of number_to_release, density, concentration or molecule_list can be set.\nCannot be set when shape is Shape.LIST.\n")
      .def_property("concentration", &ReleaseSite::get_concentration, &ReleaseSite::set_concentration, "Unit is molar (moles per liter) for volumes.\nOnly one of number_to_release, density, concentration or molecule_list can be set.\nCannot be set when shape is Shape.LIST.\n")
    ;
}

std::string GenReleaseSite::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("release_site") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("release_site")));
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
  ss << "m.ReleaseSite(" << nl;
  ss << ind << "name = " << "'" << name << "'" << "," << nl;
  if (is_set(complex)) {
    ss << ind << "complex = " << complex->export_to_python(out, ctx) << "," << nl;
  }
  if (molecule_list != std::vector<std::shared_ptr<MoleculeReleaseInfo>>() && !skip_vectors_export()) {
    ss << ind << "molecule_list = " << export_vec_molecule_list(out, ctx, exported_name) << "," << nl;
  }
  if (release_time != 0) {
    ss << ind << "release_time = " << f_to_str(release_time) << "," << nl;
  }
  if (is_set(release_pattern)) {
    ss << ind << "release_pattern = " << release_pattern->export_to_python(out, ctx) << "," << nl;
  }
  if (shape != Shape::UNSET) {
    ss << ind << "shape = " << shape << "," << nl;
  }
  if (is_set(region)) {
    ss << ind << "region = " << region->export_to_python(out, ctx) << "," << nl;
  }
  if (location != std::vector<float_t>() && !skip_vectors_export()) {
    ss << ind << "location = " << export_vec_location(out, ctx, exported_name) << "," << nl;
  }
  if (site_diameter != 0) {
    ss << ind << "site_diameter = " << f_to_str(site_diameter) << "," << nl;
  }
  if (site_radius != FLT_UNSET) {
    ss << ind << "site_radius = " << f_to_str(site_radius) << "," << nl;
  }
  if (number_to_release != FLT_UNSET) {
    ss << ind << "number_to_release = " << f_to_str(number_to_release) << "," << nl;
  }
  if (density != FLT_UNSET) {
    ss << ind << "density = " << f_to_str(density) << "," << nl;
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

std::string GenReleaseSite::export_vec_molecule_list(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < molecule_list.size(); i++) {
    const auto& item = molecule_list[i];
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

std::string GenReleaseSite::export_vec_location(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
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

