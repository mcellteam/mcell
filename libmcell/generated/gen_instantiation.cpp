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
#include "gen_instantiation.h"
#include "api/instantiation.h"
#include "api/base_chkpt_mol.h"
#include "api/geometry_object.h"
#include "api/region.h"
#include "api/release_site.h"

namespace MCell {
namespace API {

bool GenInstantiation::__eq__(const Instantiation& other) const {
  return
    vec_ptr_eq(release_sites, other.release_sites) &&
    vec_ptr_eq(geometry_objects, other.geometry_objects) &&
    vec_ptr_eq(checkpointed_molecules, other.checkpointed_molecules);
}

bool GenInstantiation::eq_nonarray_attributes(const Instantiation& other, const bool ignore_name) const {
  return
    true /*release_sites*/ &&
    true /*geometry_objects*/ &&
    true /*checkpointed_molecules*/;
}

std::string GenInstantiation::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Instantiation" << ": " <<
      "\n" << ind + "  " << "release_sites=" << vec_ptr_to_str(release_sites, ind + "  ") << ", " << "\n" << ind + "  " <<
      "geometry_objects=" << vec_ptr_to_str(geometry_objects, ind + "  ") << ", " << "\n" << ind + "  " <<
      "checkpointed_molecules=" << vec_ptr_to_str(checkpointed_molecules, ind + "  ");
  return ss.str();
}

py::class_<Instantiation> define_pybinding_Instantiation(py::module& m) {
  return py::class_<Instantiation, std::shared_ptr<Instantiation>>(m, "Instantiation")
      .def(
          py::init<
            const std::vector<std::shared_ptr<ReleaseSite>>,
            const std::vector<std::shared_ptr<GeometryObject>>,
            const std::vector<std::shared_ptr<BaseChkptMol>>
          >(),
          py::arg("release_sites") = std::vector<std::shared_ptr<ReleaseSite>>(),
          py::arg("geometry_objects") = std::vector<std::shared_ptr<GeometryObject>>(),
          py::arg("checkpointed_molecules") = std::vector<std::shared_ptr<BaseChkptMol>>()
      )
      .def("__str__", &Instantiation::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Instantiation::__eq__, py::arg("other"))
      .def("add_release_site", &Instantiation::add_release_site, py::arg("s"))
      .def("find_release_site", &Instantiation::find_release_site, py::arg("name"))
      .def("add_geometry_object", &Instantiation::add_geometry_object, py::arg("o"))
      .def("find_geometry_object", &Instantiation::find_geometry_object, py::arg("name"))
      .def("find_volume_compartment", &Instantiation::find_volume_compartment, py::arg("name"))
      .def("find_surface_compartment", &Instantiation::find_surface_compartment, py::arg("name"))
      .def("load_bngl_seed_species", &Instantiation::load_bngl_seed_species, py::arg("file_name"), py::arg("default_release_region") = nullptr, py::arg("parameter_overrides") = std::map<std::string, float_t>())
      .def("dump", &Instantiation::dump)
      .def_property("release_sites", &Instantiation::get_release_sites, &Instantiation::set_release_sites)
      .def_property("geometry_objects", &Instantiation::get_geometry_objects, &Instantiation::set_geometry_objects)
      .def_property("checkpointed_molecules", &Instantiation::get_checkpointed_molecules, &Instantiation::set_checkpointed_molecules)
    ;
}

std::string GenInstantiation::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  std::string exported_name = "instantiation";

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.Instantiation(" << nl;
  if (release_sites != std::vector<std::shared_ptr<ReleaseSite>>() && !skip_vectors_export()) {
    ss << ind << "release_sites = " << export_vec_release_sites(out, ctx, exported_name) << "," << nl;
  }
  if (geometry_objects != std::vector<std::shared_ptr<GeometryObject>>() && !skip_vectors_export()) {
    ss << ind << "geometry_objects = " << export_vec_geometry_objects(out, ctx, exported_name) << "," << nl;
  }
  if (checkpointed_molecules != std::vector<std::shared_ptr<BaseChkptMol>>() && !skip_vectors_export()) {
    ss << ind << "checkpointed_molecules = " << export_vec_checkpointed_molecules(out, ctx, exported_name) << "," << nl;
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

std::string GenInstantiation::export_vec_release_sites(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_release_sites";
  }
  else {
    exported_name = "release_sites";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < release_sites.size(); i++) {
    const auto& item = release_sites[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenInstantiation::export_vec_geometry_objects(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_geometry_objects";
  }
  else {
    exported_name = "geometry_objects";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < geometry_objects.size(); i++) {
    const auto& item = geometry_objects[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenInstantiation::export_vec_checkpointed_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_checkpointed_molecules";
  }
  else {
    exported_name = "checkpointed_molecules";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < checkpointed_molecules.size(); i++) {
    const auto& item = checkpointed_molecules[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

} // namespace API
} // namespace MCell

