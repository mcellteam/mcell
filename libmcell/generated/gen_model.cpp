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
#include <pybind11/stl.h>
#include "gen_model.h"
#include "../api/model.h"
#include "../api/config.h"
#include "../api/count.h"
#include "../api/geometry_object.h"
#include "../api/instantiation_data.h"
#include "../api/notifications.h"
#include "../api/observables.h"
#include "../api/reaction_rule.h"
#include "../api/release_site.h"
#include "../api/species.h"
#include "../api/subsystem.h"
#include "../api/surface_class.h"
#include "../api/viz_output.h"
#include "../api/warnings.h"

namespace MCell {
namespace API {

std::string GenModel::to_str(const std::string ind) const {
  #if 0 // not generated correctly yet
  std::stringstream ss;
  ss << "Model" << ": " <<
      "\n" << ind + "  " << "config=" << "(" << ((config != nullptr) ? config->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "warnings=" << "(" << ((warnings != nullptr) ? warnings->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "notifications=" << "(" << ((notifications != nullptr) ? notifications->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "species=" << vec_ptr_to_str(species, ind + "  ") << ", " << "\n" << ind + "  " <<
      "reaction_rules=" << vec_ptr_to_str(reaction_rules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_classes=" << vec_ptr_to_str(surface_classes, ind + "  ") << ", " << "\n" << ind + "  " <<
      "release_sites=" << vec_ptr_to_str(release_sites, ind + "  ") << ", " << "\n" << ind + "  " <<
      "geometry_objects=" << vec_ptr_to_str(geometry_objects, ind + "  ") << ", " << "\n" << ind + "  " <<
      "viz_outputs=" << vec_ptr_to_str(viz_outputs, ind + "  ") << ", " << "\n" << ind + "  " <<
      "counts=" << vec_ptr_to_str(counts, ind + "  ");
  return ss.str();
  #else
  return "";
  #endif
}

py::class_<Model> define_pybinding_Model(py::module& m) {
  return py::class_<Model, std::shared_ptr<Model>>(m, "Model")
      .def(
          py::init<
          >()
      )
      .def("__str__", &Model::to_str, py::arg("ind") = std::string(""))
      .def("initialize", &Model::initialize)
      .def("run_iterations", &Model::run_iterations, py::arg("iterations"))
      .def("end_simulation", &Model::end_simulation, py::arg("print_final_report") = true)
      .def("add_subsystem", &Model::add_subsystem, py::arg("subsystem"))
      .def("add_instantiation_data", &Model::add_instantiation_data, py::arg("instantiation_data"))
      .def("add_observables", &Model::add_observables, py::arg("observables"))
      .def("dump_internal_state", &Model::dump_internal_state)
      .def("export_data_model", &Model::export_data_model, py::arg("file") = STR_UNSET)
      .def("add_species", &Model::add_species, py::arg("s"))
      .def("find_species", &Model::find_species, py::arg("name"))
      .def("add_reaction_rule", &Model::add_reaction_rule, py::arg("r"))
      .def("find_reaction_rule", &Model::find_reaction_rule, py::arg("name"))
      .def("add_surface_class", &Model::add_surface_class, py::arg("sc"))
      .def("find_surface_class", &Model::find_surface_class, py::arg("name"))
      .def("add_release_site", &Model::add_release_site, py::arg("s"))
      .def("find_release_site", &Model::find_release_site, py::arg("name"))
      .def("add_geometry_object", &Model::add_geometry_object, py::arg("o"))
      .def("find_geometry_object", &Model::find_geometry_object, py::arg("name"))
      .def("add_viz_output", &Model::add_viz_output, py::arg("viz_output"))
      .def("add_count", &Model::add_count, py::arg("count"))
      .def("dump", &Model::dump)
      .def_property("config", &Model::get_config, &Model::set_config)
      .def_property("warnings", &Model::get_warnings, &Model::set_warnings)
      .def_property("notifications", &Model::get_notifications, &Model::set_notifications)
      .def_property("species", &Model::get_species, &Model::set_species)
      .def_property("reaction_rules", &Model::get_reaction_rules, &Model::set_reaction_rules)
      .def_property("surface_classes", &Model::get_surface_classes, &Model::set_surface_classes)
      .def_property("release_sites", &Model::get_release_sites, &Model::set_release_sites)
      .def_property("geometry_objects", &Model::get_geometry_objects, &Model::set_geometry_objects)
      .def_property("viz_outputs", &Model::get_viz_outputs, &Model::set_viz_outputs)
      .def_property("counts", &Model::get_counts, &Model::set_counts)
    ;
}

} // namespace API
} // namespace MCell

