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
#include "gen_single_molecule_release_info.h"
#include "../api/single_molecule_release_info.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenSingleMoleculeReleaseInfo::check_semantics() const {
  if (!is_set(species)) {
    throw ValueError("Parameter 'species' must be set.");
  }
  if (!is_set(location)) {
    throw ValueError("Parameter 'location' must be set.");
  }
}

bool GenSingleMoleculeReleaseInfo::__eq__(const GenSingleMoleculeReleaseInfo& other) const {
  return
    name == other.name &&
    (
      (species != nullptr) ?
        ( (other.species != nullptr) ?
          (species->__eq__(*other.species)) : 
          false
        ) :
        ( (other.species != nullptr) ?
          false :
          true
        )
     )  &&
    location == other.location &&
    orientation == other.orientation;
}

void GenSingleMoleculeReleaseInfo::set_initialized() {
  if (is_set(species)) {
    species->set_initialized();
  }
  initialized = true;
}

std::string GenSingleMoleculeReleaseInfo::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "location=" << location << ", " <<
      "orientation=" << orientation;
  return ss.str();
}

py::class_<SingleMoleculeReleaseInfo> define_pybinding_SingleMoleculeReleaseInfo(py::module& m) {
  return py::class_<SingleMoleculeReleaseInfo, std::shared_ptr<SingleMoleculeReleaseInfo>>(m, "SingleMoleculeReleaseInfo")
      .def(
          py::init<
            std::shared_ptr<Species>,
            const Vec3&,
            const Orientation
          >(),
          py::arg("species"),
          py::arg("location"),
          py::arg("orientation") = Orientation::NONE
      )
      .def("check_semantics", &SingleMoleculeReleaseInfo::check_semantics)
      .def("__str__", &SingleMoleculeReleaseInfo::to_str, py::arg("ind") = std::string(""))
      .def("dump", &SingleMoleculeReleaseInfo::dump)
      .def_property("species", &SingleMoleculeReleaseInfo::get_species, &SingleMoleculeReleaseInfo::set_species)
      .def_property("location", &SingleMoleculeReleaseInfo::get_location, &SingleMoleculeReleaseInfo::set_location)
      .def_property("orientation", &SingleMoleculeReleaseInfo::get_orientation, &SingleMoleculeReleaseInfo::set_orientation)
    ;
}

} // namespace API
} // namespace MCell

