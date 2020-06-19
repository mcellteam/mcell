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
#include "gen_count_term.h"
#include "../api/count_term.h"
#include "../api/count_term.h"
#include "../api/reaction_rule.h"
#include "../api/region.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenCountTerm::check_semantics() const {
}

bool GenCountTerm::__eq__(const GenCountTerm& other) const {
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
    (
      (reaction_rule != nullptr) ?
        ( (other.reaction_rule != nullptr) ?
          (reaction_rule->__eq__(*other.reaction_rule)) : 
          false
        ) :
        ( (other.reaction_rule != nullptr) ?
          false :
          true
        )
     )  &&
    (
      (region != nullptr) ?
        ( (other.region != nullptr) ?
          (region->__eq__(*other.region)) : 
          false
        ) :
        ( (other.region != nullptr) ?
          false :
          true
        )
     )  &&
    orientation == other.orientation &&
    node_type == other.node_type &&
    (
      (left_node != nullptr) ?
        ( (other.left_node != nullptr) ?
          (left_node->__eq__(*other.left_node)) : 
          false
        ) :
        ( (other.left_node != nullptr) ?
          false :
          true
        )
     )  &&
    (
      (right_node != nullptr) ?
        ( (other.right_node != nullptr) ?
          (right_node->__eq__(*other.right_node)) : 
          false
        ) :
        ( (other.right_node != nullptr) ?
          false :
          true
        )
     ) ;
}

void GenCountTerm::set_initialized() {
  if (is_set(species)) {
    species->set_initialized();
  }
  if (is_set(reaction_rule)) {
    reaction_rule->set_initialized();
  }
  if (is_set(region)) {
    region->set_initialized();
  }
  if (is_set(left_node)) {
    left_node->set_initialized();
  }
  if (is_set(right_node)) {
    right_node->set_initialized();
  }
  initialized = true;
}

void GenCountTerm::set_all_attributes_as_default_or_unset() {
  class_name = "CountTerm";
  species = nullptr;
  reaction_rule = nullptr;
  region = nullptr;
  orientation = Orientation::NOT_SET;
  node_type = ExprNodeType::LEAF;
  left_node = nullptr;
  right_node = nullptr;
}

std::string GenCountTerm::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "reaction_rule=" << "(" << ((reaction_rule != nullptr) ? reaction_rule->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "region=" << "(" << ((region != nullptr) ? region->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "orientation=" << orientation << ", " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<CountTerm> define_pybinding_CountTerm(py::module& m) {
  return py::class_<CountTerm, std::shared_ptr<CountTerm>>(m, "CountTerm")
      .def(
          py::init<
            std::shared_ptr<Species>,
            std::shared_ptr<ReactionRule>,
            std::shared_ptr<Region>,
            const Orientation,
            const ExprNodeType,
            std::shared_ptr<CountTerm>,
            std::shared_ptr<CountTerm>
          >(),
          py::arg("species") = nullptr,
          py::arg("reaction_rule") = nullptr,
          py::arg("region") = nullptr,
          py::arg("orientation") = Orientation::NOT_SET,
          py::arg("node_type") = ExprNodeType::LEAF,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &CountTerm::check_semantics)
      .def("__str__", &CountTerm::to_str, py::arg("ind") = std::string(""))
      .def("__add__", &CountTerm::__add__, py::arg("op2"))
      .def("__sub__", &CountTerm::__sub__, py::arg("op2"))
      .def("dump", &CountTerm::dump)
      .def_property("species", &CountTerm::get_species, &CountTerm::set_species)
      .def_property("reaction_rule", &CountTerm::get_reaction_rule, &CountTerm::set_reaction_rule)
      .def_property("region", &CountTerm::get_region, &CountTerm::set_region)
      .def_property("orientation", &CountTerm::get_orientation, &CountTerm::set_orientation)
      .def_property("node_type", &CountTerm::get_node_type, &CountTerm::set_node_type)
      .def_property("left_node", &CountTerm::get_left_node, &CountTerm::set_left_node)
      .def_property("right_node", &CountTerm::get_right_node, &CountTerm::set_right_node)
    ;
}

} // namespace API
} // namespace MCell

