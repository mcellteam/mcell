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
#include "gen_count_term.h"
#include "../api/count_term.h"
#include "../api/count_term.h"
#include "../api/geometry_object.h"
#include "../api/reaction_rule.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenCountTerm::check_semantics() const {
}

std::string GenCountTerm::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "reaction_rule=" << "(" << ((reaction_rule != nullptr) ? reaction_rule->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "region=" << "(" << ((region != nullptr) ? region->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

bool GenCountTerm::__eq__(const GenCountTerm& other) const {
  return
    name == other.name &&
    species->__eq__(*other.species) &&
    reaction_rule->__eq__(*other.reaction_rule) &&
    region->__eq__(*other.region) &&
    node_type == other.node_type &&
    left_node->__eq__(*other.left_node) &&
    right_node->__eq__(*other.right_node);
}

py::class_<CountTerm> define_pybinding_CountTerm(py::module& m) {
  return py::class_<CountTerm, std::shared_ptr<CountTerm>>(m, "CountTerm")
      .def(
          py::init<
            std::shared_ptr<Species>,
            std::shared_ptr<ReactionRule>,
            std::shared_ptr<GeometryObject>,
            const ExprNodeType,
            std::shared_ptr<CountTerm>,
            std::shared_ptr<CountTerm>
          >(),
          py::arg("species") = nullptr,
          py::arg("reaction_rule") = nullptr,
          py::arg("region") = nullptr,
          py::arg("node_type") = ExprNodeType::Leaf,
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
      .def_property("node_type", &CountTerm::get_node_type, &CountTerm::set_node_type)
      .def_property("left_node", &CountTerm::get_left_node, &CountTerm::set_left_node)
      .def_property("right_node", &CountTerm::get_right_node, &CountTerm::set_right_node)
    ;
}

} // namespace API
} // namespace MCell

