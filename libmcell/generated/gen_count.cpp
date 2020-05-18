/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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
#include "gen_count.h"
#include "../api/count.h"
#include "../api/count_term.h"
#include "../api/geometry_object.h"
#include "../api/reaction_rule.h"
#include "../api/species.h"

namespace MCell {
namespace API {

void GenCount::check_semantics() const {
  if (!is_set(filename)) {
    throw ValueError("Parameter 'filename' must be set.");
  }
}

std::string GenCount::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "filename=" << filename << ", " <<
      "\n" << ind + "  " << "count_expression=" << "(" << ((count_expression != nullptr) ? count_expression->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "every_n_timesteps=" << every_n_timesteps << ", " <<
      "\n" << ind + "  " << "species=" << "(" << ((species != nullptr) ? species->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "reaction_rule=" << "(" << ((reaction_rule != nullptr) ? reaction_rule->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "region=" << "(" << ((region != nullptr) ? region->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<Count> define_pybinding_Count(py::module& m) {
  return py::class_<Count, std::shared_ptr<Count>>(m, "Count")
      .def(
          py::init<
            const std::string&,
            std::shared_ptr<CountTerm>,
            const int,
            std::shared_ptr<Species>,
            std::shared_ptr<ReactionRule>,
            std::shared_ptr<GeometryObject>,
            const ExprNodeType,
            std::shared_ptr<CountTerm>,
            std::shared_ptr<CountTerm>
          >(),
          py::arg("filename"),
          py::arg("count_expression") = nullptr,
          py::arg("every_n_timesteps") = 1,
          py::arg("species") = nullptr,
          py::arg("reaction_rule") = nullptr,
          py::arg("region") = nullptr,
          py::arg("node_type") = ExprNodeType::Leaf,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &Count::check_semantics)
      .def("__str__", &Count::to_str, py::arg("ind") = std::string(""))
      .def("__add__", &Count::__add__, py::arg("op2"))
      .def("__sub__", &Count::__sub__, py::arg("op2"))
      .def("dump", &Count::dump)
      .def_property("filename", &Count::get_filename, &Count::set_filename)
      .def_property("count_expression", &Count::get_count_expression, &Count::set_count_expression)
      .def_property("every_n_timesteps", &Count::get_every_n_timesteps, &Count::set_every_n_timesteps)
      .def_property("species", &Count::get_species, &Count::set_species)
      .def_property("reaction_rule", &Count::get_reaction_rule, &Count::set_reaction_rule)
      .def_property("region", &Count::get_region, &Count::set_region)
      .def_property("node_type", &Count::get_node_type, &Count::set_node_type)
      .def_property("left_node", &Count::get_left_node, &Count::set_left_node)
      .def_property("right_node", &Count::get_right_node, &Count::set_right_node)
    ;
}

} // namespace API
} // namespace MCell

