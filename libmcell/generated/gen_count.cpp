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
#include "gen_count.h"
#include "../api/count.h"
#include "../api/complex.h"
#include "../api/count_term.h"
#include "../api/reaction_rule.h"
#include "../api/region.h"

namespace MCell {
namespace API {

void GenCount::check_semantics() const {
  if (!is_set(file_name)) {
    throw ValueError("Parameter 'file_name' must be set.");
  }
}

bool GenCount::__eq__(const GenCount& other) const {
  return
    name == other.name &&
    file_name == other.file_name &&
    (
      (count_expression != nullptr) ?
        ( (other.count_expression != nullptr) ?
          (count_expression->__eq__(*other.count_expression)) : 
          false
        ) :
        ( (other.count_expression != nullptr) ?
          false :
          true
        )
     )  &&
    multiplier == other.multiplier &&
    every_n_timesteps == other.every_n_timesteps &&
    (
      (species_pattern != nullptr) ?
        ( (other.species_pattern != nullptr) ?
          (species_pattern->__eq__(*other.species_pattern)) : 
          false
        ) :
        ( (other.species_pattern != nullptr) ?
          false :
          true
        )
     )  &&
    (
      (molecules_pattern != nullptr) ?
        ( (other.molecules_pattern != nullptr) ?
          (molecules_pattern->__eq__(*other.molecules_pattern)) : 
          false
        ) :
        ( (other.molecules_pattern != nullptr) ?
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

void GenCount::set_initialized() {
  if (is_set(count_expression)) {
    count_expression->set_initialized();
  }
  if (is_set(species_pattern)) {
    species_pattern->set_initialized();
  }
  if (is_set(molecules_pattern)) {
    molecules_pattern->set_initialized();
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

void GenCount::set_all_attributes_as_default_or_unset() {
  class_name = "Count";
  file_name = STR_UNSET;
  count_expression = nullptr;
  multiplier = 1;
  every_n_timesteps = 1;
  species_pattern = nullptr;
  molecules_pattern = nullptr;
  reaction_rule = nullptr;
  region = nullptr;
  node_type = ExprNodeType::LEAF;
  left_node = nullptr;
  right_node = nullptr;
}

std::string GenCount::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "file_name=" << file_name << ", " <<
      "\n" << ind + "  " << "count_expression=" << "(" << ((count_expression != nullptr) ? count_expression->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "multiplier=" << multiplier << ", " <<
      "every_n_timesteps=" << every_n_timesteps << ", " <<
      "\n" << ind + "  " << "species_pattern=" << "(" << ((species_pattern != nullptr) ? species_pattern->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "molecules_pattern=" << "(" << ((molecules_pattern != nullptr) ? molecules_pattern->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "reaction_rule=" << "(" << ((reaction_rule != nullptr) ? reaction_rule->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "region=" << "(" << ((region != nullptr) ? region->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")";
  return ss.str();
}

py::class_<Count> define_pybinding_Count(py::module& m) {
  return py::class_<Count, CountTerm, std::shared_ptr<Count>>(m, "Count")
      .def(
          py::init<
            const std::string&,
            std::shared_ptr<CountTerm>,
            const float_t,
            const float_t,
            std::shared_ptr<Complex>,
            std::shared_ptr<Complex>,
            std::shared_ptr<ReactionRule>,
            std::shared_ptr<Region>,
            const ExprNodeType,
            std::shared_ptr<CountTerm>,
            std::shared_ptr<CountTerm>
          >(),
          py::arg("file_name"),
          py::arg("count_expression") = nullptr,
          py::arg("multiplier") = 1,
          py::arg("every_n_timesteps") = 1,
          py::arg("species_pattern") = nullptr,
          py::arg("molecules_pattern") = nullptr,
          py::arg("reaction_rule") = nullptr,
          py::arg("region") = nullptr,
          py::arg("node_type") = ExprNodeType::LEAF,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr
      )
      .def("check_semantics", &Count::check_semantics)
      .def("__str__", &Count::to_str, py::arg("ind") = std::string(""))
      .def("__repr__", &Count::to_str, py::arg("ind") = std::string(""))
      .def("dump", &Count::dump)
      .def_property("file_name", &Count::get_file_name, &Count::set_file_name)
      .def_property("count_expression", &Count::get_count_expression, &Count::set_count_expression)
      .def_property("multiplier", &Count::get_multiplier, &Count::set_multiplier)
      .def_property("every_n_timesteps", &Count::get_every_n_timesteps, &Count::set_every_n_timesteps)
    ;
}

} // namespace API
} // namespace MCell

