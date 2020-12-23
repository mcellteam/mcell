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
#include "api/python_export_utils.h"
#include "gen_reaction_info.h"
#include "api/reaction_info.h"
#include "api/geometry_object.h"
#include "api/reaction_rule.h"

namespace MCell {
namespace API {

bool GenReactionInfo::__eq__(const ReactionInfo& other) const {
  return
    type == other.type &&
    reactant_ids == other.reactant_ids &&
    product_ids == other.product_ids &&
    (
      (is_set(reaction_rule)) ?
        (is_set(other.reaction_rule) ?
          (reaction_rule->__eq__(*other.reaction_rule)) : 
          false
        ) :
        (is_set(other.reaction_rule) ?
          false :
          true
        )
     )  &&
    time == other.time &&
    pos3d == other.pos3d &&
    (
      (is_set(geometry_object)) ?
        (is_set(other.geometry_object) ?
          (geometry_object->__eq__(*other.geometry_object)) : 
          false
        ) :
        (is_set(other.geometry_object) ?
          false :
          true
        )
     )  &&
    wall_index == other.wall_index &&
    pos2d == other.pos2d;
}

bool GenReactionInfo::eq_nonarray_attributes(const ReactionInfo& other, const bool ignore_name) const {
  return
    type == other.type &&
    true /*reactant_ids*/ &&
    true /*product_ids*/ &&
    (
      (is_set(reaction_rule)) ?
        (is_set(other.reaction_rule) ?
          (reaction_rule->__eq__(*other.reaction_rule)) : 
          false
        ) :
        (is_set(other.reaction_rule) ?
          false :
          true
        )
     )  &&
    time == other.time &&
    pos3d == other.pos3d &&
    (
      (is_set(geometry_object)) ?
        (is_set(other.geometry_object) ?
          (geometry_object->__eq__(*other.geometry_object)) : 
          false
        ) :
        (is_set(other.geometry_object) ?
          false :
          true
        )
     )  &&
    wall_index == other.wall_index &&
    pos2d == other.pos2d;
}

std::string GenReactionInfo::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "ReactionInfo" << ": " <<
      "type=" << type << ", " <<
      "reactant_ids=" << vec_nonptr_to_str(reactant_ids, ind + "  ") << ", " <<
      "product_ids=" << vec_nonptr_to_str(product_ids, ind + "  ") << ", " <<
      "\n" << ind + "  " << "reaction_rule=" << "(" << ((reaction_rule != nullptr) ? reaction_rule->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "time=" << time << ", " <<
      "pos3d=" << pos3d << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index << ", " <<
      "pos2d=" << pos2d;
  return ss.str();
}

py::class_<ReactionInfo> define_pybinding_ReactionInfo(py::module& m) {
  return py::class_<ReactionInfo, std::shared_ptr<ReactionInfo>>(m, "ReactionInfo")
      .def(
          py::init<
          >()
      )
      .def("__str__", &ReactionInfo::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ReactionInfo::__eq__, py::arg("other"))
      .def("dump", &ReactionInfo::dump)
      .def_property("type", &ReactionInfo::get_type, &ReactionInfo::set_type)
      .def_property("reactant_ids", &ReactionInfo::get_reactant_ids, &ReactionInfo::set_reactant_ids)
      .def_property("product_ids", &ReactionInfo::get_product_ids, &ReactionInfo::set_product_ids)
      .def_property("reaction_rule", &ReactionInfo::get_reaction_rule, &ReactionInfo::set_reaction_rule)
      .def_property("time", &ReactionInfo::get_time, &ReactionInfo::set_time)
      .def_property("pos3d", &ReactionInfo::get_pos3d, &ReactionInfo::set_pos3d)
      .def_property("geometry_object", &ReactionInfo::get_geometry_object, &ReactionInfo::set_geometry_object)
      .def_property("wall_index", &ReactionInfo::get_wall_index, &ReactionInfo::set_wall_index)
      .def_property("pos2d", &ReactionInfo::get_pos2d, &ReactionInfo::set_pos2d)
    ;
}

} // namespace API
} // namespace MCell

