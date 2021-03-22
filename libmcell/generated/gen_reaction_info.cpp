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
  return py::class_<ReactionInfo, std::shared_ptr<ReactionInfo>>(m, "ReactionInfo", "Data structure passed to a reaction callback.")
      .def(
          py::init<
          >()
      )
      .def("__str__", &ReactionInfo::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ReactionInfo::__eq__, py::arg("other"))
      .def("dump", &ReactionInfo::dump)
      .def_property("type", &ReactionInfo::get_type, &ReactionInfo::set_type, "Specifies whether the reaction is unimolecular or bimolecular and\nalso provides information in reactant types. \n")
      .def_property("reactant_ids", &ReactionInfo::get_reactant_ids, &ReactionInfo::set_reactant_ids, "IDs of the reacting molecules, contains 1 ID for a unimolecular or a molecule+surface class reaction , \n2 IDs for a bimolecular reaction.\nFor a bimolecular reaction, the first ID is always the molecule that was diffused and the second one \nis the molecule that was hit.\nIDs can be used to obtain location of the molecules. The position of the first molecule obtained through \nmodel.get_molecule() is the position of the diffusing molecule before the collision.\nAll the reactants are removed after return from this callback, unless they are kept by the reaction such as in A + B -> A + C.  \n")
      .def_property("product_ids", &ReactionInfo::get_product_ids, &ReactionInfo::set_product_ids, "IDs of reaction product molecules. They already exist in the simulated system together with reactants, however reactants \nwill be removed after return from this callback. \n")
      .def_property("reaction_rule", &ReactionInfo::get_reaction_rule, &ReactionInfo::set_reaction_rule, "Reaction rule of the reaction.")
      .def_property("time", &ReactionInfo::get_time, &ReactionInfo::set_time, "Time of the reaction")
      .def_property("pos3d", &ReactionInfo::get_pos3d, &ReactionInfo::set_pos3d, "Specifies where reaction occured in the 3d space, specific meaning depends on the reaction type:\n- unimolecular reaction - position of the reacting molecule,\n- volume-volume or surface-surface reaction - position of the first reactant,\n- volume-surface reaction - position where the volume molecule hit the wall with the surface molecule.\n")
      .def_property("geometry_object", &ReactionInfo::get_geometry_object, &ReactionInfo::set_geometry_object, "Set only for surface reactions or reactions with surface classes.\nObject on whose surface where the reaction occured.\n")
      .def_property("wall_index", &ReactionInfo::get_wall_index, &ReactionInfo::set_wall_index, "Set only for surface reactions or reactions with surface classes.\nIndex of wall belonging to the geometry_object where the reaction occured, \ni.e. where the volume molecule hit the wall with a surface molecule or\nwall where the diffusing surface reactant reacted.\n")
      .def_property("pos2d", &ReactionInfo::get_pos2d, &ReactionInfo::set_pos2d, "Set only for surface reactions or reactions with surface classes.\nSpecifies where reaction occured in the 2d UV coordinates defined by the wall where the reaction occured, \nspecific meaning depends on the reaction type:\n- unimolecular reaction - position of the reacting molecule,\n- volume-surface and surface-surface reaction - position of the second reactant.\n  \n  ")
    ;
}

} // namespace API
} // namespace MCell

