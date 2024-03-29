/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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

std::shared_ptr<ReactionInfo> GenReactionInfo::copy_reaction_info() const {
  std::shared_ptr<ReactionInfo> res = std::make_shared<ReactionInfo>(DefaultCtorArgType());
  res->type = type;
  res->reactant_ids = reactant_ids;
  res->product_ids = product_ids;
  res->reaction_rule = reaction_rule;
  res->time = time;
  res->pos3d = pos3d;
  res->geometry_object = geometry_object;
  res->wall_index = wall_index;
  res->pos2d = pos2d;

  return res;
}

std::shared_ptr<ReactionInfo> GenReactionInfo::deepcopy_reaction_info(py::dict) const {
  std::shared_ptr<ReactionInfo> res = std::make_shared<ReactionInfo>(DefaultCtorArgType());
  res->type = type;
  res->reactant_ids = reactant_ids;
  res->product_ids = product_ids;
  res->reaction_rule = is_set(reaction_rule) ? reaction_rule->deepcopy_reaction_rule() : nullptr;
  res->time = time;
  res->pos3d = pos3d;
  res->geometry_object = is_set(geometry_object) ? geometry_object->deepcopy_geometry_object() : nullptr;
  res->wall_index = wall_index;
  res->pos2d = pos2d;

  return res;
}

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
    true /*pos3d*/ &&
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
    true /*pos2d*/;
}

std::string GenReactionInfo::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << "ReactionInfo" << ": " <<
      "type=" << type << ", " <<
      "reactant_ids=" << vec_nonptr_to_str(reactant_ids, all_details, ind + "  ") << ", " <<
      "product_ids=" << vec_nonptr_to_str(product_ids, all_details, ind + "  ") << ", " <<
      "\n" << ind + "  " << "reaction_rule=" << "(" << ((reaction_rule != nullptr) ? reaction_rule->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "time=" << time << ", " <<
      "pos3d=" << vec_nonptr_to_str(pos3d, all_details, ind + "  ") << ", " <<
      "\n" << ind + "  " << "geometry_object=" << "(" << ((geometry_object != nullptr) ? geometry_object->to_str(all_details, ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "wall_index=" << wall_index << ", " <<
      "pos2d=" << vec_nonptr_to_str(pos2d, all_details, ind + "  ");
  return ss.str();
}

py::class_<ReactionInfo> define_pybinding_ReactionInfo(py::module& m) {
  return py::class_<ReactionInfo, std::shared_ptr<ReactionInfo>>(m, "ReactionInfo", "Data structure passed to a reaction callback registered with \nModel.register_reaction_callback.\n")
      .def(
          py::init<
          >()
      )
      .def("__copy__", &ReactionInfo::copy_reaction_info)
      .def("__deepcopy__", &ReactionInfo::deepcopy_reaction_info, py::arg("memo"))
      .def("__str__", &ReactionInfo::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &ReactionInfo::__eq__, py::arg("other"))
      .def("dump", &ReactionInfo::dump)
      .def_property("type", &ReactionInfo::get_type, &ReactionInfo::set_type, "Specifies whether the reaction is unimolecular or bimolecular and\nalso provides information on reactant types. \n")
      .def_property("reactant_ids", &ReactionInfo::get_reactant_ids, &ReactionInfo::set_reactant_ids, py::return_value_policy::reference, "IDs of the reacting molecules, contains 1 ID for a unimolecular or a molecule+surface class reaction, \n2 IDs for a bimolecular reaction.\nFor a bimolecular reaction, the first ID is always the molecule that diffused and the second one \nis the molecule that was hit.\nIDs can be used to obtain the location of the molecules. The position of the first molecule obtained through \nmodel.get_molecule() is the position of the diffusing molecule before the collision.\nAll the reactants are removed after return from this callback, unless they are kept by the reaction such as A in A + B -> A + C.  \n")
      .def_property("product_ids", &ReactionInfo::get_product_ids, &ReactionInfo::set_product_ids, py::return_value_policy::reference, "IDs of reaction product molecules. They already exist in the simulated system together with reactants; however reactants \nwill be removed after return from this callback. \n")
      .def_property("reaction_rule", &ReactionInfo::get_reaction_rule, &ReactionInfo::set_reaction_rule, "Reaction rule of the reaction that occured.")
      .def_property("time", &ReactionInfo::get_time, &ReactionInfo::set_time, "Time of the reaction.")
      .def_property("pos3d", &ReactionInfo::get_pos3d, &ReactionInfo::set_pos3d, py::return_value_policy::reference, "Specifies where reaction occurred in the 3d space, the specific meaning depends on the reaction type:\n- unimolecular reaction - position of the reacting molecule,\n- volume-volume or surface-surface reaction - position of the first reactant,\n- volume-surface reaction - position where the volume molecule hit the wall with the surface molecule.\n")
      .def_property("geometry_object", &ReactionInfo::get_geometry_object, &ReactionInfo::set_geometry_object, "The object on whose surface where the reaction occurred.\nSet only for surface reactions or reactions with surface classes.\n")
      .def_property("wall_index", &ReactionInfo::get_wall_index, &ReactionInfo::set_wall_index, "Set only for surface reactions or reactions with surface classes.\nIndex of wall belonging to the geometry_object where the reaction occured, \ni.e. wall where a volume molecule hit the surface molecule or\nwall where the diffusing surface reactant reacted.\n")
      .def_property("pos2d", &ReactionInfo::get_pos2d, &ReactionInfo::set_pos2d, py::return_value_policy::reference, "Set only for surface reactions or reactions with surface classes.\nSpecifies where reaction occurred in the 2d UV coordinates defined by the wall where the reaction occured, \nthe rspecific meaning depends on the reaction type:\n- unimolecular reaction - position of the reacting molecule,\n- volume-surface and surface-surface reaction - position of the second reactant.\n  \n  ")
    ;
}

} // namespace API
} // namespace MCell

