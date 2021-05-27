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
#include "gen_count_term.h"
#include "api/count_term.h"
#include "api/complex.h"
#include "api/count_term.h"
#include "api/reaction_rule.h"
#include "api/region.h"

namespace MCell {
namespace API {

void GenCountTerm::check_semantics() const {
}

void GenCountTerm::set_initialized() {
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

void GenCountTerm::set_all_attributes_as_default_or_unset() {
  class_name = "CountTerm";
  species_pattern = nullptr;
  molecules_pattern = nullptr;
  reaction_rule = nullptr;
  region = nullptr;
  node_type = ExprNodeType::LEAF;
  left_node = nullptr;
  right_node = nullptr;
  initial_reactions_count = 0;
}

CountTerm GenCountTerm::copy_count_term() const {
  if (initialized) {
    throw RuntimeError("Object of class CountTerm cannot be cloned with 'copy' after this object was used in model initialization.");
  }
  CountTerm res = CountTerm(DefaultCtorArgType());
  res.class_name = class_name;
  res.species_pattern = species_pattern;
  res.molecules_pattern = molecules_pattern;
  res.reaction_rule = reaction_rule;
  res.region = region;
  res.node_type = node_type;
  res.left_node = left_node;
  res.right_node = right_node;
  res.initial_reactions_count = initial_reactions_count;

  return res;
}

bool GenCountTerm::__eq__(const CountTerm& other) const {
  return
    (
      (is_set(species_pattern)) ?
        (is_set(other.species_pattern) ?
          (species_pattern->__eq__(*other.species_pattern)) : 
          false
        ) :
        (is_set(other.species_pattern) ?
          false :
          true
        )
     )  &&
    (
      (is_set(molecules_pattern)) ?
        (is_set(other.molecules_pattern) ?
          (molecules_pattern->__eq__(*other.molecules_pattern)) : 
          false
        ) :
        (is_set(other.molecules_pattern) ?
          false :
          true
        )
     )  &&
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
    (
      (is_set(region)) ?
        (is_set(other.region) ?
          (region->__eq__(*other.region)) : 
          false
        ) :
        (is_set(other.region) ?
          false :
          true
        )
     )  &&
    node_type == other.node_type &&
    (
      (is_set(left_node)) ?
        (is_set(other.left_node) ?
          (left_node->__eq__(*other.left_node)) : 
          false
        ) :
        (is_set(other.left_node) ?
          false :
          true
        )
     )  &&
    (
      (is_set(right_node)) ?
        (is_set(other.right_node) ?
          (right_node->__eq__(*other.right_node)) : 
          false
        ) :
        (is_set(other.right_node) ?
          false :
          true
        )
     )  &&
    initial_reactions_count == other.initial_reactions_count;
}

bool GenCountTerm::eq_nonarray_attributes(const CountTerm& other, const bool ignore_name) const {
  return
    (
      (is_set(species_pattern)) ?
        (is_set(other.species_pattern) ?
          (species_pattern->__eq__(*other.species_pattern)) : 
          false
        ) :
        (is_set(other.species_pattern) ?
          false :
          true
        )
     )  &&
    (
      (is_set(molecules_pattern)) ?
        (is_set(other.molecules_pattern) ?
          (molecules_pattern->__eq__(*other.molecules_pattern)) : 
          false
        ) :
        (is_set(other.molecules_pattern) ?
          false :
          true
        )
     )  &&
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
    (
      (is_set(region)) ?
        (is_set(other.region) ?
          (region->__eq__(*other.region)) : 
          false
        ) :
        (is_set(other.region) ?
          false :
          true
        )
     )  &&
    node_type == other.node_type &&
    (
      (is_set(left_node)) ?
        (is_set(other.left_node) ?
          (left_node->__eq__(*other.left_node)) : 
          false
        ) :
        (is_set(other.left_node) ?
          false :
          true
        )
     )  &&
    (
      (is_set(right_node)) ?
        (is_set(other.right_node) ?
          (right_node->__eq__(*other.right_node)) : 
          false
        ) :
        (is_set(other.right_node) ?
          false :
          true
        )
     )  &&
    initial_reactions_count == other.initial_reactions_count;
}

std::string GenCountTerm::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "\n" << ind + "  " << "species_pattern=" << "(" << ((species_pattern != nullptr) ? species_pattern->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "molecules_pattern=" << "(" << ((molecules_pattern != nullptr) ? molecules_pattern->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "reaction_rule=" << "(" << ((reaction_rule != nullptr) ? reaction_rule->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "region=" << "(" << ((region != nullptr) ? region->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "node_type=" << node_type << ", " <<
      "\n" << ind + "  " << "left_node=" << "(" << ((left_node != nullptr) ? left_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "initial_reactions_count=" << initial_reactions_count;
  return ss.str();
}

py::class_<CountTerm> define_pybinding_CountTerm(py::module& m) {
  return py::class_<CountTerm, std::shared_ptr<CountTerm>>(m, "CountTerm", "A count observable can be defined as an expression composed of addition\nor subtraction individual count terms. This class represents one count term\nin this expression.\n \n")
      .def(
          py::init<
            std::shared_ptr<Complex>,
            std::shared_ptr<Complex>,
            std::shared_ptr<ReactionRule>,
            std::shared_ptr<Region>,
            const ExprNodeType,
            std::shared_ptr<CountTerm>,
            std::shared_ptr<CountTerm>,
            const uint64_t
          >(),
          py::arg("species_pattern") = nullptr,
          py::arg("molecules_pattern") = nullptr,
          py::arg("reaction_rule") = nullptr,
          py::arg("region") = nullptr,
          py::arg("node_type") = ExprNodeType::LEAF,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr,
          py::arg("initial_reactions_count") = 0
      )
      .def("check_semantics", &CountTerm::check_semantics)
      .def("__copy__", &CountTerm::copy_count_term)
      .def("__str__", &CountTerm::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &CountTerm::__eq__, py::arg("other"))
      .def("__add__", &CountTerm::__add__, py::arg("op2"), "Create a new CountTerm that represents addition of two count terms.\nUsually used through operator '+' such as in ct1 + ct2.  \n\n- op2\n")
      .def("__sub__", &CountTerm::__sub__, py::arg("op2"), "Create a new CountTerm that represents subtraction of two count terms.\nUsually used through operator '-' such as in ct1 - ct2.  \n\n- op2\n")
      .def("dump", &CountTerm::dump)
      .def_property("species_pattern", &CountTerm::get_species_pattern, &CountTerm::set_species_pattern, "Count the number of molecules that match the given complex instance pattern.\nThis corresponds to the BNGL 'Species' specifier in the BNGL seed species section.\nCounts each molecule exactly once. \nIf the pattern has a compartment set, this specifies the counted region.\nExactly one of species_pattern, molecules_pattern, and reaction_rule must be set. \n")
      .def_property("molecules_pattern", &CountTerm::get_molecules_pattern, &CountTerm::set_molecules_pattern, "Count the number of matches of the given pattern on molecules.\nThis corresponds to the BNGL 'Molecules' specifier in the BNGL seed species section.\nThe observable will increment the count every time the pattern matches the molecule.\nFor instance, pattern A will match a complex A(a!1).B(a!1,a!2).A(b!2) twice. \nWhen the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a \nmolecule A(b.a!1).A(a!1,b!2).B(a!2) will be counted twice because the \npattern may match in two different ways. \nIf the pattern has a compartment set, the compartment is used to filter out the molecules.   \nExactly one of species_pattern, molecules_pattern, and reaction_rule must be set.\n")
      .def_property("reaction_rule", &CountTerm::get_reaction_rule, &CountTerm::set_reaction_rule, "Count the number of applications of this specific reactions that occurred since the\nstart of the simulation.\nExactly one of species_pattern, molecules_pattern, and reaction_rule must be set.\n  \n")
      .def_property("region", &CountTerm::get_region, &CountTerm::set_region, "Only a GeometryObject or SurfaceRegion can be passed as the region argument, \ncompound regions (created with +, -, *) are not supproted yet.   \nCan be combined with a compartment specified in the species_pattern or molecules_pattern.\nIf compartment in species_pattern or molecules_pattern is not specified and \nregion is left unset, counting is done in the whole world.\n")
      .def_property("node_type", &CountTerm::get_node_type, &CountTerm::set_node_type, "Internal, used to specify what type of count expression node this object represents.")
      .def_property("left_node", &CountTerm::get_left_node, &CountTerm::set_left_node, "Internal, when node_type is not Leaf, this is the left operand.")
      .def_property("right_node", &CountTerm::get_right_node, &CountTerm::set_right_node, "Internal, when node_type is not Leaf, this is the right operand.")
      .def_property("initial_reactions_count", &CountTerm::get_initial_reactions_count, &CountTerm::set_initial_reactions_count, "Used for checkpointing, allows to set initial count of reactions that occurred.\nIgnored when molecules are counted.\n")
    ;
}

std::string GenCountTerm::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "count_term_" + std::to_string(ctx.postinc_counter("count_term"));
  if (!export_even_if_already_exported()) {
    ctx.add_exported(this, exported_name);
  }

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.CountTerm(" << nl;
  if (is_set(species_pattern)) {
    ss << ind << "species_pattern = " << species_pattern->export_to_python(out, ctx) << "," << nl;
  }
  if (is_set(molecules_pattern)) {
    ss << ind << "molecules_pattern = " << molecules_pattern->export_to_python(out, ctx) << "," << nl;
  }
  if (is_set(reaction_rule)) {
    ss << ind << "reaction_rule = " << reaction_rule->export_to_python(out, ctx) << "," << nl;
  }
  if (is_set(region)) {
    ss << ind << "region = " << region->export_to_python(out, ctx) << "," << nl;
  }
  if (node_type != ExprNodeType::LEAF) {
    ss << ind << "node_type = " << node_type << "," << nl;
  }
  if (is_set(left_node)) {
    ss << ind << "left_node = " << left_node->export_to_python(out, ctx) << "," << nl;
  }
  if (is_set(right_node)) {
    ss << ind << "right_node = " << right_node->export_to_python(out, ctx) << "," << nl;
  }
  if (initial_reactions_count != 0) {
    ss << ind << "initial_reactions_count = " << initial_reactions_count << "," << nl;
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

} // namespace API
} // namespace MCell

