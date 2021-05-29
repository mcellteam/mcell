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
#include "gen_count.h"
#include "api/count.h"
#include "api/complex.h"
#include "api/count_term.h"
#include "api/reaction_rule.h"
#include "api/region.h"

namespace MCell {
namespace API {

void GenCount::check_semantics() const {
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
  name = STR_UNSET;
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
  initial_reactions_count = 0;
}

std::shared_ptr<Count> GenCount::copy_count() const {
  if (initialized) {
    throw RuntimeError("Object of class Count cannot be cloned with 'copy' after this object was used in model initialization.");
  }

  std::shared_ptr<Count> res = std::make_shared<Count>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->file_name = file_name;
  res->count_expression = count_expression;
  res->multiplier = multiplier;
  res->every_n_timesteps = every_n_timesteps;
  res->species_pattern = species_pattern;
  res->molecules_pattern = molecules_pattern;
  res->reaction_rule = reaction_rule;
  res->region = region;
  res->node_type = node_type;
  res->left_node = left_node;
  res->right_node = right_node;
  res->initial_reactions_count = initial_reactions_count;

  return res;
}

std::shared_ptr<Count> GenCount::deepcopy_count(py::dict) const {
  if (initialized) {
    throw RuntimeError("Object of class Count cannot be cloned with 'deepcopy' after this object was used in model initialization.");
  }

  std::shared_ptr<Count> res = std::make_shared<Count>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->file_name = file_name;
  res->count_expression = is_set(count_expression) ? count_expression->deepcopy_count_term() : nullptr;
  res->multiplier = multiplier;
  res->every_n_timesteps = every_n_timesteps;
  res->species_pattern = is_set(species_pattern) ? species_pattern->deepcopy_complex() : nullptr;
  res->molecules_pattern = is_set(molecules_pattern) ? molecules_pattern->deepcopy_complex() : nullptr;
  res->reaction_rule = is_set(reaction_rule) ? reaction_rule->deepcopy_reaction_rule() : nullptr;
  res->region = is_set(region) ? region->deepcopy_region() : nullptr;
  res->node_type = node_type;
  res->left_node = is_set(left_node) ? left_node->deepcopy_count_term() : nullptr;
  res->right_node = is_set(right_node) ? right_node->deepcopy_count_term() : nullptr;
  res->initial_reactions_count = initial_reactions_count;

  return res;
}

bool GenCount::__eq__(const Count& other) const {
  return
    name == other.name &&
    file_name == other.file_name &&
    (
      (is_set(count_expression)) ?
        (is_set(other.count_expression) ?
          (count_expression->__eq__(*other.count_expression)) : 
          false
        ) :
        (is_set(other.count_expression) ?
          false :
          true
        )
     )  &&
    multiplier == other.multiplier &&
    every_n_timesteps == other.every_n_timesteps &&
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

bool GenCount::eq_nonarray_attributes(const Count& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    file_name == other.file_name &&
    (
      (is_set(count_expression)) ?
        (is_set(other.count_expression) ?
          (count_expression->__eq__(*other.count_expression)) : 
          false
        ) :
        (is_set(other.count_expression) ?
          false :
          true
        )
     )  &&
    multiplier == other.multiplier &&
    every_n_timesteps == other.every_n_timesteps &&
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

std::string GenCount::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
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
      "right_node=" << "(" << ((right_node != nullptr) ? right_node->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "initial_reactions_count=" << initial_reactions_count;
  return ss.str();
}

py::class_<Count> define_pybinding_Count(py::module& m) {
  return py::class_<Count, CountTerm, std::shared_ptr<Count>>(m, "Count", "Represents a molecule or reaction count observable. \nInherits all members from class CountTerm. \nUsually just one observable is counted and a member of this class count_expression \nis not set. If an expression is needed, it is stored in the count_expression\nattribute.\n")
      .def(
          py::init<
            const std::string&,
            const std::string&,
            std::shared_ptr<CountTerm>,
            const double,
            const double,
            std::shared_ptr<Complex>,
            std::shared_ptr<Complex>,
            std::shared_ptr<ReactionRule>,
            std::shared_ptr<Region>,
            const ExprNodeType,
            std::shared_ptr<CountTerm>,
            std::shared_ptr<CountTerm>,
            const uint64_t
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("file_name") = STR_UNSET,
          py::arg("count_expression") = nullptr,
          py::arg("multiplier") = 1,
          py::arg("every_n_timesteps") = 1,
          py::arg("species_pattern") = nullptr,
          py::arg("molecules_pattern") = nullptr,
          py::arg("reaction_rule") = nullptr,
          py::arg("region") = nullptr,
          py::arg("node_type") = ExprNodeType::LEAF,
          py::arg("left_node") = nullptr,
          py::arg("right_node") = nullptr,
          py::arg("initial_reactions_count") = 0
      )
      .def("check_semantics", &Count::check_semantics)
      .def("__copy__", &Count::copy_count)
      .def("__deepcopy__", &Count::deepcopy_count, py::arg("memo"))
      .def("__str__", &Count::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Count::__eq__, py::arg("other"))
      .def("get_current_value", &Count::get_current_value, "Returns the current value for this count. Cannot be used to count reactions.\nThe model must be initialized with this Count present as one of the observables.\n")
      .def("dump", &Count::dump)
      .def_property("name", &Count::get_name, &Count::set_name, "Name of a count may be specified when one needs to search for them later. \nWhen the count is created when a BNGL file is loaded, its name is set, for instance\nwhen the following BNGL code is loaded:\n\nbegin observables\n   Molecules Acount A\nend observables\n\nthe name is set to Acount.\n")
      .def_property("file_name", &Count::get_file_name, &Count::set_file_name, "File name where this observable values will be stored.\nThe usual file_name, while using this count's name is:\nfile_name = './react_data/seed_' + str(SEED).zfill(5) + '/' + name + '.dat'.\n")
      .def_property("count_expression", &Count::get_count_expression, &Count::set_count_expression, "If this count uses more than one count term one must define a count term expression \nand set a reference to it to this attribute.  \nThe count expression must be composed only from CountTerm objects that are added or \nsubtracted.\n")
      .def_property("multiplier", &Count::get_multiplier, &Count::set_multiplier, "In some cases it might be useful to multiply the whole count by a constant to get \nfor instance concentration. The count_expression allows only addition and subtraction \nof count terms so such multiplication can be done through this attribute.\nIt can be also used to divide the resulting count by passing an inverse of the divisor (1/d).   \n")
      .def_property("every_n_timesteps", &Count::get_every_n_timesteps, &Count::set_every_n_timesteps, "Specifies periodicity of this count's output.\nValue is truncated (floored) to an integer.\nIf value is set to 0, this Count is used only on-demand through calls to its\nget_current_value method.  \n")
    ;
}

std::string GenCount::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("count") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("count")));
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
  ss << "m.Count(" << nl;
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
  if (name != STR_UNSET) {
    ss << ind << "name = " << "'" << name << "'" << "," << nl;
  }
  if (file_name != STR_UNSET) {
    ss << ind << "file_name = " << "'" << file_name << "'" << "," << nl;
  }
  if (is_set(count_expression)) {
    ss << ind << "count_expression = " << count_expression->export_to_python(out, ctx) << "," << nl;
  }
  if (multiplier != 1) {
    ss << ind << "multiplier = " << f_to_str(multiplier) << "," << nl;
  }
  if (every_n_timesteps != 1) {
    ss << ind << "every_n_timesteps = " << f_to_str(every_n_timesteps) << "," << nl;
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

