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
#include "gen_reaction_rule.h"
#include "api/reaction_rule.h"
#include "api/complex.h"

namespace MCell {
namespace API {

void GenReactionRule::check_semantics() const {
}

void GenReactionRule::set_initialized() {
  vec_set_initialized(reactants);
  vec_set_initialized(products);
  initialized = true;
}

void GenReactionRule::set_all_attributes_as_default_or_unset() {
  class_name = "ReactionRule";
  name = STR_UNSET;
  reactants = std::vector<std::shared_ptr<Complex>>();
  products = std::vector<std::shared_ptr<Complex>>();
  fwd_rate = FLT_UNSET;
  rev_name = STR_UNSET;
  rev_rate = FLT_UNSET;
  variable_rate = std::vector<std::vector<double>>();
  is_intermembrane_surface_reaction = false;
}

std::shared_ptr<ReactionRule> GenReactionRule::copy_reaction_rule() const {
  std::shared_ptr<ReactionRule> res = std::make_shared<ReactionRule>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  res->reactants = reactants;
  res->products = products;
  res->fwd_rate = fwd_rate;
  res->rev_name = rev_name;
  res->rev_rate = rev_rate;
  res->variable_rate = variable_rate;
  res->is_intermembrane_surface_reaction = is_intermembrane_surface_reaction;

  return res;
}

std::shared_ptr<ReactionRule> GenReactionRule::deepcopy_reaction_rule(py::dict) const {
  std::shared_ptr<ReactionRule> res = std::make_shared<ReactionRule>(DefaultCtorArgType());
  res->class_name = class_name;
  res->name = name;
  for (const auto& item: reactants) {
    res->reactants.push_back((is_set(item)) ? item->deepcopy_complex() : nullptr);
  }
  for (const auto& item: products) {
    res->products.push_back((is_set(item)) ? item->deepcopy_complex() : nullptr);
  }
  res->fwd_rate = fwd_rate;
  res->rev_name = rev_name;
  res->rev_rate = rev_rate;
  res->variable_rate = variable_rate;
  res->is_intermembrane_surface_reaction = is_intermembrane_surface_reaction;

  return res;
}

bool GenReactionRule::__eq__(const ReactionRule& other) const {
  return
    name == other.name &&
    vec_ptr_eq(reactants, other.reactants) &&
    vec_ptr_eq(products, other.products) &&
    fwd_rate == other.fwd_rate &&
    rev_name == other.rev_name &&
    rev_rate == other.rev_rate &&
    variable_rate == other.variable_rate &&
    is_intermembrane_surface_reaction == other.is_intermembrane_surface_reaction;
}

bool GenReactionRule::eq_nonarray_attributes(const ReactionRule& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    true /*reactants*/ &&
    true /*products*/ &&
    fwd_rate == other.fwd_rate &&
    rev_name == other.rev_name &&
    rev_rate == other.rev_rate &&
    true /*variable_rate*/ &&
    is_intermembrane_surface_reaction == other.is_intermembrane_surface_reaction;
}

std::string GenReactionRule::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "reactants=" << vec_ptr_to_str(reactants, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "products=" << vec_ptr_to_str(products, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "fwd_rate=" << fwd_rate << ", " <<
      "rev_name=" << rev_name << ", " <<
      "rev_rate=" << rev_rate << ", " <<
      "variable_rate=" << vec_nonptr_to_str(variable_rate, all_details, ind + "  ") << ", " <<
      "is_intermembrane_surface_reaction=" << is_intermembrane_surface_reaction;
  return ss.str();
}

py::class_<ReactionRule> define_pybinding_ReactionRule(py::module& m) {
  return py::class_<ReactionRule, std::shared_ptr<ReactionRule>>(m, "ReactionRule", "Represents a BioNetGen Language (BNGL) reaction rule. \nIn BNGL, a reaction is simply one or more transformations\napplied simultaneously to one or more species. The following\ntransformations (and their combinations) are allowed:\n  * Forming a bond, e.g. A(b) + B(a) -> A(b!0).B(a!0)\n  * Breaking a bond, e.g. A(b!0).B(a!0)-> A(b)+ B(a)\n  * Changing of component state, e.g. X(y~0) -> X(y~p)\n  * Creating a molecule, e.g. A(b) -> A(b) + C(d)\n  * Destroying a molecule, e.g. A(b) + B(a) -> A(b) or A -> Null \n    (Null here means that there is no product)\n  * Changing species of a bound molecule when the molecule type has the \n    same components, e.g. A(b!0).B(a!0)-> A(b!0).C(a!0)\n  \nAlso compartments may be specified in reactants (patterns) and for products.\nSpecial compartment classes supported by MCell4 are @IN and @OUT.\nThey can be used in surface reactions to constrain a reaction with a volume molecule \nhitting a surface molecule from the inside or outside of the compartment, \ne.g. A(s)@IN + S(a) -> S(a!1).A(s!1) and/or to define the location of the \nproduct, e.g. S(a!1).A(s!1) -> S(a) + A(s)@OUT.   \n")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<Complex>>,
            const std::vector<std::shared_ptr<Complex>>,
            const double,
            const std::string&,
            const double,
            const std::vector<std::vector<double>>,
            const bool
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("reactants") = std::vector<std::shared_ptr<Complex>>(),
          py::arg("products") = std::vector<std::shared_ptr<Complex>>(),
          py::arg("fwd_rate") = FLT_UNSET,
          py::arg("rev_name") = STR_UNSET,
          py::arg("rev_rate") = FLT_UNSET,
          py::arg("variable_rate") = std::vector<std::vector<double>>(),
          py::arg("is_intermembrane_surface_reaction") = false
      )
      .def("check_semantics", &ReactionRule::check_semantics)
      .def("__copy__", &ReactionRule::copy_reaction_rule)
      .def("__deepcopy__", &ReactionRule::deepcopy_reaction_rule, py::arg("memo"))
      .def("__str__", &ReactionRule::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &ReactionRule::__eq__, py::arg("other"))
      .def("to_bngl_str", &ReactionRule::to_bngl_str, "Creates a string that corresponds to the reaction rule's BNGL representation, does not contain rates.")
      .def("dump", &ReactionRule::dump)
      .def_property("name", &ReactionRule::get_name, &ReactionRule::set_name, "Name of the reaction. If this is a reversible reaction, then it is the name of the \nreaction in forward direction.\n")
      .def_property("reactants", &ReactionRule::get_reactants, &ReactionRule::set_reactants, py::return_value_policy::reference, "List of reactant patterns. Must contain one or two patterns.")
      .def_property("products", &ReactionRule::get_products, &ReactionRule::set_products, py::return_value_policy::reference, "List of reactant patterns. May be empty.")
      .def_property("fwd_rate", &ReactionRule::get_fwd_rate, &ReactionRule::set_fwd_rate, "The units of the reaction rate for uni- and bimolecular reactions are:\n  * [s^-1] for unimolecular reactions,\n  * [N^-1*s^-1] bimolecular reactions between two surface molecules on different objects \n    (this is a highly experimental feature and the unit will likely change in the future, \n     not sure if probability is computed correctly, it works the way that the surface molecule \n     is first diffused and then a potential collisions within the distance of Config.intermembrane_interaction_radius\n     are evaluated). \nOther bimolecular reaction units depend on Model.config.use_bng_units settings.\nWhen use_bng_units is False (default), traditional MCell units are used:  \n  * [M^-1*s^-1] for bimolecular reactions between either two volume molecules, a volume molecule \n                and a surface (molecule), \n  * [um^2*N^-1*s^-1] bimolecular reactions between two surface molecules on the same surface, and\nWhen use_bng_units is True, units compatible with BioNetGen's ODE, SSA, and PLA solvers are used:\n  * [um^3*N^-1*s^-1] for any bimolecular reactions, surface-surface reaction rate conversion assumes 10nm membrane thickness. \nM is the molarity of the solution and N the number of reactants.\nMay be changed after model initialization. \nSetting of value is ignored if the rate does not change. \nIf the new value differs from previous, updates all information related \nto the new rate including recomputation of reaction times for molecules if this is a\nunimolecular reaction.\n")
      .def_property("rev_name", &ReactionRule::get_rev_name, &ReactionRule::set_rev_name, "Name of the reaction in reverse direction.   \n")
      .def_property("rev_rate", &ReactionRule::get_rev_rate, &ReactionRule::set_rev_rate, "Reverse reactions rate, reaction is unidirectional when not specified.\nMay be changed after model initialization, in the case behaves the same was as for \nchanging the 'fwd_rate'. \nUses the same units as 'fwd_rate'.\n")
      .def_property("variable_rate", &ReactionRule::get_variable_rate, &ReactionRule::set_variable_rate, py::return_value_policy::reference, "The array passed as this argument must have as its items a pair of floats (time in s, rate).\nMust be sorted by time (this is not checked).      \nVariable rate is applicable only for irreversible reactions.\nWhen simulation starts and the table does not contain value for time 0, the initial fwd_rate is set to 0.\nWhen time advances after the last time in this table, the last rate is used for all subsequent iterations.   \nMembers fwd_rate and rev_rate must not be set when setting this attribute through a constructor. \nWhen this attribute is set outside of the class constructor, fwd_rate is automatically reset to an 'unset' value.\nCannot be set after model initialization. \n")
      .def_property("is_intermembrane_surface_reaction", &ReactionRule::get_is_intermembrane_surface_reaction, &ReactionRule::set_is_intermembrane_surface_reaction, "Experimental, see addintinal explanation in 'fwd' rate.\nThen set to true, this is a special type of surface-surface reaction that \nallows for two surface molecules to react when they are on different geometrical objects. \nNote: This support is limited for now, the reaction rule must be in the form of A + B -> C + D \nwhere all reactants and products must be surface molecules and their orientation must be 'any' (default). \n")
    ;
}

std::string GenReactionRule::export_to_python(std::ostream& out, PythonExportContext& ctx) {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = std::string("reaction_rule") + "_" + (is_set(name) ? fix_id(name) : std::to_string(ctx.postinc_counter("reaction_rule")));
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
  ss << "m.ReactionRule(" << nl;
  if (name != STR_UNSET) {
    ss << ind << "name = " << "'" << name << "'" << "," << nl;
  }
  if (reactants != std::vector<std::shared_ptr<Complex>>() && !skip_vectors_export()) {
    ss << ind << "reactants = " << export_vec_reactants(out, ctx, exported_name) << "," << nl;
  }
  if (products != std::vector<std::shared_ptr<Complex>>() && !skip_vectors_export()) {
    ss << ind << "products = " << export_vec_products(out, ctx, exported_name) << "," << nl;
  }
  if (fwd_rate != FLT_UNSET) {
    ss << ind << "fwd_rate = " << f_to_str(fwd_rate) << "," << nl;
  }
  if (rev_name != STR_UNSET) {
    ss << ind << "rev_name = " << "'" << rev_name << "'" << "," << nl;
  }
  if (rev_rate != FLT_UNSET) {
    ss << ind << "rev_rate = " << f_to_str(rev_rate) << "," << nl;
  }
  if (variable_rate != std::vector<std::vector<double>>() && !skip_vectors_export()) {
    ss << ind << "variable_rate = " << export_vec_variable_rate(out, ctx, exported_name) << "," << nl;
  }
  if (is_intermembrane_surface_reaction != false) {
    ss << ind << "is_intermembrane_surface_reaction = " << is_intermembrane_surface_reaction << "," << nl;
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

std::string GenReactionRule::export_vec_reactants(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < reactants.size(); i++) {
    const auto& item = reactants[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

std::string GenReactionRule::export_vec_products(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < products.size(); i++) {
    const auto& item = products[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

std::string GenReactionRule::export_vec_variable_rate(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // does not print the array itself to 'out' and returns the whole list
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < variable_rate.size(); i++) {
    const auto& item = variable_rate[i];
    if (i == 0) {
      ss << " ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    ss << "[";
    for (const auto& value: item) {
      ss << f_to_str(value) << ", ";
    }
    ss << "], ";
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

