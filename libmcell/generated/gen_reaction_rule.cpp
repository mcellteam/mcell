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
  variable_rate = std::vector<std::vector<float_t>>();
  is_intermembrane_surface_reaction = false;
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

std::string GenReactionRule::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "reactants=" << vec_ptr_to_str(reactants, ind + "  ") << ", " << "\n" << ind + "  " <<
      "products=" << vec_ptr_to_str(products, ind + "  ") << ", " << "\n" << ind + "  " <<
      "fwd_rate=" << fwd_rate << ", " <<
      "rev_name=" << rev_name << ", " <<
      "rev_rate=" << rev_rate << ", " <<
      "variable_rate=" << vec_nonptr_to_str(variable_rate, ind + "  ") << ", " <<
      "is_intermembrane_surface_reaction=" << is_intermembrane_surface_reaction;
  return ss.str();
}

py::class_<ReactionRule> define_pybinding_ReactionRule(py::module& m) {
  return py::class_<ReactionRule, std::shared_ptr<ReactionRule>>(m, "ReactionRule")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<Complex>>,
            const std::vector<std::shared_ptr<Complex>>,
            const float_t,
            const std::string&,
            const float_t,
            const std::vector<std::vector<float_t>>,
            const bool
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("reactants") = std::vector<std::shared_ptr<Complex>>(),
          py::arg("products") = std::vector<std::shared_ptr<Complex>>(),
          py::arg("fwd_rate") = FLT_UNSET,
          py::arg("rev_name") = STR_UNSET,
          py::arg("rev_rate") = FLT_UNSET,
          py::arg("variable_rate") = std::vector<std::vector<float_t>>(),
          py::arg("is_intermembrane_surface_reaction") = false
      )
      .def("check_semantics", &ReactionRule::check_semantics)
      .def("__str__", &ReactionRule::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &ReactionRule::__eq__, py::arg("other"))
      .def("to_bngl_str", &ReactionRule::to_bngl_str)
      .def("dump", &ReactionRule::dump)
      .def_property("name", &ReactionRule::get_name, &ReactionRule::set_name)
      .def_property("reactants", &ReactionRule::get_reactants, &ReactionRule::set_reactants)
      .def_property("products", &ReactionRule::get_products, &ReactionRule::set_products)
      .def_property("fwd_rate", &ReactionRule::get_fwd_rate, &ReactionRule::set_fwd_rate)
      .def_property("rev_name", &ReactionRule::get_rev_name, &ReactionRule::set_rev_name)
      .def_property("rev_rate", &ReactionRule::get_rev_rate, &ReactionRule::set_rev_rate)
      .def_property("variable_rate", &ReactionRule::get_variable_rate, &ReactionRule::set_variable_rate)
      .def_property("is_intermembrane_surface_reaction", &ReactionRule::get_is_intermembrane_surface_reaction, &ReactionRule::set_is_intermembrane_surface_reaction)
    ;
}

std::string GenReactionRule::export_to_python(std::ostream& out, PythonExportContext& ctx) const {
  if (!export_even_if_already_exported() && ctx.already_exported(this)) {
    return ctx.get_exported_name(this);
  }
  std::string exported_name = "reaction_rule_" + fix_id(name);
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
    ss << ind << "fwd_rate = " << fwd_rate << "," << nl;
  }
  if (rev_name != STR_UNSET) {
    ss << ind << "rev_name = " << "'" << rev_name << "'" << "," << nl;
  }
  if (rev_rate != FLT_UNSET) {
    ss << ind << "rev_rate = " << rev_rate << "," << nl;
  }
  if (variable_rate != std::vector<std::vector<float_t>>() && !skip_vectors_export()) {
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

std::string GenReactionRule::export_vec_reactants(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
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

std::string GenReactionRule::export_vec_products(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
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

std::string GenReactionRule::export_vec_variable_rate(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) const {
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
      ss << value << ", ";
    }
    ss << "], ";
  }
  ss << "]";
  return ss.str();
}

} // namespace API
} // namespace MCell

