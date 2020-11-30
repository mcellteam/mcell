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
}

bool GenReactionRule::__eq__(const ReactionRule& other) const {
  return
    name == other.name &&
    vec_ptr_eq(reactants, other.reactants) &&
    vec_ptr_eq(products, other.products) &&
    fwd_rate == other.fwd_rate &&
    rev_name == other.rev_name &&
    rev_rate == other.rev_rate &&
    variable_rate == other.variable_rate;
}

bool GenReactionRule::eq_nonarray_attributes(const ReactionRule& other, const bool ignore_name) const {
  return
    (ignore_name || name == other.name) &&
    true /*reactants*/ &&
    true /*products*/ &&
    fwd_rate == other.fwd_rate &&
    rev_name == other.rev_name &&
    rev_rate == other.rev_rate &&
    true /*variable_rate*/;
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
      "variable_rate=" << vec_nonptr_to_str(variable_rate, ind + "  ");
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
            const std::vector<std::vector<float_t>>
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("reactants") = std::vector<std::shared_ptr<Complex>>(),
          py::arg("products") = std::vector<std::shared_ptr<Complex>>(),
          py::arg("fwd_rate") = FLT_UNSET,
          py::arg("rev_name") = STR_UNSET,
          py::arg("rev_rate") = FLT_UNSET,
          py::arg("variable_rate") = std::vector<std::vector<float_t>>()
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
    ;
}

} // namespace API
} // namespace MCell

