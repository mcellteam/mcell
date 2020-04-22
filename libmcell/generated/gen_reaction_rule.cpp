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
#include "gen_reaction_rule.h"
#include "../api/reaction_rule.h"
#include "../api/complex_instance.h"

namespace MCell {
namespace API {

SemRes GenReactionRule::check_semantics(std::ostream& out) const {
  return SemRes::OK;
}

std::string GenReactionRule::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << get_object_name() << ": " <<
      "name=" << name << ", " <<
      "\n" << ind + "  " << "reactants=" << vec_ptr_to_str(reactants, ind + "  ") << ", " << "\n" << ind + "  " <<
      "products=" << vec_ptr_to_str(products, ind + "  ") << ", " << "\n" << ind + "  " <<
      "fwd_rate=" << fwd_rate << ", " <<
      "rev_rate=" << rev_rate;
  return ss.str();
}

py::class_<ReactionRule> define_pybinding_ReactionRule(py::module& m) {
  return py::class_<ReactionRule, std::shared_ptr<ReactionRule>>(m, "ReactionRule")
      .def(
          py::init<
            const std::string&,
            const std::vector<std::shared_ptr<ComplexInstance>>,
            const std::vector<std::shared_ptr<ComplexInstance>>,
            const float_t,
            const float_t
          >(),
          py::arg("name") = STR_UNSET,
          py::arg("reactants") = std::vector<std::shared_ptr<ComplexInstance>>(),
          py::arg("products") = std::vector<std::shared_ptr<ComplexInstance>>(),
          py::arg("fwd_rate") = FLT_UNSET,
          py::arg("rev_rate") = FLT_UNSET
        )
      .def("check_semantics", &ReactionRule::check_semantics_cerr)
      .def("__str__", &ReactionRule::to_str, py::arg("ind") = std::string(""))
      .def("dump", &ReactionRule::dump)
      .def_property("name", &ReactionRule::get_name, &ReactionRule::set_name)
      .def_property("reactants", &ReactionRule::get_reactants, &ReactionRule::set_reactants)
      .def_property("products", &ReactionRule::get_products, &ReactionRule::set_products)
      .def_property("fwd_rate", &ReactionRule::get_fwd_rate, &ReactionRule::set_fwd_rate)
      .def_property("rev_rate", &ReactionRule::get_rev_rate, &ReactionRule::set_rev_rate)
    ;
}

} // namespace API
} // namespace MCell

