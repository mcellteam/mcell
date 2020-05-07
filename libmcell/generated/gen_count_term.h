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

#ifndef API_GEN_COUNT_TERM_H
#define API_GEN_COUNT_TERM_H

#include "../api/common.h"

namespace MCell {
namespace API {

class GeometryObject;
class ReactionRule;
class Species;

#define COUNT_TERM_CTOR() \
    CountTerm( \
        std::shared_ptr<Species> species_ = nullptr, \
        std::shared_ptr<ReactionRule> reaction_rule_ = nullptr, \
        std::shared_ptr<GeometryObject> enclosed_in_object_ = nullptr \
    ) { \
      class_name = "CountTerm"; \
      species = species_; \
      reaction_rule = reaction_rule_; \
      enclosed_in_object = enclosed_in_object_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenCountTerm: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::shared_ptr<Species> species;
  virtual void set_species(std::shared_ptr<Species> new_species_) {
    species = new_species_;
  }
  virtual std::shared_ptr<Species> get_species() const {
    return species;
  }

  std::shared_ptr<ReactionRule> reaction_rule;
  virtual void set_reaction_rule(std::shared_ptr<ReactionRule> new_reaction_rule_) {
    reaction_rule = new_reaction_rule_;
  }
  virtual std::shared_ptr<ReactionRule> get_reaction_rule() const {
    return reaction_rule;
  }

  std::shared_ptr<GeometryObject> enclosed_in_object;
  virtual void set_enclosed_in_object(std::shared_ptr<GeometryObject> new_enclosed_in_object_) {
    enclosed_in_object = new_enclosed_in_object_;
  }
  virtual std::shared_ptr<GeometryObject> get_enclosed_in_object() const {
    return enclosed_in_object;
  }

  // --- methods ---
}; // GenCountTerm

class CountTerm;
py::class_<CountTerm> define_pybinding_CountTerm(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_COUNT_TERM_H
