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

#ifndef API_GEN_SUBSYSTEM_H
#define API_GEN_SUBSYSTEM_H

#include "../api/common.h"
#include "../api/base_data_class.h"

namespace MCell {
namespace API {

class ReactionRule;
class Species;

class GenSubsystem {
public:
  virtual ~GenSubsystem() {}
  // --- attributes ---
  std::vector<std::shared_ptr<ReactionRule>> reaction_rules;
  virtual void set_reaction_rules(const std::vector<std::shared_ptr<ReactionRule>> new_reaction_rules_) {
    reaction_rules = new_reaction_rules_;
  }
  virtual std::vector<std::shared_ptr<ReactionRule>> get_reaction_rules() const {
    return reaction_rules;
  }

  std::vector<std::shared_ptr<Species>> species;
  virtual void set_species(const std::vector<std::shared_ptr<Species>> new_species_) {
    species = new_species_;
  }
  virtual std::vector<std::shared_ptr<Species>> get_species() const {
    return species;
  }

  // --- methods ---
  virtual void add_species(std::shared_ptr<Species> s) = 0;
  virtual std::shared_ptr<Species> find_species(const std::string& name) = 0;
  virtual void add_reaction_rule(std::shared_ptr<ReactionRule> s) = 0;
  virtual std::shared_ptr<ReactionRule> find_reaction_rule(const std::string& name) = 0;
}; // GenSubsystem

class Subsystem;
py::class_<Subsystem> define_pybinding_Subsystem(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SUBSYSTEM_H
