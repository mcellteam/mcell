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

#ifndef API_SUBSYSTEM_H
#define API_SUBSYSTEM_H

#include "../generated/gen_subsystem.h"
#include "../api/common.h"

namespace MCell {
namespace API {

class Subsystem: public GenSubsystem {
public:
  // from generated template
  void add_species(std::shared_ptr<Species> s) override {}
  std::shared_ptr<Species> find_species(const std::string& name) override {return nullptr;}

  void add_reaction_rule(std::shared_ptr<ReactionRule> s) override {};
  std::shared_ptr<ReactionRule> find_reaction_rule(const std::string& name) override {return nullptr;}

  // added manually
  void dump() const;
};

} // namespace API
} // namespace MCell

#endif // API_SUBSYSTEM_H
