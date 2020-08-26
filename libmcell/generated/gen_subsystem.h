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

namespace MCell {
namespace API {

class ElementaryMoleculeType;
class ReactionRule;
class Species;
class SurfaceClass;

class GenSubsystem {
public:
  virtual ~GenSubsystem() {}
  std::string to_str(const std::string ind="") const ;

  // --- attributes ---
  std::vector<std::shared_ptr<Species>> species;
  virtual void set_species(const std::vector<std::shared_ptr<Species>> new_species_) {
    species = new_species_;
  }
  virtual std::vector<std::shared_ptr<Species>> get_species() const {
    return species;
  }

  std::vector<std::shared_ptr<ReactionRule>> reaction_rules;
  virtual void set_reaction_rules(const std::vector<std::shared_ptr<ReactionRule>> new_reaction_rules_) {
    reaction_rules = new_reaction_rules_;
  }
  virtual std::vector<std::shared_ptr<ReactionRule>> get_reaction_rules() const {
    return reaction_rules;
  }

  std::vector<std::shared_ptr<SurfaceClass>> surface_classes;
  virtual void set_surface_classes(const std::vector<std::shared_ptr<SurfaceClass>> new_surface_classes_) {
    surface_classes = new_surface_classes_;
  }
  virtual std::vector<std::shared_ptr<SurfaceClass>> get_surface_classes() const {
    return surface_classes;
  }

  std::vector<std::shared_ptr<ElementaryMoleculeType>> elementary_molecule_types;
  virtual void set_elementary_molecule_types(const std::vector<std::shared_ptr<ElementaryMoleculeType>> new_elementary_molecule_types_) {
    elementary_molecule_types = new_elementary_molecule_types_;
  }
  virtual std::vector<std::shared_ptr<ElementaryMoleculeType>> get_elementary_molecule_types() const {
    return elementary_molecule_types;
  }

  // --- methods ---
  virtual void add_species(std::shared_ptr<Species> s) = 0;
  virtual std::shared_ptr<Species> find_species(const std::string& name) = 0;
  virtual void add_reaction_rule(std::shared_ptr<ReactionRule> r) = 0;
  virtual std::shared_ptr<ReactionRule> find_reaction_rule(const std::string& name) = 0;
  virtual void add_surface_class(std::shared_ptr<SurfaceClass> sc) = 0;
  virtual std::shared_ptr<SurfaceClass> find_surface_class(const std::string& name) = 0;
  virtual void add_elementary_molecule_type(std::shared_ptr<ElementaryMoleculeType> mt) = 0;
  virtual std::shared_ptr<ElementaryMoleculeType> find_elementary_molecule_type(const std::string& name) = 0;
  virtual void load_bngl_molecule_types_and_reaction_rules(const std::string& filename) = 0;
}; // GenSubsystem

class Subsystem;
py::class_<Subsystem> define_pybinding_Subsystem(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SUBSYSTEM_H
