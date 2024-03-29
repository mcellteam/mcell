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

#ifndef API_GEN_SUBSYSTEM_H
#define API_GEN_SUBSYSTEM_H

#include "api/api_common.h"
#include "api/base_export_class.h"

namespace MCell {
namespace API {

class Subsystem;
class ElementaryMoleculeType;
class ReactionRule;
class Species;
class SurfaceClass;
class PythonExportContext;

#define SUBSYSTEM_CTOR() \
    Subsystem( \
        const std::vector<std::shared_ptr<Species>> species_ = std::vector<std::shared_ptr<Species>>(), \
        const std::vector<std::shared_ptr<ReactionRule>> reaction_rules_ = std::vector<std::shared_ptr<ReactionRule>>(), \
        const std::vector<std::shared_ptr<SurfaceClass>> surface_classes_ = std::vector<std::shared_ptr<SurfaceClass>>(), \
        const std::vector<std::shared_ptr<ElementaryMoleculeType>> elementary_molecule_types_ = std::vector<std::shared_ptr<ElementaryMoleculeType>>() \
    ) { \
      species = species_; \
      reaction_rules = reaction_rules_; \
      surface_classes = surface_classes_; \
      elementary_molecule_types = elementary_molecule_types_; \
    } \
    Subsystem(DefaultCtorArgType){ \
    }

class GenSubsystem: public BaseExportClass {
public:
  GenSubsystem() {
  }
  GenSubsystem(DefaultCtorArgType) {
  }
  virtual ~GenSubsystem() {}
  std::shared_ptr<Subsystem> copy_subsystem() const;
  std::shared_ptr<Subsystem> deepcopy_subsystem(py::dict = py::dict()) const;
  virtual bool __eq__(const Subsystem& other) const;
  virtual bool eq_nonarray_attributes(const Subsystem& other, const bool ignore_name = false) const;
  bool operator == (const Subsystem& other) const { return __eq__(other);}
  bool operator != (const Subsystem& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const ;

  virtual std::string export_to_python(std::ostream& out, PythonExportContext& ctx);
  virtual std::string export_vec_species(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_reaction_rules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_surface_classes(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);
  virtual std::string export_vec_elementary_molecule_types(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name);


  // --- attributes ---
  std::vector<std::shared_ptr<Species>> species;
  virtual void set_species(const std::vector<std::shared_ptr<Species>> new_species_) {
    species = new_species_;
  }
  virtual std::vector<std::shared_ptr<Species>>& get_species() {
    return species;
  }

  std::vector<std::shared_ptr<ReactionRule>> reaction_rules;
  virtual void set_reaction_rules(const std::vector<std::shared_ptr<ReactionRule>> new_reaction_rules_) {
    reaction_rules = new_reaction_rules_;
  }
  virtual std::vector<std::shared_ptr<ReactionRule>>& get_reaction_rules() {
    return reaction_rules;
  }

  std::vector<std::shared_ptr<SurfaceClass>> surface_classes;
  virtual void set_surface_classes(const std::vector<std::shared_ptr<SurfaceClass>> new_surface_classes_) {
    surface_classes = new_surface_classes_;
  }
  virtual std::vector<std::shared_ptr<SurfaceClass>>& get_surface_classes() {
    return surface_classes;
  }

  std::vector<std::shared_ptr<ElementaryMoleculeType>> elementary_molecule_types;
  virtual void set_elementary_molecule_types(const std::vector<std::shared_ptr<ElementaryMoleculeType>> new_elementary_molecule_types_) {
    elementary_molecule_types = new_elementary_molecule_types_;
  }
  virtual std::vector<std::shared_ptr<ElementaryMoleculeType>>& get_elementary_molecule_types() {
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
  virtual void load_bngl_molecule_types_and_reaction_rules(const std::string& file_name, const std::map<std::string, double>& parameter_overrides = std::map<std::string, double>()) = 0;
}; // GenSubsystem

class Subsystem;
py::class_<Subsystem> define_pybinding_Subsystem(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SUBSYSTEM_H
