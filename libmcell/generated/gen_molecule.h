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

#ifndef API_GEN_MOLECULE_H
#define API_GEN_MOLECULE_H

#include "api/common.h"
#include "api/base_introspection_class.h"

namespace MCell {
namespace API {

class Molecule;
class Species;

#define MOLECULE_CTOR_NOARGS() \
    Molecule( \
    ) { \
      class_name = "Molecule"; \
      id = MOLECULE_ID_INVALID; \
      species = nullptr; \
      pos3d = VEC3_UNSET; \
      orientation = Orientation::NOT_SET; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenMolecule: public BaseIntrospectionClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  virtual bool __eq__(const Molecule& other) const;
  virtual bool eq_nonarray_attributes(const Molecule& other, const bool ignore_name = false) const;
  bool operator == (const Molecule& other) const { return __eq__(other);}
  bool operator != (const Molecule& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  int id;
  virtual void set_id(const int new_id_) {
    if (initialized) {
      throw RuntimeError("Value 'id' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    id = new_id_;
  }
  virtual int get_id() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return id;
  }

  std::shared_ptr<Species> species;
  virtual void set_species(std::shared_ptr<Species> new_species_) {
    if (initialized) {
      throw RuntimeError("Value 'species' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    species = new_species_;
  }
  virtual std::shared_ptr<Species> get_species() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return species;
  }

  Vec3 pos3d;
  virtual void set_pos3d(const Vec3& new_pos3d_) {
    if (initialized) {
      throw RuntimeError("Value 'pos3d' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    pos3d = new_pos3d_;
  }
  virtual const Vec3& get_pos3d() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return pos3d;
  }

  Orientation orientation;
  virtual void set_orientation(const Orientation new_orientation_) {
    if (initialized) {
      throw RuntimeError("Value 'orientation' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    orientation = new_orientation_;
  }
  virtual Orientation get_orientation() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return orientation;
  }

  // --- methods ---
  virtual void remove() = 0;
}; // GenMolecule

class Molecule;
py::class_<Molecule> define_pybinding_Molecule(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MOLECULE_H
