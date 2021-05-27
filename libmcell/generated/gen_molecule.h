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

#include "api/api_common.h"
#include "api/base_introspection_class.h"

namespace MCell {
namespace API {

class Molecule;
class GeometryObject;
class PythonExportContext;

#define MOLECULE_CTOR_NOARGS() \
    Molecule( \
    ) { \
      class_name = "Molecule"; \
      id = ID_INVALID; \
      type = MoleculeType::UNSET; \
      species_id = ID_INVALID; \
      pos3d = VEC3_UNSET; \
      orientation = Orientation::NOT_SET; \
      pos2d = VEC2_UNSET; \
      geometry_object = nullptr; \
      wall_index = -1; \
      postprocess_in_ctor(); \
      check_semantics(); \
    } \
    Molecule(DefaultCtorArgType) : \
      GenMolecule(DefaultCtorArgType()) { \
      set_all_attributes_as_default_or_unset(); \
    }

class GenMolecule: public BaseIntrospectionClass {
public:
  GenMolecule() {
  }
  GenMolecule(DefaultCtorArgType) {
  }
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  void set_initialized() override;
  void set_all_attributes_as_default_or_unset() override;

  Molecule copy_molecule() const;
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

  MoleculeType type;
  virtual void set_type(const MoleculeType new_type_) {
    if (initialized) {
      throw RuntimeError("Value 'type' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    type = new_type_;
  }
  virtual MoleculeType get_type() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return type;
  }

  int species_id;
  virtual void set_species_id(const int new_species_id_) {
    if (initialized) {
      throw RuntimeError("Value 'species_id' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    species_id = new_species_id_;
  }
  virtual int get_species_id() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return species_id;
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

  Vec2 pos2d;
  virtual void set_pos2d(const Vec2& new_pos2d_) {
    if (initialized) {
      throw RuntimeError("Value 'pos2d' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    pos2d = new_pos2d_;
  }
  virtual const Vec2& get_pos2d() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return pos2d;
  }

  std::shared_ptr<GeometryObject> geometry_object;
  virtual void set_geometry_object(std::shared_ptr<GeometryObject> new_geometry_object_) {
    if (initialized) {
      throw RuntimeError("Value 'geometry_object' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    geometry_object = new_geometry_object_;
  }
  virtual std::shared_ptr<GeometryObject> get_geometry_object() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return geometry_object;
  }

  int wall_index;
  virtual void set_wall_index(const int new_wall_index_) {
    if (initialized) {
      throw RuntimeError("Value 'wall_index' of object with name " + name + " (class " + class_name + ") "
                         "cannot be set after model was initialized.");
    }
    cached_data_are_uptodate = false;
    wall_index = new_wall_index_;
  }
  virtual int get_wall_index() const {
    cached_data_are_uptodate = false; // arrays and other data can be modified through getters
    return wall_index;
  }

  // --- methods ---
  virtual void remove() = 0;
}; // GenMolecule

class Molecule;
py::class_<Molecule> define_pybinding_Molecule(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_MOLECULE_H
