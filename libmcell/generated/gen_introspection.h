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

#ifndef API_GEN_INTROSPECTION_H
#define API_GEN_INTROSPECTION_H

#include "api/api_common.h"

namespace MCell {
namespace API {

class Introspection;
class Color;
class Complex;
class GeometryObject;
class Molecule;
class Wall;
class PythonExportContext;

class GenIntrospection {
public:
  GenIntrospection() {
  }
  GenIntrospection(DefaultCtorArgType) {
  }
  virtual ~GenIntrospection() {}
  std::shared_ptr<Introspection> copy_introspection() const;
  std::shared_ptr<Introspection> deepcopy_introspection(py::dict = py::dict()) const;
  virtual bool __eq__(const Introspection& other) const;
  virtual bool eq_nonarray_attributes(const Introspection& other, const bool ignore_name = false) const;
  bool operator == (const Introspection& other) const { return __eq__(other);}
  bool operator != (const Introspection& other) const { return !__eq__(other);}
  std::string to_str(const bool all_details=false, const std::string ind="") const ;

  // --- attributes ---
  // --- methods ---
  virtual std::vector<int> get_molecule_ids(std::shared_ptr<Complex> pattern = nullptr) = 0;
  virtual std::shared_ptr<Molecule> get_molecule(const int id) = 0;
  virtual std::string get_species_name(const int species_id) = 0;
  virtual Vec3 get_vertex(std::shared_ptr<GeometryObject> object, const int vertex_index) = 0;
  virtual std::shared_ptr<Wall> get_wall(std::shared_ptr<GeometryObject> object, const int wall_index) = 0;
  virtual Vec3 get_vertex_unit_normal(std::shared_ptr<GeometryObject> object, const int vertex_index) = 0;
  virtual Vec3 get_wall_unit_normal(std::shared_ptr<GeometryObject> object, const int wall_index) = 0;
  virtual std::shared_ptr<Color> get_wall_color(std::shared_ptr<GeometryObject> object, const int wall_index) = 0;
  virtual void set_wall_color(std::shared_ptr<GeometryObject> object, const int wall_index, std::shared_ptr<Color> color) = 0;
}; // GenIntrospection

class Introspection;
py::class_<Introspection> define_pybinding_Introspection(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_INTROSPECTION_H
