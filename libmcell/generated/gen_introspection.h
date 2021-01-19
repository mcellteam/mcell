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

#ifndef API_GEN_INTROSPECTION_H
#define API_GEN_INTROSPECTION_H

#include "api/api_common.h"

namespace MCell {
namespace API {

class Introspection;
class GeometryObject;
class Molecule;
class Species;
class Wall;
class PythonExportContext;

class GenIntrospection {
public:
  virtual ~GenIntrospection() {}
  virtual bool __eq__(const Introspection& other) const;
  virtual bool eq_nonarray_attributes(const Introspection& other, const bool ignore_name = false) const;
  bool operator == (const Introspection& other) const { return __eq__(other);}
  bool operator != (const Introspection& other) const { return !__eq__(other);}
  std::string to_str(const std::string ind="") const ;

  // --- attributes ---
  // --- methods ---
  virtual std::vector<int> get_molecule_ids(std::shared_ptr<Species> species = nullptr) = 0;
  virtual std::shared_ptr<Molecule> get_molecule(const int id) = 0;
  virtual std::string get_species_name(const int species_id) = 0;
  virtual Vec3 get_vertex(std::shared_ptr<GeometryObject> object, const int vertex_index) = 0;
  virtual std::shared_ptr<Wall> get_wall(std::shared_ptr<GeometryObject> object, const int wall_index) = 0;
  virtual Vec3 get_vertex_unit_normal(std::shared_ptr<GeometryObject> object, const int vertex_index) = 0;
  virtual Vec3 get_wall_unit_normal(std::shared_ptr<GeometryObject> object, const int wall_index) = 0;
}; // GenIntrospection

class Introspection;
py::class_<Introspection> define_pybinding_Introspection(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_INTROSPECTION_H
