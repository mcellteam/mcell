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

#ifndef API_INTROSPECTION_H
#define API_INTROSPECTION_H

#include "generated/gen_introspection.h"
#include "api/api_common.h"

namespace MCell {

class World;

namespace API {

class Model;

// this class is used only as a base to Model, it is not provided through API
class Introspection: public GenIntrospection {
public:
  Introspection() :
    model_inst(nullptr),
    world(nullptr) {
  }

  void initialize_introspection(Model* model_);

  std::vector<int> get_molecule_ids(std::shared_ptr<Complex> pattern = nullptr) override;
  std::shared_ptr<Molecule> get_molecule(const int id) override;

  std::string get_species_name(const int species_id) override;

  Vec3 get_vertex(std::shared_ptr<GeometryObject> object, const int vertex_index) override;
  std::shared_ptr<Wall> get_wall(std::shared_ptr<GeometryObject> object, const int wall_index) override;

  Vec3 get_vertex_unit_normal(std::shared_ptr<GeometryObject> object, const int vertex_index) override;
  Vec3 get_wall_unit_normal(std::shared_ptr<GeometryObject> object, const int wall_index) override;

  std::shared_ptr<Color> get_wall_color(std::shared_ptr<GeometryObject> object, const int wall_index) override;
  void set_wall_color(std::shared_ptr<GeometryObject> object, const int wall_index, std::shared_ptr<Color> color) override;

  void dump() const {}

private:
  // not using name model because class Model inherits Introspection and
  // this made code a bit confusing
  Model* model_inst;

  World* world;
};

} // namespace API
} // namespace MCell

#endif // API_INTROSPECTION_H
