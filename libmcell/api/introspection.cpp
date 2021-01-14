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

#include "api/introspection.h"

#include "api/model.h"
#include "api/api_utils.h"
#include "api/molecule.h"
#include "api/wall.h"
#include "api/geometry_object.h"

#include "world.h"
#include "src4/geometry_utils.h" // TODO: we should rename the src4 directory


using namespace std;

namespace MCell {
namespace API {

void Introspection::initialize_introspection(Model* model_) {
  model_inst = model_;
  world = model_inst->get_world();
  assert(world != nullptr);
}


std::vector<int> Introspection::get_molecule_ids(std::shared_ptr<Species> species) {
  // NOTE: not very efficient
  std::vector<int> res;

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  std::vector<MCell::Molecule>& molecules = p.get_molecules();
  for (MCell::Molecule& m: molecules) {
    if (m.is_defunct()) {
      continue;
    }

    if (is_set(species) && species->species_id == m.species_id) {
      res.push_back(m.id);
    }
    else {
      res.push_back(m.id);
    }
  }

  return res;
}


std::shared_ptr<API::Molecule> Introspection::get_molecule(const int id) {
  std::shared_ptr<API::Molecule> res;
  Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  if (!p.does_molecule_exist(id)) {
    throw RuntimeError("Molecule with id " + to_string(id) + " does not exist.");
  }
  MCell::Molecule& m = p.get_m(id);
  assert(!m.is_defunct() && "is_defunct is checked already by does_molecule_exist()");

  res = make_shared<API::Molecule>();
  res->id = m.id;
  res->species_id = m.species_id;
  if (m.is_surf()) {
    const MCell::Wall& w = p.get_wall(m.s.wall_index);
    const Vec3& w_vert0 = p.get_wall_vertex(w, 0);

    res->type = MoleculeType::SURFACE;
    res->pos2d = m.s.pos * Vec2(world->config.length_unit);
    res->pos3d = GeometryUtil::uv2xyz(m.s.pos, w, w_vert0) * Vec3(world->config.length_unit);
    res->orientation = convert_orientation(m.s.orientation);

    res->geometry_object = model_inst->get_geometry_object_with_id(w.object_id);
    assert(is_set(res->geometry_object));
    res->wall_index = m.s.wall_index - res->geometry_object->first_wall_index;
  }
  else {
    res->type = MoleculeType::VOLUME;
    res->pos3d = m.v.pos * Vec3(world->config.length_unit);
    res->orientation = Orientation::NONE;
  }
  res->world = world;
  res->set_initialized();

  return res;
}


std::string Introspection::get_species_name(const int species_id) {
  if (!world->get_all_species().is_valid_id(species_id)) {
    throw RuntimeError("Species with id " + to_string(species_id) + " does not exist.");
  }
  return world->get_all_species().get(species_id).name;
}


Vec3 Introspection::get_vertex(std::shared_ptr<GeometryObject> object, const int vertex_index) {
  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  return
      p.get_geometry_vertex(object->get_partition_vertex_index(vertex_index)) *
      Vec3(world->config.length_unit);
}


std::shared_ptr<Wall> Introspection::get_wall(std::shared_ptr<GeometryObject> object, const int wall_index) {
  object->check_is_initialized();

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  const MCell::Wall& w = p.get_wall(object->get_partition_wall_index(wall_index));

  auto res = make_shared<Wall>();
  res->geometry_object = object;
  res->wall_index = wall_index;
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    res->vertices.push_back(p.get_geometry_vertex(w.vertex_indices[i]) * Vec3(world->config.length_unit));
  }
  res->area = w.area * world->config.length_unit * world->config.length_unit;
  res->unit_normal = w.normal;
  assert(cmp_eq(len3(res->unit_normal), 1));
  res->is_movable = w.is_movable;
  res->world = world;
  return res;
}


Vec3 Introspection::get_vertex_unit_normal(std::shared_ptr<GeometryObject> object, const int vertex_index) {
  object->check_is_initialized();

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  const std::vector<wall_index_t>& walls = p.get_walls_using_vertex(object->get_partition_vertex_index(vertex_index));

  if (walls.empty()) {
    throw RuntimeError("Internal error: there are no walls that use vertex with index " +
        to_string(vertex_index) + " of object " + object->name + ".");
  }

  Vec3 normals_sum = Vec3(0);
  for (wall_index_t wi: walls) {
    const MCell::Wall& w = p.get_wall(wi);
    // wall normals are already unit vectors so we can just sum them
    normals_sum = normals_sum + w.normal;
  }

  return normals_sum / Vec3(len3(normals_sum));
}


Vec3 Introspection::get_wall_unit_normal(std::shared_ptr<GeometryObject> object, const int wall_index) {
  object->check_is_initialized();

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  const MCell::Wall& w = p.get_wall(object->get_partition_wall_index(wall_index));

  // the value is is already normalized
  assert(cmp_eq(len3(w.normal), 1));
  return w.normal;
}


} // namespace API
} // namespace MCell
