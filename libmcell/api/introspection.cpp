/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "api/introspection.h"

#include "api/model.h"
#include "api/api_utils.h"
#include "api/molecule.h"
#include "api/wall.h"
#include "api/geometry_object.h"
#include "api/color.h"
#include "api/bng_converter.h"

#include "world.h"
#include "src4/geometry_utils.h" // TODO: we should rename the src4 directory
#include "bng/cplx.h"

using namespace std;

namespace MCell {
namespace API {

static void object_is_set_and_initialized(std::shared_ptr<GeometryObject> object) {
  if (!is_set(object)) {
    throw ValueError(S("Argument passed as ") + NAME_CLASS_GEOMETRY_OBJECT + " must not be None.");
  }
  object->check_is_initialized();
}

void Introspection::initialize_introspection(Model* model_) {
  model_inst = model_;
  world = model_inst->get_world();
  assert(world != nullptr);
}


std::vector<int> Introspection::get_molecule_ids(std::shared_ptr<Complex> pattern) {
  std::vector<int> res;

  // NOTE: some caching might be useful here if this function is called often
  uint_set<species_id_t> matching_species, not_matching_species;
  BNG::Cplx bng_pattern(&world->bng_engine.get_data());
  BNG::compartment_id_t primary_compartment_id = BNG::COMPARTMENT_ID_NONE;
  if (is_set(pattern)) {
    // convert to its BNG representation for matching
    BNGConverter converter(world->bng_engine.get_data(), world->bng_engine.get_config());
    bng_pattern = converter.convert_complex(*pattern, true);

    if (bng_pattern.has_compartment()) {
      // remove primary compartment were set, better explanation is in
      // MCell4Converter::convert_count_term_leaf_and_init_counting_flags
      primary_compartment_id = bng_pattern.get_primary_compartment_id();
      bng_pattern.remove_compartment_from_elem_mols(primary_compartment_id);
    }
  }

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  std::vector<MCell::Molecule>& molecules = p.get_molecules();
  for (MCell::Molecule& m: molecules) {
    if (m.is_defunct()) {
      continue;
    }
    if (is_set(pattern)) {
      // include only molecules that match the pattern
      if (matching_species.count(m.species_id) != 0) {
        res.push_back(m.id);
      }
      else if (not_matching_species.count(m.species_id) != 0) {
        // nothing to do
      }
      else {
        // not seen species
        const BNG::Species& species = world->get_all_species().get(m.species_id);
        BNG::compartment_id_t species_compartment = species.get_primary_compartment_id();
        bool match = species.matches_pattern(bng_pattern, true) &&
            (primary_compartment_id == BNG::COMPARTMENT_ID_NONE || primary_compartment_id == species_compartment);
        if (match) {
          matching_species.insert(m.species_id);
          res.push_back(m.id);
        }
        else {
          not_matching_species.insert(m.species_id);
        }
      }
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

  res = make_shared<API::Molecule>(DefaultCtorArgType());
  res->id = m.id;
  res->species_id = m.species_id;
  if (m.is_surf()) {
    const MCell::Wall& w = p.get_wall(m.s.wall_index);
    const Vec3& w_vert0 = p.get_wall_vertex(w, 0);

    res->type = MoleculeType::SURFACE;
    res->pos2d = Vec2(m.s.pos * Vec2(world->config.length_unit)).to_vec();
    res->pos3d = Vec3(GeometryUtils::uv2xyz(m.s.pos, w, w_vert0) * Vec3(world->config.length_unit)).to_vec();
    res->orientation = convert_mcell_orientation(m.s.orientation);

    res->geometry_object = model_inst->get_geometry_object_with_id(w.object_id);
    assert(is_set(res->geometry_object));
    res->wall_index = m.s.wall_index - res->geometry_object->first_wall_index;
  }
  else {
    res->type = MoleculeType::VOLUME;
    res->pos3d = Vec3(m.v.pos * Vec3(world->config.length_unit)).to_vec();
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


std::vector<double> Introspection::get_vertex(std::shared_ptr<GeometryObject> object, const int vertex_index) {
  object_is_set_and_initialized(object);

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  return
      Vec3(
          p.get_geometry_vertex(object->get_partition_vertex_index(vertex_index)) * Vec3(world->config.length_unit)
      ).to_vec();
}


std::shared_ptr<Wall> Introspection::get_wall(std::shared_ptr<GeometryObject> object, const int wall_index) {
  object_is_set_and_initialized(object);

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  const MCell::Wall& w = p.get_wall(object->get_partition_wall_index(wall_index));

  auto res = make_shared<Wall>(DefaultCtorArgType());
  res->geometry_object = object;
  res->wall_index = wall_index;
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    res->vertices.push_back(Vec3(p.get_geometry_vertex(w.vertex_indices[i]) * Vec3(world->config.length_unit)).to_vec());
  }
  res->area = w.area * world->config.length_unit * world->config.length_unit;
  assert(cmp_eq(len3(w.normal), (pos_t)1));
  res->unit_normal = w.normal.to_vec();
  res->is_movable = w.is_movable;
  res->world = world;
  return res;
}


std::vector<double> Introspection::get_vertex_unit_normal(std::shared_ptr<GeometryObject> object, const int vertex_index) {
  object_is_set_and_initialized(object);

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  // the walls that use the vertex, ascending, from the object's MeshTable
  const MCell::GeometryObject& obj = p.get_geometry_object_by_id(object->geometry_object_id);
  const MCell::MeshUtils::IndexRange walls = obj.get_walls_of_vertex(object->get_partition_vertex_index(vertex_index));

  if (walls.empty()) {
    throw RuntimeError("Internal error: there are no walls that use vertex with index " +
        to_string(vertex_index) + " of object " + object->name + ".");
  }

  Vec3 normals_sum = Vec3(0);
  for (const MCell::MeshUtils::mesh_index_t j: walls) {
    const MCell::Wall& w = p.get_wall(obj.wall_indices[j]);
    // wall normals are already unit vectors so we can just sum them
    normals_sum = normals_sum + w.normal;
  }

  return Vec3(normals_sum / Vec3(len3(normals_sum))).to_vec();
}


std::vector<double> Introspection::get_wall_unit_normal(std::shared_ptr<GeometryObject> object, const int wall_index) {
  object_is_set_and_initialized(object);

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  const MCell::Wall& w = p.get_wall(object->get_partition_wall_index(wall_index));

  // the value is is already normalized
  assert(cmp_eq(len3(w.normal), (pos_t)1));
  return w.normal.to_vec();
}


double Introspection::get_object_area(std::shared_ptr<GeometryObject> object) {
  object_is_set_and_initialized(object);

  // the area of the object's walls, kept with the checks of its mesh and
  // recomputed after its vertices move
  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  const MCell::GeometryObject& obj = p.get_geometry_object_by_id(object->geometry_object_id);
  const double length_unit = world->config.length_unit;
  return obj.get_area(p) * length_unit * length_unit;
}


double Introspection::get_object_volume(std::shared_ptr<GeometryObject> object) {
  object_is_set_and_initialized(object);

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  const MCell::GeometryObject& obj = p.get_geometry_object_by_id(object->geometry_object_id);

  // the checks stopped the run on a mesh that isn't manifold or consistently
  // oriented, so the watertight check ran; it depends only on the vertex
  // indices, so vertex moves don't change it
  const MCell::MeshUtils::MeshCheck& check = obj.get_mesh_check();
  release_assert(check.watertight_checked);
  if (!check.watertight()) {
    const size_t n = check.boundary_edges;
    throw ValueError(
        "Geometry object " + object->name + " is non-watertight (" + to_string(n) +
        (n == 1 ? " boundary edge" : " boundary edges") + "), so it has no volume; see " +
        world->config.get_geometry_report_file_name());
  }

  // the signed volume, summed over the object's shells, kept with the checks
  // and recomputed after its vertices move; computed only for a watertight
  // object
  const std::optional<double> signed_volume = obj.get_signed_volume(p);
  release_assert(signed_volume.has_value());
  const double length_unit = world->config.length_unit;
  return *signed_volume * length_unit * length_unit * length_unit;
}


std::shared_ptr<Color> Introspection::get_wall_color(std::shared_ptr<GeometryObject> object, const int wall_index) {
  object_is_set_and_initialized(object);

  wall_index_t wi = object->get_partition_wall_index(wall_index);
  const MCell::GeometryObject& obj = world->get_geometry_object(object->geometry_object_id);

  rgba_t color = obj.get_wall_color(wi);
  std::shared_ptr<Color> res = make_shared<Color>(DefaultCtorArgType());
  res->set_rgba(color); // must use setter - initializes also other attributes
  return res;
}


void Introspection::set_wall_color(
    std::shared_ptr<GeometryObject> object, const int wall_index, std::shared_ptr<Color> color) {
  object_is_set_and_initialized(object);

  object_is_set_and_initialized(object);

  wall_index_t wi = object->get_partition_wall_index(wall_index);
  MCell::GeometryObject& obj = world->get_geometry_object(object->geometry_object_id);

  obj.set_wall_color(wi, color->get_rgba());
}


} // namespace API
} // namespace MCell
