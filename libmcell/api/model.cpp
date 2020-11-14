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

#include "model.h"

#include <string>

#include "api/mcell4_converter.h"
#include "api/api_utils.h"
#include "api/molecule.h"

#include "api/species.h"
#include "api/reaction_rule.h"
#include "api/release_site.h"
#include "api/geometry_object.h"
#include "api/viz_output.h"
#include "api/count.h"
#include "api/wall.h"
#include "api/wall_wall_hit_info.h"

#include "world.h"


using namespace std;

namespace MCell {
namespace API {


Model::~Model() {
  delete world;
}


void Model::add_subsystem(std::shared_ptr<Subsystem> subsystem) {
  error_if_initialized(NAME_CLASS_SUBSYSTEM);

  append_vec_to_vec(elementary_molecule_types, subsystem->elementary_molecule_types);
  append_vec_to_vec(species, subsystem->species);
  append_vec_to_vec(surface_classes, subsystem->surface_classes);
  append_vec_to_vec(reaction_rules, subsystem->reaction_rules, false, true);
}


void Model::add_instantiation_data(std::shared_ptr<InstantiationData> instantiation_data) {
  error_if_initialized(NAME_CLASS_INSTANTIATION_DATA);

  append_vec_to_vec(release_sites, instantiation_data->release_sites);
  append_vec_to_vec(geometry_objects, instantiation_data->geometry_objects);
}


void Model::add_observables(std::shared_ptr<Observables> observables) {
  error_if_initialized(NAME_CLASS_OBSERVABLES);

  append_vec_to_vec(viz_outputs, observables->viz_outputs, true);
  append_vec_to_vec(counts, observables->counts, true);
}


void Model::initialize() {
  if (world != nullptr) {
    throw RuntimeError("Model.initialize() can be called only once");
  }

  world = new World(callbacks);

  // semantic checks are done during conversion
  MCell4Converter converter;

  converter.convert(this, world);

  // set that all used objects were initialized
  vec_set_initialized(species);
  vec_set_initialized(reaction_rules);
  vec_set_initialized(release_sites);
  vec_set_initialized(geometry_objects);
  vec_set_initialized(viz_outputs);
  vec_set_initialized(counts);

  world->init_simulation();

  initialized = true;
}


void Model::run_iterations(const float_t iterations) {
  if (world == nullptr) {
    throw RuntimeError("Model was not initialized, call Model.initialize() first");
  }
  uint output_frequency = World::determine_output_frequency(iterations);
  world->run_n_iterations(iterations, output_frequency, false);
}


void Model::end_simulation(const bool print_final_report) {
  // the first argument specifies that the last mol/rxn count and viz_output events will be run
  world->end_simulation(true, print_final_report);
}


void Model::dump_internal_state() {
  world->dump();
}


void Model::export_data_model_viz_or_full(
    const std::string& file,
    const bool only_for_visualization,
    const char* method_name) {

  if (is_set(file)) {
    world->export_data_model(file, only_for_visualization);
  }
  else {
    world->export_data_model_to_dir(get_first_viz_output_files_prefix(method_name));
  }
}


std::vector<int> Model::get_molecule_ids(std::shared_ptr<Species> species) {
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


std::shared_ptr<API::Molecule> Model::get_molecule(const int id) {
  std::shared_ptr<API::Molecule> res;
  Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  if (!p.does_molecule_exist(id)) {
    throw RuntimeError("Molecule with id " + to_string(id) + " does not exist.");
  }
  MCell::Molecule& m = p.get_m(id);
  if (m.is_defunct()) {
    throw RuntimeError("Molecule with id " + to_string(id) + " was removed.");
  }

  res = make_shared<API::Molecule>();
  res->id = m.id;
  if (m.is_surf()) {
    // TODO: res->pos3d
    res->orientation = convert_orientation(m.s.orientation);
  }
  else {
    res->pos3d = m.v.pos * Vec3(world->config.length_unit);
    res->orientation = Orientation::NONE;
  }
  res->world = world;
  res->set_initialized();

  return res;
}


Vec3 Model::get_vertex(std::shared_ptr<GeometryObject> object, const int vertex_index) {
  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  return
      p.get_geometry_vertex(object->get_partition_vertex_index(vertex_index)) *
      Vec3(world->config.length_unit);
}


std::shared_ptr<Wall> Model::get_wall(std::shared_ptr<GeometryObject> object, const int wall_index) {
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
  res->normal = w.normal; // no need to convert units here
  res->is_movable = w.is_movable;
  res->world = world;
  return res;
}


Vec3 Model::get_vertex_unit_normal(std::shared_ptr<GeometryObject> object, const int vertex_index) {
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
    normals_sum = normals_sum + w.normal;
  }

  return normals_sum / Vec3(len3(normals_sum));
}


Vec3 Model::get_wall_unit_normal(std::shared_ptr<GeometryObject> object, const int wall_index) {
  object->check_is_initialized();

  const MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  const MCell::Wall& w = p.get_wall(object->get_partition_wall_index(wall_index));

  return w.normal / Vec3(len3(w.normal));
}


void Model::add_vertex_move(
    std::shared_ptr<GeometryObject> object, const int vertex_index, const Vec3& displacement
) {
  // - currently, it is not expected that the user will have access to the scheduled vertex moves
  // - later we can use the object id to determine the partition

  object->check_is_initialized();

  release_assert(
      object->first_vertex_index + vertex_index <
      world->get_partition(PARTITION_ID_INITIAL).get_geometry_vertex_count()
  );

  vertex_moves.push_back(
      VertexMoveInfo(
          PARTITION_ID_INITIAL,
          object->geometry_object_id,
          object->get_partition_vertex_index(vertex_index),
          displacement * Vec3(world->config.rcp_length_unit) // convert to internal units
      )
  );
}


std::vector<std::shared_ptr<WallWallHitInfo>> Model::apply_vertex_moves(
    const bool collect_wall_wall_hits) {

  // run the actual vertex update
  std::set<GeometryObjectWallUnorderedPair> colliding_walls;
  world->get_partition(PARTITION_ID_INITIAL).apply_vertex_moves(vertex_moves, colliding_walls);
  vertex_moves.clear();

  std::vector<std::shared_ptr<WallWallHitInfo>> res;

  // convert information on processed hits
  if (collect_wall_wall_hits) {
    for (const auto& wall_pair: colliding_walls) {
      auto obj1 = get_geometry_object_with_id(wall_pair.geometry_object_id1);
      int obj_wall_index1 = obj1->get_object_wall_index(wall_pair.wall_index1);
      auto obj2 = get_geometry_object_with_id(wall_pair.geometry_object_id2);
      int obj_wall_index2 = obj2->get_object_wall_index(wall_pair.wall_index2);

      auto pair = make_shared<WallWallHitInfo>();
      pair->wall1 = get_wall(obj1, obj_wall_index1);
      pair->wall2 = get_wall(obj2, obj_wall_index2);

      res.push_back(pair);
    }
  }

  return res;
}


void Model::register_mol_wall_hit_callback(
    const std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)> function,
    py::object context,
    std::shared_ptr<GeometryObject> object,
    std::shared_ptr<Species> species
) {
  if (!initialized) {
    throw RuntimeError("Model must be initialized before registering callbacks");
  }

  geometry_object_id_t geometry_object_id = GEOMETRY_OBJECT_ID_INVALID;
  if (is_set(object)) {
    if (object->geometry_object_id == GEOMETRY_OBJECT_ID_INVALID) {
      throw RuntimeError("Geometry object " + object->name + " is not present in model.");
    }
    geometry_object_id = object->geometry_object_id;
  }

  species_id_t species_id = SPECIES_ID_INVALID;
  if (is_set(species)) {
    if (species->species_id == SPECIES_ID_INVALID) {
      throw RuntimeError("Species object " + object->name + " is not present in model.");
    }
    species_id = species->species_id;
  }

  callbacks.register_mol_wall_hit_callback(function, context, geometry_object_id, species_id);
}


void Model::load_bngl(
    const std::string& file_name,
    const std::string& observables_files_prefix,
    std::shared_ptr<Region> default_release_region,
    const std::map<std::string, float_t>& parameter_overrides) {

  BNG::BNGData bng_data;

  int num_errors = BNG::parse_bngl_file(file_name, bng_data, parameter_overrides);
  if (num_errors != 0) {
    throw RuntimeError("Could not parse BNGL file " + file_name + ".");
  }

  // now convert everything we parsed into the API classes so that the user can
  // inspect or manipulate it if needed
  convert_bng_data_to_subsystem_data(bng_data);

  // needs subsystem data created in the last step
  convert_bng_data_to_instantiation_data(bng_data, *this, default_release_region);

  convert_bng_data_to_observables_data(bng_data, *this, observables_files_prefix);
}


void Model::export_to_bngl(const std::string& file_name) {
  if (!initialized) {
    throw RuntimeError("Model must be initialized for BNGL export.");
  }

  string err_msg = world->export_to_bngl(file_name);
  if (err_msg != "") {
    throw RuntimeError("BNGL export failed: " + err_msg);
  }
}


// overrides from derived classes Subsystem, InstantiationData, and Observables,
// in .cpp because implementation in .h file would need too many headers to be included
void Model::add_species(std::shared_ptr<Species> s) {
  error_if_initialized(NAME_CLASS_SPECIES);
  Subsystem::add_species(s);
}

void Model::add_reaction_rule(std::shared_ptr<ReactionRule> r) {
  error_if_initialized(NAME_CLASS_REACTION_RULE);
  Subsystem::add_reaction_rule(r);
}

void Model::add_surface_class(std::shared_ptr<SurfaceClass> sc) {
  error_if_initialized(NAME_CLASS_SURFACE_CLASS);
  Subsystem::add_surface_class(sc);
}

void Model::add_release_site(std::shared_ptr<ReleaseSite> s) {
  error_if_initialized(NAME_CLASS_RELEASE_SITE);
  InstantiationData::add_release_site(s);
}

void Model::add_geometry_object(std::shared_ptr<GeometryObject> o) {
  error_if_initialized(NAME_CLASS_GEOMETRY_OBJECT);
  InstantiationData::add_geometry_object(o);
}

void Model::add_viz_output(std::shared_ptr<VizOutput> viz_output) {
  error_if_initialized(NAME_CLASS_VIZ_OUTPUT);
  Observables::add_viz_output(viz_output);
};

void Model::add_count(std::shared_ptr<Count> count) {
  error_if_initialized(NAME_CLASS_OBSERVABLES);
  Observables::add_count(count);
};

std::string Model::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Model" << ": " <<
      "\n" << ind + "  " <<
      //"config=" << "(" << ((config != nullptr) ? config->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      //"warnings=" << "(" << ((warnings != nullptr) ? warnings->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      //"notifications=" << "(" << ((notifications != nullptr) ? notifications->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "elementary_molecule_types=" << vec_ptr_to_str(elementary_molecule_types, ind + "  ") << ", " << "\n" << ind + "  " <<
      "species=" << vec_ptr_to_str(species, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_classes=" << vec_ptr_to_str(surface_classes, ind + "  ") << ", " << "\n" << ind + "  " <<
      "reaction_rules=" << vec_ptr_to_str(reaction_rules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "release_sites=" << vec_ptr_to_str(release_sites, ind + "  ") << ", " << "\n" << ind + "  " <<
      "geometry_objects=" << vec_ptr_to_str(geometry_objects, ind + "  ");
  return ss.str();
}

std::shared_ptr<GeometryObject> Model::get_geometry_object_with_id(const geometry_object_id_t id) {
  // not very efficient, we may need some caching/map later
  for (auto o: geometry_objects) {
    if (o->geometry_object_id == id) {
      return o;
    }
  }
  return std::shared_ptr<GeometryObject>(nullptr);
}


void Model::dump() const {
  cout << to_str();
}

}
}

