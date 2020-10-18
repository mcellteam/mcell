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

  world = new World();

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
    res->pos3d = m.v.pos;
    res->orientation = Orientation::NONE;
  }
  res->world = world;
  res->set_initialized();

  return res;
}


void Model::add_vertex_move(
    std::shared_ptr<GeometryObject> object, const int index, const Vec3& displacement
) {
  // currently, it is not expected that the user will have access to the scheduled vertex moves
  if (object->geometry_object_id == GEOMETRY_OBJECT_ID_INVALID) {
    throw RuntimeError("Geometry object " + object->name + " is not present in model (or model was not initialized).");
  }

  if (index < 0 || index >= (int)object->vertex_list.size()) {
    throw RuntimeError(
        "Vertex index " + to_string(index) + " is out of range for " + NAME_VERTEX_LIST + " of " + object->name + ".");
  }

  // later we can use the object id to determine the partition
  release_assert(
      object->first_vertex_index + index <
      world->get_partition(PARTITION_ID_INITIAL).get_geometry_vertex_count()
  );

  vertex_moves.push_back(
      VertexMoveInfo(
          PARTITION_ID_INITIAL,
          object->first_vertex_index + index,
          displacement * Vec3(world->config.rcp_length_unit) // convert units
      )
  );
}


void Model::apply_vertex_moves() {
  // run the actual vertex update
  world->get_partition(PARTITION_ID_INITIAL).apply_vertex_moves(vertex_moves);
  vertex_moves.clear();
}


void Model::register_wall_hit_callback(
    const std::function<void(std::shared_ptr<WallHitInfo>, py::object)> function,
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

  world->register_wall_hit_callback(function, context, geometry_object_id, species_id);
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


void Model::dump() const {
  cout << to_str();
}

}
}

