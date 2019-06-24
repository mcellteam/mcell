/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#ifndef SRC4_PARTITION_H_
#define SRC4_PARTITION_H_

#include <set>
#include <boost/container/flat_set.hpp>

#include "defines.h"
#include "molecule.h"
#include "scheduler.h"
#include "reaction.h"
#include "geometry.h"

namespace mcell {

/**
 * Class used to hold sets of ids or indices of molecules or other items
 */
// TODO: rename to subpart_mask_t
class subpartition_mask_t: public boost::container::flat_set<uint32_t> {
public:
  void set_contains_id(const uint32_t id, const bool value = true) {
    if (value) {
      assert(count(id) == 0);
      insert(id);
    }
    else {
      assert(count(id) == 1);
      erase(id);
    }
  }

  void dump();
};


/**
 * Parition class contains all molecules and other data contained in
 * one simulation block.
 */
class partition_t {
public:
  partition_t(
      const vec3_t origin_,
      const world_constants_t& world_constants_,
      simulation_stats_t& simulation_stats_
  )
    : origin_corner(origin_),
      next_molecule_id(0),
      world_constants(world_constants_),
      simulation_stats(simulation_stats_) {

    opposite_corner = origin_corner + world_constants.partition_edge_length;

    // preaallocate volume_molecules arrays and also volume_molecule_indices_per_time_step
    uint32_t num_subparts = powu(world_constants.subpartitions_per_partition_dimension, 3);
    volume_molecule_reactants_per_subpart.resize(num_subparts);
    walls_per_subpart.resize(num_subparts);

    size_t num_species = world_constants.bimolecular_reactions_map->size();
    for (auto& reactants : volume_molecule_reactants_per_subpart) {
      reactants.resize(num_species);
    }
  }


  molecule_t& get_m(const molecule_id_t idx) { // should be ID
    assert(idx < molecule_id_to_index_mapping.size());

    // code works with molecule ids, but they need to be converted to indices to the volume_molecules vector
    // because we need to defragment the contents
    uint32_t vm_vec_index = molecule_id_to_index_mapping[idx];
    assert(vm_vec_index != MOLECULE_INDEX_INVALID);
    return molecules[vm_vec_index];
  }


  molecule_id_t get_molecule_index(const molecule_t& m) {
    // simply use pointer arithmetic to compute the molecule's index
    molecule_id_t res = m.id;
    assert(res != MOLECULE_ID_INVALID);
    return res;
  }


  uint32_t get_molecule_list_index_for_time_step(const float_t time_step) {
    for (uint32_t i = 0; i < molecules_data_per_time_step_array.size(); i++) {
      if (molecules_data_per_time_step_array[i].time_step == time_step) {
        return i;
      }
    }
    return TIME_STEP_INDEX_INVALID;
  }


  uint32_t get_or_add_molecule_list_index_for_time_step(const float_t time_step) {
    uint32_t res;
    res = get_molecule_list_index_for_time_step(time_step);
    if (res == TIME_STEP_INDEX_INVALID) {
      molecules_data_per_time_step_array.push_back(
        time_step_molecules_data_t(time_step, std::vector<molecule_id_t>()));
      res = molecules_data_per_time_step_array.size() - 1;
    }
    return res;
  }


  bool in_this_partition(const vec3_t& pos) const {
    return glm::all(glm::greaterThanEqual(pos, origin_corner))
      && glm::all(glm::lessThan(pos, opposite_corner));
  }


  void get_subpart_3d_indices(const vec3_t& pos, ivec3_t& res) const {
    assert(in_this_partition(pos));
    vec3_t relative_position = pos - origin_corner;
    res = relative_position * world_constants.subpartition_edge_length_rcp;
  }


  uint32_t get_subpartition_index_from_3d_indices(const ivec3_t& indices) const {
    // example: dim: 5x5x5,  (1, 2, 3) -> 1 + 2*5 + 3*5*5 = 86
    return
        indices.x +
        indices.y * world_constants.subpartitions_per_partition_dimension +
        indices.z * world_constants.subpartitions_per_partition_dimension_squared;
  }


  uint32_t get_subpart_index_from_3d_indices(const int x, const int y, const int z) const {
    return get_subpartition_index_from_3d_indices(ivec3_t(x, y, z));
  }


  void get_subpart_3d_indices_from_index(const uint32_t index, ivec3_t& indices) const {
    uint32_t dim = world_constants.subpartitions_per_partition_dimension;
    // example: dim: 5x5x5,  86 -> (86%5, (86/5)%5, (86/(5*5))%5) = (1, 2, 3)
    indices.x = index % dim;
    indices.y = (index / dim) % dim;
    indices.z = (index / world_constants.subpartitions_per_partition_dimension_squared) % dim;
  }


  uint32_t get_subpartition_index(const vec3_t& pos) {
    ivec3_t indices;
    get_subpart_3d_indices(pos, indices);
    return get_subpartition_index_from_3d_indices(indices);
  }


  void get_subpart_llf_point(const subpart_index_t subpart_index, vec3_t& llf) const {
    ivec3_t indices;
    get_subpart_3d_indices_from_index(subpart_index, indices);
    llf = origin_corner + vec3_t(indices) * vec3_t(world_constants.subpartition_edge_length);
  }

  void get_subpart_urb_point_from_llf(const vec3_t& llf, vec3_t& urb) const {
    urb = llf + vec3_t(world_constants.subpartition_edge_length);
  }

  void change_reactants_map(molecule_t& vm, const uint32_t new_subpartition_index, bool adding, bool removing) {
    assert(world_constants.bimolecular_reactions_map->find(vm.species_id) != world_constants.bimolecular_reactions_map->end());

    // these are all the sets of indices of reactants for this particular subpartition
    species_reactants_map_t& subpart_reactants_orig_sp = volume_molecule_reactants_per_subpart[vm.v.subpart_index];
    species_reactants_map_t& subpart_reactants_new_sp = volume_molecule_reactants_per_subpart[new_subpartition_index];

    // and these are indices of possible reactants with our reactant_species_id
    const species_reaction_map_t& reactions_map = world_constants.bimolecular_reactions_map->find(vm.species_id)->second;

    // we need to set/clear flag that says that second_reactant_info.first can react with reactant_species_id
    for (const auto& second_reactant_info : reactions_map) {
      if (removing) {
        subpart_reactants_orig_sp[second_reactant_info.first].set_contains_id(vm.id, false);
      }
      if (adding) {
        subpart_reactants_new_sp[second_reactant_info.first].set_contains_id(vm.id, true);
      }
    }
  }


  void change_molecule_subpartition(molecule_t& vm, const uint32_t new_subpartition_index) {
    assert(vm.v.subpart_index < volume_molecule_reactants_per_subpart.size());
    assert(new_subpartition_index < volume_molecule_reactants_per_subpart.size());
    if (vm.v.subpart_index == new_subpartition_index) {
      return; // nothing to do
    }
#ifdef DEBUG_SUBPARTITIONS
    std::cout << "Molecule " << molecule_idx << " changed subpartition from "
        <<  vm.v.subpart_index << " to " << new_subpartition_index << ".\n";
#endif

    change_reactants_map(vm, new_subpartition_index, true, true);
    vm.v.subpart_index = new_subpartition_index;
  }


  // version that computes the right time_step_index each time it is called
  molecule_t& add_volume_molecule(const molecule_t& vm_copy, const float_t time_step) {
    uint32_t time_step_index = get_or_add_molecule_list_index_for_time_step(time_step);
    return add_volume_molecule_with_time_step_index(vm_copy, time_step_index);
  }


  // any molecule flags are set by caller after the molecule is created by this method
  molecule_t& add_volume_molecule_with_time_step_index(molecule_t vm_copy, const uint32_t time_step_index) {
    molecule_id_t molecule_id = next_molecule_id;
    next_molecule_id++;
    // and its index to the list sorted by time step
    // this is an array that changes only when molecule leaves this partition
    assert(time_step_index <= molecules_data_per_time_step_array.size());
    molecules_data_per_time_step_array[time_step_index].molecule_ids.push_back(molecule_id);

    // We always have to increase the size of the mappping array - its size is
    // large enough to hold indices for all molecules that were ever created,
    // we will need to reuse ids or compress it later
    uint32_t next_molecule_array_index = molecules.size(); // get the index of the molecule we aregoing to store
    molecule_id_to_index_mapping.push_back(next_molecule_array_index);
    assert(
        molecule_id_to_index_mapping.size() == next_molecule_id
        && "Mapping array must have value for every molecule index"
    );

    // This is the only place where we insert molecules into volume_molecules,
    // although this array size can be decreased in defragmentation
    molecules.push_back(vm_copy);
    molecule_t& new_vm = molecules.back();

    new_vm.id = molecule_id;
    new_vm.v.subpart_index = get_subpartition_index(vm_copy.v.pos);
    change_reactants_map(new_vm, new_vm.v.subpart_index, true, false);

    return new_vm;
  }


  molecule_t& add_surface_molecule(const molecule_t& sm_copy, const vec3_t& pos3d, const float_t time_step) {

    uint32_t time_step_index = get_or_add_molecule_list_index_for_time_step(time_step);

    molecule_id_t molecule_id = next_molecule_id;
    next_molecule_id++;

    // TODO: merge with surface molecules
    // and its index to the list sorted by time step
    // this is an array that changes only when molecule leaves this partition
    assert(time_step_index <= molecules_data_per_time_step_array.size());
    molecules_data_per_time_step_array[time_step_index].molecule_ids.push_back(molecule_id);

    // We always have to increase the size of the mapping array - its size is
    // large enough to hold indices for all molecules that were ever created,
    // we will need to reuse ids or compress it later
    uint32_t next_molecule_array_index = molecules.size(); // get the index of the molecule we aregoing to store
    molecule_id_to_index_mapping.push_back(next_molecule_array_index);
    assert(
        molecule_id_to_index_mapping.size() == next_molecule_id
        && "Mapping array must have value for every molecule index"
    );

    // This is the only place where we insert molecules into volume_molecules,
    // although this array size can be decreased in defragmentation
    molecules.push_back(sm_copy);
    molecule_t& new_sm = molecules.back();

    new_sm.id = molecule_id;
    // new_sm.s.subpart_index = get_subpartition_index(pos3d);
    //    change_reactants_map(new_vm, new_sm.s.subpart_index, true, false);

    return new_sm;
  }



  void set_molecule_as_defunct(molecule_t& vm) {
    // set that this molecule does not exist anymore
    vm.set_is_defunct();

    change_reactants_map(vm, 0/*unused*/, false, true);
  }


  void add_unimolecular_action(float_t diffusion_time_step, const diffuse_or_unimol_react_action_t& unimol_react_action) {
    uint32_t time_step_index = get_or_add_molecule_list_index_for_time_step(diffusion_time_step);
    molecules_data_per_time_step_array[time_step_index].calendar_for_unimol_rxs.insert(unimol_react_action, unimol_react_action.scheduled_time);
  }

  // ---------------------------------- typedefs and internal structs ----------------------------------


  // calendar type for internal action scheduling, one bucket corresponds to a single diffusion time step (diffuse_and_react_event)
  typedef calendar_t<fifo_bucket_t<diffuse_or_unimol_react_action_t>, diffuse_or_unimol_react_action_t> calendar_for_unimol_rxs_t;
  
  // this structure contains all data associated with a given diffusion time step.
  struct time_step_molecules_data_t {
    // usually initialized with empty molecule_ids_ array
    time_step_molecules_data_t(float_t time_step_, const std::vector< molecule_id_t > molecule_ids_)
      : time_step(time_step_), molecule_ids(molecule_ids_), calendar_for_unimol_rxs(time_step_) {
    }

		// diffusion time step value 
    float_t time_step;
    // molecule ids with this diffusion time step 
    std::vector<molecule_id_t> molecule_ids;
    // and unimolecular reactions scheduled while diffusing  molecules with this diffusion time step;
    // these unimolecular reactions can be associated with a molecule with a different time step, 
    // what is important is how they were created 
    calendar_for_unimol_rxs_t calendar_for_unimol_rxs;
  };

  // indexed with species_id_t
  typedef std::vector< subpartition_mask_t > species_reactants_map_t;


  // ---------------------------------- molecule getters ----------------------------------

  // usually used as constant
  std::vector< time_step_molecules_data_t >& get_molecule_data_per_time_step_array() {
    return molecules_data_per_time_step_array;
  }
  

  const std::vector<molecule_id_t>& get_molecule_ids_for_time_step_index(uint32_t time_step_index) {
    assert(time_step_index < molecules_data_per_time_step_array.size());
    return molecules_data_per_time_step_array[time_step_index].molecule_ids;
  }
  

  // calendar is usually modified by the caller
  calendar_for_unimol_rxs_t& get_unimolecular_actions_calendar_for_time_step_index(uint32_t time_step_index) {
    assert(time_step_index < molecules_data_per_time_step_array.size());
    return molecules_data_per_time_step_array[time_step_index].calendar_for_unimol_rxs;
  }
  

  const vec3_t& get_origin_corner() const {
    return origin_corner;
  }
  
  const vec3_t& get_opposite_corner() const {
    return opposite_corner;
  }

  subpartition_mask_t& get_volume_molecule_reactants(subpart_index_t subpart_index, species_id_t species_id) {
    return volume_molecule_reactants_per_subpart[subpart_index][species_id];
  }


  const std::vector<molecule_t>& get_molecules() const {
    return molecules;
  }


  std::vector<molecule_t>& get_molecules() {
    return molecules;
  }
  

  std::vector<uint32_t>& get_molecule_id_to_index_mapping() {
    return molecule_id_to_index_mapping;
  }
  
  // ---------------------------------- geometry ----------------------------------
  vertex_index_t add_geometry_vertex(const vec3_t pos) {
    vertex_index_t index = geometry_vertices.size();
    geometry_vertices.push_back(pos);
    return index;
  }

  const vec3_t get_geometry_vertex(vertex_index_t i) const {
    assert(i < geometry_vertices.size());
    return geometry_vertices[i];
  }

  // returns reference to the new wall, only sets id and index
  // register_wall must be called
  wall_t& add_uninitialized_wall(const wall_id_t id) {
    wall_index_t index = walls.size();
    walls.push_back(wall_t());
    wall_t& new_wall = walls.back();

    new_wall.id = id;
    new_wall.index = index;

    return new_wall;
  }

  // when a wall is  added with add_uninitialized_wall,
  // its type and vertices are not know yet, we must include the walls
  // into subvolumes
  // TODO: some better name?, do this per object?
  // does not have to be in .h
  void register_wall_into_supartitions(const wall_index_t wall_index) {
    wall_t& w = get_wall(wall_index);

    // also insert this triangle into walls per subparition
    subpart_indices_vector_t colliding_subparts;
    geometry::wall_subparts_collision_test(*this, w, colliding_subparts);
    for (subpart_index_t subpart_index: colliding_subparts) {
      assert(subpart_index < walls_per_subpart.size());
      walls_per_subpart[subpart_index].set_contains_id(wall_index);
    }
  }

  // returns reference to the new object, only sets id
  geometry_object_t& add_uninitialized_geometry_object(const geometry_object_id_t id) {
    geometry_object_index_t index = geometry_objects.size();
    geometry_objects.push_back(geometry_object_t());
    geometry_object_t& new_obj = geometry_objects.back();

    new_obj.id = id;

    return new_obj;
  }

  const wall_t& get_wall(wall_index_t i) const {
    assert(i < walls.size());
    return walls[i];
  }

  wall_t& get_wall(wall_index_t i) {
    assert(i < walls.size());
    return walls[i];
  }

  // maybe we will need to filter out, e.g. just reflective surfaces
  const subpartition_mask_t& get_subpart_wall_indices(subpart_index_t subpart_index) const {
    return walls_per_subpart[subpart_index];
  }

  // used to create uninitialized grid object that can be initialized with
  // grid_t::initialize, sets wall's index
  // TODO: check all this geometry creation and come up with some more systematic way
  grid_t& add_grid_for_wall(wall_t& w) {
    // TODO: check that this wall belongs to this partition
    w.grid_index = grids.size();
    grids.push_back(grid_t());
    grids.back().initialize(w.area);
    return grids.back();
  }

  grid_t& get_grid(const grid_index_t grid_index) {
    assert(grid_index < grids.size());
    return grids[grid_index];
  }

  // ---------------------------------- other ----------------------------------

  const world_constants_t& get_world_constants() const {
    return world_constants;
  }

  simulation_stats_t& get_simulation_stats() const {
    return simulation_stats;
  }

  void dump();

private:
  // left, bottom, closest (lowest z) point of the partition
  vec3_t origin_corner;
  vec3_t opposite_corner;

  // ---------------------------------- molecules ----------------------------------

  // vector containing all molecules in this partition (volume and surface)
  std::vector<molecule_t> molecules;

  // contains mapping of molecule ids to indices to the molecules array
  std::vector<uint32_t> molecule_id_to_index_mapping;

  // id of the next molecule to be created
  molecule_id_t next_molecule_id;

  // indexed by diffusion time step index
  std::vector<time_step_molecules_data_t> molecules_data_per_time_step_array;

  // indexed with subpartition index, only for vol-vol reactions
  std::vector<species_reactants_map_t> volume_molecule_reactants_per_subpart;

  // ---------------------------------- geometry objects ----------------------------------

  std::vector<vec3_t> geometry_vertices;

  // we must plan for dynamic geometry but for now its just static
  std::vector<geometry_object_t> geometry_objects;

  std::vector<wall_t> walls;
  std::vector<grid_t> grids; // not every wall has a grid

  std::vector<subpartition_mask_t> walls_per_subpart;

  const world_constants_t& world_constants; // owned by world
  simulation_stats_t& simulation_stats; // owned by world
};

} // namespace mcell

#endif // SRC4_PARTITION_H_
