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

#include "../libs/bng/rxn_container.h"
#include "defines.h"
#include "dyn_vertex_structs.h"
#include "molecule.h"
#include "scheduler.h"
#include "geometry.h"

namespace MCell {

// class used to hold potential reactants of given species
// performance critical, therefore we are using a vector for now,
// will probably need to change in the future
class SpeciesReactantsMap:
    // public std::map<species_id_t, uint_set<molecule_id_t> > {
    public std::vector<uint_set<molecule_id_t> > {

public:
  // key must exist
  void erase_existing(const species_id_t key, const molecule_id_t id) {
    //auto it = find(key);
    //assert(it != end() && "Key must exist");
    //it->second.erase_existing(id);

    assert(key < this->size());
    (*this)[key].erase_existing(id);
  }

  // key is created if it does not exist
  void insert_unique(const species_id_t key, const molecule_id_t id) {
    if (key >= this->size()) {
      // resize vector
      this->resize(key + 1);
    }
    (*this)[key].insert_unique(id);
  }

  const uint_set<molecule_id_t>& get_set(const species_id_t key) {
    if (key >= this->size()) {
      // resize vector
      return empty_set;
    }
    else {
      return (*this)[key];
    }
  }

private:
  uint_set<molecule_id_t> empty_set;
};


/**
 * Partition class contains all molecules and other data contained in
 * one simulation block.
 */
class Partition {
public:
  Partition(
      const partition_id_t id_,
      const Vec3 origin_,
      const SimulationConfig& config_,
      BNG::BNGEngine& bng_engine_,
      SimulationStats& stats_
  )
    : origin_corner(origin_),
      next_molecule_id(0),
      id(id_),
      config(config_),
      bng_engine(bng_engine_),
      stats(stats_) {

    opposite_corner = origin_corner + config.partition_edge_length;

    // pre-allocate volume_molecules arrays and also volume_molecule_indices_per_time_step
    uint32_t num_subparts = powu(config.subpartitions_per_partition_dimension, 3);
    volume_molecule_reactants_per_subpart.resize(num_subparts);
    walls_per_subpart.resize(num_subparts);
  }


  Molecule& get_m(const molecule_id_t id) {
    assert(id != MOLECULE_ID_INVALID);
    assert(id < molecule_id_to_index_mapping.size());

    // code works with molecule ids, but they need to be converted to indices to the volume_molecules vector
    // because we need to defragment the contents
    uint32_t vm_vec_index = molecule_id_to_index_mapping[id];
    assert(vm_vec_index != MOLECULE_INDEX_INVALID);
    return molecules[vm_vec_index];
  }

  const Molecule& get_m(const molecule_id_t id) const {
    assert(id != MOLECULE_ID_INVALID);
    assert(id < molecule_id_to_index_mapping.size());

    // code works with molecule ids, but they need to be converted to indices to the volume_molecules vector
    // because we need to defragment the contents
    uint32_t vm_vec_index = molecule_id_to_index_mapping[id];
    assert(vm_vec_index != MOLECULE_INDEX_INVALID);
    return molecules[vm_vec_index];
  }

  molecule_id_t get_molecule_index(const Molecule& m) {
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
        TimeStepMoleculesData(time_step, std::vector<molecule_id_t>()));
      res = molecules_data_per_time_step_array.size() - 1;
    }
    return res;
  }


  bool in_this_partition(const Vec3& pos) const {
    return glm::all(glm::greaterThanEqual(pos, origin_corner))
      && glm::all(glm::lessThan(pos, opposite_corner));
  }


  bool is_subpart_index_in_range(const int index) const {
    assert((subpart_index_t)index != SUBPART_INDEX_INVALID);
    return index >= 0 && index < (int)config.subpartitions_per_partition_dimension;
  }


  void get_subpart_3d_indices(const Vec3& pos, IVec3& res) const {
    assert(in_this_partition(pos));
    Vec3 relative_position = pos - origin_corner;
    res = relative_position * config.subpartition_edge_length_rcp;
  }

  // FIXME: use subpart_index_t, rename to subpart
  subpart_index_t get_subpart_index_from_3d_indices(const IVec3& indices) const {
    // example: dim: 5x5x5,  (1, 2, 3) -> 1 + 2*5 + 3*5*5 = 86
    return
        indices.x +
        indices.y * config.subpartitions_per_partition_dimension +
        indices.z * config.subpartitions_per_partition_dimension_squared;
  }

  subpart_index_t get_subpart_index_from_3d_indices(const int x, const int y, const int z) const {
    return get_subpart_index_from_3d_indices(IVec3(x, y, z));
  }


  void get_subpart_3d_indices_from_index(const uint32_t index, IVec3& indices) const {
    uint32_t dim = config.subpartitions_per_partition_dimension;
    // example: dim: 5x5x5,  86 -> (86%5, (86/5)%5, (86/(5*5))%5) = (1, 2, 3)
    indices.x = index % dim;
    indices.y = (index / dim) % dim;
    indices.z = (index / config.subpartitions_per_partition_dimension_squared) % dim;
  }

  subpart_index_t get_subpart_index(const Vec3& pos) const {
    IVec3 indices;
    get_subpart_3d_indices(pos, indices);
    return get_subpart_index_from_3d_indices(indices);
  }


  void get_subpart_llf_point(const subpart_index_t subpart_index, Vec3& llf) const {
    IVec3 indices;
    get_subpart_3d_indices_from_index(subpart_index, indices);
    llf = origin_corner + Vec3(indices) * Vec3(config.subpartition_edge_length);
  }

  void get_subpart_urb_point_from_llf(const Vec3& llf, Vec3& urb) const {
    urb = llf + Vec3(config.subpartition_edge_length);
  }

  void change_reactants_map(Molecule& vm, const uint32_t new_subpartition_index, bool adding, bool removing) {
    if (vm.is_surf()) {
      // nothing to do
      return;
    }

    assert(vm.is_vol() && "This function is applicable only to volume mols and ignored for surface mols");

    // and these are indices of possible reactants with our reactant_species_id
    // NOTE: this must be fast, bng engine must have this map/vector already ready
    // TODO: we can optimize this by taking just volume reactants into account
    const BNG::SpeciesRxnClassesMap* reactions_map = get_all_rxns().get_bimol_rxns_for_reactant(vm.species_id);
    if (reactions_map == nullptr) {
      // nothing to do
      return;
    }

    // these are all the sets of indices of reactants for this particular subpartition
    SpeciesReactantsMap& subpart_reactants_orig_sp = volume_molecule_reactants_per_subpart[vm.v.subpart_index];
    SpeciesReactantsMap& subpart_reactants_new_sp = volume_molecule_reactants_per_subpart[new_subpartition_index];

    // we need to set/clear flag that says that second_reactant_info.first can react with reactant_species_id
    for (const auto& second_reactant_info: *reactions_map) {
      if (second_reactant_info.second->get_num_reactions() == 0) {
        // there is a reaction class, but it has no reactions
        continue;
      }

      if (removing) {
        subpart_reactants_orig_sp.erase_existing(second_reactant_info.first, vm.id);
      }
      if (adding) {
        subpart_reactants_new_sp.insert_unique(second_reactant_info.first, vm.id);
      }
    }
  }


  void change_molecule_subpartition(Molecule& vm, const uint32_t new_subpartition_index) {
    assert(vm.v.subpart_index < volume_molecule_reactants_per_subpart.size());
    assert(new_subpartition_index < volume_molecule_reactants_per_subpart.size());
    if (vm.v.subpart_index == new_subpartition_index) {
      return; // nothing to do
    }
#ifdef DEBUG_SUBPARTITIONS
    std::cout << "Molecule " << vm.id << " changed subpartition from "
        <<  vm.v.subpart_index << " to " << new_subpartition_index << ".\n";
#endif

    change_reactants_map(vm, new_subpartition_index, true, true);
    vm.v.subpart_index = new_subpartition_index;
  }



private:
  // internal methods that sets molecule's id and
  // adds it to all relevant structures

  void add_molecule_to_diffusion_list(const Molecule& m, const uint32_t time_step_index) {

    // and its index to the list sorted by time step
    // this is an array that changes only when molecule leaves this partition or is defunct
    assert(time_step_index <= molecules_data_per_time_step_array.size());
    molecules_data_per_time_step_array[time_step_index].molecule_ids.push_back(m.id);
  }


  Molecule& add_molecule(const Molecule& vm_copy, const bool is_vol) {

    const BNG::Species& species = get_all_species().get(vm_copy.species_id);
    assert((is_vol && species.is_vol()) || (!is_vol && species.is_surf()));
    uint32_t time_step_index = get_or_add_molecule_list_index_for_time_step(species.time_step);

    molecule_id_t molecule_id = next_molecule_id;
    next_molecule_id++;

    // We always have to increase the size of the mapping array - its size is
    // large enough to hold indices for all molecules that were ever created,
    // we will need to reuse ids or compress it later
    uint32_t next_molecule_array_index = molecules.size(); // get the index of the molecule we are going to store
    molecule_id_to_index_mapping.push_back(next_molecule_array_index);
    assert(
        molecule_id_to_index_mapping.size() == next_molecule_id
        && "Mapping array must have value for every molecule index"
    );

    // This is the only place where we insert molecules into volume_molecules,
    // although this array size can be decreased in defragmentation
    molecules.push_back(vm_copy);
    Molecule& new_m = molecules.back();
    new_m.id = molecule_id;

    add_molecule_to_diffusion_list(new_m, time_step_index);

    return new_m;
  }

  // returns counted volume id for this position,
  // member function of Partition because some caching might be useful in the future
  geometry_object_id_t determine_counted_volume_id(const Vec3& pos);

  // auxiliary method for determine_counted_volume_id
  geometry_object_id_t find_smallest_counted_volume_recursively(const GeometryObject& obj, const Vec3& pos);

public:
  // any molecule flags are set by caller after the molecule is created by this method
  Molecule& add_volume_molecule(const Molecule& vm_copy) {

    Molecule& new_vm = add_molecule(vm_copy, true);

    new_vm.v.subpart_index = get_subpart_index(vm_copy.v.pos);
    change_reactants_map(new_vm, new_vm.v.subpart_index, true, false);

    // compute counted volume id for a new molecule
    new_vm.v.counted_volume_id = determine_counted_volume_id(new_vm.v.pos);

    return new_vm;
  }


  Molecule& add_surface_molecule(const Molecule& sm_copy) {

    Molecule& new_sm = add_molecule(sm_copy, false);
    return new_sm;
  }


  void set_molecule_as_defunct(Molecule& m) {
    // set that this molecule does not exist anymore
    m.set_is_defunct();

    change_reactants_map(m, 0/*unused*/, false, true);

    // remove from grid if it was not already removed
    if (m.is_surf() && m.s.grid_tile_index != TILE_INDEX_INVALID) {
      Grid& g = get_wall(m.s.wall_index).grid;
      g.reset_molecule_tile(m.s.grid_tile_index);
    }
  }


  // ---------------------------------- typedefs and internal structs ----------------------------------

  // this structure contains all data associated with a given diffusion time step.
  class TimeStepMoleculesData {
  public:
    // usually initialized with empty molecule_ids_ array
    TimeStepMoleculesData(float_t time_step_, const std::vector< molecule_id_t > molecule_ids_)
      : time_step(time_step_), molecule_ids(molecule_ids_) {
    }

		// diffusion time step value 
    float_t time_step;
    // molecule ids with this diffusion time step 
    std::vector<molecule_id_t> molecule_ids;
  };

  // ---------------------------------- molecule getters ----------------------------------

  // usually used as constant
  std::vector< TimeStepMoleculesData >& get_molecule_data_per_time_step_array() {
    return molecules_data_per_time_step_array;
  }


  const std::vector<molecule_id_t>& get_molecule_ids_for_time_step_index(uint32_t time_step_index) {
    assert(time_step_index < molecules_data_per_time_step_array.size());
    return molecules_data_per_time_step_array[time_step_index].molecule_ids;
  }


  const Vec3& get_origin_corner() const {
    return origin_corner;
  }
  
  const Vec3& get_opposite_corner() const {
    return opposite_corner;
  }

  const uint_set<molecule_id_t>& get_volume_molecule_reactants(subpart_index_t subpart_index, species_id_t species_id) {
    return volume_molecule_reactants_per_subpart[subpart_index].get_set(species_id);
  }


  const std::vector<Molecule>& get_molecules() const {
    return molecules;
  }


  std::vector<Molecule>& get_molecules() {
    return molecules;
  }
  

  std::vector<molecule_index_t>& get_molecule_id_to_index_mapping() {
    return molecule_id_to_index_mapping;
  }
  
  // ---------------------------------- geometry ----------------------------------
  vertex_index_t add_geometry_vertex(const Vec3 pos) {
    vertex_index_t index = geometry_vertices.size();
    geometry_vertices.push_back(pos);
    return index;
  }

  void remove_last_vertex(const vertex_index_t vertex_index) {
    assert(vertex_index == geometry_vertices.size() - 1 && "Check that we are removing known vertex failed");
    geometry_vertices.pop_back();
  }

  uint get_geometry_vertex_count() const {
    return geometry_vertices.size();
  }

  const Vec3& get_geometry_vertex(vertex_index_t i) const {
    assert(i < geometry_vertices.size());
    return geometry_vertices[i];
  }

  Vec3& get_geometry_vertex(vertex_index_t i) {
    assert(i < geometry_vertices.size());
    return geometry_vertices[i];
  }

  const Vec3& get_wall_vertex(const Wall& w, uint vertex_in_wall_index) const {
    assert(vertex_in_wall_index < VERTICES_IN_TRIANGLE);
    vertex_index_t i = w.vertex_indices[vertex_in_wall_index];
    assert(i < geometry_vertices.size());
    return get_geometry_vertex(i);
  }

  region_index_t add_region(const Region& reg) {
    region_index_t index = regions.size();
    regions.push_back(reg);
    return index;
  }

  // returns reference to the new wall, only sets id and index
  // register_wall must be called after the wall is initialized
  Wall& add_uninitialized_wall(const wall_id_t id) {
    wall_index_t index = walls.size();
    walls.push_back(Wall());
    Wall& new_wall = walls.back();

    new_wall.id = id;
    new_wall.index = index;

    return new_wall;
  }

  // when a wall is added with add_uninitialized_wall,
  // its type and vertices are not know yet, we must include the walls
  // into subvolumes and also for other purposes
  void finalize_wall_creation(const wall_index_t wall_index);

  // returns reference to the new object, only sets id
  GeometryObject& add_uninitialized_geometry_object(const geometry_object_id_t id) {
    geometry_object_index_t index = geometry_objects.size();
    geometry_objects.push_back(GeometryObject());
    GeometryObject& new_obj = geometry_objects.back();

    new_obj.id = id;
    new_obj.index = index;

    return new_obj;
  }

  GeometryObject& get_geometry_object(const geometry_object_index_t index) {
    assert(index < geometry_objects.size());
    return geometry_objects[index];
  }

  const GeometryObject& get_geometry_object(const geometry_object_index_t index) const {
    assert(index < geometry_objects.size());
    return geometry_objects[index];
  }

  const GeometryObjectVector& get_geometry_objects() const {
    return geometry_objects;
  }

  const Region& get_region(const region_index_t i) const {
    assert(i < regions.size());
    return regions[i];
  }

  const GeometryObject* find_geometry_object_by_name(const std::string& name) const {
    for (const GeometryObject& o: geometry_objects) {
      if (o.name == name) {
        return &o;
      }
    }
    return nullptr;
  }

  Region& get_region(const region_index_t i) {
    assert(i < regions.size());
    return regions[i];
  }

  const Region* find_region_by_name(const std::string& name) const {
    for (const Region& r: regions) {
      if (r.name == name) {
        return &r;
      }
    }
    return nullptr;
  }

  uint get_wall_count() const {
    return walls.size();
  }

  const Wall& get_wall(const wall_index_t i) const {
    assert(i < walls.size());
    const Wall& res = walls[i];
    assert(res.index == i && "Index of a wall must correspond to its position");
    return res;
  }

  Wall& get_wall(const wall_index_t i) {
    assert(i < walls.size());
    Wall& res = walls[i];
    assert(res.index == i && "Index of a wall must correspond to its position");
    return res;
  }

  Wall* get_wall_if_exists(const wall_index_t i) {
    if (i == WALL_INDEX_INVALID) {
      return nullptr;
    }
    assert(i < walls.size());
    Wall& res = walls[i];
    assert(res.index == i && "Index of a wall must correspond to its position");
    return &res;
  }

  // maybe we will need to filter out, e.g. just reflective surfaces
  const uint_set<wall_index_t>& get_subpart_wall_indices(const subpart_index_t subpart_index) const {
    return walls_per_subpart[subpart_index];
  }

  // returns nullptr if either the wall does not exist or the wall's grid was not initialized
  const Grid* get_wall_grid_if_exists(const wall_index_t wall_index) const {
    if (wall_index == WALL_INDEX_INVALID) {
      return nullptr;
    }

    const Wall& w = get_wall(wall_index);
    if (!w.has_initialized_grid()) {
      return nullptr;
    }

    return &w.grid;
  }

  const std::vector<wall_index_t>& get_walls_using_vertex(const vertex_index_t vertex_index) const {
    assert(vertex_index != VERTEX_INDEX_INVALID);
    assert(vertex_index < walls_using_vertex_mapping.size());
    return walls_using_vertex_mapping[vertex_index];
  }

  void add_child_of_directly_contained_counted_volume(
      const geometry_object_id_t parent_id,
      const geometry_object_id_t child_id) {
    // both must be counted volumes
    assert(get_geometry_object(parent_id).is_counted_volume);
    assert(get_geometry_object(child_id).is_counted_volume);
    directly_contained_counted_volume_objects[parent_id].insert(child_id);
  }

  void add_parent_that_encloses_counted_volume(
      const geometry_object_id_t child_id,
      const geometry_object_id_t parent_id) {
    // both must be counted volumes
    assert(get_geometry_object(parent_id).is_counted_volume);
    assert(get_geometry_object(child_id).is_counted_volume);
    enclosing_counted_volume_objects[child_id].insert(parent_id);
  }

  void set_that_object_has_no_enclosing_counted_volumes(
      const geometry_object_id_t child_id) {
    // create key with empty calue
    enclosing_counted_volume_objects[child_id] = uint_set<geometry_object_id_t>();
  }

  const uint_set<geometry_object_id_t>& get_enclosing_counted_volumes(
      const geometry_object_id_t child_id) const {
    auto it = enclosing_counted_volume_objects.find(child_id);
    assert(it != enclosing_counted_volume_objects.end());
    return it->second;
  }

  // ---------------------------------- dynamic vertices ----------------------------------
  // add information about a change of a specific vertex
  // order of calls is important (at least for now)
  void add_vertex_move(const vertex_index_t vertex_index, const Vec3& translation_vec) {
    scheduled_vertex_moves.push_back(VertexMoveInfo(vertex_index, translation_vec));
  }

  // do the actual changes of vertices
  void apply_vertex_moves();

private:
  void update_walls_per_subpart(const WallsWithTheirMovesMap& walls_with_their_moves, const bool insert);

  // automatically enlarges walls_using_vertex array
  void add_wall_using_vertex_mapping(vertex_index_t vertex_index, wall_index_t wall_index) {
    if (vertex_index >= walls_using_vertex_mapping.size()) {
      walls_using_vertex_mapping.resize(vertex_index + 1);
    }
    walls_using_vertex_mapping[vertex_index].push_back(wall_index);
  }
public:
  // ---------------------------------- other ----------------------------------
  BNG::SpeciesContainer& get_all_species() { return bng_engine.get_all_species(); }
  const BNG::SpeciesContainer& get_all_species() const { return bng_engine.get_all_species(); }

  BNG::RxnContainer& get_all_rxns() { return bng_engine.get_all_rxns(); }
  const BNG::RxnContainer& get_all_rxns() const { return bng_engine.get_all_rxns(); }

  void dump();

private:
  // left, bottom, closest (lowest z) point of the partition
  Vec3 origin_corner;
  Vec3 opposite_corner;

  // ---------------------------------- molecules ----------------------------------

  // vector containing all molecules in this partition (volume and surface)
  std::vector<Molecule> molecules;

  // contains mapping of molecule ids to indices to the molecules array
  std::vector<molecule_index_t> molecule_id_to_index_mapping;

  // id of the next molecule to be created
  // TODO_LATER: move to World
  molecule_id_t next_molecule_id;

  // indexed by diffusion time step index
  std::vector<TimeStepMoleculesData> molecules_data_per_time_step_array;

  // indexed with subpartition index, only for vol-vol reactions
  std::vector<SpeciesReactantsMap> volume_molecule_reactants_per_subpart;

  // ---------------------------------- geometry objects ----------------------------------

  std::vector<Vec3> geometry_vertices;

  // we must plan for dynamic geometry but for now its just static
  std::vector<GeometryObject> geometry_objects;

  std::vector<Wall> walls;

  std::vector<Region> regions;

  // indexed by vertex_index_t
  std::vector< std::vector<wall_index_t>> walls_using_vertex_mapping;

  // indexed by subpartition index, contains a set of wall indices (wall_index_t)
  std::vector< uint_set<wall_index_t> > walls_per_subpart;

  // key is object id and its values are objects ids directly contained in it,
  // defines a containment tree this way, used when placing a new molecule
  CountedVolumesMap directly_contained_counted_volume_objects;

  // key is object id and its values are objects ids that contain this volume
  // used when determining when counting the molecules
  CountedVolumesMap enclosing_counted_volume_objects;

  // ---------------------------------- dynamic vertices ----------------------------------
private:
  std::vector<VertexMoveInfo> scheduled_vertex_moves;

  // ---------------------------------- shared simulation configuration -------------------
public:
  partition_id_t id;

  // all these reference an object owned by a single World instance
  // enclose into something?
  const SimulationConfig& config;
  BNG::BNGEngine& bng_engine;
  SimulationStats& stats;
};

typedef std::vector<Partition> PartitionVector;

} // namespace mcell

#endif // SRC4_PARTITION_H_
