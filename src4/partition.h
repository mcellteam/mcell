/******************************************************************************
 *
 * Copyright (C) 2019,2020 by
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

#ifndef SRC4_PARTITION_H_
#define SRC4_PARTITION_H_

#include <set>

#include "../libs/bng/rxn_container.h"
#include "defines.h"
#include "dyn_vertex_structs.h"
#include "molecule.h"
#include "scheduler.h"
#include "geometry.h"
#include "simulation_config.h"
#include "species_flags_analyzer.h"
#include "../libmcell/api/shared_structs.h"

namespace Json {
class Value;
}

namespace MCell {


typedef std::map<counted_volume_index_t, uint> CountInGeomObjectMap;
typedef std::map<wall_index_t, uint> CountOnWallMap;
typedef uint_set<wall_index_t> WallsInSubpart; 

// class used to hold potential reactants of given species in a single subpart
// performance critical, therefore we are using a vector for now,
// we will probably need to change it in the future by deriving it from
//   std::map<species_id_t, uint_set<molecule_id_t> >,
class SpeciesReactantsMap:
    // vector is indexed by species_id
    public std::vector<uint_set<molecule_id_t> > {

public:
  // key must exist
  void erase_existing(const species_id_t key, const molecule_id_t id) {
    assert(key < this->size());
    (*this)[key].erase_existing(id);
  }

  void erase(const species_id_t key, const molecule_id_t id) {
    if (key >= this->size()) {
      return;
    }
    (*this)[key].erase(id);
  }

  // key is created if it does not exist
  void insert_unique(const species_id_t key, const molecule_id_t id) {
    if (key >= this->size()) {
      // resize vector
      this->resize(key + 1);
    }
    (*this)[key].insert_unique(id);
  }

  void insert(const species_id_t key, const molecule_id_t id) {
    if (key >= this->size()) {
      // resize vector
      this->resize(key + 1);
    }
    (*this)[key].insert(id);
  }

  bool contains(const species_id_t key, const molecule_id_t id) const {
    return (*this)[key].count(id) != 0;
  }

  const uint_set<molecule_id_t>& get_contained_set(const species_id_t key) const {
    if (key >= this->size()) {
      // not present, simply return an empty set
      return empty_set;
    }
    else {
      return (*this)[key];
    }
  }

  void clear_set(const species_id_t key) {
    if (key < this->size()) {
      (*this)[key].clear();
    }
  }

private:
  uint_set<molecule_id_t> empty_set;
};


struct Waypoint {
  Vec3 pos;
  counted_volume_index_t counted_volume_index;
};


/**
 * Partition class contains all molecules and other data contained in
 * one simulation block.
 */
class Partition {
public:
  Partition(
      const partition_id_t id_,
      const Vec3& origin_corner_,
      const SimulationConfig& config_,
      BNG::BNGEngine& bng_engine_,
      SimulationStats& stats_,
      SpeciesFlagsAnalyzer& species_flags_analyzer_
  )
    : origin_corner(origin_corner_),
      next_molecule_id(0),
      id(id_),
      config(config_),
      bng_engine(bng_engine_),
      stats(stats_),
      species_flags_analyzer(species_flags_analyzer_) {

    opposite_corner = origin_corner + config.partition_edge_length;

    // check that the subpart grid goes through (0, 0, 0),
    // (this point does not have to be contained in this partition)
    // required for correct function of raycast_with_endpoints,
    // round is required because values might be negative
    Vec3 how_many_subparts_from_000 = origin_corner/Vec3(config.subpartition_edge_length);
    release_assert(cmp_eq(round3(how_many_subparts_from_000), how_many_subparts_from_000, SQRT_EPS) &&
        "Partition is not aligned to the subpartition grid."
    );

    // pre-allocate volume_molecules arrays and also volume_molecule_indices_per_time_step
    uint32_t num_subparts = powu(config.num_subpartitions_per_partition, 3);
    volume_molecule_reactants_per_subpart.resize(num_subparts);
    walls_per_subpart.resize(num_subparts);

    // create an empty counted volume
    CountedVolume counted_volume_outside_all;
    counted_volume_index_t index = find_or_add_counted_volume(counted_volume_outside_all);
    assert(index == COUNTED_VOLUME_INDEX_OUTSIDE_ALL && "The empty counted volume must have index 0");

    rng_init(&aux_rng, 0);
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

  void get_molecules_ready_for_diffusion(MoleculeIdsVector& ready_vector) const {
    // select all molecules that are scheduled for this iteration and are not defunct
    ready_vector.clear();
    float_t time_it_end = stats.get_current_iteration() + 1;
    for (const Molecule& m: molecules) {
      // new products may have been scheduled for the previous iteration
      if (!m.is_defunct() && cmp_lt(m.diffusion_time, time_it_end, EPS)) {
        ready_vector.push_back(m.id);
      }
    }
  }

  const BNG::Species& get_species(const species_id_t species_id) const {
    return get_all_species().get(species_id);
  }

  bool in_this_partition(const Vec3& pos) const {
    return glm::all(glm::greaterThanEqual(pos, origin_corner))
      && glm::all(glm::lessThan(pos, opposite_corner));
  }

  bool is_subpart_index_in_range(const int index) const {
    return index >= 0 && index < (int)config.num_subpartitions_per_partition;
  }

  void get_subpart_3d_indices(const Vec3& pos, IVec3& res) const {
    // FIXME: make two variants of this function - with and without release check?
    assert(in_this_partition(pos) &&
        "Requested position is outside of a partition, usually a molecule diffused there. Please enlarge the partition size.");
    Vec3 relative_position = pos - origin_corner;
    res = relative_position * config.subpartition_edge_length_rcp;
  }

  subpart_index_t get_subpart_index_from_3d_indices(const IVec3& indices) const {
    // example: dim: 5x5x5,  (1, 2, 3) -> 1 + 2*5 + 3*5*5 = 86
    assert(indices.x < (int)config.num_subpartitions_per_partition);
    assert(indices.y < (int)config.num_subpartitions_per_partition);
    assert(indices.z < (int)config.num_subpartitions_per_partition);
    // this recomputation can be slow when done often, but we need subpart indices to be continuous
    return
        indices.x +
        indices.y * config.num_subpartitions_per_partition +
        indices.z * config.num_subpartitions_per_partition_squared;
  }

  subpart_index_t get_subpart_index_from_3d_indices(const int x, const int y, const int z) const {
    return get_subpart_index_from_3d_indices(IVec3(x, y, z));
  }


  void get_subpart_3d_indices_from_index(const subpart_index_t index, IVec3& indices) const {
    uint32_t dim = config.num_subpartitions_per_partition;
    // example: dim: 5x5x5,  86 -> (86%5, (86/5)%5, (86/(5*5))%5) = (1, 2, 3)
    indices.x = index % dim;
    indices.y = (index / dim) % dim;
    indices.z = (index / config.num_subpartitions_per_partition_squared) % dim;
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


  // called when a volume molecule is added and it is detected that
  // this is a new species
  void update_reactants_maps_for_new_species(
      const species_id_t new_species_id
  ) {
    // can the new species initiate a reaction?
    const BNG::Species& initiator_reactant_species = get_all_species().get(new_species_id);
    if (initiator_reactant_species.cant_initiate()) {
      // nothing to do
      return;
    }
    assert(initiator_reactant_species.is_vol());

    // TODO: get just vol reactants, have a separate cache in rxn container for it
    const BNG::ReactantRxnClassesMap* rxns_classes_map =
        get_all_rxns().get_bimol_rxns_for_reactant_any_compartment(new_species_id);
    if (rxns_classes_map == nullptr) {
      // nothing to do
      return;
    }

    // let's go through all molecules and update whether they can react with our new species
    for (const Molecule& m: molecules) {
      if (!m.is_vol() || m.is_defunct()) {
        continue;
      }

      assert(m.v.reactant_subpart_index != SUBPART_INDEX_INVALID);
      SpeciesReactantsMap& subpart_reactants_sp = volume_molecule_reactants_per_subpart[m.v.reactant_subpart_index];

      // check if any rxn matches the molecule that we are checking
      for (const auto& second_reactant_info: *rxns_classes_map) { 
      // TODO: filter first by species id
        const BNG::RxnClass* rxn_class =
            BNG::get_rxn_class_for_any_compartment(second_reactant_info.second); 
        if (rxn_class->get_num_reactions() == 0) {
          // there is a reaction class, but it has no reactions
          continue;
        }

        // find all molecules in this subpart that match the second reactant
        species_id_t second_species_id = second_reactant_info.first;
        if (m.species_id == second_species_id) {
          // this mapping may already exist because the new_species_id was known
          // when 'm' was added
          subpart_reactants_sp.insert(new_species_id, m.id);
        }
      }
    }
  }


  // - the reactant_subpart_index is index where the molecule was originally created,
  //   it might have moved to subpart_index in the meantime
  // - might invalidate Species references
  void change_vol_reactants_map_from_orig_to_current(Molecule& vm, bool adding, bool removing) {
    assert(vm.is_vol() && "This function is applicable only to volume mols and ignored for surface mols");
    assert(vm.v.subpart_index != SUBPART_INDEX_INVALID);
    assert(vm.v.reactant_subpart_index != SUBPART_INDEX_INVALID);
    assert(vm.v.subpart_index == get_subpart_index(vm.v.pos) && "Position and subpart must match all the time");
    // TODO: check !(vm.v.reactant_subpart_index == vm.v.reactant_subpart_index && adding && removing)
    const BNG::Species& vm_species = get_all_species().get(vm.species_id);
    if (!vm_species.has_bimol_vol_rxn()) {
      return;
    }

    // and these are indices of possible reactants with our reactant_species_id,
    // this might trigger creation of a new rxn class if the species were not seen before
    // TODO: get just vol reactants, have a separate cache in rxn container for it
    //      by this we can possibly remove handling of ANY compartment in rxn container
    const BNG::ReactantRxnClassesMap* rxns_classes_map =
        get_all_rxns().get_bimol_rxns_for_reactant_any_compartment(vm.species_id);
    if (rxns_classes_map == nullptr || rxns_classes_map->empty()) {
      // nothing to do, just update the index
      vm.v.reactant_subpart_index = vm.v.subpart_index;
      return;
    }

    // these are all the sets of indices of reactants for this particular subpartition
    SpeciesReactantsMap& subpart_reactants_orig_sp = volume_molecule_reactants_per_subpart[vm.v.reactant_subpart_index];
    SpeciesReactantsMap& subpart_reactants_new_sp = volume_molecule_reactants_per_subpart[vm.v.subpart_index];

    // we need to set/clear flag that says that second_reactant_info.first can react with reactant_species_id
    for (const auto& second_reactant_info: *rxns_classes_map) {
      species_id_t second_species_id = second_reactant_info.first;

      const BNG::RxnClass* rxn_class =
          BNG::get_rxn_class_for_any_compartment(second_reactant_info.second);

      if (rxn_class->get_num_reactions() == 0 ||
          !rxn_class->is_bimol_vol_rxn_class()) {
        // there is a reaction class, but it has no reactions
        continue;
      }

      // can the second reactant initiate a reaction with me?
      const BNG::Species& initiator_reactant_species = get_all_species().get(second_species_id);
      assert(initiator_reactant_species.is_vol());
      if (initiator_reactant_species.cant_initiate()) {
        // nothing to do
        continue;
      }

      if (removing) {
        subpart_reactants_orig_sp.erase_existing(second_species_id, vm.id);
      }
      if (adding) {
        subpart_reactants_new_sp.insert_unique(second_species_id, vm.id);
      }
    }

#ifdef DEBUG_EXTRA_CHECKS
    // check that we don't have any defunct mols
    if (!adding && removing) {
      for (auto subpart_vec: volume_molecule_reactants_per_subpart) {
        for (auto species_reactants: subpart_vec) {
          for (molecule_id_t id: species_reactants) {
            assert(!get_m(id).is_defunct());
          }
        }
      }
    }
#endif

    vm.v.reactant_subpart_index = vm.v.subpart_index;
  }


  void update_molecule_reactants_map(Molecule& vm) {
    assert(vm.v.subpart_index < volume_molecule_reactants_per_subpart.size());
    assert(vm.v.reactant_subpart_index < volume_molecule_reactants_per_subpart.size());
    if (vm.v.subpart_index == vm.v.reactant_subpart_index) {
      return; // nothing to do
    }
#ifdef DEBUG_SUBPARTITIONS
    std::cout << "Molecule " << vm.id << " changed subpartition from "
        <<  vm.v.reactant_subpart_index << " to " << vm.v.subpart_index << ".\n";
#endif

    change_vol_reactants_map_from_orig_to_current(vm, true, true);
  }

  molecule_id_t get_next_molecule_id_no_increment() {
    return next_molecule_id;
  }

private:
  // internal methods that sets molecule's id and
  // adds it to all relevant structures

  Molecule& add_molecule(const Molecule& vm_copy, const bool is_vol, const float_t release_delay_time) {
    const BNG::Species& species = get_all_species().get(vm_copy.species_id);
    assert((is_vol && species.is_vol()) || (!is_vol && species.is_surf()));

    // TODO: molecule IDs must be assigned by World
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

    // set to diffuse this or the next iteration,
    // for releases - it will be diffused this iteration
    // for new reaction products - it will be diffused next iteration because the DiffuseaAndReactEvent handles
    // the current iteration
    new_m.diffusion_time = stats.get_current_iteration() + release_delay_time;

    return new_m;
  }

  void update_species_for_new_molecule(const Molecule& m) {
    // make sure that the rxn for this species flags are up-to-date
    BNG::Species& sp = get_all_species().get(m.species_id);
    if (!sp.are_rxn_and_custom_flags_uptodate()) {
      sp.update_rxn_and_custom_flags(get_all_species(), get_all_rxns(), &species_flags_analyzer);
    }
    if (!sp.was_instantiated()) {
      // update rxn classes for this new species, may create new species and
      // invalidate Species reference
      get_all_rxns().get_bimol_rxns_for_reactant_any_compartment(sp.id);
    }
    // we must get a new reference
    get_all_species().get(m.species_id).inc_num_instantiations();
  }

public:
  // any molecule flags are set by caller after the molecule is created by this method
  Molecule& add_volume_molecule(const Molecule& vm_copy, const float_t release_delay_time = 0) {
    update_species_for_new_molecule(vm_copy);

    // TODO: use Species::is_instantiated instead of the known_vol_species
    if (known_vol_species.count(vm_copy.species_id) == 0) {
      // we must update reactant maps if new species were added
      update_reactants_maps_for_new_species(vm_copy.species_id);
      known_vol_species.insert(vm_copy.species_id);
    }

    // molecule must be added after the update because it was not fully set up
    Molecule& new_vm = add_molecule(vm_copy, true, release_delay_time);

    // and add this molecule to a map that tells which species can react with it
    new_vm.v.subpart_index = get_subpart_index(vm_copy.v.pos);
    new_vm.v.reactant_subpart_index = new_vm.v.subpart_index;
    // might invalidate species references
    change_vol_reactants_map_from_orig_to_current(new_vm, true, false);

    // we must refresh the species reference
    BNG::Species& sp_new_ref = get_all_species().get(vm_copy.species_id);
    // compute counted volume id for a new molecule, may define a new counted volume
    if (sp_new_ref.needs_counted_volume() && new_vm.v.counted_volume_index == COUNTED_VOLUME_INDEX_INVALID) {
      new_vm.set_counted_volume_and_compartment(*this, compute_counted_volume_using_waypoints(new_vm.v.pos));
    }
    return new_vm;
  }


  Molecule& add_surface_molecule(const Molecule& sm_copy, const float_t release_delay_time = 0) {
    update_species_for_new_molecule(sm_copy);

    Molecule& new_sm = add_molecule(sm_copy, false, release_delay_time);
    return new_sm;
  }


  void set_molecule_as_defunct(Molecule& m) {
    // set that this molecule does not exist anymore
    m.set_is_defunct();

    BNG::Species& sp = get_all_species().get(m.species_id);
    sp.dec_num_instantiations();

    if (m.is_vol()) {
      change_vol_reactants_map_from_orig_to_current(m, false, true);
    }

    // remove from grid if it was not already removed
    if (m.is_surf() && m.s.grid_tile_index != TILE_INDEX_INVALID) {
      Grid& g = get_wall(m.s.wall_index).grid;
      g.reset_molecule_tile(m.s.grid_tile_index);
    }
  }


  // ---------------------------------- molecule getters ----------------------------------

  const Vec3& get_origin_corner() const {
    return origin_corner;
  }
  
  const Vec3& get_opposite_corner() const {
    return opposite_corner;
  }

  const uint_set<molecule_id_t>& get_volume_molecule_reactants(subpart_index_t subpart_index, species_id_t species_id) const {
    assert(subpart_index < volume_molecule_reactants_per_subpart.size());
    return volume_molecule_reactants_per_subpart[subpart_index].get_contained_set(species_id);
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

  vertex_index_t add_or_find_geometry_vertex(const Vec3 pos) {
    // using exact comparison
    auto it = std::find(geometry_vertices.begin(), geometry_vertices.end(), pos);
    if (it == geometry_vertices.end()) {
      return add_geometry_vertex(pos);
    }
    else {
      return it - geometry_vertices.begin();
    }
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

  region_index_t add_region_and_set_its_index(Region& reg) {
    assert(reg.id != REGION_ID_INVALID);
    region_index_t index = regions.size();
    reg.index = index;
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

  GeometryObject& get_geometry_object_by_id(const geometry_object_id_t id) {
    GeometryObject& res = get_geometry_object((geometry_object_index_t)id);
    assert(res.id == res.index && "With a single partition, geom obj id == index");
    return res;
  }

  const GeometryObject& get_geometry_object_by_id(const geometry_object_id_t id) const {
    const GeometryObject& res = get_geometry_object((geometry_object_index_t)id);
    assert(res.id == res.index && "With a single partition, geom obj id == index");
    return res;
  }

  const GeometryObject* find_geometry_object(const std::string& name) const {
    for (auto& go: geometry_objects) {
      if (go.name == name) {
        return &go;
      }
    }
    return nullptr;
  }  

  const GeometryObjectVector& get_geometry_objects() const {
    return geometry_objects;
  }

  GeometryObjectVector& get_geometry_objects() {
    return geometry_objects;
  }

  Region& get_region(const region_index_t index) {
    assert(index < regions.size());
    return regions[index];
  }

  const Region& get_region(const region_index_t index) const {
    assert(index < regions.size());
    return regions[index];
  }

  const std::vector<Region>& get_regions() const {
    return regions;
  }

  const Region& get_region_by_id(const region_id_t id) const {
    const Region& res = get_region((region_index_t)id);
    assert(res.id == res.index && "With a single partition, region id == index");
    return res;
  }

  Region& get_region_by_id(const region_id_t id) {
    Region& res = get_region((region_index_t)id);
    assert(res.id == res.index && "With a single partition, region id == index");
    return res;
  }

  const GeometryObject* find_geometry_object_by_name(const std::string& name) const {
    for (const GeometryObject& o: geometry_objects) {
      if (o.name == name) {
        return &o;
      }
    }
    return nullptr;
  }

  const Region* find_region_by_name(const std::string& name) const {
    for (const Region& r: regions) {
      if (r.name == name) {
        return &r;
      }
    }
    return nullptr;
  }

  Region* find_region_by_name(const std::string& name) {
    for (Region& r: regions) {
      if (r.name == name) {
        return &r;
      }
    }
    return nullptr;
  }

  uint get_wall_count() const {
    return walls.size();
  }

  const std::vector<Wall>& get_walls() const {
    return walls;
  }

  std::vector<Wall>& get_walls() {
    return walls;
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

  const WallCollisionRejectionData& get_wall_collision_rejection_data(const wall_index_t i) const {
    assert(i < wall_collision_rejection_data.size());
    assert(wall_collision_rejection_data.size() == walls.size());
    const WallCollisionRejectionData& res = wall_collision_rejection_data[i];
    assert(res.has_identical_data_as_wall(walls[i]));
    return res;
  }

  void update_wall_collision_rejection_data(const Wall& w) {
    assert(w.index < wall_collision_rejection_data.size());
    assert(wall_collision_rejection_data.size() == walls.size());
    wall_collision_rejection_data[w.index] = w;
  }

  // maybe we will need to filter out, e.g. just reflective surfaces
  const WallsInSubpart& get_subpart_wall_indices(const subpart_index_t subpart_index) const {
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


  BNG::compartment_id_t get_reactant_compartment_id_for_counted_volume(
      const BNG::Species& species, const counted_volume_index_t counted_volume_index);

  // ---------------------------------- dynamic vertices ----------------------------------

  // may change the displacements in vertex_moves in some cases, do not use it afterwards
  void apply_vertex_moves(
      std::vector<VertexMoveInfo>& vertex_moves,
      std::set<GeometryObjectWallUnorderedPair>& colliding_walls);


  void move_waypoint_because_positioned_on_wall(
      const IVec3& waypoint_index, const bool reinitialize = true
  );

  // returns true if any of the walls overlap
  bool check_for_overlapped_walls(const Vec3& rand_vec) const;

private:

  void clamp_vertex_moves_to_wall_wall_collisions(
      std::vector<VertexMoveInfo>& vertex_moves,
      std::set<GeometryObjectWallUnorderedPair>& colliding_walls);

  void apply_vertex_moves_per_object(
      std::vector<VertexMoveInfo>& vertex_moves,
      std::set<GeometryObjectWallUnorderedPair>& colliding_walls);

  void update_walls_per_subpart(
      const WallsWithTheirMovesMap& walls_with_their_moves, const bool insert
  );

  // automatically enlarges walls_using_vertex array
  void add_wall_using_vertex_mapping(vertex_index_t vertex_index, wall_index_t wall_index) {
    if (vertex_index >= walls_using_vertex_mapping.size()) {
      walls_using_vertex_mapping.resize(vertex_index + 1);
    }
    walls_using_vertex_mapping[vertex_index].push_back(wall_index);
  }

  void initialize_waypoint(
      const IVec3& waypoint_index,
      const bool use_previous_waypoint,
      const IVec3& previous_waypoint_index,
      const bool keep_pos = false
  );

public:
  // ---------------------------------- other ----------------------------------
  // TODO: add get_species method
  BNG::SpeciesContainer& get_all_species() { return bng_engine.get_all_species(); }
  const BNG::SpeciesContainer& get_all_species() const { return bng_engine.get_all_species(); }

  BNG::RxnContainer& get_all_rxns() { return bng_engine.get_all_rxns(); }
  const BNG::RxnContainer& get_all_rxns() const { return bng_engine.get_all_rxns(); }

  // ---------------------------------- counting ----------------------------------

  // returns counted volume index for this position,
  counted_volume_index_t compute_counted_volume_from_scratch(const Vec3& pos);
  counted_volume_index_t compute_counted_volume_using_waypoints(const Vec3& pos);

  void initialize_all_waypoints();


  bool is_valid_waypoint_index(const IVec3& index3d) const {
    return
        index3d.x >= 0 &&
        index3d.x < (int)waypoints.size() &&
        index3d.y >= 0 &&
        index3d.y < (int)waypoints[index3d.x].size() &&
        index3d.z >= 0 &&
        index3d.z < (int)waypoints[index3d.x][index3d.y].size();
  }

  Waypoint& get_waypoint(const IVec3& index3d) {
    assert(is_valid_waypoint_index(index3d));
    return waypoints[index3d.x][index3d.y][index3d.z];
  }

  const Waypoint& get_waypoint(const IVec3& index3d) const {
    assert(is_valid_waypoint_index(index3d));
    return waypoints[index3d.x][index3d.y][index3d.z];
  }

  counted_volume_index_t find_or_add_counted_volume(const CountedVolume& cv);
  const CountedVolume& get_counted_volume(const counted_volume_index_t counted_volume_index) const {
    assert(counted_volumes_vector.size() == counted_volumes_set.size());
    assert(counted_volume_index < counted_volumes_vector.size());
    assert(counted_volumes_vector[counted_volume_index].index == counted_volume_index);
    return counted_volumes_vector[counted_volume_index];
  }

  void inc_rxn_in_volume_occured_count(
      const BNG::rxn_rule_id_t rxn_id, const counted_volume_index_t counted_volume_index) {
    assert(rxn_id != BNG::RXN_RULE_ID_INVALID);

    // counted_volume_index may be invalid when we are counting reactions in the whole world
    assert(counted_volume_index != COUNTED_VOLUME_INDEX_INVALID ||
        !get_all_rxns().get(rxn_id)->is_counted_in_volume_regions()
    );

    auto& map_for_rxn = rxn_counts_per_counted_volume[rxn_id];
    auto it_counted_volume = map_for_rxn.find(counted_volume_index);
    if (it_counted_volume == map_for_rxn.end()) {
      map_for_rxn[counted_volume_index] = 1;
    }
    else {
      it_counted_volume->second++;
    }
  }

  const CountInGeomObjectMap& get_rxn_in_volume_count_map(const BNG::rxn_rule_id_t rxn_rule_id) {
    // creates a new map if it was not created before
    return rxn_counts_per_counted_volume[rxn_rule_id];
  }

  void inc_rxn_on_surface_occured_count(const BNG::rxn_rule_id_t rxn_id, const wall_index_t wall_index) {
    assert(rxn_id != BNG::RXN_RULE_ID_INVALID);
    assert(wall_index != WALL_INDEX_INVALID);

    auto& map_for_rxn = rxn_counts_per_wall_volume[rxn_id];
    auto it_wall = map_for_rxn.find(wall_index);
    if (it_wall == map_for_rxn.end()) {
      map_for_rxn[wall_index] = 1;
    }
    else {
      it_wall->second++;
    }
  }

  const CountInGeomObjectMap& get_rxn_on_surface_count_map(const BNG::rxn_rule_id_t rxn_rule_id) {
    // creates a new map if it was not created before
    return rxn_counts_per_wall_volume[rxn_rule_id];
  }

  void remove_from_known_vol_species(const species_id_t species_id);
  void dump();
  void to_data_model(Json::Value& mcell) const;

  // ------------ Pymcell4 API ------------
  bool does_molecule_exist(const molecule_id_t id) {
    if (id == MOLECULE_ID_INVALID) {
      return false;
    }
    uint32_t vm_vec_index = molecule_id_to_index_mapping[id];
    if (vm_vec_index == MOLECULE_INDEX_INVALID) {
      return false;
    }
    return !get_m(id).is_defunct();
  }



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

  // indexed with subpartition index, only for vol-vol reactions
  std::vector<SpeciesReactantsMap> volume_molecule_reactants_per_subpart;

  // set that remembers which species we already saw, used to update
  // volume_molecule_reactants_per_subpart when needed
  std::set<species_id_t> known_vol_species;

  // ---------------------------------- geometry objects ----------------------------------

  std::vector<Vec3> geometry_vertices;

  std::vector<GeometryObject> geometry_objects;

  std::vector<Wall> walls;

  // this is a copy of normal and distance from origin for fast wall collision rejection,
  // stored in separate array to optimize cache performance in collide_wall
  std::vector<WallCollisionRejectionData> wall_collision_rejection_data;

  std::vector<Region> regions;

  // indexed by vertex_index_t
  std::vector< std::vector<wall_index_t>> walls_using_vertex_mapping;

  // indexed by subpartition index, contains a container wall indices (wall_index_t)
  std::vector< WallsInSubpart > walls_per_subpart;

  // ---------------------------------- counting ------------------------------------------
  // key is rxn rule id and its values are maps that contain current reaction counts for each
  // counted volume or wall
  // counts are reset every time MolOrRxnCountEvent is executed
  std::map< BNG::rxn_rule_id_t, CountInGeomObjectMap > rxn_counts_per_counted_volume;
  std::map< BNG::rxn_rule_id_t, CountOnWallMap > rxn_counts_per_wall_volume;

  // indexed by counted_volume_index_t
  std::vector<CountedVolume> counted_volumes_vector;
  // set for fast search
  std::set<CountedVolume> counted_volumes_set;

  // indexed by [x][y][z]
  std::vector< std::vector< std::vector< Waypoint > > > waypoints;

  // ---------------------------------- shared simulation configuration -------------------
public:
  partition_id_t id;

  // auxiliary random generator state
  mutable rng_state aux_rng;

  // all these reference an object owned by a single World instance
  // enclose into something?
  const SimulationConfig& config;
  BNG::BNGEngine& bng_engine;
  SimulationStats& stats;
  SpeciesFlagsAnalyzer& species_flags_analyzer;
};

typedef std::vector<Partition> PartitionVector;

} // namespace mcell

#endif // SRC4_PARTITION_H_
