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

#ifndef SRC4_MCELL3_WORLD_CONVERTER_H_
#define SRC4_MCELL3_WORLD_CONVERTER_H_

#include "api/callbacks.h"
#include "mcell_structs_shared.h"

#include "world.h"
#include <map>

namespace MCell {

class ReleaseEvent;
class RegionExprNode;
class ClampReleaseEvent;

class MCell3WorldConverter {
public:
  MCell3WorldConverter() :
    world(nullptr), callbacks(nullptr) {
  }

  ~MCell3WorldConverter() {
    delete world;
  }

  void reset();

  bool convert(volume* s);

private:
  bool convert_simulation_setup(volume* s);


  void create_uninitialized_walls_for_polygonal_object(const geom_object* o);

  bool convert_wall_and_update_regions(
      const wall* w, GeometryObject& object,
      const region_list* rl
  );
  bool convert_region(Partition& p, const region* r, region_index_t& region_index);
  bool convert_polygonal_object(const geom_object* o, const std::string& instantiate_name);
  bool convert_geometry_meta_object_recursively(volume* s, geom_object* meta_inst);
  bool convert_geometry_objects(volume* s);
  bool convert_species(volume* s);
  void create_clamp_release_event(
      const pathway* current_pathway, const BNG::RxnRule& rxn, const std::vector<species_id_t>& reactant_species_ids);
  bool convert_single_reaction(const rxn *rx);
  bool convert_rxns(volume* s);

  RegionExprNode* create_release_region_terms_recursively(
      release_evaluator* expr, ReleaseEvent& event_data
  );
  bool convert_release_events(volume* s);
  bool convert_viz_output_events(volume* s);
  bool convert_mol_or_rxn_count_events(volume* s);
  void update_and_schedule_concentration_clamps();

public:
  // contained world object
  World* world;

  API::Callbacks callbacks;

private:
  species_id_t get_mcell4_species_id(u_int mcell3_id) {
    auto it = mcell3_species_id_map.find(mcell3_id);
    assert(it != mcell3_species_id_map.end());
    assert(it->second != SPECIES_ID_INVALID);
    return it->second;
  }

  // mapping from mcell3 species id to mcell4 species id
  std::map<u_int, species_id_t> mcell3_species_id_map;


  void add_mcell4_vertex_index_mapping(const vector3* mcell3_vertex, PartitionVertexIndexPair pindex) {
    // check that if we are adding a vertex, it is exactly the same as there was before
    auto it = vector_ptr_to_vertex_index_map.find(mcell3_vertex);
    if (it != vector_ptr_to_vertex_index_map.end()) {
      // note: this check probably doesn't make sense because the mcell3 vertices
      // would have to change during conversion
      assert(it->second == pindex);
    }
    else {
      vector_ptr_to_vertex_index_map[mcell3_vertex] = pindex;
    }
  }

  PartitionVertexIndexPair get_mcell4_vertex_index(const vector3* mcell3_vertex) {
    auto it = vector_ptr_to_vertex_index_map.find(mcell3_vertex);
    assert(it != vector_ptr_to_vertex_index_map.end());
    return it->second;
  }

  std::map<const vector3*, PartitionVertexIndexPair> vector_ptr_to_vertex_index_map;

  void add_mcell4_wall_index_mapping(const wall* mcell3_wall, PartitionWallIndexPair pindex) {
    assert(wall_ptr_to_vertex_index_map.find(mcell3_wall) == wall_ptr_to_vertex_index_map.end() && "Wall mapping for this wall already exists");
    wall_ptr_to_vertex_index_map[mcell3_wall] = pindex;
  }

  PartitionWallIndexPair get_mcell4_wall_index(const wall* mcell3_wall) {
    auto it = wall_ptr_to_vertex_index_map.find(mcell3_wall);
    assert(it != wall_ptr_to_vertex_index_map.end());
    return it->second;
  }

  // use only through add_mcell4_wall_index_mapping, get_mcell4_wall_index
  std::map<const wall*, PartitionWallIndexPair> wall_ptr_to_vertex_index_map;


  void add_mcell4_region_index_mapping(const region* mcell3_region, PartitionWallIndexPair pindex) {
    assert(region_ptr_to_region_index_map.find(mcell3_region) == region_ptr_to_region_index_map.end() && "Region mapping for this wall already exists");
    region_ptr_to_region_index_map[mcell3_region] = pindex;
  }

  PartitionRegionIndexPair get_mcell4_region_index(const region* mcell3_region) {
    auto it = region_ptr_to_region_index_map.find(mcell3_region);
    assert(it != region_ptr_to_region_index_map.end());
    return it->second;
  }

  // use only through add_mcell4_region_index_mapping, get_mcell4_region_index
  std::map<const region*, PartitionRegionIndexPair> region_ptr_to_region_index_map;

  std::vector<ClampReleaseEvent*> concentration_clamps;
};


} // namespace mcell

#endif // SRC4_MCELL3_WORLD_CONVERTER_H_
