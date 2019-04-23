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

#include "mcell_structs.h"

#ifdef __cplusplus
extern "C"
#endif
bool mcell4_convert_mcell3_volume(struct volume* s);

#ifdef __cplusplus
extern "C"
#endif
bool mcell4_run_simulation();

#ifdef __cplusplus
extern "C"
#endif
void mcell4_delete_world();


// following code is only for C++
#ifdef __cplusplus

#include "world.h"
#include <map>

namespace mcell {

class mcell3_world_converter {
public:
  mcell3_world_converter() :
    world(nullptr) {
  }

  ~mcell3_world_converter() {
    delete world;
  }

  void reset();

  bool convert(volume* s);

  bool convert_simulation_setup(volume* s);

  bool convert_wall(
      wall* w,
      wall_t& res_wall,
      std::vector<partition_vertex_index_pair_t>& res_vertices);
  bool convert_polygonal_object(object* o);
  bool convert_geometry_objects(volume* s);

  void create_diffusion_events();
  bool convert_species_and_create_diffusion_events(volume* s);
  bool convert_single_reaction(rxn *rx);
  bool convert_reactions(volume* s);
  bool convert_release_events(volume* s);
  bool convert_viz_output_events(volume* s);

  // contained world object
  world_t* world;

private:
  species_id_t get_mcell4_species_id(u_int mcell3_id) {
    auto it = mcell3_species_id_map.find(mcell3_id);
    assert(it != mcell3_species_id_map.end());
    return it->second;
  }

  // mapping from mcell3 species id to mcell4 species id
  std::map<u_int, species_id_t> mcell3_species_id_map;

  partition_vertex_index_pair_t get_mcell4_vertex_index(vector3* mcell3_vertex) {
    auto it = vector_ptr_to_vertex_index_map.find(mcell3_vertex);
    assert(it != vector_ptr_to_vertex_index_map.end());
    return it->second;
  }

  std::map<vector3*, partition_vertex_index_pair_t> vector_ptr_to_vertex_index_map;
};


} // namespace mcell

#endif // #ifdef __cplusplus

#endif // SRC4_MCELL3_WORLD_CONVERTER_H_
