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

#include <iostream>
#include <algorithm>

#include "logging.h"

#include "partition.h"
#include "datamodel_defines.h"

#include "wall_utils.h"
#include "dyn_vertex_utils.inc"
#include "collision_utils.inc"
#include "wall_utils.inc"

using namespace std;

namespace MCell {

// when a wall is added with add_uninitialized_wall,
// its type and vertices are not know yet, we must include the walls
// into subvolumes and also for other purposes
void Partition::finalize_wall_creation(const wall_index_t wall_index) {
  Wall& w = get_wall(wall_index);

  for (vertex_index_t vi: w.vertex_indices) {
    add_wall_using_vertex_mapping(vi, wall_index);
  }

  // also insert this triangle into walls per subpartition
  SubpartIndicesVector colliding_subparts;
  GeometryUtil::wall_subparts_collision_test(*this, w, colliding_subparts);
  for (subpart_index_t subpart_index: colliding_subparts) {
    assert(subpart_index < walls_per_subpart.size());
    walls_per_subpart[subpart_index].insert_unique(wall_index);
  }
}

// remove items when 'insert' is false
void Partition::update_walls_per_subpart(const WallsWithTheirMovesMap& walls_with_their_moves, const bool insert) {
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    SubpartIndicesVector colliding_subparts;
    Wall& w = get_wall(wall_index);
    GeometryUtil::wall_subparts_collision_test(*this, w, colliding_subparts);
    for (subpart_index_t subpart_index: colliding_subparts) {
      assert(subpart_index < walls_per_subpart.size());

      if (insert) {
        walls_per_subpart[subpart_index].insert_unique(wall_index);
      }
      else {
        walls_per_subpart[subpart_index].erase_existing(wall_index);
      }
    }
  }
}


void Partition::apply_vertex_moves() {
  // 1) create a set of all affected walls with information on how much each wall moves,
  uint_set<vertex_index_t> moved_vertices_set;
  WallsWithTheirMovesMap walls_with_their_moves;
  for (const VertexMoveInfo& vertex_move_info: scheduled_vertex_moves) {

    // expecting that there we are not moving a single vertex twice
    if (moved_vertices_set.count(vertex_move_info.vertex_index) != 0) {
      mcell_error(
          "When moving dynamic vertices, each vertex may be listed just once, error for vertex with index %d.",
          vertex_move_info.vertex_index
      );
    }
    moved_vertices_set.insert(vertex_move_info.vertex_index);

    const std::vector<wall_index_t>& wall_indices = get_walls_using_vertex(vertex_move_info.vertex_index);
    for (wall_index_t wall_index: wall_indices) {
      // remember mapping wall_index -> moves
      auto it = walls_with_their_moves.find(wall_index);
      if (it == walls_with_their_moves.end()) {
        it = walls_with_their_moves.insert(make_pair(wall_index, VertexMoveInfoVector())).first;
      }
      it->second.push_back(vertex_move_info);
    }
  }

  // 2) for each wall, detect what molecules will be moved and move them right away
  //    In some cases, moving one wall might place a molecule into a path of another moved wall,
  //    however, they should be moved at the same time, so we use just the first move and skip any other further moves.
  //    Not completely sure about this, but it seems that the same behavior should be achieved when
  //    we would first collect all moves and do them later, however we are creating temporary walls,
  //    so remembering them would be more complicated.
  VolumeMoleculeMoveInfoVector volume_molecule_moves;
  uint_set<molecule_id_t> already_moved_volume_molecules;
  SurfaceMoleculeMoveInfoVector surface_molecule_moves;

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  cout << "*** Walls being moved:\n";
  for (const auto& it: walls_with_their_moves) {
    const Wall& w = get_wall(it.first);
    w.dump(*this, "  ", true);
  }
#endif

  for (const auto& it: walls_with_their_moves) {
    DynVertexUtil::collect_volume_molecules_moved_due_to_moving_wall(
        *this, it.first, it.second, already_moved_volume_molecules, volume_molecule_moves);

    DynVertexUtil::collect_surface_molecules_moved_due_to_moving_wall(
        *this, it.first, surface_molecule_moves);
  }

  // 3) get information on where these walls are and remove them
  update_walls_per_subpart(walls_with_their_moves, false);

  // 4) then we move the vertices and update relevant walls
  Geometry::update_moved_walls(*this, scheduled_vertex_moves, walls_with_their_moves);

  // 5) update subpartition info for the walls
  update_walls_per_subpart(walls_with_their_moves, true);

  // 6) move volume molecules
  for (const VolumeMoleculeMoveInfo& move_info: volume_molecule_moves) {
    DynVertexUtil::move_volume_molecule_to_closest_wall_point(*this, move_info);
  }

  // 7) move surface molecules
  // 7.1) clear grids of affected walls
  for (const auto& it: walls_with_their_moves) {
    Wall& w = get_wall(it.first);
    if (w.grid.is_initialized()) {
      w.grid.reset_all_tiles();
    }
  }

  // 7.2) do the actual movement
  // mcell3 compatibility - sort surface_molecule_moves by id in order to
  // use the same grid locations as in mcell3, not really needed for correctness
  sort( surface_molecule_moves.begin(), surface_molecule_moves.end(),
      [ ]( const auto& lhs, const auto& rhs )
      {
        return lhs.molecule_id < rhs.molecule_id;
      }
  );
  for (const SurfaceMoleculeMoveInfo& move_info: surface_molecule_moves) {
    // finally fix positions of the surface molecules
    DynVertexUtil::move_surface_molecule_to_closest_wall_point(*this, move_info);
  }

  scheduled_vertex_moves.clear();
}

void Partition::dump() {
  GeometryObject::dump_array(*this, geometry_objects);

  Region::dump_array(regions);

  for (size_t i = 0; i < walls_per_subpart.size(); i++) {
    if (!walls_per_subpart[i].empty()) {
      Vec3 llf, urb;
      get_subpart_llf_point(i, llf);
      urb = llf + Vec3(config.subpartition_edge_length);

      cout << "subpart: " << i << ", llf: " << llf << ", urb: " << urb << "\n";
      walls_per_subpart[i].dump("Indices contained in a subpartition");
    }
  }
}


void Partition::create_waypoint(const IVec3& index3d) {
  // array was allocated
  Waypoint& waypoint = get_waypoint(index3d);
  waypoint.pos =
      origin_corner +
      Vec3(config.subpartition_edge_length) * Vec3(index3d) +
      Vec3(config.subpartition_edge_length / 2);

  // figure out in which counted volumes is this waypoint present
  waypoint.counted_volume_index = compute_counted_volume_from_scratch(waypoint.pos);
}


void Partition::initialize_waypoints() {
  if (!config.has_intersecting_counted_objects) {
    // no need to initialize
    return;
  }

  // each center of a subpartition has a waypoint
  uint num_waypoints_per_dimension = config.num_subpartitions_per_partition;

  mcell_log(
      "Initializing %d waypoints... "
      "(note: the current implementation is not very efficient, if you encounter performance issues, "
      "please contact the MCell team or consider not using overlapping counted objects)",
      powu(num_waypoints_per_dimension, 3));

  waypoints.resize(num_waypoints_per_dimension);
  for (uint x = 0; x < num_waypoints_per_dimension; x++) {
    waypoints[x].resize(num_waypoints_per_dimension);
    for (uint y = 0; y < num_waypoints_per_dimension; y++) {
      waypoints[x][y].resize(num_waypoints_per_dimension);
      for (uint z = 0; z < num_waypoints_per_dimension; z++) {
        create_waypoint(IVec3(x, y, z));
      }
    }
  }
}


// method is in .cpp file because it uses inlined implementation
counted_volume_index_t Partition::compute_counted_volume_from_scratch(const Vec3& pos) {
  return CollisionUtil::compute_counted_volume_for_pos(*this, pos);
}


counted_volume_index_t Partition::compute_counted_volume_using_waypoints(const Vec3& pos) {
  return CollisionUtil::compute_counted_volume_using_waypoints(*this, pos);
}


counted_volume_index_t Partition::find_or_add_counted_volume(const CountedVolume& cv) {
  assert(counted_volumes_vector.size() == counted_volumes_set.size());

  auto it = counted_volumes_set.find(cv);
  if (it != counted_volumes_set.end()) {
    return it->index;
  }
  else {
    // we must add this new item
    CountedVolume cv_copy = cv;
    cv_copy.index = counted_volumes_set.size();
    counted_volumes_set.insert(cv_copy);
    counted_volumes_vector.push_back(cv_copy);
    return cv_copy.index;
  }
}

void Partition::to_data_model(Json::Value& mcell) const {

  Json::Value& geometrical_objects = mcell[KEY_GEOMETRICAL_OBJECTS];
  Json::Value& object_list = geometrical_objects[KEY_OBJECT_LIST];
  if (object_list.isNull()) {
    object_list = Json::Value(Json::arrayValue);
  }

  for (const GeometryObject& g: geometry_objects) {
    Json::Value object;
    g.to_data_model(*this, config, object);
    object_list.append(object);
  }

  // mod_surf_regions
  Json::Value& modify_surface_regions = mcell[KEY_MODIFY_SURFACE_REGIONS];
  DMUtil::add_version(modify_surface_regions, VER_DM_2014_10_24_1638);
  Json::Value& modify_surface_regions_list = modify_surface_regions[KEY_MODIFY_SURFACE_REGIONS_LIST];
  modify_surface_regions_list = Json::Value(Json::arrayValue);
  for (const Region& reg: regions) {
    if (reg.is_reactive()) {
      Json::Value modify_surface_region;
      reg.to_data_model(*this, modify_surface_region);
      modify_surface_regions_list.append(modify_surface_region);
    }
  }
}

} // namespace mcell
