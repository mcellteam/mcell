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

  w.present_in_subparts.clear();

  // also insert this triangle into walls per subpartition
  SubpartIndicesVector colliding_subparts;
  GeometryUtil::wall_subparts_collision_test(*this, w, colliding_subparts);
  for (subpart_index_t subpart_index: colliding_subparts) {
    assert(subpart_index < walls_per_subpart.size());

    // mapping subpart->wall
    walls_per_subpart[subpart_index].push_back(wall_index);

    // mapping wall->subpart
    w.present_in_subparts.insert(subpart_index); // TODO: use insert_unique
  }

  // make a cache-optimized copy of certain fields from Wall
  assert(wall_collision_rejection_data.size() == wall_index);
  wall_collision_rejection_data.push_back(
      WallCollisionRejectionData(w.normal, w.distance_to_origin));
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

      release_assert(false);
      /*if (insert) {
        walls_per_subpart[subpart_index].insert_unique(wall_index);
        w.present_in_subparts.insert(subpart_index);
      }
      else {
        walls_per_subpart[subpart_index].erase_existing(wall_index);
        w.present_in_subparts.erase(subpart_index);
      }*/
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
      //walls_per_subpart[i].dump("Indices contained in a subpartition");
    }
  }
}


// methods get_num_crossed_region_walls and get_num_crossed_walls_per_object (TODO)
// may fail with WALL_REDO, this means that a point is on a wall and
// therefore it cannot be safely determined where it belongs
// waypoints are used as an optimization so we can place them anywhere we like
void Partition::move_waypoint_because_positioned_on_wall(
    const IVec3& waypoint_index, const bool reinitialize) {
  Waypoint& waypoint = get_waypoint(waypoint_index);
  waypoint.pos = waypoint.pos + Vec3(SQRT_EPS);
  if (reinitialize) {
    initialize_waypoint(waypoint_index, false, IVec3(0), true);
  }
}


// keep_pos is false by default, when true, the waypoint position is reused and
// this function only updates the counted_volume_index for the new position
void Partition::initialize_waypoint(
    const IVec3& waypoint_index,
    const bool use_previous_waypoint,
    const IVec3& previous_waypoint_index,
    const bool keep_pos
) {
  // array was allocated
  Waypoint& waypoint = get_waypoint(waypoint_index);
  if (!keep_pos) {
    waypoint.pos =
        origin_corner +
        Vec3(config.subpartition_edge_length) * Vec3(waypoint_index) + // llf of a subpart
        Vec3(config.subpartition_edge_length / 2)
    ;
  }

  // it makes sense to compute counted_volume_index if there are any
  // counted objects
  if (config.has_intersecting_counted_objects) {
    // see if we can reuse the previous waypoint to accelerate computation
    if (use_previous_waypoint) {
      Waypoint& previous_waypoint = get_waypoint(previous_waypoint_index);
      map<geometry_object_index_t, uint> num_crossed_walls_per_object;

      bool must_redo_test = false;
      CollisionUtil::get_num_crossed_walls_per_object(
          *this, waypoint.pos, previous_waypoint.pos,
          num_crossed_walls_per_object, must_redo_test
      );
      if (must_redo_test) {
        // updates values referenced by waypoint
        // the redo can be also due to the previous waypoint, so we
        // will check whether the waypoint is contained without the previous waypoint
        // do not call reinitialize to avoid recursion
        move_waypoint_because_positioned_on_wall(waypoint_index, false);
      }
      // ok, test passed safely
      else if (num_crossed_walls_per_object.empty()) {
        // ok, we can simply copy counted volume from the previous waypoint
        waypoint.counted_volume_index = previous_waypoint.counted_volume_index;
        return;
      }
    }

    // figure out in which counted volumes is this waypoint present
    waypoint.counted_volume_index = compute_counted_volume_from_scratch(waypoint.pos);
  }
  else {
    waypoint.counted_volume_index = COUNTED_VOLUME_INDEX_OUTSIDE_ALL;
  }
}


void Partition::initialize_all_waypoints() {
  // we need waypoints to be initialized all the time because they
  // are used not just when dealing with counted volumes, but also
  // by regions when detecting whether a point is inside

  // each center of a subpartition has a waypoint
  uint num_waypoints_per_dimension = config.num_subpartitions_per_partition;

  mcell_log("Initializing %d waypoints... ", powu(num_waypoints_per_dimension, 3));

  bool use_previous_waypoint = false;
  IVec3 previous_waypoint_index;

  waypoints.resize(num_waypoints_per_dimension);
  for (uint x = 0; x < num_waypoints_per_dimension; x++) {
    waypoints[x].resize(num_waypoints_per_dimension);
    for (uint y = 0; y < num_waypoints_per_dimension; y++) {
      waypoints[x][y].resize(num_waypoints_per_dimension);
      for (uint z = 0; z < num_waypoints_per_dimension; z++) {
        initialize_waypoint(IVec3(x, y, z), use_previous_waypoint, previous_waypoint_index);

        use_previous_waypoint = true;
        previous_waypoint_index = IVec3(x, y, z);
      }

      // starting a new line - start from scratch
      use_previous_waypoint = false;
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
    if (reg.is_reactive() || reg.has_initial_molecules()) {
      Json::Value modify_surface_region;
      reg.to_data_model(*this, modify_surface_region);
      modify_surface_regions_list.append(modify_surface_region);
    }
  }
}

} // namespace mcell
