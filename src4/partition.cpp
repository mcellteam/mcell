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

#include <iostream>
#include <algorithm>

#include "logging.h"

#include "partition.h"
#include "datamodel_defines.h"

#include "wall_utils.h"
#include "dyn_vertex_utils.inl"
#include "collision_utils.inl"
#include "wall_utils.inl"
#include "wall_overlap.h"

using namespace std;

namespace MCell {

static void zero_out_lowest_bits_pos(pos_t* val) {
  uint tmp;
  memcpy(&tmp, val, sizeof(uint));
  tmp &= 0xFFFFF000;
  memcpy(val, &tmp, sizeof(uint));
}

Partition::Partition(
    const partition_id_t id_,
    const Vec3& origin_corner_,
    const SimulationConfig& config_,
    BNG::BNGEngine& bng_engine_,
    SimulationStats& stats_
)
  : origin_corner(origin_corner_),
    next_molecule_id(0),
    volume_molecule_reactants_per_reactant_class(config_.num_subpartitions),
    id(id_),
    config(config_),
    bng_engine(bng_engine_),
    stats(stats_) {

#if POS_T_BYTES == 4
  zero_out_lowest_bits_pos(&origin_corner.x);
  zero_out_lowest_bits_pos(&origin_corner.y);
  zero_out_lowest_bits_pos(&origin_corner.z);
#endif

  opposite_corner = origin_corner + config.partition_edge_length;

  // check that the subpart grid goes through (0, 0, 0),
  // (this point does not have to be contained in this partition)
  // required for correct function of raycast_with_endpoints,
  // round is required because values might be negative
  Vec3 how_many_subparts_from_000 = origin_corner/Vec3(config.subpartition_edge_length);
  release_assert(cmp_eq(round3(how_many_subparts_from_000), how_many_subparts_from_000, POS_SQRT_EPS) &&
      "Partition is not aligned to the subpartition grid."
  );

  // pre-allocate volume_molecules arrays and also volume_molecule_indices_per_time_step
  walls_per_subpart.resize(config.num_subpartitions);

  // create an empty counted volume
  CountedVolume counted_volume_outside_all;
  counted_volume_index_t index = find_or_add_counted_volume(counted_volume_outside_all);
  assert(index == COUNTED_VOLUME_INDEX_OUTSIDE_ALL && "The empty counted volume must have index 0");

  rng_init(&aux_rng, 0);
}

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
  GeometryUtils::wall_subparts_collision_test(*this, w, colliding_subparts);
  for (subpart_index_t subpart_index: colliding_subparts) {
    assert(subpart_index < walls_per_subpart.size());

    // mapping subpart->wall
    walls_per_subpart[subpart_index].insert(wall_index);

    // mapping wall->subpart
    w.present_in_subparts.insert(subpart_index); // TODO: use insert_unique
  }

  // make a cache-optimized copy of certain fields from Wall
  assert(wall_collision_rejection_data.size() == wall_index);
  wall_collision_rejection_data.push_back(w);
}

// remove items when 'insert' is false
void Partition::update_walls_per_subpart(const WallsWithTheirMovesMap& walls_with_their_moves, const bool insert) {
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    SubpartIndicesVector colliding_subparts;
    Wall& w = get_wall(wall_index);
    GeometryUtils::wall_subparts_collision_test(*this, w, colliding_subparts);
    assert(insert || SubpartIndicesSet(colliding_subparts.begin(), colliding_subparts.end()) == w.present_in_subparts);

    for (subpart_index_t subpart_index: colliding_subparts) {
      assert(subpart_index < walls_per_subpart.size());
      if (insert) {
        walls_per_subpart[subpart_index].insert_unique(wall_index);
        w.present_in_subparts.insert(subpart_index);
      }
      else {
        walls_per_subpart[subpart_index].erase_existing(wall_index);
        w.present_in_subparts.erase(subpart_index);
      }
    }
  }
}


void Partition::apply_vertex_moves(
    std::vector<VertexMoveInfo>& vertex_moves,
    std::set<GeometryObjectWallUnorderedPair>& colliding_walls) {
  colliding_walls.clear();

  // due to wall-wall collision detection, we must move vertices of each object separately

  // is there a single object that we are moving?
  geometry_object_id_t object_id = GEOMETRY_OBJECT_ID_INVALID;
  bool single_object = true;
  for (const VertexMoveInfo& vertex_move_info: vertex_moves) {
    if (object_id == GEOMETRY_OBJECT_ID_INVALID || object_id == vertex_move_info.geometry_object_id) {
      object_id = vertex_move_info.geometry_object_id;
    }
    else {
      single_object = false;
      break;
    }
  }

  if (single_object) {
    apply_vertex_moves_per_object(vertex_moves, colliding_walls);
  }
  else {
    // we must make a separate vector for each
    map<geometry_object_id_t, vector<VertexMoveInfo>> vertex_moves_per_object;
    for (const VertexMoveInfo& vertex_move_info: vertex_moves) {
      vertex_moves_per_object[vertex_move_info.geometry_object_id].push_back(vertex_move_info);
    }

    // and process objects one by one
    for (auto& pair_id_moves: vertex_moves_per_object) {
      apply_vertex_moves_per_object(pair_id_moves.second, colliding_walls);
    }
  }
}


void Partition::apply_vertex_moves_per_object(
    std::vector<VertexMoveInfo>& vertex_moves,
    std::set<GeometryObjectWallUnorderedPair>& colliding_walls) {

  // 0) clamp maximum movement
  clamp_vertex_moves_to_wall_wall_collisions(vertex_moves, colliding_walls);

  // 1) create a set of all affected walls with information on how much each wall moves,
  uint_set<vertex_index_t> moved_vertices_set;
  WallsWithTheirMovesMap walls_with_their_moves;
  for (VertexMoveInfo& vertex_move_info: vertex_moves) {

    // expecting that there we are not moving a single vertex twice
    if (moved_vertices_set.count(vertex_move_info.vertex_index) != 0) {
      mcell_error(
          "When moving dynamic vertices, each vertex may be listed just once, error for vertex with index %d.",
          vertex_move_info.vertex_index
      );
    }
    moved_vertices_set.insert(vertex_move_info.vertex_index);

    const std::vector<wall_index_t>& wall_indices = get_walls_using_vertex(vertex_move_info.vertex_index);

    // first check whether we can move all walls belonging to the vertex
    for (wall_index_t wall_index: wall_indices) {
      if (!get_wall(wall_index).is_movable) {
        // we must not move this vertex
        vertex_move_info.vertex_walls_are_movable = false;
        break;
      }
    }

    // if the move is applicable, create mapping in walls_with_their_moves
    if (vertex_move_info.vertex_walls_are_movable) {
      for (wall_index_t wall_index: wall_indices) {
        // remember mapping wall_index -> moves
        auto it = walls_with_their_moves.find(wall_index);
        if (it == walls_with_their_moves.end()) {
          it = walls_with_their_moves.insert(make_pair(wall_index, WallMoveInfo())).first;
        }
        it->second.vertex_moves.push_back(vertex_move_info);
      }
    }
  }

  // set whether walls change area
  for (auto& it: walls_with_their_moves) {
    // FIXME: reenable and improve check, we do not want molecules to move too far
    it.second.wall_changes_area = true; /*
        !(it.second.vertex_moves.size() == 3 &&
          it.second.vertex_moves[0].displacement == it.second.vertex_moves[1].displacement &&
          it.second.vertex_moves[1].displacement == it.second.vertex_moves[2].displacement);*/
  }


  // 2) for each wall, detect what molecules will be moved and move them right away
  //    In some cases, moving one wall might place a molecule into a path of another moved wall,
  //    however, they should be moved at the same time, so we use just the first move and skip any other further moves.
  //    Not completely sure about this, but it seems that the same behavior should be achieved when
  //    we would first collect all moves and do them later, however we are creating temporary walls,
  //    so remembering them would be more complicated.
  VolumeMoleculeMoveInfoVector volume_molecule_moves;
  MoleculeIdsSet already_moved_volume_molecules;
  SurfaceMoleculeMoveInfoVector surface_molecule_moves;

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  cout << "*** Walls being moved:\n";
  for (const auto& it: walls_with_their_moves) {
    const Wall& w = get_wall(it.first);
    w.dump(*this, "  ", true);
  }
#endif

  for (const auto& it: walls_with_their_moves) {
    DynVertexUtils::collect_volume_molecules_moved_due_to_moving_wall(
        *this, it.first, it.second, already_moved_volume_molecules, volume_molecule_moves);

    // does the area of the wall change?
    if (it.second.wall_changes_area) {
      DynVertexUtils::collect_surface_molecules_moved_due_to_moving_wall(
          *this, it.first, surface_molecule_moves);
    }
  }

  // 3) get information on where these walls are and remove them
  update_walls_per_subpart(walls_with_their_moves, false);

  // 4) then we move the vertices and update relevant walls
  Geometry::update_moved_walls(*this, vertex_moves, walls_with_their_moves);

  // 5) update subpartition info for the walls
  update_walls_per_subpart(walls_with_their_moves, true);

  // 6) move volume molecules
  for (const VolumeMoleculeMoveInfo& move_info: volume_molecule_moves) {
    DynVertexUtils::move_volume_molecule_to_closest_wall_point(*this, move_info);
  }

  // 7) move surface molecules
  // 7.1) clear grids of affected walls
  for (const auto& it: walls_with_their_moves) {
    if (it.second.wall_changes_area) {
      Wall& w = get_wall(it.first);
      if (w.grid.is_initialized()) {
        w.grid.reset_all_tiles();
      }
    }
  }

  // 7.2) do the actual movement
  // mcell3 compatibility - sort surface_molecule_moves by id in order to
  // use the same grid locations as in mcell3, not really needed for correctness
  sort( surface_molecule_moves.begin(), surface_molecule_moves.end(),
      [ ]( const SurfaceMoleculeMoveInfo& lhs, const SurfaceMoleculeMoveInfo& rhs )
      {
        return lhs.molecule_id < rhs.molecule_id;
      }
  );
  for (const SurfaceMoleculeMoveInfo& move_info: surface_molecule_moves) {
    // finally fix positions of the surface molecules (only those that are on walls that change area)
    DynVertexUtils::move_surface_molecule_to_closest_wall_point(*this, move_info);
  }
}


void Partition::clamp_vertex_moves_to_wall_wall_collisions(
    std::vector<VertexMoveInfo>& vertex_moves,
    std::set<GeometryObjectWallUnorderedPair>& colliding_walls) {

  // FIXME: the test does not handle cases when there is a sharp edge or the colliding triangle is smaller,
  // than the triangle we are moving, but for now let's keep it as it is...

  // this is a first test that checks whether our wall wont simply collide
  // if we send a ray from each of the moved
  for (VertexMoveInfo& vertex_move_info: vertex_moves) {
    assert(vertex_move_info.partition_id == id);

    const std::vector<wall_index_t>& wall_indices = get_walls_using_vertex(vertex_move_info.vertex_index);
    Vec3& displacement = vertex_move_info.displacement;
    // check that object id is consistent
    assert(!wall_indices.empty());
    assert(get_wall(wall_indices[0]).object_id == vertex_move_info.geometry_object_id);

    const Vec3& pos = get_geometry_vertex(vertex_move_info.vertex_index);
    map<geometry_object_index_t, uint> num_crossed_walls_per_object_ignored;
    bool must_redo_test;
    wall_index_t closest_hit_wall_index;

    Vec3 dir(
        (cmp_eq(displacement.x, (pos_t)0) ? 0 : ((displacement.x > 0) ? 1 : -1)),
        (cmp_eq(displacement.y, (pos_t)0) ? 0 : ((displacement.y > 0) ? 1 : -1)),
        (cmp_eq(displacement.z, (pos_t)0) ? 0 : ((displacement.z > 0) ? 1 : -1))
    );
    // extra displacement to make sure that we will end up in front of a wall,
    // not on it
    Vec3 wall_gap_displacement = dir * Vec3(MIN_WALL_GAP);

    // TODO: preferably use some function that does not collect what we hit,
    // but we do not have such
    stime_t collision_time = CollisionUtils::get_num_crossed_walls_per_object(
        *this, pos, pos + displacement + wall_gap_displacement, false,
        num_crossed_walls_per_object_ignored, must_redo_test,
        vertex_move_info.geometry_object_id,
        &closest_hit_wall_index
    );

    bool ignore_hit = false;
    if (must_redo_test) {
      mcell_log(
          "Internal warning: wall collision redo in clamp_vertex_moves_to_wall_wall_collisions is not handled yet, "
          "the vertex move that caused it won't be applied."
      );
      ignore_hit = true;
      collision_time = 0;
    }

    assert(collision_time >= 0.0 && collision_time <= 1.0);
    if (collision_time != 1.0) {

      if (collision_time < STIME_SQRT_EPS) {
        collision_time = 0.0;
      }
      else {
        collision_time -= STIME_SQRT_EPS;
      }

      // clamp the displacement (it is a reference)
      displacement = displacement * Vec3(collision_time);

      // decrement by MIN_WALL_GAP because we hit a wall and would like to stay a bit in front of it
      displacement = displacement - wall_gap_displacement;

      if (!ignore_hit) {
        assert(closest_hit_wall_index != WALL_INDEX_INVALID);

        // and remember collisions for each of the walls that are moved by this single vertex
        const Wall& hit_wall = get_wall(closest_hit_wall_index);
        for (wall_index_t wi: wall_indices) {
          colliding_walls.insert(
              GeometryObjectWallUnorderedPair(
                  vertex_move_info.geometry_object_id,
                  wi,
                  hit_wall.object_id,
                  hit_wall.index
              )
          );
        }
      }
    }
  }
}


void Partition::dump(const bool with_geometry) {
  if (with_geometry) {
    GeometryObject::dump_array(*this, geometry_objects);
  }

  Region::dump_array(regions);

#ifndef NDEBUG
  if (with_geometry) {
    for (size_t i = 0; i < walls_per_subpart.size(); i++) {
      if (!walls_per_subpart[i].empty()) {
        Vec3 llf, urb;
        get_subpart_llf_point(i, llf);
        urb = llf + Vec3(config.subpartition_edge_length);

        cout << "subpart: " << i << ", llf: " << llf << ", urb: " << urb << "\n  ";

        for (wall_index_t wi: walls_per_subpart[i]) {
          cout << wi << ", ";
        }
        cout << "\n";
      }
    }
  }
#endif
}


// methods get_num_crossed_region_walls and get_num_crossed_walls_per_object
// may fail with WALL_REDO, this means that a point is on a wall and
// therefore it cannot be safely determined where it belongs
// waypoints are used as an optimization so we can place them anywhere we like
void Partition::move_waypoint_because_positioned_on_wall(
    const IVec3& waypoint_index, const bool reinitialize) {
  Waypoint& waypoint = get_waypoint(waypoint_index);
  // rng_dbl when used in Vec3 ctor causes compiler to print warnings
  // that a constant is subtracted from uint
  double r1 = rng_dbl(&aux_rng);
  double r2 = rng_dbl(&aux_rng);
  double r3 = rng_dbl(&aux_rng);

  Vec3 random_displacement = Vec3(POS_SQRT_EPS * r1, POS_SQRT_EPS * r2, POS_SQRT_EPS * r3);

  waypoint.pos = waypoint.pos + random_displacement;
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

  // make sure that the waypoint is not on any wall, this is required for correct function
  // of regions as well
  pos_t dist2;
  do {
    wall_index_t wall_index_ignored;
    Vec2 wall_pos2d_ignored;
    dist2 = WallUtils::find_closest_wall_any_object(
        *this, waypoint.pos, POS_SQRT_EPS, false, wall_index_ignored, wall_pos2d_ignored);
    if (dist2 < POS_EPS) {
       move_waypoint_because_positioned_on_wall(waypoint_index, false);
    }
  } while (dist2 < POS_EPS);

  // another required check is that there must not be an edge between the current and
  // previous waypoint, this is required for initialization of regions
  if (use_previous_waypoint) {
    bool redo;
    int num_attempts = 0;
    Vec3& previous_waypoint_pos = get_waypoint(previous_waypoint_index).pos;
    do {
      map<geometry_object_index_t, uint> num_crossed_walls_per_object;

      CollisionUtils::get_num_crossed_walls_per_object(
          *this, waypoint.pos, previous_waypoint_pos, false, // all walls
          num_crossed_walls_per_object, redo
      );

      if (redo) {
         move_waypoint_because_positioned_on_wall(waypoint_index, false);

         if (num_attempts >= 5) {
           // give up because we are getting redos due to the previous waypoint,
           // count from scratch
           break;
         }
         num_attempts++;
      }
    } while (redo);
  }

  // it makes sense to compute counted_volume_index if there are any
  // counted objects
  if (config.has_intersecting_counted_objects) {
    // see if we can reuse the previous waypoint to accelerate computation
    if (use_previous_waypoint) {
      Waypoint& previous_waypoint = get_waypoint(previous_waypoint_index);
      map<geometry_object_index_t, uint> num_crossed_walls_per_object;

      bool must_redo_test = false;
      CollisionUtils::get_num_crossed_walls_per_object(
          *this, waypoint.pos, previous_waypoint.pos, true, // only counted objects
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
  uint num_waypoints_per_dimension = config.num_subpartitions_per_partition_edge;

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


bool Partition::check_for_overlapped_walls(const Vec3& rand_vec) const {

  typedef pair<wall_index_t, double> WallDprodPair;
  vector<WallDprodPair> wall_indices_w_dprod;
  for (const Wall& w: walls) {
    double d_prod = dot(rand_vec, w.normal);

    /* we want to place walls with opposite normals into
       neighboring positions in the sorted list */
    if (d_prod < 0) {
      d_prod = -d_prod;
    }

    wall_indices_w_dprod.push_back(make_pair(w.index, d_prod));
  }

  // sort according to dprod
  sort(wall_indices_w_dprod.begin(), wall_indices_w_dprod.end(),
      [](const WallDprodPair& a, const WallDprodPair& b) -> bool {
          return a.second < b.second;
      }
  );


  for (size_t i = 0; i < wall_indices_w_dprod.size(); i++) {
    WallDprodPair& wd = wall_indices_w_dprod[i];
    const Wall& w1 = get_wall(wd.first);

    size_t next_index = i + 1;
    while (next_index < wall_indices_w_dprod.size() &&
        (!distinguishable_f(wd.second, wall_indices_w_dprod[next_index].second, EPS))) {

      /* there may be several walls with the same (or mirror) oriented normals */
      const Wall& w2 = get_wall(wall_indices_w_dprod[next_index].first);

      if (WallOverlap::are_coplanar(*this, w1, w2, MESH_DISTINCTIVE_EPS) &&
           (WallOverlap::are_coincident(*this, w1, w2, MESH_DISTINCTIVE_EPS) ||
            WallOverlap::coplanar_walls_overlap(*this, w1, w2))
      ) {

        const string& obj1_name = get_geometry_object(w1.object_id).name;
        const string& obj2_name = get_geometry_object(w2.object_id).name;

        mcell_log(
            "Error: walls are overlapped: wall side %d from '%s' and wall "
            "side %d from '%s'.",
            w1.side, obj1_name.c_str(), w2.side, obj2_name.c_str()
        );
        return false;
      }

      next_index++;
    }
  }
  return true;
}


// method is in .cpp file because it uses inlined implementation
counted_volume_index_t Partition::compute_counted_volume_from_scratch(const Vec3& pos) {
  return CollisionUtils::compute_counted_volume_for_pos(*this, pos);
}


counted_volume_index_t Partition::compute_counted_volume_using_waypoints(const Vec3& pos) {
  return CollisionUtils::compute_counted_volume_using_waypoints(*this, pos);
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


// returns compartment ID, all compartments are represented by a counted volume
// so we can use the counted volume information to determine it
BNG::compartment_id_t Partition::get_compartment_id_for_counted_volume(const counted_volume_index_t counted_volume_index) {

  // caching, the mapping object->compartment is static and does not change even with dynamic geometry
  auto it = counted_volume_index_to_compartment_id_cache.find(counted_volume_index);
  if (it != counted_volume_index_to_compartment_id_cache.end()) {
    return it->second;
  }

  const CountedVolume& cv = get_counted_volume(counted_volume_index);

  // first collect all compartments we are in
  set<BNG::compartment_id_t> inside_of_compartments;
  for (geometry_object_index_t obj_index: cv.contained_in_objects) {
    const GeometryObject& obj = get_geometry_object(obj_index);
    if (obj.represents_compartment()) {
      inside_of_compartments.insert(obj.vol_compartment_id);
    }
  }

  if (inside_of_compartments.empty()) {
    // we are outside of any defined compartment
    counted_volume_index_to_compartment_id_cache[counted_volume_index] = BNG::COMPARTMENT_ID_NONE;
    return BNG::COMPARTMENT_ID_NONE;
  }

  if (inside_of_compartments.size() == 1) {
    // there is just one compartment
    BNG::compartment_id_t res = *inside_of_compartments.begin();
    counted_volume_index_to_compartment_id_cache[counted_volume_index] = res;
    return res;
  }

  // find the lowest level child because volume compartments are always hierarchical
  // and if we are inside of objects that correspond to EC and CP, it means that we are in CP
  BNG::compartment_id_t current = *inside_of_compartments.begin();
  bool child_found;
  do {
    const BNG::Compartment& comp = bng_engine.get_data().get_compartment(current);
    assert(comp.is_3d);

    // which child are we possibly in?
    child_found = false;
    for (BNG::compartment_id_t child: comp.children_compartments) {
      const BNG::Compartment& child_comp = bng_engine.get_data().get_compartment(child);
      if (!child_comp.is_3d) {
        // go through children of the 2d compartment (should be just one usually)
        for (BNG::compartment_id_t child3d: child_comp.children_compartments) {
          if (inside_of_compartments.count(child3d) == 1) {
            child_found = true;
            current = child3d;
            break;
          }
        }
      }
      else {
        if (inside_of_compartments.count(child) == 1) {
          child_found = true;
          current = child;
          break;
        }
      }
    }
  } while (child_found);

  counted_volume_index_to_compartment_id_cache[counted_volume_index] = current;
  return current;
}


void Partition::remove_from_known_vol_species(const species_id_t species_id) {
  if (known_vol_species.count(species_id) == 0) {
    return;
  }
  known_vol_species.erase(species_id);
}


void Partition::remove_reactant_class_usage(const BNG::reactant_class_id_t reactant_class_id) {
  volume_molecule_reactants_per_reactant_class.remove_reactant_sets_for_reactant_class(reactant_class_id);
}


void Partition::to_data_model(Json::Value& mcell, std::set<rgba_t>& used_colors) const {

  // there are two places in data model where geometry objects are
  // defined - in KEY_GEOMETRICAL_OBJECTS and KEY_MODEL_OBJECTS

  Json::Value& geometrical_objects = mcell[KEY_GEOMETRICAL_OBJECTS];
  Json::Value& object_list = geometrical_objects[KEY_OBJECT_LIST];
  if (object_list.isNull()) {
    object_list = Json::Value(Json::arrayValue);
  }

  for (const GeometryObject& g: geometry_objects) {
    Json::Value object;
    g.to_data_model_as_geometrical_object(*this, config, object, used_colors);
    object_list.append(object);
  }

  Json::Value& model_objects = mcell[KEY_MODEL_OBJECTS];
  DMUtils::add_version(model_objects, VER_DM_2018_01_11_1330);
  Json::Value& model_object_list = model_objects[KEY_MODEL_OBJECT_LIST];
  model_object_list = Json::Value(Json::arrayValue);

  for (const GeometryObject& g: geometry_objects) {
    Json::Value model_object;
    g.to_data_model_as_model_object(*this, model_object);
    model_object_list.append(model_object);
  }

  // mod_surf_regions
  Json::Value& modify_surface_regions = mcell[KEY_MODIFY_SURFACE_REGIONS];
  DMUtils::add_version(modify_surface_regions, VER_DM_2014_10_24_1638);
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


void Partition::print_periodic_stats() const {
  std::cout <<
      "Partition: molecules.size() = " << molecules.size() << "\n" <<
      "Partition: molecule_id_to_index_mapping.size() = " << molecule_id_to_index_mapping.size() << "\n" <<
      "Partition: next_molecule_id = " << next_molecule_id << "\n" <<
      "Partition: known_vol_species.size() = " << known_vol_species.size() << "\n";
}

} // namespace mcell
