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

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 */
#include <vector>

#include "logging.h"

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

#include "geometry_utils.inc"

using namespace std;

namespace MCell {

namespace CollisionUtil {

// ---------------------------------- subpartitions ----------------------------------


// move to collision utils or to partition
static vec3_t get_displacement_up_to_partition_boundary(
    const Partition& p,
    const vec3_t& pos,
    const vec3_t& displacement
    ) {

  vec3_t new_pos = pos + displacement;
  assert(p.in_this_partition(pos));
  assert(!p.in_this_partition(new_pos));

  // TODO: similar code is in collect_crossed_subparts, clean it up
  // direction of the displacement
  // TODO: get rid of the zero div check??

  vec3_t displacement_nonzero = displacement;
  debug_guard_zero_div(displacement_nonzero);
  ivec3_t dir_urb_direction = ivec3_t(glm::greaterThan(displacement_nonzero, vec3_t(0)));

  ivec3_t dir_urb_addend;
  dir_urb_addend.x = (dir_urb_direction.x == 0) ? -1 : 1;
  dir_urb_addend.y = (dir_urb_direction.y == 0) ? -1 : 1;
  dir_urb_addend.z = (dir_urb_direction.z == 0) ? -1 : 1;

  // position of edges in our direction
  vec3_t partition_edges =
      p.get_origin_corner()
      + vec3_t(dir_urb_direction) * p.config.partition_edge_length;

  vec3_t diff = partition_edges - pos;

  vec3_t displacement_rcp = 1.0/displacement_nonzero; // POSSIBLE ZERO DIV

  // time we hit a boundary
  float_t hit_time = 1;



  // first check whether we are not in fact touching one of the boundaries
  if (abs(diff.x) < EPS) {
    // only update the xyz subpartition index
    // curr_subpart_indices.x += dir_urb_addend.x;
    assert(false && "TODO");
  }
  else if (abs(diff.y) < EPS) {
    //curr_subpart_indices.y += dir_urb_addend.y;
    assert(false && "TODO");
  }
  else if (abs(diff.z) < EPS) {
    //curr_subpart_indices.z += dir_urb_addend.z;
    assert(false && "TODO");
  }
  else {
    // compute time for the next subpartition collision, let's assume that displacemnt
    // is our speed vector and the total time to travel is 1
    //
    // pos(time) = pos + displacement * time, therefore
    // time = (pos(time) - vm.v.pos) / displacement
    // =>
    // time_to_subpart_edge = (subpart_edge - vm.v.pos) / displacement_speed
    vec3_t coll_times = diff * displacement_rcp;
    assert(coll_times.x >= 0 && coll_times.y >= 0 && coll_times.z >= 0
      && "Subpartition 'edges' must be computed from direction, we cannot hit a subpart boundary that is behind us");

    // which of the times is the smallest? - i.e. which boundary we hit first
    if (coll_times.x >= 0 && coll_times.x < coll_times.y && coll_times.x <= coll_times.z) {
      // x
      hit_time = coll_times.x;
    }
    else if (coll_times.y >= 0 && coll_times.y <= coll_times.z) {
      // y
      hit_time = coll_times.y;
    }
    else if (coll_times.z >= 0) {
      // z
      hit_time = coll_times.z;
    }
    else {
      assert(false && "Collision time must not be negative");
    }
  }

  // there might be some floating point imprecisions, we want this value to be clearly in our partition,
  // so let's make the time a bit smaller
  // the displacement value is used only to find out which subpartitions we are crossing
  vec3_t new_displacement = displacement * (hit_time - EPS);
  assert(p.in_this_partition(pos + new_displacement));

  return new_displacement;
}


// This function checks if any of the neighboring subpartitions are within radius
// from pos and inserts them into crossed_subparition_indices
static void collect_neighboring_subparts(
    const Partition& p,
    const vec3_t& pos,
    const ivec3_t& subpart_indices,
    const float_t rx_radius,
    const float_t subpart_edge_len,
    subpart_indices_set_t& crossed_subpart_indices
) {
  vec3_t rel_pos = pos - p.get_origin_corner();

  // left (x)
  int x_dir_used = 0;
  float_t x_boundary = subpart_indices.x * subpart_edge_len;
  if (rel_pos.x - rx_radius < x_boundary) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x - 1, subpart_indices.y, subpart_indices.z));
    x_dir_used = -1;
  }
  // right (x)
  else if (rel_pos.x + rx_radius > x_boundary + subpart_edge_len) { // assuming that subpartitions are larger than radius
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x + 1, subpart_indices.y, subpart_indices.z));
    x_dir_used = +1;
  }

  // upper (y)
  int y_dir_used = 0;
  float_t y_boundary = subpart_indices.y * subpart_edge_len;
  if (rel_pos.y - rx_radius < y_boundary) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x, subpart_indices.y - 1, subpart_indices.z));
    y_dir_used = -1;
  }
  // right (y)
  else if (rel_pos.y + rx_radius > y_boundary + subpart_edge_len) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x, subpart_indices.y + 1, subpart_indices.z));
    y_dir_used = +1;
  }

  // front (z)
  int z_dir_used = 0;
  float_t z_boundary = subpart_indices.z * subpart_edge_len;
  if (rel_pos.z - rx_radius < z_boundary) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x, subpart_indices.y, subpart_indices.z - 1));
    z_dir_used = -1;
  }
  // back (z)
  else if (rel_pos.z + rx_radius > z_boundary + subpart_edge_len) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x, subpart_indices.y, subpart_indices.z + 1));
    z_dir_used = +1;
  }

  // we also have to count with movement in multiple dimensions
  // xy
  if (x_dir_used != 0 && y_dir_used != 0) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x + x_dir_used, subpart_indices.y + y_dir_used, subpart_indices.z));
  }

  // xz
  if (x_dir_used != 0 && z_dir_used != 0) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x + x_dir_used, subpart_indices.y, subpart_indices.z + z_dir_used));
  }

  // yz
  if (y_dir_used != 0 && z_dir_used != 0) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x, subpart_indices.y + y_dir_used, subpart_indices.z + z_dir_used));
  }

  // xyz
  if (x_dir_used != 0 && y_dir_used != 0 && z_dir_used != 0) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x + x_dir_used, subpart_indices.y + y_dir_used, subpart_indices.z + z_dir_used));
  }
}


// collect subpartition indices that we are crossing and that are within radius
// of vm that moves by displacement
static inline void __attribute__((always_inline)) collect_crossed_subparts(
  const Partition& p,
  const Molecule& vm, // molecule that we are diffusing
  const vec3_t& displacement,
  const float_t rx_radius,
  const float_t sp_edge_length,
  const bool collect_for_walls,
  SubpartIndicesVector& crossed_subparts_for_walls, // crossed subparts considered for wall collision
  subpart_indices_set_t& crossed_subparts_for_molecules // crossed subparts considered for molecule collisions
) {
  if (collect_for_walls) {
    assert(crossed_subparts_for_walls.empty());
  }
  assert(crossed_subparts_for_molecules.empty());

  // remember the starting subpartition
  if (collect_for_walls) {
    crossed_subparts_for_walls.push_back(vm.v.subpart_index);
  }
  crossed_subparts_for_molecules.insert(vm.v.subpart_index);

  // destination
  vec3_t dest_pos = vm.v.pos + displacement;

  // urb - upper, right, bottom
  vec3_t displacement_nonzero = displacement;
  debug_guard_zero_div(displacement_nonzero);
  ivec3_t dir_urb_direction = ivec3_t(glm::greaterThan(displacement_nonzero, vec3_t(0)));
  assert(dir_urb_direction.x == 0 || dir_urb_direction.x == 1);
  assert(dir_urb_direction.y == 0 || dir_urb_direction.y == 1);
  assert(dir_urb_direction.z == 0 || dir_urb_direction.z == 1);

  // get 3d indices of start and end subpartitions
  ivec3_t src_subpart_indices, dest_subpart_indices;
  assert(vm.v.subpart_index != SUBPART_INDEX_INVALID);
  p.get_subpart_3d_indices_from_index(vm.v.subpart_index, src_subpart_indices);
  p.get_subpart_3d_indices(dest_pos, dest_subpart_indices);


  // the radius must be adjusted so that it takes into account distance in the 3D space,
  // taking the radius in just one dimension is not sufficient as shown here:
  //
  //     a
  //    /|
  // __/_|__
  //  b  |x
  //
  // distance of the points a and b from x is higher than radius, but
  // distance of x to any point on the line might be less than radius
  //
  // mcell3 extends the bounding box in all dimensions and is overly pessimistic
  //
  // the (possibly) logic behind this is that the case has something to do with
  // the diameter of a square, computing precise length would be most probably worse perfrmance-wise
  // TODO: explain why this should work
  // a more efficient variant if definitely possible
  float_t rx_radius_for_neighbors = rx_radius * SQRT2;

  // first check what's around the starting point
  collect_neighboring_subparts(
      p, vm.v.pos, src_subpart_indices, rx_radius_for_neighbors, sp_edge_length,
      crossed_subparts_for_molecules
  );

  // collect subpartitions on the way by always finding the point where a subpartition boundary is hit
  // we must do it even when we are crossing just one subpartition because we might hit others while
  // moving along them
  if ( !glm::all( glm::equal(dest_subpart_indices, src_subpart_indices) ) ) {

    subpart_index_t dest_subpart_index = p.get_subpart_index_from_3d_indices(dest_subpart_indices);

    ivec3_t dir_urb_addend;
    dir_urb_addend.x = (dir_urb_direction.x == 0) ? -1 : 1;
    dir_urb_addend.y = (dir_urb_direction.y == 0) ? -1 : 1;
    dir_urb_addend.z = (dir_urb_direction.z == 0) ? -1 : 1;

    vec3_t curr_pos = vm.v.pos;
    ivec3_t curr_subpart_indices = src_subpart_indices;

    subpart_index_t curr_subpart_index;

    vec3_t displacement_rcp = 1.0/displacement_nonzero; // POSSIBLE ZERO DIV

    do {
      // subpartition edges
      // = origin + subparition index * length + is_urb * length

      // NOTE: some of these computation can be moved out of the loop, not sure if it will improve performance
      vec3_t sp_len_as_vec3 = vec3_t(sp_edge_length);
      vec3_t sp_edges =
          p.get_origin_corner()
          +  vec3_t(curr_subpart_indices) * sp_len_as_vec3 // llf edge
          + vec3_t(dir_urb_direction) * sp_len_as_vec3; // move if we go urb

      vec3_t diff = sp_edges - curr_pos;

      // first check whether we are not in fact touching one of the boundaries
      if (abs(diff.x) < EPS) {
        // only update the xyz subpartition index
        curr_subpart_indices.x += dir_urb_addend.x;\
        // in some cases, we can run out of partition
        // this means that we missed the destination partition which is fine since collection of subparts
        // is only an optimization, but still we must terminate
        if (!p.is_subpart_index_in_range(curr_subpart_indices.x)) {
          break;
        }
      }
      else if (abs(diff.y) < EPS) {
        curr_subpart_indices.y += dir_urb_addend.y;
        if (!p.is_subpart_index_in_range(curr_subpart_indices.y)) {
          break;
        }
      }
      else if (abs(diff.z) < EPS) {
        curr_subpart_indices.z += dir_urb_addend.z;
        if (!p.is_subpart_index_in_range(curr_subpart_indices.z)) {
          break;
        }
      }
      else {
        // compute time for the next subpartition collision, let's assume that displacemnt
        // is our speed vector and the total time to travel is 1
        //
        // pos(time) = pos + displacement * time, therefore
        // time = (pos(time) - vm.v.pos) / displacement
        // =>
        // time_to_subpart_edge = (subpart_edge - vm.v.pos) / displacement_speed
        vec3_t coll_times = diff * displacement_rcp;
        assert(coll_times.x >= 0 && coll_times.y >= 0 && coll_times.z >= 0
            && "Subpartition 'edges' must be computed from direction, we cannot hit a subpart boundary that is behind us");

        // which of the times is the smallest? - i.e. which boundary we hit first
        if (coll_times.x >= 0 && coll_times.x < coll_times.y && coll_times.x <= coll_times.z) {
          // new position on the edge of the subpartition
          curr_pos += displacement * coll_times.x;
          // and also update the xyz subpartition index
          curr_subpart_indices.x += dir_urb_addend.x;
          if (!p.is_subpart_index_in_range(curr_subpart_indices.x)) {
            break;
          }
        }
        else if (coll_times.y >= 0 && coll_times.y <= coll_times.z) {
          curr_pos += displacement * coll_times.y;
          curr_subpart_indices.y += dir_urb_addend.y;
          if (!p.is_subpart_index_in_range(curr_subpart_indices.y)) {
            break;
          }
        }
        else if (coll_times.z >= 0) {
          curr_pos += displacement * coll_times.z;
          curr_subpart_indices.z += dir_urb_addend.z;
          if (!p.is_subpart_index_in_range(curr_subpart_indices.z)) {
            break;
          }
        }
        else {
          break;
        }
      }

      curr_subpart_index = p.get_subpart_index_from_3d_indices(curr_subpart_indices);
      if (collect_for_walls) {
        crossed_subparts_for_walls.push_back(curr_subpart_index);
      }
      crossed_subparts_for_molecules.insert(curr_subpart_index);

      // also neighbors
      collect_neighboring_subparts(
          p, curr_pos, curr_subpart_indices, rx_radius_for_neighbors, sp_edge_length,
          crossed_subparts_for_molecules
      );

    } while (curr_subpart_index != dest_subpart_index);
  }

  // finally check also neighbors in destination
  collect_neighboring_subparts(
      p, dest_pos, dest_subpart_indices, rx_radius_for_neighbors, sp_edge_length,
      crossed_subparts_for_molecules
  );
}

static inline void get_corner_points_for_subpart_colection(
    const vec3_t& pos,
    const vec3_t& move,
    const float_t rx_radius,
    small_vector<vec3_t>& pts
) {
  assert(pts.size() >= 4);
  assert(!cmp_eq(move, 0));
  /* A) compute 4 points where the traces start
          p0--p1
          |\  /|
          | pos|
          |/  \|
          p2--p3
   */

   /* 1) first vector perpendicular to the move id computed this way
       v = (vx, vy, vz), move = (mx, my, mx)

       dot product:
       mx*vx + my*vy + mz*vz == 0

       then we set 2 values to be 1:
       mx*1 + my*1 + mz*vz == 0

       assuming that mz is not 0, then:
       v = (1, 1, (mx+my / mz))

       precision of the initial points is not so important
   */
  uint largest_dim = get_largest_abs_dim_index(move);
  vec3_t v;
  switch(largest_dim) {
    case 0:
      v = vec3_t( (move.y + move.z)/move.x, 1, 1);
      break;
    case 1:
      v = vec3_t( 1, (move.x + move.z)/move.y, 1);
      break;
    case 2:
      v = vec3_t( 1, 1, (move.x + move.y)/move.z);
      break;
    default:
      assert(false);
  }

  // 2) get a vector perpendicular to v and move
  vec3_t vp = cross(v, move);

  // 3) and compute the vectors v1-v4 of length rx_radius * SQRT2
  float_t ratio = 1/len3(v) * rx_radius * SQRT2;
  float_t ratiop = 1/len3(vp) * rx_radius * SQRT2;

  vec3_t v0 = v * vec3_t(ratio);
  vec3_t v1 = vp * vec3_t(ratiop);
  vec3_t v2 = -v1;
  vec3_t v3 = -v2;

  pts[0] = pos + v0;
  pts[1] = pos + v1;
  pts[2] = pos + v2;
  pts[3] = pos + v3;

#ifdef DEBUG_SUBPARTITIONS
  std::cout << "get_corner_points_for_subpart_colection: pos: " << pos
      << ", 0:" << pts[0] << ", 1:" << pts[1] << ", 2:" << pts[2] << ", 3:" << pts[3] << "\n";
#endif
}

// a potentially more efficient variant of collect_crossed_subparts
// 1 variant for molecules,
// another for walls (direct crosses, no radius)
static inline void __attribute__((always_inline)) collect_crossed_subparts2(
  const Partition& p,
  const Molecule& vm, // molecule that we are diffusing
  const vec3_t& displacement,
  const float_t rx_radius,
  const float_t sp_edge_length,
  const bool collect_for_walls,
  SubpartIndicesVector& crossed_subparts_for_walls, // crossed subparts considered for wall collision
  subpart_indices_set_t& crossed_subparts_for_molecules // crossed subparts considered for molecule collisions
) {
  if (collect_for_walls) {
    assert(crossed_subparts_for_walls.empty());
  }
  assert(crossed_subparts_for_molecules.empty());

  if (cmp_eq(displacement, 0)) {
    // we are practically not moving (add some nearby positions?)
    crossed_subparts_for_walls.push_back(vm.v.subpart_index);
    crossed_subparts_for_molecules.insert(vm.v.subpart_index);
    return;
  }

  const int NUM_CORNER_POINTS = 4;
  const int MOL_POS_POINT_INDEX = 4;

  // The idea here is to have 4 lines that represent a moving square instead of a circle,
  // the 4 lines are then checked for collisions with subpartitions. The rx radius must be smaller
  // than the subpart size, therefore the actual subpartition through which the ray from the
  // molecule's position goes are included automatically.
  //
  // However, when we are collecting walls, we want just the exct subpartitions and not neighbors.
  //
  small_vector<vec3_t> curr_positions;
  if (collect_for_walls) {
    curr_positions.resize(NUM_CORNER_POINTS + 1);
    curr_positions[MOL_POS_POINT_INDEX] = vm.v.pos;
  }
  else {
    curr_positions.resize(NUM_CORNER_POINTS);
  }

  uint num_points = curr_positions.size();

  small_vector<bool> reached_last_subpart;
  reached_last_subpart.resize(num_points, false);
  uint num_reached_last_subpart = 0;

  // set pts 0 - 3
  get_corner_points_for_subpart_colection(vm.v.pos, displacement, rx_radius, curr_positions);

  // compute how often a subpart plane is hit along the way for x, y, and z dimensions

  // remember the starting subpartitions
  small_vector<subpart_index_t> start_subpart_indices;
  start_subpart_indices.resize(num_points);
  for (uint i = 0; i < NUM_CORNER_POINTS; i++) {
    // the starting points may be out of our partition - this is not supported yet
    if (p.in_this_partition(curr_positions[i])) {
      start_subpart_indices[i] = p.get_subpart_index(curr_positions[i]);
    }
    else {
      std::cout << "Point " << i << ", curr_positions[i]: " << curr_positions[i] << "\n";
      mcell_error(
          "collect_crossed_subparts: Subpart search starting points are outside of the partition. Not supported yet.");
    }
  }
  if (collect_for_walls) {
    start_subpart_indices[MOL_POS_POINT_INDEX] = vm.v.subpart_index;
  }

  // urb - upper, right, bottom
  vec3_t displacement_nonzero = displacement;
  debug_guard_zero_div(displacement_nonzero);
  ivec3_t dir_urb_direction = ivec3_t(glm::greaterThan(displacement_nonzero, vec3_t(0)));
  assert(dir_urb_direction.x == 0 || dir_urb_direction.x == 1);
  assert(dir_urb_direction.y == 0 || dir_urb_direction.y == 1);
  assert(dir_urb_direction.z == 0 || dir_urb_direction.z == 1);

  // get 3d indices of start and end subpartitions
  small_vector<ivec3_t> curr_subpart_indices3d;
  small_vector<ivec3_t> dest_subpart_indices3d;
  small_vector<subpart_index_t> dest_subpart_indices;

  curr_subpart_indices3d.resize(num_points);
  dest_subpart_indices3d.resize(num_points);
  dest_subpart_indices.resize(num_points);

  for (uint i = 0; i < num_points; i++) {
    crossed_subparts_for_molecules.insert(start_subpart_indices[i]);

    p.get_subpart_3d_indices_from_index( start_subpart_indices[i], curr_subpart_indices3d[i]);

    vec3_t dest_pos = curr_positions[i] + displacement;

    // truncate if needed, does not happen often
    if (!p.in_this_partition(dest_pos)) {
      dest_pos = curr_positions[i] + get_displacement_up_to_partition_boundary(p, curr_positions[i], displacement);
    }

    p.get_subpart_3d_indices(dest_pos, dest_subpart_indices3d[i]);

    #ifdef DEBUG_SUBPARTITIONS
    DUMP_CONDITION4P(
      std::cout << "Point " << i << ", start subpart indices: " <<  curr_subpart_indices3d[i] << "\n";
      std::cout << "Point " << i << ", dest  subpart indices: " <<  dest_subpart_indices3d[i] << "\n";
    );
    #endif

    dest_subpart_indices[i] = p.get_subpart_index_from_3d_indices(dest_subpart_indices3d[i]);

    if (start_subpart_indices[i] == dest_subpart_indices[i]) {
      reached_last_subpart[i] = true;
      num_reached_last_subpart++;
    }
  }

  // for each dimension, process each trace and each plane
  // collect subpartitions on the way by always finding the point where a subpartition boundary is hit
  // we must do it even when we are crossing just one subpartition because we might hit others while
  // moving along them


  ivec3_t dir_urb_addend;
  dir_urb_addend.x = (dir_urb_direction.x == 0) ? -1 : 1;
  dir_urb_addend.y = (dir_urb_direction.y == 0) ? -1 : 1;
  dir_urb_addend.z = (dir_urb_direction.z == 0) ? -1 : 1;

  vec3_t displacement_rcp = 1.0/displacement_nonzero; // POSSIBLE ZERO DIV
  vec3_t sp_len_as_vec3 = vec3_t(sp_edge_length);

  while(num_reached_last_subpart != num_points) {

    for (uint i = 0; i < num_points; i++){
      if (reached_last_subpart[i]) {
        continue;
      }

      // subpartition edges
      // = origin + subparition index * length + is_urb * length
      vec3_t sp_edges =
        p.get_origin_corner()
        +  vec3_t(curr_subpart_indices3d[i]) * sp_len_as_vec3 // llf edge
        + vec3_t(dir_urb_direction) * sp_len_as_vec3; // move if we go urb

      vec3_t& curr_pos = curr_positions[i];
      vec3_t diff = sp_edges - curr_pos;

      int changed_subpart_index;

      // first check whether we are not in fact touching one of the boundaries
      if (abs(diff.x) < EPS) {
        // only update the xyz subpartition index
        curr_subpart_indices3d[i].x += dir_urb_addend.x;
        changed_subpart_index = curr_subpart_indices3d[i].x;
      }
      else if (abs(diff.y) < EPS) {
        curr_subpart_indices3d[i].y += dir_urb_addend.y;
        changed_subpart_index = curr_subpart_indices3d[i].y;
      }
      else if (abs(diff.z) < EPS) {
        curr_subpart_indices3d[i].z += dir_urb_addend.z;
        changed_subpart_index = curr_subpart_indices3d[i].z;
      }
      else {
        // compute time for the next subpartition collision, let's assume that displacemnt
        // is our speed vector and the total time to travel is 1
        //
        // pos(time) = pos + displacement * time, therefore
        // time = (pos(time) - vm.v.pos) / displacement
        // =>
        // time_to_subpart_edge = (subpart_edge - vm.v.pos) / displacement_speed
        vec3_t coll_times = diff * displacement_rcp;

        assert(coll_times.x >= 0 && coll_times.y >= 0 && coll_times.z >= 0
            && "Subpartition 'edges' for the main point must be computed from direction, we cannot hit a subpart boundary that is behind us");

        // which of the times is the smallest? - i.e. which boundary we hit first
        if (coll_times.x >= 0 && coll_times.x < coll_times.y && coll_times.x <= coll_times.z) {
          // new position on the edge of the subpartition
          curr_pos += displacement * coll_times.x;
          // and also update the xyz subpartition index
          curr_subpart_indices3d[i].x += dir_urb_addend.x;
          changed_subpart_index = curr_subpart_indices3d[i].x;
        }
        else if (coll_times.y >= 0 && coll_times.y <= coll_times.z) {
          curr_pos += displacement * coll_times.y;
          curr_subpart_indices3d[i].y += dir_urb_addend.y;
          changed_subpart_index = curr_subpart_indices3d[i].y;
        }
        else if (coll_times.z >= 0) {
          curr_pos += displacement * coll_times.z;
          curr_subpart_indices3d[i].z += dir_urb_addend.z;
          changed_subpart_index = curr_subpart_indices3d[i].z;
        }
        else {
          assert(false);
          reached_last_subpart[i] = true;
          num_reached_last_subpart++;
          continue;
        }
      }

      // in some cases, we can run out of partition
      // this means that we missed the destination partition which is fine since collection of subparts
      // is only an optimization, but still we must terminate this search
      if (!p.is_subpart_index_in_range(changed_subpart_index)) {
        reached_last_subpart[i] = true;
        num_reached_last_subpart++;
        continue;
      }


      subpart_index_t curr_subpart_index = p.get_subpart_index_from_3d_indices(curr_subpart_indices3d[i]);
      #ifdef DEBUG_SUBPARTITIONS
      DUMP_CONDITION4P(
        std::cout << "Point " << i << ", new subpart indices: " <<  curr_subpart_indices3d[i] << "\n";
      );
      #endif

      if (collect_for_walls && i == MOL_POS_POINT_INDEX) {
        crossed_subparts_for_walls.push_back(curr_subpart_index);
      }
      crossed_subparts_for_molecules.insert(curr_subpart_index);

      if (curr_subpart_index == dest_subpart_indices[i]) {
        // do we have to reach the last subpart at the same time?
        // ... expecting that yes
        reached_last_subpart[i] = true;
        num_reached_last_subpart++;
      }

    } // for
  } // while
}


#if 0
// one more option is to use point-line distance and
// somehow go along the trace and check neighboring points
//
static inline void __attribute__((always_inline)) collect_crossed_subparts3(
  const Partition& p,
  const Molecule& vm, // molecule that we are diffusing
  const vec3_t& displacement,
  const float_t rx_radius,
  const float_t sp_edge_length,
  const bool collect_for_walls,
  SubpartIndicesVector& crossed_subparts_for_walls, // crossed subparts considered for wall collision
  subpart_indices_set_t& crossed_subparts_for_molecules, // crossed subparts considered for molecule collisions
  subpart_index_t& last_subpart_index
) {
  // starting point


  // reuse collide_mol


  // 1) check corners
}
#endif

#if 0
// last option is to fall back to the mcell3 implementation (maybe try this first of optimization is needed)
static inline void __attribute__((always_inline)) collect_crossed_subparts4(
#endif


// ---------------------------------- molecule collisions ----------------------------------

// check whether diffused_vm molecule collision that moves by displacement can collide
// with colliding_vm; returns true if there can be a collision and returns relative collision
// time and relative position
static bool collide_mol(
    const Molecule& diffused_vm,
    const vec3_t& displacement,
    const Molecule& colliding_vm,
    const float_t rx_radius_3d,
    float_t& rel_collision_time,
    vec3_t& rel_collision_pos
) {

  const vec3_t& pos = colliding_vm.v.pos; /* Position of target molecule */
  vec3_t dir = pos - diffused_vm.v.pos;  /* From starting point of moving molecule to target */

  float_t d = glm::dot((glm_vec3_t)dir, (glm_vec3_t)displacement);        /* Dot product of movement vector and vector to target */

  /* Miss the molecule if it's behind us */
  if (d < 0) {
    return false;
  }

  float_t movelen2 = glm::dot((glm_vec3_t)displacement, (glm_vec3_t)displacement); /* Square of distance the moving molecule travels */

  /* check whether the test molecule is further than the displacement. */
  if (d > movelen2) {
    return false;
  }

  /* check whether the moving molecule will miss interaction disk of the
     test molecule.*/
  float_t dirlen2 = glm::dot((glm_vec3_t)dir, (glm_vec3_t)dir);
  float_t sigma2 = rx_radius_3d * rx_radius_3d;   /* Square of interaction radius */
  if (movelen2 * dirlen2 - d * d > movelen2 * sigma2) {
    return false;
  }

  /* reject collisions with itself */
  if (diffused_vm.id == colliding_vm.id) {
    return false;
  }

  /* defunct - not probable */
  if (colliding_vm.is_defunct()) {
    return false;
  }

  rel_collision_time = d / movelen2;

  rel_collision_pos = diffused_vm.v.pos + rel_collision_time * displacement;
  return true;
}


// body of the collision detection loop
// made into separate function to be possibly able to make some optimizations over it in the future
static void collide_mol_loop_body(
    Partition& p,
    const Molecule& vm,
    const molecule_id_t colliding_vm_id,
    const vec3_t& remaining_displacement,
    const float_t radius,
    collision_vector_t& molecule_collisions
) {

  Molecule& colliding_vm = p.get_m(colliding_vm_id);

  // we would like to compute everything that's needed just once
  float_t time;
  vec3_t position;
  // collide_mol must be inlined because many things are computed all over there
  if (collide_mol(vm, remaining_displacement, colliding_vm, radius, time, position)) {
    const RxnClass* rxn_class = p.all_reactions.get_specific_reaction_class(vm, colliding_vm);
    assert(rxn_class != nullptr);
    molecule_collisions.push_back(
        Collision(CollisionType::VOLMOL_VOLMOL, &p, vm.id, time, position, colliding_vm.id, rxn_class)
    );
  }
}


// ---------------------------------- wall collisions ----------------------------------


/***************************************************************************
jump_away_line:
  In: starting coordinate
      vector we were going to move along and need to change
      fraction of way we moved before noticing we were hitting a edge
      location of the first vertex of the edge
      location of the second vertex of the edge
      normal vector to the surface containing our edge
  Out: No return value.  Movement vector is slightly changed.
***************************************************************************/
static void jump_away_line(
    const vec3_t& p,
    const float_t k, const vec3_t& A, const vec3_t& B, const vec3_t& n, rng_state& rng,
    vec3_t& v /*inout*/
) {
  vec3_t e, f;
  float_t le_1, tiny;

  e = B - A;
  float_t elen2 = glm::dot((glm_vec3_t)e, (glm_vec3_t)e);
  le_1 = 1.0 / sqrt(elen2);

  e = e * vec3_t(le_1);

  f.x = n.y * e.z - n.z * e.y;
  f.y = n.z * e.x - n.x * e.z;
  f.z = n.x * e.y - n.y * e.x;

  tiny = EPS * (MCell::abs_max_2vec(p, v) + 1.0) /
         (k * max3d(glm::abs((glm_vec3_t)f)));
  if ((rng_uint(&rng) & 1) == 0) {
    tiny = -tiny;
  }
  v.x -= tiny * f.x;
  v.y -= tiny * f.y;
  v.z -= tiny * f.z;
}


/***************************************************************************
collide_wall:
  In: point: starting coordinate
      move: vector to move along
      face: wall we're checking for a collision
      t: double to store time of collision
      hitpt: vector to store the location of the collision
      update_move: flag to signal whether we should modify the movement vector
        in an ambiguous case (i.e. if we hit an edge or corner); if not, any
        ambiguous cases are treated as a miss.
  Out: Integer value indicating what happened
         COLLIDE_MISS  missed
         COLLIDE_FRONT hit the front face (face normal points out of)
         COLLIDE_BACK  hit the back face
         COLLIDE_REDO  hit an edge and modified movement vector; redo
  Note: t and/or hitpt may be modified even if there is no collision
        Not highly optimized yet.  May want to project to Cartesian
        coordinates for speed (as MCell2 did, and Rex implemented
        in pre-40308 backups in vol_utils.c).  When reflecting, use
        the value of t returned, not hitpt (reflections happen slightly
        early to avoid rounding errors with superimposed planes).
***************************************************************************/
static CollisionType collide_wall(
    const Partition& p,
    const vec3_t& pos, const Wall& face,
    rng_state &rng,
    const bool update_move,
    vec3_t& move,
    float_t& collision_time, vec3_t& collision_pos
) {
  p.stats.inc_ray_polygon_tests();

  float_t dp, dv, dd;
  float_t d_eps;

  const vec3_t& normal = face.normal;

  dp = dot(normal, pos);
  dv = dot(normal, move);
  dd = dp - face.distance_to_origin;

  if (dd > 0) {
    d_eps = EPS;
    if (dd < d_eps)
      d_eps = 0.5 * dd;

    /* Start & end above plane */
    if (dd + dv > d_eps) {
      return CollisionType::WALL_MISS;
    }
  }
  else {
    d_eps = -EPS;
    if (dd > d_eps)
      d_eps = 0.5 * dd;

    /* Start & end below plane */
    if (dd < 0.0 && dd + dv < d_eps) {
      return CollisionType::WALL_MISS;
    }
  }

  float_t a;

  if (dd == 0.0) {
    /* Start beside plane, end above or below */
    if (dv != 0.0)
      return CollisionType::WALL_MISS;

    // in case that the trajectory is parallel to the wall?
    // update the displacement a bit
    if (update_move) {
      a = (MCell::abs_max_2vec(pos, move) + 1.0) * EPS_C;
      if ((rng_uint(&rng) & 1) == 0)
        a = -a;
      if (dd == 0.0) {
        move = move - vec3_t(a) * normal;
      }
      else {
        move = move * vec3_t(1.0 - a);
      }
      return CollisionType::WALL_REDO;
    }
    else {
      return CollisionType::WALL_MISS;
    }
  }

  a = 1.0 / dv;
  a *= -dd; /* Time we actually hit */
  collision_time = a;

  collision_pos = pos + a * move;

  vec3_t face_vert[VERTICES_IN_TRIANGLE];
  face_vert[0] = p.get_geometry_vertex(face.vertex_indices[0]);
  face_vert[1] = p.get_geometry_vertex(face.vertex_indices[1]);
  face_vert[2] = p.get_geometry_vertex(face.vertex_indices[2]);

  vec3_t local = collision_pos - face_vert[0];

  float_t b = dot(local, face.unit_u);
  float_t c = dot(local, face.unit_v);

  float_t f;
  if (face.uv_vert2.v < 0.0) {
    c = -c;
    f = -face.uv_vert2.v;
  }
  else {
    f = face.uv_vert2.v;
  }

  if (c > 0) {
    float_t g, h;
    g = b * f;
    h = c * face.uv_vert2.u;
    if (g > h) {
      if (c * face.uv_vert1_u + g < h + face.uv_vert1_u * face.uv_vert2.v) {
        if (dv > 0) {
          return CollisionType::WALL_BACK;
        }
        else {
          return CollisionType::WALL_FRONT;
        }
      }
      else if ((!distinguishable_f(
          c * face.uv_vert1_u + g,
          h + face.uv_vert1_u * face.uv_vert2.v,
          EPS_C))) {
        if (update_move) {
          jump_away_line(pos, a, face_vert[1], face_vert[2], face.normal, rng, move);
          return CollisionType::WALL_REDO;
        }
        else {
          return CollisionType::WALL_MISS;
        }
      }
      else {
        return CollisionType::WALL_MISS;
      }
    }
    else if (!distinguishable_f(g, h, EPS)) {
      if (update_move) {
        jump_away_line(pos, a, face_vert[2], face_vert[0], face.normal, rng, move);
        return CollisionType::WALL_REDO;
      }
      else {
        return CollisionType::WALL_MISS;
      }
    }
    else {
      return CollisionType::WALL_MISS;
    }
  }
  else if (!distinguishable_f(c, 0, EPS)) /* Hit first edge! */
  {
    if (update_move) {
      jump_away_line(pos, a, face_vert[0], face_vert[1], face.normal, rng, move);
      return CollisionType::WALL_REDO;
    }
    else {
      return CollisionType::WALL_MISS;
    }
  }
  else {
    return CollisionType::WALL_MISS;
  }
}


static bool is_immediate_collision(const float_t time) {
  return time < EPS;
}

// called only from ray_trace_vol
static void collect_wall_collisions(
    Partition& p,
    const Molecule& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const subpart_index_t subpart_index,
    const wall_index_t previous_reflected_wall,
    rng_state& rng,
    // displacement can be changed in case we needed to 'REDO' the collision, also,
    // if there was a hit, is changed to the closest displacement
    vec3_t& displacement,
    vec3_t& displacement_up_to_wall_collision, // overwritten only when there is a wall collision
    collision_vector_t& wall_collisions
) {

  // remember which was the closest hit to update displacement
  float_t closest_hit_time = TIME_FOREVER;
  vec3_t closest_hit_pos;

  // check each wall in this subpartition
  const uint_set<wall_index_t>& wall_indices = p.get_subpart_wall_indices(subpart_index);

  for (auto it = wall_indices.begin(); it != wall_indices.end(); it++) {
    wall_index_t wall_index = *it;

    if (wall_index == previous_reflected_wall)
      continue;

    const Wall& w = p.get_wall(wall_index);
    float_t collision_time;
    vec3_t collision_pos;

#ifdef DEBUG_COLLISIONS_WALL_EXTRA
    SimulationStats* world = &p.stats;
    // just faking the name for the dump condition macro - FIXME - use better name
    DUMP_CONDITION4(
        std::cout << "Checking wall:\n";
        w.dump(p, "", true);
    );
#endif

    CollisionType collision_type =
        collide_wall(p, vm.v.pos, w, rng, true, displacement, collision_time, collision_pos);

#ifdef DEBUG_COLLISIONS_WALL_EXTRA
    DUMP_CONDITION4(
        if (collision_type == CollisionType::WALL_REDO || collision_type == CollisionType::WALL_FRONT || collision_type == CollisionType::WALL_BACK) {
          cout << "Collide wall: vm pos: " << vm.v.pos << ", displacement: " << displacement << "\n";
          w.dump(p, "", true);
          cout << "collision time: " << collision_time << ", collision pos: " << collision_pos << "\n";
        }
    );
#endif

    if (collision_type == CollisionType::WALL_REDO) {
      // molecule was 'jumped' and we need to run collision detection over all walls again
      wall_collisions.clear();
      it = wall_indices.begin();
    }
    else if (collision_type != CollisionType::WALL_MISS) {
      // remember the closest hit position, but ignore immediate
      if (collision_time < closest_hit_time && !is_immediate_collision(collision_time)) {
        closest_hit_pos = collision_pos;
      }

      p.stats.inc_ray_polygon_colls();
      wall_collisions.push_back(
          Collision(collision_type, &p, vm.id, collision_time, collision_pos, wall_index)
      );
    }
  }

  if (!wall_collisions.empty()) {
    // update displacement up to the wall collision
    displacement_up_to_wall_collision = closest_hit_pos - vm.v.pos;
  }
}


// returns true if there was a collision with a wall,
// HIT REDO is not supported
static bool collide_wall_test(
    const Partition& p,
    const vec3_t& pos, const Wall& face,
    const vec3_t& move
) {
  rng_state unused_rng_state;
  float_t ignored_collision_time;
  vec3_t ignored_collision_pos;

  vec3_t tmp_move = move;
  CollisionType res;
#ifndef NDEBUG
  uint orig_randcnt = unused_rng_state.randcnt;
#endif
  res = CollisionUtil::collide_wall(
      p, pos, face, unused_rng_state, false, tmp_move,
      ignored_collision_time, ignored_collision_pos
  );
#ifndef NDEBUG
  assert(orig_randcnt == unused_rng_state.randcnt && "collide_wall_test should not trigger usage of rng");
#endif


  switch (res) {
    case CollisionType::WALL_MISS:
      return false;
    case CollisionType::WALL_FRONT:
    case CollisionType::WALL_BACK:
      #ifdef DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS
        cout << "# Detecting collision for molecule with original pos " << pos <<
          " at " << ignored_collision_pos << " with wall:\n";
        face.dump(p, "", true);
      #endif
      return true;
      break;
    case CollisionType::WALL_REDO:
      mcell_error("Collision REDO is not handled yet in dynamic vertices.");
      return false;
    default:
      assert(false);
      return false;
  }
}

/**

 Detection of whether two lines on the same plane (not checked) cross.

 From http://mathworld.wolfram.com/Line-LineIntersection.html

  The intersection of two lines containing the points
    x_1=(x_1,y_1,z_1) and x_2=(x_2,y_2,z_2), and
    x_3=(x_3,y_3,z_3) and x_4=(x_4,y_4,z_4), respectively,

  can also be found directly by simultaneously solving
  x = x_1+(x_2-x_1)s     (17)
  x = x_3+(x_4-x_3)t     (18)

 together with the condition that the four points be coplanar (i.e., the lines are not skew),

  his set of equations can be solved for s to yield

   s=((cxb)·(axb))/(|axb|^2)   (20)

  where
  a = x_2-x_1   (21)
  b = x_4-x_3   (22)
  c = x_3-x_1   (23)

  the point of intersection can then be immediately found by plugging back in for s to obtain
    x=x_1+a((cxb)·(axb))/(|axb|^2) (24)
*/
static bool collide_line_and_line_test(
    const vec3_t& e, const vec3_t& f, const vec3_t& vfe /* == f - e*/,
    const vec3_t& o, const vec3_t& p
) {
  // rename arguments to get the names used in comment
  const vec3_t& x1 = e;
  const vec3_t& x2 = f;
  const vec3_t& x3 = o;
  const vec3_t& x4 = p;

  const vec3_t& a = vfe;
  vec3_t b = x4 - x3;
  vec3_t c = x3 - x1;

  vec3_t a_x_b = cross(a, b);
  vec3_t c_x_b = cross(c, b);

  float_t len_squared_a_x_b = len3_squared(a_x_b);
  if (cmp_eq(len_squared_a_x_b, 0)) {
    // s would be too large if the divisor is close to 0
    return false;
  }

  // s=((cxb)·(axb))/(|axb|^2)
  float_t s = dot(c_x_b, a_x_b) / len_squared_a_x_b;

  // check whether we are in segment e-f
  if (s < 0.0 || s > 1.0) {
    return false;
  }

  // check whether we are in segment o-p
  // a' -> x4 - x3 = b
  // b' -> x2 - x1 = a
  // c' -> x1 - x3 = d
  vec3_t d = x1 - x3;

  // t = ((c'xb')·(a'xb'))/(|a'xb'|^2)
  // t = ((dxa)·(bxa))/(|bxa|^2)
  vec3_t d_x_a = cross(d, a);
  vec3_t b_x_a = cross(b, a);

  float_t len_squared_b_x_a = len3_squared(b_x_a);
  if (cmp_eq(len_squared_b_x_a, 0)) {
    assert(false && "Should not happen anymore");
    // t would be too large if the divisor is close to 0
    return false;
  }

  float_t t = dot(d_x_a, b_x_a) / len_squared_b_x_a;

  if (t < 0.0 || t > 1.0) {
    return false;
  }

  return true;
}


/**
 collide_moving_line_and_static_line_test

 Detection of ray casts for situations there at least 2 of the vertices of a wall were moved:

    f
     \
 k..->\.o..->..l
 |     \|      |
 |      x      |
 |      |\     |
 m..->..p.\->..n
           \
            e


 Original edge of the moved wall: k m
 New edge of the moved wall: l n
 Percentage of the path traveled by the wall: t  (must be in range 0..1)
 Position of our molecule: e
 Destination position after ray trace: f
 Point where the ray hits the moving edge: x

 Point on the path from k to l: o = k + t*(l-k)
 Point on the path from m to n: p = m + t*(n - m)

 Normal vector of plane Po defined by e, f, o: no = cross(e-f, o-f)
 A vector on plane Pp defined by e, f, p: v = p - e

 First step in determining x is to find value of t that then tells us the coordinates of o and p.

 For this, we must find value of t where Po and Pp are the same.

 Value t can be computed from dot(no, v) == 0
 Analytical solutions are too complicated (hundreds of multiplications), therefore a
 numerical solution is used instead. Allowed range of t == 0 .. 1 allows us to quickly
 throw out cases where our ray does not cross the moving edge.

 f(t) = dot(
    cross((e - f), (-f + k + (-k + l)t) ),
    (-e + m + (-m + n)t)
 )

 f'(t) =
   dot(
    cross((e - f), (-k + l)),
    (-e + m - m t + n t)
   ) +
   dot(
    cross((e - f), (-f + k - k t + l t),
    (-m + n)
   )


 We need to find t such that f(t) == 0.

 Mathematica code used to obtain the equations above:

 o = k + t*(l - k)
 p = m + t*(n - m)

 no = Cross[e - f, o - f]
 fun[] := Dot[no, p - e]
 dfun[] := Simplify[D[fun[], t]]

*/

namespace Local {

static float_t compute_f(
    const vec3_t& e, const vec3_t& f,
    const vec3_t& k, const vec3_t& l,
    const vec3_t& m, const vec3_t& n,
    const float_t t
) {
  return
      dot(
          cross((e - f), (-f + k + (-k + l) * vec3_t(t)) ),
          (-e + m + (-m + n) * vec3_t(t))
      );
}

static float_t compute_df(
    const vec3_t& e, const vec3_t& f,
    const vec3_t& k, const vec3_t& l,
    const vec3_t& m, const vec3_t& n,
    const float_t t
) {
  return
      dot(
          cross(
              (e - f),
              (-k + l)
          ),
          (-e + m + (- m + n) * vec3_t(t))
      )
      +
      dot(
          cross(
              (e - f),
              (-f + k + (- k + l) * vec3_t(t))
          ),
          (-m + n)
      );
}

} // namespace Local


// returns true if plane that contains all points efop was found and also crosses line
// segments kl and mn
// sets points o and p
// returns false is such plane does not exist
static bool find_plane_crossing_efop(
    const vec3_t& e, const vec3_t& f, const vec3_t& move /* f = e + move */,
    const vec3_t& k, const vec3_t& l,
    const vec3_t& m, const vec3_t& n,
    vec3_t& o, vec3_t& p) {

  assert(cmp_eq(f, e + move));

  // first determine t, it must be in the range 0..1
  // f is continuous and (should?) have just one solution

  // starting from 0
  float_t t = 0;
  float_t t_previous = GIGANTIC4;

  // use Newton's Method to find solution to 't' (https://en.wikipedia.org/wiki/Newton's_method)
  // using high precision, we should converge quickly
  bool dft_is_zero = false;
  bool t_out_of_range = false;

  float_t ft;
  while (!cmp_eq(t, t_previous, EPS)) {
    ft = Local::compute_f(e, f, k, l, m, n, t);
    float_t dft = Local::compute_df(e, f, k, l, m, n, t);

    if (cmp_eq(dft, 0, EPS)) {
      dft_is_zero = true;
      break;
    }

    t_previous = t;
    t = t_previous - ft / dft;

    if (t < 0.0 || t > 1.0) {
      t_out_of_range = true;
      break;
    }
  }

  if (t_out_of_range) {
    return false;
  }

  bool ft_is_zero = (cmp_eq(ft, 0, EPS));

  if (dft_is_zero && !ft_is_zero) {
    // we found a minimum that is not zero however
    // there is no plane that could connect all 4 points efop
    return false;
  }

  // ok, we found our plane and we know the value of 't',
  // we can continue with figuring out out whether we really cross
  // the object defined by the moving edge
  o = k + (l - k) * vec3_t(t);
  p = m + (n - m) * vec3_t(t);
  return true;
}

static bool collide_moving_line_and_static_line_test(
    const Partition& part,
    const vec3_t& e, const vec3_t& move /* f = e + move */,
    const vec3_t& k, const vec3_t& l,
    const vec3_t& m, const vec3_t& n,
    const bool edge_moved,
    const Wall* wall_if_edge_defines_triangle /* might be nullptr */
) {
  assert(e != e + move && "The static line must not be a point");
  assert(k != m && "The source moving line must not be a point");
  assert(l != n && "The destination moving line must not be a point");

  if (!edge_moved) {
    return false;
  }

  assert(!(k == l && m == n) && "The line must move");

  // 1) if klmn creates a triangle, one of the arguments kmml or kmn is not nullptr
  if (wall_if_edge_defines_triangle) {
    return collide_wall_test(part, e, *wall_if_edge_defines_triangle, move);
  }

  // 2) if klmn is more complex and possibly these points do not lie on the same plane, therefore simple triangle
  //    ray trace cannot be used
  vec3_t o, p;
  vec3_t f = e + move;
  bool found = find_plane_crossing_efop(e, f, move, k, l, m, n, o, p);
  if (!found) {
    return false;
  }

  // now we must find whether the line segments ef and op intersect
  bool collides = collide_line_and_line_test(e, f, move, o, p);


#ifdef DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS
  if (collides) {
    cout << "# Detecting collision for molecule with original pos " << e <<
      " with a moving edge with points klmn " << k << ", " << l  << ", " << m  << ", " << n << "\n";
  }
#endif

  return collides;
}


// ---------------------------------- other detections ----------------------------------

/*************************************************************************
 is_point_inside_release_region:
    Simplified against MCell3 variant - we are dealing with just one
    region

    Does not deal with subpartitions yet

    Expects that the region is an enclosed volume

*************************************************************************/
static bool is_point_inside_region(const Partition& p, const vec3_t pos, const Region& reg) {

  // cast ray along the whole partition
  vec3_t move(0, 0, p.config.partition_edge_length);

  // and count intersects
  int num_hits = 0;
  for (const auto& wall_edges_pair: reg.walls_and_edges) {

    // check collision with each wall
    const Wall& w = p.get_wall(wall_edges_pair.first);
    bool collides = CollisionUtil::collide_wall_test(p, pos, w, move);
    if (collides) {
      num_hits++;
    }
  }

  // odd number of hits means that we are inside
  return num_hits % 2 == 1;
}

// ---------------------------------- reflections ----------------------------------

/******************************************************************************
 *
 * the reflect_or_periodic_bc helper function is used in diffuse_3D to handle
 * either reflections or periodic boundary conditions for a diffusing molecule
 * encountering a wall
 *
 * Return values:
 *
 *  0 : indicates that the molecule reflected off a wall
 *  1 : indicates that the molecule hit a periodic box and was moved to a
 *      position in the neighboring image
 *
 ******************************************************************************/
static int reflect_or_periodic_bc(
    const Partition& p,
    const Collision& collision,
    Molecule& vm, // moves vm to the reflection point
    vec3_t& displacement,
    float_t& remaining_time_step, // same as t_steps
    wall_index_t& reflected_wall_index
) {

  const Wall& w = p.get_wall(collision.colliding_wall_index);
  wall_index_t reflect_w = collision.colliding_wall_index;
  float_t reflect_t = collision.time;

  /* Update molecule location to the point of reflection (originally in register_hits) */
  vm.v.pos = collision.pos;

  /* Reduce our remaining available time. */
  remaining_time_step *= (1.0 - reflect_t);

  reflected_wall_index = reflect_w;

  float_t reflect_factor = -2.0 * glm::dot((glm_vec3_t)displacement, (glm_vec3_t)w.normal);

  // Set displacement for remainder of step length
  // No PBCs or non-traditional PBCs
  displacement = (displacement + vec3_t(reflect_factor) * w.normal) * vec3_t(1.0 - reflect_t);

  return 0;
}


} // namespace coll_util

} // namespace mcell
