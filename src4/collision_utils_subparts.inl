/******************************************************************************
 *
 * Copyright (C) 2019,2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_COLLISION_UTILS_SUBPARTS_INC_
#define SRC4_COLLISION_UTILS_SUBPARTS_INC_

#include <vector>

#include "logging.h"

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

#include "geometry_utils.inl"

using namespace std;


namespace MCell {

namespace CollisionUtils {

// This function checks if any of the neighboring subpartitions are within radius
// from pos and inserts them into crossed_subparition_indices
static inline void __attribute__((always_inline)) collect_neighboring_subparts(
    const Partition& p,
    const Vec3& pos,
    const IVec3& subpart_indices,
    const pos_t rxn_radius,
    const pos_t subpart_edge_len,
    SubpartIndicesSet& crossed_subpart_indices
) {
  const pos_t part_len = p.config.partition_edge_length;

  Vec3 rel_pos = pos - p.get_origin_corner();

  Vec3 rxn_radius3(rxn_radius);
  Vec3 rel_pos_plus_radius = rel_pos + rxn_radius3;
  Vec3 rel_pos_minus_radius = rel_pos - rxn_radius3;

  Vec3 boundary = Vec3(subpart_indices) * Vec3(subpart_edge_len);

  // left (x)
  int x_dir_used = 0;
  if (rel_pos_minus_radius.x < boundary.x && rel_pos_minus_radius.x > 0.0) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x - 1, subpart_indices.y, subpart_indices.z));
    x_dir_used = -1;
  }
  // right (x)
  else if (rel_pos_plus_radius.x > boundary.x + subpart_edge_len && rel_pos_plus_radius.x < part_len) { // assuming that subpartitions are larger than radius
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x + 1, subpart_indices.y, subpart_indices.z));
    x_dir_used = +1;
  }

  // upper (y)
  int y_dir_used = 0;
  if (rel_pos_minus_radius.y < boundary.y && rel_pos_minus_radius.y > 0.0) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x, subpart_indices.y - 1, subpart_indices.z));
    y_dir_used = -1;
  }
  // right (y)
  else if (rel_pos_plus_radius.y > boundary.y + subpart_edge_len && rel_pos_plus_radius.y < part_len) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x, subpart_indices.y + 1, subpart_indices.z));
    y_dir_used = +1;
  }

  // front (z)
  int z_dir_used = 0;
  if (rel_pos_minus_radius.z < boundary.z && rel_pos_minus_radius.z > 0.0) {
    crossed_subpart_indices.insert(
        p.get_subpart_index_from_3d_indices(subpart_indices.x, subpart_indices.y, subpart_indices.z - 1));
    z_dir_used = -1;
  }
  // back (z)
  else if (rel_pos_plus_radius.z > boundary.z + subpart_edge_len && rel_pos_plus_radius.z < part_len) {
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
static inline subpart_index_t __attribute__((always_inline)) collect_crossed_subparts(
  const Partition& p,
  const Molecule& vm, // molecule that we are diffusing
  const Vec3& displacement,
  const pos_t rxn_radius,
  const pos_t sp_edge_length,
  const bool collect_for_molecules,
  const bool collect_for_walls,
  SubpartIndicesVector& crossed_subparts_for_walls, // crossed subparts considered for wall collision
  SubpartIndicesSet& crossed_subparts_for_molecules // crossed subparts considered for molecule collisions
) {
  subpart_index_t current_subpart_index = vm.v.subpart_index;

  if (collect_for_walls) {
    assert(crossed_subparts_for_walls.empty());
  }
  assert(crossed_subparts_for_molecules.empty());

  // remember the starting subpartition
  if (collect_for_walls) {
    crossed_subparts_for_walls.push_back(current_subpart_index);
  }
  if (collect_for_molecules) {
    crossed_subparts_for_molecules.insert(current_subpart_index);
  }

  // destination
  Vec3 dest_pos = vm.v.pos + displacement;

  // urb - upper, right, bottom
  Vec3 displacement_nonzero = displacement;
  guard_zero_div(displacement_nonzero);
  IVec3 dir_urb_direction = IVec3(glm::greaterThan(displacement_nonzero, Vec3(0)));
  assert(dir_urb_direction.x == 0 || dir_urb_direction.x == 1);
  assert(dir_urb_direction.y == 0 || dir_urb_direction.y == 1);
  assert(dir_urb_direction.z == 0 || dir_urb_direction.z == 1);

  // get 3d indices of start and end subpartitions
  IVec3 src_subpart_indices, dest_subpart_indices;
  assert(current_subpart_index != SUBPART_INDEX_INVALID);
  p.get_subpart_3d_indices_from_index(current_subpart_index, src_subpart_indices);
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
  pos_t rxn_radius_for_neighbors = rxn_radius * POS_SQRT2;

  // first check what's around the starting point
  if (collect_for_molecules && p.config.use_expanded_list) {
    collect_neighboring_subparts(
        p, vm.v.pos, src_subpart_indices, rxn_radius_for_neighbors, sp_edge_length,
        crossed_subparts_for_molecules
    );
  }

  // collect subpartitions on the way by always finding the point where a subpartition boundary is hit
  // we must do it even when we are crossing just one subpartition because we might hit others while
  // moving along them
  subpart_index_t dest_subpart_index = p.get_subpart_index_from_3d_indices_allow_outside(dest_subpart_indices);
  if ( current_subpart_index != dest_subpart_index) {

    IVec3 dir_urb_addend(
        (dir_urb_direction.x == 0) ? -1 : 1,
        (dir_urb_direction.y == 0) ? -1 : 1,
        (dir_urb_direction.z == 0) ? -1 : 1
    );

    Vec3 curr_pos = vm.v.pos;
    IVec3 curr_subpart_indices = src_subpart_indices;

    subpart_index_t curr_subpart_index;

    Vec3 displacement_rcp = Vec3(1.0)/displacement_nonzero;

    do {
      // subpartition edges
      // = origin + subparition index * length + is_urb * length

      // NOTE: some of these computation can be moved out of the loop, not sure if it will improve performance
      Vec3 sp_len_as_vec3(sp_edge_length);
      Vec3 sp_edges =
          p.get_origin_corner()
          + Vec3(curr_subpart_indices) * sp_len_as_vec3 // llf edge
          + Vec3(dir_urb_direction) * sp_len_as_vec3; // move if we go urb

      Vec3 diff = sp_edges - curr_pos;

#if POS_T_BYTES == 4
      pos_t dist = POS_SQRT_EPS;
#else
      pos_t dist = POS_EPS;
#endif	  

      // compute time for the next subpartition collision, let's assume that displacemnt
      // is our speed vector and the total time to travel is 1
      //
      // pos(time) = pos + displacement * time, therefore
      // time = (pos(time) - vm.v.pos) / displacement
      // =>
      // time_to_subpart_edge = (subpart_edge - vm.v.pos) / displacement_speed
      Vec3 coll_times = diff * displacement_rcp;
      //note: one of coll_times may be negative due to float imprecisions

      // which of the times is the smallest? - i.e. which boundary we hit first (values may be negative due to flt32 imprecisions)
      if (coll_times.x < coll_times.y && coll_times.x <= coll_times.z) {
        // new position on the edge of the subpartition
        curr_pos += displacement * coll_times.x;
        // and also update the xyz subpartition index
        curr_subpart_indices.x += dir_urb_addend.x;
        if (!p.is_subpart_index_in_range(curr_subpart_indices.x)) {
          break;
        }
      }
      else if (coll_times.y <= coll_times.z) {
        curr_pos += displacement * coll_times.y;
        curr_subpart_indices.y += dir_urb_addend.y;
        if (!p.is_subpart_index_in_range(curr_subpart_indices.y)) {
          break;
        }
      }
      else  {
        curr_pos += displacement * coll_times.z;
        curr_subpart_indices.z += dir_urb_addend.z;
        if (!p.is_subpart_index_in_range(curr_subpart_indices.z)) {
          break;
        }
      }


      curr_subpart_index = p.get_subpart_index_from_3d_indices(curr_subpart_indices);
      if (collect_for_walls) {
        crossed_subparts_for_walls.push_back(curr_subpart_index);
      }
      if (collect_for_molecules) {
        crossed_subparts_for_molecules.insert(curr_subpart_index);
      }

      // also neighbors
      if (collect_for_molecules && p.config.use_expanded_list) {
        collect_neighboring_subparts(
            p, curr_pos, curr_subpart_indices, rxn_radius_for_neighbors, sp_edge_length,
            crossed_subparts_for_molecules
        );
      }

    } while (curr_subpart_index != dest_subpart_index);
  }

  // finally check also neighbors in destination
  if (collect_for_molecules && p.config.use_expanded_list) {
    collect_neighboring_subparts(
        p, dest_pos, dest_subpart_indices, rxn_radius_for_neighbors, sp_edge_length,
        crossed_subparts_for_molecules
    );
  }

  return dest_subpart_index;
}

} // namespace CollisionUtil
} // namespace MCell

#endif // SRC4_COLLISION_UTILS_SUBPARTS_INC_
