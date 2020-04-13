
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
//#if 0 // legacy implementation, kept in case it was needed for reference
// This function checks if any of the neighboring subpartitions are within radius
// from pos and inserts them into crossed_subparition_indices
static void collect_neighboring_subparts(
    const Partition& p,
    const Vec3& pos,
    const IVec3& subpart_indices,
    const float_t rx_radius,
    const float_t subpart_edge_len,
    subpart_indices_set_t& crossed_subpart_indices
) {
  Vec3 rel_pos = pos - p.get_origin_corner();

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
  const Vec3& displacement,
  const float_t rx_radius,
  const float_t sp_edge_length,
  const bool collect_for_walls,
  SubpartIndicesVector& crossed_subparts_for_walls, // crossed subparts considered for wall collision
  subpart_indices_set_t& crossed_subparts_for_molecules // crossed subparts considered for molecule collisions
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
  crossed_subparts_for_molecules.insert(current_subpart_index);

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

    IVec3 dir_urb_addend;
    dir_urb_addend.x = (dir_urb_direction.x == 0) ? -1 : 1;
    dir_urb_addend.y = (dir_urb_direction.y == 0) ? -1 : 1;
    dir_urb_addend.z = (dir_urb_direction.z == 0) ? -1 : 1;

    Vec3 curr_pos = vm.v.pos;
    IVec3 curr_subpart_indices = src_subpart_indices;

    subpart_index_t curr_subpart_index;

    Vec3 displacement_rcp = 1.0/displacement_nonzero; // POSSIBLE ZERO DIV

    do {
      // subpartition edges
      // = origin + subparition index * length + is_urb * length

      // NOTE: some of these computation can be moved out of the loop, not sure if it will improve performance
      Vec3 sp_len_as_vec3 = Vec3(sp_edge_length);
      Vec3 sp_edges =
          p.get_origin_corner()
          +  Vec3(curr_subpart_indices) * sp_len_as_vec3 // llf edge
          + Vec3(dir_urb_direction) * sp_len_as_vec3; // move if we go urb

      Vec3 diff = sp_edges - curr_pos;

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
        Vec3 coll_times = diff * displacement_rcp;
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

} // namespace CollisionUtil
} // namespace MCell
