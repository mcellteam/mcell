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

#ifndef SRC4_COLLISION_UTILS_INC_
#define SRC4_COLLISION_UTILS_INC_

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

#ifdef COLLECT_SUBPARTS_LEGACY
#include "collision_utils_legacy.inc"
#endif


namespace MCell {

namespace CollisionUtil {

// ---------------------------------- subpartitions ----------------------------------


// move to collision utils or to partition
static Vec3 get_displacement_up_to_partition_boundary(
    const Partition& p,
    const Vec3& pos,
    const Vec3& displacement
    ) {
#ifndef NDEBUG
  Vec3 new_pos = pos + displacement;
  assert(p.in_this_partition(pos));
  assert(!p.in_this_partition(new_pos));
#endif

  Vec3 displacement_nonzero = displacement;
  guard_zero_div(displacement_nonzero);
  IVec3 dir_urb_direction = IVec3(glm::greaterThan(displacement_nonzero, Vec3(0)));

  // position of edges in our direction
  Vec3 partition_edges =
      p.get_origin_corner()
      + Vec3(dir_urb_direction) * p.config.partition_edge_length;

  Vec3 diff = partition_edges - pos;

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
    Vec3 coll_times = diff / displacement_nonzero;
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
  Vec3 new_displacement = displacement * (hit_time - EPS);
  assert(p.in_this_partition(pos + new_displacement));

  return new_displacement;
}


static inline void get_corner_points_for_subpart_colection(
    const Vec3& pos,
    const Vec3& move,
    const float_t rx_radius,
    small_vector<Vec3>& pts
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

       vz = (-mx-my)/mz

       assuming that mz is not 0, then:
       v = (1, 1, -(mx+my / mz))

       precision of the initial points is not so important
   */
  uint largest_dim = get_largest_abs_dim_index(move);
  Vec3 v;
  switch(largest_dim) {
    case 0:
      v = Vec3( -(move.y + move.z)/move.x, 1, 1);
      break;
    case 1:
      v = Vec3( 1, -(move.x + move.z)/move.y, 1);
      break;
    case 2:
      v = Vec3( 1, 1, -(move.x + move.y)/move.z);
      break;
    default:
      assert(false);
  }

  // 2) get a vector perpendicular to v and move
  Vec3 vp = cross(v, move);

  // 3) and compute the vectors v1-v4 of length rx_radius * SQRT2
  float_t radius = rx_radius * SQRT2 * RX_RADIUS_MULTIPLIER;
  float_t ratio = 1/len3(v) * radius;
  float_t ratiop = 1/len3(vp) * radius;

  Vec3 v0 = v * Vec3(ratio);
  Vec3 v1 = vp * Vec3(ratiop);
  Vec3 v2 = -v1;
  Vec3 v3 = -v0;

  assert(
      cmp_eq(len3(v0), radius, SQRT_EPS) && cmp_eq(len3(v1), radius, SQRT_EPS) &&
      cmp_eq(len3(v2), radius, SQRT_EPS) && cmp_eq(len3(v3), radius, SQRT_EPS));
  assert(
      cmp_eq(dot(v0, move), 0, SQRT_EPS) && cmp_eq(dot(v1, move), 0, SQRT_EPS) &&
      cmp_eq(dot(v2, move), 0, SQRT_EPS) && cmp_eq(dot(v3, move), 0, SQRT_EPS));

  pts[0] = pos + v0;
  pts[1] = pos + v1;
  pts[2] = pos + v2;
  pts[3] = pos + v3;

#ifdef DEBUG_SUBPARTITIONS
  std::cout << "get_corner_points_for_subpart_colection: pos: " << pos << ", \n"
      << "0:" << pts[0] << "\n1:" << pts[1] << "\n2:" << pts[2] << "\n3:" << pts[3] << "\n";
#endif
}


/*
 Cast a ray through a volume by specifying the start and end positions

 base implementation from https://bitbucket.org/volumesoffun/polyvox/src/9a71004b1e72d6cf92c41da8995e21b652e6b836/include/PolyVox/Raycast.inl?at=develop&fileviewer=file-view-default

 The MIT License (MIT)

 Copyright (c) 2015 David Williams and Matthew Williams

 This function is based on Christer Ericson's code and description of the 'Uniform Grid Intersection Test' in
 'Real Time Collision Detection'. The following information from the errata on the book website is also relevent:

  pages 326-327. In the function VisitCellsOverlapped() the two lines calculating tx and ty are incorrect.
  The less-than sign in each line should be a greater-than sign. That is, the two lines should read:

  float tx = ((x1 > x2) ? (x1 - minx) : (maxx - x1)) / Abs(x2 - x1);
  float ty = ((y1 > y2) ? (y1 - miny) : (maxy - y1)) / Abs(y2 - y1);

  Thanks to Jetro Lauha of Fathammer in Helsinki, Finland for reporting this error.
*/

static void raycast_with_endpoints(
    const Partition& p, const Vec3& pt1, const Vec3& pt2,
    const uint pt_index /*only for debug*/, const bool collect_for_walls, /* wall subparts are collected for the internal ray */
    const IVec3& dir, const Vec3& abs_d_rcp, const Vec3& deltat,
    SubpartIndicesVector& crossed_subparts_for_walls,
    subpart_indices_set_t& crossed_subparts_for_molecules
)
{
  // FIXME: used lower case, these are not really contants...
  const float_t cell_side = p.config.subpartition_edge_length;
  const float_t cell_size_rcp = p.config.subpartition_edge_length_rcp;

  // computing subpart index using get_subpart_index_from_3d_indices every time is rather expensive
  const int subpart_x_addend = dir.x;
  const int subpart_y_addend = dir.y * p.config.num_subpartitions_per_partition;
  const int subpart_z_addend = dir.z * p.config.num_subpartitions_per_partition_squared;

  IVec3 indices;
  p.get_subpart_3d_indices(pt1, indices);
  subpart_index_t subpart_index = p.get_subpart_index_from_3d_indices(indices);

  IVec3 end_indices;
  p.get_subpart_3d_indices(pt2, end_indices);

  // TODO: optimize when start and end are the same? check with benchmark

  Vec3 min = floor3(pt1 * Vec3(cell_size_rcp)) * Vec3(cell_side);
  Vec3 max = min + Vec3(cell_side);

  Vec3 t;
  t.x = ((dir.x == -1) ? (pt1.x - min.x) : (max.x - pt1.x)) * abs_d_rcp.x;
  t.y = ((dir.y == -1) ? (pt1.y - min.y) : (max.y - pt1.y)) * abs_d_rcp.y;
  t.z = ((dir.z == -1) ? (pt1.z - min.z) : (max.z - pt1.z)) * abs_d_rcp.z;

#ifdef DEBUG_SUBPARTITIONS
  std::cout << "Point " << pt_index << ", start subpart indices " << indices << "\n";
#endif

  do {
    if (!collect_for_walls) {
      crossed_subparts_for_molecules.insert(subpart_index);
    }
    if (collect_for_walls) {
      crossed_subparts_for_walls.push_back(subpart_index);
    }

    if (t.x <= t.y && t.x <= t.z)
    {
      if (indices.x == end_indices.x) {
        break;
      }
      t.x += deltat.x;
      indices.x += dir.x;
      subpart_index += subpart_x_addend;
    }
    else if (t.y <= t.z)
    {
      if (indices.y == end_indices.y) {
        break;
      }
      t.y += deltat.y;
      indices.y += dir.y;
      subpart_index += subpart_y_addend;
    }
    else
    {
      if (indices.z == end_indices.z) {
        break;
      }
      t.z += deltat.z;
      indices.z += dir.z;
      subpart_index += subpart_z_addend;
    }

    #ifdef DEBUG_SUBPARTITIONS
      std::cout << "Point " << pt_index << ", new subpart indices: " << indices << "\n";
    #endif
  }
  while (true);

#ifdef DEBUG_SUBPARTITIONS
  std::cout << "Point " << pt_index << ", last subpart indices: " << indices << "\n";
#endif

  // note: indices != end_indices at this point, this is ok (at least all tests pass)
}


static inline void compute_inputs_for_raycast_with_endpoints(
    const float_t subpartition_edge_length, const Vec3& displacement,
    Vec3& dir, Vec3& abs_d_rcp, Vec3& deltat
) {
  const int di = ((displacement.x > 0) ? 1 : ((displacement.x < 0) ? -1 : 0)); // displacement direction
  const int dj = ((displacement.y > 0) ? 1 : ((displacement.y < 0) ? -1 : 0));
  const int dk = ((displacement.z > 0) ? 1 : ((displacement.z < 0) ? -1 : 0));
  dir = Vec3(di, dj, dk);

  // corner points
  Vec3 abs_d = abs3(displacement);
  abs_d.x = (abs_d.x == 0) ? EPS : abs_d.x;
  abs_d.y = (abs_d.y == 0) ? EPS : abs_d.y;
  abs_d.z = (abs_d.z == 0) ? EPS : abs_d.z;

  abs_d_rcp = float_t(1.0)/abs_d;
  deltat = Vec3(subpartition_edge_length) * abs_d_rcp;
}


#ifndef COLLECT_SUBPARTS_LEGACY
static inline void collect_crossed_subparts(
  const Partition& p,
  const Molecule& vm, // molecule that we are diffusing
  const Vec3& displacement,
  const float_t rx_radius,
  const float_t sp_edge_length,
  SubpartIndicesVector& crossed_subparts_for_walls, // crossed subparts considered for wall collision
  subpart_indices_set_t& crossed_subparts_for_molecules // crossed subparts considered for molecule collisions
) {
  assert(p.in_this_partition(vm.v.pos));
  assert(crossed_subparts_for_walls.empty());
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
  small_vector<Vec3> start_positions;
  start_positions.resize(NUM_CORNER_POINTS + 1);

  uint num_points = start_positions.size();

  // set pts 0 - 3
  get_corner_points_for_subpart_colection(vm.v.pos, displacement, rx_radius, start_positions);

  // we need to move the points a bit backwards and also forwards
  float_t displacement_length = len3(displacement);
  Vec3 displacement_unit = displacement/Vec3(displacement_length);
  Vec3 displacement_of_radius_length = displacement_unit * Vec3(rx_radius * SQRT2 * RX_RADIUS_MULTIPLIER);

  // move molecule collision detection points a bit back
  small_vector<subpart_index_t> start_subpart_indices;
  for (uint i = 0; i < NUM_CORNER_POINTS; i++) {
    Vec3 start_pos = start_positions[i] - displacement_of_radius_length;

    if (!p.in_this_partition(start_pos)) {
      // we got out of the partition, move the corner point just a little so that we still fit
      // because we are sending negated displacement, we also get negated result
      Vec3 new_displacement_neg = get_displacement_up_to_partition_boundary(p, start_positions[i], -displacement_of_radius_length);
      start_pos = start_positions[i] + new_displacement_neg;
      assert(p.in_this_partition(start_pos));
    }

    start_positions[i] = start_pos;
  }
  // move wall detection point only a tiny bit back to deal with precision issues
  start_positions[MOL_POS_POINT_INDEX] = vm.v.pos - displacement_unit * Vec3(SQRT_EPS);

#ifdef DEBUG_SUBPARTITIONS
  std::cout << "Corrected corner points for subpart colection:\n"
      << "0:" << start_positions[0] << "\n1:" << start_positions[1] << "\n2:" << start_positions[2] << "\n3:" << start_positions[3] << "\n";
#endif

  // we already subtracted displacement_of_radius_length from the current positions
  Vec3 extended_displacement = displacement + Vec3(2)*displacement_of_radius_length;

  small_vector<Vec3> dest_positions;
  dest_positions.resize(num_points);

  for (uint i = 0; i < num_points; i++) {
    Vec3 dest_pos;

    if (i == MOL_POS_POINT_INDEX) {
      // similarly as we moved the start, let's move the move wall detection end point a bit further
      dest_pos = start_positions[i] + displacement + displacement_unit * Vec3(SQRT_EPS);;
    }
    else {
      dest_pos = start_positions[i] + extended_displacement;
    }

    // truncate if needed, does not happen often
    if (!p.in_this_partition(dest_pos)) {
      dest_pos = start_positions[i] + get_displacement_up_to_partition_boundary(p, start_positions[i], extended_displacement);
    }

    dest_positions[i] = dest_pos;

    #ifdef DEBUG_SUBPARTITIONS
    DUMP_CONDITION4P(
      std::cout << "Point " << i << ", start pos: " << start_positions[i] << "\n";
      std::cout << "Point " << i << ", dest  pos: " << dest_positions[i] << "\n";
    );
    #endif
  }

  // for each dimension, process each trace and each plane
  // collect subpartitions on the way by always finding the point where a subpartition boundary is hit
  // we must do it even when we are crossing just one subpartition because we might hit others while
  // moving along them
  Vec3 dir;
  Vec3 abs_d_rcp;
  Vec3 deltat;
  compute_inputs_for_raycast_with_endpoints(
      p.config.subpartition_edge_length, extended_displacement,
      dir, abs_d_rcp, deltat
  );

  for (uint i = 0; i < start_positions.size(); i++){
    raycast_with_endpoints(
        p, start_positions[i], dest_positions[i],
        i, i == MOL_POS_POINT_INDEX,
        dir, abs_d_rcp, deltat,
        crossed_subparts_for_walls, crossed_subparts_for_molecules
    );
  }

  // first for walls might not be present either, should be the first item
  assert(!crossed_subparts_for_walls.empty());
  if (crossed_subparts_for_walls.front() != vm.v.subpart_index) {
    crossed_subparts_for_walls.push_back(vm.v.subpart_index);
  }

  // some subparts can be missed when checking just the corner points
  for (subpart_index_t index: crossed_subparts_for_walls) {
    crossed_subparts_for_molecules.insert(index);
  }
}

#endif // COLLECT_SUBPARTS_LEGACY

// ---------------------------------- molecule collisions ----------------------------------

// check whether diffused_vm molecule collision that moves by displacement can collide
// with colliding_vm; returns true if there can be a collision and returns relative collision
// time and relative position
static bool collide_mol(
    const Molecule& diffused_vm,
    const Vec3& displacement,
    const Molecule& colliding_vm,
    const float_t rx_radius_3d,
    float_t& rel_collision_time,
    Vec3& rel_collision_pos
) {

  const Vec3& pos = colliding_vm.v.pos; /* Position of target molecule */
  Vec3 dir = pos - diffused_vm.v.pos;  /* From starting point of moving molecule to target */

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
    const Vec3& remaining_displacement,
    const float_t radius,
    collision_vector_t& molecule_collisions
) {

  Molecule& colliding_vm = p.get_m(colliding_vm_id);

  // we would like to compute everything that's needed just once
  float_t time;
  Vec3 position;
  // collide_mol must be inlined because many things are computed all over there
  if (collide_mol(vm, remaining_displacement, colliding_vm, radius, time, position)) {

    BNG::RxnClass* rxn_class =
        p.get_all_rxns().get_bimol_rxn_class(vm.species_id, colliding_vm.species_id);

    assert(rxn_class != nullptr && "If a collision of two molecules was detected, there must be a reaction for them.");

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
    const Vec3& p,
    const float_t k, const Vec3& A, const Vec3& B, const Vec3& n, rng_state& rng,
    Vec3& v /*inout*/
) {
  Vec3 e, f;
  float_t le_1, tiny;

  e = B - A;
  float_t elen2 = glm::dot((glm_vec3_t)e, (glm_vec3_t)e);
  le_1 = 1.0 / sqrt(elen2);

  e = e * Vec3(le_1);

  f.x = n.y * e.z - n.z * e.y;
  f.y = n.z * e.x - n.x * e.z;
  f.z = n.x * e.y - n.y * e.x;

  tiny = EPS * (MCell::abs_max_2vec(p, v) + 1.0) /
         (k * max3(glm::abs((glm_vec3_t)f)));
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
    const Vec3& pos, const Wall& face,
    rng_state &rng,
    const bool update_move,
    Vec3& move,
    float_t& collision_time, Vec3& collision_pos
) {
  p.stats.inc_ray_polygon_tests();

  float_t dp, dv, dd;
  float_t d_eps;

  const Vec3& normal = face.normal;

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
      a = (MCell::abs_max_2vec(pos, move) + 1.0) * EPS;
      if ((rng_uint(&rng) & 1) == 0)
        a = -a;
      if (dd == 0.0) {
        move = move - Vec3(a) * normal;
      }
      else {
        move = move * Vec3(1.0 - a);
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

  Vec3 face_vert[VERTICES_IN_TRIANGLE];
  face_vert[0] = p.get_geometry_vertex(face.vertex_indices[0]);
  face_vert[1] = p.get_geometry_vertex(face.vertex_indices[1]);
  face_vert[2] = p.get_geometry_vertex(face.vertex_indices[2]);

  Vec3 local = collision_pos - face_vert[0];

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
          EPS))) {
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
    const subpart_index_t subpart_index, // the wall that we hit must be inside of this partition
    const wall_index_t last_hit_wall_index,
    rng_state& rng,
    // displacement can be changed in case we needed to 'REDO' the collision, also,
    // if there was a hit, is changed to the closest displacement
    Vec3& displacement,
    Vec3& displacement_up_to_wall_collision, // overwritten only when there is a wall collision
    collision_vector_t& wall_collisions
) {

  // remember which was the closest hit to update displacement
  float_t closest_hit_time = TIME_FOREVER;
  Vec3 closest_hit_pos(0);

  // check each wall in this subpartition
  const uint_set<wall_index_t>& wall_indices = p.get_subpart_wall_indices(subpart_index);

  for (auto it = wall_indices.begin(); it != wall_indices.end(); it++) {
    wall_index_t wall_index = *it;

    if (wall_index == last_hit_wall_index) {
      continue;
    }

    const Wall& w = p.get_wall(wall_index);
    float_t collision_time;
    Vec3 collision_pos;

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

      // the hit must be inside this subpartition otherwise an optimization in ray_trace()
      // that tells to end searching for hits in other subparts once a hit was found won't work
      // (we might have walls that are in the current subpart but a hit will actually occur later
      //  than a hit in another subpart)
      if (p.get_subpart_index(collision_pos) != subpart_index) {
        continue;
      }

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
    const Vec3& pos, const Wall& face,
    const Vec3& move
) {
  rng_state unused_rng_state;
  float_t ignored_collision_time;
  Vec3 ignored_collision_pos;

  Vec3 tmp_move = move;
  CollisionType res;
#ifndef NDEBUG
  unused_rng_state.randcnt = 0;
  uint orig_randcnt = unused_rng_state.randcnt;
#endif
  res = CollisionUtil::collide_wall(
      p, pos, face, unused_rng_state, false, tmp_move,
      ignored_collision_time, ignored_collision_pos
  );
#ifndef NDEBUG
  assert(orig_randcnt == unused_rng_state.randcnt && "collide_wall_test should not trigger usage of rng");
#endif

#ifdef DEBUG_COLLISIONS_WALL_EXTRA
  Vec3 face_vert[VERTICES_IN_TRIANGLE];
  face_vert[0] = p.get_geometry_vertex(face.vertex_indices[0]);
  face_vert[1] = p.get_geometry_vertex(face.vertex_indices[1]);
  face_vert[2] = p.get_geometry_vertex(face.vertex_indices[2]);

  cout <<
      "Testing wall collision: pos: " << pos << ", move: " << move <<
      ", v1:" << face_vert[0] << ", v2:" << face_vert[0] << ", v3:" << face_vert[0] << "\n  result: " ;

  switch (res) {
    case CollisionType::WALL_MISS: cout << "MISS"; break;
    case CollisionType::WALL_FRONT: cout << "FRONT"; break;
    case CollisionType::WALL_BACK: cout << "BACK"; break;
    case CollisionType::WALL_REDO: cout << "REDO"; break;
    default: cout << "UNEXPECTED"; break;
  }
  cout << "\n";
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
    const Vec3& e, const Vec3& f, const Vec3& vfe /* == f - e*/,
    const Vec3& o, const Vec3& p
) {
  // rename arguments to get the names used in comment
  const Vec3& x1 = e;
  const Vec3& x2 = f;
  const Vec3& x3 = o;
  const Vec3& x4 = p;

  const Vec3& a = vfe;
  Vec3 b = x4 - x3;
  Vec3 c = x3 - x1;

  Vec3 a_x_b = cross(a, b);
  Vec3 c_x_b = cross(c, b);

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
  Vec3 d = x1 - x3;

  // t = ((c'xb')·(a'xb'))/(|a'xb'|^2)
  // t = ((dxa)·(bxa))/(|bxa|^2)
  Vec3 d_x_a = cross(d, a);
  Vec3 b_x_a = cross(b, a);

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
    const Vec3& e, const Vec3& f,
    const Vec3& k, const Vec3& l,
    const Vec3& m, const Vec3& n,
    const float_t t
) {
  return
      dot(
          cross((e - f), (-f + k + (-k + l) * Vec3(t)) ),
          (-e + m + (-m + n) * Vec3(t))
      );
}

static float_t compute_df(
    const Vec3& e, const Vec3& f,
    const Vec3& k, const Vec3& l,
    const Vec3& m, const Vec3& n,
    const float_t t
) {
  return
      dot(
          cross(
              (e - f),
              (-k + l)
          ),
          (-e + m + (- m + n) * Vec3(t))
      )
      +
      dot(
          cross(
              (e - f),
              (-f + k + (- k + l) * Vec3(t))
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
    const Vec3& e, const Vec3& f, const Vec3& move /* f = e + move */,
    const Vec3& k, const Vec3& l,
    const Vec3& m, const Vec3& n,
    Vec3& o, Vec3& p) {

  assert(cmp_eq(f, e + move));

  // first determine t, it must be in the range 0..1
  // f is continuous and (should?) have just one solution

  // starting from 0
  float_t t = 0;
  float_t t_previous = FLT_GIGANTIC;

  // use Newton's Method to find solution to 't' (https://en.wikipedia.org/wiki/Newton's_method)
  // using high precision, we should converge quickly
  bool dft_is_zero = false;
  bool t_out_of_range = false;

  float_t ft;
  while (!cmp_eq(t, t_previous, SQRT_EPS)) {
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
  o = k + (l - k) * Vec3(t);
  p = m + (n - m) * Vec3(t);
  return true;
}

static bool collide_moving_line_and_static_line_test(
    const Partition& part,
    const Vec3& e, const Vec3& move /* f = e + move */,
    const Vec3& k, const Vec3& l,
    const Vec3& m, const Vec3& n,
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
  Vec3 o, p;
  Vec3 f = e + move;
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
static bool is_point_inside_region(const Partition& p, const Vec3& pos, const Region& reg) {

  // cast ray along the whole partition
  Vec3 move(EPS, EPS, p.config.partition_edge_length);

  // and count intersects
  int num_hits = 0;
  for (const auto& wall_edges_pair: reg.walls_and_edges) {

    // check collision with each wall
    // TODO LATER: optimize with subparts
    // FIXME: the collide_Wall_test is not reliable in all cases
    const Wall& w = p.get_wall(wall_edges_pair.first);
    bool collides = CollisionUtil::collide_wall_test(p, pos, w, move);
    if (collides) {
      num_hits++;
    }
  }

  // odd number of hits means that we are inside
  return num_hits % 2 == 1;
}

// ---------------------------------- counting ----------------------------------

// does not set res.index, the index needs to be obtained from partition
static counted_volume_index_t compute_counted_volume_for_pos(Partition& p, const Vec3& pos) {

  CountedVolume cv;
#if 1
  // send a ray trace and collect all walls we crossed,
  // then if we hit wall of an object odd number of times, we are inside of an object

  // assuming that the whole object fits into our partition
  // can be optimized by choosing some close object
  Vec3 dst = p.get_origin_corner();

  // move a bit along by a magic number the side of the partition so that the ray trace hits correctly
  // this is a temporary measure to minimize REDOs with wall collision
  dst.x += p.config.subpartition_edge_length / 10;
  dst.y += p.config.subpartition_edge_length / 20;

  Vec3 displacement = dst - pos - Vec3(EPS); // from here up to a corner of the partition

  // collect which subpartitions we crossed
  Vec3 dir;
  Vec3 abs_d_rcp;
  Vec3 deltat;
  compute_inputs_for_raycast_with_endpoints(
      p.config.subpartition_edge_length, displacement,
      dir, abs_d_rcp, deltat
  );

  SubpartIndicesVector crossed_subparts_for_walls;
  subpart_indices_set_t crossed_subparts_for_molecules_ignored;
  raycast_with_endpoints(
      p, pos, dst,
      0, true,
      dir, abs_d_rcp, deltat,
      crossed_subparts_for_walls, crossed_subparts_for_molecules_ignored
  );

  // now check which walls we hit
  map<geometry_object_index_t, uint> crossed_wall_per_object;
  uint_set<wall_index_t> already_checked_walls;
  for (subpart_index_t subpart_w_walls_index: crossed_subparts_for_walls) {

    // check each wall in this subpartition
    const uint_set<wall_index_t>& wall_indices = p.get_subpart_wall_indices(subpart_w_walls_index);

    for (auto it = wall_indices.begin(); it != wall_indices.end(); it++) {
      wall_index_t wall_index = *it;

      if (already_checked_walls.count(wall_index) != 0) {
        // we already checked this wall - the same wall can be in multiple subparts
        continue;
      }
      already_checked_walls.insert(wall_index);

      const Wall& w = p.get_wall(wall_index);

      // do we care about this object?
      if (!p.get_geometry_object(w.object_index).is_counted_volume) {
        continue;
      }

      float_t collision_time_ignored;
      Vec3 collision_pos_ignored;

      rng_state rng; // using a separate random generator for jump_away lines
      rng_init(&rng, 0);

      CollisionType collision_type;
      uint num_tries = 0;
      do {
        collision_type = collide_wall(
              p, pos, w, rng, true, displacement,
              collision_time_ignored, collision_pos_ignored);

        num_tries++;
        release_assert(num_tries < 100 && "Cannot determine molecule or waypoint position.");
        // FIXME: we really do not want REDO, it is highly unreliable
        assert(collision_type != CollisionType::WALL_REDO);
      } while (collision_type == CollisionType::WALL_REDO);

      if (collision_type == CollisionType::WALL_FRONT || collision_type == CollisionType::WALL_BACK) {
        auto it = crossed_wall_per_object.find(w.object_index);
        if (it == crossed_wall_per_object.end()) {
          crossed_wall_per_object[w.object_index] = 1;
        }
        else {
          it->second++;
        }
      }
    }
  }

  // finally construct the set of object indices that encompass the position
  cv.contained_in_objects.clear();
  for (auto it: crossed_wall_per_object) {
    assert(it.second != 0);
    // hit object odd number of times - pos is inside
    if (it.second % 2 == 1) {
      // inside of this object
      cv.contained_in_objects.insert(it.first);
    }
  }
#else
  // the implementation above does not work correctly in all cases,
  // the WALL_REDO handling is wrong and not sure how it can be fixed
  // using VTK to make sure that we get correct result
  // maybe we can use the implementation above to optimize the counting,
  // but for now let's keep it simple

  cv.contained_in_objects.clear();
  for (GeometryObject& obj: p.get_geometry_objects()) {
    if (obj.is_counted_volume) {
      if (CountedVolumesUtil::is_point_inside_counted_volume(obj, pos)) {
        cv.contained_in_objects.insert(obj.index);
      }
    }
  }
#endif

  return p.find_or_add_counted_volume(cv);
}


// TODO: optimization - search in a direction where a molecule is moving,
// we must probably just crossed a wall
static counted_volume_index_t compute_counted_volume_using_waypoints(
    Partition& p,
    const Vec3& pos) {

  bool wall_found = false;
  const Waypoint* waypoint = nullptr;

  if (p.config.has_intersecting_counted_objects) {

    // subpartition origin point
    IVec3 index3d;
    subpart_index_t subpart_index = p.get_subpart_index(pos);
    p.get_subpart_3d_indices_from_index(subpart_index, index3d);

    // waypoint is always in the center of a subpartition
    waypoint = &p.get_waypoint(index3d);

    // there must not be a counted wall between our point and the waypoint
    Vec3 displacement = waypoint->pos - pos;
    const uint_set<wall_index_t>& wall_indices = p.get_subpart_wall_indices(subpart_index);
    for (auto it = wall_indices.begin(); it != wall_indices.end(); it++) {
      const Wall& w = p.get_wall(*it);

      // do we care about this object?
      if (!p.get_geometry_object(w.object_index).is_counted_volume) {
        continue;
      }

      float_t collision_time_ignored;
      Vec3 collision_pos_ignored;

      rng_state rng; // using a separate random generator for jump_away lines
      rng_init(&rng, 0);

      CollisionType collision_type;
      collision_type = collide_wall(
            p, pos, w, rng, true, displacement,
            collision_time_ignored, collision_pos_ignored
      );

      if (collision_type == CollisionType::WALL_FRONT ||
          collision_type == CollisionType::WALL_BACK ||
          collision_type == CollisionType::WALL_REDO) {
        wall_found = true;
        break;
      }
    }
  } // check with

  if (p.config.has_intersecting_counted_objects && !wall_found) {
    assert(waypoint != nullptr);
    // nothing is obstructing,
    // we can simply copy counted volume index from the waypoint
    p.stats.inc_num_waypoints_used();
    return waypoint->counted_volume_index;
  }
  else {
    // either we have no waypoints or we need to recompute the position
    // because a wall is obstructing
    p.stats.inc_recomputations_of_counted_volume();
    return CollisionUtil::compute_counted_volume_for_pos(p, pos);
  }
}


static void update_counted_volume_id_when_crossing_wall(
    Partition& p,
    const Wall& w,
    const Collision& collision,
    const Vec3& displacement, // TODO: not sure if really needed
    Molecule& vm
) {
  assert(vm.is_vol());

  const GeometryObject& obj = p.get_geometry_object(w.object_id);
  if (obj.is_counted_volume) {
    
    assert(obj.counted_volume_index_inside != COUNTED_VOLUME_INDEX_INVALID);
    assert(obj.counted_volume_index_outside != COUNTED_VOLUME_INDEX_INVALID);

    // which direction?
    if (collision.type == CollisionType::WALL_BACK) {
      // going outside
      assert(obj.counted_volume_index_inside == COUNTED_VOLUME_INDEX_INTERSECTS ||
          obj.counted_volume_index_inside == vm.v.counted_volume_index);

      if (obj.counted_volume_index_outside != COUNTED_VOLUME_INDEX_INTERSECTS) {
        // no intersect - it is clear what is outside
        vm.v.counted_volume_index = obj.counted_volume_index_outside;
      }
      else {

        float_t bump = EPS; // hit from back - we go in the direction of the normal
        Vec3 displacement = Vec3(2 * bump) * w.normal;

        // intersect - need to check waypoints or simply recompute the counted volumes
        vm.v.counted_volume_index =
            compute_counted_volume_using_waypoints(p, vm.v.pos + displacement);
      }
    }
    else if (collision.type == CollisionType::WALL_FRONT) {
      // going inside
      assert(obj.counted_volume_index_outside == COUNTED_VOLUME_INDEX_INTERSECTS ||
          obj.counted_volume_index_outside == vm.v.counted_volume_index);

      if (obj.counted_volume_index_outside != COUNTED_VOLUME_INDEX_INTERSECTS) {
        // no intersect - it is clear what is inside
        vm.v.counted_volume_index = obj.counted_volume_index_inside;
      }
      else {
        float_t bump = -EPS; // hit from front - we go against the direction of the normal
        Vec3 displacement = Vec3(2 * bump) * w.normal;

        vm.v.counted_volume_index =
            compute_counted_volume_using_waypoints(p, vm.v.pos + displacement);
      }
    }
    else {
      assert(false);
    }
  }
}

// ---------------------------------- reflections and other wall interactions ----------------------------------

static void cross_transparent_wall(
    Partition& p,
    const Collision& collision,
    const Vec3& full_displacement, // whole displacement computed for this molecule // TODO: not used
    Molecule& vm, // moves vm to the reflection point
    Vec3& remaining_displacement,
    float_t& t_steps,
    float_t& elapsed_molecule_time,
    wall_index_t& last_hit_wall_index
) {
#ifdef DEBUG_COUNTED_VOLUMES
  vm.dump(p, "- Before crossing: ");
#endif

  const Wall& w = p.get_wall(collision.colliding_wall_index);

  // Update molecule location to the point of wall crossing
  vm.v.pos = collision.pos;
  vm.v.subpart_index = p.get_subpart_index(vm.v.pos);

  update_counted_volume_id_when_crossing_wall(p, w, collision, remaining_displacement, vm);

#ifdef DEBUG_TRANSPARENT_SURFACES
  std::cout << "Crossed a transparent wall, side: " << w.side << "\n";
#endif

  // ignore the collision time, it is a bit earlier and does not fit for multiple collisions in the
  // same time step
  float_t t_smash = collision.time;

  remaining_displacement = remaining_displacement * Vec3(1.0 - t_smash);
  elapsed_molecule_time += t_steps * t_smash;

  t_steps *= (1.0 - t_smash);
  if (t_steps < EPS) {
    t_steps = EPS;
  }

  last_hit_wall_index = w.index;

#ifdef DEBUG_COUNTED_VOLUMES
  vm.dump(p, "- After crossing: ");
#endif
}

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
    Vec3& displacement,
    float_t& remaining_time_step, // same as t_steps
    wall_index_t& last_hit_wall_index
) {

  const Wall& w = p.get_wall(collision.colliding_wall_index);
  wall_index_t reflect_w = collision.colliding_wall_index;
  float_t t_reflect = collision.time;

  /* Update molecule location to the point of reflection (originally in register_hits) */
  vm.v.pos = collision.pos;
  vm.v.subpart_index = p.get_subpart_index(vm.v.pos);

  /* Reduce our remaining available time. */
  remaining_time_step *= (1.0 - t_reflect);

  last_hit_wall_index = reflect_w;

  float_t reflect_factor = -2.0 * glm::dot((glm_vec3_t)displacement, (glm_vec3_t)w.normal);

  // Set displacement for remainder of step length
  // No PBCs or non-traditional PBCs
  displacement = (displacement + Vec3(reflect_factor) * w.normal) * Vec3(1.0 - t_reflect);

  return 0;
}

} // namespace CollisionUtil

} // namespace MCell

#endif // SRC4_COLLISION_UTILS_INC_
