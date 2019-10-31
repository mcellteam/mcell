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
  subpart_indices_set_t& crossed_subparts_for_molecules, // crossed subparts considered for molecule collisions
  subpart_index_t& last_subpart_index
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
  p.get_subpart_3d_indices_from_index(vm.v.subpart_index, src_subpart_indices);
  p.get_subpart_3d_indices(dest_pos, dest_subpart_indices);

  // first check what's around the starting point
  collect_neighboring_subparts(
      p, vm.v.pos, src_subpart_indices, rx_radius, sp_edge_length,
      crossed_subparts_for_molecules
  );

  // collect subpartitions on the way by always finding the point where a subpartition boundary is hit
  // we must do it eve when we are crossing just one subpartition because we might hit others while
  // moving along them
  if ( !glm::all( glm::equal(dest_subpart_indices, src_subpart_indices) ) ) {

    subpart_index_t dest_subpart_index = p.get_subpartition_index_from_3d_indices(dest_subpart_indices);
    last_subpart_index = dest_subpart_index;

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
        curr_subpart_indices.x += dir_urb_addend.x;
      }
      else if (abs(diff.y) < EPS) {
        curr_subpart_indices.y += dir_urb_addend.y;
      }
      else if (abs(diff.z) < EPS) {
        curr_subpart_indices.z += dir_urb_addend.z;
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
        assert(coll_times.x >= 0 && coll_times.y >= 0 && coll_times.z >= 0 && "Edges must be computed from direction");

        // which of the times is the smallest? - i.e. which boundary we hit first
        if (coll_times.x >= 0 && coll_times.x < coll_times.y && coll_times.x <= coll_times.z) {
          // new position on the edge of the subpartition
          curr_pos += displacement * coll_times.x;
          // and also update the xyz subpartition index
          curr_subpart_indices.x += dir_urb_addend.x;
        }
        else if (coll_times.y >= 0 && coll_times.y <= coll_times.z) {
          curr_pos += displacement * coll_times.y;
          curr_subpart_indices.y += dir_urb_addend.y;
        }
        else if (coll_times.z >= 0) {
          curr_pos += displacement * coll_times.z;
          curr_subpart_indices.z += dir_urb_addend.z;
        }
        else {
          break;
        }
      }

      curr_subpart_index = p.get_subpartition_index_from_3d_indices(curr_subpart_indices);
      if (collect_for_walls) {
        crossed_subparts_for_walls.push_back(curr_subpart_index);
      }
      crossed_subparts_for_molecules.insert(curr_subpart_index);

      // also neighbors
      collect_neighboring_subparts(
          p, curr_pos, curr_subpart_indices, rx_radius, sp_edge_length,
          crossed_subparts_for_molecules
      );

    } while (curr_subpart_index != dest_subpart_index);
  }
  else {
    // subpartition index did not change
    last_subpart_index = vm.v.subpart_index;
  }

  // finally check also neighbors in destination
  collect_neighboring_subparts(
      p, dest_pos, dest_subpart_indices, rx_radius, sp_edge_length,
      crossed_subparts_for_molecules
  );
}


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
    const WorldConstants& world_constants,
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
    const Reaction* rx = world_constants.get_reaction(vm, colliding_vm);
    assert(rx != nullptr);
    molecule_collisions.push_back(
        Collision(CollisionType::VOLMOL_VOLMOL, &p, vm.id, time, position, colliding_vm.id, rx)
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
  p.get_simulation_stats().inc_ray_polygon_tests();

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

    // NOTE: what is this doing?
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

static void collect_wall_collisions(
    Partition& p,
    const Molecule& vm, // molecule that we are diffusing, we are changing its pos  and possibly also subvolume
    const subpart_index_t subpart_index,
    const wall_index_t previous_reflected_wall,
    rng_state& rng,
    vec3_t& displacement, // can be changed
    collision_vector_t& wall_collisions
) {
  // check each wall in this subpartition
  const UintSet& wall_indices = p.get_subpart_wall_indices(subpart_index);

  for (UintSet::const_iterator it = wall_indices.begin(); it != wall_indices.end(); it++) {
    wall_index_t wall_index = *it;

    if (wall_index == previous_reflected_wall)
      continue;

    const Wall& w = p.get_wall(wall_index);
    float_t collision_time;
    vec3_t collision_pos;

#ifdef DEBUG_COLLISIONS_WALL_EXTRA
    SimulationStats* world = &p.get_simulation_stats();
    // just faking the name for the dump condition macro
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
      p.get_simulation_stats().inc_ray_polygon_colls();
      wall_collisions.push_back(
          Collision(collision_type, &p, vm.id, collision_time, collision_pos, wall_index)
      );
    }
  }
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
