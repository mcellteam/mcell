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

#include "geometry_utils.inc"
#include "collision_utils.inc"
#include "diffusion_utils.inc"
#include "dyn_vertex_structs.h"

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


// returns the closest distance
// wall must belong to the same object and region
float_t Partition::find_closest_wall(
    const vec3_t& pos, const wall_index_t wall_that_moved_molecule,
    wall_index_t& best_wall_index,
    vec2_t& best_wall_pos2d
) {

  const Wall& moved_wall = get_wall(wall_that_moved_molecule);

  best_wall_index = WALL_INDEX_INVALID;
  float_t best_d2 = GIGANTIC4;
  best_wall_pos2d = vec2_t(0);

  // TODO: use wall_that_moved_molecule to optimize the search,
  // however for now we are using the same approach as in mcell3 because
  // we need to match the behavior
  const GeometryObject& obj = get_geometry_object(moved_wall.object_index);

  for (const wall_index_t wall_index: obj.wall_indices) {
    const Wall& w = get_wall(wall_index);

    // TODO: regions - we are checking only object id here
    if (!moved_wall.is_same_region(w)) {
      continue;
    }

    vec2_t wall_pos2d;
    float_t d2 = GeometryUtil::closest_interior_point(*this, pos, w, wall_pos2d);

    if (d2 <= best_d2) { // the <= is to emulate behavior of mcell3 that goes through the walls in opposite order
      best_d2 = d2;
      best_wall_index = w.index;
      best_wall_pos2d = wall_pos2d;
    }
  }

  return best_d2;
}

// based on place_mol_relative_to_mesh
void Partition::move_volume_molecule_to_closest_wall_point(const VolumeMoleculeMoveInfo& molecule_move_info) {

  Molecule& vm = get_m(molecule_move_info.molecule_id);
  assert(vm.is_vol());

  // 1) check all walls to get a reference hopefully
  //  then try to limit only to wall's neighbors - we would really like to avoid any regions (maybe...?)
  wall_index_t best_wall_index;
  vec2_t best_wall_pos2d;
  find_closest_wall(vm.v.pos, molecule_move_info.wall_index, best_wall_index, best_wall_pos2d);

  if (best_wall_index == WALL_INDEX_INVALID) {
    mcell_error("Could not find a wall close to volume molecule with id %d.\n", vm.id);
  }
  Wall& wall = get_wall(best_wall_index);

  vec3_t new_pos3d = GeometryUtil::uv2xyz(best_wall_pos2d, wall, get_wall_vertex(wall, 0));

#ifdef DEBUG_DYNAMIC_GEOMETRY
  vm.dump(*this, "", "Moving vm towards new wall: ", simulation_stats.get_current_iteration(), 0);
  wall.dump(*this, "", true);
#endif

  // displacement
  float_t bump = (molecule_move_info.place_above) ? EPS : -EPS;
  vec3_t displacement(2 * bump * wall.normal.x, 2 * bump * wall.normal.y, 2 * bump * wall.normal.z);

  // move the molecule a bit so that it ends up at the correct side of the wall
  vm.v.pos = new_pos3d;
  vec3_t new_pos_after_diffuse;
  diffusion_util::tiny_diffuse_3D(*this, vm, displacement, wall.index, new_pos_after_diffuse);


  // TODO:
  // Make sure we didn't end up on a neighbor's wall, which is kind of easy to
  // do with, for example, a shrinking box/cuboid.
  // - see place_mol_relative_to_mesh


  // move the molecule and also update the subpartition information
  vm.v.pos = new_pos_after_diffuse;
  subpart_index_t new_subpart = get_subpartition_index(vm.v.pos);
  change_molecule_subpartition(vm, new_subpart);


#ifdef DEBUG_DYNAMIC_GEOMETRY
  vm.dump(*this, "", "Vm after being moved: ", simulation_stats.get_current_iteration /*iteration*/, 0);
#endif
}


bool is_point_above_plane_defined_by_wall(const Partition& p, const Wall& w, const vec3_t& pos) {

  // make vector pointing from any point to our position
  vec3_t w0_pos = pos - p.get_wall_vertex(w, 0);

  // dot product with normal gives ||a|| * ||b|| * cos(phi)
  float_t dot_prod = dot(w0_pos, w.normal);
  assert(!cmp_eq(dot_prod, 0, EPS) && "Checked point is on the plane");
  return dot_prod > 0;
}


/*************************************************************************
nearest_free:
  In: a surface grid
      a vector in u,v coordinates on that surface
      the maximum distance we can search for free spots
  Out: integer containing the index of the closest unoccupied grid point
       to the vector, or -1 if no unoccupied points are found in range
  Note: we assume you've already checked the grid element that contains
        the point, so we don't bother looking there first.
  Note: if no unoccupied tile is found, found_dist2 contains distance to
        closest occupied tile.
*************************************************************************/
// TODO: cleanup, move to grid_util
tile_index_t nearest_free(
    const Wall& wall, const vec2_t& v, const float_t max_d2,
    float_t& found_dist2) {

  const Grid& g = wall.grid;

  tile_index_t h;
  int i, j;
  uint k;
  int span;
  int can_flip;
  tile_index_t tile_index;
  float_t d2;
  float_t f, ff, fff;
  float_t over3n = 0.333333333333333 / (float_t)(g.num_tiles_along_axis);


  /* check whether the grid is fully occupied */
  if (g.is_full()) {
    found_dist2 = 0;
    return TILE_INDEX_INVALID;
  }

  tile_index = TILE_INDEX_INVALID;
  d2 = 2 * max_d2 + 1.0;

  // this seems quite inefficient

  for (k = 0; k < g.num_tiles_along_axis; k++) {
    f = v.v - ((float_t)(3 * k + 1)) * over3n * wall.uv_vert2.v;
    ff = f - over3n * wall.uv_vert2.v;
    ff *= ff;
    f *= f;
    if (f > max_d2 && ff > max_d2) {
      continue; /* Entire strip is too far away */
    }

    span = (g.num_tiles_along_axis - k);
    for (j = 0; j < span; j++) {
      can_flip = (j != span - 1);
      for (i = 0; i <= can_flip; i++) {
        fff =
            v.u - over3n * ((float_t)(3 * j + i + 1) * wall.uv_vert1_u +
                             (float_t)(3 * k + i + 1) * wall.uv_vert2.u);
        fff *= fff;
        if (i) {
          fff += ff;
        }
        else {
          fff += f;
        }

        if (fff < max_d2 && (tile_index == TILE_INDEX_INVALID || fff < d2)) {
          h = (g.num_tiles_along_axis - k) - 1;
          h = h * h + 2 * j + i;

          if (g.get_molecule_on_tile(h) == MOLECULE_ID_INVALID) {
            tile_index = h;
            d2 = fff;
          }
          else if (tile_index == TILE_INDEX_INVALID) {
            if (fff < d2) {
              d2 = fff;
            }
          }
        }
      }
    }
  }

  found_dist2 = d2;

  return tile_index;
}

/*************************************************************************
search_nbhd_for_free:
  In: the wall that we ought to be in
      a vector in u,v coordinates on that surface where we should go
      the maximum distance we can search for free spots
      a place to store the index of our free slot
      a function that we'll call to make sure a wall is okay
      context for that function passed in by whatever called us
  Out: pointer to the wall that has the free slot, or NULL if no wall
       exist in range.
  Note: usually the calling function will create a grid if needed and
        check the grid element at u,v; if that is not done this function
        will return the correct result but not efficiently.
  Note: This is not recursive.  It should be made recursive.
*************************************************************************/
void search_nbhd_for_free(
    Partition& p,
    const wall_index_t origin_wall_index, const vec2_t& closest_pos2d, const float_t max_search_d2,
    wall_index_t& res_wall_index, tile_index_t& res_tile_index
  ) {

  wall_index_t best_wall_index = origin_wall_index;
  tile_index_t best_tile_index = TILE_INDEX_INVALID;
  vec2_t best_pos = closest_pos2d;

  const Wall& origin_wall = p.get_wall(origin_wall_index);
  assert(origin_wall.grid.is_initialized());

  // Find index and distance of nearest free grid element on origin wall,
  // returns TILE_INDEX_INVALID when there is not space left
  float_t d2_unused;
  best_tile_index = nearest_free(origin_wall, best_pos, max_search_d2, d2_unused);

  if (best_tile_index != TILE_INDEX_INVALID) {
    res_wall_index = best_wall_index;
    res_tile_index = best_tile_index;
    return;
  }

  float_t best_d2 = 2.0 * max_search_d2 + 1.0;
  const vec2_t& point = closest_pos2d;

  /* if there are no free slots on the origin wall - look around */
  /* Check for closer free grid elements on neighboring walls */
  for (edge_index_t j = 0; j < EDGES_IN_TRIANGLE; j++) {
    if (origin_wall.edges[j].backward_index == WALL_INDEX_INVALID)
      continue;

    wall_index_t there_wall_index;
    if (origin_wall.edges[j].forward_index == origin_wall_index) {
      there_wall_index = origin_wall.edges[j].backward_index;
    }
    else {
      there_wall_index = origin_wall.edges[j].forward_index;
    }

    Wall& there_wall = p.get_wall(there_wall_index);

    // TODO - wall regions - we are checking only geom object for now
    if (!origin_wall.is_same_region(there_wall)) {
      continue;
    }

    /* check whether there are any available spots on the neighbor wall */

    if (there_wall.grid.is_initialized()) {
      if (there_wall.grid.is_full()) {
        continue;
      }
    }

    /* Calculate distance between point and edge j of origin wall */
    vec2_t vurt0, vurt1;
    switch (j) {
    case 0:
      vurt0 = vec2_t(0);
      vurt1.u = origin_wall.uv_vert1_u;
      vurt1.v = 0;
      break;
    case 1:
      vurt0.u = origin_wall.uv_vert1_u;
      vurt0.v = 0;
      vurt1 = origin_wall.uv_vert2;
      break;
    case 2:
      vurt0 = origin_wall.uv_vert2;
      vurt1 = vec2_t(0);
      break;
    default:
      /* default case should not occur since 0<=j<=2 */
      assert(false);
    }

    vec2_t pt, ed;
    ed = vurt1 - vurt0;
    pt = point - vurt0;

    float_t d2;
    d2 = dot2(pt, ed);
    d2 = len2_squared(pt) -
         d2 * d2 / len2_squared(ed); /* Distance squared to line */

    /* Check for free grid element on neighbor if point to edge distance is
     * closer than best_d2  */
    if (d2 < best_d2) {

      if (!there_wall.grid.is_initialized()) {
        there_wall.grid.initialize(p, there_wall);
      }

      GeometryUtil::traverse_surface(origin_wall, point, j, pt);
      tile_index_t i = nearest_free(there_wall, pt, max_search_d2, d2);

      if (i != TILE_INDEX_INVALID && d2 < best_d2) {
        best_tile_index = i;
        best_d2 = d2;
        best_wall_index = there_wall_index;
      }
    }
  }

  res_wall_index = best_wall_index;
  res_tile_index = best_tile_index;
}


void find_closest_tile_on_wall(
    Partition& p,
    const wall_index_t closest_wall_index, const vec2_t& closest_pos2d,
    const float_t closest_d2, const float_t search_d2,
    wall_index_t& found_wall_index, tile_index_t& found_tile_index, vec2_t& found_pos2d
) {

  tile_index_t closest_tile_index;

  const Wall& w = p.get_wall(closest_wall_index);
  closest_tile_index = GridUtil::uv2grid_tile_index(closest_pos2d, w);

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  cout << "find_closest_tile_on_wall: closest_wall_index: " << closest_wall_index <<
      ", closest_tile_index: " << closest_tile_index << "\n";
  w.grid.dump();
#endif

  molecule_id_t mol_on_tile = w.grid.get_molecule_on_tile(closest_tile_index);


  if (mol_on_tile == MOLECULE_ID_INVALID) {
    // ok, tile is empty
    found_wall_index = closest_wall_index;
    found_tile_index = closest_tile_index;
    found_pos2d = closest_pos2d;
  }
  else {
    // need to find a tile that is close

    // squared distance
    float_t max_search_d2 = search_d2 - closest_d2;
    if (max_search_d2 <= EPS_C * EPS_C) {
      mcell_error("Search distance for find_closest_tile_on_wall is too small");
    }

    wall_index_t new_wall_index;
    wall_index_t new_tile_index;
    search_nbhd_for_free(
        p, closest_wall_index, closest_pos2d, max_search_d2,
        new_wall_index, new_tile_index
    );
    assert(new_wall_index != TILE_INDEX_INVALID);

    vec2_t new_pos2d;
    const Wall& w = p.get_wall(new_wall_index);
    if (p.get_world_constants().randomize_smol_pos) {
      // molecules are processed in different order than in mcell3...
      assert(false && "TODO - need RNG for dynamic vertices and surface molecules");
      // new_pos2d = GridUtil::grid2uv_random(w, new_tile_index, rng); // XXX
    }
    else {
      new_pos2d = GridUtil::grid2uv(w, new_tile_index);
    }

    found_wall_index = new_wall_index;
    found_tile_index = new_tile_index; //GridUtil::uv2grid_tile_index(closest_pos2d, w);
    found_pos2d = new_pos2d;
  }
}

// insert_surface_molecule & place_surface_molecule in mcell3
void Partition::move_surface_molecule_to_closest_wall_point(
    const SurfaceMoleculeMoveInfo& molecule_move_info) {

  Molecule& sm = get_m(molecule_move_info.molecule_id);
  assert(sm.is_surf());

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  sm.dump(*this, "", "Moving sm towards new wall: ", simulation_stats.get_current_iteration, 0);
#endif

  // 1) check all walls to get a reference hopefully
  //  then try to limit only to wall's neighbors - we would really like to avoid any regions (maybe...?)
  wall_index_t best_wall_index;
  vec2_t best_wall_pos2d;
  float_t best_d2 = find_closest_wall(molecule_move_info.pos3d, molecule_move_info.wall_index, best_wall_index, best_wall_pos2d);

  if (best_wall_index == WALL_INDEX_INVALID) {
    mcell_error("Could not find a wall close to volume molecule with id %d.\n", sm.id);
  }

  // place molecule onto the found wall, all the remaining information about the molecule stays the same
  wall_index_t found_wall_index;
  tile_index_t found_tile_index;
  vec2_t found_pos2d;

  find_closest_tile_on_wall(
      *this, best_wall_index, best_wall_pos2d, best_d2, get_world_constants().vacancy_search_dist2,
      found_wall_index, found_tile_index, found_pos2d
  );
  assert(found_wall_index != WALL_INDEX_INVALID);
  assert(found_tile_index != TILE_INDEX_INVALID);

  sm.s.wall_index = found_wall_index;
  sm.s.grid_tile_index = found_tile_index;
  sm.s.pos = found_pos2d; // this will need to be fixed
  Wall& w = get_wall(found_wall_index);
  w.grid.set_molecule_tile(found_tile_index, sm.id);

  // TODO reschedule unimolar

#ifdef DEBUG_DYNAMIC_GEOMETRY
  sm.dump(*this, "", "Sm after being moved: ", simulation_stats.get_current_iteration /*iteration*/, 0);
#endif
}


// TODO: move to dyn_vertex_utils? -> probably yes
void Partition::collect_volume_molecules_moved_due_to_moving_wall(
    const wall_index_t moved_wall_index, const VertexMoveInfoVector& move_infos,
    // TODO: maybe remove this set already_moved_molecules since we already have the molecules in the vector molecule_moves?
    UintSet& already_moved_molecules,
    VolumeMoleculeMoveInfoVector& molecule_moves
) {

  Wall& orig_wall = get_wall(moved_wall_index);

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  cout << "*** Moving vertices of a wall:\n";
  orig_wall.dump(*this, "", false);
  for (const VertexMoveInfo& info: move_infos) {
    cout << "vertex index: " << info.vertex_index << "\n";
    cout << "original position: " << get_geometry_vertex(info.vertex_index) << "\n";
    cout << "translation: " << info.translation_vec << "\n";
  }
  cout << "***\n";
#endif

  assert(move_infos.size() > 0 && move_infos.size() <= 3 && "Move infos belong to the wall that is being moved");

  // construct a virtual space where we have all the walls with new and old
  // positions

  const vertex_index_t* orig_indices = orig_wall.vertex_indices;

  // create 3 new temporary vertices and insert them into the partition,
  // we will create them and then erase, so we will create a new vertex even if
  // it is the same one
  // we will need them for our arbitrary walls
  vertex_index_t new_indices[VERTICES_IN_TRIANGLE];
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {

    // copy position of the original vertex
    vec3_t new_vertex = get_geometry_vertex(orig_indices[i]);

    // should we move it?
    bool moved = false;
    for (const VertexMoveInfo& move_info: move_infos) {
      if (move_info.vertex_index == orig_indices[i]) {
        assert(!moved);
        new_vertex = new_vertex + move_info.translation_vec;
        #ifndef NDEBUG
          break;
        #endif
      }
    }

    new_indices[i] = add_geometry_vertex(new_vertex);
  }

  const vec3_t* o[VERTICES_IN_TRIANGLE] = {
      &get_geometry_vertex(orig_indices[0]),
      &get_geometry_vertex(orig_indices[1]),
      &get_geometry_vertex(orig_indices[2])
  };

  const vec3_t* n[VERTICES_IN_TRIANGLE] = {
      &get_geometry_vertex(new_indices[0]),
      &get_geometry_vertex(new_indices[1]),
      &get_geometry_vertex(new_indices[2])
  };


  // opposite triangle (wall after being moved)
  // first true argument - we need the wall constants to be precomputed,
  // second false argument - we do not care about edge information
  Wall new_wall = Wall(*this, new_indices[0], new_indices[1], new_indices[2], true, false);


  // TODO: code below can be moved to a separate function to detect whether a point is in a mesh
#ifdef DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS
  // TODO: move to some 'dump utils'
  // script mcell/utils/blender_debug_scripts/dyn_vertex_check.py can be used to visualize the
  // collision detection
  cout << "Constructed object from initial and final triangle:\n";

  // dump the object in a form that can be imported to blender
  map<vertex_index_t, uint> map_assigned_index;
  cout << "mesh_verts = [\n";
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    cout << "  " << get_geometry_vertex(orig_indices[i]) << ", #" << i << "\n";
    map_assigned_index[orig_indices[i]] = i;
  }
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    cout << "  " << get_geometry_vertex(new_indices[i]) << ", #" << i + VERTICES_IN_TRIANGLE << "\n";
    map_assigned_index[new_indices[i]] = i + VERTICES_IN_TRIANGLE;
  }
  cout << "]\n";

  // dumping only the orig and new wall
  cout << "mesh_faces = [\n";
  cout << "  (";
  for (uint k = 0; k < VERTICES_IN_TRIANGLE; k++) {
    cout << map_assigned_index[ orig_wall.vertex_indices[k] ] << ", ";
  }
  cout << "),\n";
  cout << "  (";
  for (uint k = 0; k < VERTICES_IN_TRIANGLE; k++) {
    cout << map_assigned_index[new_wall.vertex_indices[k] ] << ", ";
  }
  cout << "),\n";


  // FIXME": how to display the interconnections as well?
  /*cout << "mesh_faces = [\n";
  for (uint i = 0; i < TRIANGLES_IN_MOVED_TRIANGLE_MESH; i++) {
    cout << "  (";
    for (uint k = 0; k < VERTICES_IN_TRIANGLE; k++) {
      cout << map_assigned_index[ moved_triangle_walls[i]->vertex_indices[k] ] << ", ";
    }
    cout << "),\n";
  }
  */
  cout << "]\n";
  cout << "add_mesh(mesh_verts, mesh_faces, \"my_mesh\")\n";
#endif


  // if moving by one edge creates just a triangle, store this information
  bool egde_moved[EDGES_IN_TRIANGLE];
  Wall* wall_if_edge_defines_triangle[EDGES_IN_TRIANGLE];

  for (uint i1 = 0; i1 < EDGES_IN_TRIANGLE; i1++) {
    uint i2 = (i1 + 1) % EDGES_IN_TRIANGLE;

    bool v1_same = *o[i1] == *n[i1];
    bool v2_same = *o[i2] == *n[i2];

    if (v1_same && v2_same) {
      egde_moved[i1] = false;
      wall_if_edge_defines_triangle[i1] = nullptr;
    }
    else if (v1_same) {
      egde_moved[i1] = true;
      wall_if_edge_defines_triangle[i1] = new Wall(*this, orig_indices[i1], orig_indices[i2], new_indices[i2], true, false);
    }
    else if (v2_same) {
      egde_moved[i1] = true;
      wall_if_edge_defines_triangle[i1] = new Wall(*this, orig_indices[i1], orig_indices[i2], new_indices[i1], true, false);
    }
    else {
      egde_moved[i1] = true;
      wall_if_edge_defines_triangle[i1] = nullptr;
    }
  }

  // TODO: optimize only for molecules in relevant subpartitions
  for (Molecule& m: molecules) {

    if (m.is_surf()) {
      continue;
    }

    if (already_moved_molecules.count(m.id) == 1) {
      //assert(false && "this should not happen anymore?");
      continue;
    }

#ifdef DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS
    // cout << "# Detecting collision for molecule with id " << m.id << " at " << m.v.pos << "\n";
    // TODO: move to some 'dump utils'
    cout << "mol" << m.id << " = " << m.v.pos << "\n";
    cout << "add_point(mol" << m.id << ", \"molecule" << m.id << "\")\n";
#endif

    // now check a single projection against all walls created by these two triangles
    // if the number of crosses is odd, then we are inside
    int num_hits = 0;

    // cast ray along the whole partition
    vec3_t move(0, 0, get_world_constants().partition_edge_length);

    // check collision with the original wall
    bool collides = CollisionUtil::collide_wall_test(*this, m.v.pos, orig_wall, move);
    if (collides) { num_hits++; }

    // and with the new wall as well
    collides = CollisionUtil::collide_wall_test(*this, m.v.pos, new_wall, move);
    if (collides) { num_hits++; }

    // now, let's deal with the areas that are 'drawn' by the moving edges
    // NOTE: many of the values can be pre-computed, but let's keep it simple for now
    for (uint i1 = 0; i1 < EDGES_IN_TRIANGLE; i1++) {
      uint i2 = (i1 + 1) % EDGES_IN_TRIANGLE;

      collides = CollisionUtil::collide_moving_line_and_static_line_test(
          *this,
          m.v.pos, move,
          *o[i1], *n[i1], *o[i2], *n[i2],
          egde_moved[i1], wall_if_edge_defines_triangle[i1]
      );
      if (collides) {
        num_hits++;
      }
    }


    if (num_hits % 2 == 1) {

      // we are moving with a single wall here.
      // different from MCell3 behavior where it tries to find the closest point on the mesh with a given name
      // this might be super expensive if the mesh is large

      // first we need to figure out on which side of the new wall we should place the molecule
      // with regards to its normal
      bool place_above = is_point_above_plane_defined_by_wall(*this, orig_wall, m.v.pos);

      molecule_moves.push_back(VolumeMoleculeMoveInfo(m.id, orig_wall.index, place_above));

      // and remember that we must not be moving it anymore
      // should work even without it...
      already_moved_molecules.insert_unique(m.id);

    }
  }

  // free added walls (wall with index 0 existed before)
  for (uint i = 1; i < EDGES_IN_TRIANGLE; i++) {
    if (wall_if_edge_defines_triangle[i] != nullptr) {
      delete wall_if_edge_defines_triangle[i];
    }
  }

  // and remove added vertices (checking that we are removing the right ones)
  for (int i = VERTICES_IN_TRIANGLE - 1; i >= 0; i--) {
    remove_last_vertex(new_indices[i]);
  }

}


// TODO: move to dyn_vertex_utils? -> probably yes
void Partition::collect_surface_molecules_moved_due_to_moving_wall(
    const wall_index_t moved_wall_index,
    SurfaceMoleculeMoveInfoVector& molecule_moves) {

  // get all surface molecules that belong to a given wall

  // We have no direct set that would say which molecules belong to a given wall,
  // attempts were done, but even maintaining a simple set in grid costs ~5% of runtime.
  // With the data that we have available, the simples option is to
  // go through the array in grid and check all existing molecules.
  // There might be other solutions that use subpartitions but let's keep it simple for now.
  const Grid* g = get_wall_grid_if_exists(moved_wall_index);
  if (g == nullptr) {
    return;
  }
  small_vector<molecule_id_t> molecule_ids;
  g->get_contained_molecules(molecule_ids);

  const Wall& w = get_wall(moved_wall_index);
  const vec3_t& vert0 = get_geometry_vertex(w.vertex_indices[0]);

  for (molecule_id_t id: molecule_ids) {

    vec3_t pos3d = GeometryUtil::uv2xyz(get_m(id).s.pos, w, vert0);

    molecule_moves.push_back(
        SurfaceMoleculeMoveInfo(id, moved_wall_index, pos3d)
    );
  }
}



void Partition::apply_vertex_moves() {
  // 1) create a set of all affected walls with information on how much each wall moves,
  UintSet moved_vertices_set;
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
  UintSet already_moved_volume_molecules;
  SurfaceMoleculeMoveInfoVector surface_molecule_moves;

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  cout << "*** Walls being moved:\n";
  for (const auto& it: walls_with_their_moves) {
    const Wall& w = get_wall(it.first);
    w.dump(*this, "  ", true);
  }
#endif

  for (const auto& it: walls_with_their_moves) {
    collect_volume_molecules_moved_due_to_moving_wall(it.first, it.second, already_moved_volume_molecules, volume_molecule_moves);

    collect_surface_molecules_moved_due_to_moving_wall(it.first, surface_molecule_moves);
  }

  // 3) get information on where these walls are and remove them
  update_walls_per_subpart(walls_with_their_moves, false);

  // 4) then we move the vertices and update relevant walls
  Geometry::update_moved_walls(*this, scheduled_vertex_moves, walls_with_their_moves);

  // 5) update subpartition info for the walls
  update_walls_per_subpart(walls_with_their_moves, true);

  // 6) move volume molecules
  for (const VolumeMoleculeMoveInfo& move_info: volume_molecule_moves) {
    move_volume_molecule_to_closest_wall_point(move_info);
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
    // we need original xyz position of the surface mols
    // TODO: the new wall must belong to the same object
    move_surface_molecule_to_closest_wall_point(move_info);
  }

  scheduled_vertex_moves.clear();
}


void Partition::dump() {
  for (GeometryObject& obj:geometry_objects) {
    obj.dump(*this, "  ");
  }

  for (size_t i = 0; i < walls_per_subpart.size(); i++) {
    if (!walls_per_subpart[i].empty()) {
      vec3_t llf, urb;
      get_subpart_llf_point(i, llf);
      urb = llf + vec3_t(world_constants.subpartition_edge_length);

      cout << "subpart: " << i << ", llf: " << llf << ", urb: " << urb << "\n";
      walls_per_subpart[i].dump();
    }
  }
}

} // namespace mcell
