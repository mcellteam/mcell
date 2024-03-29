/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_DYN_VERTEX_UTILS_INC_
#define SRC4_DYN_VERTEX_UTILS_INC_

#include <iostream>

#include "logging.h"

#include "partition.h"
#include "dyn_vertex_structs.h"
#include "geometry.h"
#include "grid_utils.inl"
#include "geometry_utils.inl"
#include "diffusion_utils.inl"
#include "collision_utils.inl"
#include "wall_utils.inl"

using namespace std;

namespace MCell {

namespace DynVertexUtils {


// based on place_mol_relative_to_mesh
static void move_volume_molecule_to_closest_wall_point(Partition& p, const VolumeMoleculeMoveInfo& molecule_move_info) {

  Molecule& vm = p.get_m(molecule_move_info.molecule_id);
  assert(vm.is_vol());

  // 1) check all walls to get a reference hopefully
  //  then try to limit only to wall's neighbors - we would really like to avoid any regions (maybe...?)
  wall_index_t best_wall_index;
  Vec2 best_wall_pos2d;
  WallUtils::find_closest_wall(p, vm.v.pos, molecule_move_info.wall_index, false, best_wall_index, best_wall_pos2d);

  if (best_wall_index == WALL_INDEX_INVALID) {
    mcell_error("Could not find a wall close to volume molecule with id %d.\n", vm.id);
  }
  Wall& wall = p.get_wall(best_wall_index);

  Vec3 new_pos3d = GeometryUtils::uv2xyz(best_wall_pos2d, wall, p.get_wall_vertex(wall, 0));

#ifdef DEBUG_DYNAMIC_GEOMETRY
  vm.dump(p, "", "Moving vm towards new wall: ", p.stats.get_current_iteration(), 0);
  wall.dump(p, "", true);
#endif

  // displacement
  pos_t bump = (molecule_move_info.place_above) ? POS_EPS : -POS_EPS;
  Vec3 displacement(2 * bump * wall.normal.x, 2 * bump * wall.normal.y, 2 * bump * wall.normal.z);

  // move the molecule a bit so that it ends up at the correct side of the wall
  vm.v.pos = new_pos3d;
  vm.v.subpart_index = p.get_subpart_index(vm.v.pos);
  Vec3 new_pos_after_diffuse;
  DiffusionUtils::tiny_diffuse_3D(p, vm, displacement, wall.index, new_pos_after_diffuse);

  // TODO:
  // Make sure we didn't end up on a neighbor's wall, which is kind of easy to
  // do with, for example, a shrinking box/cuboid.
  // - see place_mol_relative_to_mesh

  // move the molecule and also update the information on subpartition reactants
  vm.v.pos = new_pos_after_diffuse;
  vm.v.subpart_index = p.get_subpart_index(vm.v.pos);
  p.update_molecule_reactants_map(vm);

#ifdef DEBUG_DYNAMIC_GEOMETRY
  vm.dump(p, "", "Vm after being moved: ", p.stats.get_current_iteration() /*iteration*/, 0);
#endif
}


// insert_surface_molecule & place_surface_molecule in mcell3
// TODO: quite similar code as in WallUtils::place_surface_molecule,
// can be it be somehow merged?
static void move_surface_molecule_to_closest_wall_point(
    Partition& p,
    const SurfaceMoleculeMoveInfo& molecule_move_info) {

  Molecule& sm = p.get_m(molecule_move_info.molecule_id);
  assert(sm.is_surf());

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  sm.dump(p, "", "Moving sm towards new wall: ", p.stats.get_current_iteration(), 0);
#endif

  // 1) check all walls to get a reference hopefully
  //  then try to limit only to wall's neighbors - we would really like to avoid any regions (maybe...?)
  wall_index_t best_wall_index;
  Vec2 best_wall_pos2d;
  pos_t best_d2 = WallUtils::find_closest_wall(
      p, molecule_move_info.pos3d, molecule_move_info.wall_index, true, best_wall_index, best_wall_pos2d);

  if (best_wall_index == WALL_INDEX_INVALID) {
    mcell_error("Could not find a wall close to volume molecule with id %d.\n", sm.id);
  }

  // place molecule onto the found wall, all the remaining information about the molecule stays the same
  wall_index_t found_wall_index;
  tile_index_t found_tile_index;
  Vec2 found_pos2d(FLT_INVALID);

  GridUtils::find_closest_tile_on_wall(
      p, best_wall_index, best_wall_pos2d, best_d2, p.config.vacancy_search_dist2,
      found_wall_index, found_tile_index, found_pos2d
  );
  assert(found_wall_index != WALL_INDEX_INVALID);
  assert(found_tile_index != TILE_INDEX_INVALID);
  assert(found_pos2d != Vec2(FLT_INVALID));

  sm.s.wall_index = found_wall_index;
  sm.s.grid_tile_index = found_tile_index;
  sm.s.pos = found_pos2d;
  Wall& w = p.get_wall(found_wall_index);
  w.grid.set_molecule_tile(found_tile_index, sm.id);

  // TODO: reschedule unimolar

#ifdef DEBUG_DYNAMIC_GEOMETRY
  sm.dump(p, "", "Sm after being moved: ", p.stats.get_current_iteration() /*iteration*/, 0);
#endif
}


namespace Local {
static void dump_blender_display_code(
    const Partition& p,
    const Wall& orig_wall, const vertex_index_t* orig_indices,
    const Wall& new_wall, const vertex_index_t* new_indices) {

  // NOTE: move to some 'dump utils'?
  // script mcell/utils/blender_debug_scripts/dyn_vertex_check.py can be used to visualize the
  // collision detection
  cout << "Constructed object from initial and final triangle:\n";

  // dump the object in a form that can be imported to blender
  map<vertex_index_t, uint> map_assigned_index;
  cout << "mesh_verts = [\n";
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    cout << "  " << p.get_geometry_vertex(orig_indices[i]) << ", #" << i << "\n";
    map_assigned_index[orig_indices[i]] = i;
  }
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    cout << "  " << p.get_geometry_vertex(new_indices[i]) << ", #" << i + VERTICES_IN_TRIANGLE << "\n";
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

  // NOTE: how to display the interconnections as well?
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
}
} // namespace Local


static void collect_volume_molecules_moved_due_to_moving_wall(
    Partition& p,
    const wall_index_t moved_wall_index, const WallMoveInfo& move_info,
    MoleculeIdsSet& already_moved_molecules,
    VolumeMoleculeMoveInfoVector& molecule_moves
) {

  const Wall& orig_wall = p.get_wall(moved_wall_index);

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  cout << "*** Moving vertices of a wall:\n";
  orig_wall.dump(p, "", false);
  for (const VertexMoveInfo& info: move_info) {
    cout << "vertex index: " << info.vertex_index << "\n";
    cout << "original position: " << p.get_geometry_vertex(info.vertex_index) << "\n";
    cout << "translation: " << info.displacement << "\n";
  }
  cout << "***\n";
#endif

  assert(move_info.vertex_moves.size() > 0 && move_info.vertex_moves.size() <= 3 && "Move infos belong to the wall that is being moved");

  // construct a virtual space where we have all the walls with new and old
  // positions

  const vertex_index_t* orig_indices = orig_wall.vertex_indices;

  // create 3 new temporary vertices
  Vec3 n[VERTICES_IN_TRIANGLE];
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {

    // copy position of the original vertex
    Vec3 new_vertex = p.get_geometry_vertex(orig_indices[i]);

    // should we move it?
    bool moved = false;
    for (const VertexMoveInfo* vertex_move: move_info.vertex_moves) {
      if (vertex_move->vertex_index == orig_indices[i]) {
        assert(!moved);
        new_vertex = new_vertex + vertex_move->displacement;
        #ifndef NDEBUG
          break;
        #endif
      }
    }

    n[i] = new_vertex;
  }

  const Vec3* o[VERTICES_IN_TRIANGLE] = {
      &p.get_geometry_vertex(orig_indices[0]),
      &p.get_geometry_vertex(orig_indices[1]),
      &p.get_geometry_vertex(orig_indices[2])
  };

  // opposite triangle (wall after being moved)
  // first true argument - we need the wall constants to be precomputed,
  // second false argument - we do not care about edge information
  WallWithVertices new_wall = WallWithVertices(p, n[0], n[1], n[2], true);

  #ifdef DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS
    Local::dump_blender_display_code(p, orig_wall, orig_indices, new_wall, new_indices);
  #endif

  // if moving by one edge creates just a triangle, store this information
  bool egde_moved[EDGES_IN_TRIANGLE];
  Wall* wall_if_edge_defines_triangle[EDGES_IN_TRIANGLE];

  for (uint i1 = 0; i1 < EDGES_IN_TRIANGLE; i1++) {
    uint i2 = (i1 + 1) % EDGES_IN_TRIANGLE;

    bool v1_same = *o[i1] == n[i1];
    bool v2_same = *o[i2] == n[i2];

    if (v1_same && v2_same) {
      egde_moved[i1] = false;
      wall_if_edge_defines_triangle[i1] = nullptr;
    }
    else if (v1_same) {
      egde_moved[i1] = true;
      wall_if_edge_defines_triangle[i1] = new WallWithVertices(p, *o[i1], *o[i2], n[i2], true);
    }
    else if (v2_same) {
      egde_moved[i1] = true;
      wall_if_edge_defines_triangle[i1] = new WallWithVertices(p, *o[i1], *o[i2], n[i1], true);
    }
    else {
      egde_moved[i1] = true;
      wall_if_edge_defines_triangle[i1] = nullptr;
    }
  }

  // TODO DYN GEOM: optimize only for molecules in relevant subpartitions
  std::vector<Molecule>& molecules = p.get_molecules();
  for (Molecule& m: molecules) {
    if (m.is_defunct()) {
      continue;
    }
    if (m.is_surf()) {
      continue;
    }

    if (already_moved_molecules.count(m.id) == 1) {
      //assert(false && "this should not happen anymore?");
      continue;
    }

#ifdef DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS
    cout << "mol" << m.id << " = " << m.v.pos << "\n";
    cout << "add_point(mol" << m.id << ", \"molecule" << m.id << "\")\n";
#endif

    // now check a single projection against all walls created by these two triangles
    // if the number of crosses is odd, then we are inside
    int num_hits = 0;

    // cast ray along the whole partition
    Vec3 move(0, 0, p.config.partition_edge_length);

    // check collision with the original wall
    bool collides = CollisionUtils::collide_wall_test(p, orig_wall, m.v.pos, move);
    if (collides) { num_hits++; }

    // and with the new wall as well, the wall does not exist in the partition
    collides = CollisionUtils::collide_wall_test(p, new_wall, m.v.pos, move);
    if (collides) { num_hits++; }

    // now, let's deal with the areas that are 'drawn' by the moving edges
    // NOTE: many of the values can be pre-computed, but let's keep it simple for now
    for (uint i1 = 0; i1 < EDGES_IN_TRIANGLE; i1++) {
      uint i2 = (i1 + 1) % EDGES_IN_TRIANGLE;

      collides = CollisionUtils::collide_moving_line_and_static_line_test(
          p,
          m.v.pos, move,
          *o[i1], n[i1], *o[i2], n[i2],
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
      bool place_above = GeometryUtils::is_point_above_plane_defined_by_wall(p, orig_wall, m.v.pos);

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
}


void collect_surface_molecules_moved_due_to_moving_wall(
    const Partition& p,
    const wall_index_t moved_wall_index,
    SurfaceMoleculeMoveInfoVector& molecule_moves,
    MoleculeIdsVector& paired_molecules) {

  // We have no direct set that would say which molecules belong to a given wall,
  // attempts were done, but even maintaining a simple set in grid costs ~5% of runtime.
  // With the data that we have available, the simples option is to
  // go through the array in grid and check all existing molecules.
  // There might be other solutions that use subpartitions but let's keep it simple for now.
  const Grid* g = p.get_wall_grid_if_exists(moved_wall_index);
  if (g == nullptr) {
    return;
  }
  small_vector<molecule_id_t> molecule_ids;
  g->get_contained_molecules(molecule_ids);

  const Wall& w = p.get_wall(moved_wall_index);
  const Vec3& vert0 = p.get_geometry_vertex(w.vertex_indices[0]);

  for (molecule_id_t id: molecule_ids) {
    const Molecule& sm = p.get_m(id);
    assert(sm.is_surf());
    assert(!sm.is_defunct());
    Vec3 pos3d = GeometryUtils::uv2xyz(sm.s.pos, w, vert0);

    molecule_moves.push_back(
        SurfaceMoleculeMoveInfo(id, moved_wall_index, pos3d)
    );

    if (p.get_paired_molecule(id) != MOLECULE_ID_INVALID) {
      paired_molecules.push_back(id);
    }
  }
}

} // namespace DynVertexUtil

} // namespace MCell

#endif // SRC4_DYN_VERTEX_UTILS_INC_
