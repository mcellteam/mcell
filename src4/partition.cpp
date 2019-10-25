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

#include "partition.h"

#include "dyn_vertex_utils.h"
#include "geometry_utils.inc"
#include "collision_utils.inc"

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


// TODO: move to dyn_vertex_utils? -> probably yes
void Partition::move_molecules_due_to_moving_wall(const wall_index_t moved_wall_index, const VertexMoveInfoVector& move_infos) {

  // construct a virtual space where we have all the walls with new and old
  // positions

  Wall& orig_wall = get_wall(moved_wall_index);
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

  // create new arbitrary walls and initialize them,
  // we do not need to insert them into the partition
  const uint TRIANGLES_IN_MOVED_TRIANGLE_MESH = 8;
  Wall* moved_triangle_walls[TRIANGLES_IN_MOVED_TRIANGLE_MESH];
  moved_triangle_walls[0] = &orig_wall;

  vertex_index_t o0 = orig_indices[0], o1 = orig_indices[1], o2 = orig_indices[2];
  vertex_index_t n0 = new_indices[0], n1 = new_indices[1], n2 = new_indices[2];

  // opposite triangle (wall after being moved)
  // first true argument - we need the wall constants to be precomputed,
  // second false argument - we do not care about edge information
  moved_triangle_walls[1] = new Wall(*this, n0, n1, n2, true, false);

  // triangles that connect the orig and new wall
  moved_triangle_walls[2] = new Wall(*this, o0, o1, n0, true, false);
  moved_triangle_walls[3] = new Wall(*this, n0, o1, n1, true, false);

  moved_triangle_walls[4] = new Wall(*this, o1, o2, n1, true, false);
  moved_triangle_walls[5] = new Wall(*this, n1, o2, n2, true, false);

  moved_triangle_walls[6] = new Wall(*this, o2, o0, n2, true, false);
  moved_triangle_walls[7] = new Wall(*this, n2, o0, n0, true, false);

#ifdef DEBUG_DYNAMIC_GEOMETRY
  // TODO: move to some 'dump utils'
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
  cout << "mesh_faces = [\n";
  for (uint i = 0; i < TRIANGLES_IN_MOVED_TRIANGLE_MESH; i++) {
    cout << "  (";
    for (uint k = 0; k < VERTICES_IN_TRIANGLE; k++) {
      cout << map_assigned_index[ moved_triangle_walls[i]->vertex_indices[k] ] << ", ";
    }
    cout << "),\n";
  }
  cout << "]\n";
  cout << "add_mesh(mesh_verts, mesh_faces, \"my_mesh\")\n";
#endif

  // TODO: optimize only for molecules in relevant subpartitions
  for (const Molecule& m: molecules) {

    if (m.is_surf()) {
      continue;
    }

#ifdef DEBUG_DYNAMIC_GEOMETRY
    // cout << "# Detecting collision for molecule with id " << m.id << " at " << m.v.pos << "\n";
    // TODO: move to some 'dump utils'
    cout << "mol" << m.id << " = " << m.v.pos << "\n";
    cout << "add_point(mol" << m.id << ", \"molecule" << m.id << "\")\n";
#endif

    // now check a single projection against all walls created by these two triangles
    // if the number of crosses is odd, then we are inside
    int num_hits = 0;
    rng_state unused_rng_state;
    float_t ignored_collision_time;
    vec3_t ignored_collision_pos;

    // cast ray along the whole partition
    vec3_t move(0, 0, get_world_constants().partition_edge_length);

    for (uint i = 0; i < TRIANGLES_IN_MOVED_TRIANGLE_MESH; i++) {

      // TODO: skip degenerate triangles
      CollisionType res = CollisionUtil::collide_wall(
          *this, m.v.pos, *moved_triangle_walls[i],
          unused_rng_state, false,
          move,
          ignored_collision_time, ignored_collision_pos
      );

      switch (res) {
        case CollisionType::WALL_MISS:
          break;
        case CollisionType::WALL_FRONT:
        case CollisionType::WALL_BACK:
          #ifdef DEBUG_DYNAMIC_GEOMETRY
            cout << "# Detecting collision for molecule with id " << m.id <<
              " at " << ignored_collision_pos << " with wall " << i <<"\n";

          #endif
          num_hits++;
          break;
        case CollisionType::WALL_REDO:
          mcell_error("Collision REDO is not handled yet in dynamic vertices.");
          break;
        default:
          assert(false);
      }
    }

    if (num_hits % 2 == 1) {
      // TODO: move molecule to a new position close to the
      mcell_error("Collision detected!:)");
    }
  }

  // free added walls (wall with index 0 existed before)
  for (uint i = 1; i < TRIANGLES_IN_MOVED_TRIANGLE_MESH; i++) {
    delete moved_triangle_walls[i];
  }

  // and remove added vertices (checking that we are removing the right ones)
  for (int i = VERTICES_IN_TRIANGLE - 1; i >= 0; i--) {
    remove_last_vertex(new_indices[i]);
  }

}


void Partition::apply_vertex_moves() {

  // TODO: expecting that there we are not moving a single vertex twice, add a check for it

  // 1) create a set of all affected walls with information on how much each wall moves,
  WallsWithTheirMovesMap walls_with_their_moves;
  for (const VertexMoveInfo& move_info: scheduled_vertex_moves) {
    const std::vector<wall_index_t>& wall_indices = get_walls_using_vertex(move_info.vertex_index);

    for (wall_index_t wall_index: wall_indices) {
      // remember mapping wall_index -> moves
      auto it = walls_with_their_moves.find(wall_index);
      if (it == walls_with_their_moves.end()) {
        it = walls_with_their_moves.insert(make_pair(wall_index, VertexMoveInfoVector())).first;
      }
      it->second.push_back(move_info);
    }
  }

  // 2) for each wall, detect what molecules will be moved and move them right away
  for (auto it : walls_with_their_moves) {
    move_molecules_due_to_moving_wall(it.first, it.second);
  }

  // 3) get information on where these walls are and remove them
  update_walls_per_subpart(walls_with_their_moves, false);

  // 4) then we move the vertices and update relevant walls
  DynVertexUtils::move_vertices_and_update_walls(*this, scheduled_vertex_moves, walls_with_their_moves);

  // 5) and update subpartition info for the walls
  update_walls_per_subpart(walls_with_their_moves, true);


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
