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

#include "dyn_vertex_utils.inc"

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


static void dump_counted_volumes_map(const std::string name, const CountedVolumesMap& map) {
  cout << name << ": ";
  if (map.empty()) {
    cout << "empty\n";
  }
  else {
    cout << "\n";
  }

  for (auto it: map) {
    cout << "  " << it.first << " -> {";
    it.second.dump();
    cout << "}\n";
  }
}

void Partition::dump() {
  GeometryObject::dump_array(*this, geometry_objects);

  // dump
  dump_counted_volumes_map(
      "directly_contained_counted_volume_objects",
      directly_contained_counted_volume_objects
  );
  dump_counted_volumes_map(
      "enclosing_counted_volume_objects",
      enclosing_counted_volume_objects
  );

  Region::dump_array(regions);

  for (size_t i = 0; i < walls_per_subpart.size(); i++) {
    if (!walls_per_subpart[i].empty()) {
      Vec3 llf, urb;
      get_subpart_llf_point(i, llf);
      urb = llf + Vec3(config.subpartition_edge_length);

      cout << "subpart: " << i << ", llf: " << llf << ", urb: " << urb << "\n";
      walls_per_subpart[i].dump("Indices contained in a subpartition");
    }
  }
}


geometry_object_id_t Partition::find_smallest_counted_volume_recursively(const GeometryObject& obj, const Vec3& pos) {

  assert(obj.is_counted_volume);
  auto it = directly_contained_counted_volume_objects.find(obj.id);
  if (it == directly_contained_counted_volume_objects.end()) {
    // there are no contained objects
    return obj.id;
  }

  const uint_set<geometry_object_id_t>& subobject_ids = it->second;

  for (geometry_object_id_t subobj_id: subobject_ids) {
    assert(subobj_id != obj.id);
    GeometryObject& subobj = get_geometry_object(subobj_id);
    if (CountedVolumesUtil::is_point_inside_counted_volume(subobj, pos)) {
      // the hierarchy must be a tree, so we directly find the path to the leaf of the object 'containment' tree
      return find_smallest_counted_volume_recursively(subobj, pos);
    }
  }

  // none of the contained objects contains point at 'pos'
  return obj.id;
}


geometry_object_id_t Partition::determine_counted_volume_id(const Vec3& pos) {
  // find the first object we are in
  for (GeometryObject& obj: geometry_objects) {
    if (obj.is_counted_volume && CountedVolumesUtil::is_point_inside_counted_volume(obj, pos)) {

      // follow the containment graph to determine the smallest counted volume
      geometry_object_id_t res = find_smallest_counted_volume_recursively(obj, pos);

#ifdef DEBUG_COUNTED_VOLUMES
      cout <<
          "Counted volumes: assigned obj " << get_geometry_object(res).name << " (id:" << res << ")" <<
          " for pos " << pos << "\n";
#endif

			return res;
    }
  }

#ifdef DEBUG_COUNTED_VOLUMES
      cout << "Counted volumes: assigned 'outside' for pos " << pos << "\n";
#endif

  // nothing found
  return COUNTED_VOLUME_ID_OUTSIDE_ALL;
}


void Partition::to_data_model(Json::Value& mcell) const {

  Json::Value& geometrical_objects = mcell[KEY_GEOMETRICAL_OBJECTS];
  Json::Value& object_list = geometrical_objects[KEY_OBJECT_LIST];

  for (const GeometryObject& g: geometry_objects) {
    Json::Value object;
    g.to_data_model(object, *this, config);
    object_list.append(object);
  }
}

} // namespace mcell
