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

#include "bng/bng.h"

#include "rng.h" // MCell 3
#include "isaac64.h"
#include "mcell_structs.h"
#include "logging.h"
#include "wall_util.h"

#include "partition.h"
#include "geometry.h"
#include "release_event.h"
#include "datamodel_defines.h"

#include "geometry_utils.h"
#include "geometry_utils.inc" // uses get_wall_bounding_box, maybe not include this file
#include "collision_utils.inc"
#include "dump_state.h"

using namespace std;

namespace MCell {

// TODO: replace int with wall_index_t where possible

/***************************************************************************
compatible_edges:
  In: array of pointers to walls
      index of first wall
      index of edge in first wall
      index of second wall
      index of edge in second wall
  Out: 1 if the edge joins the two walls
       0 if not (i.e. the wall doesn't contain the edge or the edge is
       traversed in the same direction in each or the two walls are
       actually the same wall)
***************************************************************************/
static int compatible_edges(
    const Partition& p,
    const vector<wall_index_t>& faces, int wA, int eA, int wB, int eB) {

  const Vec3 *vA0, *vA1, *vA2, *vB0, *vB1, *vB2;

  const Wall& wall_a = p.get_wall(faces[wA]);
  const Wall& wall_b = p.get_wall(faces[wB]);

  if ((wA < 0) || (eA < 0) || (wB < 0) || (eB < 0))
    return 0;

  vA0 = &p.get_wall_vertex(wall_a, eA);
  if (eA == 2) {
    vA1 = &p.get_wall_vertex(wall_a, 0);
  }
  else {
    vA1 = &p.get_wall_vertex(wall_a, eA + 1);
  }

  if (eA == 0) {
    vA2 = &p.get_wall_vertex(wall_a, 2);
  }
  else {
    vA2 = &p.get_wall_vertex(wall_a, eA - 1);
  }

  vB0 = &p.get_wall_vertex(wall_b, eB);
  if (eB == 2) {
    vB1 = &p.get_wall_vertex(wall_b, 0);
  }
  else {
    vB1 = &p.get_wall_vertex(wall_b, eB + 1);
  }

  if (eB == 0) {
    vB2 = &p.get_wall_vertex(wall_b, 2);
  }
  else {
    vB2 = &p.get_wall_vertex(wall_b, eB - 1);
  }

  return ((vA0 == vB1 && vA1 == vB0 && vA2 != vB2) ||
          (vA0->x == vB1->x && vA0->y == vB1->y && vA0->z == vB1->z &&
           vA1->x == vB0->x && vA1->y == vB0->y && vA1->z == vB0->z &&
           !(vA2->x == vB2->x && vA2->y == vB2->y && vA2->z == vB2->z)));
}


/*****************************************************************************
 have_common_region checks if wall1 and wall2 located on the (same) object
 are part of a common region or not
******************************************************************************/
static bool have_common_region(
    const Partition& p, const GeometryObject& obj, int wall1, int wall2) {

  const Wall& w1 = p.get_wall(obj.wall_indices[wall1]);
  const Wall& w2 = p.get_wall(obj.wall_indices[wall2]);

  for (region_index_t index: w1.regions) {
    if (w2.regions.count(index) == 1) {
      return true;
    }
  }
  return false;
}


/***************************************************************************
refine_edge_pairs:
  In: the head of a linked list of shared edges
      array of pointers to walls
  Out: No return value.  The best-matching pair of edges percolates up
       to be first in the list.  "Best-matching" means that the edge
       is traversed in different directions by each face, and that the
       normals of the two faces are as divergent as possible.
***************************************************************************/
static void refine_edge_pairs(
    const Partition& p, const GeometryObject& obj, poly_edge *pe, const vector<wall_index_t>& faces) {

#define TSWAP(x, y) temp = (x); (x) = (y); (y) = temp

  int temp;

  float_t best_align = 2;
  bool share_region = false;
  poly_edge* best_p1 = pe;
  poly_edge* best_p2 = pe;
  int best_n1 = 1;
  int best_n2 = 2;

  poly_edge* p1 = pe;
  int n1 = 1;
  while (p1 != NULL && p1->n >= n1) {
    int wA, eA;
    if (n1 == 1) {
      wA = p1->face[0];
      eA = p1->edge[0];
    } else {
      wA = p1->face[1];
      eA = p1->edge[1];
    }

    poly_edge* p2;
    int n2;
    if (n1 == 1) {
      n2 = n1 + 1;
      p2 = p1;
    } else {
      n2 = 1;
      p2 = p1->next;
    }
    while (p2 != NULL && p2->n >= n2) {
      int wB, eB;
      if (n2 == 1) {
        wB = p2->face[0];
        eB = p2->edge[0];
      } else {
        wB = p2->face[1];
        eB = p2->edge[1];
      }

      // as soon as we hit an incompatible edge we can break out of the p2 loop
      // and continue scanning the next p1
      if (compatible_edges(p, faces, wA, eA, wB, eB)) {
        const Wall& wall_a = p.get_wall(wA);
        const Wall& wall_b = p.get_wall(wB);
        assert(wall_a.wall_constants_precomputed);
        assert(wall_b.wall_constants_precomputed);

        float_t align = dot(wall_a.normal, wall_b.normal);

        // as soon as two walls have a common region we only consider walls who
        // share (any) region. We need to reset the best_align to make sure we
        // don't pick any wall that don't share a region discovered previously
        bool common_region = have_common_region(p, obj, wA, wB);
        if (common_region) {
          if (!share_region) {
            best_align = 2;
          }
          share_region = true;
        }

        if (common_region || !share_region) {
          if (align < best_align) {
            best_p1 = p1;
            best_p2 = p2;
            best_n1 = n1;
            best_n2 = n2;
            best_align = align;
          }
        }
      } else {
        break;
      }

      if (n2 == 1)
        n2++;
      else {
        p2 = p2->next;
        n2 = 1;
      }
    }

    if (n1 == 1)
      n1++;
    else {
      p1 = p1->next;
      n1 = 1;
    }
  }

  /* swap best match into top spot */
  if (best_align > 1.0)
    return; /* No good pairs. */

  TSWAP(best_p1->face[best_n1-1], pe->face[0]);
  TSWAP(best_p1->edge[best_n1-1], pe->edge[0]);
  TSWAP(best_p2->face[best_n2-1], pe->face[1]);
  TSWAP(best_p2->edge[best_n2-1], pe->edge[1]);

#undef TSWAP
}

/***************************************************************************
surface_net:
  In: array of pointers to walls
      integer length of array
  Out: -1 if the surface is a manifold, 0 if it is not, 1 on malloc failure
       Walls end up connected across their edges.
  Note: Two edges must have their vertices listed in opposite order (i.e.
        connect two faces pointing the same way) to be linked.  If more than
        two faces share the same edge and can be linked, the faces with
        normals closest to each other will be linked.  We do not assume that
        the object is connected.  All pieces must be a manifold, however,
        for the entire object to be a manifold.  (That is, there must not
        be any free edges anywhere.)  It is possible to build weird, twisty
        self-intersecting things.  The behavior of these things during a
        simulation is not guaranteed to be well-defined.
***************************************************************************/
static int surface_net(Partition& p, GeometryObject& obj) {

  uint nfaces = obj.wall_indices.size();
  vector<wall_index_t> facelist = obj.wall_indices;

  Edge *e;
  int is_closed = 1;

  struct edge_hashtable eht;
  int nkeys = (3 * nfaces) / 2;
  if (ehtable_init(&eht, nkeys))
    return 1;

  for (uint face_index = 0; face_index < nfaces; face_index++) {

    Wall& w = p.get_wall(facelist[face_index]);

    int k;
    int nedge = 3;

    for (int j = 0; j < nedge; j++) {

      if (j + 1 < nedge)
        k = j + 1;
      else
        k = 0;

      poly_edge pe;
      const Vec3& vert_j = p.get_wall_vertex(w, j);
      pe.v1x = vert_j.x;
      pe.v1y = vert_j.y;
      pe.v1z = vert_j.z;
      const Vec3& vert_k = p.get_wall_vertex(w, k);
      pe.v2x = vert_k.x;
      pe.v2y = vert_k.y;
      pe.v2z = vert_k.z;
      pe.face[0] = face_index;
      pe.edge[0] = j;

      if (ehtable_add(&eht, &pe))
        return 1;
    }
  }

  for (int i = 0; i < nkeys; i++) {
    poly_edge *pep = (eht.data + i);
#ifdef DEBUG_GEOM_OBJ_INITIALIZATION
    dump_poly_edge(i, pep);
#endif
    while (pep != NULL) {
      if (pep->n > 2) {
        refine_edge_pairs(p, obj, pep, facelist);
      }
      if (pep->n >= 2) {
        if (pep->face[0] != -1 && pep->face[1] != -1) {
          if (compatible_edges(p, facelist, pep->face[0], pep->edge[0], pep->face[1], pep->edge[1])) {

            Wall& face0 = p.get_wall(facelist[pep->face[0]]);
            Wall& face1 = p.get_wall(facelist[pep->face[1]]);

            assert(pep->face[0] < (int)facelist.size());
            assert(pep->face[1] < (int)facelist.size());
            face0.nb_walls[pep->edge[0]] = facelist[pep->face[1]];
            face1.nb_walls[pep->edge[1]] = facelist[pep->face[0]];

            Edge e;
            e.forward_index = facelist[pep->face[0]];
            e.backward_index = facelist[pep->face[1]];
            e.edge_num_used_for_init = pep->edge[0];

            assert(face0.wall_constants_precomputed);
            assert(face1.wall_constants_precomputed);
            e.reinit_edge_constants(p);

            face0.edges[pep->edge[0]] = e;
            face1.edges[pep->edge[1]] = e;
          }

        }
        else {
          is_closed = 0;
        }
      } else if (pep->n == 1) {
        is_closed = 0;
        Edge e;

        e.forward_index = facelist[pep->face[0]];
        e.backward_index = WALL_INDEX_INVALID;
        /* Don't call reinit_edge_constants unless both edges are set */
      }
      pep = pep->next;
    }
  }

  ehtable_kill(&eht);
  return -is_closed; /* We use 1 to indicate malloc failure so return 0/-1 */
}


void GeometryObject::initialize_neighboring_walls_and_their_edges(Partition& p) {
  surface_net(p, *this);
}


void GeometryObject::dump(const Partition& p, const std::string ind) const {
  cout << ind <<
      "GeometryObject: id:" << id << ", name:" << name <<
      ", is_counted_volume " << is_counted_volume << "\n";

  for (wall_index_t i: wall_indices) {
    cout << ind << "  " << i << ": \n";
    p.get_wall(i).dump(p, ind + "    ");
  }
}


void GeometryObject::dump_array(const Partition& p, const std::vector<GeometryObject>& vec) {
  cout << "GeometryObject array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump(p, "  ");
  }
}

void GeometryObject::to_data_model(
    const Partition& p, const SimulationConfig& config, Json::Value& object) const {

  vector<uint> vertex_order;

  object[KEY_NAME] = DMUtil::remove_obj_name_prefix(parent_name, name);
  object[KEY_MATERIAL_NAMES].append(Json::Value(KEY_VALUE_MEMBRANE));


  bool first = true; // to indicate when to use a comma

  // first note which vertices this object uses
  uint_set<vertex_index_t> used_vertex_indices;
  for (wall_index_t wall_index: wall_indices) {
    const Wall& w = p.get_wall(wall_index);

    for (vertex_index_t vertex_index: w.vertex_indices) {
      used_vertex_indices.insert(vertex_index);
    }
  }

  // then generate vertices and remember mapping
  Json::Value& vertex_list = object[KEY_VERTEX_LIST];
  map<vertex_index_t, uint> map_vertex_index_to_vertex_list_index;
  uint current_index_in_vertex_list = 0;
  for (vertex_index_t i = 0; i < p.get_geometry_vertex_count(); i++) {

    if (used_vertex_indices.count(i) == 1) {
      // define mapping vertex_index -> index in vertex array
      map_vertex_index_to_vertex_list_index[i] = current_index_in_vertex_list;
      current_index_in_vertex_list++;

      // append triple x, y, z
      Vec3 pos = p.get_geometry_vertex(i) * Vec3(p.config.length_unit);
      DMUtil::append_triplet(vertex_list, pos.x, pos.y, pos.z);
    }
  }

  // element connections - they correspond to the ordering of walls
  Json::Value& element_connections = object[KEY_ELEMENT_CONNECTIONS];
  for (wall_index_t wall_index: wall_indices) {
    const Wall& w = p.get_wall(wall_index);

    Json::Value vertex_indices;
    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      assert(map_vertex_index_to_vertex_list_index.count(w.vertex_indices[i]) == 1);
      vertex_indices.append(map_vertex_index_to_vertex_list_index[w.vertex_indices[i]]);
    }

    element_connections.append(vertex_indices);
  }

  // surface regions
  insertion_ordered_set<region_index_t> region_indices;
  map<wall_index_t, uint> map_wall_index_to_order_index_in_object;
  for (size_t i = 0; i < wall_indices.size(); i++) {

    map_wall_index_to_order_index_in_object[wall_indices[i]] = i;

    const Wall& w = p.get_wall(wall_indices[i]);
    for (region_index_t region_index: w.regions) {
      region_indices.insert_ordered(region_index);
    }
  }

  vector<Json::Value> surface_regions_vec;
  for (region_index_t region_index: region_indices.get_as_vector()) {
    const Region& reg = p.get_region(region_index);

    if (reg.name_has_all_suffix()) {
      continue;
    }

    Json::Value surface_region;

    string name = DMUtil::get_surface_region_name(reg.name);
    surface_region[KEY_NAME] = name;

    Json::Value include_elements;

    for (auto it: reg.walls_and_edges) {
       include_elements.append(map_wall_index_to_order_index_in_object[it.first]);
    }
    surface_region[KEY_INCLUDE_ELEMENTS] = include_elements;
    surface_regions_vec.push_back(surface_region);
  }

  // we may create the define_surface_regions only when there is at least one region
  if (!surface_regions_vec.empty()) {
    Json::Value& define_surface_regions = object[KEY_DEFINE_SURFACE_REGIONS];

    for (auto& surface_region: surface_regions_vec) {
      define_surface_regions.append(surface_region);
    }
  }
}


void InitialRegionMolecules::to_data_model(
    const BNG::SpeciesContainer& all_species,
    Json::Value& initial_region_molecules
) const {
  initial_region_molecules[KEY_MOLECULE] = all_species.get(species_id).name;
  initial_region_molecules[KEY_ORIENT] = DMUtil::orientation_to_str(orientation);
  if (const_num_not_density) {
    initial_region_molecules[KEY_MOLECULE_NUMBER] = to_string(release_num);
  }
  else {
    initial_region_molecules[KEY_MOLECULE_DENSITY] = DMUtil::f_to_string(release_density);
  }
}


void Region::init_from_whole_geom_object(const GeometryObject& obj) {
  name = obj.name + REGION_ALL_SUFFIX_W_COMMA;
  geometry_object_id = obj.id;

  // simply add all walls
  for (const wall_index_t& wi: obj.wall_indices) {
    add_wall_to_walls_and_edges(wi, false);
  }
}


void Region::init_surface_region_edges(const Partition& p) {
  assert(!walls_and_edges.empty());

  // must be run after edge must be initializaion,
  // however, not all edges need to be initialized
  for (auto& it: walls_and_edges) {
    const Wall& w = p.get_wall(it.first);
    for (edge_index_t ei = 0; ei < EDGES_IN_TRIANGLE; ei++) {
      const Edge& edge = w.edges[ei];
      // uninitialized edges are considered to be bordering
      if (!edge.is_initialized()) {
        it.second.insert(ei);
      }
      else {
        // is one of the forward or backward walls in our region?
        // if not, then the edge is on the border of this region
        if (walls_and_edges.count(edge.forward_index) == 0 ||
            walls_and_edges.count(edge.backward_index) == 0) {
          it.second.insert(ei);
        }
      }
    }
  }
}


/* tetrahedralVol returns the (signed) volume of the tetrahedron spanned by
 * the vertices a, b, c, and d.
 * The formula was taken from "Computational Geometry" (2nd Ed) by J. O'Rourke
 */
static float_t tetrahedral_volume(
    const Vec3& a, const Vec3& b,
    const Vec3& c, const Vec3& d
) {
  // TODO: maybe use vector operators
  return 1.0 / 6.0 *
      (-1 * (a.z - d.z) * (b.y - d.y) * (c.x - d.x) +
      (a.y - d.y) * (b.z - d.z) * (c.x - d.x) +
      (a.z - d.z) * (b.x - d.x) * (c.y - d.y) -
      (a.x - d.x) * (b.z - d.z) * (c.y - d.y) -
      (a.y - d.y) * (b.x - d.x) * (c.z - d.z) +
      (a.x - d.x) * (b.y - d.y) * (c.z - d.z));
}


/***************************************************************************
  Based on is_manifold:
  In: r: A region. This region must already be painted on walls. The edges must
         have already been added to the object (i.e. sharpened).
      count_regions_flag: This is usually set, unless we are only checking
                          volumes for dynamic geometries.
  Out: 1 if the region is a manifold, 0 otherwise.
  Note: by "manifold" we mean "orientable compact two-dimensional
        manifold without boundaries embedded in R3"
***************************************************************************/
// TODO: call with initialization?
void Region::initialize_volume_info_if_needed(const Partition& p) {
  if (volume_info_initialized) {
    return;
  }

  if (walls_and_edges.empty()) {
    mcell_internal_error("Region '%s' has NULL wall array!", name.c_str());
  }

  // initialize bounding box
  Geometry::compute_region_bounding_box(p, *this, bounding_box_llf, bounding_box_urb);

  // use the center of the region bounding box as reference point for
  // computing the volume
  Vec3 d = (bounding_box_llf + bounding_box_urb) * Vec3(0.5);

  perf() << "Computing volume of region " << name << "\n";

  volume = 0;
  float_t current_volume = 0;
  for (auto it: walls_and_edges) {

    if (!it.second.empty()) {
      // does this wall represent a region border?
      region_is_manifold = false;
      volume_info_initialized = true;
    }

    const Wall& w = p.get_wall(it.first);

    // compute volume of tetrahedron with w as its face
    current_volume += tetrahedral_volume(
        p.get_wall_vertex(w, 0), p.get_wall_vertex(w, 1), p.get_wall_vertex(w, 2), d);
  }

  perf() << "  - volume computed\n";

  volume = current_volume;
  region_is_manifold = true;
  volume_info_initialized = true;
}


void Region::initialize_wall_subpart_mapping_if_needed(const Partition& p) {
  if (walls_per_subpart_initialized) {
    return;
  }

  perf() << "Preparing wall subpart mapping for region " << name << "\n";

  // first initialize region wall subpart mapping
  walls_per_subpart.clear();

  // FIXME: update when a wall of this region is moved
  for (const auto& wall_and_edges_pair: walls_and_edges) {
    const Wall& w = p.get_wall(wall_and_edges_pair.first);

    for (subpart_index_t si: w.present_in_subparts) {
      walls_per_subpart[si].push_back(w.index);
    }
  }

  perf() << "  - wall subpart mapping prepared\n";

  walls_per_subpart_initialized = true;
}

// returns true if the waypoint with current_waypoint_index is in this region
// similar code as in Partition::create_waypoint
bool Region::initialize_region_waypoint(
    const Partition& p,
    const IVec3& current_waypoint_index,
    const bool use_previous_waypoint,
    const IVec3& previous_waypoint_index,
    const bool previous_waypoint_present_in_region
) {

  // waypoint is always in the center of a subpartition
  // thanks to the extra margin on each side, we might be out of the partition
  if (p.is_valid_waypoint_index(current_waypoint_index)) {
    const Waypoint& waypoint = p.get_waypoint(current_waypoint_index);

    if (use_previous_waypoint) {
      const Waypoint& previous_waypoint = p.get_waypoint(previous_waypoint_index);

      uint num_crossed = CollisionUtil::get_num_crossed_region_walls(
          p, waypoint.pos, previous_waypoint.pos,
          *this
      );

      if ((num_crossed % 2 == 0 && previous_waypoint_present_in_region) ||
          (num_crossed % 2 == 1 && !previous_waypoint_present_in_region)) {

        // still in the same region
        waypoints_in_this_region.insert(current_waypoint_index);
        return true;
      }
      else {
        return false;
      }

    }

    if (CollisionUtil::is_point_inside_region_no_waypoints(p, waypoint.pos, *this)) {
      waypoints_in_this_region.insert(current_waypoint_index);
      return true;
    }
  }

  return false;
}

void Region::initialize_region_waypoints_if_needed(const Partition& p) {

  if (region_waypoints_initialized) {
    return;
  }

  perf() << "Initializing waypoints for region " << name << "\n";

  assert(walls_per_subpart_initialized);

  // get initial waypoint
  waypoints_in_this_region.clear();

  IVec3 llf_waypoint_index;
  subpart_index_t subpart_index = p.get_subpart_index(bounding_box_llf);
  p.get_subpart_3d_indices_from_index(subpart_index, llf_waypoint_index);

  // then compute how many waypoints in each dimension we need to check
  Vec3 region_dims = bounding_box_urb - bounding_box_llf;
  // num_waypoints need to be incremented by 2 - for each of the side regions
  IVec3 num_waypoints = region_dims / Vec3(p.config.subpartition_edge_length) + Vec3(2);

  bool previous_waypoint_present = false;
  bool use_previous_waypoint = false;
  IVec3 previous_waypoint_index;

  for (int x = 0; x < num_waypoints.x; x++) {
    for (int y = 0; y < num_waypoints.y; y++) {
      for (int z = 0; z < num_waypoints.z; z++) {

        IVec3 current_waypoint_index( llf_waypoint_index + IVec3(x, y, z) );

        previous_waypoint_present = initialize_region_waypoint(
            p, current_waypoint_index,
            use_previous_waypoint, previous_waypoint_index, previous_waypoint_present
        );

        use_previous_waypoint = true;
        previous_waypoint_index = current_waypoint_index;
      }

      // starting a new line - start from scratch
      // NOTE: maybe we do not need to restart this every line
      use_previous_waypoint = false;
      previous_waypoint_present = false;
    }
  }

  region_waypoints_initialized = true;

  perf() << "  - region waypoints initialized\n";
}



bool Region::is_point_inside(const Partition& p, const Vec3& pos) {

  initialize_volume_info_if_needed(p);
  initialize_wall_subpart_mapping_if_needed(p);
  initialize_region_waypoints_if_needed(p);

  // as a fallback for debugging one can call directly
  //return CollisionUtil::is_point_inside_region_no_waypoints(p, pos, *this);

  assert(is_manifold() && "Not sure what to do here yet");


  // get a waypoint close to this position
  IVec3 waypoint_index;
  subpart_index_t subpart_index = p.get_subpart_index(pos); // TODO: get indices directly
  p.get_subpart_3d_indices_from_index(subpart_index, waypoint_index);
  const Waypoint& waypoint = p.get_waypoint(waypoint_index);

  map<geometry_object_index_t, uint> num_crossed_walls_per_object;
  uint num_crossed = CollisionUtil::get_num_crossed_region_walls(
      p, pos, waypoint.pos,
      *this
  );

  bool waypoint_is_inside_this_region = waypoints_in_this_region.count(waypoint_index) != 0;
  if (waypoint_is_inside_this_region) {
    return num_crossed % 2 == 0;
  }
  else {
    return num_crossed % 2 == 1;
  }
}


void Region::dump(const std::string ind) const {
  cout << ind << "Region : " <<
      "name:" << name <<
      ", species_id: " << ((species_id == SPECIES_ID_INVALID) ? string("invalid") : to_string(species_id)) <<
      "\n";

  for (auto& wall_it: walls_and_edges) {
    cout << ind << "  " << "wall " << wall_it.first << ", region edges: {";
    for (auto& edge: wall_it.second) {
      cout << edge << ", ";
    }
    cout << "}\n";
  }
}


void Region::dump_array(const std::vector<Region>& vec) {
  cout << "Region array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump("  ");
  }
}


void Region::to_data_model(const Partition& p, Json::Value& modify_surface_region) const {

  if (initial_region_molecules.empty()) {
    DMUtil::add_version(modify_surface_region, VER_DM_2018_01_11_1330);
  }
  else {
    DMUtil::add_version(modify_surface_region, VER_DM_2020_07_12_1600);
  }

  modify_surface_region[KEY_DESCRIPTION] = "";
  modify_surface_region[KEY_OBJECT_NAME] =
      DMUtil::remove_obj_name_prefix(p.get_geometry_object(geometry_object_id).name);

  string region_name = DMUtil::get_region_name(name);
  if (region_name == "ALL") {
    modify_surface_region[KEY_REGION_SELECTION] = VALUE_ALL;
    modify_surface_region[KEY_REGION_NAME] = "";
  }
  else {
    modify_surface_region[KEY_REGION_SELECTION] = VALUE_SEL;
    modify_surface_region[KEY_REGION_NAME] = region_name;
  }

  modify_surface_region[KEY_NAME] = ""; // don't care

  if (species_id != SPECIES_ID_INVALID) {
    modify_surface_region[KEY_SURF_CLASS_NAME] = p.get_all_species().get(species_id).name;
  }

  if (!initial_region_molecules.empty()) {
    Json::Value& initial_region_molecules_list = modify_surface_region[KEY_INITIAL_REGION_MOLECULES_LIST];
    for (const InitialRegionMolecules& info: initial_region_molecules) {
      Json::Value initial_region_molecules_item;
      info.to_data_model(p.get_all_species(), initial_region_molecules_item);
      initial_region_molecules_list.append(initial_region_molecules_item);
    }
  }
}



// may be also used for reinitialization
void Grid::initialize(const Partition& p, const Wall& w) {

  if (is_initialized()) {
    // keep the same number of items
    #ifndef NDEBUG
      small_vector<molecule_id_t> molecule_ids;
      get_contained_molecules(molecule_ids);
      assert(molecule_ids.size() == num_occupied);
    #endif
  }
  else {
    num_occupied = 0;
  }

  num_tiles_along_axis = (int)ceil_f(sqrt_f(w.area));
  if (num_tiles_along_axis < 1) {
    num_tiles_along_axis = 1;
  }

  num_tiles = num_tiles_along_axis * num_tiles_along_axis;

  molecules_per_tile.resize(num_tiles, MOLECULE_ID_INVALID);

  strip_width_rcp = 1.0 / (w.uv_vert2.v / ((float_t)num_tiles_along_axis));
  vert2_slope = w.uv_vert2.u / w.uv_vert2.v;
  fullslope = w.uv_vert1_u / w.uv_vert2.v;

  binding_factor = ((float_t)num_tiles) / w.area;

  const Vec3& vert0_tmp = p.get_wall_vertex(w, 0);

  vert0.u = dot(vert0_tmp, w.unit_u);
  vert0.v = dot(vert0_tmp, w.unit_v);
}


// populates array molecules with ids of molecules belonging to this grid
void Grid::get_contained_molecules(
    small_vector<molecule_id_t>& molecule_ids
) const {
  // might be optimized by retaining a vector/set that gets invalidated if anything changes
  molecule_ids.clear();
  for (molecule_id_t id: molecules_per_tile) {
    if (id != MOLECULE_ID_INVALID) {
      molecule_ids.push_back(id);
    }
  }
}


void Grid::dump() const {
  // dumping just occupied locations and base info for now
  cout << "Grid: num_tiles: " << num_tiles << ", num_occupied: " << num_occupied << "\n";
  for (uint i = 0; i < molecules_per_tile.size(); i++) {
    molecule_id_t id = molecules_per_tile[i];
    if (id != MOLECULE_ID_INVALID) {
      cout << "[" << i << "]" << id << "\n";
    }
  }
}


/***************************************************************************
init_edge_transform
  In: e: pointer to an edge
      edgenum: integer telling which edge (0-2) of the "forward" face we are
  Out: No return value.  Coordinate transform in edge struct is set.
  Note: Don't call this on a non-shared edge.
***************************************************************************/
void Edge::reinit_edge_constants(const Partition& p) {
  assert(is_shared_edge());

  if (!is_initialized()) {
    return;
  }

  Wall wf = p.get_wall(forward_index);
  Wall wb = p.get_wall(backward_index);

#ifdef DEBUG_EDGE_INITIALIZATION
  std::cout << "Edge initialization, edgenum: " << edge_num_used_for_init << "\n";
  wf.dump(p, "", true);
  wb.dump(p, "", true);
#endif

  edge_index_t i = edge_num_used_for_init;
  assert(i < VERTICES_IN_TRIANGLE);
  edge_index_t j = i + 1;
  if (j == VERTICES_IN_TRIANGLE)
    j = 0;

  const Vec3& wf_vert_0 = p.get_geometry_vertex(wf.vertex_indices[0]);
  const Vec3& wf_vert_i = p.get_geometry_vertex(wf.vertex_indices[i]);
  const Vec3& wf_vert_j = p.get_geometry_vertex(wf.vertex_indices[j]);

  const Vec3& wb_vert_0 = p.get_geometry_vertex(wb.vertex_indices[0]);


  /* Intermediate basis from the perspective of the forward frame */
  Vec3 diff_i_0 = wf_vert_i - wf_vert_0;

  Vec2 O_f;
  O_f.u = dot(diff_i_0, wf.unit_u); // should be called after wall init
  O_f.v = dot(diff_i_0, wf.unit_v); /* Origin */

  Vec3 diff_j_0 = wf_vert_j - wf_vert_0;
  Vec2 temp_ff;
  temp_ff.u = dot(diff_j_0, wf.unit_u) - O_f.u;
  temp_ff.v = dot(diff_j_0, wf.unit_v) - O_f.v; /* Far side of e */

  assert(!cmp_eq(len2_squared(temp_ff), 0.0, EPS));
  float_t d_f = 1.0 / len2(temp_ff);

  Vec2 ehat_f, fhat_f;
  ehat_f = temp_ff * d_f; /* ehat along edge */
  fhat_f.u = -ehat_f.v;
  fhat_f.v = ehat_f.u; /* fhat 90 degrees CCW */

  /* Intermediate basis from the perspective of the backward frame */
  Vec3 diff_i_b0 = wf_vert_i - wb_vert_0;
  Vec2 O_b;
  O_b.u = dot(diff_i_b0, wb.unit_u);
  O_b.v = dot(diff_i_b0, wb.unit_v); /* Origin */

  Vec3 diff_j_b0 = wf_vert_j - wb_vert_0;
  Vec2 temp_fb;
  temp_fb.u = dot(diff_j_b0, wb.unit_u) - O_b.u;
  temp_fb.v = dot(diff_j_b0, wb.unit_v) - O_b.v; /* Far side of e */

  assert(!cmp_eq(len2_squared(temp_fb), 0.0, EPS));
  float_t d_b = 1.0 / len2(temp_fb);

  Vec2 ehat_b, fhat_b;

  ehat_b = temp_fb * d_b; /* ehat along edge */
  fhat_b.u = -ehat_b.v;
  fhat_b.v = ehat_b.u; /* fhat 90 degrees CCW */

  /* Calculate transformation matrix */

  float_t mtx[2][2];
  mtx[0][0] = ehat_f.u * ehat_b.u + fhat_f.u * fhat_b.u;
  mtx[0][1] = ehat_f.v * ehat_b.u + fhat_f.v * fhat_b.u;
  mtx[1][0] = ehat_f.u * ehat_b.v + fhat_f.u * fhat_b.v;
  mtx[1][1] = ehat_f.v * ehat_b.v + fhat_f.v * fhat_b.v;

  /* Calculate translation vector */

  Vec2 q;
  q = O_b;

  q.u -= mtx[0][0] * O_f.u + mtx[0][1] * O_f.v;
  q.v -= mtx[1][0] * O_f.u + mtx[1][1] * O_f.v;

  /* Store the results */
  cos_theta = mtx[0][0];
  sin_theta = mtx[0][1];
  translate = q;

#ifdef DEBUG_EDGE_INITIALIZATION
  dump();
#endif
}

void Edge::debug_check_values_are_uptodate(const Partition& p) {
  if (!is_initialized()) {
    return;
  }
#ifndef NDEBUG
  Vec2 orig_translate = translate;
#endif
  float_t orig_cos_theta = cos_theta;
  float_t orig_sin_theta = sin_theta;
#ifdef DEBUG_EDGE_INITIALIZATION
  dump();
#endif
  reinit_edge_constants(p);
#ifdef DEBUG_EDGE_INITIALIZATION
  dump();
#endif
  assert(cmp_eq(orig_translate, translate));
  assert(cmp_eq(orig_cos_theta, cos_theta));
  assert(cmp_eq(orig_sin_theta, sin_theta));
}


void Edge::dump(const std::string ind) const {
  Vec2 translate_rounded;
  translate_rounded.x = (translate.x < EPS) ? 0.0 : translate.x;
  translate_rounded.y = (translate.y < EPS) ? 0.0 : translate.y;

  cout << ind <<
      "Edge: "
      "forward_index: " << forward_index <<
      ", backward_index: " << backward_index <<
      ", translate: " << translate_rounded <<
      ", cos_theta: " << ((cos_theta < EPS) ? 0.0 : cos_theta) <<
      ", sin_theta: " << ((sin_theta < EPS) ? 0.0 : sin_theta) <<
      ", edge_num_used_for_init: " << edge_num_used_for_init << "\n";
}


void Wall::precompute_wall_constants(const Partition& p) {

  const Vec3& v0 = p.get_geometry_vertex(vertex_indices[0]);
  const Vec3& v1 = p.get_geometry_vertex(vertex_indices[1]);
  const Vec3& v2 = p.get_geometry_vertex(vertex_indices[2]);

  Vec3 vA, vB, vX;
  vA = v1 - v0;
  vB = v2 - v0;
  vX = cross(vA, vB);

  area = 0.5 * len3(vX);

  if (!distinguishable_f(area, 0, EPS)) {
    /* this is a degenerate polygon.
    * perform initialization and quit. */
    normal = Vec3(0);
    unit_u = Vec3(0);
    unit_v = Vec3(0);
    uv_vert1_u = 0;
    uv_vert2 = Vec2(0);
    distance_to_origin = 0;
    return;
  }

  Vec3 f1 = v1 - v0;
  float_t f1_len_squared = len3_squared(f1);
  assert(f1_len_squared != 0);
  float_t inv_f1_len = 1 / sqrt_f(f1_len_squared);

  unit_u = f1 * Vec3(inv_f1_len);

  Vec3 f2 = v2 - v0;
  normal = cross(unit_u, f2);

  float_t norm_len_squared = len3_squared(normal);
  assert(norm_len_squared != 0);
  float_t inv_norm_len = 1 / sqrt_f(norm_len_squared);
  normal = normal * Vec3(inv_norm_len);

  unit_v = cross(normal, unit_u);

  distance_to_origin =  dot(v0, normal);

  uv_vert1_u = dot(f1, unit_u);
  uv_vert2.u = dot(f2, unit_u);
  uv_vert2.v = dot(f2, unit_v);

  wall_constants_precomputed = true;
}

void Wall::reinit_edge_constants(const Partition& p) {
  // edges, uses info set by init_tri_wall
  for (uint edge_index = 0; edge_index < EDGES_IN_TRIANGLE; edge_index++) {
    if (edges[edge_index].is_shared_edge()) {
      edges[edge_index].reinit_edge_constants(p);
    }
  }
}


void Wall::dump(const Partition& p, const std::string ind, const bool for_diff) const {
  if (for_diff) {
    cout << "wall[side: " << side << "]: ";
    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      vertex_index_t vertex_index = vertex_indices[i];
      Vec3 pos = p.get_geometry_vertex(vertex_index);
      cout << pos;
      if (i != VERTICES_IN_TRIANGLE - 1) {
        cout << ", ";
      }
    }
    cout << "\n";
  }
  else {
    cout << ind <<
        "id: " << id << ", index: " << index << ", side: " << side <<
        ", object_id: " << object_id << ", object_index: " << object_index << "\n";

    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      vertex_index_t vertex_index = vertex_indices[i];
      Vec3 pos = p.get_geometry_vertex(vertex_index);
      cout << ind << "vertex_index: " << vertex_index << ":" << pos << "\n";
    }

    cout << ind << "edges:\n";
    for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
      cout << ind << i << ":\n";
      edges[i].dump(ind + "  ");
    }

    cout << ind;
    for (uint i = 0; i < EDGES_IN_TRIANGLE; i++) {
      cout << "nb_walls[" << i << "]: " << nb_walls[i] << ", ";
    }
    cout << "\n";

    cout << ind << "regions: ";
    for (region_index_t i: regions) {
      cout << i << ", ";
    }
    cout << "\n";

    cout << ind
      << "uv_vert1_u: " << uv_vert1_u
      << ", uv_vert2: " << uv_vert2
      << ", area: " << area
      << ", normal: " << normal << "\n";

    cout << ind
      << "unit_u: " << unit_u
      << ", unit_v: " << unit_v
      << ", distance_to_origin: " << distance_to_origin
      << "\n";
  }
}


namespace Geometry {


/***************************************************************************
create_region_bbox:
  In: a region
  Out: pointer to a 2-element array contining the LLF and URB corners of
       a bounding box around the region, or NULL if out of memory.
***************************************************************************/
void compute_region_bounding_box(
    const Partition& p, const Region& r,
    Vec3& llf, Vec3& urb
)
{
  bool found_first_wall = false;
  for (auto it: r.walls_and_edges) {
    const Wall& w = p.get_wall(it.first);
    const Vec3* verts[VERTICES_IN_TRIANGLE];
    verts[0] = &p.get_wall_vertex(w, 0);
    verts[1] = &p.get_wall_vertex(w, 1);
    verts[2] = &p.get_wall_vertex(w, 2);

    if (!found_first_wall) {
      llf = *verts[0];
      urb = *verts[0];
      found_first_wall = true;
    }

    for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
      const Vec3* v = verts[i];
      if (llf.x > v->x) {
        llf.x = v->x;
      }
      else if (urb.x < v->x) {
        urb.x = v->x;
      }

      if (llf.y > v->y) {
        llf.y = v->y;
      }
      else if (urb.y < v->y) {
        urb.y = v->y;
      }

      if (llf.z > v->z) {
        llf.z = v->z;
      }
      else if (urb.z < v->z) {
        urb.z = v->z;
      }
    }
  }
}


/***************************************************************************
eval_rel_region_bbox:
  In: release expression for a 3D region release
      place to store LLF corner of the bounding box for the release
      place to store URB corner
  Out: true on success, false on failure.  Bounding box is set based on release
       expression (based boolean intersection of bounding boxes for each
       region).  The function reports failure if any region is unclosed.
***************************************************************************/
// TODO: not checking whether object is closed - use some library for that
bool compute_region_expr_bounding_box(
    World* world, const RegionExprNode* expr,
    Vec3& llf, Vec3& urb
) {
  int count_regions_flag = 1;
  Partition& p = world->get_partition(0);

  if (expr->op == RegionExprOperator::Leaf) {
    Region& reg = p.get_region_by_id(expr->region_id);
    compute_region_bounding_box(p, reg, llf, urb);
    return true;
  }
  else {
    Vec3 llf_left, urb_left;
    compute_region_expr_bounding_box(world, expr->left, llf_left, urb_left);

    Vec3 llf_right, urb_right;
    compute_region_expr_bounding_box(world, expr->right, llf_right, urb_right);

    llf = llf_left;
    urb = urb_left;

    if (expr->op == RegionExprOperator::Union) {
      // TODO: make some function for it
      if (llf.x > llf_right.x) {
        llf.x = llf_right.x;
      }
      if (llf.y > llf_right.y) {
        llf.y = llf_right.y;
      }
      if (llf.z > llf_right.z) {
        llf.z = llf_right.z;
      }
      if (urb.x < urb_right.x) {
        urb.x = urb_right.x;
      }
      if (urb.y < urb_right.y) {
        urb.y = urb_right.y;
      }
      if (urb.z < urb_right.z) {
        urb.z = urb_right.z;
      }
    }
    else if (expr->op == RegionExprOperator::Difference) {
      // for difference/subtraction the MCell3 implementation returns
      // the llf_left/urb_left
    }
    else if (expr->op == RegionExprOperator::Intersect) {
      if (llf.x < llf_right.x) {
        llf.x = llf_right.x;
      }
      if (llf.y < llf_right.y) {
        llf.y = llf_right.y;
      }
      if (llf.z < llf_right.z) {
        llf.z = llf_right.z;
      }
      if (urb.x > urb_right.x) {
        urb.x = urb_right.x;
      }
      if (urb.y > urb_right.y) {
        urb.y = urb_right.y;
      }
      if (urb.z > urb_right.z) {
        urb.z = urb_right.z;
      }
    }
    else {
      assert(false);
      return false;
    }
    return true;
  }
}



// this is the entry point called from Partition class
void update_moved_walls(
    Partition& p,
    const std::vector<VertexMoveInfo>& scheduled_vertex_moves,
    // we can compute all the information already from scheduled_vertex_moves,
    // but the keys of the map walls_with_their_moves are the walls that we need to update
    const WallsWithTheirMovesMap& walls_with_their_moves
) {

  // move all vertices
  for (const VertexMoveInfo& move_info: scheduled_vertex_moves) {
    Vec3& vertex_ref = p.get_geometry_vertex(move_info.vertex_index);
    vertex_ref = vertex_ref + move_info.translation_vec;
    if (! p.in_this_partition(vertex_ref) ) {
      mcell_log("Error: Crossing partitions is not supported yet.\n");
      exit(1);
    }
  }

  // update walls
  // FIXME: this should be placed in geometry.cpp
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    Wall& w = p.get_wall(wall_index);

    // first we need to update all wall constants
    w.precompute_wall_constants(p);

    // reinitialize grid
    if (w.grid.is_initialized()) {
      w.grid.initialize(p, w);
    }
  }


  // edges need to be fixed after all wall have been moved
  // otherwise the edge initialization would be using
  // inconsistent data
  // we need to update also edges of neighboring walls
  uint_set<wall_index_t> walls_to_be_updated;
  // NOTE: we might consider sharing edges in the same way as in MCell3
  for (auto it: walls_with_their_moves) {
    wall_index_t wall_index = it.first;
    Wall& w = p.get_wall(wall_index);

    walls_to_be_updated.insert(wall_index);
    for (uint n = 0; n < EDGES_IN_TRIANGLE; n++) {
      walls_to_be_updated.insert(w.nb_walls[n]);
    }
  }
  for (wall_index_t wall_index: walls_to_be_updated) {
    if (wall_index == WALL_INDEX_INVALID) {
      continue;
    }
    Wall& w = p.get_wall(wall_index);
    w.reinit_edge_constants(p);
  }
}



/*****************************************************************
check_for_overlapped_walls:
  In: rng: random number generator
      n_subvols: number of subvolumes
      subvol: a subvolume
  Out: 0 if no errors, the world geometry is successfully checked for
       overlapped walls.
       1 if there are any overlapped walls.
******************************************************************/
int check_for_overlapped_walls(World* world) {

  /* pick up a random vector */
  Vec3 rand_vector;
  rand_vector.x = rng_dbl(&world->rng);
  rand_vector.y = rng_dbl(&world->rng);
  rand_vector.z = rng_dbl(&world->rng);

  // TODO, it was included for now only because the rand gen was called
#if 0
  for (int i = 0; i < n_subvols; i++) {
    struct subvolume *sv = &(subvol[i]);
    struct wall_aux_list *head = NULL;

    for (struct wall_list *wlp = sv->wall_head; wlp != NULL; wlp = wlp->next) {
      double d_prod = dot_prod(&rand_vector, &(wlp->this_wall->normal));
      /* we want to place walls with opposite normals into
         neighboring positions in the sorted linked list */
      if (d_prod < 0)
        d_prod = -d_prod;

      struct wall_aux_list *newNode;
      newNode = CHECKED_MALLOC_STRUCT(struct wall_aux_list, "wall_aux_list");
      newNode->this_wall = wlp->this_wall;
      newNode->d_prod = d_prod;

      sorted_insert_wall_aux_list(&head, newNode);
    }

    for (struct wall_aux_list *curr = head; curr != NULL; curr = curr->next) {
      struct wall *w1 = curr->this_wall;

      struct wall_aux_list *next_curr = curr->next;
      while ((next_curr != NULL) &&
             (!distinguishable(curr->d_prod, next_curr->d_prod, EPS_C))) {
        /* there may be several walls with the same (or mirror)
           oriented normals */
        struct wall *w2 = next_curr->this_wall;

        if (are_walls_coplanar(w1, w2, MESH_DISTINCTIVE)) {
          if ((are_walls_coincident(w1, w2, MESH_DISTINCTIVE) ||
               coplanar_tri_overlap(w1, w2))) {
            mcell_error(
                "walls are overlapped: wall %d from '%s' and wall "
                "%d from '%s'.",
                w1->side, w1->parent_object->sym->name, w2->side,
                w2->parent_object->sym->name);
          }
        }
        next_curr = next_curr->next;
      }
    }
    /* free memory */
    if (head != NULL)
      delete_wall_aux_list(head);
  }
#endif
  return 0;
}


} /* namespace Geometry */

} /* namespace MCell */
