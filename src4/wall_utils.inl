/******************************************************************************
 *
 * Copyright (C) 2006-2017,2019-2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_WALL_UTILS_INC_
#define SRC4_WALL_UTILS_INC_

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to give the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 */
#include <vector>
#include <limits.h>

#include "bng/bng.h"

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

#include "geometry_utils.h"

#include "geometry_utils.inl"
#include "collision_utils.inl"
#include "rxn_utils.inl"

namespace MCell {

namespace WallUtils {

/***********************************************************************
walls_share_full_edge:
  In: two walls
  Out: 1 if the walls share a full edge, 0 - otherwise.
       Here by "full" we mean that the shared edge has two endpoints
       that are the vertices of both walls w1 and w2.
************************************************************************/
static bool walls_share_full_edge(const Partition& p, const Wall& w1, const Wall& w2) {
  uint i, k;
  uint count = 0;

  /* count number of shared vertices between two walls */
  for (i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    for (k = 0; k < VERTICES_IN_TRIANGLE; k++) {

      if (!distinguishable_vec3( p.get_wall_vertex(w1, i), p.get_wall_vertex(w2, k), POS_EPS) ) {
        count++;
      }
    }
  }

  return count == 2;
}


/*************************************************************************
find_nbr_walls_shared_one_vertex:
   In: the origin wall
       array with information about which vertices of the origin wall
          are shared with neighbor wall (they are indices in the
          global "world->walls_using_vertex" array).
   Out: linked list of the neighbor walls that have only one common
        vertex with the origin wall (not edge-to-edge walls, but
        vertex-to-vertex walls).
   Note: the "origin" wall is not included in the list
**************************************************************************/
static void find_nbr_walls_shared_one_vertex(
    const Partition& p,
    const Wall& origin_wall,
    vertex_index_t shared_verts[VERTICES_IN_TRIANGLE],
    WallIndicesVector& neighboring_walls
)
{
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    if (shared_verts[i] != VERTEX_INDEX_INVALID) {

      const std::vector<wall_index_t>& wall_indices = p.get_walls_using_vertex(shared_verts[i]);
      for (wall_index_t wi: wall_indices) {
        const Wall& w = p.get_wall(wi);

        if (w.id == origin_wall.id) {
          // we do not care about current wall
          continue;
        }

        if (!walls_share_full_edge(p, origin_wall, w)) {
          neighboring_walls.push_back(w.id);
        }
      }
    }
  }
}

#if 0 // unused for now
/***************************************************************************
wall_contains_both_vertices:
  In: wall
      two vertices
  Out: Returns 1 if the wall contains both vertices above, and 0 otherwise.

  Different from MCell3 implementation, we are just comparing vertex indices,
  this is used for neighboring edge detection and we assume that the
  neighbor is in the same mesh, therefore the vertices must have the same
  index.

***************************************************************************/
static bool wall_contains_both_vertices(
    const Wall& w, const vertex_index_t vert_A, const vertex_index_t vert_B) {

  assert(vert_A != vert_B && "Vertices must be different");

  uint count = 0;

  for (uint i = 0; i < 3; i++) {
    vertex_index_t vi = w.vertex_indices[i];

    if (vi == vert_A || vi == vert_B) {
      count++;
    }
  }

  return count == 2;
}
#endif

/***********************************************************************
find_shared_vertices_for_neighbor_walls:
   In: original wall
       neighbor wall
       index of the neighbor wall vertex that is
           shared with original wall (return value)
       index of the neighbor wall vertex that is
           shared with original wall (return value)
   Out: neighbor wall shared vertices indices are set up.
***********************************************************************/
static void find_shared_vertices_for_neighbor_walls(const Wall& orig_wall,
                                             const Wall& nb_wall,
                                             int& shared_vert_1,
                                             int& shared_vert_2) {

  shared_vert_1 = -1;
  shared_vert_2 = -1;

  if (nb_wall.vertex_indices[0] == orig_wall.vertex_indices[0] ||
      nb_wall.vertex_indices[0] == orig_wall.vertex_indices[1] ||
      nb_wall.vertex_indices[0] == orig_wall.vertex_indices[2]) {
    shared_vert_1 = 0;
  }

  if (nb_wall.vertex_indices[1] == orig_wall.vertex_indices[0] ||
      nb_wall.vertex_indices[1] == orig_wall.vertex_indices[1] ||
      nb_wall.vertex_indices[1] == orig_wall.vertex_indices[2]) {
    if (shared_vert_1 < 0) {
      shared_vert_1 = 1;
    }
    else {
      shared_vert_2 = 1;
    }
  }

  if (nb_wall.vertex_indices[2] == orig_wall.vertex_indices[0] ||
      nb_wall.vertex_indices[2] == orig_wall.vertex_indices[1] ||
      nb_wall.vertex_indices[2] == orig_wall.vertex_indices[2]) {
    if (shared_vert_1 < 0) {
      shared_vert_1 = 2;
    }
    else {
      shared_vert_2 = 2;
    }
  }

  assert(shared_vert_1 != -1);
  assert(shared_vert_2 != -1);
}

/*************************************************************************
find_shared_edge_index_of_neighbor_wall:
  In: original wall
      neighbor wall
  Out: index of the shared edge in the coordinate system of neighbor wall.

**************************************************************************/
static edge_index_t find_shared_edge_index_of_neighbor_wall(
    const Wall& orig_wall,
    const Wall& nbr_wall) {
  edge_index_t nbr_edge_ind;
  int shared_vert_ind_1, shared_vert_ind_2;

  find_shared_vertices_for_neighbor_walls(
      orig_wall, nbr_wall, shared_vert_ind_1, shared_vert_ind_2);

  if ((shared_vert_ind_1 + shared_vert_ind_2) == 1) {
    nbr_edge_ind = EDGE_INDEX_0;
  }
  else if ((shared_vert_ind_1 + shared_vert_ind_2) == 2) {
    nbr_edge_ind = EDGE_INDEX_2;
  }
  else if ((shared_vert_ind_1 + shared_vert_ind_2) == 3) {
    nbr_edge_ind = EDGE_INDEX_1;
  }
  else {
    assert(false);
    nbr_edge_ind = EDGE_INDEX_INVALID;
  }

  return nbr_edge_ind;
}

#if 0
/****************************************************************************
find_neighbor_wall_and_edge:
  In: orig_wall: wall
      orig_edge_ind: wall edge index (in the coordinate system of "wall")
      nbr_wall: neighbor wall (return value)
      nbr_edge_ind: index of the edge in the coordinate system of "neighbor
                    wall" that is shared with "wall" and coincides with the
                    edge with "wall edge index" (return value)

****************************************************************************/
static void find_neighbor_wall_and_edge(
    const Partition& p,
    const Wall& orig_wall, const edge_index_t orig_edge_index,
    wall_index_t& nbr_wall_index, edge_index_t& nbr_edge_index) {

  struct wall *w;
  vertex_index_t vert_A, vert_B;

  switch (orig_edge_index) {
  case EDGE_INDEX_0:
    vert_A = orig_wall.vertex_indices[0];
    vert_B = orig_wall.vertex_indices[1];
    break;
  case EDGE_INDEX_1:
    vert_A = orig_wall.vertex_indices[1];
    vert_B = orig_wall.vertex_indices[2];
    break;
  case EDGE_INDEX_2:
    vert_A = orig_wall.vertex_indices[2];
    vert_B = orig_wall.vertex_indices[0];
    break;
  default:
    assert(false && "Invalid edge index");
  }

  for (uint i = 0; i < 3; i++) {
    wall_index_t wi = orig_wall.nb_walls[i];
    if (wi == WALL_INDEX_INVALID) {
      continue;
    }
    const Wall& w = p.get_wall(wi);
    // NOTE: wall_contains_both_vertices and find_shared_edge_index_of_neighbor_wall can be probably optimized
    // into a single call
    if (wall_contains_both_vertices(w, vert_A, vert_B)) {
      nbr_wall_index = wi;
      nbr_edge_index = find_shared_edge_index_of_neighbor_wall(orig_wall, w);
      return;
    }
  }
}
#endif // unused for now


/***********************************************************************
is_wall_edge_region_border:
  In: wall
      wall's edge
  Out: 1 if the edge is a region's border, and 0 - otherwise.
  Note: we do not specify any particular region here, any region will
        suffice
************************************************************************/
/***********************************************************************
is_wall_edge_restricted_region_border: - when region_must_be_reactive is false
  In: wall
      wall's edge
      surface molecule
  Out: 1 if the edge is a restricted region's border for above surface molecule
       0 - otherwise.
  Note: we do not specify any particular region here, any region will
        suffice for which special reactions (REFL/ABSORB) are defined.
************************************************************************/
static bool is_wall_edge_region_border(
    const Partition& p,
    const Wall& this_wall,
    const edge_index_t edge_index,
    const bool region_must_be_reactive
) {

  for (region_index_t region_index: this_wall.regions) {
    const Region& reg = p.get_region(region_index);

    if (region_must_be_reactive && !reg.is_reactive()) {
      continue;
    }

    bool is_edge = reg.is_edge(this_wall.index, edge_index);
    if (is_edge) {
      return true;
    }
  }

  return false;
}


/***************************************************************************
wall_in_box:
  In: array of pointers to vertices for wall (should be 3)
      normal vector for wall
      distance from wall to origin (point normal form)
      first corner of bounding box
      opposite corner of bounding box
  Out: nonzero if the wall intersects the box.  0 otherwise.
***************************************************************************/
static int wall_in_box(
    const Partition& p,
    const Wall& w,
    const Vec3& llf, const Vec3& urb) {

  const Vec3* vert[VERTICES_IN_TRIANGLE] = {
      &p.get_wall_vertex(w, 0),
      &p.get_wall_vertex(w, 1),
      &p.get_wall_vertex(w, 2)
  };

  /* Check if any vertex of the wall is in the box. */
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    if (point_in_box(*vert[i], llf, urb)) {
      return 1;
    }
  }

  /* Check if any wall edge intersects any face of the box */
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    pos_t r, a3, a4;

    const Vec3 *v1, *v2;
    v2 = vert[i];
    v1 = (i == 0) ? vert[VERTICES_IN_TRIANGLE - 1] : vert[i - 1];

		// NOTE: already tried to simplify the code below, but left it like it in the end
    /* x-faces */
    if ((v1->x <= llf.x && llf.x < v2->x) ||
        (v1->x > llf.x && llf.x >= v2->x)) {
      r = (llf.x - v1->x) / (v2->x - v1->x);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->z + r * (v2->z - v1->z);
      if (llf.y <= a3 && a3 <= urb.y && llf.z <= a4 && a4 <= urb.z)
        return 2;
    }
    if ((v1->x <= urb.x && urb.x < v2->x) ||
        (v1->x > urb.x && urb.x >= v2->x)) {
      r = (urb.x - v1->x) / (v2->x - v1->x);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->z + r * (v2->z - v1->z);
      if (llf.y <= a3 && a3 <= urb.y && llf.z <= a4 && a4 <= urb.z)
        return 3;
    }

    /* y-faces */
    if ((v1->y <= llf.y && llf.y < v2->y) ||
        (v1->y > llf.y && llf.y >= v2->y)) {
      r = (llf.y - v1->y) / (v2->y - v1->y);
      a3 = v1->x + r * (v2->x - v1->x);
      a4 = v1->z + r * (v2->z - v1->z);
      if (llf.x <= a3 && a3 <= urb.x && llf.z <= a4 && a4 <= urb.z)
        return 4;
    }
    if ((v1->y <= urb.y && urb.y < v2->y) ||
        (v1->y > urb.y && urb.y >= v2->y)) {
      r = (urb.y - v1->y) / (v2->y - v1->y);
      a3 = v1->x + r * (v2->x - v1->x);
      a4 = v1->z + r * (v2->z - v1->z);
      if (llf.x <= a3 && a3 <= urb.x && llf.z <= a4 && a4 <= urb.z)
        return 5;
    }

    /* z-faces */
    if ((v1->z <= llf.z && llf.z < v2->z) ||
        (v1->z > llf.z && llf.z >= v2->z)) {
      r = (llf.z - v1->z) / (v2->z - v1->z);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->x + r * (v2->x - v1->x);
      if (llf.y <= a3 && a3 <= urb.y && llf.x <= a4 && a4 <= urb.x)
        return 6;
    }
    if ((v1->z <= urb.z && urb.z < v2->z) ||
        (v1->z > urb.z && urb.z >= v2->z)) {
      r = (urb.z - v1->z) / (v2->z - v1->z);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->x + r * (v2->x - v1->x);
      if (llf.y <= a3 && a3 <= urb.y && llf.x <= a4 && a4 <= urb.x)
        return 7;
    }
  }

  /* Lookup table for vertex-edge mapping for a cube */
  int which_x1[12] = { 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1 };
  int which_y1[12] = { 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0 };
  int which_z1[12] = { 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0 };
  int which_x2[12] = { 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0 };
  int which_y2[12] = { 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1 };
  int which_z2[12] = { 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0 };

  int edge1_vt[12] = { 0, 1, 3, 2, 6, 7, 5, 4, 0, 1, 3, 4 };
  int edge2_vt[12] = { 1, 3, 2, 6, 7, 5, 4, 0, 2, 5, 7, 2 };

  const Vec3& normal = w.normal;
  pos_t d = w.distance_to_origin;


  /* Check if any box edge intersects the wall */

  int n_opposite = 0;
  pos_t vu_[VERTICES_IN_TRIANGLE * 2]; /* Assume wall has 3 vertices */
  pos_t* vv_ = &(vu_[VERTICES_IN_TRIANGLE]);

  /* Wall coordinate system n,u,v */
  Vec3 n = normal;
  Vec3 u = *vert[1] - *vert[0];
  pos_t len_u_squared = len3_squared(u);
  assert(len_u_squared != 0);
  pos_t r_u = 1 / sqrt_p(len_u_squared);

  u = u * Vec3(r_u);
  Vec3 v = cross(n, u);
  for (uint j = 0; j < VERTICES_IN_TRIANGLE; j++) {
    vu_[j] = dot(*vert[j], u);
    vv_[j] = dot(*vert[j], v);
  }

  /* Test every edge. */
  Vec3 bb = llf;
  pos_t d_box[8];
  d_box[0] = dot(bb, n);;

  for (uint i = 0; i < 12; i++) {
    pos_t a1, a2;
    Vec3 ba;

    if (i < 7) /* Visiting new vertices in order */
    {
      ba = bb;
      bb.x = (which_x2[i]) ? urb.x : llf.x;
      bb.y = (which_y2[i]) ? urb.y : llf.y;
      bb.z = (which_z2[i]) ? urb.z : llf.z;
      a2 = d_box[edge2_vt[i]] = dot(bb, n);
      a1 = d_box[edge1_vt[i]];

      if ((a1 - d < 0 && a2 - d < 0) || (a1 - d > 0 && a2 - d > 0)) {
        continue;
      }
      else {
        n_opposite++;
      }
    }
    else { /* Revisiting old vertices out of order */
      /*      if (!n_opposite) return 0; */
      a1 = d_box[edge1_vt[i]];
      a2 = d_box[edge2_vt[i]];

      if ((a1 - d < 0 && a2 - d < 0) || (a1 - d > 0 && a2 - d > 0))
        continue;

      n_opposite++;
      ba.x = (which_x1[i]) ? urb.x : llf.x;
      ba.y = (which_y1[i]) ? urb.y : llf.y;
      ba.z = (which_z1[i]) ? urb.z : llf.z;
      bb.x = (which_x2[i]) ? urb.x : llf.x;
      bb.y = (which_y2[i]) ? urb.y : llf.y;
      bb.z = (which_z2[i]) ? urb.z : llf.z;
    }
    /* Now ba,bb = box edge endpoints ; a1,a2 = distances along wall normal */
    pos_t r = (d - a1) / (a2 - a1);
    Vec3 c = ba + Vec3(r) * (bb - ba);
    pos_t cu = dot(c, u);
    pos_t cv = dot(c, v);
    /* Test for internal intersection point in wall coordinate space */
    int temp = 0;
    for (uint j = 0; j < VERTICES_IN_TRIANGLE; j++) {
      uint k = (j == 0) ? VERTICES_IN_TRIANGLE - 1 : j - 1;
      if ((vu_[k] < cu && cu <= vu_[j]) || (vu_[k] >= cu && cu > vu_[j])) {
        r = (cu - vu_[k]) / (vu_[j] - vu_[k]);
        if ((vv_[k] + r * (vv_[j] - vv_[k])) > cv)
          temp++;
      }
    }
    if (temp & 1)
      return 8 + i;
  }

  return 0;
}



/***********************************************************************
find_restricted_regions_by_wall:
  In: wall
      surface molecule
  Out: an object's region list if the wall belongs to the region
          that is restrictive (REFL/ABSORB) to the surface molecule
       NULL - if no such regions found
  Note: regions called "ALL" or the ones that have ALL_ELEMENTS are not
        included in the return "region list".
************************************************************************/
static void find_restricted_regions_by_wall(
    Partition& p,
    const Wall& this_wall,
    const Molecule& sm,
	RegionIndicesSet& regions) {

  assert(regions.empty() && "Function does not clean the resulting array");
  assert(sm.is_surf());

  const BNG::Species& species = p.get_species(sm.species_id);

  if (!species.can_interact_with_border()) {
    // no restricted region
    return;
  }

  // get all possible molecule-surface reactions
  BNG::RxnClassesVector matching_rxn_classes;
  RxnUtils::trigger_intersect(p, sm, sm.s.orientation, this_wall, true, matching_rxn_classes);

  // collect species that can react
  uint_set<species_id_t> restricted_surf_classes;
  for (const BNG::RxnClass* rxn_class: matching_rxn_classes) {
    release_assert(rxn_class->is_simple() && "BNG not supported here yet");
    if (rxn_class->is_reflect_type() || rxn_class->is_absorb_region_border_type_incl_all_molecules()) {
      assert(rxn_class->is_bimol());

      // rxn class may use ALL_MOLECULES or ALL_SURFACE_MOLECULES therefore we cannot use 
      // rxn_class->get_second_species_id
      species_id_t second_reactant_id = rxn_class->get_reactive_surface_reactant_species_id();
      assert(p.get_species(second_reactant_id).is_reactive_surface());

      restricted_surf_classes.insert(second_reactant_id);
    }
  }

  // and store all regions that match the reacting species
  for (region_index_t region_index: this_wall.regions) {
    const Region& reg = p.get_region(region_index);
    if (restricted_surf_classes.count(reg.species_id) != 0) {
      regions.insert(region_index);
    }
  }
}


/*****************************************************************
wall_belongs_to_all_regions_in_region_list:
  In: wall
      region_list
  Out: 1 if wall belongs to all regions in the region list
       0 otherwise.
  Note: Wall can belong to several regions simultaneously.
******************************************************************/
static bool wall_belongs_to_all_regions_in_region_list(
    const Wall& this_wall, const RegionIndicesSet& regions) {

  if (regions.empty()) {
    return false;
  }

  for (region_index_t region_index: regions) {
    if (this_wall.regions.count(region_index) == 0) {
      return false;
    }
  }

  return true;
}


/*****************************************************************
wall_belongs_to_any_region_in_region_list:
  In: wall
      region_list
  Out: 1 if wall belongs to any region in the region list
       0 otherwise.
  Note: Wall can be belong to several regions simultaneously.
  Note: It is assumed that both wall and region list are defined for
        the same object.
******************************************************************/
static bool wall_belongs_to_any_region_in_region_list(
		const Wall& this_wall, const RegionIndicesSet& regions) {
	
  for (region_index_t region_index: regions) {
    if (this_wall.regions.count(region_index) != 0) {
      return true;
    }
  }  

  return false;
}


/*****************************************************************
walls_belong_to_at_least_one_different_restricted_region:
  In: wall and surface molecule on it
      wall and surface molecule on it
  Out: 1 if both walls belong to at least one different restricted region
       relative to the properties of surface molecule,
       0 otherwise.
  Note: Wall can be belong to several regions simultaneously.
        Restricted region is one for the boundary of which the reactions
        REFL/ABSORB are declared.
******************************************************************/
static bool walls_belong_to_at_least_one_different_restricted_region(
    Partition& p,
    const Wall& w1, const Molecule& sm1,
    const Wall& w2, const Molecule& sm2) {
  
  RegionIndicesSet regions1;
  find_restricted_regions_by_wall(p, w1, sm1, regions1);
  RegionIndicesSet regions2;
  find_restricted_regions_by_wall(p, w2, sm2, regions2);

  if (regions1.empty() && regions2.empty()) {
    return false;
  }

  bool res = false;
  if (regions1.empty()) {
    /* Is wall 1 part of all restricted regions rl_2, then these restricted
     * regions just encompass wall 1 */
    if (wall_belongs_to_all_regions_in_region_list(w1, regions2)) {
      return false;
    }
    else {
      return true;
    }

    return res;
  }

  if (regions2.empty()) {
    /* Is wall 2 part of all restricted regions rl_1, then these restricted
     * regions just encompass wall 2 */
    if (wall_belongs_to_all_regions_in_region_list(w2, regions1)){
      return false;
    }
    else {
      return true;
    }
  }

  // return true if the 2nd set of regions does not contain all the regions that 
  // the 1st set contains
  for (region_index_t region_index: regions1) {
    if (regions2.count(region_index) == 0) {
      return true;
    }
  }

  return false;
}



// returns square of the closest distance
// wall must belong to the same object and region
static pos_t find_closest_wall(
    Partition& p,
    const Vec3& pos, const wall_index_t wall_that_moved_molecule,
    const bool must_have_at_least_one_tile,
    wall_index_t& best_wall_index,
    Vec2& best_wall_pos2d
) {

  const Wall& moved_wall = p.get_wall(wall_that_moved_molecule);

  best_wall_index = WALL_INDEX_INVALID;
  pos_t best_d2 = DBL_GIGANTIC;
  best_wall_pos2d = Vec2(0);

  // NOTE: use wall_that_moved_molecule to optimize the search,
  // however for now we are using the same approach as in mcell3 because
  // we need to match the behavior
  const GeometryObject& obj = p.get_geometry_object(moved_wall.object_index);

  for (const wall_index_t wall_index: obj.wall_indices) {
    Wall& w = p.get_wall(wall_index);

    if (!moved_wall.is_same_region(w)) {
      continue;
    }

    if (!w.has_initialized_grid()) {
      w.initialize_grid(p);
    }

    if (must_have_at_least_one_tile && w.grid.num_tiles == 0) {
      continue;
    }

    Vec2 wall_pos2d;
    pos_t d2 = GeometryUtils::closest_interior_point(p, pos, w, wall_pos2d);

    if (d2 <= best_d2) { // the <= is to emulate behavior of mcell3 that goes through the walls in opposite order
      best_d2 = d2;
      best_wall_index = w.index;
      best_wall_pos2d = wall_pos2d;
    }
  }

  return best_d2;
}


// returns square of the closest distance
static pos_t find_closest_wall_any_object(
    Partition& p,
    const Vec3& pos,
    const pos_t search_d2, // squared search diameter
    const bool must_have_at_least_one_tile,
    wall_index_t& best_wall_index,
    Vec2& best_wall_pos2d
) {

  best_wall_index = WALL_INDEX_INVALID;
  pos_t best_d2 = DBL_GIGANTIC;
  best_wall_pos2d = Vec2(0);

  // find which subpartitions we need to check
  pos_t search_d = sqrt_p(search_d2);
  IVec3 subpart_indices;
  p.get_subpart_3d_indices(pos, subpart_indices);
  SubpartIndicesSet crossed_subpart_indices;
  crossed_subpart_indices.insert(p.get_subpart_index_from_3d_indices(subpart_indices));
  CollisionUtils::collect_neighboring_subparts(
      p, pos, subpart_indices, search_d, p.config.subpart_edge_length,
      crossed_subpart_indices
  );

  for (subpart_index_t subpart_index: crossed_subpart_indices) {
    const WallsInSubpart& walls_in_subpart = p.get_subpart_wall_indices(subpart_index);
    for (wall_index_t wi: walls_in_subpart) {
      Wall& w = p.get_wall(wi);

      if (must_have_at_least_one_tile && !w.has_initialized_grid()) {
        w.initialize_grid(p);
      }

      if (must_have_at_least_one_tile && w.grid.num_tiles == 0) {
        continue;
      }

      Vec2 wall_pos2d;
      pos_t d2 = GeometryUtils::closest_interior_point(p, pos, w, wall_pos2d);

      if (d2 <= best_d2 && d2 <= search_d2) { // the <= is to emulate behavior of mcell3 that goes through the walls in opposite order
        best_d2 = d2;
        best_wall_index = w.index;
        best_wall_pos2d = wall_pos2d;
      }
    }
  }

  return best_d2;
}


} // namespace WallUtil

} // namespace MCell

#endif // SRC4_WALL_UTILS_INC_
