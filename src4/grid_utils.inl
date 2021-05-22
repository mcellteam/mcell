/******************************************************************************
 *
 * Copyright (C) 2006-2017,2020-2021 by
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

#ifndef SRC4_GRID_UTILS_INC_
#define SRC4_GRID_UTILS_INC_

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

#include "wall_utils.inl"
#include "geometry_utils.inl"

namespace MCell {

namespace GridUtils {



/*************************************************************************
xyz2grid and uv2grid:
  In: a vector and a surface grid
  Out: int containing the index on the grid of that vector
  Note: xyz2grid just does a dot-product to uv coordinates first.
        Error checking for a valid point is done.
*************************************************************************/
static tile_index_t xyz2grid_tile_index(
    const Partition& p,
    const Vec3& v,
    const Wall& w
) {
  const Grid& g = w.grid;
  assert(g.is_initialized());
  const Vec3& unit_u = w.unit_u;
  const Vec3& unit_v = w.unit_v;

  if (g.num_tiles == 1) {
    return 0;
  }

  tile_index_t tile_idx_0, tile_idx_mid, tile_idx_last;
  /* find tile indices of the corner tiles */
  tile_idx_0 = 0;
  /* see function "move_strip_up()" */
  tile_idx_mid = g.num_tiles - 2 * (g.num_tiles_along_axis) + 1;
  tile_idx_last = g.num_tiles - 1;

  const Vec3& vert_0 = p.get_wall_vertex(w, 0);
  const Vec3& vert_1 = p.get_wall_vertex(w, 1);
  const Vec3& vert_2 = p.get_wall_vertex(w, 2);

  if (!(distinguishable_vec3(v, vert_0, POS_EPS))) {
    return tile_idx_mid;
  }
  if (!(distinguishable_vec3(v, vert_1, POS_EPS))) {
    return tile_idx_last;
  }
  if (!(distinguishable_vec3(v, vert_2, POS_EPS))) {
    return tile_idx_0;
  }

  // check
  if (!(GeometryUtils::point_in_triangle(v, vert_0, vert_1, vert_2))) {
    mcell_internal_error("Error in function 'uv2grid()': point is outside wall.");
  }

  pos_t i = dot(v, unit_u) - g.vert0.u;
  pos_t j = dot(v, unit_v) - g.vert0.v;

  int strip, stripe, flip;
  pos_t striploc, striprem, stripeloc, striperem;
  striploc = j * g.strip_width_rcp;
  strip = (int)striploc;
  striprem = striploc - strip;

  strip = g.num_tiles_along_axis - strip - 1;

  pos_t u0, u1_u0;
  u0 = j * g.vert2_slope;
  u1_u0 = w.uv_vert1_u - j * g.fullslope;

  stripeloc = ((i - u0) / u1_u0) * (strip + (1 - striprem));
  stripe = (int)(stripeloc);
  striperem = stripeloc - stripe;

  flip = (striperem < 1 - striprem) ? 0 : 1;

  tile_index_t idx = strip * strip + 2 * stripe + flip;

  if ((u_int)idx >= g.num_tiles) {
    mcell_internal_error("Error in function 'xyz2grid()': returning tile index "
                         "%d while wall has %u tiles",
                         idx, g.num_tiles);
  }

  return idx;
}

static tile_index_t uv2grid_tile_index(
    const Vec2& v,
    const Wall& w
) {
  const Grid& g = w.grid;
  assert(g.is_initialized());
  assert(g.num_tiles > 0);
  pos_t i, j;
  pos_t u0, u1_u0;
  pos_t striploc, striprem, stripeloc, striperem;
  int strip, stripe, flip, idx;
  Vec2 vert_0, vert_1;
  int tile_idx_0, tile_idx_mid, tile_idx_last;

  if (g.num_tiles == 1) {
    return 0;
  }

  /* find tile indices of the corner tiles */
  tile_idx_0 = 0;
  /* see function "move_strip_up()" */
  tile_idx_mid = g.num_tiles - 2 * g.num_tiles_along_axis + 1;
  tile_idx_last = g.num_tiles - 1;

  vert_0.u = vert_0.v = 0;
  vert_1.u = w.uv_vert1_u;
  vert_1.v = 0;

  if (!distinguishable_vec2(v, vert_0, POS_EPS)) {
    return tile_idx_mid;
  }
  if (!distinguishable_vec2(v, vert_1, POS_EPS)) {
    return tile_idx_0;
  }
  if (!distinguishable_vec2(v, w.uv_vert2, POS_EPS)) {
    return tile_idx_last;
  }

  // check
  if (!(GeometryUtils::point_in_triangle_2D(v, vert_0, vert_1, w.uv_vert2))) {
    mcell_internal_error("Error in function 'uv2grid()': point is outside wall.");
  }

  i = v.u;
  j = v.v;

  striploc = j * g.strip_width_rcp;
  strip = (int)striploc;
  striprem = striploc - strip;

  strip = g.num_tiles_along_axis - strip - 1;

  u0 = j * g.vert2_slope;
  u1_u0 = w.uv_vert1_u - j * g.fullslope;

  stripeloc = ((i - u0) / u1_u0) * (strip + (1 - striprem));
  stripe = (int)(stripeloc);
  striperem = stripeloc - stripe;

  flip = (striperem < 1 - striprem) ? 0 : 1;
  idx = strip * strip + 2 * stripe + flip;

  if ((u_int)idx >= g.num_tiles) {
    mcell_internal_error(
        "Error in function 'uv2grid_tile_index()': returning tile index "
        "%d while wall has %u tiles",
        idx, g.num_tiles
    );
  }

  return idx;
}


/*************************************************************************
grid2xyz and grid2uv and grid2uv_random
  In: a surface grid
      index of a tile on that grid
      vector to store the results
  Out: vector contains the coordinates of the center of that tile, or
       a random coordinate within that tile.
       WARNING: no error checking--index assumed to be valid!
  Note: grid2xyz just multiplies by uv unit vectors at the end.
*************************************************************************/

static Vec3 grid2xyz(Partition& p, const Wall& w, tile_index_t index) {
  const Grid& g = w.grid;
  const Vec3& unit_u = w.unit_u;
  const Vec3& unit_v = w.unit_v;
  int root;
  int rootrem;
  int k, j, i;
  pos_t ucoef, vcoef, over3n;

  root = (int)(sqrt_p(index));
  rootrem = index - root * root;
  k = g.num_tiles_along_axis - root - 1;
  j = rootrem / 2;
  i = rootrem - 2 * j;

  over3n = 1 / (pos_t)(3 * g.num_tiles_along_axis);

  ucoef = ((pos_t)(3 * j + i + 1)) * over3n * w.uv_vert1_u +
          ((pos_t)(3 * k + i + 1)) * over3n * w.uv_vert2.u;
  vcoef = ((pos_t)(3 * k + i + 1)) * over3n * w.uv_vert2.v;

  const Vec3& vert0 = p.get_geometry_vertex(w.vertex_indices[0]);
  Vec3 res;
  res = Vec3(ucoef) * unit_u  + Vec3(vcoef) * unit_v  + vert0;
  return res;
}


static Vec2 grid2uv(const Wall& w, tile_index_t index) {
  const Grid& g = w.grid;
  Vec2 res;
  int root;
  int rootrem;
  int k, j, i;
  pos_t over3n;

  root = (int)(sqrt_p(index));
  rootrem = index - root * root;
  k = g.num_tiles_along_axis - root - 1;
  j = rootrem / 2; // integer division
  i = rootrem - 2 * j;

  over3n = 1 / (pos_t)(3 * g.num_tiles_along_axis);

  res.u = ((pos_t)(3 * j + i + 1)) * over3n * w.uv_vert1_u +
         ((pos_t)(3 * k + i + 1)) * over3n * w.uv_vert2.u;
  res.v = ((pos_t)(3 * k + i + 1)) * over3n * w.uv_vert2.v;
  return res;
}


static Vec2 grid2uv_random(
    const Wall& w, const tile_index_t tile_index,
    rng_state& rng
) {
  const Grid& g = w.grid;
  int root;
  int rootrem;
  int k, j, i;
  pos_t over_n;
  pos_t u_ran, v_ran;

  root = (int)(sqrt_p(tile_index));
  rootrem = tile_index - root * root;
  k = g.num_tiles_along_axis - root - 1;
  j = rootrem / 2;
  i = rootrem - 2 * j;

  over_n = 1 / (pos_t)(g.num_tiles_along_axis);

  u_ran = rng_dbl(&rng);
  v_ran = 1 - sqrt_p(rng_dbl(&rng));

  Vec2 res;
  res.u =
      ((pos_t)(j + i) + (1 - 2 * i) * (1 - v_ran) * u_ran) * over_n * w.uv_vert1_u +
      ((pos_t)(k + i) + (1 - 2 * i) * v_ran) * over_n * w.uv_vert2.u;
  res.v =
      ((pos_t)(k + i) + (1 - 2 * i) * v_ran) * over_n * w.uv_vert2.v;

  return res;
}

/*****************************************************************************
is_inner_tile:
   In: Surface grid
       Index of the tile on that grid
   Out: Returns 1 if the tile is an inner tile
        (not on the border with the neighbor walls).
        Returns 0 if the tile is on the border with the neighbor wall.
*****************************************************************************/
static bool is_inner_tile(const Grid& g, tile_index_t index) {
  int root, rootrem, strip, stripe, flip;

  root = (int)(sqrt_p(index));
  rootrem = index - root * root;
  strip = g.num_tiles_along_axis - root - 1;
  stripe = rootrem / 2;
  flip = rootrem - 2 * stripe;

  if (strip == 0|| stripe == 0) {
    return false;
  }

  if (strip + stripe == (int)g.num_tiles_along_axis - 1) {
    return false;
  }

  if (strip + stripe == (int)g.num_tiles_along_axis - 2 && flip == 1) {
    return false;
  }

  return true;
}


/*****************************************************************************
is_corner_tile:
   In: Surface grid
       Index of the tile on that grid
   Out: Returns 1 if the tile is a corner tile
        (there are only three corners on each wall).
        Returns 0 if the tile is not a corner tile.
*****************************************************************************/
static bool is_corner_tile(const Grid& g, tile_index_t index) {
  if (index == 0 || index == g.num_tiles - 1) {
    return true;
  }

  /* tile index at the wall corner with vertex 0 */
  tile_index_t tile_idx_mid = g.num_tiles - 2 * (g.num_tiles_along_axis) + 1;
  if (index == tile_idx_mid) {
    return true;
  }

  return false;
}


// auxiliary function for grid_neighbors
static void get_grid_neighbors_single_grid_and_index(
    Partition& p,
    const Wall& wall,
    tile_index_t tile_index,
    wall_index_t nb_wall_index,
    wall_index_t &res_nb_wall_index,
    tile_index_t& nb_index
) {
  const Grid& grid = wall.grid;

  if (wall.nb_walls[nb_wall_index] == WALL_INDEX_INVALID) {
    res_nb_wall_index = WALL_INDEX_INVALID;
  }
  else if (!p.get_wall(wall.nb_walls[nb_wall_index]).has_initialized_grid()) {
    res_nb_wall_index = WALL_INDEX_INVALID;
  }
  else {
    molecule_id_t mid = grid.get_molecule_on_tile(tile_index);

    // NOTE: this seems to be the same in all calls
    Vec3 loc3d;
    if (mid != MOLECULE_ID_INVALID) {
      const Molecule& sm = p.get_m(mid);
      const Vec3& wall_vert0 = p.get_geometry_vertex(wall.vertex_indices[0]);
      loc3d = GeometryUtils::uv2xyz(sm.s.pos, wall, wall_vert0);
    }
    else {
      loc3d = grid2xyz(p, wall, tile_index);
    }

    const Wall& nb_wall = p.get_wall(wall.nb_walls[nb_wall_index]);
    Vec2 near2d;
    pos_t d = GeometryUtils::closest_interior_point(p, loc3d, nb_wall, near2d);

    if (!distinguishable_f(d, POS_GIGANTIC, POS_EPS)) {
      res_nb_wall_index = WALL_INDEX_INVALID;
    }
    else {
      res_nb_wall_index = nb_wall.index;
      nb_index = uv2grid_tile_index(near2d, nb_wall);
    }
  }
}

/*************************************************************************
grid_neighbors:
  In: a surface grid
      an index on that grid
      flag that tells whether we have to create a grid on the
          neighbor wall if there is no grid there
      an array[3] of pointers to be filled in with neighboring grid(s)
      an array[3] of pointers to be filled in with neighboring indices
  Out: no return value.  The three nearest neighbors are returned,
       which may be on a neighboring grid if the supplied index is
       at an edge.  If there is no neighbor in one of the three
       directions, the neighboring grid pointer is set to NULL.
  Note: the three neighbors are returned in the same order as the
        edges, i.e. the 0th will be the nearest neighbor in the
        direction of the 0th edge, and so on.
  Note: If this code is used to find neighboring molecules,
        the "create_grid_flag" should be set to zero.
        In such case if a nearby wall exists but has no grid placed
        on it, this function returns NULL for that grid, even though
        there is space there (just no molecules).
        If this code is used to find free spots, the "create_grid_flag"
        should be set to 1 (or any positive value) and the function
        returns newly created grid for this wall.
*************************************************************************/
// original argument create_grid_flag - assumed to be always 0
static void grid_neighbors(
    Partition& p,
    const Wall& wall,
    tile_index_t tile_index,
    // output
    wall_index_t nb_wall[EDGES_IN_TRIANGLE],
    tile_index_t nb_idx[EDGES_IN_TRIANGLE]
) {
  const Grid& grid = wall.grid;

  int i, j, k, root, rootrem;
  Vec3 loc_3d;
  Vec2 near_2d;
  pos_t d;

  /* Calculate strip (k), stripe (j), and flip (i) indices from idx */
  root = (int)(sqrt_p(tile_index));
  rootrem = tile_index - root * root;
  k = root;
  j = rootrem / 2; // integer division
  i = rootrem - 2 * j;

  /* First look "left" (towards edge 2) */
  if (j > 0 || i > 0) /* all tiles except upright tiles in stripe 0 */
  {
    nb_wall[2] = wall.index;
    nb_idx[2] = tile_index - 1;
  }
  else /* upright tiles in stripe 0 */
  {
    get_grid_neighbors_single_grid_and_index(p, wall, tile_index, 2, nb_wall[2], nb_idx[2]);
  }

  /* Then "right" (towards edge 1) */
  if (j < k) /* all tiles except upright tiles in last stripe */
  {
    nb_wall[1] = wall.index;
    nb_idx[1] = tile_index + 1;
  }
  else /* upright tiles in last stripe */
  {
    get_grid_neighbors_single_grid_and_index(p, wall, tile_index, 1, nb_wall[1], nb_idx[1]);
  }

  /* Finally "up/down" (towards edge 0 if not flipped) */
  if (i || k + 1 < (int)grid.num_tiles_along_axis) /* all tiles except upright tiles in last strip */
  {
    nb_wall[0] = wall.index;
    if (i) {
      nb_idx[0] = 2 * j + (k - 1) * (k - 1); /* unflip and goto previous strtile_index*/
    }
    else {
      nb_idx[0] = 1 + 2 * j + (k + 1) * (k + 1); /* flip and goto next strip */
    }
  }
  else /* upright tiles in last strip */
  {
    get_grid_neighbors_single_grid_and_index(p, wall, tile_index, 0, nb_wall[0], nb_idx[0]);
  }
}


/*************************************************************************
tile_orientation:
  In: a vector to the point on the grid and a surface grid
  Out: 0 if the triangle containg the point  is upright,
       and 1 if it is inverted.
       WARNING: no error checking--point assumed to be valid.
*************************************************************************/
static int tile_orientation(const Vec2 v, const Wall& w) {
  const Grid& g = w.grid;
  pos_t i, j;
  pos_t u0, u1_u0;
  pos_t striploc, striprem, stripeloc, striperem;
  int strip, stripe, flip;

  i = v.u;
  j = v.v;

  striploc = j * g.strip_width_rcp;
  strip = (int)striploc;
  striprem = striploc - strip;

  strip = g.num_tiles_along_axis - strip - 1;

  u0 = j * g.vert2_slope;
  u1_u0 = w.uv_vert1_u - j * g.fullslope;

  stripeloc = ((i - u0) / u1_u0) * (((pos_t)strip) + (1 - striprem));
  stripe = (int)(stripeloc);
  striperem = stripeloc - stripe;

  flip = (striperem < 1 - striprem) ? 0 : 1;

  return flip;
}


/*****************************************************************************
move_strip_up:
   In: Wall surface grid and tile index on that grid
   Out: Returns index of the tile that is exactly above the tile
        with index "idx".  Both "idx" and return index should be
        on the same surface grid.
        Returns (-1) if there is no such tile that satisfies the above
        condition (e.g. when the "idx" tile is in the (strip == 0).)
   Note: The direction "up" means moving to the strip where tile indices
         are greater than index "idx".
*****************************************************************************/
static tile_index_t move_strip_up(const Grid& grid, tile_index_t index) {
  int root;
  int tile_up_idx; /* return value */

  root = (int)(sqrt_p(index)) + 1;

  if ((int)grid.num_tiles_along_axis == root) {
    /* tile above is on another wall */
    tile_up_idx = -1;
  } else {
    tile_up_idx = index + 2 * root;
  }

  return tile_up_idx;
}

/*****************************************************************************
move_strip_down:
   In: Wall surface grid and tile index on that grid
   Out: Returns index of the tile that is exactly below the tile
        with index "idx".  Both "idx" and return index should be
        on the same surface grid. Both "idx" and "tile down"
        should share the side.
        Returns (-1) if there is no such tile that satisfies the above
        conditions.
   Note: The direction "down" means moving to the strip where tile indices
         are less than index "idx".
*****************************************************************************/
static tile_index_t move_strip_down(const Grid& grid, tile_index_t index) {
  int root, rootrem, strip, stripe, flip;
  int num_tiles_per_strip;
  tile_index_t tile_down_index; /* return value */

  /* find internal coordinates (strip, stripe, flip) */
  root = (int)(sqrt_p(index));
  rootrem = index - root * root;
  strip = grid.num_tiles_along_axis - root - 1;
  stripe = rootrem / 2;
  flip = rootrem - 2 * stripe;

  num_tiles_per_strip = 2 * (grid.num_tiles_along_axis) - 2 * strip - 1;

  if (is_inner_tile(grid, index)) {
    tile_down_index = index - num_tiles_per_strip + 1;
  } else {
    if ((strip == 0) && (stripe > 0)) {
      if (index == grid.num_tiles - 1) {
        tile_down_index = TILE_INDEX_INVALID;
      } else {
        tile_down_index = index - num_tiles_per_strip + 1;
      }
    } else {
      /* here we are at the left or right border layers */
      if (flip == 0) {
        tile_down_index = TILE_INDEX_INVALID;
      } else {
        tile_down_index = index - num_tiles_per_strip + 1;
      }
    }
  }

  return tile_down_index;
}


/*****************************************************************************
find_shared_vertices_corner_tile_parent_wall:
   In: Surface grid
       Index of the tile on that grid
       3-member array of indices (in the global array "walls_using_vertex"
            of parent wall vertices that are shared with other walls
            (return value)
   Out: Returns 3-member array of the indices of parent wall vertices
        in the global "world->walls_using_vertex" array
        that are shared with other neighbor walls.
        If the wall vertex is not shared the corresponding value
        in the return array is not set.
   Note: Used only for corner tiles.  Here some of the tile vertices
         coincide with the wall vertices which in turn may be shared
         with the neighbor walls.
*****************************************************************************/
// NOTE: can be optimized (somehow) if needed
static bool neighboring_wall_uses_this_vertex(
    const Partition& p,
    const Wall& w,
    const vertex_index_t vi
) {
  for (wall_index_t i = 0; i < EDGES_IN_TRIANGLE; i++) {
    if (w.nb_walls[i] != WALL_INDEX_INVALID) {
      const Wall& nw = p.get_wall(w.nb_walls[i]);
      for (uint svi = 0; svi < VERTICES_IN_TRIANGLE; svi++) {
        if (vi == nw.vertex_indices[svi]) {
          return true;
        }
      }
    }
  }
  return false;
}


static void find_shared_vertices_corner_tile_parent_wall(
    const Partition& p,
    const Wall& w,
    const Grid& g,
    const tile_index_t tile_index,
    vertex_index_t shared_verts[VERTICES_IN_TRIANGLE] // res
) {
  vertex_index_t global_vert_index;

  /* check if we are at vertex 0 */
  if (tile_index == (g.num_tiles - 2 * g.num_tiles_along_axis + 1)) {

    // how to figure out that the vertex is shared? - not sure if the approach is correct one
    // there must be neighboring walls?
    vertex_index_t vi = w.vertex_indices[0];
    if (neighboring_wall_uses_this_vertex(p, w, vi)) {
      shared_verts[0] = vi;
    }
    else {
      shared_verts[0] = VERTEX_INDEX_INVALID;
    }
  }

  /* check if we are at vertex 1 */
  if (tile_index == g.num_tiles - 1) {

    vertex_index_t vi = w.vertex_indices[1];
    if (neighboring_wall_uses_this_vertex(p, w, vi)) {
      shared_verts[1] = vi;
    }
    else {
      shared_verts[1] = VERTEX_INDEX_INVALID;
    }
  }

  /* check if we are at vertex 2 */
  if (tile_index == 0) {
    vertex_index_t vi = w.vertex_indices[2];
    if (neighboring_wall_uses_this_vertex(p, w, vi)) {
      shared_verts[2] = vi;
    }
    else {
      shared_verts[2] = VERTEX_INDEX_INVALID;
    }
  }
}


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
static void find_shared_vertices_for_neighbor_walls(
    const Partition& p,
    const Wall& orig_wall,
    const Wall& nb_wall,
    int& shared_vert_1,
    int& shared_vert_2
) {
  const Vec3& orig_wall_vert0 = p.get_wall_vertex(orig_wall, 0);
  const Vec3& orig_wall_vert1 = p.get_wall_vertex(orig_wall, 1);
  const Vec3& orig_wall_vert2 = p.get_wall_vertex(orig_wall, 2);

  const Vec3& nb_wall_vert0 = p.get_wall_vertex(nb_wall, 0);
  const Vec3& nb_wall_vert1 = p.get_wall_vertex(nb_wall, 1);
  const Vec3& nb_wall_vert2 = p.get_wall_vertex(nb_wall, 2);

  if ((!distinguishable_vec3(nb_wall_vert0, orig_wall_vert0, POS_EPS)) ||
      (!distinguishable_vec3(nb_wall_vert0, orig_wall_vert1, POS_EPS)) ||
      (!distinguishable_vec3(nb_wall_vert0, orig_wall_vert2, POS_EPS))) {
     shared_vert_1 = 0;
  }

  if ((!distinguishable_vec3(nb_wall_vert1, orig_wall_vert0, POS_EPS)) ||
      (!distinguishable_vec3(nb_wall_vert1, orig_wall_vert1, POS_EPS)) ||
      (!distinguishable_vec3(nb_wall_vert1, orig_wall_vert2, POS_EPS))) {
    if (shared_vert_1 < 0) {
      shared_vert_1 = 1;
    }
    else {
      shared_vert_2 = 1;
    }
  }

  if ((!distinguishable_vec3(nb_wall_vert2, orig_wall_vert0, POS_EPS)) ||
      (!distinguishable_vec3(nb_wall_vert2, orig_wall_vert1, POS_EPS)) ||
      (!distinguishable_vec3(nb_wall_vert2, orig_wall_vert2, POS_EPS))) {
    if (shared_vert_1 < 0) {
      shared_vert_1 = 2;
    }
    else {
      shared_vert_2 = 2;
    }
  }
}


/**************************************************************************
grid_all_neighbors_across_walls_through_vertices:
  In: a surface molecule
      linked list of the neighbor walls that share one vertex only
      surface grid of thew all the molecule sits on or where the hit happens
      flag that tells whether we need to create a grid on a neighbor wall
      flag that tells whether we are searching for reactant
          (value = 1) or doing product placement (value = 0)
      a linked list of  neighbor tiles (return value)
      a length of the linked list above (return value)
  Out: The list of nearest neighbors are returned,
       Neighbors should share common vertex.
  Note: This version allows looking for the neighbors at the neighbor walls
       that are connected to the start wall through vertices only. Also
       the function takes care of REFLECTIVE/ABSORPTIVE region borders.
****************************************************************************/
static void grid_all_neighbors_across_walls_through_vertices(
    Partition& p,
    const Molecule* sm, // may be nullptr, not used currently but should be
    const WallIndicesVector& neighboring_walls,
    const Wall& wall,
    bool create_grid_flag,
    bool search_for_reactant,
    TileNeighborVector& neighbors
) {
  const Grid& grid = wall.grid;
  /* index of the neighbor wall vertex in "world->all_vertices" array  that coincides with tile vertex */
  vertex_index_t nbr_wall_vertex_index;
  /* index of the neighbor tile */
  tile_index_t nbr_tile_index;

  /* check for possible reflection (absorption) from the wall edges
     that may be region borders.  This is INSIDE_OUT check against
     molecule's own wall */
  uint_set<region_index_t> restricted_regions;
  if (sm != nullptr && search_for_reactant) {
    const BNG::Species& species = p.get_species(sm->species_id);
    if (species.can_interact_with_border()) {
      WallUtils::find_restricted_regions_by_wall(p, wall, *sm, restricted_regions);
    }
  }

  /* only one corner tile from each neighbor wall
     can be a neighbor to our start tile */

  /* since the neighbor walls are connected to the origin wall by just
     one vertex, and this code is valid only for the corner tile
     on the origin wall, from each neighbor wall we will pick up
     only one corner tile that shares a vertex with the origin wall */
  for (wall_index_t wi: neighboring_walls) {
    Wall& neighboring_wall = p.get_wall(wi);

    if (!neighboring_wall.has_initialized_grid()) {
      if (create_grid_flag) {
        neighboring_wall.initialize_grid(p);
      }
      else {
        continue;
      }
    }

    /* if there is a restricted region list for own wall
       and neighbor wall does NOT belong to all regions in this list -
       we assume that neighbor wall lies outside the
       restricted region boundary and we DO NOT add
       tile on such wall to the list of neighbor tiles  */
    if (search_for_reactant &&
        !restricted_regions.empty() &&
        !WallUtils::wall_belongs_to_all_regions_in_region_list(neighboring_wall, restricted_regions)) {
        continue;
    }

    /* Similar test done OUTSIDE-IN */
    // TODO!
#if 0
    if (sm != NULL) {
      if (search_for_reactant && (sm->properties->flags & SPECIES_FLAG_CAN_REGION_BORDER)) {
        rlp_head_nbr_wall = find_restricted_regions_by_wall(world, neighboring_wall, sm);

        if (rlp_head_nbr_wall != NULL) {
          if (!wall_belongs_to_all_regions_in_region_list(sm->grid->surface,
                                                          rlp_head_nbr_wall)) {
            delete_void_list((struct void_list *)rlp_head_nbr_wall);
            continue;
          } else {
            /* we will add the tile on this wall to the list */
            delete_void_list((struct void_list *)rlp_head_nbr_wall);
          }
        }
      }
    }
#endif

    Grid& neighboring_grid = neighboring_wall.grid;

    /* find the index of the neighbor tile */
    if (grid.num_tiles == 1) {
      // how often is this branch really used?
      nbr_tile_index = 0;
    }
    else {
      nbr_wall_vertex_index = VERTEX_INDEX_INVALID;
      nbr_tile_index = TILE_INDEX_INVALID;

      for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
        for (uint k = 0; k < VERTICES_IN_TRIANGLE; k++) {
          if (wall.vertex_indices[i] == neighboring_wall.vertex_indices[k]) {
            nbr_wall_vertex_index = neighboring_wall.vertex_indices[k];
            break; // shouldn't we terminate the loop completely?
          }
        }
      }

      if (nbr_wall_vertex_index == VERTEX_INDEX_INVALID) {
        mcell_internal_error("Error identifying tile on the neighbor wall.");
      }

      /* find the index of the neighbor tile */
      if (nbr_wall_vertex_index == neighboring_wall.vertex_indices[0]) {
        nbr_tile_index = neighboring_grid.num_tiles - 2 * (neighboring_grid.num_tiles_along_axis) + 1;
      }
      else if (nbr_wall_vertex_index == neighboring_wall.vertex_indices[1]) {
        nbr_tile_index = neighboring_grid.num_tiles - 1;
      }
      else if (nbr_wall_vertex_index == neighboring_wall.vertex_indices[2]) {
        nbr_tile_index = 0;
      }
      if (nbr_tile_index == TILE_INDEX_INVALID) {
        mcell_internal_error("Error identifying tile on the neighbor wall.");
      }
    }

    neighbors.push_front(WallTileIndexPair(neighboring_wall.index, nbr_tile_index));
  }
}


/*************************************************************************
bisect:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      pos_t we are using to bisect the array
  Out: index of the largest element in the array smaller than the bisector
*************************************************************************/
static int bisect(const std::vector<pos_t>& list, int n, pos_t val) {
  int lo = 0;
  int hi = n;
  int mid = 0;
  while (hi - lo > 1) {
    mid = (hi + lo) / 2;
    if (list[mid] > val) {
      hi = mid;
    }
    else {
      lo = mid;
    }
  }
  return lo;
}


/*************************************************************************
bisect_high:
  In: array of floats, sorted low to high
      int saying how many floats there are
      pos_t we are using to bisect the array
  Out: index of the smallest element in the array larger than the bisector
*************************************************************************/
static int bisect_high(const std::vector<pos_t>& list, int n, pos_t val) {
  int lo = 0;
  int hi = n - 1;
  int mid = 0;
  while (hi - lo > 1) {
    mid = (hi + lo) / 2;
    if (list[mid] > val) {
      hi = mid;
    }
    else {
      lo = mid;
    }
  }
  if (list[lo] > val) {
    return lo;
  }
  else {
    return hi;
  }
}


static void push_front_neighbor_if_not_present(
    TileNeighborVector& neighbors,
    const WallTileIndexPair& info
) {
  auto it = std::find(neighbors.begin(), neighbors.end(), info);
  if (it == neighbors.end()) {
    neighbors.push_front(info);
  }
}


/**************************************************************************
add_more_tile_neighbors_to_list_fast:
  In: a linked list of  neighbor tiles
      a surface grid of the start tile
      (strip, stripe, flip) coordinates of the start tile
      3D coordinates of the shared edge
          (edge has a direction from "start" to "end")
      index of the shared edge looking from the original orid
          (0 - for edge between vert[0] and vert[1],
           1 - for edge between vert[1] and vert[2],
           2 - for edge between vert[0] and vert[2])
      a surface grid of the neighbor wall
  Out: Returns number of new neighbor tiles added to the original
       linked list of neighbor tiles.
       Neighbors should share either common edge or common vertex.
  Note:  The function looks only for the tiles that are the neighbors of
        "orig_grid/orig_idx" and reside on the neighbor "new_grid".
***************************************************************************/
static int add_more_tile_neighbors_to_list_fast(
    const Partition& p,
    const Wall& orig_wall,
    int orig_strip,
    uint orig_stripe,
    int orig_flip,
    const Vec3& start, const Vec3& end,
    int edge_index,
    const Wall& new_wall,
    TileNeighborVector& neighbors
) {
  assert(orig_wall.index != new_wall.index
      && "Function 'add_more_tile_neighbors_to_list()' should be called for different grids (walls) only");

  const Grid& orig_grid = orig_wall.grid;
  const Grid& new_grid = new_wall.grid;

  bool invert_orig_pos = false; /* flag */
  bool check_side_flag;     /* flag */
  pos_t edge_length;      /* length of the shared edge */
  /* positions of the shared vertices along the shared edge,
     measured relative to the shared edge length for the original tile */
  pos_t orig_pos_1 = -1, orig_pos_2 = -1;
  /* number of tile vertices on the common edge */
  const int new_pos_size = new_grid.num_tiles_along_axis + 1;
  /* array of the positions of tile vertices on the common edge */
  std::vector<pos_t> new_pos(new_pos_size);

  /* each tile vertex on the common shared edge is connected to
     3 tiles (the end points of the shared edge are connected
     to 1 tile). */
  /* 2-dimensional array of the tile indices */
  std::vector<std::vector<tile_index_t>> new_tile_indices(new_pos_size);
  for (auto& elem: new_tile_indices) {
    elem.resize(3);
  }

  int i, k;

  /* what indices of the shared vertices refer to the vertex "start"
     and "end" for "original" wall */
  int new_start_index, new_end_index;

  int tiles_added = 0; /* counter of added tiles */

  /* find out relative positions of vertices of original tile
     on the common edge */
  edge_length = distance3(start, end);
  if (orig_stripe == 0) {
    if (orig_strip > 0) {
      if (orig_flip == 0) {
        orig_pos_1 = orig_strip * edge_length / (orig_grid.num_tiles_along_axis);
        orig_pos_2 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
      }
      else { /* (orig_flip == 1) */
        orig_pos_1 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
      }
    }
    else {
      /* find out common edge refers to what side of the original wall */
      if (edge_index == 0) {
        if (orig_flip == 0) {
          orig_pos_1 = orig_stripe * edge_length / (orig_grid.num_tiles_along_axis);
          orig_pos_2 = (orig_stripe + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }
        else { /* (orig_flip == 1) */
          orig_pos_1 = (orig_stripe + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }
      }
      else if (edge_index == 1) {
        if (orig_flip == 0) {
          orig_pos_1 = (orig_strip) * edge_length / (orig_grid.num_tiles_along_axis);
        }
        else { /* (orig_flip == 1) */
          orig_pos_1 = orig_strip * edge_length / (orig_grid.num_tiles_along_axis);
          orig_pos_2 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }
      }
      else if (edge_index == 2) {
        if (orig_flip == 0) {
          orig_pos_1 = orig_strip * edge_length / (orig_grid.num_tiles_along_axis);
          orig_pos_2 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }
        else { /* (orig_flip == 1) */
          orig_pos_1 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }

      }
      else {
        assert(false && "Error in the function 'add_more_tile_neighbors_to_list_fast()'.");
      }
    }
  }

  check_side_flag = 0;
  if ((orig_strip == 0) && (orig_stripe > 0)) {
    if (orig_stripe == orig_grid.num_tiles_along_axis - 1) {
      check_side_flag = 1;
    }
    if ((orig_stripe == orig_grid.num_tiles_along_axis - 2) && (orig_flip == 1)) {
      check_side_flag = 1;
    }
    if (!check_side_flag) {
      if (orig_flip == 0) {
        orig_pos_1 = orig_stripe * edge_length / (orig_grid.num_tiles_along_axis);
        orig_pos_2 = (orig_stripe + 1) * edge_length / (orig_grid.num_tiles_along_axis);
      }
      else { /* (orig_flip == 1) */
        orig_pos_1 = (orig_stripe + 1) * edge_length / (orig_grid.num_tiles_along_axis);
      }
    }
    else {
      /* find out common edge refers to what side of the original wall */
      if (edge_index == 0) {
        if (orig_flip == 0) {
          orig_pos_1 = orig_stripe * edge_length / (orig_grid.num_tiles_along_axis);
          orig_pos_2 = (orig_stripe + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }
        else { /* (orig_flip == 1) */
          orig_pos_1 = (orig_stripe + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }
      }
      else if (edge_index == 1) {
        if (orig_flip == 0) {
          orig_pos_1 = orig_strip * edge_length / (orig_grid.num_tiles_along_axis);
          orig_pos_2 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }
        else { /* (orig_flip == 1) */
          orig_pos_1 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
        }
      }
      else {
        assert(false && "Error in the function 'add_more_tile_neighbors_to_list_fast()'.");
      }
    }
  }

  if ((orig_strip > 0) && (orig_stripe > 0)) {
    if (orig_flip == 0) {
      orig_pos_1 = orig_strip * edge_length / (orig_grid.num_tiles_along_axis);
      orig_pos_2 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
    }
    else { /* (orig_flip == 1) */
      orig_pos_1 = (orig_strip + 1) * edge_length / (orig_grid.num_tiles_along_axis);
    }
  }

  /* what vertices of new wall are shared with original wall */
  int shared_vert_1 = -1, shared_vert_2 = -1;

  find_shared_vertices_for_neighbor_walls(
      p, orig_wall, new_wall, shared_vert_1, shared_vert_2);

  /* set the value of 'invert_orig_pos' flag */
  if (!distinguishable_vec3(start, p.get_wall_vertex(new_wall, shared_vert_1), POS_EPS)) {
    new_start_index = shared_vert_1;
    new_end_index = shared_vert_2;
  }
  else {
    new_start_index = shared_vert_2;
    new_end_index = shared_vert_1;
  }

  if (new_start_index > new_end_index) {
    invert_orig_pos = 1;
  }

  /* our coordinate system here is effectively one-dimensional
     shared edge. For the neighbor grid it's orientation is opposite
     to the one we used for the original grid. Let's recalculate
     the positions of the common points on the shared edge depending
     on the neigbor grid shared edge direction */
  if (invert_orig_pos) {
    orig_pos_1 = edge_length - orig_pos_1;
    if (orig_pos_2 > 0) {
      orig_pos_2 = edge_length - orig_pos_2;
    }
  }

  /* find out relative positions of vertices of tile structure
     on the common edge for "new_grid" */
  for (i = 0; i < new_pos_size; i++) {
    new_pos[i] = i * edge_length / (new_grid.num_tiles_along_axis);
  }

  /* index of the shared edge in terms of the "new_grid".
   Here the edge between vert[0] and vert[1] has index 0,
   the edge between vert[1] and vert[2] has index 1,
   and the edge between vert[0] and vert[2] has index 2. */
  int new_edge_index = 0;
  if ((shared_vert_1 + shared_vert_2) == 1) {
    new_edge_index = 0;
  }
  else if ((shared_vert_1 + shared_vert_2) == 2) {
    new_edge_index = 2;
  }
  else if ((shared_vert_1 + shared_vert_2) == 3) {
    new_edge_index = 1;
  }
  else {
    assert(false && "Error in the function 'add_more_tile_neighbors_to_list_fast()'.");
  }

  /* fill out the array with tile indices for the border layer
     adjacent to the shared edge */
  int last_value;
  if (new_edge_index == 0) {
    new_tile_indices[0][0] = -1;
    new_tile_indices[0][1] = -1;
    new_tile_indices[0][2] = new_grid.num_tiles - 2 * (new_grid.num_tiles_along_axis) + 1;
    last_value = new_tile_indices[0][2];

    for (i = 1; i < new_pos_size - 1; i++) {
      for (k = 0; k < 3; k++) {
        new_tile_indices[i][k] = last_value + k;
      }
      last_value = new_tile_indices[i][2];
    }

    new_tile_indices[new_pos_size - 1][0] = last_value;
    new_tile_indices[new_pos_size - 1][1] = -1;
    new_tile_indices[new_pos_size - 1][2] = -1;

  }
  else if (new_edge_index == 1) {
    new_tile_indices[0][0] = -1;
    new_tile_indices[0][1] = -1;
    new_tile_indices[0][2] = new_grid.num_tiles - 1;
    last_value = new_tile_indices[0][2];

    for (i = 1; i < new_pos_size - 1; i++) {
      for (k = 0; k < 2; k++) {
        new_tile_indices[i][k] = last_value - k;
      }
      last_value = new_tile_indices[i][1];
      new_tile_indices[i][2] = move_strip_down(new_grid, last_value);
      last_value = new_tile_indices[i][2];
    }
    new_tile_indices[new_pos_size - 1][0] = last_value;
    new_tile_indices[new_pos_size - 1][1] = -1;
    new_tile_indices[new_pos_size - 1][2] = -1;

  }
  else { /* (new_edge_index == 2) */
    new_tile_indices[0][0] = -1;
    new_tile_indices[0][1] = -1;
    new_tile_indices[0][2] = new_grid.num_tiles - 2 * (new_grid.num_tiles_along_axis) + 1;
    last_value = new_tile_indices[0][2];

    for (i = 1; i < new_pos_size - 1; i++) {
      for (k = 0; k < 2; k++) {
        new_tile_indices[i][k] = last_value + k;
      }
      last_value = new_tile_indices[i][1];
      new_tile_indices[i][2] = move_strip_down(new_grid, last_value);
      last_value = new_tile_indices[i][2];
    }
    new_tile_indices[new_pos_size - 1][0] = last_value;
    new_tile_indices[new_pos_size - 1][1] = -1;
    new_tile_indices[new_pos_size - 1][2] = -1;
  }

  int ind_high, ind_low = -1;
  if (orig_pos_1 > orig_pos_2) {
    ind_high = bisect_high(new_pos, new_pos_size, orig_pos_1); // mcell3 function
    if (orig_pos_2 > 0) {
      ind_low = bisect(new_pos, new_pos_size, orig_pos_2); // mcell3 function
    }

  }
  else {
    ind_high = bisect_high(new_pos, new_pos_size, orig_pos_2); // mcell3 function

    if (orig_pos_1 > 0) {
      ind_low = bisect(new_pos, new_pos_size, orig_pos_1); // mcell3 function
    }
  }

  if (ind_low >= 0) {
    for (i = ind_low + 1; i < ind_high; i++) {
      for (k = 0; k < 3; k++) {
        push_front_neighbor_if_not_present(
            neighbors, WallTileIndexPair(new_wall.index, new_tile_indices[i][k]));
      }
    }

  }
  else {
    push_front_neighbor_if_not_present(
        neighbors, WallTileIndexPair(new_wall.index, new_tile_indices[ind_high][0]));
  }

  return tiles_added;
}


// w might be nullptr
static bool wall_has_initialized_grid(const Wall* w) {
  return w != nullptr && w->has_initialized_grid();
}


static bool check_if_can_move_through_border(
    Partition& p,
    const Molecule& sm,
    const Wall& wall,
    const uint_set<region_index_t>& this_wall_restricted_regions,
    const Wall* neighbor_wall) {

  if (neighbor_wall != NULL) {
    if (!WallUtils::wall_belongs_to_all_regions_in_region_list(*neighbor_wall, this_wall_restricted_regions)) {
      return false;
    }

    uint_set<region_index_t> nb_regions;
    WallUtils::find_restricted_regions_by_wall(p, *neighbor_wall, sm, nb_regions);
    if (!WallUtils::wall_belongs_to_all_regions_in_region_list(*neighbor_wall, nb_regions)) {
      return false;
    }
  }

  return true;
}

/**************************************************************************
grid_all_neighbors_across_walls_through_edges:
  In: a surface molecule
      surface grid of the wall where molecule sits on or where hit happens
      index of the tile for the condition above
      flag that tells whether we need to create a grid on a neighbor wall
      flag that tells whether we are searching for reactant
          (value = 1) or doing product placement (value = 0)
      a linked list of  neighbor tiles (return value)
      a length of the linked list above (return value)
  Out: The list of nearest neighbors are returned,
       Neighbors should share common edge.
  Note: This version allows looking for the neighbors at the neighbor walls
        that are connected to the start wall through edges only.
****************************************************************************/
static void grid_all_neighbors_across_walls_through_edges(
    Partition& p,
    const Molecule* sm, // may be nullptr
    const Wall& wall,
    const tile_index_t tile_index,
    bool create_grid_flag,
    bool search_for_reactant,
    TileNeighborVector& neighbors
) {
  const Grid& grid = wall.grid;
  assert(tile_index < grid.num_tiles && "Surface molecule tile index is out of bounds for this grid");


  int root, rootrem, strip, stripe, flip;
  tile_index_t temp_index;

  /* flags */
  bool move_thru_border_0 = true;
  bool move_thru_border_1 = true;
  bool move_thru_border_2 = true;

  // the wall pointers might be nullptrs
  Wall* nb_walls[EDGES_IN_TRIANGLE];
  nb_walls[0] = p.get_wall_if_exists(wall.nb_walls[0]);
  nb_walls[1] = p.get_wall_if_exists(wall.nb_walls[1]);
  nb_walls[2] = p.get_wall_if_exists(wall.nb_walls[2]);

  /* check for possible reflection (absorption) from the wall edges
     that may be region borders.  These are INSIDE_OUT and OUTSIDE-IN
     checks against molecule's own wall and neighbor wall */
  if (sm != nullptr && search_for_reactant) {
    const BNG::Species& species = p.get_species(sm->species_id);
    if (species.can_interact_with_border()) {
      assert(sm->s.wall_index == wall.index);
  
      uint_set<region_index_t> regions;
      WallUtils::find_restricted_regions_by_wall(p, wall, *sm, regions);
  
      move_thru_border_0 = check_if_can_move_through_border(p, *sm, wall, regions, nb_walls[0]);
      move_thru_border_1 = check_if_can_move_through_border(p, *sm, wall, regions, nb_walls[1]);
      move_thru_border_2 = check_if_can_move_through_border(p, *sm, wall, regions, nb_walls[2]);
    }
  }
  
  /* find (strip, stripe, flip) coordinates of the tile */
  root = (int)(sqrt_p(tile_index));
  rootrem = tile_index - root * root;
  strip = grid.num_tiles_along_axis - root - 1;
  stripe = rootrem / 2;
  flip = rootrem - 2 * stripe;

  const Vec3& wall_vert0 = p.get_wall_vertex(wall, 0);
  const Vec3& wall_vert1 = p.get_wall_vertex(wall, 1);
  const Vec3& wall_vert2 = p.get_wall_vertex(wall, 2);

  if (create_grid_flag) {
    for (uint kk = 0; kk < EDGES_IN_TRIANGLE; kk++) {
      if (nb_walls[kk] != nullptr && !nb_walls[kk]->has_initialized_grid()) {
        nb_walls[kk]->initialize_grid(p);
      }
    }
  }

  wall_index_t wall_index = wall.index;
  if (stripe == 0) {
    if (flip > 0) /* inverted tile */
    {
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 1));
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 1));
      if (strip < (int)grid.num_tiles_along_axis - 2) {
        neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 2));
      }

      /* put in the list tiles that are on the row below the start tile
         but on the same grid */
      temp_index = move_strip_down(grid, tile_index);
      assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

      neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
      if (strip < (int)grid.num_tiles_along_axis - 2) {
        neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 1));
        neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 2));
      }

      if (strip > 0) {
        /* put in the list tiles that are on the row above the start tile
           but on the same grid */
        temp_index = move_strip_up(grid, tile_index);
        assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

        neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
        neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 1));
        neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 1));
      }

      /* get the neighbors from the neighbor walls */
      if (wall_has_initialized_grid(nb_walls[2]) && move_thru_border_2) {
        add_more_tile_neighbors_to_list_fast(
            p,
            wall,
            strip, stripe, flip, wall_vert0, wall_vert2, 2,
            *nb_walls[2],
            neighbors
        );
      }

      if (strip == 0 && wall_has_initialized_grid(nb_walls[0]) && move_thru_border_0) {
        add_more_tile_neighbors_to_list_fast(
            p,
            wall,
            strip, stripe, flip, wall_vert0, wall_vert1, 0,
            *nb_walls[0],
            neighbors
        );
      }

      if (strip == (int)grid.num_tiles_along_axis - 2 && wall_has_initialized_grid(nb_walls[1]) && move_thru_border_1) {
        add_more_tile_neighbors_to_list_fast(
            p,
            wall,
            strip, stripe, flip, wall_vert1, wall_vert2, 0,
            *nb_walls[1],
            neighbors
        );
      }

    }
    else {        /* upright tile (flip == 0)  */
      if (tile_index == 0) /* it is a special case */
      {
        if (grid.num_tiles > 1) {
          /* put in the list tiles that are on the row above the start tile */
          temp_index = move_strip_up(grid, tile_index);
          assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

          neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
          neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 1));
          neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 1));
        }
        else {
          if (wall_has_initialized_grid(nb_walls[0]) && move_thru_border_0) {
            add_more_tile_neighbors_to_list_fast(
                p,
                wall,
                strip, stripe, flip, wall_vert0, wall_vert1, 0,
                *nb_walls[0],
                neighbors
            );
          }
        }

        if (wall_has_initialized_grid(nb_walls[1]) && move_thru_border_1) {
          add_more_tile_neighbors_to_list_fast(
              p,
              wall,
              strip, stripe, flip, wall_vert1, wall_vert2, 1,
              *nb_walls[1],
              neighbors
          );
        }

        if (wall_has_initialized_grid(nb_walls[2]) && move_thru_border_2) {
          add_more_tile_neighbors_to_list_fast(
              p,
              wall,
              strip, stripe, flip, wall_vert0, wall_vert2, 2,
              *nb_walls[2],
              neighbors
          );
        }

      }
      else { /* if (tile_index != 0) */

        neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 1));
        neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 2));

        /* put in the list tiles that are on the row below the start tile */
        temp_index = move_strip_down(grid, tile_index + 1);
        assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

        neighbors.push_front(WallTileIndexPair(wall_index, temp_index));

        if (strip > 0) {
          /* put in the list tiles that are on the row above the start tile
             but on the same grid */
          temp_index = move_strip_up(grid, tile_index);
          assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

          neighbors.push_front(WallTileIndexPair(wall_index,  temp_index));
          neighbors.push_front(WallTileIndexPair(wall_index,  temp_index - 1));
          neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 1));
          neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 2));
        }
        else { /* strip == 0 */
          /* put in the list tiles that are on the row above the start tile
          but on the different grid */
          /* it is the top left corner - special case */
          if (wall_has_initialized_grid(nb_walls[0]) && move_thru_border_0) {
            add_more_tile_neighbors_to_list_fast(
                p,
                wall,
                strip, stripe, flip, wall_vert0, wall_vert1, 0,
                *nb_walls[0],
                neighbors
            );
          }

          if (wall_has_initialized_grid(nb_walls[2]) && move_thru_border_2) {
            add_more_tile_neighbors_to_list_fast(
                p,
                wall,
                strip, stripe, flip, wall_vert0, wall_vert2, 2,
                *nb_walls[2],
                neighbors
            );
          }
        }

      } /* end if-else (tile_index == 0) */

    } /* end upright or inverted tile */

  } /* end if (stripe == 0) */

  if ((strip == 0) && (stripe > 0)) {
    /* put in the list tiles that are on the same row */
    neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 1));
    neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 2));

    if ((stripe < (int)grid.num_tiles_along_axis - 2) || ((stripe == (int)grid.num_tiles_along_axis - 2) && (flip == 0))) {
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 1));
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 2));
    }
    else if ((stripe == (int)grid.num_tiles_along_axis - 2) && (flip == 1)) {
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 1));
    }

    /* put in the list tiles that are on the row below */
    if (flip > 0) {
      temp_index = move_strip_down(grid, tile_index);
      assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

      neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 1));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 2));
      if (stripe < (int)grid.num_tiles_along_axis - 2) {
        neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 1));
        neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 2));
      }
    }
    else { /* (flip == 0) */
      if (tile_index < grid.num_tiles - 1) {
        temp_index = move_strip_down(grid, tile_index);
        assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

        neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
        neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 1));
        neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 1));
      }
      else {
        /* this is a corner tile */
        temp_index = move_strip_down(grid, tile_index - 1);
        assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

        neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
      }
    }

    /* put in the list tiles that are on the row above */
    if (wall_has_initialized_grid(nb_walls[0]) && move_thru_border_0) {
      add_more_tile_neighbors_to_list_fast(
          p,
          wall,
          strip, stripe, flip, wall_vert0, wall_vert1, 0,
          *nb_walls[0],
          neighbors
      );
    }

    /* put in the list tiles that are on the side */
    if ((tile_index == grid.num_tiles - 1) ||
        (tile_index == grid.num_tiles - 2)) {

      if (wall_has_initialized_grid(nb_walls[1]) && move_thru_border_1) {
        add_more_tile_neighbors_to_list_fast(
            p,
            wall,
            strip, stripe, flip, wall_vert1, wall_vert2, 1,
            *nb_walls[1],
            neighbors
        );
      }
    }

  } /* end if ((strip == 0) && (stripe > 0)) */

  if ((strip > 0) && (stripe > 0)) {
    /* We are guaranteed to be in the right layer of the wall.
       Note that no extra checking is done here for that
       assumption besides calling the function "is_inner_tile()"
       before calling this function. */
    if (flip > 0) {
      /* put in the list tiles that are on the same row */
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 1));
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 2));
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 1));

      /* put in the list tiles that are above the current row */
      temp_index = move_strip_up(grid, tile_index);
      assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

      neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 1));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 1));
      /* put in the list tiles that are below the current row */
      temp_index = move_strip_down(grid, tile_index);
      assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

      neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 1));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 2));

    } else { /* (flip == 0) */

      /* put in the list tiles that are on the same row */
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 1));
      neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 2));

      /* put in the list tiles that are above the current row */
      temp_index = move_strip_up(grid, tile_index);
      assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

      neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 1));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index - 2));
      neighbors.push_front(WallTileIndexPair(wall_index, temp_index + 1));
      /* put in the list tiles that are below the current row */
      temp_index = move_strip_down(grid, tile_index - 1);
      assert(temp_index != TILE_INDEX_INVALID && "Error in navigating on the grid");

      neighbors.push_front(WallTileIndexPair(wall_index, temp_index));
    }
    /* put in the list tiles that are on the side */
    if (wall_has_initialized_grid(nb_walls[1]) && move_thru_border_1) {
      add_more_tile_neighbors_to_list_fast(
          p,
          wall,
          strip, stripe, flip, wall_vert1, wall_vert2, 1,
          *nb_walls[1],
          neighbors
      );
    }
  } /* end if ((strip > 0) && (stripe > 0)) */

}



/*****************************************************************************
grid_all_neighbors_for_inner_tile:
  In: a surface grid
      an index on that grid
      a point on that tile
      a linked list of  neighbor tiles (return value)
      a length of the linked list above (return value)
  Out: The list and list length of nearest neighbors
       are returned.  Neighbors should share either common edge or
       common vertice.
  Note: The code below is valid only for the inner tile - the one that has 12
        neighbor tiles all belonging to the same grid as the start tile.
*****************************************************************************/
static void grid_all_neighbors_for_inner_tile(
    Partition& p,
    const Wall& wall,
    tile_index_t tile_index,
    const Vec2& pos,
    TileNeighborVector& neighbors
) {
  const Grid& grid = wall.grid;

  assert(tile_index < grid.num_tiles);

  // Neighboring surface grids (edge-to-edge)
  wall_index_t sw[EDGES_IN_TRIANGLE] = {WALL_INDEX_INVALID, WALL_INDEX_INVALID, WALL_INDEX_INVALID};
  // Indices on those grids (edge-to-edge) of neighbor molecules
  tile_index_t si[EDGES_IN_TRIANGLE] = {TILE_INDEX_INVALID, TILE_INDEX_INVALID, TILE_INDEX_INVALID};

  /* find neighbors to react with */
  grid_neighbors(p, wall, tile_index, sw, si);

  assert(!(wall.index != sw[0] || wall.index != sw[1] || wall.index != sw[2])
      && "Function is called for the tile %d that is not an inner tile."
  );

  int vert_nbr_index = -1;
  for (uint kk = 0; kk < 3; kk++) {
    if ((si[kk] != tile_index - 1) && (si[kk] != tile_index + 1)) {
      vert_nbr_index = si[kk];
      break;
    }
  }

  /* The tile has 2 neighbors to the left and 2 neighbors to the right */
  // original code appends to front/head,
  // OPTIMIZE: can we push the data from the end rather?
  wall_index_t wall_index = wall.index;
  neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 1));
  neighbors.push_front(WallTileIndexPair(wall_index, tile_index - 2));
  neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 1));
  neighbors.push_front(WallTileIndexPair(wall_index, tile_index + 2));

  /* find the orientation of the tile */
  int tile_orient = tile_orientation(pos, wall);

  tile_index_t temp_ind;
  if (tile_orient == 0) {
    /* upright tile has 5 neighbors in the row above it */
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index));
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index - 1));
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index - 2));
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index + 1));
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index + 2));

    /* upright tile has 3 neighbors in the row below it */
    temp_ind = move_strip_down(grid, tile_index);
    assert(temp_ind != TILE_INDEX_INVALID && "Function is called for a tile that is not an inner tile.");

    neighbors.push_front(WallTileIndexPair(wall_index, temp_ind));
    neighbors.push_front(WallTileIndexPair(wall_index, temp_ind - 1));
    neighbors.push_front(WallTileIndexPair(wall_index, temp_ind + 1));
  }
  else {
    /* inverted tile has 3 neighbors in the row above it  */
    temp_ind = move_strip_up(grid, tile_index);
    assert(temp_ind != TILE_INDEX_INVALID && "Function is called for a tile that is not an inner tile.");

    neighbors.push_front(WallTileIndexPair(wall_index, temp_ind));
    neighbors.push_front(WallTileIndexPair(wall_index, temp_ind - 1));
    neighbors.push_front(WallTileIndexPair(wall_index, temp_ind + 1));

    /*   inverted tile has 5 neighbors in the row below it  */
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index));
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index - 1));
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index - 2));
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index + 1));
    neighbors.push_front(WallTileIndexPair(wall_index, vert_nbr_index + 2));
  }

  assert(neighbors.size() == 12 && "Function is called for a tile that is not an inner tile.");
}


/**************************************************************************
find_neighbor_tiles:
  In: a surface molecule
      wall where hit happens, or surface molecule is located
      index of the tile where hit happens, or surface molecule is located
      flag that tells whether we need to create a grid on a neighbor wall
      flag that tells whether we are searching for reactant
          (value = 1) or doing product placement (value = 0)
      a linked list of  neighbor tiles (return value)
      a length of the linked list above (return value)
  Out:
      The list of nearest neighbors
      Neighbors should share either common edge or common vertex.
  Note: This version allows looking for the neighbors at the neighbor walls
       that are connected to the start wall through vertices only.
****************************************************************************/
static void find_neighbor_tiles(
    Partition& p,
    const Molecule* sm, // may be nullptr
    const Wall& wall,
    tile_index_t tile_index,
    bool create_grid_flag,
    bool search_for_reactant,
    TileNeighborVector& neighbors
) {
  const Grid& grid = wall.grid;
  assert(tile_index != TILE_INDEX_INVALID);

  if (is_inner_tile(grid, tile_index)) {

    Vec2 pos = grid2uv(wall, tile_index);
    grid_all_neighbors_for_inner_tile(p, wall, tile_index, pos, neighbors);
  }
  else {
    if (is_corner_tile(grid, tile_index)) {

      /* find tile vertices that are shared with the parent wall */
      vertex_index_t shared_verts[VERTICES_IN_TRIANGLE] = {VERTEX_INDEX_INVALID, VERTEX_INDEX_INVALID, VERTEX_INDEX_INVALID};
      find_shared_vertices_corner_tile_parent_wall(p, wall, grid, tile_index, shared_verts);

      /* create list of neighbor walls that share one vertex
         with the start tile  (not edge-to-edge neighbor walls) */
      WallIndicesVector neighboring_walls;
      WallUtils::find_nbr_walls_shared_one_vertex(p, wall, shared_verts, neighboring_walls);

      if (!neighboring_walls.empty()) {
        grid_all_neighbors_across_walls_through_vertices(
            p, sm, neighboring_walls, wall, create_grid_flag, search_for_reactant, neighbors);
      }

      grid_all_neighbors_across_walls_through_edges(
          p, sm, wall, tile_index, create_grid_flag, search_for_reactant, neighbors);

    }
    else {
      grid_all_neighbors_across_walls_through_edges(
          p, sm, wall, tile_index, create_grid_flag, search_for_reactant, neighbors);
    }
  }

#ifdef DEBUG_GRIDS
  neighbors.dump(__FUNCTION__, "  ");
#endif
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
static tile_index_t nearest_free(
    const Wall& wall, const Vec2& v, const pos_t max_d2,
    pos_t& found_dist2) {

  const Grid& g = wall.grid;

  tile_index_t tile_index;
  pos_t d2;
  const pos_t over3n = 0.333333333333333 / (pos_t)(g.num_tiles_along_axis);

  /* check whether the grid is fully occupied */
  if (g.is_full()) {
    found_dist2 = 0;
    return TILE_INDEX_INVALID;
  }

  tile_index = TILE_INDEX_INVALID;
  d2 = 2 * max_d2 + 1;

  // this seems quite inefficient

  for (uint k = 0; k < g.num_tiles_along_axis; k++) {
    pos_t f = v.v - ((pos_t)(3 * k + 1)) * over3n * wall.uv_vert2.v;
    pos_t ff = f - over3n * wall.uv_vert2.v;
    ff *= ff;
    f *= f;
    if (f > max_d2 && ff > max_d2) {
      continue; /* Entire strip is too far away */
    }

    uint span = (g.num_tiles_along_axis - k);
    for (uint j = 0; j < span; j++) {
      uint can_flip = (j != span - 1);
      for (uint i = 0; i <= can_flip; i++) {
        pos_t fff =
            v.u - over3n * ((pos_t)(3 * j + i + 1) * wall.uv_vert1_u +
                             (pos_t)(3 * k + i + 1) * wall.uv_vert2.u);
        fff *= fff;
        if (i) {
          fff += ff;
        }
        else {
          fff += f;
        }

        if (fff < max_d2 && (tile_index == TILE_INDEX_INVALID || fff < d2)) {
          tile_index_t h = (g.num_tiles_along_axis - k) - 1;
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
*************************************************************************/
static void search_nbhd_for_free(
    Partition& p,
    const wall_index_t origin_wall_index, const Vec2& closest_pos2d, const pos_t max_search_d2,
    wall_index_t& res_wall_index, tile_index_t& res_tile_index,
    set<wall_index_t>& visited_walls,
    int remaining_attempts = 100
  ) {

  release_assert(remaining_attempts > 0
      && "Reached threshold for neighbor search, could not find an empty slot on a wall.");

  wall_index_t best_wall_index = origin_wall_index;
  tile_index_t best_tile_index = TILE_INDEX_INVALID;
  Vec2 best_pos = closest_pos2d;

  const Wall& origin_wall = p.get_wall(origin_wall_index);
  assert(origin_wall.grid.is_initialized());

  // Find index and distance of nearest free grid element on origin wall,
  // returns TILE_INDEX_INVALID when there is not space left
  pos_t d2_unused;
  best_tile_index = GridUtils::nearest_free(origin_wall, best_pos, max_search_d2, d2_unused);

  if (best_tile_index != TILE_INDEX_INVALID) {
    res_wall_index = best_wall_index;
    res_tile_index = best_tile_index;
    return;
  }

  pos_t best_d2 = 2 * max_search_d2 + 1;
  const Vec2& point = closest_pos2d;

  vector<wall_index_t> walls_to_try_if_fail;

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

    if (!origin_wall.is_same_region(there_wall)) {
      continue;
    }

    /* check whether there are any available spots on the neighbor wall */
    if (there_wall.grid.is_initialized()) {
      if (there_wall.grid.is_full()) {
        if (visited_walls.count(there_wall_index) == 0) {
          walls_to_try_if_fail.push_back(there_wall_index);
        }
        continue;
      }
    }

    /* Calculate distance between point and edge j of origin wall */
    Vec2 vurt0, vurt1;
    switch (j) {
      case 0:
        vurt0 = Vec2(0);
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
        vurt1 = Vec2(0);
        break;
      default:
        /* default case should not occur since 0<=j<=2 */
        assert(false);
    }

    Vec2 pt, ed;
    ed = vurt1 - vurt0;
    pt = point - vurt0;

    pos_t d2;
    d2 = dot2(pt, ed);
    d2 = len2_squared(pt) - d2 * d2 / len2_squared(ed); /* Distance squared to line */

    /* Check for free grid element on neighbor if point to edge distance is
     * closer than best_d2  */
    if (d2 < best_d2) {
      if (!there_wall.grid.is_initialized()) {
        there_wall.grid.initialize(p, there_wall);
      }

      GeometryUtils::traverse_surface(origin_wall, point, j, pt);
      tile_index_t i = GridUtils::nearest_free(there_wall, pt, max_search_d2, d2);

      if (i != TILE_INDEX_INVALID && d2 < best_d2) {
        best_tile_index = i;
        best_d2 = d2;
        best_wall_index = there_wall_index;
      }
    }
  }

  if (best_tile_index == TILE_INDEX_INVALID && remaining_attempts != 0) {
    // TODO: choose random direction
    // TODO: it would be better to check the closest walls first
    visited_walls.insert(origin_wall_index);
    for (wall_index_t wi: walls_to_try_if_fail) {
      search_nbhd_for_free(
          p, wi, closest_pos2d, max_search_d2,
          res_wall_index, res_tile_index,
          visited_walls,
          remaining_attempts-1
      );
      if (res_tile_index != TILE_INDEX_INVALID) {
        assert(res_wall_index != WALL_INDEX_INVALID);
        // found a slot
        return;
      }
    }
  }

  res_wall_index = best_wall_index;
  res_tile_index = best_tile_index;
}


static void find_closest_tile_on_wall(
    Partition& p,
    const wall_index_t closest_wall_index, const Vec2& closest_pos2d,
    const pos_t closest_d2, const double search_d2,
    wall_index_t& found_wall_index, tile_index_t& found_tile_index, Vec2& found_pos2d
) {
  tile_index_t closest_tile_index;

  const Wall& w = p.get_wall(closest_wall_index);
  closest_tile_index = GridUtils::uv2grid_tile_index(closest_pos2d, w);

#ifdef DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
  std::cout << "find_closest_tile_on_wall: closest_wall_index: " << closest_wall_index <<
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
    pos_t max_search_d2 = search_d2 - closest_d2;
    if (max_search_d2 <= POS_EPS * POS_EPS) {
      mcell_error("Search distance for find_closest_tile_on_wall is too small");
    }

    wall_index_t new_wall_index;
    wall_index_t new_tile_index;
    set<wall_index_t> visited_walls_ignored;
    search_nbhd_for_free(
        p, closest_wall_index, closest_pos2d, max_search_d2,
        new_wall_index, new_tile_index,
        visited_walls_ignored
    );
    assert(new_wall_index != TILE_INDEX_INVALID);
    assert(new_tile_index != TILE_INDEX_INVALID);

    Vec2 new_pos2d;
    const Wall& w = p.get_wall(new_wall_index);
    if (p.config.randomize_smol_pos) {
      // molecules are processed in different order than in mcell3, this might break
      // compatibility, but the behavior is the same here
      new_pos2d = GridUtils::grid2uv_random(w, new_tile_index, p.aux_rng);
    }
    else {
      new_pos2d = GridUtils::grid2uv(w, new_tile_index);
    }

    found_wall_index = new_wall_index;
    found_tile_index = new_tile_index;
    found_pos2d = new_pos2d;
  }
}


static molecule_id_t place_single_molecule_onto_grid(
    Partition& p,
    rng_state& rng,
    Wall& wall,
    const tile_index_t tile_index,
    const bool override_pos_on_wall,
    const Vec2& pos_on_wall_override_value,
    const species_id_t species_id,
    const orientation_t orientation,
    const double current_time,
    const double release_delay_time

) {

  Vec2 pos_on_wall;
  if (override_pos_on_wall) {
    // MCell3 does not recompute the position when surface molecules are released
    // in the LIST mode/shape
    pos_on_wall = pos_on_wall_override_value;
  }
  else {
    if (p.config.randomize_smol_pos) {
      pos_on_wall = GridUtils::grid2uv_random(wall, tile_index, rng);
    }
    else {
      pos_on_wall = GridUtils::grid2uv(wall, tile_index);
    }
  }

  Molecule sm_to_add(MOLECULE_ID_INVALID, species_id, pos_on_wall, current_time);
  sm_to_add.s.wall_index = wall.index;
  if (orientation == ORIENTATION_NONE) {
    sm_to_add.s.orientation = (rng_uint(&rng) & 1) ? 1 : -1;
  }
  else {
    sm_to_add.s.orientation = orientation;
  }
  sm_to_add.s.grid_tile_index = tile_index;
  sm_to_add.set_flag(MOLECULE_FLAG_SURF);
  sm_to_add.set_flag(MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN);

  // returned molecule has its id set
  Molecule& new_sm = p.add_surface_molecule(
      sm_to_add,
      release_delay_time
  );

  wall.grid.set_molecule_tile(tile_index, new_sm.id);

  return new_sm.id;
}


/*************************************************************************
place_surface_molecule
  In: species for the new molecule
      3D location of the new molecule
      orientation of the new molecule
      diameter to search for a free surface spot
      schedule time for the new molecule
  Out: pointer to the new molecule, or NULL if no free spot was found.
  Note: This function halts the program if it runs out of memory.
        This function is similar to insert_surface_molecule, but it does
        not schedule the molecule or add it to the count.  This is done
        to simplify the logic when placing a surface macromolecule.
        (i.e. place all molecules, and once we're sure we've succeeded,
        schedule them all and count them all.)
 *************************************************************************/
static molecule_id_t place_surface_molecule_to_closest_pos(
    Partition& p,
    rng_state& rng,
    const Vec3& pos,
    const species_id_t species_id,
    const orientation_t orientation,
    const double search_diam,
    const double current_time,
    const double release_delay_time
) {
  int grid_index = 0;
  int *grid_index_p = &grid_index;

  pos_t search_d2;
  if (search_diam <= POS_EPS) {
    search_d2 = POS_EPS * POS_EPS;
  }
  else {
    search_d2 = search_diam * search_diam;
  }

  wall_index_t best_wall_index;
  Vec2 best_wall_pos2d;
  pos_t best_d2 = WallUtils::find_closest_wall_any_object(
      p, pos, search_d2, true, best_wall_index, best_wall_pos2d);

  if (best_wall_index == WALL_INDEX_INVALID) {
    return MOLECULE_ID_INVALID;
  }
  // best wall might be full
  Wall& best_w = p.get_wall(best_wall_index);

  // place molecule onto the found wall, all the remaining information about the molecule stays the same
  wall_index_t found_wall_index;
  tile_index_t found_tile_index;

  if (!best_w.has_initialized_grid()) {
    best_w.initialize_grid(p);
  }

  found_tile_index = GridUtils::uv2grid_tile_index(best_wall_pos2d, best_w);

  if (best_w.grid.get_molecule_on_tile(found_tile_index) != MOLECULE_ID_INVALID) {
    pos_t d2 = search_d2 - best_d2;

    if (d2 <= POS_EPS * POS_EPS) {
      // no free tile close enough was found
      return MOLECULE_ID_INVALID;
    }
    else {
      set<wall_index_t> visited_walls;
      search_nbhd_for_free(
          p, best_w.index, best_wall_pos2d, d2,
          found_wall_index, found_tile_index,
          visited_walls
      );

      if (found_wall_index == WALL_INDEX_INVALID) {
        return MOLECULE_ID_INVALID;
      }
    }
  }
  else {
    found_wall_index = best_wall_index;
  }

  assert(found_wall_index != WALL_INDEX_INVALID);
  assert(found_tile_index != TILE_INDEX_INVALID);

  Wall& found_w = p.get_wall(found_wall_index);
  return place_single_molecule_onto_grid(
      p, rng, found_w, found_tile_index, true, best_wall_pos2d,
      species_id, orientation, current_time, release_delay_time);
}

} // namespace GridUtil
} // namespace MCell

#endif // SRC4_GRID_UTILS_INC_
