/******************************************************************************
 *
 * Copyright (C) 2006-2017,2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

/* The triangle overlap code in this file is based on the triangle/triangle
 * intersection test routine by by Tomas Moller, 1997.
 *
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * In contrast to Moller's general triangle-triangle intersection routine
 * our code only tests for area overlap of coplanar triangles.
 * In particular it will treat triangles which share an edge or vertex as
 * non-overlapping. Only triangles with a bona-fide area overlap will be
 * flagged. 
 *
 * The code returns 1 if the provided triangles have an area overlap and 0
 * otherwise.
 */

#include "geometry.h"
#include "partition.h"

using namespace std;

namespace MCell {
namespace WallOverlap {

#include "wall_overlap.h"


/******************************************************************
are_walls_coplanar:
  In: first wall
      second wall
      accuracy of the comparison
  Out: 1 if the walls are coplanar
       0 if walls are not coplanar
  Note: see "Real-time rendering" 2nd Ed., by Tomas Akenine-Moller and
        Eric Haines, pp. 590-591
******************************************************************/
static bool are_coplanar(const Partition& p, const Wall& w1, const Wall& w2, const pos_t eps) {

  /* find the plane equation of the second wall in the form (n*x + d2 = 0) */

  pos_t d2, d1_0, d1_1, d1_2;

  const Vec3 w2_vert0 = p.get_wall_vertex(w2, 0);

  const Vec3 w1_vert0 = p.get_wall_vertex(w1, 0);
  const Vec3 w1_vert1 = p.get_wall_vertex(w1, 1);
  const Vec3 w1_vert2 = p.get_wall_vertex(w1, 2);

  d2 = -dot(w2.normal, w2_vert0);

  /* check whether all vertices of the first wall satisfy
     plane equation of the second wall */
  d1_0 = dot(w2.normal, w1_vert0) + d2;
  d1_1 = dot(w2.normal, w1_vert1) + d2;
  d1_2 = dot(w2.normal, w1_vert2) + d2;

  if ((!distinguishable_f(d1_0, 0, eps)) && (!distinguishable_f(d1_1, 0, eps)) &&
      (!distinguishable_f(d1_2, 0, eps))) {
    return true;
  }

  return false;
}


/*******************************************************************
are_walls_coincident:
  In: first wall
      second wall
      accuracy of the comparison
  Out: 0 if the walls are not coincident
       1 if the walls are coincident
*******************************************************************/
static bool are_coincident(const Partition& p, const Wall& w1, const Wall& w2, const pos_t eps) {

  int count = 0;

  for (uint w1_vert = 0; w1_vert < VERTICES_IN_TRIANGLE; w1_vert++) {
    for (uint w2_vert = 0; w2_vert < VERTICES_IN_TRIANGLE; w2_vert++) {
      if (!distinguishable_vec3(
          p.get_wall_vertex(w1, w1_vert), p.get_wall_vertex(w2, w2_vert), eps)) {
        count++;
      }
    }
  }

  if (count >= 3)
    return 1;

  return 0;
}


/* this edge to edge test is based on Franlin Antonio's gem:
 * "Faster Line Segment Intersection", in Graphics Gems III,
 * pp. 199-202 */
static bool edge_edge_test(
    const Vec3& v0, const Vec3& u0, const Vec3& u1,
    uint i0, uint i1,
    pos_t Ax, pos_t Ay) {

  pos_t Bx = u0[i0] - u1[i0];
  pos_t By = u0[i1] - u1[i1];
  pos_t Cx = v0[i0] - u0[i0];
  pos_t Cy = v0[i1] - u0[i1];
  pos_t f = Ay * Bx - Ax * By;
  pos_t d = By * Cx - Bx * Cy;
  if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
    pos_t e = Ax * Cy - Ay * Cx;

    // ignore edge or vertex overlaps
    bool dz = !distinguishable_f(d, 0.0, POS_EPS);
    bool ez = !distinguishable_f(e, 0.0, POS_EPS);
    bool df = !distinguishable_f(d, f, POS_EPS);
    bool ef = !distinguishable_f(e, f, POS_EPS);
    if ((dz && ez) || (dz && ef) || (df && ez) || (df && ef)) {
      return false;
    }

    if (f > 0) {
      if (e >= 0 && e <= f) {
        return true;
      }
    }
    else {
      if (e <= 0 && e >= f) {
        return true;
      }
    }
  }
  return false;
}


/* edge_against_tri_edges tests if edge v0v1 intersects with any edge of
 * triangle u0u1u2 */
static inline bool edge_against_tri_edges(
    const Vec3& v0, const Vec3& v1, const Vec3& u0,
    const Vec3& u1, const Vec3& u2, uint i0,
    uint i1) {

  pos_t Ax = v1[i0] - v0[i0];
  pos_t Ay = v1[i1] - v0[i1];
  /* test edge u0,u1 against v0,v1 */
  if (edge_edge_test(v0, u0, u1, i0, i1, Ax, Ay)) {
    return true;
  }
  /* test edge u1,u2 against v0,v1 */
  if (edge_edge_test(v0, u1, u2, i0, i1, Ax, Ay)) {
    return true;
  };
  /* test edge u2,u1 against v0,v1 */
  if (edge_edge_test(v0, u2, u0, i0, i1, Ax, Ay)) {
    return true;
  }
  return false;
}


/* point_in_tri tests if point v0 is contained in triangle u0u1u2 */
static bool point_in_tri(const Vec3& v0, const Vec3& u0, const Vec3& u1, const Vec3& u2,
                               uint i0, uint i1) {
  /* is T1 completly inside T2? */
  /* check if v0 is inside tri(u0,u1,u2) */
  pos_t a = u1[i1] - u0[i1];
  pos_t b = -(u1[i0] - u0[i0]);
  pos_t c = -a * u0[i0] - b * u0[i1];
  pos_t d0 = a * v0[i0] + b * v0[i1] + c;

  a = u2[i1] - u1[i1];
  b = -(u2[i0] - u1[i0]);
  c = -a * u1[i0] - b * u1[i1];
  pos_t d1 = a * v0[i0] + b * v0[i1] + c;

  a = u0[i1] - u2[i1];
  b = -(u0[i0] - u2[i0]);
  c = -a * u2[i0] - b * u2[i1];
  pos_t d2 = a * v0[i0] + b * v0[i1] + c;
  if (d0 * d1 > 0.0) {
    if (d0 * d2 > 0.0) return true;
  }
  return false;
}


/* coplanar_tri_tri checks if triangles v0v1v2 and u0u1u2 have area overlap. */
static bool coplanar_tri_tri(
    const Vec3& n,
    const Vec3& v0, const Vec3& v1, const Vec3& v2,
    const Vec3& u0, const Vec3& u1, const Vec3& u2) {

  uint i0, i1;
  /* first project onto an axis-aligned plane, that maximizes the area */
  /* of the triangles, compute indices: i0,i1. */
  Vec3 a = abs3(n);
  if (a.x > a.y) {
    if (a.x > a.z) {
      i0 = 1; /* a.x is greatest */
      i1 = 2;
    } else {
      i0 = 0; /* a.z is greatest */
      i1 = 1;
    }
  } else /* a.x<=a.y */
  {
    if (a.z > a.y) {
      i0 = 0; /* a.z is greatest */
      i1 = 1;
    } else {
      i0 = 0; /* a.y is greatest */
      i1 = 2;
    }
  }

  /* test all edges of triangle 1 against the edges of triangle 2 */
  if (edge_against_tri_edges(v0, v1, u0, u1, u2, i0, i1) ||
      edge_against_tri_edges(v1, v2, u0, u1, u2, i0, i1) ||
      edge_against_tri_edges(v2, v0, u0, u1, u2, i0, i1)) {
    return true;
  }

  /* finally, test if tri1 is totally contained in tri2 or vice versa */
  if (point_in_tri(v0, u0, u1, u2, i0, i1) ||
      point_in_tri(u0, v0, v1, v2, i0, i1)) {
    return true;
  }
  return false;
}


/* coplanar_tri_operlap tests if w1 and w2 overlap and returns 1 if yes
 * and 0 otherwise.
 *
 * NOTE 1: w1 and w2 are assumed to be coplanar.
 * NOTE 2: this function will only return 1 if w1 and w2 have a bona
 *         fide area overlap. If w1 and w2 share an edge or vertex
 *         this will *not* be treated as an overlap.
 */
static bool coplanar_walls_overlap(const Partition& p, const Wall& w1, const Wall& w2) {
  assert(are_coplanar(p, w1, w2, POS_EPS));

  const Vec3& v0 = p.get_wall_vertex(w1, 0);
  const Vec3& v1 = p.get_wall_vertex(w1, 1);
  const Vec3& v2 = p.get_wall_vertex(w1, 2);
  const Vec3& u0 = p.get_wall_vertex(w2, 0);
  const Vec3& u1 = p.get_wall_vertex(w2, 1);
  const Vec3& u2 = p.get_wall_vertex(w2, 2);

  /* plane equation 1: N1.X+d1=0 */
  return coplanar_tri_tri(w1.normal, v0, v1, v2, u0, u1, u2);
}


static uint get_num_same_vertex_coords(const Partition& p, const Wall& w1, const Wall& w2) {

  // assuming no vertices for one wall are identical
  uint num_same = 0;
  for (uint i = 0; i < VERTICES_IN_TRIANGLE; i++) {
    const Vec3& v1 = p.get_wall_vertex(w1, i);
    for (uint k = 0; k < VERTICES_IN_TRIANGLE; k++) {
      const Vec3& v2 = p.get_wall_vertex(w2, k);
      if (v1 == v2) {
        num_same++;
      }
    }
  }
  return num_same;
}


bool check_for_overlapped_walls(Partition& p, const Vec3& rand_vec) {

  typedef pair<wall_index_t, double> WallDprodPair;
  vector<WallDprodPair> wall_indices_w_dprod;
  for (const Wall& w: p.get_walls()) {
    double d_prod = dot(rand_vec, w.normal);

    /* we want to place walls with opposite normals into
       neighboring positions in the sorted list */
    if (d_prod < 0) {
      d_prod = -d_prod;
    }

    wall_indices_w_dprod.push_back(make_pair(w.index, d_prod));
  }

  // sort according to dprod
  sort(wall_indices_w_dprod.begin(), wall_indices_w_dprod.end(),
      [](const WallDprodPair& a, const WallDprodPair& b) -> bool {
          return a.second < b.second;
      }
  );


  for (size_t i = 0; i < wall_indices_w_dprod.size(); i++) {
    WallDprodPair& wd = wall_indices_w_dprod[i];
    const Wall& w1 = p.get_wall(wd.first);

    size_t next_index = i + 1;
    while (next_index < wall_indices_w_dprod.size() &&
        (!distinguishable_f(wd.second, wall_indices_w_dprod[next_index].second, EPS))) {

      /* there may be several walls with the same (or mirror) oriented normals */
      const Wall& w2 = p.get_wall(wall_indices_w_dprod[next_index].first);

      const string& obj1_name = p.get_geometry_object(w1.object_id).name;
      const string& obj2_name = p.get_geometry_object(w2.object_id).name;

      uint num_same = get_num_same_vertex_coords(p, w1, w2);
      if (num_same == 3) {
        errs() << "wall overlap: wall side " << w1.side << " from '" << obj1_name <<
            "' overlaps wall side " << w2.side << " from '" << obj2_name <<
            "'.\n";
      }

      if (WallOverlap::are_coplanar(p, w1, w2, MESH_DISTINCTIVE_EPS) &&
           (WallOverlap::are_coincident(p, w1, w2, MESH_DISTINCTIVE_EPS) ||
            WallOverlap::coplanar_walls_overlap(p, w1, w2))
      ) {
        // check for a shared wall
        uint num_same = get_num_same_vertex_coords(p, w1, w2);
        if (num_same == 2) {
          // ok, allowed
          errs() << "edge overlap: wall side " << w1.side << " from '" << obj1_name <<
              "' overlaps wall side " << w2.side << " from '" << obj2_name <<
              "'.\n";
        }
        else if (num_same == 3) {
          // 1) set the primary wall -> first found
          // 2) unify vertex indices

          //merge_wall(w1, w2);
          cout << "wall overlap TODO\n";
          return false;
        }
        else {
          errs() << "walls are overlapped: wall side " << w1.side << " from '" << obj1_name <<
              "' overlaps wall side " << w2.side << " from '" << obj2_name <<
              "', the only overlapping walls that are allowed are those that have the same vertex coordinates.\n";
          return false;
        }
      }

      next_index++;
    }
  }
  return true;
}

} // namespace WallOverlap
} // namespace MCell
