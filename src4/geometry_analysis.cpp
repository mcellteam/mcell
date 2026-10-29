/******************************************************************************
 *
 * Copyright (C) 2026 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

// Built as C++23 with floating-point contraction off (see src4/CMakeLists.txt):
// the volume formula must give the same bits on every platform.

#include "geometry_analysis.h"

#include <algorithm>
#include <bit>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "predicates.h"

namespace MCell {
namespace GeometryAnalysis {

namespace {

// ---------------------------------------------------------------------------
// vector arithmetic, used for numbers only, never for decisions

Point sub(const Point& a, const Point& b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double dot(const Point& a, const Point& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// ---------------------------------------------------------------------------
// exact signs

int sign(double v) {
  return (v > 0.0) - (v < 0.0);
}

// sign of (a - d) . ((b - d) x (c - d)); positive when d lies below the plane
// through a, b and c, seen with a, b, c counterclockwise from above
int orient3d(const Point& a, const Point& b, const Point& c, const Point& d) {
  return sign(Predicates::orient3d(a.data(), b.data(), c.data(), d.data()));
}

using Point2 = std::array<double, 2>;

// sign of (a - c) x (b - c); positive when a, b and c are counterclockwise
int orient2d(const Point2& a, const Point2& b, const Point2& c) {
  return sign(Predicates::orient2d(a.data(), b.data(), c.data()));
}

// coordinates in the plane perpendicular to the axis, in cyclic order so that
// the 2D orientation has the sign of that component of the 3D normal
Point2 project(const Point& p, int axis) {
  switch (axis) {
    case 0:
      return {p[1], p[2]};
    case 1:
      return {p[2], p[0]};
    default:
      return {p[0], p[1]};
  }
}

// the axis along which triangle abc projects to the largest area, or -1 when
// a, b and c lie on one line; any axis with a nonzero area gives exact results
int projection_axis(const Point& a, const Point& b, const Point& c) {
  int best = -1;
  double best_area = 0.0;
  for (int axis = 0; axis < 3; axis++) {
    const Point2 pa = project(a, axis), pb = project(b, axis), pc = project(c, axis);
    const double area = std::fabs(Predicates::orient2d(pa.data(), pb.data(), pc.data()));
    if (area > best_area) {
      best_area = area;
      best = axis;
    }
  }
  return best;
}

bool collinear(const Point& a, const Point& b, const Point& c) {
  return projection_axis(a, b, c) < 0;
}

// p lies in the plane of the triangle abc (not collinear); is it in the
// closed triangle?
bool coplanar_point_in_triangle(const Point& p, const Point& a, const Point& b, const Point& c) {
  const int axis = projection_axis(a, b, c);
  const Point2 pp = project(p, axis), pa = project(a, axis), pb = project(b, axis), pc = project(c, axis);
  const int o = orient2d(pa, pb, pc);
  const int d1 = orient2d(pa, pb, pp), d2 = orient2d(pb, pc, pp), d3 = orient2d(pc, pa, pp);
  return (d1 == 0 || d1 == o) && (d2 == 0 || d2 == o) && (d3 == 0 || d3 == o);
}

bool point_on_triangle(const Point& p, const Point& a, const Point& b, const Point& c) {
  return orient3d(a, b, c, p) == 0 && coplanar_point_in_triangle(p, a, b, c);
}

// r lies on the line through p and q; is it on the closed segment?
bool collinear_point_on_segment(const Point& p, const Point& q, const Point& r) {
  for (int k = 0; k < 3; k++) {
    if (r[k] < std::min(p[k], q[k]) || r[k] > std::max(p[k], q[k])) {
      return false;
    }
  }
  return true;
}

bool point_on_segment(const Point& p, const Point& q, const Point& r) {
  return collinear(p, q, r) && collinear_point_on_segment(p, q, r);
}

// closed 2D segments pq and rs share a point
bool segments_intersect_2d(const Point2& p, const Point2& q, const Point2& r, const Point2& s) {
  const auto on_segment = [](const Point2& a, const Point2& b, const Point2& c) {
    return std::min(a[0], b[0]) <= c[0] && c[0] <= std::max(a[0], b[0]) &&
           std::min(a[1], b[1]) <= c[1] && c[1] <= std::max(a[1], b[1]);
  };
  const int o1 = orient2d(p, q, r), o2 = orient2d(p, q, s);
  const int o3 = orient2d(r, s, p), o4 = orient2d(r, s, q);
  if (o1 * o2 < 0 && o3 * o4 < 0) {
    return true;
  }
  return (o1 == 0 && on_segment(p, q, r)) || (o2 == 0 && on_segment(p, q, s)) ||
         (o3 == 0 && on_segment(r, s, p)) || (o4 == 0 && on_segment(r, s, q));
}

// closed 3D segments pq and rs share a point; either may have zero length
bool segments_intersect(const Point& p, const Point& q, const Point& r, const Point& s) {
  if (p == q) {
    return point_on_segment(r, s, p);
  }
  if (r == s) {
    return point_on_segment(p, q, r);
  }
  if (orient3d(p, q, r, s) != 0) {
    return false; // not coplanar
  }
  int axis = projection_axis(p, q, r);
  if (axis < 0) {
    axis = projection_axis(p, q, s);
  }
  if (axis < 0) {
    // all four points lie on one line
    return collinear_point_on_segment(p, q, r) || collinear_point_on_segment(p, q, s) ||
           collinear_point_on_segment(r, s, p) || collinear_point_on_segment(r, s, q);
  }
  return segments_intersect_2d(project(p, axis), project(q, axis), project(r, axis), project(s, axis));
}

// the closed segment pq (possibly of zero length) and the closed triangle abc
// (not collinear) share a point
bool segment_triangle_intersect(const Point& p, const Point& q, const Point& a, const Point& b, const Point& c) {
  const int sp = orient3d(a, b, c, p), sq = orient3d(a, b, c, q);
  if (sp * sq > 0) {
    return false;
  }
  if (sp == 0 && sq == 0) {
    // in the plane of the triangle
    if (coplanar_point_in_triangle(p, a, b, c) || coplanar_point_in_triangle(q, a, b, c)) {
      return true;
    }
    const int axis = projection_axis(a, b, c);
    const Point2 pp = project(p, axis), pq = project(q, axis);
    const Point2 pa = project(a, axis), pb = project(b, axis), pc = project(c, axis);
    return segments_intersect_2d(pp, pq, pa, pb) || segments_intersect_2d(pp, pq, pb, pc) ||
           segments_intersect_2d(pp, pq, pc, pa);
  }
  if (sp == 0) {
    return coplanar_point_in_triangle(p, a, b, c);
  }
  if (sq == 0) {
    return coplanar_point_in_triangle(q, a, b, c);
  }
  // pq crosses the plane; the line through p and q meets the closed triangle
  // when it passes on the same side of all three edges
  const int e1 = orient3d(p, q, a, b), e2 = orient3d(p, q, b, c), e3 = orient3d(p, q, c, a);
  return !((e1 > 0 || e2 > 0 || e3 > 0) && (e1 < 0 || e2 < 0 || e3 < 0));
}

// the two extreme points of three collinear points
std::pair<Point, Point> collinear_extremes(const Point& a, const Point& b, const Point& c) {
  int axis = 0;
  while (axis < 2 && a[axis] == b[axis] && b[axis] == c[axis]) {
    axis++;
  }
  std::array<Point, 3> pts{a, b, c};
  std::ranges::sort(pts, {}, [axis](const Point& p) { return p[axis]; });
  return {pts[0], pts[2]};
}

// ---------------------------------------------------------------------------
// overlap of two closed, non-degenerate triangles: Guigue and Devillers, "Fast
// and robust triangle-triangle overlap test using orientation predicates",
// Journal of Graphics Tools 8(1), 2003, with every sign computed exactly

bool gd_vertex_2d(const Point2& p1, const Point2& q1, const Point2& r1,
                  const Point2& p2, const Point2& q2, const Point2& r2) {
  if (orient2d(r2, p2, q1) >= 0) {
    if (orient2d(r2, q2, q1) <= 0) {
      if (orient2d(p1, p2, q1) > 0) {
        return orient2d(p1, q2, q1) <= 0;
      }
      if (orient2d(p1, p2, r1) >= 0) {
        return orient2d(q1, r1, p2) >= 0;
      }
      return false;
    }
    if (orient2d(p1, q2, q1) <= 0) {
      if (orient2d(r2, q2, r1) <= 0) {
        return orient2d(q1, r1, q2) >= 0;
      }
      return false;
    }
    return false;
  }
  if (orient2d(r2, p2, r1) >= 0) {
    if (orient2d(q1, r1, r2) >= 0) {
      return orient2d(p1, p2, r1) >= 0;
    }
    if (orient2d(q1, r1, q2) >= 0) {
      return orient2d(r2, r1, q2) >= 0;
    }
    return false;
  }
  return false;
}

bool gd_edge_2d(const Point2& p1, const Point2& q1, const Point2& r1,
                const Point2& p2, const Point2& r2) {
  if (orient2d(r2, p2, q1) >= 0) {
    if (orient2d(p1, p2, q1) >= 0) {
      return orient2d(p1, q1, r2) >= 0;
    }
    if (orient2d(q1, r1, p2) >= 0) {
      return orient2d(r1, p1, p2) >= 0;
    }
    return false;
  }
  if (orient2d(r2, p2, r1) >= 0) {
    if (orient2d(p1, p2, r1) >= 0) {
      if (orient2d(p1, r1, r2) >= 0) {
        return true;
      }
      return orient2d(q1, r1, r2) >= 0;
    }
    return false;
  }
  return false;
}

// both triangles counterclockwise
bool gd_ccw_2d(const Point2& p1, const Point2& q1, const Point2& r1,
               const Point2& p2, const Point2& q2, const Point2& r2) {
  if (orient2d(p2, q2, p1) >= 0) {
    if (orient2d(q2, r2, p1) >= 0) {
      if (orient2d(r2, p2, p1) >= 0) {
        return true;
      }
      return gd_edge_2d(p1, q1, r1, p2, r2);
    }
    if (orient2d(r2, p2, p1) >= 0) {
      return gd_edge_2d(p1, q1, r1, r2, q2);
    }
    return gd_vertex_2d(p1, q1, r1, p2, q2, r2);
  }
  if (orient2d(q2, r2, p1) >= 0) {
    if (orient2d(r2, p2, p1) >= 0) {
      return gd_edge_2d(p1, q1, r1, q2, p2);
    }
    return gd_vertex_2d(p1, q1, r1, q2, r2, p2);
  }
  return gd_vertex_2d(p1, q1, r1, r2, p2, q2);
}

bool gd_2d(const Point2& p1, const Point2& q1, const Point2& r1,
           const Point2& p2, const Point2& q2, const Point2& r2) {
  if (orient2d(p1, q1, r1) < 0) {
    if (orient2d(p2, q2, r2) < 0) {
      return gd_ccw_2d(p1, r1, q1, p2, r2, q2);
    }
    return gd_ccw_2d(p1, r1, q1, p2, q2, r2);
  }
  if (orient2d(p2, q2, r2) < 0) {
    return gd_ccw_2d(p1, q1, r1, p2, r2, q2);
  }
  return gd_ccw_2d(p1, q1, r1, p2, q2, r2);
}

bool gd_coplanar(const Point& p1, const Point& q1, const Point& r1,
                 const Point& p2, const Point& q2, const Point& r2) {
  // both triangles lie in one plane, which the projection along this axis maps
  // one to one onto the coordinate plane
  const int axis = projection_axis(p1, q1, r1);
  return gd_2d(project(p1, axis), project(q1, axis), project(r1, axis),
               project(p2, axis), project(q2, axis), project(r2, axis));
}

bool gd_check_min_max(const Point& p1, const Point& q1, const Point& r1,
                      const Point& p2, const Point& q2, const Point& r2) {
  if (orient3d(q2, p2, p1, q1) > 0) {
    return false;
  }
  if (orient3d(r2, p2, r1, p1) > 0) {
    return false;
  }
  return true;
}

// triangle 1 is rotated so that p1 is alone on its side of the plane of
// triangle 2; dp2, dq2 and dr2 are the sides of triangle 2's points relative
// to the plane of triangle 1
bool gd_3d(const Point& p1, const Point& q1, const Point& r1,
           const Point& p2, const Point& q2, const Point& r2,
           int dp2, int dq2, int dr2) {
  if (dp2 > 0) {
    if (dq2 > 0) {
      return gd_check_min_max(p1, r1, q1, r2, p2, q2);
    }
    if (dr2 > 0) {
      return gd_check_min_max(p1, r1, q1, q2, r2, p2);
    }
    return gd_check_min_max(p1, q1, r1, p2, q2, r2);
  }
  if (dp2 < 0) {
    if (dq2 < 0) {
      return gd_check_min_max(p1, q1, r1, r2, p2, q2);
    }
    if (dr2 < 0) {
      return gd_check_min_max(p1, q1, r1, q2, r2, p2);
    }
    return gd_check_min_max(p1, r1, q1, p2, q2, r2);
  }
  if (dq2 < 0) {
    if (dr2 >= 0) {
      return gd_check_min_max(p1, r1, q1, q2, r2, p2);
    }
    return gd_check_min_max(p1, q1, r1, p2, q2, r2);
  }
  if (dq2 > 0) {
    if (dr2 > 0) {
      return gd_check_min_max(p1, r1, q1, p2, q2, r2);
    }
    return gd_check_min_max(p1, q1, r1, q2, r2, p2);
  }
  if (dr2 > 0) {
    return gd_check_min_max(p1, q1, r1, r2, p2, q2);
  }
  if (dr2 < 0) {
    return gd_check_min_max(p1, r1, q1, r2, p2, q2);
  }
  return gd_coplanar(p1, q1, r1, p2, q2, r2);
}

bool gd_triangles_intersect(const Point& p1, const Point& q1, const Point& r1,
                            const Point& p2, const Point& q2, const Point& r2) {
  // sides of triangle 1's points relative to the plane of triangle 2
  const int dp1 = orient3d(p1, p2, q2, r2), dq1 = orient3d(q1, p2, q2, r2), dr1 = orient3d(r1, p2, q2, r2);
  if (dp1 * dq1 > 0 && dp1 * dr1 > 0) {
    return false;
  }
  // sides of triangle 2's points relative to the plane of triangle 1
  const int dp2 = orient3d(p2, q1, r1, p1), dq2 = orient3d(q2, q1, r1, p1), dr2 = orient3d(r2, q1, r1, p1);
  if (dp2 * dq2 > 0 && dp2 * dr2 > 0) {
    return false;
  }
  if (dp1 > 0) {
    if (dq1 > 0) {
      return gd_3d(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2);
    }
    if (dr1 > 0) {
      return gd_3d(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2);
    }
    return gd_3d(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2);
  }
  if (dp1 < 0) {
    if (dq1 < 0) {
      return gd_3d(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2);
    }
    if (dr1 < 0) {
      return gd_3d(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2);
    }
    return gd_3d(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2);
  }
  if (dq1 < 0) {
    if (dr1 >= 0) {
      return gd_3d(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2);
    }
    return gd_3d(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2);
  }
  if (dq1 > 0) {
    if (dr1 > 0) {
      return gd_3d(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2);
    }
    return gd_3d(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2);
  }
  if (dr1 > 0) {
    return gd_3d(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2);
  }
  if (dr1 < 0) {
    return gd_3d(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2);
  }
  return gd_coplanar(p1, q1, r1, p2, q2, r2);
}

// ---------------------------------------------------------------------------
// cleaning, edges, orientation

// exact equality of points as a hash key; vtkMergePoints compares with ==, so
// 0.0 and -0.0 are the same point
struct PointKey {
  std::uint64_t x, y, z;
  bool operator==(const PointKey&) const = default;
};

std::uint64_t key_bits(double d) {
  return std::bit_cast<std::uint64_t>(d == 0.0 ? 0.0 : d);
}

struct PointKeyHash {
  std::size_t operator()(const PointKey& k) const {
    std::uint64_t h = k.x * 0x9E3779B97F4A7C15ull;
    h ^= k.y + 0x632BE59BD9B4E019ull + (h << 6) + (h >> 2);
    h ^= k.z + 0x94D049BB133111EBull + (h << 6) + (h >> 2);
    return static_cast<std::size_t>(h ^ (h >> 31));
  }
};

// undirected edges, one entry per use, sorted
std::vector<std::uint64_t> sorted_edges(const Mesh& mesh) {
  std::vector<std::uint64_t> edges;
  edges.reserve(mesh.triangles.size() * 3);
  for (const Triangle& t : mesh.triangles) {
    for (int i = 0; i < 3; i++) {
      const std::uint32_t u = t[i], v = t[(i + 1) % 3];
      const auto [a, b] = std::minmax(u, v);
      edges.push_back((std::uint64_t{a} << 32) | b);
    }
  }
  std::ranges::sort(edges);
  return edges;
}

// true when the predicate holds for the number of uses of every edge
template<typename UsesOk>
bool all_edge_uses(const Mesh& mesh, UsesOk uses_ok) {
  const std::vector<std::uint64_t> edges = sorted_edges(mesh);
  for (std::size_t i = 0; i < edges.size();) {
    std::size_t j = i;
    while (j < edges.size() && edges[j] == edges[i]) {
      j++;
    }
    if (!uses_ok(j - i)) {
      return false;
    }
    i = j;
  }
  return true;
}

// the triangle order that vtkPolyDataNormals (consistency on, no splitting)
// produces before vtkMassProperties sees the mesh: breadth-first waves from
// the lowest unvisited triangle, reversing neighbors whose order disagrees
std::vector<Triangle> orient_consistently(const Mesh& mesh) {
  const int n = static_cast<int>(mesh.triangles.size());
  std::vector<Triangle> out = mesh.triangles;

  // point -> triangles in ascending order, as vtkPolyData::BuildLinks
  std::vector<int> start(mesh.points.size() + 1, 0);
  for (const Triangle& t : mesh.triangles) {
    for (int v : t) {
      start[v + 1]++;
    }
  }
  std::partial_sum(start.begin(), start.end(), start.begin());
  std::vector<int> links(start.back());
  std::vector<int> fill(start.begin(), start.end() - 1);
  for (int c = 0; c < n; c++) {
    for (int v : mesh.triangles[c]) {
      links[fill[v]++] = c;
    }
  }

  std::vector<char> visited(n, 0);
  std::vector<int> wave, next_wave, neighbors;
  for (int seed = 0; seed < n; seed++) {
    if (visited[seed]) {
      continue;
    }
    visited[seed] = 1;
    wave.assign(1, seed);
    while (!wave.empty()) {
      for (int c : wave) {
        const Triangle pts = out[c];
        for (int j = 0; j < 3; j++) {
          const int j1 = (j + 1) % 3;
          // vtkPolyData::GetCellEdgeNeighbors(c, pts[j], pts[j1])
          neighbors.clear();
          for (int k = start[pts[j]]; k < start[pts[j] + 1]; k++) {
            const int cell = links[k];
            if (cell == c) {
              continue;
            }
            for (int l = start[pts[j1]]; l < start[pts[j1] + 1]; l++) {
              if (links[l] == cell) {
                neighbors.push_back(cell);
                break;
              }
            }
          }
          for (int nb : neighbors) {
            if (visited[nb]) {
              continue;
            }
            Triangle& np = out[nb];
            int l = 0;
            while (l < 3 && np[l] != pts[j1]) {
              l++;
            }
            if (np[(l + 1) % 3] != pts[j]) {
              std::ranges::reverse(np);
            }
            visited[nb] = 1;
            next_wave.push_back(nb);
          }
        }
      }
      std::swap(wave, next_wave);
      next_wave.clear();
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// bounding volume hierarchy

struct Box {
  Point lo, hi;
};

bool boxes_overlap(const Box& a, const Box& b) {
  return a.lo[0] <= b.hi[0] && b.lo[0] <= a.hi[0] &&
         a.lo[1] <= b.hi[1] && b.lo[1] <= a.hi[1] &&
         a.lo[2] <= b.hi[2] && b.lo[2] <= a.hi[2];
}

bool box_contains(const Box& b, const Point& p) {
  return b.lo[0] <= p[0] && p[0] <= b.hi[0] &&
         b.lo[1] <= p[1] && p[1] <= b.hi[1] &&
         b.lo[2] <= p[2] && p[2] <= b.hi[2];
}

Box triangle_box(const std::array<Point, 3>& t) {
  Box b;
  for (int k = 0; k < 3; k++) {
    b.lo[k] = std::min({t[0][k], t[1][k], t[2][k]});
    b.hi[k] = std::max({t[0][k], t[1][k], t[2][k]});
  }
  return b;
}

struct Node {
  Box box;
  int first; // inner node: index of the left child, the right one follows; leaf: first triangle
  int count; // triangles in a leaf, 0 for an inner node
  int size;  // triangles in the subtree
};

constexpr int LEAF_SIZE = 4;

// might the segment from p to p + dir touch the box? Never answers no for a
// box that the segment touches: the box is padded far beyond rounding errors
bool segment_may_touch_box(const Point& p, const Point& dir, const Box& box) {
  double tmin = 0.0, tmax = 1.0;
  for (int k = 0; k < 3; k++) {
    const double pad =
        1e-9 * (std::fabs(box.lo[k]) + std::fabs(box.hi[k]) + std::fabs(p[k]) + std::fabs(dir[k])) +
        std::numeric_limits<double>::min();
    const double lo = box.lo[k] - pad, hi = box.hi[k] + pad;
    if (dir[k] == 0.0) {
      if (p[k] < lo || p[k] > hi) {
        return false;
      }
      continue;
    }
    const double t_lo = (lo - p[k]) / dir[k], t_hi = (hi - p[k]) / dir[k];
    tmin = std::max(tmin, std::min(t_lo, t_hi));
    tmax = std::min(tmax, std::max(t_lo, t_hi));
    if (tmin > tmax) {
      return false;
    }
  }
  return true;
}

// fixed pseudo-random unit vectors, the same on every platform
const std::vector<Point>& ray_directions() {
  static const std::vector<Point> dirs = [] {
    constexpr int NUM_DIRECTIONS = 64;
    std::vector<Point> res;
    std::uint64_t s = 0x2545F4914F6CDD1Dull;
    const auto rnd = [&s] {
      s ^= s >> 12;
      s ^= s << 25;
      s ^= s >> 27;
      return static_cast<double>((s * 0x2545F4914F6CDD1Dull) >> 11) * 0x1p-53 * 2.0 - 1.0;
    };
    while (res.size() < NUM_DIRECTIONS) {
      const Point d{rnd(), rnd(), rnd()};
      const double l = std::sqrt(dot(d, d));
      if (l > 0.2 && l <= 1.0) {
        res.push_back({d[0] / l, d[1] / l, d[2] / l});
      }
    }
    return res;
  }();
  return dirs;
}

enum class RayHit { Miss, Cross, Degenerate };

// the segment pq against the triangle abc (not collinear), where p does not lie
// on the surface and q lies outside the mesh's bounding box
RayHit ray_hit(const Point& p, const Point& q, const Point& a, const Point& b, const Point& c) {
  const int sp = orient3d(a, b, c, p), sq = orient3d(a, b, c, q);
  if (sp == 0 && sq == 0) {
    return RayHit::Degenerate; // the segment runs in the plane of the triangle
  }
  if (sp == 0 || sq == 0 || sp == sq) {
    return RayHit::Miss; // p and q aren't on the triangle, so meeting the plane only at one of them misses it
  }
  const int e1 = orient3d(p, q, a, b), e2 = orient3d(p, q, b, c), e3 = orient3d(p, q, c, a);
  if ((e1 > 0 || e2 > 0 || e3 > 0) && (e1 < 0 || e2 < 0 || e3 < 0)) {
    return RayHit::Miss;
  }
  if (e1 == 0 || e2 == 0 || e3 == 0) {
    return RayHit::Degenerate; // through an edge or a vertex
  }
  return RayHit::Cross;
}

} // namespace

// ---------------------------------------------------------------------------
// Solid

struct Solid::Impl {
  Mesh mesh;
  std::vector<Node> nodes;
  std::vector<std::array<Point, 3>> tris; // triangle coordinates in leaf order
  std::vector<char> degenerate;           // per entry of tris: its points lie on one line
  std::vector<std::vector<int>> components; // points of each connected component

  explicit Impl(Mesh m);

  bool touches(const Impl& other) const;
  bool on_surface(const Point& p) const;
  Side classify(const Point& p) const;
  bool any_component_inside(const Impl& other) const;
};

Solid::Impl::Impl(Mesh m)
  : mesh(std::move(m)) {

  const int n = static_cast<int>(mesh.triangles.size());

  // connected components, ordered by their smallest point index
  std::vector<int> parent(mesh.points.size());
  std::iota(parent.begin(), parent.end(), 0);
  const auto find = [&parent](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  for (const Triangle& t : mesh.triangles) {
    for (int i = 1; i < 3; i++) {
      const int a = find(t[0]), b = find(t[i]);
      if (a != b) {
        parent[std::max(a, b)] = std::min(a, b);
      }
    }
  }
  std::vector<int> component_of_root(mesh.points.size(), -1);
  for (int i = 0; i < static_cast<int>(mesh.points.size()); i++) {
    const int r = find(i);
    if (component_of_root[r] < 0) {
      component_of_root[r] = static_cast<int>(components.size());
      components.emplace_back();
    }
    components[component_of_root[r]].push_back(i);
  }

  if (n == 0) {
    return;
  }

  // top-down median split along the widest spread of triangle box centers
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::vector<Box> boxes(n);
  std::vector<Point> centers(n);
  for (int t = 0; t < n; t++) {
    const Triangle& tri = mesh.triangles[t];
    boxes[t] = triangle_box({mesh.points[tri[0]], mesh.points[tri[1]], mesh.points[tri[2]]});
    for (int k = 0; k < 3; k++) {
      centers[t][k] = (boxes[t].lo[k] + boxes[t].hi[k]) * 0.5;
    }
  }

  struct Range {
    int node, begin, end;
  };
  nodes.reserve(2 * (n / LEAF_SIZE) + 2);
  nodes.push_back({});
  std::vector<Range> stack{{0, 0, n}};
  while (!stack.empty()) {
    const Range r = stack.back();
    stack.pop_back();
    Node node;
    node.box = boxes[order[r.begin]];
    Box spread{centers[order[r.begin]], centers[order[r.begin]]};
    for (int i = r.begin + 1; i < r.end; i++) {
      const int t = order[i];
      for (int k = 0; k < 3; k++) {
        node.box.lo[k] = std::min(node.box.lo[k], boxes[t].lo[k]);
        node.box.hi[k] = std::max(node.box.hi[k], boxes[t].hi[k]);
        spread.lo[k] = std::min(spread.lo[k], centers[t][k]);
        spread.hi[k] = std::max(spread.hi[k], centers[t][k]);
      }
    }
    node.size = r.end - r.begin;
    if (node.size <= LEAF_SIZE) {
      node.first = r.begin;
      node.count = node.size;
      nodes[r.node] = node;
      continue;
    }
    int axis = 0;
    for (int k = 1; k < 3; k++) {
      if (spread.hi[k] - spread.lo[k] > spread.hi[axis] - spread.lo[axis]) {
        axis = k;
      }
    }
    const int mid = (r.begin + r.end) / 2;
    std::nth_element(order.begin() + r.begin, order.begin() + mid, order.begin() + r.end,
        [&centers, axis](int a, int b) { return centers[a][axis] < centers[b][axis]; });
    node.first = static_cast<int>(nodes.size());
    node.count = 0;
    nodes[r.node] = node;
    nodes.push_back({});
    nodes.push_back({});
    stack.push_back({node.first, r.begin, mid});
    stack.push_back({node.first + 1, mid, r.end});
  }

  tris.resize(n);
  degenerate.resize(n);
  for (int i = 0; i < n; i++) {
    const Triangle& t = mesh.triangles[order[i]];
    tris[i] = {mesh.points[t[0]], mesh.points[t[1]], mesh.points[t[2]]};
    degenerate[i] = collinear(tris[i][0], tris[i][1], tris[i][2]);
  }
}

bool Solid::Impl::touches(const Impl& other) const {
  if (nodes.empty() || other.nodes.empty()) {
    return false;
  }
  std::vector<std::pair<int, int>> stack{{0, 0}};
  while (!stack.empty()) {
    const auto [ia, ib] = stack.back();
    stack.pop_back();
    const Node& na = nodes[ia];
    const Node& nb = other.nodes[ib];
    if (!boxes_overlap(na.box, nb.box)) {
      continue;
    }
    if (na.count && nb.count) {
      for (int i = na.first; i < na.first + na.count; i++) {
        const auto& ta = tris[i];
        const Box box_a = triangle_box(ta);
        for (int j = nb.first; j < nb.first + nb.count; j++) {
          const auto& tb = other.tris[j];
          if (boxes_overlap(box_a, triangle_box(tb)) &&
              triangles_intersect(ta[0], ta[1], ta[2], tb[0], tb[1], tb[2])) {
            return true;
          }
        }
      }
    }
    else if (nb.count || (!na.count && na.size >= nb.size)) {
      stack.push_back({na.first, ib});
      stack.push_back({na.first + 1, ib});
    }
    else {
      stack.push_back({ia, nb.first});
      stack.push_back({ia, nb.first + 1});
    }
  }
  return false;
}

bool Solid::Impl::on_surface(const Point& p) const {
  // collinear triangles are skipped: in a closed mesh their points lie on the
  // edges of their neighbors
  std::vector<int> stack{0};
  while (!stack.empty()) {
    const Node& node = nodes[stack.back()];
    stack.pop_back();
    if (!box_contains(node.box, p)) {
      continue;
    }
    if (node.count) {
      for (int i = node.first; i < node.first + node.count; i++) {
        if (!degenerate[i] && point_on_triangle(p, tris[i][0], tris[i][1], tris[i][2])) {
          return true;
        }
      }
    }
    else {
      stack.push_back(node.first);
      stack.push_back(node.first + 1);
    }
  }
  return false;
}

Side Solid::Impl::classify(const Point& p) const {
  if (nodes.empty() || !box_contains({mesh.lo, mesh.hi}, p)) {
    return Side::Outside;
  }
  if (on_surface(p)) {
    return Side::OnSurface;
  }

  // parity of the crossings of a segment from p to beyond the bounding box;
  // a segment that grazes an edge, a vertex or a plane is replaced by one in
  // another direction
  const Point diagonal = sub(mesh.hi, mesh.lo);
  const double length = 2.0 * std::sqrt(dot(diagonal, diagonal)) + 1.0;
  std::vector<int> stack;
  for (const Point& d : ray_directions()) {
    const Point q{p[0] + d[0] * length, p[1] + d[1] * length, p[2] + d[2] * length};
    const Point dir = sub(q, p);
    int crossings = 0;
    bool degenerate_hit = false;
    stack.assign(1, 0);
    while (!stack.empty() && !degenerate_hit) {
      const Node& node = nodes[stack.back()];
      stack.pop_back();
      if (!segment_may_touch_box(p, dir, node.box)) {
        continue;
      }
      if (node.count) {
        for (int i = node.first; i < node.first + node.count; i++) {
          if (degenerate[i]) {
            continue;
          }
          const RayHit hit = ray_hit(p, q, tris[i][0], tris[i][1], tris[i][2]);
          if (hit == RayHit::Cross) {
            crossings++;
          }
          else if (hit == RayHit::Degenerate) {
            degenerate_hit = true;
            break;
          }
        }
      }
      else {
        stack.push_back(node.first);
        stack.push_back(node.first + 1);
      }
    }
    if (!degenerate_hit) {
      return (crossings % 2 == 1) ? Side::Inside : Side::Outside;
    }
  }
  return Side::Undecided;
}

// valid when the surfaces don't touch: then each connected component lies
// wholly inside or wholly outside the other solid, and one point decides
bool Solid::Impl::any_component_inside(const Impl& other) const {
  for (const std::vector<int>& component : components) {
    for (int point : component) {
      const Side side = other.classify(mesh.points[point]);
      if (side == Side::Inside) {
        return true;
      }
      if (side == Side::Outside) {
        break;
      }
      // undecided: try another point of the same component
    }
  }
  return false;
}

Solid::Solid(Mesh mesh)
  : pimpl(std::make_unique<Impl>(std::move(mesh))) {
}

Solid::~Solid() = default;
Solid::Solid(Solid&&) noexcept = default;
Solid& Solid::operator=(Solid&&) noexcept = default;

const Mesh& Solid::mesh() const {
  return pimpl->mesh;
}

// ---------------------------------------------------------------------------
// public functions

Mesh clean_mesh(const std::vector<Point>& points, const std::vector<Triangle>& triangles) {
  Mesh res;
  std::vector<int> new_index(points.size(), -1);
  std::unordered_map<PointKey, int, PointKeyHash> unique;
  unique.reserve(points.size() * 2 + 1);
  res.triangles.reserve(triangles.size());

  for (const Triangle& t : triangles) {
    Triangle ids;
    for (int i = 0; i < 3; i++) {
      if (t[i] < 0 || static_cast<std::size_t>(t[i]) >= points.size()) {
        throw std::out_of_range("Triangle point index " + std::to_string(t[i]) + " is out of range.");
      }
      int& id = new_index[t[i]];
      if (id < 0) {
        const Point& p = points[t[i]];
        const auto [it, inserted] = unique.emplace(
            PointKey{key_bits(p[0]), key_bits(p[1]), key_bits(p[2])}, static_cast<int>(res.points.size()));
        if (inserted) {
          res.points.push_back(p);
        }
        id = it->second;
      }
      ids[i] = id;
    }
    if (ids[0] != ids[1] && ids[1] != ids[2] && ids[2] != ids[0]) {
      res.triangles.push_back(ids);
    }
  }

  if (!res.points.empty()) {
    res.lo = res.hi = res.points[0];
    for (const Point& p : res.points) {
      for (int k = 0; k < 3; k++) {
        res.lo[k] = std::min(res.lo[k], p[k]);
        res.hi[k] = std::max(res.hi[k], p[k]);
      }
    }
  }
  return res;
}

bool is_closed(const Mesh& mesh) {
  return all_edge_uses(mesh, [](std::size_t uses) { return uses == 2; });
}

bool is_watertight(const Mesh& mesh) {
  return all_edge_uses(mesh, [](std::size_t uses) { return uses != 1; });
}

VolumeArea volume_and_area_vtk(const Mesh& mesh) {
  // the arithmetic follows vtkMassProperties::RequestData expression by
  // expression, so that the results are bit-identical to VTK's
  const std::vector<Triangle> tris = orient_consistently(mesh);
  const double num_cells = static_cast<double>(tris.size());
  double munc[3] = {0.0, 0.0, 0.0}, vol[3] = {0.0, 0.0, 0.0}, kxyz[3];
  double wxyz = 0.0, wxy = 0.0, wxz = 0.0, wyz = 0.0;
  double surface_area = 0.0;

  for (const Triangle& t : tris) {
    double x[3], y[3], z[3];
    for (int idx = 0; idx < 3; idx++) {
      const Point& p = mesh.points[t[idx]];
      x[idx] = p[0];
      y[idx] = p[1];
      z[idx] = p[2];
    }
    double i[3], j[3], k[3], u[3], absu[3], ii[3], jj[3], kk[3];
    i[0] = (x[1] - x[0]);
    j[0] = (y[1] - y[0]);
    k[0] = (z[1] - z[0]);
    i[1] = (x[2] - x[0]);
    j[1] = (y[2] - y[0]);
    k[1] = (z[2] - z[0]);
    i[2] = (x[2] - x[1]);
    j[2] = (y[2] - y[1]);
    k[2] = (z[2] - z[1]);
    u[0] = (j[0] * k[1] - k[0] * j[1]);
    u[1] = (k[0] * i[1] - i[0] * k[1]);
    u[2] = (i[0] * j[1] - j[0] * i[1]);
    const double length = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    if (length != 0.0) {
      u[0] /= length;
      u[1] /= length;
      u[2] /= length;
    }
    else {
      u[0] = u[1] = u[2] = 0.0;
    }
    absu[0] = std::fabs(u[0]);
    absu[1] = std::fabs(u[1]);
    absu[2] = std::fabs(u[2]);
    if ((absu[0] > absu[1]) && (absu[0] > absu[2])) {
      munc[0]++;
    }
    else if ((absu[1] > absu[0]) && (absu[1] > absu[2])) {
      munc[1]++;
    }
    else if ((absu[2] > absu[0]) && (absu[2] > absu[1])) {
      munc[2]++;
    }
    else if ((absu[0] == absu[1]) && (absu[0] == absu[2])) {
      wxyz++;
    }
    else if ((absu[0] == absu[1]) && (absu[0] > absu[2])) {
      wxy++;
    }
    else if ((absu[0] == absu[2]) && (absu[0] > absu[1])) {
      wxz++;
    }
    else if ((absu[1] == absu[2]) && (absu[0] < absu[2])) {
      wyz++;
    }
    else {
      // vtkMassProperties reports "Unpredicted situation"
      const double nan = std::numeric_limits<double>::quiet_NaN();
      return {nan, nan};
    }
    ii[0] = i[0] * i[0];
    ii[1] = i[1] * i[1];
    ii[2] = i[2] * i[2];
    jj[0] = j[0] * j[0];
    jj[1] = j[1] * j[1];
    jj[2] = j[2] * j[2];
    kk[0] = k[0] * k[0];
    kk[1] = k[1] * k[1];
    kk[2] = k[2] * k[2];
    const double a = std::sqrt(ii[1] + jj[1] + kk[1]);
    const double b = std::sqrt(ii[0] + jj[0] + kk[0]);
    const double c = std::sqrt(ii[2] + jj[2] + kk[2]);
    const double s = 0.5 * (a + b + c);
    const double area = std::sqrt(std::fabs(s * (s - a) * (s - b) * (s - c)));
    surface_area += area;
    const double zavg = (z[0] + z[1] + z[2]) / 3.0;
    const double yavg = (y[0] + y[1] + y[2]) / 3.0;
    const double xavg = (x[0] + x[1] + x[2]) / 3.0;
    vol[2] += (area * u[2] * zavg);
    vol[1] += (area * u[1] * yavg);
    vol[0] += (area * u[0] * xavg);
  }
  kxyz[0] = (munc[0] + (wxyz / 3.0) + ((wxy + wxz) / 2.0)) / num_cells;
  kxyz[1] = (munc[1] + (wxyz / 3.0) + ((wxy + wyz) / 2.0)) / num_cells;
  kxyz[2] = (munc[2] + (wxyz / 3.0) + ((wxz + wyz) / 2.0)) / num_cells;
  return {std::fabs(kxyz[0] * vol[0] + kxyz[1] * vol[1] + kxyz[2] * vol[2]), surface_area};
}

bool triangles_intersect(
    const Point& a0, const Point& a1, const Point& a2,
    const Point& b0, const Point& b1, const Point& b2) {
  const bool a_flat = collinear(a0, a1, a2), b_flat = collinear(b0, b1, b2);
  if (!a_flat && !b_flat) {
    return gd_triangles_intersect(a0, a1, a2, b0, b1, b2);
  }
  // a triangle whose points lie on one line is the segment between its
  // extreme points
  if (a_flat && b_flat) {
    const auto [p, q] = collinear_extremes(a0, a1, a2);
    const auto [r, s] = collinear_extremes(b0, b1, b2);
    return segments_intersect(p, q, r, s);
  }
  if (a_flat) {
    const auto [p, q] = collinear_extremes(a0, a1, a2);
    return segment_triangle_intersect(p, q, b0, b1, b2);
  }
  const auto [p, q] = collinear_extremes(b0, b1, b2);
  return segment_triangle_intersect(p, q, a0, a1, a2);
}

Side point_in_solid(const Solid& solid, const Point& point) {
  return solid.pimpl->classify(point);
}

Containment containment(const Solid& first, const Solid& second) {
  const Solid::Impl& a = *first.pimpl;
  const Solid::Impl& b = *second.pimpl;
  if (a.nodes.empty() || b.nodes.empty() || !boxes_overlap({a.mesh.lo, a.mesh.hi}, {b.mesh.lo, b.mesh.hi})) {
    return Containment::Disjoint;
  }
  if (a.touches(b)) {
    return (a.mesh.points == b.mesh.points) ? Containment::Identical : Containment::Intersect;
  }
  if (a.any_component_inside(b)) {
    return Containment::FirstInSecond;
  }
  if (b.any_component_inside(a)) {
    return Containment::SecondInFirst;
  }
  return Containment::Disjoint;
}

std::string number_to_string(double value) {
  // double-conversion's ToShortest with UNIQUE_ZERO and
  // EMIT_POSITIVE_EXPONENT_SIGN, in decimal notation for exponents from -6 to
  // 20, as vtkNumberToString configures it
  if (std::isnan(value)) {
    return "NaN";
  }
  if (std::isinf(value)) {
    return value < 0 ? "-Infinity" : "Infinity";
  }
  if (value == 0.0) {
    return "0";
  }
  char buf[64];
  const auto r = std::to_chars(buf, buf + sizeof(buf), value, std::chars_format::scientific);
  const char* p = buf;
  std::string res;
  if (*p == '-') {
    res += '-';
    p++;
  }
  std::string digits;
  while (p < r.ptr && *p != 'e') {
    if (*p != '.') {
      digits += *p;
    }
    p++;
  }
  int exponent = 0;
  std::from_chars(p + 1 + (p[1] == '+' ? 1 : 0), r.ptr, exponent);
  const int length = static_cast<int>(digits.size());
  const int decimal_point = exponent + 1;
  if (-6 <= exponent && exponent < 21) {
    if (decimal_point <= 0) {
      res += "0.";
      res.append(-decimal_point, '0');
      res += digits;
    }
    else if (decimal_point >= length) {
      res += digits;
      res.append(decimal_point - length, '0');
    }
    else {
      res.append(digits, 0, decimal_point);
      res += '.';
      res.append(digits, decimal_point, std::string::npos);
    }
  }
  else {
    res += digits[0];
    if (length > 1) {
      res += '.';
      res.append(digits, 1, std::string::npos);
    }
    res += 'e';
    res += exponent < 0 ? '-' : '+';
    res += std::to_string(std::abs(exponent));
  }
  return res;
}

bool write_obj(const std::string& file_prefix, const std::vector<NamedMesh>& meshes) {
  if (meshes.empty()) {
    return false; // vtkOBJExporter writes nothing without objects
  }
  // text mode, like vtkOBJExporter
  std::ofstream obj(file_prefix + ".obj");
  std::ofstream mtl(file_prefix + ".mtl");
  if (!obj || !mtl) {
    return false;
  }
  // the material file's name without its directory, as vtksys::SystemTools::GetFilenameName
#ifdef _WIN32
  const char* separators = "/\\";
#else
  const char* separators = "/";
#endif
  const std::string mtl_path = file_prefix + ".mtl";
  const std::size_t slash = mtl_path.find_last_of(separators);
  const std::string mtl_name = (slash == std::string::npos) ? mtl_path : mtl_path.substr(slash + 1);

  obj << "#  wavefront obj file written by the visualization toolkit\n\n";
  obj << "mtllib " << mtl_name << "\n\n";
  mtl << "# wavefront mtl file written by the visualization toolkit\n\n";
  int id_start = 1;
  for (const NamedMesh& m : meshes) {
    mtl << "newmtl mtl_" << m.name << "\n"
        << "Ka 0 0 0\nKd 1 1 1\nKs 0 0 0\nNs 1\nTr 0.25\nd 0.25\nillum 3\n";
    obj << "\no " << m.name << "\n";
    for (const Point& p : m.mesh.points) {
      obj << "v " << number_to_string(p[0]) << " " << number_to_string(p[1]) << " " << number_to_string(p[2]) << "\n";
    }
    obj << "\nusemtl mtl_" << m.name << "\n";
    for (const Triangle& t : m.mesh.triangles) {
      obj << "f " << t[0] + id_start << " " << t[1] + id_start << " " << t[2] + id_start << " \n";
    }
    id_start += static_cast<int>(m.mesh.points.size());
  }
  obj.close();
  mtl.close();
  return obj.good() && mtl.good();
}

} // namespace GeometryAnalysis
} // namespace MCell
