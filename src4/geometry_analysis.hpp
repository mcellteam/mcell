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

#ifndef SRC4_GEOMETRY_ANALYSIS_HPP_
#define SRC4_GEOMETRY_ANALYSIS_HPP_

// Analysis of triangle meshes: volume and area, the parts that face the wrong
// way, containment of one mesh in another, and OBJ export. It replaced MCell's
// use of VTK. Whether a mesh is watertight comes from its
// vertex indices (mesh_table.hpp). Decisions (do two meshes touch, is a point
// inside) are exact: every sign comes from Shewchuk's predicates. Volume and area
// follow MCell's Region and Wall code expression by expression
// (volume_and_area). Like the predicates, the
// decisions are exact as long as no intermediate product underflows or
// overflows, i.e. for coordinate differences roughly between 1e-80 and 1e80
// (or zero), far beyond what MCell models use.
//
// This header is valid C++17 so that the rest of MCell can include it; the
// implementation is C++23.

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace MCell {
namespace GeometryAnalysis {

using Point = std::array<double, 3>;
using Triangle = std::array<int, 3>; // indices into Mesh::points

struct Mesh {
  std::vector<Point> points;
  std::vector<Triangle> triangles;
  Point lo{}, hi{}; // bounding box of the points, zeros for an empty mesh
};

// the triangles as given, with their points numbered in the order the triangles
// first use them; nothing is merged or dropped: points with equal coordinates
// stay apart, and every triangle is kept; unused points are left out
Mesh unmerged_mesh(const std::vector<Point>& points, const std::vector<Triangle>& triangles);

struct InsideOutParts {
  std::size_t components = 0; // connected components that face the wrong way
  double signed_volume = 0.0; // the sum of their signed volumes
};

// the parts of a mesh that face the wrong way. The parts are the connected
// components (triangles joined through edges); a component must face outward
// (a positive signed volume) when an even number of the mesh's other
// components enclose it and inward (a cavity) when an odd number does. Only
// closed, consistently oriented components are judged, and only those whose
// signed volume is large enough for its sign to be certain (a flat component
// is never judged); a component on the surface of another can't be judged
// either
InsideOutParts find_inside_out_parts(const Mesh& mesh);

struct PartFacing {
  double signed_volume = 0.0; // of the part's triangles, for a closed, consistently oriented part
  bool judged = false;        // see judge_parts
  bool wrong = false;         // judged, and facing the wrong way
};

// which way each part of a mesh faces, by the rule of find_inside_out_parts.
// part_of[t] is the part of triangle t, from 0 to closed.size() - 1; closed[p]:
// every edge of part p has exactly two of its triangles; consistent[p]: they
// also run along each edge in opposite directions. A part is judged when it is
// closed and consistently oriented, its signed volume (volume_and_area) is
// large enough for its sign to be certain, and for each other closed part whose
// bounding box contains its own, one of its points lies strictly inside or
// outside that part
std::vector<PartFacing> judge_parts(
    const std::vector<Point>& points, const std::vector<Triangle>& triangles,
    const std::vector<int>& part_of, const std::vector<char>& closed, const std::vector<char>& consistent);

struct VolumeArea {
  double volume;
  double area;
};

// signed volume and surface area as MCell's simulation computes them, and the
// function it uses for the objects' areas and volumes (GeometryObject::get_area,
// get_signed_volume): the volume sums the signed volumes of the tetrahedra that
// join each triangle to the center of the triangles' bounding box, so it is
// negative when the triangles face inward; the area sums 0.5 |(b - a) x (c - a)|
// over the triangles (Wall::initialize_wall_constants); the triangles are used as
// given and in order, without cleaning or reorientation.
// Only for watertight triangles: otherwise the sum of the tetrahedra depends on
// the reference point and is no volume, so it is not computed (surface_area)
VolumeArea volume_and_area(const std::vector<Point>& points, const std::vector<Triangle>& triangles);

// the area of volume_and_area alone, the same terms in the same order, so with
// the same bits; for triangles that aren't watertight
double surface_area(const std::vector<Point>& points, const std::vector<Triangle>& triangles);

// closed triangles a and b share at least one point, exactly
bool triangles_intersect(
    const Point& a0, const Point& a1, const Point& a2,
    const Point& b0, const Point& b1, const Point& b2);

enum class Side { Outside, Inside, OnSurface, Undecided };

enum class Containment { Disjoint, Identical, Intersect, FirstInSecond, SecondInFirst };

class Solid;

// where the point lies with respect to the closed surface; Undecided only if
// every one of many ray directions grazes an edge or a vertex
Side point_in_solid(const Solid& solid, const Point& point);

// relation of two closed solids, with the rules of MCell's former VTK code:
// surfaces that touch or cross are Identical when their point lists are equal
// and Intersect otherwise; else the first solid is inside the second when any
// of its connected components is inside, then the same the other way round;
// else Disjoint
Containment containment(const Solid& first, const Solid& second);

// a mesh prepared for containment tests; the mesh should be closed
class Solid {
public:
  explicit Solid(Mesh mesh);
  ~Solid();
  Solid(Solid&&) noexcept;
  Solid& operator=(Solid&&) noexcept;

  const Mesh& mesh() const;

private:
  struct Impl;
  std::unique_ptr<Impl> pimpl;

  friend Side point_in_solid(const Solid& solid, const Point& point);
  friend Containment containment(const Solid& first, const Solid& second);
};

// shortest text that reads back as the same double, in the notation of
// vtkNumberToString (decimal for exponents -6 to 20, e.g. 1e+21, 1e-07)
std::string number_to_string(double value);

// one object of an OBJ file: its name, its points and its triangles
struct NamedMesh {
  std::string name;
  std::vector<Point> points;
  std::vector<Triangle> triangles; // indices into points
};

// writes <file_prefix>.obj and <file_prefix>.mtl in the layout of MCell's
// patched vtkOBJExporter, which wrote them before: one named object per mesh,
// with its points and triangles in order, each with a white material of
// opacity 0.25; the first line of each file names MCell. Returns false if a
// file can't be written
bool write_obj(const std::string& file_prefix, const std::vector<NamedMesh>& meshes);

} // namespace GeometryAnalysis
} // namespace MCell

#endif // SRC4_GEOMETRY_ANALYSIS_HPP_
