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

#ifndef SRC4_GEOMETRY_ANALYSIS_H_
#define SRC4_GEOMETRY_ANALYSIS_H_

// Analysis of closed triangle meshes: cleaning, closedness, volume and area,
// containment of one mesh in another, and OBJ export. It replaces MCell's use
// of VTK. Decisions (is a mesh closed, do two meshes touch, is a point inside)
// are exact: every sign comes from Shewchuk's predicates. Volume and area
// follow vtkMassProperties expression by expression. Like the predicates, the
// decisions are exact as long as no intermediate product underflows or
// overflows, i.e. for coordinate differences roughly between 1e-80 and 1e80
// (or zero), far beyond what MCell models use.
//
// This header is valid C++17 so that the rest of MCell can include it; the
// implementation is C++23.

#include <array>
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

// cleans triangles the way vtkTriangleFilter + vtkCleanPolyData (tolerance 0)
// did: exactly equal points are merged (0.0 equals -0.0), points are numbered
// in the order the triangles first use them, unused points are dropped, and
// triangles that use a point twice are dropped (VTK turned those into lines)
Mesh clean_mesh(const std::vector<Point>& points, const std::vector<Triangle>& triangles);

// every edge belongs to exactly two triangles
// (vtkSelectEnclosedPoints::IsSurfaceClosed)
bool is_closed(const Mesh& mesh);

// no edge belongs to only one triangle; edges with three or more triangles are
// allowed (vtkFeatureEdges, boundary edges only)
bool is_watertight(const Mesh& mesh);

struct VolumeArea {
  double volume;
  double area;
};

// unsigned volume and surface area, computed as vtkPolyDataNormals
// (consistency on, splitting off) followed by vtkMassProperties computes them;
// both are NaN in the case vtkMassProperties calls an "unpredicted situation"
VolumeArea volume_and_area_vtk(const Mesh& mesh);

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

struct NamedMesh {
  std::string name;
  Mesh mesh;
};

// writes <file_prefix>.obj and <file_prefix>.mtl with the same bytes as
// MCell's patched vtkOBJExporter: one named object per mesh, each with a white
// material of opacity 0.25; returns false if a file can't be written
bool write_obj(const std::string& file_prefix, const std::vector<NamedMesh>& meshes);

} // namespace GeometryAnalysis
} // namespace MCell

#endif // SRC4_GEOMETRY_ANALYSIS_H_
