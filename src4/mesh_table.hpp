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

#ifndef SRC4_MESH_TABLE_HPP_
#define SRC4_MESH_TABLE_HPP_

// The mesh of one geometry object as its file gives it: the connectivity comes
// only from the vertex indices of the walls (MDL ELEMENT_CONNECTIONS, the API's
// wall_list), never from coordinates. Nothing is merged, reordered or repaired;
// problems are only found and reported.
//
// The MeshTable holds the connectivity as a struct of arrays built in three
// passes over the walls' sides; the checks and the geometry report are built on
// it. Ids are 0-based positions: vertex v is the object's v-th vertex, wall w
// its w-th wall. Entry i = 3w + k is side k of wall w, the walk from the wall's
// vertex k to its vertex (k + 1) % 3, the side across which MCell's
// Wall::nb_walls[k] lies. An edge is a pair of vertices (lo, hi) with lo < hi;
// edges are numbered in the order the entries first use them.
//
// This is MCell's side of coboundary.py: the table gives the sparse matrices d0
// (edges x vertices), d1 (walls x edges) and A20 (walls x vertices) of that code
// as index arrays, with the same numbering. The checks, in the order an object
// is checked:
//   pass 0         indices out of range, a wall that uses a vertex twice, walls
//                  with the same three vertices; for the Python API's lists also
//                  vertices and walls without three values, and surface-region
//                  wall indices out of range
//   zero area      walls whose area MCell's Wall code takes as zero
//   equal coordinates  vertices with exactly the same coordinates (0.0 equals
//                  -0.0); not an error, since connectivity comes from indices
//   1.1 manifold   no edge with three or more walls, and no other non-manifold
//                  vertex (two or more open fans: a bow-tie; a closed fan and
//                  another fan: a pinch)
//   1.2 consistently oriented (if 1.1 holds): the two walls of every edge run
//                  along it in opposite directions
//   1.3 watertight (if 1.1 and 1.2 hold): no edge of only one wall
//   shells         walls joined through edges of two walls; the signed volume of
//                  each watertight shell, and the shells that face the wrong way
//
// Part of the geometry library, which depends on nothing else in MCell, so that
// the test program mesh_utils_test can link it alone. The converters run the
// checks on a model's geometry objects, check_lists for the Python API and
// check_parsed_walls for MDL, and MeshUtils::GeometryReport (mesh_utils.hpp)
// writes the report and the terminal lines.
//
// This header is valid C++17 so that the rest of MCell can include it; the
// implementation is C++23.

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "geometry_analysis.hpp"

namespace MCell {
namespace MeshUtils {

using Point = GeometryAnalysis::Point;
// a wall's three vertex indices, as the model gives them
using WallVertices = GeometryAnalysis::Triangle;

// ids of the vertices, walls, entries and edges of one object
using mesh_index_t = std::uint32_t;
inline constexpr mesh_index_t MESH_INDEX_NONE = std::numeric_limits<mesh_index_t>::max();

// ---------------------------------------------------------------------------
// pass 0 without the table

struct IndexProblems {
  // more vertices, or more than a third as many walls, as mesh_index_t can
  // number; nothing else is checked then
  bool too_large = false;
  // vertices without three coordinates and walls without three vertex indices,
  // ascending; only the Python API's lists can hold them (check_lists)
  std::vector<mesh_index_t> malformed_vertices;
  std::vector<mesh_index_t> malformed_walls;
  // walls with an index below 0 or not below the number of vertices, ascending
  std::vector<mesh_index_t> out_of_range;
  // walls that use one vertex twice, ascending
  std::vector<mesh_index_t> repeated_vertex;

  bool empty() const {
    return !too_large && malformed_vertices.empty() && malformed_walls.empty() && out_of_range.empty() &&
        repeated_vertex.empty();
  }
};

IndexProblems check_indices(std::size_t n_vertices, const std::vector<WallVertices>& walls);

// ---------------------------------------------------------------------------
// the table

struct MeshTable {
  // per entry i = 3w + k
  std::vector<mesh_index_t> edge_of; // the entry's edge
  std::vector<std::int8_t> sign_of;  // +1 when the entry runs from lo to hi, else -1

  // per edge
  std::vector<mesh_index_t> lo, hi;     // its vertices, lo < hi
  std::vector<mesh_index_t> wall_count; // the number of entries on the edge
  std::vector<mesh_index_t> entry0;     // the first of them
  std::vector<mesh_index_t> entry1;     // the second, or MESH_INDEX_NONE

  // per vertex
  std::vector<mesh_index_t> walls_per_vertex; // walls that use the vertex
  std::vector<mesh_index_t> edges_per_vertex; // edges that end at the vertex

  std::size_t n_vertices() const {
    return walls_per_vertex.size();
  }
  std::size_t n_walls() const {
    return edge_of.size() / 3;
  }
  std::size_t n_edges() const {
    return lo.size();
  }

  // the entry on the other side of entry i (wall j / 3, its side j % 3), when
  // i's edge has exactly two walls; MESH_INDEX_NONE otherwise
  mesh_index_t neighbor(const mesh_index_t i) const {
    const mesh_index_t e = edge_of[i];
    if (wall_count[e] != 2) {
      return MESH_INDEX_NONE;
    }
    return entry0[e] == i ? entry1[e] : entry0[e];
  }

  // the vertex where entry i starts: vertex i % 3 of wall i / 3
  mesh_index_t vertex(const mesh_index_t i) const {
    const mesh_index_t e = edge_of[i];
    return sign_of[i] > 0 ? lo[e] : hi[e];
  }
};

// builds the table (passes 0-2); the walls must have passed check_indices
MeshTable build_mesh_table(std::size_t n_vertices, const std::vector<WallVertices>& walls);

// ---------------------------------------------------------------------------
// checks on the table

// groups of walls with the same three vertices, in any order: each group
// ascending, the groups in the order of their first walls
std::vector<std::vector<mesh_index_t>> find_duplicate_walls(
    const MeshTable& table, const std::vector<WallVertices>& walls);

// 1.1; each vertex is listed under the first case that applies to it, so a vertex
// on an edge of three or more walls is not tested for the other two cases. The
// cases are kept apart for the test against coboundary.py only; the report and the
// terminal show the vertices as one list, unclassified (vertices())
struct NonManifold {
  std::vector<mesh_index_t> edges;            // edges with three or more walls, ascending
  std::vector<mesh_index_t> edge_vertices;    // the ends of those edges
  std::vector<mesh_index_t> bowtie_vertices;  // two or more open fans
  std::vector<mesh_index_t> pinched_vertices; // a closed fan and two or more fans in all

  bool empty() const {
    return edges.empty() && bowtie_vertices.empty() && pinched_vertices.empty();
  }

  // every non-manifold vertex once, ascending
  std::vector<mesh_index_t> vertices() const;
};

NonManifold find_non_manifold(const MeshTable& table, const std::vector<WallVertices>& walls);

// the walls of each of the given edges, ascending; the table keeps only the first
// two walls of an edge
std::vector<std::vector<mesh_index_t>> walls_of_edges(
    const MeshTable& table, const std::vector<mesh_index_t>& edges);

// 1.2: edges of two walls that both run along the edge in the same direction,
// ascending
std::vector<mesh_index_t> find_inconsistent_edges(const MeshTable& table);

// 1.3: edges of one wall, ascending
std::vector<mesh_index_t> find_boundary_edges(const MeshTable& table);

// walls joined through edges of exactly two walls; shells are numbered in the
// order of their lowest walls
struct Shells {
  std::vector<mesh_index_t> shell_of; // per wall
  std::size_t count = 0;
};

Shells find_shells(const MeshTable& table);

// the items of one vertex's list
struct IndexRange {
  const mesh_index_t* first = nullptr;
  const mesh_index_t* last = nullptr;

  const mesh_index_t* begin() const {
    return first;
  }
  const mesh_index_t* end() const {
    return last;
  }
  std::size_t size() const {
    return static_cast<std::size_t>(last - first);
  }
  bool empty() const {
    return first == last;
  }
  mesh_index_t operator[](const std::size_t i) const {
    return first[i];
  }
};

// a list per vertex: the items of vertex v are items[begin[v]] up to, not
// including, items[begin[v + 1]]
struct VertexLists {
  std::vector<mesh_index_t> begin;
  std::vector<mesh_index_t> items;

  IndexRange of(const mesh_index_t v) const {
    return {items.data() + begin[v], items.data() + begin[v + 1]};
  }
};

// the walls that use each vertex, ascending
VertexLists walls_of_vertices(const MeshTable& table, const std::vector<WallVertices>& walls);

// the same, with each wall's vertices taken from the table (MeshTable::vertex)
VertexLists walls_of_vertices(const MeshTable& table);

// the 1-ring of each vertex: the other ends of its edges, in the order of the edges
VertexLists neighbors_of_vertices(const MeshTable& table);

// ---------------------------------------------------------------------------
// checks on coordinates

// half the length of (v1 - v0) x (v2 - v0), with the expression of
// Wall::initialize_wall_constants, so that the bits are the same
double wall_area(const Point& v0, const Point& v1, const Point& v2);

// MCell's test for a zero area, !distinguishable_f(area, 0, EPS) with EPS 1e-12,
// as in Wall::initialize_wall_constants: the area is at most 1e-12 (internal
// units squared), or not a number
bool is_zero_area(double area);

// walls whose area is zero, ascending; the points are in internal units
std::vector<mesh_index_t> find_zero_area_walls(
    const std::vector<Point>& points, const std::vector<WallVertices>& walls);

// the orphaned vertices, which no wall of the table uses, ascending. An edge
// can't be orphaned: the table makes each edge from a side of a wall
std::vector<mesh_index_t> find_orphaned_vertices(const MeshTable& table);

// groups of vertices whose coordinates are exactly equal (0.0 equals -0.0): each
// group ascending, the groups in the order of their first vertices
std::vector<std::vector<mesh_index_t>> find_equal_coordinates(const std::vector<Point>& points);

// ---------------------------------------------------------------------------
// all checks of one object, and its block of the geometry report

struct ShellFacts {
  std::size_t walls = 0;      // the number of walls
  mesh_index_t first_wall = 0;
  double area = 0.0;          // internal units squared
  bool watertight = false;    // no edge of the shell has only one wall
  // internal units cubed; computed only when the shell is watertight, since
  // otherwise the sum of the tetrahedra depends on the reference point
  std::optional<double> signed_volume;
  bool facing_judged = false; // see GeometryAnalysis::judge_parts
  bool wrong_facing = false;
};

struct EdgeWalls {
  mesh_index_t lo = 0, hi = 0;
  std::vector<mesh_index_t> walls;
};

// a surface region's wall indices: as the model gives them, or those out of range
struct SurfaceRegionWalls {
  std::string name;
  std::vector<int> walls;
};

struct MeshCheck {
  std::size_t n_vertices = 0;
  std::size_t n_walls = 0; // the walls checked

  // MDL's REMOVE_ELEMENTS: the walls are numbered as in the file, which has
  // more walls than were checked. The file's id of each checked wall, and the
  // file's ids of the walls left out; both empty when no wall was left out
  std::vector<mesh_index_t> file_wall_ids;
  std::vector<mesh_index_t> removed_walls;

  IndexProblems indices;
  // surface regions with wall indices out of range: each such index once, in
  // the order the region gives them; they don't stop the other checks
  std::vector<SurfaceRegionWalls> bad_regions;
  // the rest needs the table, built only if indices.empty()
  bool table_built = false;
  MeshTable table;

  // vertices that no checked wall uses; with MDL's REMOVE_ELEMENTS also those
  // that only removed walls use. An error that doesn't stop the other checks
  std::vector<mesh_index_t> orphaned_vertices;
  std::vector<std::vector<mesh_index_t>> duplicate_walls;
  std::vector<mesh_index_t> zero_area_walls;
  // checked whatever the indices, unless a vertex doesn't have three
  // coordinates; the coordinates of each group, as given
  std::vector<std::vector<mesh_index_t>> equal_coordinates;
  std::vector<Point> equal_coordinates_at;

  NonManifold non_manifold;
  std::vector<EdgeWalls> non_manifold_edges; // the walls of non_manifold.edges

  bool orientation_checked = false;          // 1.1 holds
  std::vector<EdgeWalls> inconsistent_edges;

  bool watertight_checked = false;           // 1.1 and 1.2 hold
  std::size_t boundary_edges = 0;

  // with watertight_checked
  Shells shells;
  std::vector<ShellFacts> shell_facts;
  double area = 0.0;          // of all walls, internal units squared
  // of all walls, internal units cubed; computed only when watertight
  std::optional<double> signed_volume;

  bool manifold() const {
    return table_built && non_manifold.empty();
  }
  bool consistently_oriented() const {
    return orientation_checked && inconsistent_edges.empty();
  }
  bool watertight() const {
    return watertight_checked && boundary_edges == 0;
  }
  std::size_t wrong_facing_shells() const;
  // pass 0 (bad regions included), orphaned vertices, walls with the same three
  // vertices, zero area, 1.1 or 1.2 failed
  bool has_errors() const;
};

// checks one object. points: the vertices in internal units, as MCell simulates
// with them; walls: the vertex indices as given; coordinates: the vertices as the
// model gives them, for the equal-coordinate check and its report (the API's
// vertex_list in microns; may be the same as points)
MeshCheck check_mesh(
    const std::vector<Point>& points, const std::vector<WallVertices>& walls,
    const std::vector<Point>& coordinates);

// checks an object of the Python API from its lists as given: vertex_list in
// microns, wall_list, and the wall indices of its surface regions. Their items
// may hold any number of values; with a vertex that doesn't have three
// coordinates or a wall that doesn't have three vertex indices there is no mesh
// to check, so only pass 0 runs, and equal coordinates when every vertex has
// three. points: vertex_list in internal units, as MCell converts it, read only
// when every vertex has three coordinates
MeshCheck check_lists(
    const std::vector<std::vector<double>>& vertex_list, const std::vector<std::vector<int>>& wall_list,
    const std::vector<Point>& points, const std::vector<SurfaceRegionWalls>& regions);

// checks an MDL object from the walls MCell3 parsed (polygon_object::element).
// removed: the walls MCell3 removed (polygon_object::side_removed), both those of
// REMOVE_ELEMENTS and the zero-area walls it dropped; only the area tells them
// apart. A removed wall with zero area is checked, so that it is reported, and so
// is one with a vertex index out of range, which pass 0 reports; the others are
// taken as REMOVE_ELEMENTS and left out, and the report gives the file's wall ids
MeshCheck check_parsed_walls(
    const std::vector<Point>& points, const std::vector<WallVertices>& file_walls,
    const std::vector<char>& removed, const std::vector<Point>& coordinates);

// writes the object's block of the geometry report, which starts with
// "Object <name>:"; every list of ids gives its count and says what the ids are,
// and walls have the file's ids. Areas and volumes are converted with
// length_unit (microns per internal unit), coordinates are written as given
void write_report_block(
    std::ostream& out, const std::string& name, const MeshCheck& check, double length_unit);

// ---------------------------------------------------------------------------
// the terminal lines

// a problem of one object, as its terminal line gives it
struct ProblemLine {
  bool error = true; // false: a warning, after which the run goes on
  std::string code;  // e.g. MESH_ZERO_AREA_WALLS
  std::string text;  // e.g. "1 wall with zero area"
};

// the object's problems in the order of its report block, errors first. Shells
// that face the wrong way have no line; the report lists them
std::vector<ProblemLine> problem_lines(const MeshCheck& check);

// "Error: [<code>] <object>: <text>; see <report>", or "Warning: ..."
std::string format_problem_line(const ProblemLine& line, const std::string& object, const std::string& report);

// the one stop once every object is checked, e.g. "2 geometry objects have
// errors: Scene.fin, Scene.flipped; see reports/geometry_report_00001.txt"
std::string stop_message(const std::vector<std::string>& objects, const std::string& report);

// the one stop after the usage errors (MeshUtils::UsageErrors), e.g. "2
// geometry objects have usage errors: Scene.box, CP; see
// reports/geometry_report_00001.txt"
std::string usage_stop_message(const std::vector<std::string>& objects, const std::string& report);

} // namespace MeshUtils
} // namespace MCell

#endif // SRC4_MESH_TABLE_HPP_
