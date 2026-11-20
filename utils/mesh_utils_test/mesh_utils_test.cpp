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

// Built as C++23 (see CMakeLists.txt next to this file).
//
// mesh_utils_test <data dir>
//
// Checks the mesh table and the mesh checks of src4/mesh_table.hpp against the
// expected results stored in <data dir> (mcell_tests/tests/utils_tests/mesh_utils):
//
//   connectivity/<name>.walls, <name>.expected
//       the table and the checks that need only the walls, on the 131 meshes of
//       possible_implementation/check_edge_table_proto.py; the expected results
//       are those of coboundary.py, SciPy and NumPy (see the files' headers)
//   reports/<name>.mesh, <name>.report, <name>.lines
//       the checks of an object as the converters run them (check_lists for the
//       Python API, check_parsed_walls for MDL), its report block
//       (write_report_block) and its terminal lines (problem_lines), on small
//       meshes with coordinates, one per error or warning; the expected blocks
//       and lines are written by hand
//
// The singular and plural form of every terminal line, and the message of the
// one stop, are checked against texts written here.
//
// It prints one line per difference and a summary, and exits with 0 when
// everything agrees, 1 otherwise. It reads only the files given and runs nothing
// else.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "mesh_table.hpp"

namespace fs = std::filesystem;
using namespace MCell;
using namespace MCell::MeshUtils;

namespace {

int n_failures = 0;
int n_checks = 0;

void fail(const std::string& where, const std::string& what) {
  std::cout << "FAIL " << where << ": " << what << '\n';
  n_failures++;
}

template <class A, class B>
void expect_equal(const std::string& where, const std::string& what, const A& got, const B& want) {
  n_checks++;
  if (!(got.size() == want.size() && std::equal(got.begin(), got.end(), want.begin(),
        [](const auto& x, const auto& y) { return static_cast<long long>(x) == static_cast<long long>(y); }))) {
    std::ostringstream s;
    s << what << " differs: got " << got.size() << " values, expected " << want.size();
    for (std::size_t i = 0; i < std::min(got.size(), want.size()); i++) {
      if (static_cast<long long>(got[i]) != static_cast<long long>(want[i])) {
        s << "; first difference at " << i << ": got " << static_cast<long long>(got[i])
          << ", expected " << static_cast<long long>(want[i]);
        break;
      }
    }
    fail(where, s.str());
  }
}

void expect_true(const std::string& where, const std::string& what, const bool ok) {
  n_checks++;
  if (!ok) {
    fail(where, what);
  }
}

// the non-comment lines of a file
std::vector<std::string> data_lines(const fs::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot read " + path.string());
  }
  std::vector<std::string> res;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] != '#') {
      res.push_back(line);
    }
  }
  return res;
}

// "<key> <count>" at lines[at]
std::size_t read_count(const std::vector<std::string>& lines, std::size_t at, const std::string& key) {
  std::istringstream s(lines.at(at));
  std::string k;
  std::size_t n = 0;
  if (!(s >> k >> n) || k != key) {
    throw std::runtime_error("expected '" + key + " <count>' in line: " + lines.at(at));
  }
  return n;
}

struct WallsFile {
  std::size_t n_vertices = 0;
  std::vector<WallVertices> walls;
};

void read_walls(const std::vector<std::string>& lines, std::size_t& at, std::vector<WallVertices>& walls) {
  const std::size_t n = read_count(lines, at++, "walls");
  for (std::size_t w = 0; w < n; w++) {
    std::istringstream s(lines.at(at++));
    WallVertices v;
    if (!(s >> v[0] >> v[1] >> v[2])) {
      throw std::runtime_error("bad wall line: " + lines.at(at - 1));
    }
    walls.push_back(v);
  }
}

WallsFile read_walls_file(const fs::path& path) {
  const std::vector<std::string> lines = data_lines(path);
  WallsFile res;
  std::size_t at = 0;
  res.n_vertices = read_count(lines, at++, "vertices");
  read_walls(lines, at, res.walls);
  return res;
}

// an object as a model gives it: the lines of vertices and walls may hold any
// number of values, as the Python API's lists may; then optional lines
//   region <name> <wall indices>   a surface region (the Python API)
//   removed <wall ids>             the object is MDL's, as MCell3 parsed it, and
//                                  these are the walls MCell3 removed
struct MeshFile {
  double length_unit = 1.0;
  std::vector<std::vector<double>> vertex_list; // in microns
  std::vector<std::vector<int>> wall_list;
  std::vector<SurfaceRegionWalls> regions;
  bool parsed = false; // a "removed" line was given
  std::vector<char> removed;
};

template <class T>
std::vector<T> values_of(const std::string& line) {
  std::istringstream s(line);
  std::vector<T> res;
  T x;
  while (s >> x) {
    res.push_back(x);
  }
  if (!s.eof()) {
    throw std::runtime_error("bad line: " + line);
  }
  return res;
}

MeshFile read_mesh_file(const fs::path& path) {
  const std::vector<std::string> lines = data_lines(path);
  MeshFile res;
  std::size_t at = 0;
  {
    std::istringstream s(lines.at(at++));
    std::string k;
    if (!(s >> k >> res.length_unit) || k != "length_unit") {
      throw std::runtime_error("expected 'length_unit <value>' in " + path.string());
    }
  }
  const std::size_t nv = read_count(lines, at++, "vertices");
  for (std::size_t v = 0; v < nv; v++) {
    res.vertex_list.push_back(values_of<double>(lines.at(at++)));
  }
  const std::size_t nw = read_count(lines, at++, "walls");
  for (std::size_t w = 0; w < nw; w++) {
    res.wall_list.push_back(values_of<int>(lines.at(at++)));
  }
  for (; at < lines.size(); at++) {
    std::istringstream s(lines[at]);
    std::string key;
    s >> key;
    std::string rest;
    std::getline(s, rest);
    if (key == "region") {
      std::istringstream r(rest);
      SurfaceRegionWalls region;
      r >> region.name;
      std::string indices;
      std::getline(r, indices);
      region.walls = values_of<int>(indices);
      res.regions.push_back(region);
    }
    else if (key == "removed") {
      res.parsed = true;
      res.removed.assign(res.wall_list.size(), 0);
      for (const int w : values_of<int>(rest)) {
        res.removed.at(static_cast<std::size_t>(w)) = 1;
      }
    }
    else {
      throw std::runtime_error("unknown line in " + path.string() + ": " + lines[at]);
    }
  }
  return res;
}

// "<key> <values>" lines; sign_of is a string of + and -, duplicate_walls has
// groups separated by |
struct Expected {
  std::map<std::string, std::vector<long long>> values;
  std::vector<int> signs;
  std::vector<std::vector<long long>> duplicate_groups;

  const std::vector<long long>& at(const std::string& key) const {
    const auto it = values.find(key);
    if (it == values.end()) {
      throw std::runtime_error("no " + key + " in the expected results");
    }
    return it->second;
  }
};

Expected read_expected(const fs::path& path) {
  Expected res;
  for (const std::string& line : data_lines(path)) {
    std::istringstream s(line);
    std::string key;
    s >> key;
    if (key == "sign_of") {
      std::string signs;
      s >> signs;
      for (const char c : signs) {
        res.signs.push_back(c == '+' ? 1 : -1);
      }
    }
    else if (key == "duplicate_walls") {
      std::string token;
      std::vector<long long> group;
      while (s >> token) {
        if (token == "|") {
          res.duplicate_groups.push_back(group);
          group.clear();
        }
        else {
          group.push_back(std::stoll(token));
        }
      }
      if (!group.empty()) {
        res.duplicate_groups.push_back(group);
      }
    }
    else {
      std::vector<long long> v;
      long long x;
      while (s >> x) {
        v.push_back(x);
      }
      res.values[key] = v;
    }
  }
  return res;
}

// ---------------------------------------------------------------------------
// connectivity: the table and the checks on the walls, against coboundary.py

void check_connectivity(const std::string& name, const WallsFile& m, const Expected& x) {
  const std::string where = "connectivity/" + name;
  const std::vector<WallVertices>& walls = m.walls;
  const std::size_t n_entries = 3 * walls.size();

  expect_true(where, "the indices are valid", check_indices(m.n_vertices, walls).empty());
  const MeshTable t = build_mesh_table(m.n_vertices, walls);

  // edges: d1 in CSR form (edge_of, sign_of) and the edges' vertex pairs
  const std::size_t n_edges = static_cast<std::size_t>(x.at("edges").at(0));
  expect_true(where, "the number of edges is " + std::to_string(n_edges), t.n_edges() == n_edges);
  std::vector<long long> lo_hi;
  for (std::size_t e = 0; e < t.n_edges(); e++) {
    lo_hi.push_back(t.lo[e]);
    lo_hi.push_back(t.hi[e]);
  }
  expect_equal(where, "the edges' vertices", lo_hi, x.at("edge_vertices"));
  expect_equal(where, "edge_of", t.edge_of, x.at("edge_of"));
  expect_equal(where, "sign_of", t.sign_of, x.signs);

  // what follows from them: d1's column counts and first two rows, and the
  // column counts of A20 and d0
  std::vector<long long> count(n_edges, 0), entry0(n_edges, -1), entry1(n_edges, -1);
  const std::vector<long long>& edge_of = x.at("edge_of");
  for (std::size_t i = 0; i < edge_of.size() && i < n_entries; i++) {
    const std::size_t e = static_cast<std::size_t>(edge_of[i]);
    if (e >= n_edges) {
      continue;
    }
    if (count[e] == 0) {
      entry0[e] = static_cast<long long>(i);
    }
    else if (count[e] == 1) {
      entry1[e] = static_cast<long long>(i);
    }
    count[e]++;
  }
  std::vector<long long> got_entry1;
  for (const mesh_index_t i : t.entry1) {
    got_entry1.push_back(i == MESH_INDEX_NONE ? -1 : static_cast<long long>(i));
  }
  expect_equal(where, "wall_count", t.wall_count, count);
  expect_equal(where, "entry0", t.entry0, entry0);
  expect_equal(where, "entry1", got_entry1, entry1);
  std::vector<long long> walls_per_vertex(m.n_vertices, 0), edges_per_vertex(m.n_vertices, 0);
  for (const WallVertices& w : walls) {
    for (const int v : w) {
      walls_per_vertex[v]++;
    }
  }
  for (const long long v : x.at("edge_vertices")) {
    edges_per_vertex[v]++;
  }
  expect_equal(where, "walls_per_vertex", t.walls_per_vertex, walls_per_vertex);
  expect_equal(where, "edges_per_vertex", t.edges_per_vertex, edges_per_vertex);

  // the neighbor across each side: the other entry of an edge of two walls
  bool neighbors_ok = true;
  for (mesh_index_t i = 0; i < n_entries; i++) {
    const mesh_index_t j = t.neighbor(i);
    if (t.wall_count[t.edge_of[i]] == 2) {
      neighbors_ok = neighbors_ok && j != MESH_INDEX_NONE && j != i && t.neighbor(j) == i &&
          t.edge_of[j] == t.edge_of[i];
    }
    else {
      neighbors_ok = neighbors_ok && j == MESH_INDEX_NONE;
    }
  }
  expect_true(where, "the neighbors across sides are mutual and on the same edge", neighbors_ok);

  // the checks
  const NonManifold nm = find_non_manifold(t, walls);
  expect_equal(where, "non-manifold edges", nm.edges, x.at("non_manifold_edges"));
  expect_equal(where, "vertices of non-manifold edges (case 1)", nm.edge_vertices, x.at("case1_vertices"));
  expect_equal(where, "bow-tie vertices (case 2)", nm.bowtie_vertices, x.at("case2_vertices"));
  expect_equal(where, "pinched vertices (case 3)", nm.pinched_vertices, x.at("case3_vertices"));
  expect_equal(where, "boundary edges", find_boundary_edges(t), x.at("boundary_edges"));
  expect_equal(where, "inconsistent edges", find_inconsistent_edges(t), x.at("inconsistent_edges"));

  const Shells shells = find_shells(t);
  expect_true(where, "the number of shells", static_cast<long long>(shells.count) == x.at("shells").at(0));
  expect_equal(where, "shell_of", shells.shell_of, x.at("shell_of"));

  const std::vector<std::vector<mesh_index_t>> dup = find_duplicate_walls(t, walls);
  bool dup_ok = dup.size() == x.duplicate_groups.size();
  for (std::size_t g = 0; dup_ok && g < dup.size(); g++) {
    dup_ok = std::equal(dup[g].begin(), dup[g].end(), x.duplicate_groups[g].begin(), x.duplicate_groups[g].end(),
                        [](mesh_index_t a, long long b) { return static_cast<long long>(a) == b; });
  }
  expect_true(where, "walls with the same three vertices", dup_ok);

  // the walls of the non-manifold edges: every wall that uses the edge, ascending
  const std::vector<std::vector<mesh_index_t>> on_edges = walls_of_edges(t, nm.edges);
  bool on_edges_ok = on_edges.size() == nm.edges.size();
  for (std::size_t j = 0; on_edges_ok && j < nm.edges.size(); j++) {
    const mesh_index_t e = nm.edges[j];
    std::vector<mesh_index_t> want;
    for (mesh_index_t w = 0; w < walls.size(); w++) {
      const bool has_lo = std::find(walls[w].begin(), walls[w].end(), static_cast<int>(t.lo[e])) != walls[w].end();
      const bool has_hi = std::find(walls[w].begin(), walls[w].end(), static_cast<int>(t.hi[e])) != walls[w].end();
      if (has_lo && has_hi) {
        want.push_back(w);
      }
    }
    on_edges_ok = on_edges[j] == want;
  }
  expect_true(where, "the walls of the non-manifold edges", on_edges_ok);

  // the lists per vertex, by their definition: the walls that use the vertex,
  // ascending (A20's columns); the other ends of the vertex's edges, by edge (d0's)
  const VertexLists vw = walls_of_vertices(t, walls);
  std::vector<std::vector<mesh_index_t>> want_walls(m.n_vertices);
  for (mesh_index_t w = 0; w < walls.size(); w++) {
    for (const int v : walls[w]) {
      want_walls[v].push_back(w);
    }
  }
  const VertexLists vn = neighbors_of_vertices(t);
  std::vector<std::vector<mesh_index_t>> want_neighbors(m.n_vertices);
  const std::vector<long long>& ev = x.at("edge_vertices");
  for (std::size_t e = 0; 2 * e + 1 < ev.size(); e++) {
    want_neighbors[ev[2 * e]].push_back(static_cast<mesh_index_t>(ev[2 * e + 1]));
    want_neighbors[ev[2 * e + 1]].push_back(static_cast<mesh_index_t>(ev[2 * e]));
  }
  bool lists_ok = vw.begin.size() == m.n_vertices + 1 && vn.begin.size() == m.n_vertices + 1;
  for (std::size_t v = 0; lists_ok && v < m.n_vertices; v++) {
    lists_ok = std::vector<mesh_index_t>(vw.items.begin() + vw.begin[v], vw.items.begin() + vw.begin[v + 1]) ==
        want_walls[v] &&
        std::vector<mesh_index_t>(vn.items.begin() + vn.begin[v], vn.items.begin() + vn.begin[v + 1]) ==
        want_neighbors[v];
  }
  expect_true(where, "the walls and the 1-ring of each vertex", lists_ok);

  // the same lists with each wall's vertices taken from the table, which the
  // simulation keeps without the walls; one vertex's list through of()
  bool entries_ok = true;
  for (mesh_index_t w = 0; entries_ok && w < walls.size(); w++) {
    for (mesh_index_t k = 0; k < 3; k++) {
      entries_ok = entries_ok && t.vertex(3 * w + k) == static_cast<mesh_index_t>(walls[w][k]);
    }
  }
  expect_true(where, "the vertex where each entry starts", entries_ok);
  const VertexLists vt = walls_of_vertices(t);
  bool from_table_ok = vt.begin == vw.begin && vt.items == vw.items;
  for (std::size_t v = 0; from_table_ok && v < m.n_vertices; v++) {
    const IndexRange r = vt.of(static_cast<mesh_index_t>(v));
    from_table_ok = std::vector<mesh_index_t>(r.begin(), r.end()) == want_walls[v] && r.size() == want_walls[v].size();
  }
  expect_true(where, "the walls of each vertex from the table alone", from_table_ok);
}

// ---------------------------------------------------------------------------
// reports: check_mesh and write_report_block on meshes with coordinates

std::string read_text(const fs::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot read " + path.string());
  }
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

const char* const REPORT = "reports/geometry_report_00001.txt";

std::string terminal_lines(const MeshCheck& check, const std::string& name) {
  std::string res;
  for (const ProblemLine& line : problem_lines(check)) {
    res += format_problem_line(line, name, REPORT) + '\n';
  }
  return res;
}

void check_report(
    const std::string& name, const MeshFile& m, const std::string& expected_block,
    const std::string& expected_lines) {
  const std::string where = "reports/" + name;

  // internal units, as the API converts vertex_list (vertex * rcp_length_unit)
  const double rcp_length_unit = 1.0 / m.length_unit;
  const bool three_coordinates = std::all_of(m.vertex_list.begin(), m.vertex_list.end(),
      [](const std::vector<double>& v) { return v.size() == 3; });
  std::vector<Point> coordinates, points;
  if (three_coordinates) {
    for (const std::vector<double>& v : m.vertex_list) {
      coordinates.push_back({v[0], v[1], v[2]});
      points.push_back({v[0] * rcp_length_unit, v[1] * rcp_length_unit, v[2] * rcp_length_unit});
    }
  }

  // the checks as the converters run them; walls: the walls checked
  MeshCheck check;
  std::vector<WallVertices> walls;
  if (m.parsed) {
    // MDL: MCell3 parsed every wall as three indices
    std::vector<WallVertices> file_walls;
    for (const std::vector<int>& w : m.wall_list) {
      file_walls.push_back({w.at(0), w.at(1), w.at(2)});
    }
    check = check_parsed_walls(points, file_walls, m.removed, coordinates);
    for (std::size_t w = 0; w < file_walls.size(); w++) {
      if (check.file_wall_ids.empty() ||
          std::find(check.file_wall_ids.begin(), check.file_wall_ids.end(), w) != check.file_wall_ids.end()) {
        walls.push_back(file_walls[w]);
      }
    }
  }
  else {
    check = check_lists(m.vertex_list, m.wall_list, points, m.regions);
    if (check.table_built) {
      for (const std::vector<int>& w : m.wall_list) {
        walls.push_back({w[0], w[1], w[2]});
      }
    }
  }

  std::ostringstream out;
  write_report_block(out, name, check, m.length_unit);
  n_checks++;
  if (out.str() != expected_block) {
    fail(where, "the report block differs; expected:\n" + expected_block + "got:\n" + out.str());
  }
  const std::string lines = terminal_lines(check, name);
  n_checks++;
  if (lines != expected_lines) {
    fail(where, "the terminal lines differ; expected:\n" + expected_lines + "got:\n" + lines);
  }

  if (!check.watertight_checked) {
    return;
  }
  // the object's area is the sum of the walls' areas in order, with the bits of
  // Wall::initialize_wall_constants
  double area = 0.0;
  for (const WallVertices& w : walls) {
    area += wall_area(points[w[0]], points[w[1]], points[w[2]]);
  }
  expect_true(where, "the area is the sum of wall_area over the walls", area == check.area);

  // the shells that face the wrong way are the parts find_inside_out_parts
  // finds in the mesh as given (nothing merged)
  const GeometryAnalysis::InsideOutParts parts =
      GeometryAnalysis::find_inside_out_parts(GeometryAnalysis::unmerged_mesh(points, walls));
  double volume = 0.0;
  for (const ShellFacts& f : check.shell_facts) {
    if (f.wrong_facing) {
      volume += f.signed_volume.value();
    }
  }
  expect_true(where, "find_inside_out_parts agrees on the shells that face the wrong way",
              parts.components == check.wrong_facing_shells() && parts.signed_volume == volume);
}

// ---------------------------------------------------------------------------
// MCell's zero-area test and the wall area

void check_area_functions() {
  const std::string where = "area functions";
  expect_true(where, "zero is a zero area", is_zero_area(0.0));
  expect_true(where, "1e-12 is a zero area", is_zero_area(1e-12));
  expect_true(where, "1.1e-12 is not", !is_zero_area(1.1e-12));
  expect_true(where, "1 is not", !is_zero_area(1.0));
  expect_true(where, "1e300 is not", !is_zero_area(1e300));
  expect_true(where, "NaN is a zero area", is_zero_area(std::numeric_limits<double>::quiet_NaN()));
  expect_true(where, "a right triangle with legs 1 has area 0.5",
              wall_area({0, 0, 0}, {1, 0, 0}, {0, 1, 0}) == 0.5);
  expect_true(where, "collinear points have area 0", wall_area({0, 0, 0}, {1, 1, 1}, {2, 2, 2}) == 0.0);
}

// ---------------------------------------------------------------------------
// the singular and plural form of every terminal line, and the one stop

void expect_lines(const std::string& what, const MeshCheck& check, const std::vector<std::string>& texts) {
  std::string want;
  for (const std::string& t : texts) {
    want += t + "; see " + REPORT + '\n';
  }
  const std::string got = terminal_lines(check, "x");
  n_checks++;
  if (got != want) {
    fail("terminal lines", what + ": expected:\n" + want + "got:\n" + got);
  }
}

void check_line_texts() {
  // pass 0, in the order of the report, with one of each and with two
  MeshCheck one;
  one.indices.malformed_vertices = {2};
  one.indices.malformed_walls = {1};
  one.indices.out_of_range = {3};
  one.indices.repeated_vertex = {4};
  one.bad_regions = {{"r", {9}}};
  expect_lines("pass 0, one of each", one,
      {"Error: [MESH_INVALID_INDICES] x: 1 vertex that doesn't have three coordinates, "
       "1 wall that doesn't have three vertex indices, 1 wall with a vertex index out of range, "
       "1 wall that uses a vertex twice, 1 surface region with a wall index out of range"});

  MeshCheck two;
  two.indices.malformed_vertices = {2, 5};
  two.indices.malformed_walls = {1, 6};
  two.indices.out_of_range = {3, 7};
  two.indices.repeated_vertex = {4, 8};
  two.bad_regions = {{"r", {9}}, {"s", {-1}}};
  expect_lines("pass 0, two of each", two,
      {"Error: [MESH_INVALID_INDICES] x: 2 vertices that don't have three coordinates, "
       "2 walls that don't have three vertex indices, 2 walls with a vertex index out of range, "
       "2 walls that use a vertex twice, 2 surface regions with a wall index out of range"});

  MeshCheck large;
  large.indices.too_large = true;
  expect_lines("too large", large, {"Error: [MESH_INVALID_INDICES] x: too many vertices or walls to check"});

  // the other errors, then the warnings
  MeshCheck c1;
  c1.table_built = true;
  c1.duplicate_walls = {{0, 1}};
  c1.zero_area_walls = {2};
  c1.non_manifold.bowtie_vertices = {0};
  c1.inconsistent_edges = {{0, 1, {0, 1}}};
  c1.equal_coordinates = {{0, 3}};
  c1.watertight_checked = true;
  c1.boundary_edges = 1;
  expect_lines("one of each", c1,
      {"Error: [MESH_DUPLICATE_WALLS] x: 1 group of walls with the same three vertices",
       "Error: [MESH_ZERO_AREA_WALLS] x: 1 wall with zero area",
       "Error: [MESH_NOT_MANIFOLD] x: 1 non-manifold vertex, where groups of walls touch only at a vertex",
       "Error: [MESH_INCONSISTENT_ORIENTATION] x: 1 pair of walls runs along its shared edge in the same direction",
       "Warning: [MESH_DUPLICATE_VERTICES] x: 1 group of vertices with equal coordinates",
       "Warning: [MESH_NOT_WATERTIGHT] x: 1 boundary edge"});

  MeshCheck c2;
  c2.table_built = true;
  c2.duplicate_walls = {{0, 1}, {2, 3}};
  c2.zero_area_walls = {2, 4};
  c2.non_manifold.bowtie_vertices = {0};
  c2.non_manifold.pinched_vertices = {5};
  c2.inconsistent_edges = {{0, 1, {0, 1}}, {1, 2, {1, 2}}, {2, 3, {2, 3}}};
  c2.equal_coordinates = {{0, 3}, {1, 4}};
  c2.watertight_checked = true;
  c2.boundary_edges = 4;
  expect_lines("more of each", c2,
      {"Error: [MESH_DUPLICATE_WALLS] x: 2 groups of walls with the same three vertices",
       "Error: [MESH_ZERO_AREA_WALLS] x: 2 walls with zero area",
       "Error: [MESH_NOT_MANIFOLD] x: 2 non-manifold vertices, where groups of walls touch only at a vertex",
       "Error: [MESH_INCONSISTENT_ORIENTATION] x: 3 pairs of walls run along their shared edge in the same direction",
       "Warning: [MESH_DUPLICATE_VERTICES] x: 2 groups of vertices with equal coordinates",
       "Warning: [MESH_NOT_WATERTIGHT] x: 4 boundary edges"});

  // non-manifold edges: their ends are non-manifold vertices too
  MeshCheck e1;
  e1.table_built = true;
  e1.non_manifold.edges = {0};
  e1.non_manifold.edge_vertices = {0, 1};
  e1.non_manifold_edges = {{0, 1, {0, 1, 2}}};
  expect_lines("one edge", e1,
      {"Error: [MESH_NOT_MANIFOLD] x: 1 edge with 3 or more walls, 2 non-manifold vertices"});
  MeshCheck e2 = e1;
  e2.non_manifold.edges = {0, 4};
  e2.non_manifold.edge_vertices = {0, 1, 2, 3};
  e2.non_manifold_edges = {{0, 1, {0, 1, 2}}, {2, 3, {3, 4, 5}}};
  expect_lines("two edges", e2,
      {"Error: [MESH_NOT_MANIFOLD] x: 2 edges with 3 or more walls, 4 non-manifold vertices"});

  // the stop
  n_checks++;
  if (stop_message({"Scene.fin"}, REPORT) != std::string("1 geometry object has errors: Scene.fin; see ") + REPORT) {
    fail("stop", "one object: got " + stop_message({"Scene.fin"}, REPORT));
  }
  n_checks++;
  if (stop_message({"Scene.fin", "Scene.flipped"}, REPORT) !=
      std::string("2 geometry objects have errors: Scene.fin, Scene.flipped; see ") + REPORT) {
    fail("stop", "two objects: got " + stop_message({"Scene.fin", "Scene.flipped"}, REPORT));
  }
}

std::vector<fs::path> files_with(const fs::path& dir, const std::string& extension) {
  std::vector<fs::path> res;
  for (const fs::directory_entry& e : fs::directory_iterator(dir)) {
    if (e.is_regular_file() && e.path().extension() == extension) {
      res.push_back(e.path());
    }
  }
  std::sort(res.begin(), res.end());
  return res;
}

} // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: mesh_utils_test <data dir>\n";
    return 1;
  }
  const fs::path dir = argv[1];
  try {
    check_area_functions();
    check_line_texts();

    const std::vector<fs::path> walls_files = files_with(dir / "connectivity", ".walls");
    for (const fs::path& p : walls_files) {
      const std::string name = p.stem().string();
      check_connectivity(name, read_walls_file(p), read_expected(p.parent_path() / (name + ".expected")));
    }

    const std::vector<fs::path> mesh_files = files_with(dir / "reports", ".mesh");
    for (const fs::path& p : mesh_files) {
      const std::string name = p.stem().string();
      check_report(name, read_mesh_file(p), read_text(p.parent_path() / (name + ".report")),
                   read_text(p.parent_path() / (name + ".lines")));
    }

    if (walls_files.empty() || mesh_files.empty()) {
      fail(dir.string(), "no test data found");
    }
    std::cout << "mesh_utils_test: " << walls_files.size() << " connectivity meshes, " << mesh_files.size()
              << " reports, " << n_checks << " checks, " << n_failures << " failed\n";
  }
  catch (const std::exception& e) {
    std::cout << "FAIL: " << e.what() << '\n';
    return 1;
  }
  return n_failures == 0 ? 0 : 1;
}
