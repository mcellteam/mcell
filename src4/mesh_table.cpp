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
// wall areas must have the bits of Wall::initialize_wall_constants, and volumes
// those of GeometryAnalysis::volume_and_area.
//
// The table follows possible_implementation/edge_table_cpp.md of the build
// root (passes 0-2 and the checks of its section 5); its expected results come
// from coboundary.py (see mcell_tests/tests/utils_tests/mesh_utils).

#include "mesh_table.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <ostream>
#include <tuple>

namespace MCell {
namespace MeshUtils {

namespace {

// the vertex where entry i starts, and the one where it ends
mesh_index_t start_vertex(const std::vector<WallVertices>& walls, const mesh_index_t i) {
  return static_cast<mesh_index_t>(walls[i / 3][i % 3]);
}

mesh_index_t end_vertex(const std::vector<WallVertices>& walls, const mesh_index_t i) {
  return static_cast<mesh_index_t>(walls[i / 3][(i % 3 + 1) % 3]);
}

// union-find with path halving
mesh_index_t find_root(std::vector<mesh_index_t>& parent, mesh_index_t x) {
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}

// the corner 3w + slot of vertex v in the wall of entry i; v is an end of entry i
mesh_index_t corner(const std::vector<WallVertices>& walls, const mesh_index_t i, const mesh_index_t v) {
  const mesh_index_t first = i - i % 3;
  return start_vertex(walls, i) == v ? i : first + (i % 3 + 1) % 3;
}

// the three vertices of a wall, ascending
std::array<int, 3> sorted_vertices(const WallVertices& w) {
  std::array<int, 3> s{w[0], w[1], w[2]};
  std::sort(s.begin(), s.end());
  return s;
}

// a total order of coordinates in which exactly equal values, and only those,
// are equivalent; -0.0 is taken as 0.0
std::uint64_t coordinate_key(const double x) {
  return std::bit_cast<std::uint64_t>(x == 0.0 ? 0.0 : x);
}

// areas and volumes, to 10 significant digits; -0 is written as 0
std::string format_number(const double x) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.10g", x + 0.0);
  return buf;
}

std::string counted(const std::size_t n, const char* one, const char* many) {
  return std::to_string(n) + " " + (n == 1 ? one : many);
}

void write_ids(std::ostream& out, const std::vector<mesh_index_t>& ids) {
  for (const mesh_index_t id : ids) {
    out << ' ' << id;
  }
}

// a list of ids with its count and what the ids are, e.g. " 2 (walls 3 7)"
void write_counted_ids(
    std::ostream& out, const std::vector<mesh_index_t>& ids, const char* one, const char* many) {
  out << ' ' << ids.size() << " (" << (ids.size() == 1 ? one : many);
  write_ids(out, ids);
  out << ')';
}

} // namespace

// ---------------------------------------------------------------------------
// pass 0 without the table

IndexProblems check_indices(const std::size_t n_vertices, const std::vector<WallVertices>& walls) {
  IndexProblems res;
  if (n_vertices >= MESH_INDEX_NONE || walls.size() >= MESH_INDEX_NONE / 3) {
    res.too_large = true;
    return res;
  }
  for (std::size_t w = 0; w < walls.size(); w++) {
    const WallVertices& v = walls[w];
    const bool in_range = std::all_of(v.begin(), v.end(), [n_vertices](const int i) {
      return i >= 0 && static_cast<std::size_t>(i) < n_vertices;
    });
    if (!in_range) {
      res.out_of_range.push_back(static_cast<mesh_index_t>(w));
    }
    else if (v[0] == v[1] || v[1] == v[2] || v[2] == v[0]) {
      res.repeated_vertex.push_back(static_cast<mesh_index_t>(w));
    }
  }
  return res;
}

// ---------------------------------------------------------------------------
// the table

MeshTable build_mesh_table(const std::size_t n_vertices, const std::vector<WallVertices>& walls) {
  const mesh_index_t n_vert = static_cast<mesh_index_t>(n_vertices);
  const mesh_index_t n = static_cast<mesh_index_t>(3 * walls.size()); // entries

  MeshTable t;
  t.walls_per_vertex.assign(n_vert, 0);
  t.edges_per_vertex.assign(n_vert, 0);
  for (const WallVertices& w : walls) {
    for (const int v : w) {
      t.walls_per_vertex[v]++;
    }
  }

  // pass 1: every entry goes to the bucket of the end of its edge with fewer
  // walls (ties: lo), which keeps buckets short; count the entries per bucket
  std::vector<mesh_index_t> key(n), other(n);
  std::vector<mesh_index_t> begin(static_cast<std::size_t>(n_vert) + 1, 0);
  for (mesh_index_t i = 0; i < n; i++) {
    const mesh_index_t u = start_vertex(walls, i), w = end_vertex(walls, i);
    const mesh_index_t lo = std::min(u, w), hi = std::max(u, w);
    const mesh_index_t k = t.walls_per_vertex[hi] < t.walls_per_vertex[lo] ? hi : lo;
    key[i] = k;
    other[i] = k == lo ? hi : lo;
    begin[k + 1]++;
  }
  for (mesh_index_t k = 0; k < n_vert; k++) {
    begin[k + 1] += begin[k]; // bucket k: slots begin[k] up to begin[k + 1]
  }
  std::vector<mesh_index_t> end(begin.begin(), begin.end() - 1); // every bucket starts empty

  // pass 2: number the edges in the order of first appearance, filling the
  // table on the way
  std::vector<mesh_index_t> slot_other(n), slot_edge(n);
  t.edge_of.resize(n);
  t.sign_of.resize(n);
  for (mesh_index_t i = 0; i < n; i++) {
    const mesh_index_t k = key[i], p = other[i];
    mesh_index_t j = begin[k];
    const mesh_index_t stop = end[k];
    while (j < stop && slot_other[j] != p) {
      j++;
    }
    mesh_index_t e;
    if (j == stop) {
      e = static_cast<mesh_index_t>(t.lo.size());
      slot_other[j] = p;
      slot_edge[j] = e;
      end[k] = stop + 1;
      const mesh_index_t lo = std::min(k, p), hi = std::max(k, p);
      t.lo.push_back(lo);
      t.hi.push_back(hi);
      t.wall_count.push_back(1);
      t.entry0.push_back(i);
      t.entry1.push_back(MESH_INDEX_NONE);
      t.edges_per_vertex[lo]++;
      t.edges_per_vertex[hi]++;
    }
    else {
      e = slot_edge[j];
      if (t.wall_count[e] == 1) {
        t.entry1[e] = i;
      }
      t.wall_count[e]++;
    }
    t.edge_of[i] = e;
    t.sign_of[i] = start_vertex(walls, i) < end_vertex(walls, i) ? 1 : -1;
  }
  return t;
}

// ---------------------------------------------------------------------------
// checks on the table

std::vector<std::vector<mesh_index_t>> walls_of_edges(
    const MeshTable& table, const std::vector<mesh_index_t>& edges) {

  std::vector<std::vector<mesh_index_t>> res(edges.size());
  std::vector<mesh_index_t> position(table.n_edges(), MESH_INDEX_NONE);
  for (std::size_t j = 0; j < edges.size(); j++) {
    position[edges[j]] = static_cast<mesh_index_t>(j);
  }
  // a wall uses three different edges, so each wall comes once, ascending
  for (mesh_index_t i = 0; i < table.edge_of.size(); i++) {
    const mesh_index_t j = position[table.edge_of[i]];
    if (j != MESH_INDEX_NONE) {
      res[j].push_back(i / 3);
    }
  }
  return res;
}

std::vector<std::vector<mesh_index_t>> find_duplicate_walls(
    const MeshTable& table, const std::vector<WallVertices>& walls) {

  // walls with the same three vertices share all their edges, so the walls on
  // a wall's side 0 are the candidates; the table keeps only two walls per edge,
  // the others are listed for the edges with more
  std::vector<mesh_index_t> crowded;
  for (mesh_index_t e = 0; e < table.n_edges(); e++) {
    if (table.wall_count[e] >= 3) {
      crowded.push_back(e);
    }
  }
  const std::vector<std::vector<mesh_index_t>> crowded_walls = walls_of_edges(table, crowded);
  std::vector<mesh_index_t> crowded_index(table.n_edges(), MESH_INDEX_NONE);
  for (std::size_t j = 0; j < crowded.size(); j++) {
    crowded_index[crowded[j]] = static_cast<mesh_index_t>(j);
  }

  // the group of each wall, by its first wall
  std::vector<mesh_index_t> group_of(walls.size(), MESH_INDEX_NONE);
  std::vector<std::vector<mesh_index_t>> groups;
  for (mesh_index_t w = 0; w < walls.size(); w++) {
    const mesh_index_t e = table.edge_of[3 * w];
    const std::array<int, 3> s = sorted_vertices(walls[w]);
    mesh_index_t first = MESH_INDEX_NONE; // the lowest earlier wall with the same vertices
    if (table.wall_count[e] == 2) {
      const mesh_index_t other = table.neighbor(3 * w) / 3;
      if (other < w && sorted_vertices(walls[other]) == s) {
        first = other;
      }
    }
    else if (table.wall_count[e] >= 3) {
      for (const mesh_index_t other : crowded_walls[crowded_index[e]]) {
        if (other >= w) {
          break;
        }
        if (sorted_vertices(walls[other]) == s) {
          first = other;
          break;
        }
      }
    }
    if (first == MESH_INDEX_NONE) {
      continue;
    }
    if (group_of[first] == MESH_INDEX_NONE) {
      group_of[first] = static_cast<mesh_index_t>(groups.size());
      groups.push_back({first});
    }
    group_of[w] = group_of[first];
    groups[group_of[first]].push_back(w);
  }
  std::sort(groups.begin(), groups.end());
  return groups;
}

NonManifold find_non_manifold(const MeshTable& table, const std::vector<WallVertices>& walls) {
  NonManifold res;
  const mesh_index_t n_vert = static_cast<mesh_index_t>(table.n_vertices());

  // case 1: the ends of edges with three or more walls
  std::vector<char> case1(n_vert, 0);
  for (mesh_index_t e = 0; e < table.n_edges(); e++) {
    if (table.wall_count[e] >= 3) {
      res.edges.push_back(e);
      case1[table.lo[e]] = 1;
      case1[table.hi[e]] = 1;
    }
  }

  // closed fans: a union-find over the walls' corners, joined at each edge of
  // two walls; joining two corners that are already connected closes a fan
  std::vector<mesh_index_t> parent(table.edge_of.size());
  std::iota(parent.begin(), parent.end(), 0);
  std::vector<mesh_index_t> closed(n_vert, 0);
  for (mesh_index_t e = 0; e < table.n_edges(); e++) {
    if (table.wall_count[e] != 2) {
      continue;
    }
    for (const mesh_index_t v : {table.lo[e], table.hi[e]}) {
      const mesh_index_t x = find_root(parent, corner(walls, table.entry0[e], v));
      const mesh_index_t y = find_root(parent, corner(walls, table.entry1[e], v));
      if (x == y) {
        closed[v]++;
      }
      else {
        parent[std::max(x, y)] = std::min(x, y);
      }
    }
  }

  for (mesh_index_t v = 0; v < n_vert; v++) {
    if (case1[v]) {
      res.edge_vertices.push_back(v);
      continue;
    }
    // without an edge of three or more walls, the open fans at v are the edges
    // at v minus the walls at v
    const std::int64_t open_fans =
        static_cast<std::int64_t>(table.edges_per_vertex[v]) - static_cast<std::int64_t>(table.walls_per_vertex[v]);
    if (open_fans >= 2) {
      res.bowtie_vertices.push_back(v); // case 2
    }
    else if (closed[v] >= 1 && closed[v] + open_fans >= 2) {
      res.pinched_vertices.push_back(v); // case 3
    }
  }
  return res;
}

std::vector<mesh_index_t> NonManifold::vertices() const {
  // the three lists are disjoint, so no vertex repeats
  std::vector<mesh_index_t> res = edge_vertices;
  res.insert(res.end(), bowtie_vertices.begin(), bowtie_vertices.end());
  res.insert(res.end(), pinched_vertices.begin(), pinched_vertices.end());
  std::sort(res.begin(), res.end());
  return res;
}

std::vector<mesh_index_t> find_inconsistent_edges(const MeshTable& table) {
  std::vector<mesh_index_t> res;
  for (mesh_index_t e = 0; e < table.n_edges(); e++) {
    if (table.wall_count[e] == 2 && table.sign_of[table.entry0[e]] == table.sign_of[table.entry1[e]]) {
      res.push_back(e);
    }
  }
  return res;
}

std::vector<mesh_index_t> find_boundary_edges(const MeshTable& table) {
  std::vector<mesh_index_t> res;
  for (mesh_index_t e = 0; e < table.n_edges(); e++) {
    if (table.wall_count[e] == 1) {
      res.push_back(e);
    }
  }
  return res;
}

Shells find_shells(const MeshTable& table) {
  const mesh_index_t n_walls = static_cast<mesh_index_t>(table.n_walls());
  std::vector<mesh_index_t> parent(n_walls);
  std::iota(parent.begin(), parent.end(), 0);
  for (mesh_index_t e = 0; e < table.n_edges(); e++) {
    if (table.wall_count[e] == 2) {
      const mesh_index_t x = find_root(parent, table.entry0[e] / 3);
      const mesh_index_t y = find_root(parent, table.entry1[e] / 3);
      if (x != y) {
        parent[std::max(x, y)] = std::min(x, y);
      }
    }
  }
  Shells res;
  res.shell_of.resize(n_walls);
  std::vector<mesh_index_t> shell_of_root(n_walls, MESH_INDEX_NONE);
  for (mesh_index_t w = 0; w < n_walls; w++) {
    const mesh_index_t r = find_root(parent, w);
    if (shell_of_root[r] == MESH_INDEX_NONE) {
      shell_of_root[r] = static_cast<mesh_index_t>(res.count++);
    }
    res.shell_of[w] = shell_of_root[r];
  }
  return res;
}

VertexLists walls_of_vertices(const MeshTable& table, const std::vector<WallVertices>& walls) {
  VertexLists res;
  const std::size_t n_vert = table.n_vertices();
  res.begin.assign(n_vert + 1, 0);
  for (std::size_t v = 0; v < n_vert; v++) {
    res.begin[v + 1] = res.begin[v] + table.walls_per_vertex[v];
  }
  res.items.resize(res.begin[n_vert]);
  std::vector<mesh_index_t> next(res.begin.begin(), res.begin.end() - 1);
  for (mesh_index_t w = 0; w < walls.size(); w++) {
    for (const int v : walls[w]) {
      res.items[next[v]++] = w;
    }
  }
  return res;
}

VertexLists neighbors_of_vertices(const MeshTable& table) {
  VertexLists res;
  const std::size_t n_vert = table.n_vertices();
  res.begin.assign(n_vert + 1, 0);
  for (std::size_t v = 0; v < n_vert; v++) {
    res.begin[v + 1] = res.begin[v] + table.edges_per_vertex[v];
  }
  res.items.resize(res.begin[n_vert]);
  std::vector<mesh_index_t> next(res.begin.begin(), res.begin.end() - 1);
  for (mesh_index_t e = 0; e < table.n_edges(); e++) {
    res.items[next[table.lo[e]]++] = table.hi[e];
    res.items[next[table.hi[e]]++] = table.lo[e];
  }
  return res;
}

// ---------------------------------------------------------------------------
// checks on coordinates

double wall_area(const Point& v0, const Point& v1, const Point& v2) {
  // Wall::initialize_wall_constants: vA = v1 - v0, vB = v2 - v0,
  // area = 0.5 * len3(cross(vA, vB)) with glm's cross product
  const double ax = v1[0] - v0[0], ay = v1[1] - v0[1], az = v1[2] - v0[2];
  const double bx = v2[0] - v0[0], by = v2[1] - v0[1], bz = v2[2] - v0[2];
  const double cx = ay * bz - by * az;
  const double cy = az * bx - bz * ax;
  const double cz = ax * by - bx * ay;
  return 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
}

bool is_zero_area(const double area) {
  // distinguishable_f(area, 0, EPS) of libbng/bng/shared_defines.h, negated
  const double c = std::fabs(area - 0.0);
  double a = std::fabs(area);
  if (a < 1) {
    a = 1;
  }
  const double b = 0.0;
  double eps = 1e-12; // EPS
  if (b < a) {
    eps *= a;
  }
  else {
    eps *= b;
  }
  return !(c > eps);
}

std::vector<mesh_index_t> find_zero_area_walls(
    const std::vector<Point>& points, const std::vector<WallVertices>& walls) {

  std::vector<mesh_index_t> res;
  for (mesh_index_t w = 0; w < walls.size(); w++) {
    const WallVertices& v = walls[w];
    if (is_zero_area(wall_area(points[v[0]], points[v[1]], points[v[2]]))) {
      res.push_back(w);
    }
  }
  return res;
}

std::vector<std::vector<mesh_index_t>> find_equal_coordinates(const std::vector<Point>& points) {
  using Key = std::tuple<std::uint64_t, std::uint64_t, std::uint64_t>;
  std::vector<std::pair<Key, mesh_index_t>> keyed;
  keyed.reserve(points.size());
  for (mesh_index_t v = 0; v < points.size(); v++) {
    const Point& p = points[v];
    if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2])) {
      continue; // not equal to anything
    }
    keyed.push_back({Key{coordinate_key(p[0]), coordinate_key(p[1]), coordinate_key(p[2])}, v});
  }
  std::sort(keyed.begin(), keyed.end()); // equal keys by ascending vertex

  std::vector<std::vector<mesh_index_t>> groups;
  for (std::size_t i = 0; i < keyed.size();) {
    std::size_t j = i + 1;
    while (j < keyed.size() && keyed[j].first == keyed[i].first) {
      j++;
    }
    if (j - i >= 2) {
      std::vector<mesh_index_t> group;
      for (std::size_t k = i; k < j; k++) {
        group.push_back(keyed[k].second);
      }
      groups.push_back(std::move(group));
    }
    i = j;
  }
  std::sort(groups.begin(), groups.end());
  return groups;
}

// ---------------------------------------------------------------------------
// all checks of one object, and its block of the geometry report

std::size_t MeshCheck::wrong_facing_shells() const {
  return static_cast<std::size_t>(std::count_if(
      shell_facts.begin(), shell_facts.end(), [](const ShellFacts& s) { return s.wrong_facing; }));
}

bool MeshCheck::has_errors() const {
  return !indices.empty() || !duplicate_walls.empty() || !zero_area_walls.empty() ||
      (table_built && !non_manifold.empty()) || !inconsistent_edges.empty();
}

MeshCheck check_mesh(
    const std::vector<Point>& points, const std::vector<WallVertices>& walls,
    const std::vector<Point>& coordinates) {

  MeshCheck c;
  c.n_vertices = points.size();
  c.n_walls = walls.size();
  c.equal_coordinates = find_equal_coordinates(coordinates);

  // pass 0; without valid indices there is no table
  c.indices = check_indices(points.size(), walls);
  if (!c.indices.empty()) {
    return c;
  }
  c.zero_area_walls = find_zero_area_walls(points, walls);
  c.table = build_mesh_table(points.size(), walls);
  c.table_built = true;
  c.duplicate_walls = find_duplicate_walls(c.table, walls);

  // 1.1
  c.non_manifold = find_non_manifold(c.table, walls);
  if (!c.non_manifold.empty()) {
    const std::vector<std::vector<mesh_index_t>> lists = walls_of_edges(c.table, c.non_manifold.edges);
    for (std::size_t j = 0; j < lists.size(); j++) {
      const mesh_index_t e = c.non_manifold.edges[j];
      c.non_manifold_edges.push_back({c.table.lo[e], c.table.hi[e], lists[j]});
    }
    return c;
  }

  // 1.2
  c.orientation_checked = true;
  for (const mesh_index_t e : find_inconsistent_edges(c.table)) {
    c.inconsistent_edges.push_back(
        {c.table.lo[e], c.table.hi[e], {c.table.entry0[e] / 3, c.table.entry1[e] / 3}});
  }
  if (!c.inconsistent_edges.empty()) {
    return c;
  }

  // 1.3
  c.watertight_checked = true;
  const std::vector<mesh_index_t> boundary = find_boundary_edges(c.table);
  c.boundary_edges = boundary.size();

  // shells
  c.shells = find_shells(c.table);
  const std::size_t n_shells = c.shells.count;
  std::vector<std::vector<WallVertices>> shell_walls(n_shells);
  c.shell_facts.resize(n_shells);
  for (mesh_index_t w = 0; w < walls.size(); w++) {
    const mesh_index_t s = c.shells.shell_of[w];
    if (shell_walls[s].empty()) {
      c.shell_facts[s].first_wall = w;
    }
    shell_walls[s].push_back(walls[w]);
  }
  std::vector<char> closed(n_shells, 1);
  for (const mesh_index_t e : boundary) {
    closed[c.shells.shell_of[c.table.entry0[e] / 3]] = 0;
  }
  std::vector<int> part_of(walls.size());
  for (std::size_t w = 0; w < walls.size(); w++) {
    part_of[w] = static_cast<int>(c.shells.shell_of[w]);
  }
  const std::vector<GeometryAnalysis::PartFacing> facing = GeometryAnalysis::judge_parts(
      points, walls, part_of, closed, std::vector<char>(n_shells, 1));
  for (std::size_t s = 0; s < n_shells; s++) {
    ShellFacts& f = c.shell_facts[s];
    const GeometryAnalysis::VolumeArea va = GeometryAnalysis::volume_and_area(points, shell_walls[s]);
    f.walls = shell_walls[s].size();
    f.area = va.area;
    f.watertight = closed[s];
    f.signed_volume = closed[s] ? va.volume : 0.0;
    f.facing_judged = facing[s].judged;
    f.wrong_facing = facing[s].wrong;
  }

  // the whole object, as Region computes its volume and area
  const GeometryAnalysis::VolumeArea va = GeometryAnalysis::volume_and_area(points, walls);
  c.area = va.area;
  c.signed_volume = c.watertight() ? va.volume : 0.0;
  return c;
}

void write_report_block(
    std::ostream& out, const std::string& name, const MeshCheck& check,
    const std::vector<Point>& coordinates, const double length_unit) {

  const double area_unit = length_unit * length_unit;
  const double volume_unit = length_unit * length_unit * length_unit;

  out << "Object " << name << ": " << counted(check.n_vertices, "vertex", "vertices") << ", "
      << counted(check.n_walls, "wall", "walls");
  if (check.table_built) {
    out << ", " << counted(check.table.n_edges(), "edge", "edges");
  }
  out << '\n';

  // pass 0
  if (check.indices.too_large) {
    out << "  Too many vertices or walls to check\n";
  }
  if (!check.indices.out_of_range.empty()) {
    out << "  Walls with a vertex index out of range:";
    write_counted_ids(out, check.indices.out_of_range, "wall", "walls");
    out << '\n';
  }
  if (!check.indices.repeated_vertex.empty()) {
    out << "  Walls that use a vertex twice:";
    write_counted_ids(out, check.indices.repeated_vertex, "wall", "walls");
    out << '\n';
  }
  if (!check.duplicate_walls.empty()) {
    out << "  Walls with the same three vertices: " << counted(check.duplicate_walls.size(), "group", "groups") << '\n';
    for (const std::vector<mesh_index_t>& g : check.duplicate_walls) {
      out << "    walls";
      write_ids(out, g);
      out << '\n';
    }
  }

  if (!check.zero_area_walls.empty()) {
    out << "  Zero-area walls:";
    write_counted_ids(out, check.zero_area_walls, "wall", "walls");
    out << '\n';
  }

  if (!check.equal_coordinates.empty()) {
    out << "  Vertices with equal coordinates: "
        << counted(check.equal_coordinates.size(), "group", "groups") << '\n';
    for (const std::vector<mesh_index_t>& g : check.equal_coordinates) {
      const Point& p = coordinates[g[0]];
      out << "    vertices";
      write_ids(out, g);
      out << " at (" << GeometryAnalysis::number_to_string(p[0]) << ", "
          << GeometryAnalysis::number_to_string(p[1]) << ", "
          << GeometryAnalysis::number_to_string(p[2]) << ")\n";
    }
  }

  // 1.1
  if (!check.table_built) {
    out << "  Manifold mesh: not checked\n";
  }
  else if (check.non_manifold.empty()) {
    out << "  Manifold mesh: True\n";
  }
  else {
    // the edges' walls, then every non-manifold vertex in one list, unclassified
    out << "  Manifold mesh: False\n";
    if (!check.non_manifold_edges.empty()) {
      out << "    edges with 3 or more walls: " << check.non_manifold_edges.size() << '\n';
      for (const EdgeWalls& e : check.non_manifold_edges) {
        out << "      edge " << e.lo << '-' << e.hi << ": walls";
        write_ids(out, e.walls);
        out << '\n';
      }
    }
    out << "    non-manifold vertices:";
    write_counted_ids(out, check.non_manifold.vertices(), "vertex", "vertices");
    out << '\n';
  }

  // 1.2
  if (!check.orientation_checked) {
    out << "  Consistently oriented: not checked\n";
  }
  else if (check.inconsistent_edges.empty()) {
    out << "  Consistently oriented: True\n";
  }
  else {
    out << "  Consistently oriented: False\n";
    for (const EdgeWalls& e : check.inconsistent_edges) {
      out << "    walls " << e.walls[0] << " and " << e.walls[1] << ", edge " << e.lo << '-' << e.hi << '\n';
    }
  }

  // 1.3
  if (!check.watertight_checked) {
    out << "  Watertight: not checked\n";
    return;
  }
  if (check.boundary_edges == 0) {
    out << "  Watertight: True\n";
  }
  else {
    out << "  Watertight: False (" << counted(check.boundary_edges, "boundary edge", "boundary edges") << ")\n";
  }

  // shells
  out << "  Shells: " << check.shells.count << ", area " << format_number(check.area * area_unit) << " um^2";
  if (check.watertight()) {
    out << ", signed volume " << format_number(check.signed_volume * volume_unit) << " um^3";
  }
  out << '\n';
  const std::size_t wrong = check.wrong_facing_shells();
  if (wrong != 0) {
    out << "  Shells facing the wrong way: " << wrong << '\n';
    for (std::size_t s = 0; s < check.shell_facts.size(); s++) {
      const ShellFacts& f = check.shell_facts[s];
      if (f.wrong_facing) {
        out << "    shell " << s << ": " << counted(f.walls, "wall", "walls") << " from wall " << f.first_wall
            << ", signed volume " << format_number(f.signed_volume * volume_unit) << " um^3\n";
      }
    }
  }
}

} // namespace MeshUtils
} // namespace MCell
