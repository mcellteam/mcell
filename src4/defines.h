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

#ifndef SRC4_DEFINES_H_
#define SRC4_DEFINES_H_

#include <stdint.h>
#include <vector>
#include <set>
#include <string>
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>

#include "mcell_structs.h"
#include "debug_config.h"

// warning: do not use floating point types from directly,
// we need to be able to control the precision
#include "../libs/glm/glm.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include "../libs/glm/gtx/component_wise.hpp"

#include "../libs/bng/defines_shared.h" // may define INDEXER_WA

#if defined(NDEBUG) && defined(INDEXER_WA)
#warning "INDEXER_WA is enabled and this will lead to lower performance"
#endif

#ifndef INDEXER_WA
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>
#endif

// this file must not depend on any other from mcell4 otherwise there
// might be some nasty cyclic include dependencies

/**
 * Usage of IDs and indexes:
 *
 * ID - world-unique identifier of an object
 * Index -partition-unique identifier of an object
 *
 * Why use indices instead of pointers:
 *   - Having a value that is independent on actual position in memory can be beneficial,
 *     this way, we can migrate the data to a different piece of memory and the index will be still valid
 *
 *   - Do we really need to move data? - Yes
 *      - Growing arrays get reallocated
 *      - Defragmentation of molecule arrays
 *
 * Use only IDs?
 *   - This would mean that for every access, we need to use indirection, might be cheap, but also might
 *     block some superscalar executions...
 *   - Also, every type of object would have its own mapping array for all the possible IDs
 *     in each partition
 *
 * What about cases when a wall changes its partition?
 *   - If I hold just its index, I cannot know that it has moved.
 *
 * Go back to single partition? Do not think about multicore execution for now?
 *   - Single partition system won't need to use indices, just IDs
 *   - However, ids do not allow for defragmentation without fixing up all existing references
 *
 * Keep the implementation it as it is for now? - Yes
 *   - molecules use ids, everything else uses indices
 *   - refactor once we will have a good reason to do it - e.g. initially we might now want to allow walls to cross partitions
 *
 * Another approach:
 *   - Id + index for everything that can change?
 *   - Id for everything that does not change?
 *
 */

namespace MCell {

// import shared declarations from BNG into our own namespace
using BNGCommon::float_t;
using BNGCommon::EPS;
using BNGCommon::SQRT_EPS;
using BNGCommon::FLT_GIGANTIC;
using BNGCommon::fabs_f;
using BNGCommon::cmp_eq;
using BNGCommon::distinguishable_f;
using BNGCommon::sqrt_f;
using BNGCommon::pow_f;
using BNGCommon::floor_f;
using BNGCommon::round_f;

// ---------------------------------- optimization macros ----------------------------------
#if defined(likely) or defined(unlikely)
#error "Macros 'likely' or 'unlikely' are already defined"
#endif

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

// ---------------------------------- float types ----------------------------------

const float_t SCHEDULER_COMPARISON_EPS = 1e-10;

// ---------------------------------- configurable constants----------------------------------

#define ENABLE_LEGACY_CALLBACKS

//#define ENABLE_SORT_MOLS_BY_SUBPART // improves cache performance but produces different result
const uint SORT_MOLS_BY_SUBPART_PERIODICITY = 20;

const uint DEFRAGMENTATION_PERIODICITY = 100;
const float_t PARTITION_EDGE_LENGTH_DEFAULT_UM = 10; // large for now because we have just one partition
const float_t PARTITION_EDGE_EXTRA_MARGIN_UM = 0.01;
const float_t SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT = 1;

const uint MAX_SUBPARTS_PER_PARTITION = 300;

// ---------------------------------- fixed constants and specific typedefs -------------------
const float_t POS_INVALID = FLT_MAX; // cannot be NAN because we cannot do any comparison with NANs
const float_t LENGTH_INVALID = FLT_MAX;

const float_t TIME_INVALID = -256;
const float_t TIME_FOREVER = 1e20; // based on MCell3
const float_t TIME_SIMULATION_START = 0;

const float_t SQRT2 = 1.41421356238;
const float_t RX_RADIUS_MULTIPLIER = 1.3; // TEMPORARY - we should figure out why some collisions with subparts are missed..., but maybe it won't have big perf impact...

const uint INT_INVALID = INT32_MAX;
const uint ID_INVALID = UINT32_MAX; // general invalid id, should not be used when a definition for a specific type is available
const uint ID_INVALID2 = UINT32_MAX - 1; // second general invalid id, might be used for special values
const uint INDEX_INVALID = UINT32_MAX; // general invalid index, should not be used when a definition for a specific type is available
const uint INDEX_INVALID2 = UINT32_MAX - 2; // second general invalid index, might be used for special values

// molecule id is a unique identifier of a molecule,
// no 2 molecules may have the same ID in the course of a simulation (at least for now)
typedef uint molecule_id_t;
const molecule_id_t MOLECULE_ID_INVALID = ID_INVALID;

// molecule index is index into partition's molecules array, indices and ids are
// identical until the first defragmentation that shuffles molecules in the molecules array
typedef uint molecule_index_t;
const molecule_index_t MOLECULE_INDEX_INVALID = INDEX_INVALID;

// for now, this is the partition that contains point 0,0,0 in its center
typedef uint partition_id_t;
const partition_id_t PARTITION_ID_INITIAL = 0;
const partition_id_t PARTITION_ID_INVALID = INDEX_INVALID;

typedef uint subpart_index_t;
const subpart_index_t SUBPART_INDEX_INVALID = INDEX_INVALID;

// time step is used in partition to make sets of molecules that can be diffused with
// different periodicity
const uint TIME_STEP_INDEX_INVALID = INDEX_INVALID;

const char* const NAME_INVALID = "name_invalid";
const char* const NAME_NOT_SET = "name_not_set";

const uint64_t BUCKET_INDEX_INVALID = UINT64_MAX;

const uint VERTICES_IN_TRIANGLE = 3;
const uint EDGES_IN_TRIANGLE = VERTICES_IN_TRIANGLE; // same of course as above, but different name to specify what we are counting

typedef uint vertex_index_t; // index in partition's vertices
const vertex_index_t VERTEX_INDEX_INVALID = INDEX_INVALID;

typedef uint wall_index_t; // index in partition's walls
const wall_index_t WALL_INDEX_INVALID = INDEX_INVALID;

typedef uint tile_index_t; // index of a tile in a grid
const tile_index_t TILE_INDEX_INVALID = INDEX_INVALID;

typedef uint edge_index_t; // index of an edge in a wall, must be in range 0..2
const edge_index_t EDGE_INDEX_0 = 0;
const edge_index_t EDGE_INDEX_1 = 1;
const edge_index_t EDGE_INDEX_2 = 2;
const edge_index_t EDGE_INDEX_WITHIN_WALL = 3; // used in find_edge_point
const edge_index_t EDGE_INDEX_CANNOT_TELL = 4;
const edge_index_t EDGE_INDEX_INVALID = INDEX_INVALID;


typedef uint count_buffer_id_t; // index of a tile in a grid
const count_buffer_id_t COUNT_BUFFER_ID_INVALID = INDEX_INVALID;


/* contains information about the neighbors of the tile */
class WallTileIndexPair {
public:
  WallTileIndexPair()
    : wall_index(WALL_INDEX_INVALID), tile_index(TILE_INDEX_INVALID)
    {
  }

  WallTileIndexPair(const wall_index_t wall_index_, const tile_index_t tile_index_)
    : wall_index(wall_index_), tile_index(tile_index_)
    {
  }

  wall_index_t wall_index;  /* surface grid the tile is on */
  tile_index_t tile_index;  /* index on that tile */
};

typedef uint wall_id_t; // world-unique wall id
const wall_id_t WALL_ID_INVALID = ID_INVALID;
const wall_id_t WALL_ID_NOT_IN_PARTITION = ID_INVALID2;

typedef uint region_id_t; // world-unique region id
const region_id_t REGION_ID_INVALID = ID_INVALID;

typedef uint region_index_t; // index in partition's regions
const region_index_t REGION_INDEX_INVALID = INDEX_INVALID;

typedef uint geometry_object_index_t;
const geometry_object_index_t GEOMETRY_OBJECT_INDEX_INVALID = INDEX_INVALID;

typedef uint geometry_object_id_t; // world-unique unique geometry object id
const geometry_object_id_t GEOMETRY_OBJECT_ID_INVALID = ID_INVALID;

// each partition has a list of volumes that are counted
// each counted volume is identified by a set of geometry objects it is in
// completely
typedef uint counted_volume_index_t;
// volume outside of all objects, i.e. WORLD-all objects,
// this is a special counted volume that has no geometry objects
const counted_volume_index_t COUNTED_VOLUME_INDEX_OUTSIDE_ALL = 0;

// special value that tells that the specific index must be computed
// based on the specific position
const counted_volume_index_t COUNTED_VOLUME_INDEX_INTERSECTS = INDEX_INVALID2;

const counted_volume_index_t COUNTED_VOLUME_INDEX_INVALID = INDEX_INVALID;


typedef std::pair<partition_id_t, wall_index_t> PartitionWallIndexPair;
typedef std::pair<partition_id_t, region_index_t> PartitionRegionIndexPair;
typedef std::pair<partition_id_t, vertex_index_t> PartitionVertexIndexPair;

typedef std::pair<float_t, PartitionWallIndexPair> CummAreaPWallIndexPair;

typedef std::map<geometry_object_id_t, uint_set<geometry_object_id_t>> CountedVolumesMap;

const int BASE_CONTAINER_ALLOC = 16;

#ifndef INDEXER_WA
// WARNING: std::set_intersection and possibly other algorithms do not work correctly with dense_hash_set
typedef uint_dense_hash_set<subpart_index_t> SubpartIndicesSet; // boost-based uint_set is worse on average
typedef boost::container::small_vector<subpart_index_t, BASE_CONTAINER_ALLOC>  SubpartIndicesVector;
#else
typedef std::vector<subpart_index_t> SubpartIndicesVector;
typedef std::set<subpart_index_t> SubpartIndicesSet;
#endif

template<typename T>
class insertion_ordered_set {
public:
  void insert_ordered(const T& value) {
    if (s.count(value) == 0) {
      s.insert(value);
      v.push_back(value);
    }
  }

  const T& operator [] (const size_t i) const {
    assert(i < size());
    return v[i];
  }

  size_t size() const {
    assert(v.size() == s.size());
    return v.size();
  }

  const std::vector<T>& get_as_vector() const {
    return v;
  }

private:
  std::vector<T> v;
  std::set<T> s;
};

// ---------------------------------- vector types ----------------------------------

#if FLOAT_T_BYTES == 8
typedef glm::dvec3 glm_vec3_t;
typedef glm::dvec2 glm_vec2_t;
typedef glm::dmat4x4 mat4x4;
#else
typedef glm::fvec3 glm_vec3_t;
typedef glm::fvec2 glm_vec2_t;
typedef glm::fmat4x4 mat4x4;
#endif

typedef glm::ivec3 glm_ivec3_t;
typedef glm::uvec3 UVec3;
typedef glm::bvec3 BVec3;

struct IVec3: public glm_ivec3_t {
  IVec3() = default;
  IVec3(const glm_vec3_t& a) { x = a.x; y = a.y; z = a.z; }
  IVec3(const IVec3& a) : glm_ivec3_t(a.x, a.y, a.z) { }
  IVec3(const int x_, const int y_, const int z_) { x = x_; y = y_; z = z_; }
  IVec3(const int xyz) { x = xyz; y = xyz; z = xyz; }
  IVec3(const std::vector<int>& xyz) { assert(xyz.size() == 3); x = xyz[0]; y = xyz[1]; z = xyz[2]; }
  IVec3(const vector3& v3) { x = v3.x; y = v3.y; z = v3.z; }

  // arbitrary ordering in order to use IVec3 as keys in sets and maps
  bool operator < (const IVec3& other) const {
    if (x < other.x) {
      return true;
    }
    else if (x == other.x) {
      if (y < other.y) {
        return true;
      }
      else if (y == other.y) {
        if (z < other.z) {
          return true;
        }
        else {
          return false;
        }
      }
      else {
        return false;
      }
    }
    else {
      return false;
    }
  }
};

struct Vec3: public glm_vec3_t {
  Vec3() = default;
  Vec3(const glm_vec3_t& a) { x = a.x; y = a.y; z = a.z; }
  Vec3(const Vec3& a) : glm_vec3_t(a.x, a.y, a.z) { }
  Vec3(const IVec3& a) : glm_vec3_t(a.x, a.y, a.z) { }
  Vec3(const vector3& a) { x = a.x; y = a.y; z = a.z; }
  Vec3(const float_t x_, const float_t y_, const float_t z_) { x = x_; y = y_; z = z_; }
  Vec3(const float_t xyz) { x = xyz; y = xyz; z = xyz; }
  Vec3(const std::vector<float_t>& xyz) { assert(xyz.size() == 3); x = xyz[0]; y = xyz[1]; z = xyz[2]; }

  bool is_valid() const { return !(x == POS_INVALID || y == POS_INVALID || z == POS_INVALID); }

  std::string to_string() const;
  void dump(const std::string extra_comment, const std::string ind) const;
};

// usually are .u and .v used to access contained values
struct Vec2: public glm_vec2_t {
  Vec2() = default;
  Vec2(const glm_vec2_t& a) { x = a.x; y = a.y; }
  Vec2(const Vec2& a) : glm_vec2_t(a.x, a.y) { }
  Vec2(const vector2& a) { x = a.u; y = a.v; }
  Vec2(const float_t x_, const float_t y_) { x = x_; y = y_; }
  Vec2(const float_t xy) { x = xy; y = xy; }
  Vec2(const std::vector<float_t>& xy) { assert(xy.size() == 2); x = xy[0]; y = xy[1]; }

  bool is_valid() const { return !(x == POS_INVALID || y == POS_INVALID); }

  std::string to_string() const;
  void dump(const std::string extra_comment, const std::string ind) const;
};



static inline std::ostream & operator<<(std::ostream &out, const Vec3& a) {
  out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return out;
}
static inline std::ostream & operator<<(std::ostream &out, const Vec2& a) {
  out << "(" << a.u << ", " << a.v << ")";
  return out;
}


// ---------------------------------- auxiliary functions ----------------------------------

static inline float_t log_f(const float_t x) {
  assert(x != 0);
#if FLOAT_T_BYTES == 8
  return log(x);
#else
  return logf(x);
#endif
}

static inline float_t ceil_f(const float_t x) {
#if FLOAT_T_BYTES == 8
  return ceil(x);
#else
  return ceilf(x);
#endif
}

static inline float_t floor_to_multiple(const float_t val, float_t multiple) {
  assert(val >= 0);
  return (float_t)((int)((val + EPS)/ multiple)) * multiple;
}

static inline float_t floor_to_multiple_allow_negative(const float_t val, float_t multiple) {
  if (val >= 0) {
    return (float_t)((int)((val + EPS)/ multiple)) * multiple;
  }
  else {
    // we need to floor towards the lower negative value, the code above would
    // ceil the value for negative inputs
    return (float_t)((int)((val + EPS - multiple)/ multiple)) * multiple;
  }
}

static inline Vec3 floor_to_multiple(const Vec3& val, float_t multiple) {
  assert(val.x >= 0 && val.y >=0 && val.z >= 0);
  return (Vec3)((glm::ivec3)((val + Vec3(EPS))/ multiple)) * multiple;
}

static inline Vec3 floor_to_multiple_allow_negative(const Vec3& val, float_t multiple) {
  Vec3 res;
  res.x = floor_to_multiple_allow_negative(val.x, multiple);
  res.y = floor_to_multiple_allow_negative(val.y, multiple);
  res.z = floor_to_multiple_allow_negative(val.z, multiple);
  return res;
}

static inline float_t ceil_to_multiple(const float_t val, float_t multiple) {
  assert(val >= 0);
  float_t res = floor_to_multiple(val, multiple);
  if (!cmp_eq(val, res)) {
    res += multiple;
  }
  return res;
}

static inline Vec3 ceil_to_multiple(const Vec3& val, float_t multiple) {
  Vec3 res;
  res.x = ceil_to_multiple(val.x, multiple);
  res.y = ceil_to_multiple(val.y, multiple);
  res.z = ceil_to_multiple(val.z, multiple);
  return res;
}

static inline bool cmp_eq(const Vec3& a, const Vec3& b, const float_t eps) {
  return cmp_eq(a.x, b.x, eps) && cmp_eq(a.y, b.y, eps) && cmp_eq(a.z, b.z, eps);
}

static inline bool cmp_eq(const Vec3& a, const Vec3& b) {
  return cmp_eq(a, b, EPS);
}

static inline bool cmp_eq(const Vec2& a, const Vec2& b, const float_t eps) {
  return cmp_eq(a.x, b.x, eps) && cmp_eq(a.y, b.y, eps);
}

static inline bool cmp_eq(const Vec2& a, const Vec2& b) {
  return cmp_eq(a, b, EPS);
}

static inline bool cmp_lt(const float_t a, const float_t b, const float_t eps) {
  return a < b && !cmp_eq(a, b, eps);
}

static inline bool cmp_gt(const float_t a, const float_t b, const float_t eps) {
  return a > b && !cmp_eq(a, b, eps);
}


static inline uint powu(const uint a, const uint n) {
  uint res = a;
  for (uint i = 1; i < n; i++) {
    res *= a;
  }
  return res;
}


static inline float_t max3(const Vec3& v) {
  return glm::compMax((glm_vec3_t)v);
}

static inline float_t min3(const Vec3& v) {
  return glm::compMin((glm_vec3_t)v);
}

static inline int max3_i(const IVec3& v) {
  return glm::compMax((glm_ivec3_t)v);
}

static inline int min3_i(const IVec3& v) {
  return glm::compMin((glm_ivec3_t)v);
}

static inline Vec3 abs3(const Vec3& v) {
  return glm::abs((glm_vec3_t)v);
}

static inline Vec3 floor3(const Vec3& v) {
  return glm::floor((glm_vec3_t)v);
}

static inline Vec3 round3(const Vec3& v) {
  return glm::round((glm_vec3_t)v);
}

/* abs_max_2vec picks out the largest (absolute) value found among two vectors
 * (useful for properly handling floating-point rounding error). */
static inline float_t abs_max_2vec(const Vec3& v1, const Vec3& v2) {
  glm_vec3_t v1abs = abs3(v1);
  glm_vec3_t v2abs = abs3(v2);
  Vec3 vmax = glm::max(v1abs, v2abs);
  return MCell::max3(vmax);
}

static inline float_t determinant2(const Vec2& v1, const Vec2& v2) {
  return v1.u * v2.v - v1.v * v2.u;
}

static inline float_t dot2(const Vec2& v1, const Vec2& v2) {
  return glm::dot((glm_vec2_t)v1, (glm_vec2_t)v2);
}

static inline float_t len2_squared(const Vec2& v1) {
  return v1.u * v1.u + v1.v * v1.v;
}

// FIXME: use places where to use this function
static inline float_t len2(const Vec2& v1) {
  return sqrt_f(len2_squared(v1));
}

static inline float_t dot(const Vec3& v1, const Vec3& v2) {
  return glm::dot((glm_vec3_t)v1, (glm_vec3_t)v2);
}

static inline float_t len3_squared(const Vec3& v1) {
  return v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
}

// FIXME: use places where to use this function
static inline float_t len3(const Vec3& v1) {
  return sqrt_f(len3_squared(v1));
}

static inline float_t distance3(const Vec3& v1, const Vec3& v2) {
  return sqrt_f( len3_squared(v1 - v2) );
}

static inline uint get_largest_abs_dim_index(const Vec3& v) {
  Vec3 a = glm::abs(glm_vec3_t(v));
  if (a.x > a.y) {
    if (a.x > a.z) {
      return 0; // x
    }
    else {
      return 2; // z
    }
  }
  else {
    if (a.y > a.z) {
      return 1; // y
    }
    else {
      return 2; // z
    }

  }
}

/***************************************************************************
point_in_box:
  In:  pt - we want to find out if this point is in the box
       llf - lower left front corner of the box
       urb - upper right back corner of the box

  Out: Returns true if point is in box, 0 otherwise
***************************************************************************/
static inline bool point_in_box(const Vec3& pt, const Vec3& llf, const Vec3& urb) {
  return
      pt.x >= llf.x && pt.x <= urb.x &&
      pt.y >= llf.y && pt.y <= urb.y &&
      pt.z >= llf.z && pt.z <= urb.z;
}

/**
 * Performs vector cross product.
 * Computes the cross product of two vector3's v1 and v2 storing the result
 * in vector3 v3.
 */
static inline Vec3 cross(const Vec3& v1, const Vec3& v2) {
  return glm::cross((glm_vec3_t)v1, (glm_vec3_t)v2);
}


static inline int distinguishable_vec2(const Vec2& a, const Vec2& b, float_t eps) {
  float_t c, cc, d;

  /* Find largest coordinate */
  c = fabs_f(a.u);

  d = fabs_f(a.v);
  if (d > c)
    c = d;

  d = fabs_f(b.u);
  if (d > c)
    c = d;

  d = fabs_f(b.v);
  if (d > c)
    c = d;

  /* Find largest difference */
  cc = fabs_f(a.u - b.u);

  d = fabs_f(a.v - b.v);
  if (d > cc)
    cc = d;

  /* Make sure fractional difference is at least eps and absolute difference is
   * at least (eps*eps) */
  if (c < eps)
    c = eps;
  return (c * eps < cc);
}


static inline bool distinguishable_vec3(const Vec3& a, const Vec3& b, float_t eps) {
  float_t c, cc, d;

  /* Find largest coordinate */
  c = fabs_f(a.x);

  d = fabs_f(a.y);
  if (d > c)
    c = d;

  d = fabs_f(a.z);
  if (d > c)
    c = d;

  d = fabs_f(b.x);
  if (d > c)
    c = d;

  d = fabs_f(b.y);
  if (d > c)
    c = d;

  d = fabs_f(b.z);
  if (d > c)
    c = d;

  /* Find largest difference */
  cc = fabs_f(a.x - b.x);

  d = fabs_f(a.y - b.y);
  if (d > cc)
    cc = d;

  d = fabs_f(a.z - b.z);
  if (d > cc)
    cc = d;

  /* Make sure fractional difference is at least eps and absolute difference is
   * at least (eps*eps) */
  if (c < eps)
    c = eps;
  return (c * eps < cc);
}


static inline void guard_zero_div(Vec3& val) {
  if (val.x == 0) {
    val.x = FLT_MIN;
  }
  if (val.y == 0) {
    val.y = FLT_MIN;
  }
  if (val.z == 0) {
    val.z = FLT_MIN;
  }
}

/*
 * Stats collected during simulation, contains also number of the current iteration
 */
class SimulationStats {
public:
  SimulationStats() {
    reset();
  }
  void inc_ray_voxel_tests() {
    ray_voxel_tests++;
  }
  void inc_ray_polygon_tests() {
    ray_polygon_tests++;
  }
  void inc_ray_polygon_colls() {
    ray_polygon_colls++;
  }

  void inc_mol_wall_reflections() {
    mol_wall_reflections++;
  }

  void inc_vol_mol_vol_mol_collisions() {
    vol_mol_vol_mol_collisions++;
  }

  // new mcell4 stats
  void inc_mol_moves_between_walls() {
    mol_moves_between_walls++;
  }

  void inc_recomputations_of_counted_volume() {
    recomputations_of_counted_volume++;
  }

  void inc_num_waypoints_used() {
    num_waypoints_used++;
  }
  void inc_diffuse_3d_calls() {
    diffuse_3d_calls++;
  }

  void dump();

  const uint64_t& get_current_iteration() const {
    return current_iteration;
  }

  uint64_t& get_current_iteration() {
    return current_iteration;
  }

  void reset() {
    current_iteration = 0;
    ray_voxel_tests = 0;
    ray_polygon_tests = 0;
    ray_polygon_colls = 0;
    mol_wall_reflections = 0;
    vol_mol_vol_mol_collisions = 0;
    mol_moves_between_walls = 0;
    num_waypoints_used = 0;
    recomputations_of_counted_volume = 0;
    diffuse_3d_calls = 0;
  }

private:
  uint64_t current_iteration;

  uint64_t ray_voxel_tests;
  uint64_t ray_polygon_tests;
  uint64_t ray_polygon_colls;

  uint64_t mol_wall_reflections;
  uint64_t vol_mol_vol_mol_collisions;

  uint64_t mol_moves_between_walls;

  uint64_t num_waypoints_used;
  uint64_t recomputations_of_counted_volume;

  uint64_t diffuse_3d_calls;
};

/*
 * Constant data set in initialization useful for all classes, single object is owned by world
 */
// TODO: cleanup all unnecessary argument passing, e.g. in diffuse_react_event.cpp
class SimulationConfig: public BNG::BNGConfig {
public:
  // TODO: add ctor

  // configuration

  float_t vacancy_search_dist2; /* Square of distance to search for free grid
                                  location to place surface product */

  Vec3 partition0_llf;

  float_t partition_edge_length; // TODO: rename to side
  uint num_subpartitions_per_partition;
  uint num_subpartitions_per_partition_squared;
  float_t subpartition_edge_length; // == partition_edge_length / subpartitions_per_partition_dimension
  float_t subpartition_edge_length_rcp; // == 1/subpartition_edge_length

  // other options
  bool use_expanded_list; /* If set, check neighboring subvolumes for mol-mol
                            interactions */
  bool randomize_smol_pos; /* If set, always place surface molecule at random
                             location instead of center of grid */

  // initialized in World::init_counted_volumes
  // also tells whether waypoints in a partition were initialized
  bool has_intersecting_counted_objects;


  void init() {
    BNGConfig::init();
    init_subpartition_edge_length();
  }

  void dump();

private:
  void init_subpartition_edge_length() {
    // probably not needed anymore
    //assert(num_subpartitions_per_partition % 2 == 0
    //    && "Implementation of raycast_with_endpoints requires that central subparts are aligned with the axes and not shifted");
    release_assert(partition_edge_length > 0);
    subpartition_edge_length = partition_edge_length / (float_t)num_subpartitions_per_partition;
    subpartition_edge_length_rcp = 1.0/subpartition_edge_length;
    num_subpartitions_per_partition_squared = powu(num_subpartitions_per_partition, 2);
  }

};
} // namespace mcell

#endif // SRC4_DEFINES_H_
