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

#ifndef SWIG
#include <stdint.h>
#include <cassert>
#include <cmath>


#include "../libs/bng/common_defines.h"
#include "../libs/bng/bng_defines.h"

#include "mcell_structs.h"
#include "debug_config.h"

// warning: do not use floating point types from directly,
// we need to be able to control the precision
#include "../libs/glm/glm.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include "../libs/glm/gtx/component_wise.hpp"

// this file must not depend on any other from mcell4 otherwise there
// might be some nasty cyclic include dependencies
#endif

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

typedef Common::float_t float_t;

// ---------------------------------- optimization macros ----------------------------------
#if defined(likely) or defined(unlikely)
#error "Macros 'likely' or 'unlikely' are already defined"
#endif

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)



const float_t SCHEDULER_COMPARISON_EPS = 1e-10;

// ---------------------------------- configurable constants----------------------------------

const uint DEFRAGMENTATION_PERIODICITY = 500;
const float_t PARTITION_EDGE_LENGTH_DEFAULT = 10 * 100 /*100 = 1/length unit*/; // large for now because we have just one partition
const float_t SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT = 1;


// ---------------------------------- fixed constants and specific typedefs -------------------

const float_t DIFFUSION_CONSTANT_ZER0 = 0;

// unique species id
typedef uint species_id_t;
const species_id_t SPECIES_ID_INVALID = ID_INVALID;


// molecule index is index into partition's molecules array, indices and ids are
// identical until the first defragmentation that shuffles molecules in the molecules array
typedef uint molecule_index_t;
const molecule_index_t MOLECULE_INDEX_INVALID = INDEX_INVALID;

// for now, this is the partition that contains point 0,0,0 in its center
typedef uint partition_index_t;
const partition_index_t PARTITION_INDEX_INITIAL = 0;
const partition_index_t PARTITION_INDEX_INVALID = INDEX_INVALID;

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

typedef uint region_index_t; // index in partition's regions
const region_index_t REGION_INDEX_INVALID = INDEX_INVALID;

typedef uint geometry_object_index_t;
const geometry_object_index_t GEOMETRY_OBJECT_INDEX_INVALID = INDEX_INVALID;

typedef uint geometry_object_id_t; // world-unique unique geometry object id
const geometry_object_id_t GEOMETRY_OBJECT_ID_INVALID = ID_INVALID;

typedef int32_t orientation_t;
const orientation_t ORIENTATION_DOWN = -1;
const orientation_t ORIENTATION_NONE = 0;
const orientation_t ORIENTATION_UP = 1;

typedef std::pair<partition_index_t, wall_index_t> PartitionWallIndexPair;
typedef std::pair<partition_index_t, region_index_t> PartitionRegionIndexPair;
typedef std::pair<partition_index_t, vertex_index_t> PartitionVertexIndexPair;

typedef std::pair<float_t, PartitionWallIndexPair> CummAreaPWallIndexPair;


class Reaction;

#ifndef INDEXER_WA
typedef boost::container::small_vector<subpart_index_t, 8>  SubpartIndicesVector;
typedef boost::container::small_vector<const Reaction*, 8>  ReactionsVector;
#else
typedef std::vector<subpart_index_t> SubpartIndicesVector;
typedef std::vector<const Reaction*> ReactionsVector;
#endif

// ---------------------------------- vector types ----------------------------------

typedef glm::dvec3 glm_vec3_t;
typedef glm::dvec2 glm_vec2_t;
typedef glm::ivec3 ivec3_t;
typedef glm::uvec3 uvec3_t;
typedef glm::bvec3 bvec3_t;
typedef glm::dmat4x4 mat4x4;

// NOTE: rename to Vec3? - not sure, this is a really simple object...
struct vec3_t: public glm_vec3_t {
  vec3_t() = default;
  vec3_t(const glm_vec3_t& a) { x = a.x; y = a.y; z = a.z; }
  vec3_t(const vec3_t& a) : glm_vec3_t(a.x, a.y, a.z) { }
  vec3_t(const ivec3_t& a) : glm_vec3_t(a.x, a.y, a.z) { }
  vec3_t(const vector3& a) { x = a.x; y = a.y; z = a.z; }
  vec3_t(const float_t x_, const float_t y_, const float_t z_) { x = x_; y = y_; z = z_; }
  vec3_t(const float_t xyz) { x = xyz; y = xyz; z = xyz; }

  bool is_valid() const { return !(x == POS_INVALID || y == POS_INVALID || z == POS_INVALID); }

  std::string to_string() const;
  void dump(const std::string extra_comment, const std::string ind) const;
};

// usually are .u and .v used to access contained values
struct vec2_t: public glm_vec2_t {
  vec2_t() = default;
  vec2_t(const glm_vec2_t& a) { x = a.x; y = a.y; }
  vec2_t(const vec2_t& a) : glm_vec2_t(a.x, a.y) { }
  vec2_t(const vector2& a) { x = a.u; y = a.v; }
  vec2_t(const float_t x_, const float_t y_) { x = x_; y = y_; }
  vec2_t(const float_t xy) { x = xy; y = xy; }

  bool is_valid() const { return !(x == POS_INVALID || y == POS_INVALID); }

  std::string to_string() const;
  void dump(const std::string extra_comment, const std::string ind) const;
};

std::ostream & operator<<(std::ostream &out, const vec3_t &a);
std::ostream & operator<<(std::ostream &out, const vec2_t &a);


// ---------------------------------- auxiliary functions ----------------------------------

static inline float_t sqrt_f(const float_t x) {
#if FLOAT_T_BYTES == 8
  return sqrt(x);
#else
  return sqrtf(x);
#endif
}

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

static inline float_t fabs_f(const float_t x) {
#if FLOAT_T_BYTES == 8
  return fabs(x);
#else
  return fabsf(x);
#endif
}


static inline float_t floor_to_multiple(const float_t val, float_t multiple) {
  return (float_t)((int)(val / multiple)) * multiple;
}

static inline vec3_t floor_to_multiple(const vec3_t& val, float_t multiple) {
  return (vec3_t)((glm::ivec3)(val / multiple)) * multiple;
}

static inline bool cmp_eq(const float_t a, const float_t b, const float_t eps) {
  return fabs_f(a - b) < eps;
}

static inline bool cmp_eq(float_t a, float_t b) {
  return cmp_eq(a, b, EPS);
}

static inline bool cmp_eq(const vec3_t& a, const vec3_t& b, const float_t eps) {
  return cmp_eq(a.x, b.x, eps) && cmp_eq(a.y, b.y, eps) && cmp_eq(a.z, b.z, eps);
}

static inline bool cmp_eq(const vec3_t& a, const vec3_t& b) {
  return cmp_eq(a, b, EPS);
}

static inline bool cmp_eq(const vec2_t& a, const vec2_t& b, const float_t eps) {
  return cmp_eq(a.x, b.x, eps) && cmp_eq(a.y, b.y, eps);
}

static inline bool cmp_eq(const vec2_t& a, const vec2_t& b) {
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

static inline float_t max3d(const vec3_t& v) {
  return glm::compMax((glm_vec3_t)v);
}

/* abs_max_2vec picks out the largest (absolute) value found among two vectors
 * (useful for properly handling floating-point rounding error). */
static inline float_t abs_max_2vec(const vec3_t& v1, const vec3_t& v2) {
  glm_vec3_t v1abs = glm::abs((glm_vec3_t)v1);
  glm_vec3_t v2abs = glm::abs((glm_vec3_t)v2);
  vec3_t vmax = glm::max(v1abs, v2abs);
  return MCell::max3d(vmax);
}

static inline float_t determinant2(const vec2_t& v1, const vec2_t& v2) {
  return v1.u * v2.v - v1.v * v2.u;
}

static inline float_t dot2(const vec2_t& v1, const vec2_t& v2) {
  return glm::dot((glm_vec2_t)v1, (glm_vec2_t)v2);
}

static inline float_t len2_squared(const vec2_t& v1) {
  return v1.u * v1.u + v1.v * v1.v;
}

// FIXME: use places where to use this function
static inline float_t len2(const vec2_t& v1) {
  return sqrt_f(len2_squared(v1));
}

static inline float_t dot(const vec3_t& v1, const vec3_t& v2) {
  return glm::dot((glm_vec3_t)v1, (glm_vec3_t)v2);
}

static inline float_t len3_squared(const vec3_t& v1) {
  return v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
}

// FIXME: use places where to use this function
static inline float_t len3(const vec3_t& v1) {
  return sqrt_f(len3_squared(v1));
}

static inline float_t distance3(const vec3_t& v1, const vec3_t& v2) {
  return sqrt_f( len3_squared(v1 - v2) );
}

/***************************************************************************
point_in_box:
  In:  pt - we want to find out if this point is in the box
       llf - lower left front corner of the box
       urb - upper right back corner of the box

  Out: Returns true if point is in box, 0 otherwise
***************************************************************************/
static inline bool point_in_box(const vec3_t& pt, const vec3_t& llf, const vec3_t& urb) {
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
static inline vec3_t cross(const vec3_t& v1, const vec3_t& v2) {
  return glm::cross((glm_vec3_t)v1, (glm_vec3_t)v2);
}

// returns true when whether two values are measurably different
inline bool distinguishable_f(float_t a, float_t b, float_t eps) {
  float_t c = fabs_f(a - b);
  a = fabs_f(a);
  if (a < 1) {
    a = 1;
  }
  b = fabs_f(b);

  if (b < a) {
    eps *= a;
  } else {
    eps *= b;
  }
  return (c > eps);
}



static inline int distinguishable_vec2(const vec2_t& a, const vec2_t& b, float_t eps) {
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


static inline bool distinguishable_vec3(const vec3_t& a, const vec3_t& b, float_t eps) {
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


static inline void debug_guard_zero_div(float_t& val) {
#ifndef NDEBUG
  // if we divide by such a small number, result is practically the same as
  // if we would return inf during division
  if (val == 0) {
    val = FLT_MIN;
  }
#endif
}

static inline void debug_guard_zero_div(vec3_t& val) {
#ifndef NDEBUG
  if (val.x == 0) {
    val.x = FLT_MIN;
  }
  if (val.y == 0) {
    val.y = FLT_MIN;
  }
  if (val.z == 0) {
    val.z = FLT_MIN;
  }
#endif
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

  // new mcell4 stats
  void inc_mol_moves_between_walls() {
    mol_moves_between_walls++;
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
    mol_moves_between_walls = 0;
  }

private:
  uint64_t current_iteration;

  uint64_t ray_voxel_tests;
  uint64_t ray_polygon_tests;
  uint64_t ray_polygon_colls;
  uint64_t mol_moves_between_walls;
};

/*
 * Constant data set in initialization useful for all classes, single object is owned by world
 */
// TODO: cleanup all unnecessary argument passing, e.g. in diffuse_react_event.cpp
class SimulationConfig {
public:
  // configuration
  float_t time_unit;
  float_t length_unit;
  float_t rx_radius_3d;
  float_t vacancy_search_dist2;

  float_t partition_edge_length;
  uint subpartitions_per_partition_dimension;
  uint subpartitions_per_partition_dimension_squared;
  float_t subpartition_edge_length; // == partition_edge_length / subpartitions_per_partition_dimension
  float_t subpartition_edge_length_rcp; // == 1/subpartition_edge_length

  // other options
  bool use_expanded_list;
  bool randomize_smol_pos;

  void init() {
    init_subpartition_edge_length();
  }

  void dump();

private:
  void init_subpartition_edge_length() {
    if (partition_edge_length != 0) {
      subpartition_edge_length = partition_edge_length / (float_t)subpartitions_per_partition_dimension;
      subpartition_edge_length_rcp = 1.0/subpartition_edge_length;
    }
    subpartitions_per_partition_dimension_squared = powu(subpartitions_per_partition_dimension, 2);
  }

};
} // namespace mcell

#endif // SRC4_DEFINES_H_
