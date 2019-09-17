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

#ifndef NDEBUG
#define INDEXER_WA // Don't know yet how to convince Eclipse to correctly index boost containers
#endif

#if defined(NDEBUG) && defined(INDEXER_WA)
#warning "INDEXER_WA is enabled and this will lead to lower performance"
#endif

#include <stdint.h>
#include <vector>
#include <string>
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <unordered_map>
#include "../libs/boost/container/small_vector.hpp"
#include "../libs/boost/container/flat_set.hpp"

#include "mcell_structs.h"
#include "debug_config.h"

// warning: do not use floating point types from directly,
// we need to be able to control the precision
#include "../libs/glm/glm.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include "../libs/glm/gtx/component_wise.hpp"

// this file must not depend on any other from mcell4 otherwise there
// might be some nasty cyclic include dependencies


namespace MCell {

// ---------------------------------- optimization macros ----------------------------------
#if defined(likely) or defined(unlikely)
#error "Macros 'likely' or 'unlikely' are already defined"
#endif

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

// ---------------------------------- float types ----------------------------------

typedef double float_t; // will be changed to float
#define FLOAT_T_BYTES 8

#if FLOAT_T_BYTES == 8
const float_t EPS = 1e-12; // same as EPS_C
const float_t SQRT_EPS = 1e-6;
const float_t GIGANTIC4 = 1e140;
#else
#error "Base type float32 is not supported yet"
#endif

const float_t SCHEDULER_COMPARISON_EPS = 1e-10;

// ---------------------------------- configurable constants----------------------------------

const uint DEFRAGMENTATION_PERIODICITY = 500;
const float_t PARTITION_EDGE_LENGTH_DEFAULT = 10 * 100 /*100 = 1/length unit*/; // large for now because we have just one partition
const float_t SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT = 1;


// ---------------------------------- fixed costants and specific typedefs -------------------
const float_t POS_INVALID = FLT_MAX; // cannot be NAN because we cannot do any comparison with NANs

const float_t TIME_INVALID = -1;
const float_t TIME_FOREVER = FLT_MAX; // this max is sufficient for both float and double
const float_t TIME_SIMULATION_START = 0;

const float_t DIFFUSION_CONSTANT_ZER0 = 0;

const uint INDEX_INVALID = UINT32_MAX;

// unique species id
typedef uint species_id_t;
const species_id_t SPECIES_ID_INVALID = UINT32_MAX;

// molecule id is a unique identifier of a molecule,
// no 2 molecules may have the same ID in the course of a simulation (at least for now)
typedef uint molecule_id_t;
const molecule_id_t MOLECULE_ID_INVALID = UINT32_MAX;

// molecule index is index into partition's molecules array, indices and ids are
// identical until the first defragmentation that shuffles molecules in the molecules array
typedef uint molecule_index_t;
const molecule_index_t MOLECULE_INDEX_INVALID = UINT32_MAX;

// for now, this is the partition that contains point 0,0,0 in its center
typedef uint partition_index_t;
const partition_index_t PARTITION_INDEX_INITIAL = 0;
const partition_index_t PARTITION_INDEX_INVALID = UINT32_MAX;

typedef uint subpart_index_t;
const subpart_index_t SUBPART_INDEX_INVALID = UINT32_MAX;

// time step is used in partition to make sets of molecules that can be diffused with
// different periodicity
const uint TIME_STEP_INDEX_INVALID = UINT32_MAX;

const char* const NAME_INVALID = "name_invalid";
const char* const NAME_NOT_SET = "name_not_set";

const uint64_t BUCKET_INDEX_INVALID = UINT64_MAX;

const uint VERTICES_IN_TRIANGLE = 3;
const uint EDGES_IN_TRIANGLE = VERTICES_IN_TRIANGLE; // same of course as above, but different name to specify what we are counting

typedef uint vertex_index_t; // index in partition's vertices
const vertex_index_t VERTEX_INDEX_INVALID = UINT32_MAX;

//typedef uint grid_index_t; // index in partition's walls
//const grid_index_t GRID_INDEX_INVALID = UINT32_MAX;

typedef uint tile_index_t; // index of a tile in a grid
const tile_index_t TILE_INDEX_INVALID = UINT32_MAX;

typedef uint wall_index_t; // index in partition's walls
const wall_index_t WALL_INDEX_INVALID = UINT32_MAX;

typedef uint wall_id_t; // world-unique wall id
const wall_id_t WALL_ID_INVALID = UINT32_MAX;

//typedef uint wall_class_index_t; // index in world's wall classes

typedef uint geometry_object_index_t;
const geometry_object_index_t GEOMETRY_OBJECT_INDEX_INVALID = UINT32_MAX;

typedef uint geometry_object_id_t; // world-unique unique geometry object id
const geometry_object_id_t GEOMETRY_OBJECT_ID_INVALID = UINT32_MAX;

typedef int32_t orientation_t;
const orientation_t ORIENTATION_DOWN = -1;
const orientation_t ORIENTATION_NONE = 0;
const orientation_t ORIENTATION_UP = 1;

typedef std::pair<partition_index_t, wall_index_t> partition_wall_index_pair_t;
typedef std::pair<partition_index_t, vertex_index_t> partition_vertex_index_pair_t;

typedef std::pair<float_t, partition_wall_index_pair_t> cum_area_pwall_index_pair_t;


class reaction_t;
#ifndef INDEXER_WA
template<class T, class Allocator=boost::container::new_allocator<T>>
using small_vector = boost::container::small_vector<T, 8, Allocator>;

typedef boost::container::small_vector<subpart_index_t, 8>  subpart_indices_vector_t;
typedef boost::container::small_vector<const reaction_t*, 8>  reactions_vector_t;
#else
template<typename T, typename _Alloc = std::allocator<T>  >
using small_vector = std::vector<T, _Alloc>;

typedef std::vector<subpart_index_t> subpart_indices_vector_t;
typedef std::vector<const reaction_t*> reactions_vector_t;
#endif


/**
 * Class used to hold sets of ids or indices of molecules or other items
 */
class uint_set_t: public boost::container::flat_set<uint> {
public:
  void set_contains_id(const uint id, const bool value = true) {
    if (value) {
      assert(count(id) == 0);
      insert(id);
    }
    else {
      assert(count(id) == 1);
      erase(id);
    }
  }

  void dump();
};

// ---------------------------------- vector types ----------------------------------

typedef glm::dvec3 glm_vec3_t;
typedef glm::dvec2 glm_vec2_t;
typedef glm::ivec3 ivec3_t;
typedef glm::uvec3 uvec3_t;
typedef glm::bvec3 bvec3_t;
typedef glm::dmat4x4 mat4x4;

struct vec3_t: public glm_vec3_t{
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
struct vec2_t: public glm_vec2_t{
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

inline bool cmp_eq(float_t a, float_t b) {
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

static inline float_t dot(const vec3_t& v1, const vec3_t& v2) {
  return glm::dot((glm_vec3_t)v1, (glm_vec3_t)v2);
}

static inline float_t len3_squared(const vec3_t& v1) {
  return v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
}

static inline float_t distance3(const vec3_t& v1, const vec3_t& v2) {
  return sqrt_f( len3_squared(v1 - v2) );
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


class reaction_t;
#ifndef INDEXER_WA
typedef std::unordered_map<species_id_t, reaction_t*> species_reaction_map_t;
typedef species_reaction_map_t unimolecular_reactions_map_t;
typedef std::unordered_map< species_id_t, species_reaction_map_t > bimolecular_reactions_map_t;
#else
typedef std::map<species_id_t, reaction_t*> species_reaction_map_t;
typedef species_reaction_map_t unimolecular_reactions_map_t;
typedef std::map<species_id_t, species_reaction_map_t> bimolecular_reactions_map_t;
#endif

/*
 * Constant data set in initialization useful for all classes, single object is owned by world
 */
struct simulation_stats_t {
  simulation_stats_t()
    : ray_voxel_tests(0), ray_polygon_tests(0), ray_polygon_colls(0) {
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

  void dump();
private:
  uint64_t ray_voxel_tests;
  uint64_t ray_polygon_tests;
  uint64_t ray_polygon_colls;
};

} // namespace mcell

#endif // SRC4_DEFINES_H_
