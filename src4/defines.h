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

// NOTES:
// rename partition to part?

#ifndef SRC4_DEFINES_H_
#define SRC4_DEFINES_H_

#include <stdint.h>

#include <vector>
#include <string>
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <unordered_map>

#include "mcell_structs.h"

#include "debug_config.h"

//#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES
//#define GLM_FORCE_AVX2
//#define GLM_FORCE_PURE

// warning: do not use directly, we need to be able to control the precision
#include "../libs/glm/glm.hpp"

// this file must not depend on any other from mcell4 otherwise there
// might be some nasty include dependencies



namespace mcell {

typedef double float_t; // soon to be changed to float
#define FLOAT_T_BYTES 8

#if FLOAT_T_BYTES == 8
#define EPS 1e-12 // same as EPS_C
#else
#error "TODO: float32"
#endif

const float_t TIME_INVALID = NAN;
const float_t TIME_SIMULATION_START = 0;

const uint32_t DEFRAGMENTATION_PERIODICITY = 500;

const float_t PARTITION_EDGE_LENGTH_DEFAULT = 10 * 100 /*1/length unit*/; // maybe too large for now

//const float_t SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT = 20; // mcell3 uses logarithmic scaling, this is not useful here
const float_t SUBPARTITIONS_PER_PARTITION_DIMENSION_DEFAULT = 1;

const float_t SCHEDULER_COMPARISON_EPS = 1e-10;

typedef glm::dvec3 glm_vec3_t;
typedef glm::ivec3 ivec3_t;
typedef glm::bvec3 bvec3_t;

struct vec3_t: public glm_vec3_t{

	vec3_t() : glm::dvec3(0) {}
	vec3_t(const glm::dvec3& a) { x = a.x; y = a.y; z = a.z; }
	vec3_t(const vec3_t& a) : glm::dvec3(a.x, a.y, a.z) { /*x = a.x; y = a.y; z = a.z;*/ }
	vec3_t(const vector3& a) { x = a.x; y = a.y; z = a.z; }
	vec3_t(const float_t x_, const float_t y_, const float_t z_) { x = x_; y = y_; z = z_; }
	vec3_t(const float_t xyz) { x = xyz; y = xyz; z = xyz; }


	/*ivec3_t get_direction_vec() {
		ivec3_t res;
		// maybe I do not need the 0 value
		res.x = (x == 0) ? 0 : ((x > 0) ? 1 : -1);
		res.y = (y == 0) ? 0 : ((y > 0) ? 1 : -1);
		res.z = (z == 0) ? 0 : ((z > 0) ? 1 : -1);
		return res;
	}*/

	// exact match
	//bool operator == (const vec3_t& a) { return x == a.x && y == a.y && z == a.z; }
	void dump(const std::string extra_comment, const std::string ind);
};


std::ostream & operator<<(std::ostream &out, const vec3_t &a);

//const int MAX_MOLECULES_PER_PARTITION = 32*32*32 /*32k*/; //temporary, must work dynamically




//typedef uint16_t species_id_t; - slower
typedef uint32_t species_id_t;
const int SPECIES_ID_INVALID = USHRT_MAX;

typedef uint32_t molecule_idx_t;
const uint32_t MOLECULE_IDX_INVALID = UINT32_MAX;

const uint32_t MOLECULE_INDEX_INVALID = UINT32_MAX;

// for now, this is the partition that contains point 0,0,0
// in its center
const uint32_t PARTITION_INDEX_INITIAL = 0;

const uint32_t PARTITION_INDEX_INVALID = UINT32_MAX;
const uint32_t SUBPARTITION_INDEX_INVALID = UINT32_MAX;

const uint32_t TIME_STEP_INDEX_INVALID = UINT32_MAX;

typedef glm::dmat4x4 mat4x4;

const char* const NAME_INVALID = "invalid_name";

static inline float_t floor_to_multiple(const float_t val, float_t multiple) {
	return (float_t)((int)(val / multiple)) * multiple;
}

static inline vec3_t floor_to_multiple(const vec3_t& val, float_t multiple) {
	return (vec3_t)((glm::ivec3)(val / multiple)) * multiple;
}

static inline bool cmp_eq(const float_t a, const float_t b, const float_t eps) {
	return fabs(a - b) < eps;
}

static inline bool cmp_lt(const float_t a, const float_t b, const float_t eps) {
	return a < b && !cmp_eq(a, b, eps);
}

static inline uint32_t powu(const uint32_t a, const uint32_t n) {
	uint32_t res = a;
	for (uint32_t i = 1; i < n; i++) {
		res *= a;
	}
	return res;
}

class reaction_t;
typedef std::unordered_map<species_id_t, reaction_t*> species_reaction_map_t;
typedef std::unordered_map< species_id_t, species_reaction_map_t > bimolecular_reactions_map_t;

// constants useful for all classes, single objectis owned by world
struct world_constants_t {
  float_t time_unit;
  float_t length_unit;

  float_t rx_radius_3d;

  float_t partition_edge_length;
  uint32_t subpartitions_per_partition_dimension;
  uint32_t subpartitions_per_partition_dimension_squared;
  float_t subpartition_edge_length; // == partition_edge_length / subpartitions_per_partition_dimension
  float_t subpartition_edge_length_rcp; // == 1/subpartition_edge_length

  const bimolecular_reactions_map_t* bimolecular_reactions_map; // owned by world


  void init_subpartition_edge_length() {
  	if (partition_edge_length != 0) {
			subpartition_edge_length = partition_edge_length / (float_t)subpartitions_per_partition_dimension;
			subpartition_edge_length_rcp = 1.0/subpartition_edge_length;
  	}
  	subpartitions_per_partition_dimension_squared = powu(subpartitions_per_partition_dimension, 2);
  }

  void init(bimolecular_reactions_map_t* bimolecular_reactions_map_) {
  	bimolecular_reactions_map = bimolecular_reactions_map_;
  	init_subpartition_edge_length();
  }

  void dump();
};


} /* namespace mcell */

#endif // SRC4_DEFINES_H_
