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
#include <string>
#include <bitset>
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>

#include "mcell_structs.h"

// warning: do not use directly, we need to be able to control the precision
#include "../libs/glm/glm.hpp"


// diverse debug macros
//#define DEBUG_SCHEDULER


namespace mcell {

typedef double float_t; // soon to be changed to float
#define FLOAT_T_BYTES 8

const float_t TIME_INVALID = NAN;
const float_t TIME_SIMULATION_START = 0;

const float_t PARTITION_EDGE_LENGTH_DEFAULT = 10; // maybe too large for now

const float_t SCHEDULER_COMPARISON_EPS = 1e-10;

struct vec3_t: public glm::dvec3 {

	vec3_t() : glm::dvec3(0) {}
	vec3_t(const glm::dvec3& a) { x = a.x; y = a.y; z = a.z; }
	vec3_t(const vec3_t& a) : glm::dvec3(a.x, a.y, a.z) { /*x = a.x; y = a.y; z = a.z;*/ }
	vec3_t(const vector3& a) { x = a.x; y = a.y; z = a.z; }
	vec3_t(const float_t x_, const float_t y_, const float_t z_) { x = x_; y = y_; z = z_; }
	vec3_t(const float_t xyz) { x = xyz; y = xyz; z = xyz; }

	// exact match
	//bool operator == (const vec3_t& a) { return x == a.x && y == a.y && z == a.z; }
};

std::ostream & operator<<(std::ostream &out, const vec3_t &a);


const int MAX_MOLECULES_PER_PARTITION = 32*32*32 /*32k*/; //temporary, must work dynamically

//typedef std::bitset<MAX_MOLECULES_PER_PARTITION> subpartition_mask_t;


typedef uint16_t species_id_t;
const int SPECIES_ID_INVALID = USHRT_MAX;

typedef uint32_t molecule_index_t;

// for now, this is the partition that contains point 0,0,0
// in its center
const uint32_t PARTITION_INDEX_INITIAL = 0;

const uint32_t PARTITION_INDEX_INVALID = UINT32_MAX;

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

} /* namespace mcell */

#endif // SRC4_DEFINES_H_
