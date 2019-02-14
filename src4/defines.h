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

namespace mcell {

typedef double float_t; // soon to be changed to float
const float_t TIME_INVALID = NAN;

struct vec3_t {
	float_t x;
	float_t y;
	float_t z;

	vec3_t() : x(0), y(0), z(0) {}
	vec3_t(const vec3_t& a) { x = a.x; y = a.y; z = a.z; }
	vec3_t(const vector3& a) { x = a.x; y = a.y; z = a.z; }
	vec3_t(const float_t x_, const float_t y_, const float_t z_) { x = x_; y = y_; z = z_; }

	// exact match
	bool operator == (const vec3_t& a) { return x == a.x && y == a.y && z == a.z; }
};

std::ostream & operator<<(std::ostream &out, const vec3_t &a);


const int MAX_MOLECULES_PER_PARTITION = 32*32*32 /*32k*/; //temporary, must work dynamically

typedef std::bitset<MAX_MOLECULES_PER_PARTITION> subpartition_mask_t;


typedef uint16_t species_id_t;
const int SPECIES_ID_INVALID = USHRT_MAX;


typedef glm::dmat4x4 mat4x4;

} /* namespace mcell */

#endif // SRC4_DEFINES_H_
