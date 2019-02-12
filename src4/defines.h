#ifndef SRC4_DEFINES_H_
#define SRC4_DEFINES_H_

#include <stdint.h>

#include <vector>
#include <string>
#include <bitset>
#include <cassert>

namespace mcell {

typedef double float_t; // soon to be changed to float

struct vec3_t {
	float_t x;
	float_t y;
	float_t z;
};

const int MAX_MOLECULES_PER_PARTITION = 32*32*32 /*32k*/; //temporary, must work dynamically

typedef std::bitset<MAX_MOLECULES_PER_PARTITION> subpartition_mask_t;

} /* namespace mcell */

#endif // SRC4_DEFINES_H_
