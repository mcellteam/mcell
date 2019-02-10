#ifndef SRC4_DEFINES_H_
#define SRC4_DEFINES_H_

#include <stdint.h>

namespace mcell {

typedef double float_t; // soon to be changed to float

struct vec3_t {
	float_t x;
	float_t y;
	float_t z;
};

} /* namespace mcell */

#endif // SRC4_DEFINES_H_
