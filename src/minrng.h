#ifndef MINRNG_H
#define MINRNG_H

#include <inttypes.h>

typedef uint32_t  ub4;

#define rot(x,k) (((x)<<(k))|((x)>>(32-(k))))

struct mrng_state
{
  ub4 a;
  ub4 b;
  ub4 c;
  ub4 d;
};

ub4 mrng_generate( struct mrng_state *x );
void mrng_init( struct mrng_state *x, ub4 seed );
#define mrng_uint32(rng) \
      ( mrng_generate(rng) )

#define mrng_dbl32(rng) \
   ( DBL32 * (double) mrng_uint32(rng) )

#endif
