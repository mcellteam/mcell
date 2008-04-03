#include "minrng.h"

ub4 mrng_generate( struct mrng_state *x )
{
  ub4 e = x->a - rot(x->b, 27);
  x->a = x->b ^ rot(x->c, 17);
  x->b = x->c + x->d;
  x->c = x->d + e;
  x->d = e + x->a;
  return x->d;
}

void mrng_init( struct mrng_state *x, ub4 seed )
{
  ub4 i;
  x->a = 0xf1ea5eed, x->b = x->c = x->d = seed;
  for (i=0; i<20; ++i) {
    (void) mrng_generate(x);
  }
}
