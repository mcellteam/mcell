#include "rng.h"
#include <math.h>

double rng_gauss(struct rng_state *rng)
{
  double r,v1,v2;

  v1 = rng_dbl(rng) - 1.0;
  v2 = rng_dbl(rng) - 1.0;
  r = v1*v1 + v2*v2;
  
  return v1*sqrt(-2.0*log(r)/r);
}

