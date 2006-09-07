#include "rng.h"
#include <math.h>
#include "mcell_structs.h"

extern struct volume *world; 

double rng_gauss(struct rng_state *rng)
{
  double r,v1,v2;

  v1 = rng_dbl(rng) - 1.0;
  if(world->notify->final_summary == NOTIFY_FULL){
      world->random_number_use++;
  }
  v2 = rng_dbl(rng) - 1.0;
  if(world->notify->final_summary == NOTIFY_FULL){
      world->random_number_use++;
  }
  r = v1*v1 + v2*v2;
  
  return v1*sqrt(-2.0*log(r)/r);
}

