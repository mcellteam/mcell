#ifndef RAN4_H
#define RAN4_H

struct ran4_state
{
  unsigned int i_state;
  unsigned int e_state;
};

void ran4_init(struct ran4_state *rng,int seed);
unsigned int ran4_uint32(struct ran4_state *rng);
double ran4_dbl32(struct ran4_state *rng);

#endif
