#ifndef MCELL_COUNT_UTIL
#define MCELL_COUNT_UTIL

#include "mcell_structs.h"

void count_hit(struct molecule *m,struct wall *w,int direction);
void count_react(struct rxn *rx,int path,double timestep);
void count_crossings(struct molecule *m,struct subvolume *sv,
                     struct vector3 *move);

#endif
