#ifndef MCELL_DIFFUSE
#define MCELL_DIFFUSE

#include "mcell_structs.h"

void pick_displacement(struct vector3 *v,double scale);
struct collision* ray_trace(struct molecule *m, struct collision *c,
                            struct subvolume *sv, struct vector3 *v);
int diffuse_3D(struct molecule *m,double target_time);


#endif