#ifndef MCELL_DIFFUSE
#define MCELL_DIFFUSE

#include "mcell_structs.h"

void pick_displacement(struct vector3 *v,double scale);
struct collision* ray_trace(struct molecule *m, struct collision *c,
                            struct subvolume *sv, struct vector3 *v);
struct molecule* diffuse_3D(struct molecule *m,double max_time,int inert,double newbie_time);
void run_timestep(struct storage *local,double release_time,double checkpt_time);


#endif
