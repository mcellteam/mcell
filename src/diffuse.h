#ifndef MCELL_DIFFUSE
#define MCELL_DIFFUSE

#include "mcell_structs.h"

struct wall* ray_trace_2d(struct grid_molecule *g, struct vector2 *disp, struct vector2 *loc);
struct collision* ray_trace(struct storage *local,
                            struct volume_molecule *m,
                            struct collision *c,
                            struct subvolume *sv,
                            struct vector3 *v,
                            struct wall *reflectee);
struct sp_collision* ray_trace_trimol(struct storage *local,
                                      struct volume_molecule *m, 
                                      struct sp_collision *c,
                                      struct subvolume *sv,
                                      struct vector3 *v,
                                      struct wall *reflectee,
                                      double walk_start_time);
void run_timestep(struct storage *local, double release_time,double checkpt_time);
void run_concentration_clamp(double t_now);


#endif
