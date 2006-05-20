#ifndef MCELL_DIFFUSE
#define MCELL_DIFFUSE

#include "mcell_structs.h"

void pick_displacement(struct vector3 *v,double scale);
void pick_clamped_displacement(struct vector3 *v,struct molecule *m);
struct wall* ray_trace_2d(struct grid_molecule *g, struct vector2 *disp, struct vector2 *loc);
struct collision* ray_trace(struct molecule *m, struct collision *c,
                            struct subvolume *sv, struct vector3 *v,
			    struct wall *reflectee);
struct molecule* diffuse_3D(struct molecule *m,double max_time,int inert);
struct grid_molecule* diffuse_2D(struct grid_molecule *g,double max_time);
struct grid_molecule* react_2D(struct grid_molecule *g,double t);
void run_timestep(struct storage *local,double release_time,double checkpt_time);
void run_concentration_clamp(double t_now);


#endif
