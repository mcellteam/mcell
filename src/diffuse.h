#ifndef MCELL_DIFFUSE
#define MCELL_DIFFUSE

#include "mcell_structs.h"

void pick_displacement(struct vector3 *v,double scale);
void pick_2d_displacement(struct vector2 *v, double scale);
void pick_release_displacement(struct vector3 *in_disk,struct vector3 *away,double scale);
void pick_clamped_displacement(struct vector3 *v,struct volume_molecule *m);
struct wall* ray_trace_2d(struct grid_molecule *g, struct vector2 *disp, struct vector2 *loc);
struct collision* ray_trace(struct volume_molecule *m, struct collision *c,
                            struct subvolume *sv, struct vector3 *v,
			    struct wall *reflectee);
struct sp_collision* ray_trace_trimol(struct volume_molecule *m, 
                            struct sp_collision *c,
                            struct subvolume *sv, struct vector3 *v,
			    struct wall *reflectee, double walk_start_time);
struct volume_molecule* diffuse_3D(struct volume_molecule *m,double max_time,int inert);
struct volume_molecule* diffuse_3D_big_list(struct volume_molecule *m,double max_time,int inert);
struct grid_molecule* diffuse_2D(struct grid_molecule *g,double max_time, double *advance_time);
struct grid_molecule* react_2D(struct grid_molecule *g,double t);
struct grid_molecule* react_2D_all_neighbors(struct grid_molecule *g,double t);
struct grid_molecule* react_2D_trimol_all_neighbors(struct grid_molecule *g,double t);
void run_timestep(struct storage *local,double release_time,double checkpt_time);
void run_concentration_clamp(double t_now);


#endif
