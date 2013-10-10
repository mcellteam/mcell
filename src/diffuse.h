#ifndef MCELL_DIFFUSE
#define MCELL_DIFFUSE

#include "mcell_structs.h"

void pick_displacement(struct vector3 *v,double scale, 
    struct rng_state *rng);

void pick_2d_displacement(struct vector2 *v, double scale, 
    struct rng_state *rng);

void pick_release_displacement(struct vector3 *in_disk, struct vector3 *away,
    double scale, double *r_step_release, double *d_step, 
    u_int radial_subdivisions, int directions_mask, u_int num_directions, 
    double rx_radius_3d, struct rng_state *rng);

void pick_clamped_displacement(struct vector3 *v,struct volume_molecule *m,
    double *r_step_surfce, struct rng_state *rng, u_int radial_subdivision);

struct wall* ray_trace_2d(struct grid_molecule *g, struct vector2 *disp, 
    struct vector2 *loc, int *kill_me, struct rxn **rxp, 
    struct hit_data **hd_info, struct vector3 *all_vertices);

struct collision* ray_trace(struct volume *world, struct volume_molecule *m, 
    struct collision *c, struct subvolume *sv, struct vector3 *v, 
    struct wall *reflectee);

struct sp_collision* ray_trace_trimol(struct volume *world, 
    struct volume_molecule *m, struct sp_collision *c, struct subvolume *sv, 
    struct vector3 *v, struct wall *reflectee, double walk_start_time);

struct volume_molecule* diffuse_3D(struct volume *world, 
    struct volume_molecule *m, double max_time, int inert);

struct volume_molecule* diffuse_3D_big_list(struct volume *world,
    struct volume_molecule *m, double max_time,int inert);

struct grid_molecule* diffuse_2D(struct volume *world, 
    struct grid_molecule *g, double max_time, double *advance_time);

struct grid_molecule* react_2D(struct volume *world, struct grid_molecule *g, 
    double t, enum notify_level_t molecule_collision_report, 
    int grid_grid_reaction_flag, long long *grid_grid_colls);

struct grid_molecule* react_2D_all_neighbors(struct volume *world,
    struct grid_molecule *g, double t, 
    enum notify_level_t molecule_collision_report, 
    int grid_grid_reaction_flag, long long *grid_grid_colls);

struct grid_molecule* react_2D_trimol_all_neighbors(struct volume *world,
    struct grid_molecule *g, double t, 
    enum notify_level_t molecule_collision_report, 
    enum notify_level_t final_summary, int grid_grid_reaction_flag, 
    long long *grid_grid_colls);

void run_timestep(struct volume *world, struct storage *local,
    double release_time, double checkpt_time);

void run_concentration_clamp(struct volume *world, double t_now);


#endif
