#ifndef MCELL_COUNT_UTIL
#define MCELL_COUNT_UTIL

#include "mcell_structs.h"

int region_listed(struct region_list *rl,struct region *r);

void count_region_update(struct volume *world, struct species *sp,
  struct region_list *rl, int direction, int crossed, struct vector3 *loc,
  double t);

void count_region_border_update(struct species *sp, struct hit_data *hd_info,
    int count_hashmask, struct counter **count_hash);

void count_region_from_scratch(struct volume *world, 
    struct abstract_molecule *am, struct rxn_pathname *rxpn, int n,
    struct vector3 *loc,struct wall *my_wall,double t);

void count_moved_grid_mol(struct volume *world, struct grid_molecule *g,
    struct surface_grid *sg, struct vector2 *loc, int count_hashmask, 
    struct counter **count_hash, long long *ray_polygon_colls); 

void fire_count_event(struct counter *event,int n,struct vector3 *where,
    byte what);

int place_waypoints(struct volume *world);

int prepare_counters(struct volume *world); 

int check_counter_geometry(int count_hashmask, struct counter **count_hash,
    byte *place_waypoints_flag);

int expand_object_output(struct output_request *request,struct object *obj);
int object_has_geometry(struct object *obj);

/* hit data for region borders */
void update_hit_data(struct hit_data** hd_head,
                     struct wall* current,
                     struct wall* target,
                     struct grid_molecule *g,
                     struct vector2 boundary_pos,
                     int direction,
                     int crossed);

/************************************************************
 * Complex counting
 ************************************************************/
int count_complex(struct volume *world,
                  struct volume_molecule *cmplex,
                  struct volume_molecule *replaced_subunit,
                  int replaced_subunit_idx);
int count_complex_surface(struct grid_molecule *cmplex,
                          struct grid_molecule *replaced_subunit,
                          int replaced_subunit_idx);
int count_complex_surface_new(struct grid_molecule *cmplex);

#endif
