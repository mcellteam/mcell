#ifndef MCELL_COUNT_UTIL
#define MCELL_COUNT_UTIL

#include "mcell_structs.h"

int region_listed(struct region_list *rl,struct region *r);
void count_region_update(struct species *sp,
                         struct region_list *rl,
                         int direction,
                         int crossed,
                         struct vector3 *loc,
                         double t);
void count_region_border_update(struct species *sp, struct hit_data *hd_info);
void count_region_from_scratch(struct abstract_molecule *am,struct rxn_pathname *rxpn,int n,struct vector3 *loc,struct wall *my_wall,double t);
void count_moved_grid_mol(struct grid_molecule *g,struct surface_grid *sg,struct vector2 *loc);
void fire_count_event(struct counter *event,int n,struct vector3 *where,byte what);

int place_waypoints(void);
int prepare_counters(void);
int check_counter_geometry(void);
int expand_object_output(struct output_request *request,struct object *obj);
int object_has_geometry(struct object *obj);

/************************************************************
 * Complex counting
 ************************************************************/
int count_complex(struct volume_molecule *cmplex,
                  struct volume_molecule *replaced_subunit,
                  int replaced_subunit_idx);
int count_complex_surface(struct grid_molecule *cmplex,
                          struct grid_molecule *replaced_subunit,
                          int replaced_subunit_idx);
int count_complex_surface_new(struct grid_molecule *cmplex);

#endif
