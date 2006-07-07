#ifndef MCELL_COUNT_UTIL
#define MCELL_COUNT_UTIL

#include "mcell_structs.h"

void update_collision_count(struct species *sp,struct region_list *rl,int direction,int crossed,double factor,struct vector3 *loc,double t);
int region_listed(struct region_list *rl,struct region *r);
void count_me_by_region(struct abstract_molecule *me,int n,struct rxn_pathname *rxp,double t);
int place_waypoints();
int prepare_counters();
int expand_object_output(struct output_request *request,struct object *obj);
int object_has_geometry(struct object *obj);
int instantiate_request(struct output_request *request);
struct counter* create_new_counter(struct region *where,void *who,byte what);
void fire_count_event(struct counter *event,int n,struct grid_molecule *g,struct vector3 *where,byte what);


#endif
