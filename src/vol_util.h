#ifndef MCELL_VOL_UTIL
#define MCELL_VOL_UTIL

#include "mcell_structs.h"

int inside_subvolume(struct vector3 *point,struct subvolume *subvol);
struct subvolume* find_course_subvol(struct vector3 *loc);
struct subvolume* traverse_subvol(struct subvolume *here,struct vector3 *point,int which);
struct subvolume* next_subvol(struct vector3 *here,struct vector3 *move,struct subvolume *sv);
struct subvolume* find_subvolume(struct vector3 *loc,struct subvolume *guess);
double collide_sv_time(struct vector3 *point,struct vector3 *move,struct subvolume *sv);

struct molecule* insert_molecule(struct molecule *m,struct molecule *guess);
void excert_molecule(struct molecule *m);
int insert_molecule_list(struct molecule *m);
struct molecule* migrate_molecule(struct molecule *m,struct subvolume *new_sv);

int release_inside_regions(struct release_site_obj *rso,struct molecule *m,int n);
int release_molecules(struct release_event_queue *req);

int set_partitions();
double distance_point_line(struct vector3 *q, struct vector3 *v0, struct vector3 *v1);
int navigate_world(int curr_index, int direction);
int navigate_world_by_edge(int curr_index, int direction1, int direction2);
int navigate_world_by_corner(int curr_index, int direction1, int direction2, int direction3);

void path_bounding_box(struct vector3 *loc, struct vector3 *displacement, struct vector3 *llf, struct vector3 *urb);
#endif
