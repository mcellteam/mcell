#ifndef MCELL_VOL_UTIL
#define MCELL_VOL_UTIL

#include "mcell_structs.h"

int bisect(double *list,int n,double val);
int bisect_near(double *list,int n,double val);

int inside_subvolume(struct vector3 *point,struct subvolume *subvol);
struct subvolume* find_course_subvol(struct vector3 *loc);
struct subvolume* traverse_subvol(struct subvolume *here,struct vector3 *point,int which);
struct subvolume* next_subvol(struct vector3 *here,struct vector3 *move,struct subvolume *sv);
struct subvolume* find_subvolume(struct vector3 *loc,struct subvolume *guess);

struct molecule* insert_molecule(struct molecule *m,struct molecule *guess);
void excert_molecule(struct molecule *m);
void insert_molecule_list(struct molecule *m);
struct molecule* migrate_molecule(struct molecule *m,struct subvolume *new_sv);

void release_molecules(struct release_event_queue *req);

void set_partitions();

#endif
