#ifndef MCELL_VOL_UTIL
#define MCELL_VOL_UTIL

#include "mcell_structs.h"

int bisect(double *list,int n,double val);
int inside_subvolume(struct vector3 *point,struct subvolume *subvol);
struct subvolume* find_course_subvol(struct vector3 *loc);
struct subvolume* traverse_subvol(struct subvolume *here,struct vector3 *point,int which);
struct subvolume* find_subvolume(struct vector3 *loc,struct subvolume *guess);

struct molecule* insert_molecule(struct molecule *m,struct molecule *guess);
void excert_molecule(struct molecule *m);
void insert_molecule_list(struct molecule *m);

#endif
