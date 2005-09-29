#ifndef MCELL_GRID
#define MCELL_GRID

#include "mcell_structs.h"

void xyz2uv(struct vector3 *a,struct wall *w,struct vector2 *b);
void uv2xyz(struct vector2 *a,struct wall *w,struct vector3 *b);

int xyz2grid(struct vector3 *v,struct surface_grid *g);
void grid2xyz(struct surface_grid *g,int index,struct vector3 *v);

int uv2grid(struct vector2 *v,struct surface_grid *g);
void grid2uv(struct surface_grid *g,int index,struct vector2 *v);
void grid2uv_random(struct surface_grid *g,int index,struct vector2 *v);

void init_grid_geometry(struct surface_grid *g);
int create_grid(struct wall *w,struct subvolume *guess);

int nearest_free(struct surface_grid *g,struct vector2 *v,double max_d2,double *found_dist2);
struct wall *search_nbhd_for_free(struct wall *origin,struct vector2 *point,double max_d2,int *found_idx,
                             int (*ok)(void*,struct wall*),void *context);

#endif
