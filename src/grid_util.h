#ifndef MCELL_GRID
#define MCELL_GRID

#include "mcell_structs.h"

int xyz2grid(struct vector3 *v,struct surface_grid *g);
void grid2xyz(struct surface_grid *g,int index,struct vector3 *v);
void init_grid_geometry(struct surface_grid *g);

#endif
