/* To compile: gcc -lm validate_grid_util.c grid_util.c rng.c */
/* To run: ./a.out */
/* Generates a bunch of random points inside a triangle, finds
   their index, then finds the center of that triangle.  If you
   plot starting and ending location in Matlab, you can visually
   check to make sure it's working right.
   
   Output format is xi yi zi idx xf yf zf where xi yi zi are the
   location of the random point, idx is the index into the grid,
   and xf yf zf are the coordinates of the center of that tile.
*/


#include <stdio.h>

#include "rng.h"
#include "mcell_structs.h"
#include "grid_util.h"

int main()
{
  int seed = 1;
  struct surface_grid g;
  struct wall w;
  struct vector3 v,vcent;
  int i,idx;
  
  ran4_init(&seed);
  
  g.n = 5;
  g.surface = &w;
  
  w.uv_vert1_u = 3.0;
  w.uv_vert2.u = 1.0;
  w.uv_vert2.v = 1.0;
  
  w.unit_u.x = 1.0;
  w.unit_u.y = 0.0;
  w.unit_u.z = 0.0;
  
  w.unit_v.x = 0.0;
  w.unit_v.y = 1.0;
  w.unit_v.z = 0.0;
  
  init_grid_geometry(&g);

  for (i=0;i<10000;i++)
  {
    v.x = 3.0*rng_double(3*i);
    v.y = 1.0*rng_double(3*i+1);
    v.z = 2.0*rng_double(3*i+2)-1.0;
    
    if (v.x > v.y && v.x < 3.0 - 2.0*v.y && v.x > 0.0 && v.y > 0.0)
    {
      idx = xyz2grid(&v , &g);
      grid2xyz(&g , idx , &vcent);
      
      printf("%6.4f %6.4f %6.4f %3d %6.4f %6.4f %6.4f\n",
             v.x,v.y,v.z,idx,vcent.x,vcent.y,vcent.z);
    }
  }
}
