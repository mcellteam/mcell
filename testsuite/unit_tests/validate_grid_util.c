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

void init_wall(struct wall *w,struct vector3 *v0,struct vector3 *v1,struct vector3 *v2)
{
  double f,fx,fy,fz;
  struct vector3 vA,vB,vX;

  w->vert[0] = v0;
  w->vert[1] = v1;
  w->vert[2] = v2;
  w->edges[0] = NULL;
  w->edges[1] = NULL;
  w->edges[2] = NULL;
  w->nb_walls[0] = NULL;
  w->nb_walls[1] = NULL;
  w->nb_walls[2] = NULL;
  
  vectorize(v0, v1, &vA);
  vectorize(v0, v2, &vB);
  cross_prod(&vA , &vB , &vX);
  w->area = 0.5 * vect_length(&vX);

/*  
  w->area = sqrt(0.5 * ( ((v1->x - v0->x)*(v1->x - v0->x)+
                          (v1->y - v0->y)*(v1->y - v0->y)+
                          (v1->z - v0->z)*(v1->z - v0->z)) *
                         ((v2->x - v0->x)*(v2->x - v0->x)+
                          (v2->y - v0->y)*(v2->y - v0->y)+
                          (v2->z - v0->z)*(v2->z - v0->z)) +
                         ((v1->x - v0->x)*(v2->x - v0->x)+
                          (v1->y - v0->y)*(v2->y - v0->y)+
                          (v1->z - v0->z)*(v2->z - v0->z)) *
                         ((v1->x - v0->x)*(v2->x - v0->x)+
                          (v1->y - v0->y)*(v2->y - v0->y)+
                          (v1->z - v0->z)*(v2->z - v0->z)) ) );
*/

  fx = (v1->x - v0->x);
  fy = (v1->y - v0->y);
  fz = (v1->z - v0->z);
  f = 1 / sqrt( fx*fx + fy*fy + fz*fz );
  
  w->unit_u.x = fx * f;
  w->unit_u.y = fy * f;
  w->unit_u.z = fz * f;
  
  fx = (v2->x - v0->x);
  fy = (v2->y - v0->y);
  fz = (v2->z - v0->z);

  w->normal.x = w->unit_u.y * fz - w->unit_u.z * fy;
  w->normal.y = w->unit_u.z * fx - w->unit_u.x * fz;
  w->normal.z = w->unit_u.x * fy - w->unit_u.y * fx;
  f = 1 / sqrt( w->normal.x*w->normal.x + w->normal.y*w->normal.y + w->normal.z*w->normal.z );
  w->normal.x *= f;
  w->normal.y *= f;
  w->normal.z *= f;
  w->unit_v.x = w->normal.y * w->unit_u.z - w->normal.z * w->unit_u.y;
  w->unit_v.y = w->normal.z * w->unit_u.x - w->normal.x * w->unit_u.z;
  w->unit_v.z = w->normal.x * w->unit_u.y - w->normal.y * w->unit_u.x;
  w->d = v0->x * w->normal.x + v0->y * w->normal.y + v0->z * w->normal.z;
  
  w->uv_vert1_u = (w->vert[1]->x - w->vert[0]->x)*w->unit_u.x + 
                  (w->vert[1]->y - w->vert[0]->y)*w->unit_u.y +
                  (w->vert[1]->z - w->vert[0]->z)*w->unit_u.z;
  w->uv_vert2.u = (w->vert[2]->x - w->vert[0]->x)*w->unit_u.x + 
                  (w->vert[2]->y - w->vert[0]->y)*w->unit_u.y +
                  (w->vert[2]->z - w->vert[0]->z)*w->unit_u.z;
  w->uv_vert2.v = (w->vert[2]->x - w->vert[0]->x)*w->unit_v.x + 
                  (w->vert[2]->y - w->vert[0]->y)*w->unit_v.y +
                  (w->vert[2]->z - w->vert[0]->z)*w->unit_v.z;
  
  w->mol = NULL;
  w->mol_count = 0;
  w->effectors = NULL;
  w->viz_state = EXCLUDE_OBJ; 
}

int main()
{
  int seed = 1;
  struct surface_grid g;
  struct wall w;
  struct vector3 v,vcent,wv0,wv1,wv2;
  int i,idx;
  
  ran4_init(&seed);
  
  g.n = 5;
  g.surface = &w;
  
  wv0.x = 0;
  wv0.y = 0;
  wv0.z = 0;
  wv1.x = 3;
  wv1.y = 0;
  wv1.z = 0;
  wv2.x = -0.5;
  wv2.y = 1;
  wv2.z = 0;
  
  init_wall(&w,&wv0,&wv1,&wv2);
  
  init_grid_geometry(&g);

  for (i=0;i<10000;i++)
  {
    v.x = 3.5*rng_double(3*i)-0.5;
    v.y = 1.0*rng_double(3*i+1);
    v.z = 2.0*rng_double(3*i+2)-1.0;
    
    if (v.y > 0 && 3.5*v.y + v.x < 3.0 && v.y+2.0*v.x > 0.0)
    {
      idx = xyz2grid(&v , &g);
      grid2xyz(&g , idx , &vcent);
      
      printf("%6.4f %6.4f %6.4f %3d %6.4f %6.4f %6.4f\n",
             v.x,v.y,v.z,idx,vcent.x,vcent.y,vcent.z);
    }
  }
}
