/**************************************************************************\
** File: diffuse.c                                                        **
**                                                                        **
** Purpose: Moves molecules around the world with reactions and collisions**
**                                                                        **
** Testing status: compiles.                                              **
\**************************************************************************/



#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"
#include "mem_util.h"
#include "sched_util.h"
#include "util.h"

#include "mcell_structs.h"
#include "count_util.h"
#include "grid_util.h"
#include "vol_util.h"
#include "wall_util.h"
#include "react.h"
#include "react_output.h"


#define MULTISTEP_WORTHWHILE 2
#define MULTISTEP_PERCENTILE 0.99
#define MULTISTEP_FRACTION 0.9
#define MAX_UNI_TIMESKIP 5000


/* EXD_TIME_CALC is a local #define in exact_disk */
/* EXD_SPAN_CALC is a local #define in exact_disk */

/* CLEAN_AND_RETURN(x) is a local #define in diffuse_3D */
/* ERROR_AND_QUIT is a local #define in diffuse_3D */


extern struct volume *world;
extern double scaling_sum, scaling_count;



/*************************************************************************
pick_2d_displacement:
  In: vector2 to store the new displacement
      scale factor to apply to the displacement
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         2D molecule, scaled by the scaling factor.
*************************************************************************/

void pick_2d_displacement(struct vector2 *v,double scale)
{
  struct vector2 a;
  double f;
  
  /* TODO: this is exact case only--see if lookup is faster. */
  do
  {
    a.u = 2.0*rng_dbl(world->rng)-1.0;
    a.v = 2.0*rng_dbl(world->rng)-1.0;
    f = a.u*a.u + a.v*a.v;
  } while (f<0.01 || f>1.0);
  
  f = (1.0/f) * sqrt(log( 1/(1-rng_dbl(world->rng)) )) * scale;
  
  v->u = (a.u*a.u-a.v*a.v)*f;
  v->v = (2.0*a.u*a.v)*f;
  
//  printf("Disp %.8f %.8f %.8f\n",v->u,v->v,scale);
}


/*************************************************************************
pick_clamped_displacement:
  In: vector3 to store the new displacement
      molecule that just came through the surface
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         3D molecule that has come through a surface from a uniform
         concentration on the other side.
  Note: m->previous_wall points to the wall we're coming from, and
        m->index is the orientation we came off with
*************************************************************************/

void pick_clamped_displacement(struct vector3 *v,struct molecule *m)
{
  double p;
  double r_n;
  struct vector2 r_uv;
  struct wall *w = m->previous_wall;
  
  p = rng_dbl(world->rng);
  
  /* Correct distribution along normal from surface (from lookup table) */
  r_n = m->index * m->properties->space_step * 
        world->r_step_surface[ rng_uint(world->rng) & (world->radial_subdivisions-1) ];
  
  /* This is just a guess at an appropriate approximate planar distribution */
  pick_2d_displacement(&r_uv,sqrt(p)*m->properties->space_step); 
  
  v->x = r_n*w->normal.x + r_uv.u*w->unit_u.x + r_uv.v*w->unit_v.x;
  v->y = r_n*w->normal.y + r_uv.u*w->unit_u.y + r_uv.v*w->unit_v.y;
  v->z = r_n*w->normal.z + r_uv.u*w->unit_u.z + r_uv.v*w->unit_v.z;
}



/*************************************************************************
pick_displacement:
  In: vector3 to store the new displacement
      scale factor to apply to the displacement
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         3D molecule, scaled by the scaling factor.
*************************************************************************/

void pick_displacement(struct vector3 *v,double scale)
{
  double r;
  
  if (world->fully_random)
  {
    double h,r_sin_phi,theta;
    
    r = scale * world->r_step[ rng_uint(world->rng) & (world->radial_subdivisions-1) ];
    h = 2.0*rng_dbl(world->rng) - 1.0;
    theta = 2.0*MY_PI*rng_dbl(world->rng);
    
    r_sin_phi = r * sqrt(1.0 - h*h);
    
    v->x = r_sin_phi * cos(theta);
    v->y = r_sin_phi * sin(theta);
    v->z = r*h;
  }
  else
  {
    u_int x_bit,y_bit,z_bit;
    u_int r_bit,thetaphi_bit;
    u_int bits;
    u_int idx;
    
    bits = rng_uint(world->rng);
    
    x_bit =        (bits & 0x80000000);
    y_bit =        (bits & 0x40000000);
    z_bit =        (bits & 0x20000000);
    thetaphi_bit = (bits & 0x1FFFF000) >> 12;
    r_bit =        (bits & 0x00000FFF);
    
    r = scale * world->r_step[ r_bit & (world->radial_subdivisions - 1) ];
    
    idx = (thetaphi_bit & world->directions_mask);
    while ( idx >= world->num_directions)
    {
      idx = ( rng_uint(world->rng) & world->directions_mask);
    }
    
    idx *= 3;
    
    if (x_bit) v->x = r * world->d_step[idx];
    else v->x = -r * world->d_step[idx];
    
    if (y_bit) v->y = r * world->d_step[idx+1];
    else v->y = -r * world->d_step[idx+1];

    if (z_bit) v->z = r * world->d_step[idx+2];
    else v->z = -r * world->d_step[idx+2];
  }
}


/*************************************************************************
ray_trace_2d:
  In: molecule that is moving
      displacement vector from current to new location
      place to store new coordinate (in coord system of new wall)
  Out: wall at endpoint of movement vector, plus location of that endpoint
       in the coordinate system of the new wall.
*************************************************************************/

struct wall* ray_trace_2d(struct grid_molecule *g,struct vector2 *disp,struct vector2 *pos)
{
  struct vector2 first_pos,old_pos,boundary_pos;
  struct vector2 this_pos,this_disp;
  struct vector2 new_disp,reflector;
  struct wall *this_wall,*target_wall;
  struct rxn *rx;
  int new_wall_index;
  double f;
  
  this_wall = g->grid->surface;
  
  first_pos.u = g->s_pos.u;
  first_pos.v = g->s_pos.v;
  
  this_pos.u = g->s_pos.u;
  this_pos.v = g->s_pos.v;
  this_disp.u = disp->u;
  this_disp.v = disp->v;
  
  while (1) /* Will break out with return or break when we're done traversing walls */
  {
    new_wall_index = find_edge_point(this_wall,&this_pos,&this_disp,&boundary_pos);
    
    if (new_wall_index==-2) /* Ambiguous edge collision--just give up */
    {
      g->s_pos.u = first_pos.u;
      g->s_pos.v = first_pos.v;
      return NULL;
    }
  
    if (new_wall_index==-1) /* We didn't hit the edge.  Stay inside this wall. */
    {
      pos->u = this_pos.u + this_disp.u;
      pos->v = this_pos.v + this_disp.v;
      
      g->s_pos.u = first_pos.u;
      g->s_pos.v = first_pos.v;
      return this_wall;
    }
    
    old_pos.u = this_pos.u;
    old_pos.v = this_pos.v;
    target_wall = traverse_surface(this_wall,&old_pos,new_wall_index,&this_pos);
  
    if (target_wall!=NULL)
    {
      rx = trigger_intersect(g->properties->hashval,(struct abstract_molecule*)g,g->orient,target_wall);
      if (rx==NULL || rx->n_pathways!=RX_REFLEC)
      {
	this_disp.u = old_pos.u + this_disp.u;
	this_disp.v = old_pos.v + this_disp.v;
	traverse_surface(this_wall,&this_disp,new_wall_index,&new_disp);
	this_disp.u = new_disp.u - this_pos.u;
	this_disp.v = new_disp.v - this_pos.v;
	this_wall = target_wall;
	
	continue;
      }
    }
    
    /* If we reach this point, assume we reflect off edge */
    /* Note that this_pos has been corrupted by traverse_surface; use old_pos */
    new_disp.u = this_disp.u - (boundary_pos.u - old_pos.u);
    new_disp.v = this_disp.v - (boundary_pos.v - old_pos.v);
    switch (new_wall_index)
    {
      case 0:
	new_disp.v *= -1.0;
	break;
      case 1:
	reflector.u = -this_wall->uv_vert2.v;
	reflector.v = this_wall->uv_vert2.u-this_wall->uv_vert1_u;
	f = 1.0/sqrt(reflector.u*reflector.u + reflector.v*reflector.v);
	reflector.u *= f;
	reflector.v *= f;
	f = 2.0 * (new_disp.u*reflector.u + new_disp.v*reflector.v);
	new_disp.u -= f*reflector.u;
	new_disp.v -= f*reflector.v;
	break;
      case 2:
	reflector.u = this_wall->uv_vert2.v;
	reflector.v = -this_wall->uv_vert2.u;
	f = 1.0/sqrt(reflector.u*reflector.u + reflector.v*reflector.v);
	reflector.u *= f;
	reflector.v *= f;
	f = 2.0 * (new_disp.u*reflector.u + new_disp.v*reflector.v);
	new_disp.u -= f*reflector.u;
	new_disp.v -= f*reflector.v;
	break;
    }
    
    this_pos.u = boundary_pos.u;
    this_pos.v = boundary_pos.v;
    this_disp.u = new_disp.u;
    this_disp.v = new_disp.v;
  }
  
  g->s_pos.u = first_pos.u;
  g->s_pos.v = first_pos.v;
  return NULL;
}


/*************************************************************************
ray_trace:
  In: molecule that is moving
      linked list of potential collisions with molecules (we could react)
      subvolume that we start in
      displacement vector from current to new location
      wall we have reflected off of and should not hit again
  Out: collision list of walls and molecules we intersected along our ray
       (current subvolume only), plus the subvolume wall.  Will always
       return at least the subvolume wall--NULL indicates an out of
       memory error.
*************************************************************************/

struct collision* ray_trace(struct molecule *m, struct collision *c,
                            struct subvolume *sv, struct vector3 *v,
			    struct wall *reflectee)
{
  struct collision *smash,*shead;
  struct abstract_molecule *a;
  struct wall_list *wlp;
  struct wall_list fake_wlp;
  double dx,dy,dz;
  /* time, in units of of the molecule's time step, at which molecule
     will cross the x,y,z partitions, respectively. */ 
  double tx,ty,tz;
  int i,j,k;
  
  shead = NULL;
  smash = (struct collision*) mem_get(sv->local_storage->coll);
  if(smash == NULL) return NULL;

  fake_wlp.next = sv->wall_head;
    
  for (wlp = sv->wall_head ; wlp != NULL; wlp = wlp->next)
  {
    if (wlp->this_wall==reflectee) continue;
    
    i = collide_wall(&(m->pos),v,wlp->this_wall,&(smash->t),&(smash->loc));

    if (i==COLLIDE_REDO)
    {
      if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);
      shead = NULL;
      wlp = &fake_wlp;
      continue;
    }
    else if (i!=COLLIDE_MISS)
    {
      smash->what = COLLIDE_WALL + i;
      smash->target = (void*) wlp->this_wall;
      smash->next = shead;
      shead = smash;
      smash = (struct collision*) mem_get(sv->local_storage->coll);
      if (smash==NULL)
      {
	if (shead!=NULL) mem_put_list(sv->local_storage->coll,shead);
	return NULL;
      }
    }
  }

  dx=dy=dz=0.0;
  i=-10;
  if (v->x < 0.0)
  {
    dx = world->x_fineparts[ sv->llf.x ] - m->pos.x;
    i = 0;
  }
  else if (v->x > 0.0)
  {
    dx = world->x_fineparts[ sv->urb.x ] - m->pos.x;
    i = 1;
  }

  j=-10;
  if (v->y < 0.0)
  {
    dy = world->y_fineparts[ sv->llf.y ] - m->pos.y;
    j = 0;
  }
  else if (v->y > 0.0)
  {
    dy = world->y_fineparts[ sv->urb.y ] - m->pos.y;
    j = 1;
  }

  k=-10;
  if (v->z < 0.0)
  {
    dz = world->z_fineparts[ sv->llf.z ] - m->pos.z;
    k = 0;
  }
  else if (v->z > 0.0)
  {
    dz = world->z_fineparts[ sv->urb.z ] - m->pos.z;
    k = 1;
  }
  
  if (i+j+k < -15) /* Two or three vectors are zero */
  {
    if (i >= 0) /* X is the nonzero one */
    {
      smash->t = dx / v->x;
      smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NX + i;
    }
    else if (j>=0) /* Y is nonzero */
    {
      smash->t = dy / v->y;
      smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NY + j;
    }
    else if (k>=0) /* Z is nonzero */
    {
      smash->t = dz / v->z;
      smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NZ + k;
    }
    else
    {
      smash->t = FOREVER;
      smash->what = COLLIDE_SUBVOL; /*Wrong, but we'll never hit it, so it's ok*/
    }
  }
  else /* Zero or one vectors zero--use alternate method */
  {
    tx = dx * v->y * v->z;
    ty = v->x * dy * v->z;
    tz = v->x * v->y * dz;
    if (tx<0.0) tx = -tx;
    if (ty<0.0) ty = -ty;
    if (tz<0.0) tz = -tz;
  
    if (tx < ty)
    {
      if (tx < tz)
      {
        smash->t = dx / v->x;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NX + i;
      }
      else
      {
        smash->t = dz / v->z;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NZ + k;
      }
    }
    else
    {
      if (ty < tz)
      {
        smash->t = dy / v->y;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NY + j;
      }
      else
      {
        smash->t = dz / v->z;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NZ + k;
      }
    }
  }
  
  smash->loc.x = m->pos.x + smash->t * v->x;
  smash->loc.y = m->pos.y + smash->t * v->y;
  smash->loc.z = m->pos.z + smash->t * v->z;
  
  smash->target = sv;
  smash->next = shead;
  shead = smash;

  for ( ; c!=NULL ; c = c->next)
  {
    a = (struct abstract_molecule*)c->target;
    if (a->properties==NULL) continue;
    
    i = collide_mol(&(m->pos),v,a,&(c->t),&(c->loc));
    if (i != COLLIDE_MISS)
    {
      smash = (struct collision*) mem_get(sv->local_storage->coll);
      if (smash==NULL)
      {
	mem_put_list(sv->local_storage->coll,shead);
	return NULL;
      }
      memcpy(smash,c,sizeof(struct collision));
      
      smash->what = COLLIDE_MOL + i;

      smash->next = shead;
      shead = smash;
    }
  }
  
  return shead;
}



/*************************************************************************
estimate_disk:
  In: location of moving molecule at time of collision
      movement vector for moving molecule
      interaction radius
      subvolume the moving molecule is in
      the moving molecule
      the target molecule at time of collision
  Out: The fraction of a full interaction disk that is actually
       accessible to the moving molecule, estimated using Monte Carlo
       integration, or TARGET_OCCLUDED if the target is not accessible.
*************************************************************************/

double estimate_disk(struct vector3 *loc,struct vector3 *mv,double R,struct subvolume *sv,struct molecule *moving,struct molecule *target)
{
  int rpt,idx,bits;
  double area;
  struct vector3 u,v,loc_to_targ;
  double d2_mv_i,a,b,t;
  double upperU;
  double upperV;
  double lowerU;
  double lowerV;
  struct wall_list *wl;
  struct rxn *rx;
  
  area = 0;
  d2_mv_i = 1.0/(mv->x*mv->x + mv->y*mv->y + mv->z*mv->z);
  
  loc_to_targ.x = target->pos.x - loc->x;
  loc_to_targ.y = target->pos.y - loc->y;
  loc_to_targ.z = target->pos.z - loc->z;
  
  for (rpt = 0; rpt < 1 ; rpt++)
  {
  
  upperU = lowerU = upperV = lowerV = 1.0;
  
  do
  {
    bits = rng_uint(world->rng);
    idx = bits & world->directions_mask;
  } while (idx >= world->num_directions);
  
  idx *= 3;
  if (bits&0x80000000) u.x = world->d_step[idx]; else u.x = -world->d_step[idx];
  if (bits&0x40000000) u.y = world->d_step[idx+1]; else u.y = -world->d_step[idx+1];
  if (bits&0x20000000) u.z = world->d_step[idx+2]; else u.z = -world->d_step[idx+2];
  
  a = (u.x*mv->x + u.y*mv->y + u.z*mv->z);
  
  if (a*a*d2_mv_i < 0.9)  /* Vectors too closely aligned */
  {
    rpt--;
    continue;
  }
  
  a *= d2_mv_i;
  u.x = u.x - a*mv->x;
  u.y = u.y - a*mv->y;
  u.z = u.z - a*mv->z;
  b = R/sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
  u.x *= b;
  u.y *= b;
  u.z *= b;
  a = sqrt(d2_mv_i);
  v.x = a*(mv->y*u.z - mv->z*u.y);
  v.y = a*(-mv->x*u.z + mv->z*u.x);
  v.z = a*(mv->x*u.y - mv->y*u.x);
  
  for (wl = sv->wall_head ; wl!=NULL ; wl = wl->next)
  {
    if ( (moving->properties->flags && CAN_MOLWALL) != 0 )
    {
      rx = trigger_intersect(moving->properties->hashval,(struct abstract_molecule*)moving,0,wl->this_wall);
      if (rx != NULL && (rx->n_pathways==RX_TRANSP))
      {
	continue; /* We can move through this wall! */
      }
    }
    
    t = touch_wall(loc,&u,wl->this_wall);
    if (t>0.0 && t<upperU) upperU = t;
    if (t<0.0 && -t<lowerU) lowerU = -t;
    
    t = touch_wall(loc,&v,wl->this_wall);
    if (t>0.0 && t<upperV) upperV = t;
    if (t<0.0 && -t<lowerV) lowerV = -t;
    
    if (rpt==0)
    {
      t = touch_wall(loc,&loc_to_targ,wl->this_wall);
      if (t>0 && t<1) return TARGET_OCCLUDED;  /* This wall blocked us! */
    }
  }

  if (u.x > EPS_C)
  {
    u.x = 1/u.x;
    t = (world->x_fineparts[sv->urb.x] - loc->x)*u.x;
    if (t < upperU) upperU = t;
    t = (loc->x - world->x_fineparts[sv->llf.x])*u.x;
    if (t < lowerU) lowerU = t;
  }
  else if (u.x < -EPS_C)
  {
    u.x = 1/u.x;
    t = (world->x_fineparts[sv->llf.x] - loc->x)*u.x;
    if (t < upperU) upperU = t;
    t = (loc->x - world->x_fineparts[sv->urb.x])*u.x;
    if (t < lowerU) lowerU = t;
  }
  if (u.y > EPS_C)
  {
    u.y = 1/u.y;
    t = (world->y_fineparts[sv->urb.y] - loc->y)*u.y;
    if (t < upperU) upperU = t;
    t = (loc->y - world->y_fineparts[sv->llf.y])*u.y;
    if (t < lowerU) lowerU = t;
  }
  else if (u.y < -EPS_C)
  {
    u.y = 1/u.y;
    t = (world->y_fineparts[sv->llf.y] - loc->y)*u.y;
    if (t < upperU) upperU = t;
    t = (loc->y - world->y_fineparts[sv->urb.y])*u.y;
    if (t < lowerU) lowerU = t;
  }
  if (u.z > EPS_C)
  {
    u.z = 1/u.z;
    t = (world->z_fineparts[sv->urb.z] - loc->z)*u.z;
    if (t < upperU) upperU = t;
    t = (loc->z - world->z_fineparts[sv->llf.z])*u.z;
    if (t < lowerU) lowerU = t;
  }
  else if (u.z < -EPS_C)
  {
    u.z = 1/u.z;
    t = (world->z_fineparts[sv->llf.z] - loc->z)*u.z;
    if (t < upperU) upperU = t;
    t = (loc->z - world->z_fineparts[sv->urb.z])*u.z;
    if (t < lowerU) lowerU = t;
  }

  if (v.x > EPS_C)
  {
    v.x = 1/v.x;
    t = (world->x_fineparts[sv->urb.x] - loc->x)*v.x;
    if (t < upperV) upperV = t;
    t = (loc->x - world->x_fineparts[sv->llf.x])*v.x;
    if (t < lowerV) lowerV = t;
  }
  else if (v.x < -EPS_C)
  {
    v.x = 1/v.x;
    t = (world->x_fineparts[sv->llf.x] - loc->x)*v.x;
    if (t < upperV) upperV = t;
    t = (loc->x - world->x_fineparts[sv->urb.x])*v.x;
    if (t < lowerV) lowerV = t;
  }
  if (v.y > EPS_C)
  {
    v.y = 1/v.y;
    t = (world->y_fineparts[sv->urb.y] - loc->y)*v.y;
    if (t < upperV) upperV = t;
    t = (loc->y - world->y_fineparts[sv->llf.y])*v.y;
    if (t < lowerV) lowerV = t;
  }
  else if (v.y < -EPS_C)
  {
    v.y = 1/v.y;
    t = (world->y_fineparts[sv->llf.y] - loc->y)*v.y;
    if (t < upperV) upperV = t;
    t = (loc->y - world->y_fineparts[sv->urb.y])*v.y;
    if (t < lowerV) lowerV = t;
  }
  if (v.z > EPS_C)
  {
    v.z = 1/v.z;
    t = (world->z_fineparts[sv->urb.z] - loc->z)*v.z;
    if (t < upperV) upperV = t;
    t = (loc->z - world->z_fineparts[sv->llf.z])*v.z;
    if (t < lowerV) lowerV = t;
  }
  else if (v.z < -EPS_C)
  {
    v.z = 1/v.z;
    t = (world->z_fineparts[sv->llf.z] - loc->z)*v.z;
    if (t < upperV) upperV = t;
    t = (loc->z - world->z_fineparts[sv->urb.z])*v.z;
    if (t < lowerV) lowerV = t;
  }

  if (upperU < 0 || upperU > 1.0 ||
      lowerU < 0 || lowerU > 1.0 ||
      upperV < 0 || upperV > 1.0 ||
      lowerV < 0 || lowerV > 1.0)
  {
    fprintf(world->log_file, "File '%s', Line %ld: MCell should not get to this point.  Please report this message.\n", __FILE__, (long)__LINE__);
  }

  area += upperU*upperU + lowerU*lowerU + upperV*upperV + lowerV*lowerV;
/*
  if (a > 1.1) printf("Correction factor %.2f\n",a);
  if (a < 1.0-EPS_C) printf("MUDDY BLURDER! a=%.2f R=%.2f u=[%.2f %.2f %.2f] %.2f %.2f %.2f %.2f\n",
                            a,R,u.x,u.y,u.z,upperU,lowerU,upperV,lowerV);
*/
  }
  
  if (rpt==0) return 1.0;
  return area/(4.0*rpt);
}




/******************************/
/** exact_disk stuff follows **/
/******************************/

/* This is the same as the normal vector3, but the coordinates have different names */
struct exd_vector3
{
  double m,u,v;
};


/****************************************************************
exd_zetize:
In: y coordinate (as in atan2)
    x coordinate (as in atan2)
Out: Zeta value corresponding to (y,x), in the range 0 to 4.
     Zeta is a substitute for the angle theta, and this function
     is a substitute for atan2(y,x) which returns theta. Like
     theta, zeta increases throughout the unit circle, but it
     only has 8-fold symmetry instead of perfect symmetry.  Zeta
     values 0-1 are the first quadrant, 1-2 the second, and so
     on.  Zeta is a monotonically increasing function of theta,
     but requires ~9x less computation time--valuable for when
     you need to sort by angle but don't need the angle itself.
Note: This is a utility finction in 'exact_disk()'.
****************************************************************/
/* Speed: 9ns (compare with 84ns for atan2) */
/* Added extra computations--speed not retested yet */
double exd_zetize(double y,double x)
{
  if (y>=0.0)
  {  
    if (x>=0)
    {
      if (x<y) return 1.0-0.5*x/y;
      else return 0.5*y/x;
    }
    else
    {
      if (-x<y) return 1.0-0.5*x/y;
      else return 2.0+0.5*y/x;
    }
  }
  else
  {
    if (x<=0)
    {
      if (y<x) return 3.0-0.5*x/y;
      else return 2.0+0.5*y/x;
    }
    else
    {
      if (x<-y) return 3.0-0.5*x/y;
      else return 4.0+0.5*y/x;
    }
  }   
}  


/*********************************************************************
exd_coordize:
In: movement vector
    place to store unit movement vector (first basis vector)
    place to store second basis vector
    place to store third basis vector
Out: No return value.  Unit vectors m,u,v are set such that vector m
     is in the direction of vector mv, and vectors u and v are
     orthogonal to m.  The vectors m,u,v, form a right-handed
     coordinate system.
Note: This is a utility function for 'exact_disk()'.
*********************************************************************/ 
/* Speed: 86ns on azzuri (as marked + 6ns function call overhead) */
void exd_coordize(struct vector3 *mv,struct vector3 *m,struct vector3 *u,struct vector3 *v)
{
  double a;
  
  /* Normalize input vector -- 27ns */
  a = 1.0/sqrt(mv->x*mv->x + mv->y*mv->y + mv->z*mv->z);
  m->x = a*mv->x; m->y = a*mv->y; m->z = a*mv->z;
  
  /* Find orthogonal vectors -- 21ns */ 
  if (m->x*m->x > m->y*m->y)
  {
    if (m->x*m->x > m->z*m->z)
    {
      if (m->y*m->y > m->z*m->z)
      {
        u->x = m->y; u->y = -m->x; u->z = 0.0; a = 1.0 - m->z*m->z;
        v->x = m->z*m->x; v->y = m->z*m->y; v->z = -a;
      }
      else
      {   
        u->x = m->z; u->y = 0.0; u->z = -m->x; a = 1.0 - m->y*m->y;
        v->x = -m->y*m->x; v->y = a; v->z = -m->y*m->z;
      }
    }  
    else
    {   
      u->x = -m->z; u->y = 0.0; u->z = m->x; a = 1.0 - m->y*m->y;
      v->x = m->y*m->x; v->y = -a; v->z = m->y*m->z;
    }
  }  
  else
  {   
    if (m->y*m->y > m->z*m->z)
    {
      if (m->x*m->x > m->z*m->z)
      {
        u->x = -m->y; u->y = m->x; u->z = 0.0; a = 1.0 - m->z*m->z;
        v->x = -m->z*m->x; v->y = -m->z*m->y; v->z = a;
      }
      else
      {   
        u->x = 0.0; u->y = m->z; u->z = -m->y; a = 1.0 - m->x*m->x;
        v->x = -a; v->y = m->x*m->y; v->z = m->x*m->z;
      }
    }  
    else
    {   
      u->x = 0.0; u->y = -m->z; u->z = m->y; a = 1.0 - m->x*m->x;
      v->x = a; v->y = -m->x*m->y; v->z = -m->x*m->z;
    }
  }  
  
  /* Normalize orthogonal vectors -- 32ns */
  a = 1/sqrt(a);
  u->x *= a;
  u->y *= a;
  u->z *= a;
  v->x *= a;
  v->y *= a;
  v->z *= a;
}


/*************************************************************************
exact_disk:
  In: location of moving molecule at time of collision
      movement vector for moving molecule
      interaction radius
      subvolume the moving molecule is in
      the moving molecule
      the target molecule at time of collision
  Out: The fraction of a full interaction disk that is actually
       accessible to the moving molecule, computed exactly from the
       geometry, or TARGET_OCCLUDED if the path to the target molecule is
       blocked.  If there is a memory error, it returns EXD_OUT_OF_MEMORY.
*************************************************************************/

double exact_disk(struct vector3 *loc,struct vector3 *mv,double R,struct subvolume *sv,struct molecule *moving,struct molecule *target)
{
#define EXD_SPAN_CALC(v1,v2,p) ((v1)->u - (p)->u)*((v2)->v - (p)->v)  -  ((v2)->u - (p)->u)*((v1)->v - (p)->v)
#define EXD_TIME_CALC(v1,v2,p) ((p)->u*(v1)->v - (p)->v*(v1)->u) / ((p)->v*((v2)->u-(v1)->u) - (p)->u*((v2)->v-(v1)->v))
  struct wall_list *wl;
  struct wall *w;
  struct vector3 llf,urb;
  
  struct exd_vector3 v0muv,v1muv,v2muv;
  struct exd_vertex pa,pb;
  struct exd_vertex *ppa,*ppb,*pqa,*pqb,*vertex_head,*vp,*vq,*vr,*vs;
  struct rxn *rx;
  double pa_pb;
  int n_verts,n_edges;
  int p_flags;
  
  double R2;
  int uncoordinated;
  struct vector3 m,u,v;
  struct exd_vector3 Lmuv;
  struct exd_vertex g;
  double m2_i;
  double l_n,m_n;
  double a,b,c,d,r,s,t,A,zeta,last_zeta;
  int i;
  
  /* Initialize */
  vertex_head = NULL;
  n_verts = 0;
  n_edges = 0;
  
  /* Partially set up coordinate systems for first pass */
  R2 = R*R;
  m2_i = 1.0/(mv->x*mv->x + mv->y*mv->y + mv->z*mv->z);
  uncoordinated = 1;
  Lmuv.m = Lmuv.u = Lmuv.v = 0.0;  /* Keep compiler happy */
  g.u = g.v = g.r2 = g.zeta = 0.0; /* More compiler happiness */

  /* Find walls that occlude the interaction disk (or block the reaction) */
  for (wl=sv->wall_head;wl!=NULL;wl=wl->next)
  {
    w = wl->this_wall;
    
    /* Ignore this wall if it is too far away! */
    
    /* Find distance from plane of wall to molecule */    
    l_n = loc->x*w->normal.x + loc->y*w->normal.y + loc->z*w->normal.z;
    d = w->d - l_n;
    
    /* See if we're within interaction distance of wall */
    m_n = mv->x*w->normal.x + mv->y*w->normal.y + mv->z*w->normal.z;
    
    if ( d*d >= R2*(1 - m2_i*m_n*m_n) ) continue;

    
    /* Ignore this wall if no overlap between wall & disk bounding boxes */
    
    /* Find wall bounding box */
    urb.x = llf.x = w->vert[0]->x;
    if (w->vert[1]->x < llf.x) llf.x = w->vert[1]->x;
    else urb.x = w->vert[1]->x;
    if (w->vert[2]->x < llf.x) llf.x = w->vert[2]->x;
    else if (w->vert[2]->x > urb.x) urb.x = w->vert[2]->x;

    urb.y = llf.y = w->vert[0]->y;
    if (w->vert[1]->y < llf.y) llf.y = w->vert[1]->y;
    else urb.y = w->vert[1]->y;
    if (w->vert[2]->y < llf.y) llf.y = w->vert[2]->y;
    else if (w->vert[2]->y > urb.y) urb.y = w->vert[2]->y;

    urb.z = llf.z = w->vert[0]->z;
    if (w->vert[1]->z < llf.z) llf.z = w->vert[1]->z;
    else urb.z = w->vert[1]->z;
    if (w->vert[2]->z < llf.z) llf.z = w->vert[2]->z;
    else if (w->vert[2]->z > urb.z) urb.z = w->vert[2]->z;
    
    /* Reject those without overlapping bounding boxes */
    b = R2*(1.0-mv->x*mv->x*m2_i);
    a = llf.x - loc->x; if (a>0 && a*a >= b) continue;
    a = loc->x - urb.x; if (a>0 && a*a >= b) continue;

    b = R2*(1.0-mv->y*mv->y*m2_i);
    a = llf.y - loc->y; if (a>0 && a*a >= b) continue;
    a = loc->y - urb.y; if (a>0 && a*a >= b) continue;

    b = R2*(1.0-mv->z*mv->z*m2_i);
    a = llf.z - loc->z; if (a>0 && a*a >= b) continue;
    a = loc->z - urb.z; if (a>0 && a*a >= b) continue;

    
    /* Ignore this wall if moving molecule can travel through it */
    
    /* Reject those that the moving particle can travel through */
    if ( (moving->properties->flags && CAN_MOLWALL) != 0 )
    {
      rx = trigger_intersect(moving->properties->hashval,(struct abstract_molecule*)moving,0,w);
      if (rx != NULL && (rx->n_pathways==RX_TRANSP))
      {
	continue;
      }
    }
    
    
    /* Find line of intersection between wall and disk */
    
    /* Set up coordinate system and convert vertices */
    if (uncoordinated)
    {
      exd_coordize(mv,&m,&u,&v);
      
      Lmuv.m = loc->x*m.x + loc->y*m.y + loc->z*m.z;    
      Lmuv.u = loc->x*u.x + loc->y*u.y + loc->z*u.z;    
      Lmuv.v = loc->x*v.x + loc->y*v.y + loc->z*v.z;
      
      g.u = (target->pos.x - loc->x)*u.x + (target->pos.y - loc->y)*u.y + (target->pos.z - loc->z)*u.z;
      g.v = (target->pos.x - loc->x)*v.x + (target->pos.y - loc->y)*v.y + (target->pos.z - loc->z)*v.z;
      g.r2 = g.u*g.u+g.v*g.v;
      g.zeta = exd_zetize(g.v,g.u);
      
      uncoordinated=0;
    }
    
    v0muv.m = w->vert[0]->x*m.x + w->vert[0]->y*m.y + w->vert[0]->z*m.z - Lmuv.m;
    v0muv.u = w->vert[0]->x*u.x + w->vert[0]->y*u.y + w->vert[0]->z*u.z - Lmuv.u;
    v0muv.v = w->vert[0]->x*v.x + w->vert[0]->y*v.y + w->vert[0]->z*v.z - Lmuv.v;
    
    v1muv.m = w->vert[1]->x*m.x + w->vert[1]->y*m.y + w->vert[1]->z*m.z - Lmuv.m;
    v1muv.u = w->vert[1]->x*u.x + w->vert[1]->y*u.y + w->vert[1]->z*u.z - Lmuv.u;
    v1muv.v = w->vert[1]->x*v.x + w->vert[1]->y*v.y + w->vert[1]->z*v.z - Lmuv.v;

    v2muv.m = w->vert[2]->x*m.x + w->vert[2]->y*m.y + w->vert[2]->z*m.z - Lmuv.m;
    v2muv.u = w->vert[2]->x*u.x + w->vert[2]->y*u.y + w->vert[2]->z*u.z - Lmuv.u;
    v2muv.v = w->vert[2]->x*v.x + w->vert[2]->y*v.y + w->vert[2]->z*v.z - Lmuv.v;
    
    /* Draw lines between points and pick intersections with plane of m=0 */
    if ( (v0muv.m < 0) == (v1muv.m < 0) )  /* v0,v1 on same side */
    {
      if ( (v2muv.m < 0) == (v1muv.m < 0) ) continue;
      t = v0muv.m/(v0muv.m-v2muv.m);
      pa.u = v0muv.u + (v2muv.u-v0muv.u)*t;
      pa.v = v0muv.v + (v2muv.v-v0muv.v)*t;
      t = v1muv.m/(v1muv.m-v2muv.m);
      pb.u = v1muv.u + (v2muv.u-v1muv.u)*t;
      pb.v = v1muv.v + (v2muv.v-v1muv.v)*t;
    }
    else if ( (v0muv.m<0) == (v2muv.m<0) ) /* v0,v2 on same side */
    {
      t = v0muv.m/(v0muv.m-v1muv.m);
      pa.u = v0muv.u + (v1muv.u-v0muv.u)*t;
      pa.v = v0muv.v + (v1muv.v-v0muv.v)*t;
      t = v2muv.m/(v2muv.m-v1muv.m);
      pb.u = v2muv.u + (v1muv.u-v2muv.u)*t;
      pb.v = v2muv.v + (v1muv.v-v2muv.v)*t;
    }
    else /* v1, v2 on same side */
    {
      t = v1muv.m/(v1muv.m-v0muv.m);
      pa.u = v1muv.u + (v0muv.u-v1muv.u)*t;
      pa.v = v1muv.v + (v0muv.v-v1muv.v)*t;
      t = v2muv.m/(v2muv.m-v0muv.m);
      pb.u = v2muv.u + (v0muv.u-v2muv.u)*t;
      pb.v = v2muv.v + (v0muv.v-v2muv.v)*t;
    }
    
    /* Intersect line with circle; skip this wall if no intersection */
    pa.r2 = pa.u*pa.u + pa.v*pa.v;
    pb.r2 = pb.u*pb.u + pb.v*pb.v;
    t=0; s=1;
    
    if (pa.r2 > R2 || pb.r2 > R2)
    {
      pa_pb = pa.u*pb.u + pa.v*pb.v;
      a = 1.0/(pa.r2 + pb.r2 - 2*pa_pb);
      b = (pa_pb - pa.r2)*a;
      c = (R2 - pa.r2)*a;
      d = b*b+c;
      if (d<=0) continue;
      d = sqrt(d);
      t = -b-d;
      if (t>=1) continue;
      if (t<0) t=0;
      s = -b+d;
      if (s<=0) continue;
      if (s>1) s=1;
    }
    
    
    /* Add this edge to the growing list, or return -1 if edge blocks target */
    
    /* Construct final endpoints and prepare to store them */
    ppa = (struct exd_vertex*)mem_get( sv->local_storage->exdv );
    ppb = (struct exd_vertex*)mem_get( sv->local_storage->exdv );
    if (ppa==NULL || ppb==NULL) return EXD_OUT_OF_MEMORY;
    if (t>0)
    {
      ppa->u = pa.u + t*(pb.u-pa.u);
      ppa->v = pa.v + t*(pb.v-pa.v);
      ppa->r2 = ppa->u*ppa->u + ppa->v*ppa->v;
      ppa->zeta = exd_zetize(ppa->v,ppa->u);
    }
    else
    {
      ppa->u = pa.u; ppa->v = pa.v; ppa->r2 = pa.r2;
      ppa->zeta = exd_zetize(pa.v,pa.u);
    }
    if (s<1)
    {
      ppb->u = pa.u + s*(pb.u-pa.u);
      ppb->v = pa.v + s*(pb.v-pa.v);
      ppb->r2 = ppb->u*ppb->u + ppb->v*ppb->v;
      ppb->zeta = exd_zetize(ppb->v,ppb->u);
    }
    else
    {
      ppb->u = pb.u; ppb->v = pb.v; ppb->r2 = pb.r2;
      ppb->zeta = exd_zetize(pb.v,pb.u);
    }
    
    /* It's convenient if ppa is earlier, ccw, than ppb */
    a = (ppb->zeta - ppa->zeta);
    if (a<0) a+=4.0;
    if (a >= 2.0)
    {
      vp = ppb; ppb=ppa; ppa=vp;
      a=4.0-a;
    }
    
    /* Detect a blocked reaction: line is between origin and target */
    b = (g.zeta - ppa->zeta);
    if (b<0) b+=4.0;
    
    if (b<a)
    {
      c = (ppa->u-g.u)*(ppb->v-g.v)-(ppa->v-g.v)*(ppb->u-g.u);
      if (c<=0) /* Blocked!! */
      {
	ppa->next=ppb; ppb->next=vertex_head;
	mem_put_list( sv->local_storage->exdv , ppa );
	return TARGET_OCCLUDED;
      }
    }
    
    ppa->role = EXD_HEAD;
    ppb->role = EXD_TAIL;
    ppa->e = ppb;
    ppb->e = NULL;
    
    ppb->next = vertex_head;
    ppa->next = ppb;
    vertex_head = ppa;
    n_verts += 2;
    n_edges++;
  }

  
  /* Find partition boundaries that occlude the interaction disk */
  if (!world->use_expanded_list) /* We'll hit partitions */
  {  
    /* First see if any overlap */
    p_flags = 0;
    
    d = loc->x - world->x_fineparts[ sv->llf.x ];
    if (d<R)
    {
      c = R2*(mv->y*mv->y + mv->z*mv->z)*m2_i;
      if (d*d<c) p_flags |= X_NEG_BIT;
      d = world->x_fineparts[ sv->urb.x ] - loc->x;
      if (d*d<c) p_flags |= X_POS_BIT;
    }
    else
    {
      d = world->x_fineparts[ sv->urb.x ] - loc->x;
      if (d<R && d*d<R2*(mv->y*mv->y + mv->z*mv->z)*m2_i) p_flags |= X_POS_BIT;
    }
  
    d = loc->y - world->y_fineparts[ sv->llf.y ];
    if (d<R)
    {
      c = R2*(mv->x*mv->x + mv->z*mv->z)*m2_i;
      if (d*d<c) p_flags |= Y_NEG_BIT;
      d = world->y_fineparts[ sv->urb.y ] - loc->y;
      if (d*d<c) p_flags |= Y_POS_BIT;
    }
    else
    {
      d = world->y_fineparts[ sv->urb.y ] - loc->y;
      if (d<R && d*d<R2*(mv->x*mv->x + mv->z*mv->z)*m2_i) p_flags |= Y_POS_BIT;
    }
  
    d = loc->z - world->z_fineparts[ sv->llf.z ];
    if (d<R)
    {
      c = R2*(mv->y*mv->y + mv->x*mv->x)*m2_i;
      if (d*d<c) p_flags |= Z_NEG_BIT;
      d = world->z_fineparts[ sv->urb.z ] - loc->z;
      if (d*d<c) p_flags |= Z_POS_BIT;
    }
    else
    {
      d = world->z_fineparts[ sv->urb.z ] - loc->z;
      if (d<R && d*d<R2*(mv->y*mv->y + mv->x*mv->x)*m2_i) p_flags |= Z_POS_BIT;
    }
  
    /* Now find the lines created by any that do overlap */
    if (p_flags)
    {
      if (uncoordinated) exd_coordize(mv,&m,&u,&v);
      uncoordinated = 0;
      
      for (i=1;i<=p_flags;i*=2)
      {
	if ((i & p_flags)!=0)
	{
	  /* Load up the relevant variables */
	  switch (i)
	  {
	    case X_NEG_BIT:
	      d = world->x_fineparts[ sv->llf.x ] - loc->x;
	      a = u.x; b = v.x;
	      break;
	    case X_POS_BIT:
	      d = world->x_fineparts[ sv->urb.x ] - loc->x;
	      a = u.x; b = v.x;
	      break;
	    case Y_NEG_BIT:
	      d = world->y_fineparts[ sv->llf.y ] - loc->y;
	      a = u.y; b = v.y;
	      break;
	    case Y_POS_BIT:
	      d = world->y_fineparts[ sv->urb.y ] - loc->y;
	      a = u.y; b = v.y;
	      break;
	    case Z_NEG_BIT:
	      d = world->z_fineparts[ sv->llf.z ] - loc->z;
	      a = u.z; b = v.z;
	      break;
	    case Z_POS_BIT:
	      d = world->z_fineparts[ sv->urb.z ] - loc->z;
	      a = u.z; b = v.z;
	      break;
	    default:
	      continue;
	  }
	  
	  if (a==0)
	  {
	    s = d/b;
	    if (s*s>R2)
	    {
	      fprintf(world->log_file, "File '%s', Line %ld: MCell should not come to this point.  Please report this message.\n", __FILE__, (long)__LINE__); 
	      continue;
	    }
	    t = sqrt(R2-s*s);
	    pa.u = t; pa.v = s;
	    pb.u = -t; pb.v = s;
	  }
	  else if (b==0)
	  {
	    t = d/b;
	    if (t*t>R2)
	    {
	      fprintf(world->log_file, "File '%s', Line %ld: MCell should not come to this point.  Please report this message.\n", __FILE__, (long)__LINE__); 
	      continue;
	    }
	    s = sqrt(R2-t*t);
	    pa.u = t; pa.v = s;
	    pb.u = t; pb.v = -s;
	  }
	  else
	  {
	    c = a*a+b*b;
	    s = d*b;
	    if (d*d>R2*c)
	    {
	      fprintf(world->log_file, "File '%s', Line %ld: MCell should not come to this point.  Please report this message.\n", __FILE__, (long)__LINE__); 
	      continue;
	    }
	    t = sqrt(R2*c-d*d);
	    c = 1.0/c;
	    r = 1.0/a;
	    pa.v = c*(s+t*a);
	    pa.u = (d-b*pa.v)*r;
	    pb.v = c*(s-t*a);
	    pb.u = (d-b*pb.v)*r;
	  }
	  
	  /* Create memory for the pair of vertices */
	  ppa = (struct exd_vertex*)mem_get( sv->local_storage->exdv );
	  ppb = (struct exd_vertex*)mem_get( sv->local_storage->exdv );
	  if (ppa==NULL || ppb==NULL) return EXD_OUT_OF_MEMORY;
	  
	  a = exd_zetize(pa.v,pa.u);
	  b = exd_zetize(pb.v,pb.u);
	  c = b-a;
	  if (c<0) c += 4;
	  if (c<2)
	  {
	    ppa->u = pa.u; ppa->v = pa.v; ppa->r2 = pa.u*pa.u+pa.v*pa.v; ppa->zeta = a;
	    ppb->u = pb.u; ppb->v = pb.v; ppb->r2 = pb.u*pb.u+pb.v*pb.v; ppb->zeta = b;
	  }
	  else
	  {
	    ppb->u = pa.u; ppb->v = pa.v; ppb->r2 = pa.u*pa.u+pa.v*pa.v; ppb->zeta = a;
	    ppa->u = pb.u; ppa->v = pb.v; ppa->r2 = pb.u*pb.u+pb.v*pb.v; ppa->zeta = b;	  
	  }
	  
	  ppa->role = EXD_HEAD;
	  ppb->role = EXD_TAIL;
	  ppa->e = ppb;
	  ppb->e = NULL;
	  
	  ppb->next = vertex_head;
	  ppa->next = ppb;
	  vertex_head = ppa;
	  n_verts += 2;
	  n_edges++;	
	}
      }
    }
  }
  
  
  /* Now that we have everything, see if we can perform simple calculations */
  
  /* Did we even find anything?  If not, return full area */
  if (n_edges==0)
  {
    return 1.0;
  }
  /* If there is only one edge, just calculate it */
  else if (n_edges==1)
  {
    ppa = vertex_head;
    ppb = ppa->e;
    
    a = ppa->u*ppb->u+ppa->v*ppb->v;
    b = ppa->u*ppb->v-ppa->v*ppb->u;
    if (a<=0) /* Angle > pi/2 */
    {
      s = atan(-a/b)+0.5*MY_PI;
    }
    else
    {
      s = atan(b/a);
    }
    A = (0.5*b + R2*(MY_PI-0.5*s))/(MY_PI*R2);
    
    mem_put_list( sv->local_storage->exdv , vertex_head );
    return A;
  }
  
  
  /* If there are multiple edges, calculating area is more complex. */
  
  /* Insertion sort the multiple edges */
  vp=vertex_head->next;
  ppa = ppb = vertex_head;
  ppa->next = NULL;
  ppa->span = NULL;
  while (vp!=NULL)
  {
    /* Snip off one item from old list to add */
    vp->span = NULL;
    vq = vp->next;
    
    /* Add it to list with ppa as head and ppb as tail */
    if (vp->zeta < ppa->zeta)
    {
      vp->next = ppa;
      ppa = vp;
    }
    else
    {
      for (pqa=ppa;pqa->next!=NULL;pqa=pqa->next)
      {
	if (vp->zeta < pqa->next->zeta) break;
      }
      vp->next = pqa->next;
      pqa->next = vp;
      if (vp->next==NULL) ppb = vp;
    }
    
    /* Repeat for remainder of old list */
    vp = vq;
  }
  
  /* Close circular list */ 
  vertex_head = ppa;
  ppb->next = ppa;
  
  
  /* Walk around the circle, inserting points where lines cross */
  ppb=NULL;
  for (ppa=vertex_head ; ppa!=vertex_head || ppb==NULL; ppa=ppa->next)
  {
    if (ppa->role != EXD_HEAD) continue;
    ppb=ppa->e;
    
    for (pqa=ppa->next;pqa!=ppb;pqa=pqa->next)
    {
      if (pqa->role != EXD_HEAD) continue;
      pqb=pqa->e;
      
      /* Create displacement vectors */
      pa.u = ppb->u-ppa->u;
      pa.v = ppb->v-ppa->v;
      pb.u = pqb->u-pqa->u;
      pb.v = pqb->v-pqa->v;
      r = pb.u*pa.v - pa.u*pb.v;
      
      /* Check if lines are parallel--combine if so */
      if (r*r<EPS_C*(pa.u*pa.u+pa.v*pa.v)*(pb.u*pb.u+pb.v*pb.v))
      {
	pqa->e = NULL;
	pqa->role = EXD_OTHER;
	
	a = pqb->zeta - ppb->zeta;
	if (a<0) a += 4.0;
	
	if (a>2) /* Other line is completely contained inside us */
	{
	  pqb->role = EXD_OTHER;
	}
	else  /* We have a new endpoint, so we need to check crosses again */
	{
	  ppa->e = pqb;
	  ppb->role = EXD_OTHER;
	  ppb=pqb;
	  pqa=ppa;
	}
	continue;
      }
      
      /* Check if these lines cross and find times at which they do */
      s = (ppa->u-pqa->u)*pa.v - (ppa->v-pqa->v)*pa.u;
      if (s*r <= EPS_C*R2*R2) continue;
      t = s/r;
      if (t>=1-EPS_C) continue;
      if (pa.u*pa.u>pa.v*pa.v)
      {
	s = (pqa->u-ppa->u+t*pb.u)*pa.u;
	if (s <= EPS_C*R2 || s >= pa.u*pa.u*(1.0-EPS_C)) continue;
      }
      else
      {
	s = (pqa->v-ppa->v+t*pb.v)*pa.v;
	if (s <= EPS_C*R2 || s >= pa.v*pa.v*(1.0-EPS_C)) continue;
      }
      
      /* Create intersection point */
      vq = (struct exd_vertex*)mem_get( sv->local_storage->exdv );
      if (vq==NULL) return EXD_OUT_OF_MEMORY;
      vq->u = pqa->u + t*pb.u;
      vq->v = pqa->v + t*pb.v;
      vq->r2 = vq->u*vq->u + vq->v*vq->v;
      vq->zeta = exd_zetize(vq->v,vq->u);
      vq->e = ppb;
      vq->span = NULL;
      vq->role = EXD_CROSS;
      
      /* Insert new point into the list */
      for (vp=ppa;vp!=ppb;vp=vp->next)
      {
	a = vq->zeta-vp->next->zeta;
	if (a>2.0) a -= 4.0;
	else if (a<-2.0) a += 4.0;
	
	if (a<0) break;
      }
      
      vq->next = vp->next;
      vp->next = vq;
      if (vq->zeta < vertex_head->zeta) vertex_head = vq;
    }
  }
  
  /* Collapse nearby points in zeta and R */
  for (vp=vertex_head,vq=NULL ; vq!=vertex_head ; vp=vq)
  {
    for (vq=vp->next;vq!=vertex_head;vq=vq->next)
    {
      if (vq->zeta-vp->zeta < EPS_C)
      {
	vq->zeta = vp->zeta;
	if (-EPS_C < vq->r2-vp->r2 && EPS_C > vq->r2-vp->r2)
	{
	  vq->r2=vp->r2;
	  /* Mark crosses that occur multiple times--only need one */
//	  if (vq->role==EXD_CROSS && vp->role != EXD_OTHER) vq->role = EXD_OTHER;
//	  else if (vp->role==EXD_CROSS && vq->role != EXD_OTHER) vp->role = EXD_OTHER;
	}
      }
      else break;
    }
  }
  
  /* Register all spanning line segments */
  vq=NULL;
  for (vp=vertex_head ; vp!=vertex_head || vq==NULL ; vp=vp->next)
  {
    if (vp->role != EXD_HEAD) continue;
    
    for (vq=vp->next ; vq!=vp->e ; vq=vq->next)
    {
      if (vq->zeta==vp->zeta) continue;
      if (vq->zeta==vp->e->zeta) break;
      if (vq->role==EXD_OTHER) continue;
      
      vr = (struct exd_vertex*) mem_get( sv->local_storage->exdv );
      if (vr==NULL) return EXD_OUT_OF_MEMORY;
      vr->next = vq->span;
      vq->span = vr;
      vr->e = vp;
      vr->zeta = vq->zeta;
      vr->role = EXD_SPAN;
    }
  }

  /* Now we finally walk through and calculate the area */  
  A=0.0;
  zeta=0.0;
  last_zeta=-1;
  vr = vs = NULL;
  for (vp=vertex_head ; zeta < 4.0-EPS_C ; vp=vp->next)
  {
    if (vp->role == EXD_OTHER) continue;
    if (vp->zeta==last_zeta) continue;
    last_zeta = vp->zeta;
    
    /* Store data for the next tentatively approved point */
    if (vs==&pa) vr=&pb;
    else vr=&pa;
    vr->u = vp->u;
    vr->v = vp->v;
    vr->zeta = vp->zeta;
    if (vp->role == EXD_TAIL)
    {
      vr->r2 = R2*(1.0+EPS_C);
      vr->e = NULL;
    }
    else
    {
      vr->r2 = vp->r2;
      vr->e = vp->e;
    }
    
    /* Check head points at same place to see if they're closer */
    for (vq=vp->next ; vq->zeta==last_zeta ; vq=vq->next)
    {
      if (vq->role==EXD_HEAD)
      {
	if (vq->r2 < vp->r2 || vr->e==NULL)
	{
	  vr->u = vq->u;
	  vr->v = vq->v;
	  vr->r2 = vq->r2;
	  vr->e = vq->e;
	}
	else if (vq->r2 == vr->r2)
	{
	  b = EXD_SPAN_CALC(vr,vr->e,vq->e);
	  if (b>0) vr->e = vq->e;
	}
      }
    }
	
    
    /* Check each span to see if anything is closer than our approval point */
    for (vq=vp->span ; vq!=NULL ; vq=vq->next)
    {
      ppa = vq->e;
      ppb = ppa->e;
      b = EXD_SPAN_CALC(ppa,ppb,vr);
      c = b*b;
      if (c < R2*R2*EPS_C)  /* Span crosses the point */
      {
	if (vr->e==NULL)
	{
	  vr->r2 = vr->u*vr->u + vr->v*vr->v;
	  vr->e = ppb;
	}
	else
	{
	  b = EXD_SPAN_CALC(vr,vr->e,ppb);
	  if (b>0) vr->e = ppb;
	}
      }
      else if (b<0 || vr->e==NULL)  /* Span is inside the point or spans tail */
      {
	t = EXD_TIME_CALC(ppa,ppb,vp);
	vr->u = ppa->u + t*(ppb->u - ppa->u);
	vr->v = ppa->v + t*(ppb->v - ppa->v);
	vr->r2 = vr->u*vr->u + vr->v*vr->v;
	vr->e = ppb;
      }
    }
    
    /* Should have an approved point in vr */
    if (vs==NULL) /* No angle traversed yet */
    {
      vs = vr;
    }
    else
    {
      c = vr->zeta - vs->zeta;
      if (c<0) c+=4.0;
      if (/*vr->e != vs->e || c+zeta >= 4.0-EPS_C*/ c>EPS_C)
      {
	zeta += c;
	if (vs->e == NULL || (vs->e->zeta-vs->zeta)*(vs->e->zeta-vs->zeta) < EPS_C*EPS_C)
	{
	  if (c>=2.0) /* More than pi */
	  {
	    vs->u=-vs->u;
	    vs->v=-vs->v;
	    A += 0.5*MY_PI*R2;
	  }
	  a = vs->u*vr->u + vs->v*vr->v;
	  b = vs->u*vr->v - vs->v*vr->u;
	  if (a<=0) /* More than pi/2 */
	  {
    	    s = atan( -a / b )+0.5*MY_PI;
	  }
	  else
	  {
	    s = atan( b / a );
	  }
	  A += 0.5*s*R2;
	}
	else
	{
	  if (vs->e->zeta == vr->zeta)
	  {
	    A += 0.5*( vs->u*vs->e->v - vs->v*vs->e->u );
	  }
	  else
	  {
	    t = EXD_TIME_CALC(vs,vs->e,vr);
	    b = vs->u+(vs->e->u-vs->u)*t;
	    c = vs->v+(vs->e->v-vs->v)*t;
	    A += 0.5*(vs->u*c - vs->v*b);
	  }
	}
	vs=vr;
      }
      else
      { 
	if (vr->e!=NULL) vs=vr; 
      }
    }
  }
  
  
  /* Finally, let's clean up the mess we made! */
  
  /* Deallocate lists */
  ppa = vertex_head->next;
  vertex_head->next = NULL;
  mem_put_list( sv->local_storage->exdv , ppa );
  
  /* Return fractional area */
  
  return A/(MY_PI*R2);
  
#undef EXD_TIME_CALC
#undef EXD_SPAN_CALC 
}

/**************************/
/** done with exact_disk **/
/**************************/


#ifdef DEBUG
/* Debugging function searching for misplaced molecules in the Min simulation */
void scream_if_misplaced(struct molecule *m)
{
  if (m->pos.x*world->length_unit > 4.0 || m->pos.x < 0.0)
  {
    printf("Out of X bounds.\n");
  }
  if (fabs(m->pos.y*world->length_unit) > 0.5)
  {
    printf("Out of Y bounds.\n");
  }
  if (fabs(m->pos.z*world->length_unit) > 0.5)
  {
    printf("Out of Z bounds.\n");
  }  
}

/* Debugging function: print a string and some details about a molecule. */
void tell_loc(struct molecule *m,char *s)
{
  if (0 || s[0] == '\0')
  printf("%sMy name is %x and I live at %.3f,%.3f,%.3f\n",
         s,(int)m,m->pos.x*world->length_unit,m->pos.y*world->length_unit,m->pos.z*world->length_unit);
}

/* Debugging function: search a schedule for a specific item */
int search_schedule_for_me(struct schedule_helper *sch,struct abstract_element *ae)
{
  struct abstract_element *aep;
  int i;
  
  for ( aep = sch->current ; aep != NULL ; aep = aep->next )
  {
    if (aep == ae) return 1;
  }
  for (i=0;i<sch->buf_len;i++)
  {
    for ( aep = sch->circ_buf_head[i] ; aep!=NULL ; aep = aep->next )
    {
      if (aep == ae) return 1;
    }
  }
  
  if (sch->next_scale == NULL) return 0;
  else return search_schedule_for_me(sch->next_scale , ae);
}

/* Debugging function: see if an object was allocated by a given mem_helper */
int search_memory_for_me(struct mem_helper *mh,struct abstract_list *al)
{
  int i;
  
  for (i=0;i<mh->buf_len;i++)
  {
    if (al == (struct abstract_list*)(mh->heap_array + mh->record_size*i)) return 1;
  }
  
  if (mh->next_helper == NULL) return 0;
  else return search_memory_for_me(mh->next_helper,al);
}

/* Debugging function: see if we got a circular molecule list inside a SV */
int test_subvol_for_circular(struct subvolume *sv)
{
  struct molecule *mp,*smp,*psmp;
  int warned = 0;
  
  psmp = NULL;
  mp = smp = sv->mol_head;
  do
  {
    if (!warned && smp->subvol != sv)
    {
      printf("Occupancy leak from %x to %x through %x to %x\n",(int)sv,(int)smp->subvol,(int)psmp,(int)smp);
      warned = 1;
    }
    psmp = smp;
    smp = smp->next_v;
    mp = mp->next_v;
    if (mp!=NULL) mp = mp->next_v;
  } while (mp != NULL && smp != NULL && mp != smp);
  
  if (mp != NULL) return 1;
  return 0;
}
#endif



/*************************************************************************
gather_walls_first:
  In: A list of collisions
      Tolerance below which two collisions are considered simultaneous
  Out: The same list of collisions with simultaneous ones sorted to put
       walls first among equals.
  Note: This isn't very efficient for lots of coincident objects since
        it uses a bubble-sort-like algorithm.  The list that is passed
	in should have already been mergesorted, though, so it shouldn't
	be too bad.
*************************************************************************/

struct collision* gather_walls_first(struct collision *shead,double tol)
{
  struct collision *cp,*co,*ci,*cf,*ct;
  
  co = NULL;
  cp = shead;
  while (cp->next != NULL)
  {
    if (cp->next->t - cp->t > tol || (cp->what&COLLIDE_WALL) != 0 )
    {
      co = cp;
      cp = cp->next;
    }
    else
    {
      for (ct=cp; ; ct=ct->next)  /* Find any wall */
      {
        if (ct->next==NULL) return shead;
        
        if (ct->next->t - cp->t > tol)
        {
          co = ct;
          cp = ct->next;
          break;
        }
        if ((ct->next->what&COLLIDE_WALL) != 0)
        {
          ci = ct->next;
          for (cf=ci ; cf->next!=NULL ; cf=cf->next)  /* Find last wall */
          {
            if (cf->next->t - cp->t > tol || (cf->next->what&COLLIDE_WALL)==0) break;
          }
          
          if (co==NULL) shead = ci;
	  else co->next=ci;
	  
          ct->next = cf->next;
          cf->next = cp;
          co = cf;
          
          break;
        }
      }
    }
  }
  return shead;
}


#if 0
/* Under development. */
int can_hit_target(struct molecule *m,struct molecule *targ,struct subvolume *sv)
{
  const double TOL = 10.0*EPS_C;   /* Two walls are coincident if this close */
  struct vector3 to_target;                /* Vector from molecule to target */
  struct collision *list;         /* List of collisions between mol & target */
  struct *cp;                         /* Primary iterator for collision list */
  struct *x;                        /* Iterator for coincident walls in list */
  struct rxn *rx;
  int k;
  
  to_target.x = targ->pos.x - m->pos.x;
  to_target.y = targ->pos.y - m->pos.y;
  to_target.z = targ->pos.z - m->pos.z;
  
  list = ray_trace(m,NULL,sv,&to_target);
  if (list == NULL) return RX_NO_MEM;

  if (list->next == NULL) return RX_A_OK;
  
  list = (struct collision*)ae_list_sort((struct abstract_element*)list);
  
  for (cp = list ; cp != NULL ; cp = cp->next)
  {
    if ((cp->what & COLLIDE_WALL) == 0) return RX_A_OK;
    
    if (cp->next != NULL && (cp->next->what&COLLIDE_WALL)!=0 && (cp->next->t - cp->t) < TOL)
    {
      if ( (m->properties->flags&CAN_MOLWALL) == 0 ) rx = NULL;
      else
      {
	if ((cp->what & COLLIDE_MASK) == COLLIDE_FRONT) k = 1;
	else k = -1;
      
	rx = trigger_intersect(
                       m->properties->hashval,(struct abstract_molecule*)m,k,
                       (struct *wall)cp->target
                     );
      }
    }
    else
    {
      if ( (m->properties->flags&CAN_MOLWALL) == 0 ) rx = NULL;
      else
      {
	if ((cp->what & COLLIDE_MASK) == COLLIDE_FRONT) k = 1;
	else k = -1;
      
	rx = trigger_intersect(
                       m->properties->hashval,(struct abstract_molecule*)m,k,
                       (struct *wall)cp->target
                     );
      }
      
      if (rx==NULL) return RX_BLOCKED;
      if (rx->n_pathways <= RX_SPECIAL) continue;
      else return RX_BLOCKED;
    }
  }
}
#endif


/****************************************************************************
safe_diffusion_step:
  In: molecule that is moving
      linked list of potential collisions with molecules from the 
                starting subvolume
  Out: The estimated number of diffusion steps this molecule can take before
       something interesting might happen to it, or 1.0 if something might
       happen within one timestep.
  Note: Each molecule uses its own timestep.  Only molecules that the moving
  	molecule can react with directly are counted (secondary reaction
	products are ignored, so might be skipped).  "Might happen" is to
	the 99% confidence level (i.e. the distance you'd have to go before
	1% of the molecules will have gotten far enough to have a chance of
	reacting, although those 1% will probably not go in the right
	direction).  This doesn't take into account the diffusion of other
	target molecules, so it may introduce errors for clouds of molecules
	diffusing into each other from a distance.
	*FIXME*: Add a flag to make this be very conservative or to turn
	this off entirely, aside from the TIME_STEP_MAX= directive.
****************************************************************************/
double safe_diffusion_step(struct molecule *m,struct collision *shead)
{
  double d2;
  double d2_nearmax;
  double d2min = GIGANTIC;
  struct subvolume *sv = m->subvol;
  struct wall *w;
  struct wall_list *wl;
  struct collision *smash;
  double steps;
  struct molecule *mp;
  
  d2_nearmax = m->properties->space_step * world->r_step[ (int)(world->radial_subdivisions * MULTISTEP_PERCENTILE) ];
  d2_nearmax *= d2_nearmax;

  if ( (m->properties->flags & (CAN_MOLMOL|CANT_INITIATE)) == CAN_MOLMOL )
  {
    for (smash = shead ; smash != NULL ; smash = smash->next)
    {
      mp = (struct molecule*)smash->target;
      d2 = (m->pos.x - mp->pos.x)*(m->pos.x - mp->pos.x) +
	   (m->pos.y - mp->pos.y)*(m->pos.y - mp->pos.y) +
	   (m->pos.z - mp->pos.z)*(m->pos.z - mp->pos.z);
      if (d2 < d2min) d2min = d2;
    }
  }
  for (wl = sv->wall_head ; wl!=NULL ; wl = wl->next)
  {
    w = wl->this_wall;
    d2 = (w->normal.x*m->pos.x + w->normal.y*m->pos.y + w->normal.z*m->pos.z) - w->d;
    d2 *= d2;
    if (d2 < d2min) d2min = d2;
  }

  d2 = (m->pos.x - world->x_fineparts[ sv->llf.x ]);
  d2 *= d2;
  if (d2 < d2min) d2min = d2;
  
  d2 = (m->pos.x - world->x_fineparts[ sv->urb.x ]);
  d2 *= d2;
  if (d2 < d2min) d2min = d2;

  d2 = (m->pos.y - world->y_fineparts[ sv->llf.y ]);
  d2 *= d2;
  if (d2 < d2min) d2min = d2;
  
  d2 = (m->pos.y - world->y_fineparts[ sv->urb.y ]);
  d2 *= d2;
  if (d2 < d2min) d2min = d2;

  d2 = (m->pos.z - world->z_fineparts[ sv->llf.z ]);
  d2 *= d2;
  if (d2 < d2min) d2min = d2;
  
  d2 = (m->pos.z - world->z_fineparts[ sv->urb.z ]);
  d2 *= d2;
  if (d2 < d2min) d2min = d2;
  
  if (d2min < d2_nearmax) steps = 1.0;
  else steps = d2min / d2_nearmax;
	
  if (steps < MULTISTEP_WORTHWHILE) steps = 1.0;
  
  return steps;
}
 
/****************************************************************************
expand_collision_list:
  In: molecule that is moving
      displacement to the new location
      subvolume that we start in
      linked list of potential collisions with molecules from the 
                starting subvolume
  Out: To the current collision list of molecules we may react with in our 
       subvolume the molecules from the neighbor subvolumes are added.
       The molecules are added only when the molecule displacement bounding box 
       intersects with the subvolume bounding box.
       Returns new list of collisions.
	
****************************************************************************/
struct collision* expand_collision_list(struct molecule *m, struct vector3 *mv, struct subvolume *sv, struct collision *shead1)
{
  struct collision *smash;
  struct molecule *mp;
  /* neighbors of the current subvolume */
  struct subvolume *new_sv;
  struct rxn *rx;
  /* lower left and upper_right corners of the molecule path
     bounding box expanded by R. */
  struct vector3 path_llf, path_urb;
  /* lower left and upper right corners of the subvolume */
  struct vector3 new_sv_llf, new_sv_urb;
  struct vector3 v0, v1;
  /* index of the neighbor subvolume */
  int neighbor_index;
  int i;
  double d, d1, d2;  /* distance */
  double R, R2;  /* molecule interaction radius */
  
 
  R = (world->rx_radius_3d); 
  R2 = R*R;

  /* find the molecule path bounding box. */
  path_bounding_box(&(m->pos), mv, &path_llf, &path_urb);
  
     new_sv = (struct subvolume *)(sv->neighbor[X_POS]);
     if(new_sv != NULL)
     {
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue; 
                /* if molecule is farther than R from the shared SV plane */ 
                d = mp->pos.x - world->x_fineparts[sv->urb.x];
                if(d > R) continue; 
                     
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: Out of memory,  trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
     }
     new_sv = (struct subvolume *)(sv->neighbor[X_NEG]);
     if(new_sv != NULL)
     {
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV plane */ 
                d = world->x_fineparts[sv->llf.x] - mp->pos.x;
                if(d > R) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: Out of memory, trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}								
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
     }
     new_sv = (struct subvolume *)(sv->neighbor[Z_POS]);
     if(new_sv != NULL)
     {
       new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
       new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
       new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
       new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
       new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
       new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
       if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
         /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV plane */ 
                d = mp->pos.z - world->z_fineparts[sv->urb.z];
                if(d > R) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: Out of memory, trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


       }
     }

     new_sv = (struct subvolume *)(sv->neighbor[Z_NEG]);
     if(new_sv != NULL)
     {
       new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
       new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
       new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
       new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
       new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
       new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
       if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV plane */ 
                d = world->z_fineparts[sv->llf.z] - mp->pos.z;
                if(d > R) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


       }
     }

     new_sv = (struct subvolume *)(sv->neighbor[Y_POS]);
     if(new_sv != NULL)
     {
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV plane */ 
                d = mp->pos.y - world->y_fineparts[sv->urb.y];
                if(d > R) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
     }

     new_sv = (struct subvolume *)(sv->neighbor[Y_NEG]);
     if(new_sv != NULL)
     {
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV plane */ 
                d = world->y_fineparts[sv->llf.y] - mp->pos.y;
                if(d > R) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
     }

   /* let's reach subvolumes located "edge-to-edge" from the top face
      of the current subvolume. */
   /* go (-X and +Z) */
     	neighbor_index = sv->index -(world->nz_parts-1)*(world->ny_parts-1) + 1;
        if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
        {
   	  new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
          
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = mp->pos.z - world->z_fineparts[sv->urb.z];
                d2 = world->x_fineparts[sv->llf.x] - mp->pos.x;
                if((d1 > R) && (d2 > R)) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
   	}
   /* go (+Y and +Z) */
      neighbor_index = sv->index  + (world->nz_parts-1) + 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = mp->pos.z - world->z_fineparts[sv->urb.z];
                d2 = mp->pos.y -world->y_fineparts[sv->urb.y];                
                if((d1 > R) && (d2 > R)) continue;
  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
      } 
   
   /* go (+X and +Z) */
      neighbor_index = sv->index  + (world->nz_parts-1)*(world->ny_parts-1) + 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = mp->pos.z - world->z_fineparts[sv->urb.z];
                d2 = mp->pos.x - world->x_fineparts[sv->urb.x]; 
                if((d1 > R) && (d2 > R)) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* go (-Y and +Z) */
      neighbor_index = sv->index  - (world->nz_parts-1) + 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue; 
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = mp->pos.z - world->z_fineparts[sv->urb.z];
                d2 = world->y_fineparts[sv->llf.y] - mp->pos.y;
                if((d1 > R) && (d2 > R)) continue;
                
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* let's reach subvolumes located "edge-to-edge" from the bottom face
      of the current subvolume. */
   /* go (-X and -Z) */
      neighbor_index = sv->index - (world->nz_parts-1)*(world->ny_parts-1) - 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue; 
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = world->z_fineparts[sv->llf.z] - mp->pos.z;
                d2 = world->x_fineparts[sv->llf.x] - mp->pos.x;
                if((d1 > R) && (d2 > R)) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* go (+Y and -Z) */
      neighbor_index = sv->index  + (world->nz_parts-1) - 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = world->z_fineparts[sv->llf.z] - mp->pos.z;
                d2 = mp->pos.y - world->y_fineparts[sv->urb.y];
                if((d1 > R) && (d2 > R)) continue; 
                
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* go (+X and -Z) */
      neighbor_index = sv->index  + (world->nz_parts-1)*(world->ny_parts-1) - 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = world->z_fineparts[sv->llf.z] - mp->pos.z;
                d2 = mp->pos.x - world->x_fineparts[sv->urb.x];
                if((d1 > R) && (d2 > R)) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* go (-Y and -Z) */
      neighbor_index = sv->index  - (world->nz_parts-1) - 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = world->z_fineparts[sv->llf.z] - mp->pos.z;
                d2 = world->y_fineparts[sv->llf.y] - mp->pos.y;
                if((d1 > R) && (d2 > R)) continue;
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: Out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* let's reach subvolumes located "edge-to-edge" from the vertical edges
      of the current subvolume. */
   /* go (-Y and -X) */
      neighbor_index = sv->index  - (world->nz_parts-1) - (world->nz_parts-1)*(world->ny_parts-1);
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue; 
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = world->y_fineparts[sv->llf.y] - mp->pos.y;
                d2 = world->x_fineparts[sv->llf.x] - mp->pos.x;
                if((d1 > R) && (d2 > R)) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* go (+Y and -X) */
      neighbor_index = sv->index  + (world->nz_parts-1) - (world->nz_parts-1)*(world->ny_parts-1);
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = mp->pos.y - world->y_fineparts[sv->urb.y];
                d2 = world->x_fineparts[sv->llf.x] - mp->pos.x;  
                if((d1 > R) && (d2 > R)) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* go (+Y and +X) */
      neighbor_index = sv->index  + (world->nz_parts-1) + (world->nz_parts-1)*(world->ny_parts-1);
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue; 
                /* if molecule is farther than R from the shared SV edge */
                d1 = mp->pos.y - world->y_fineparts[sv->urb.y];
                d2 = mp->pos.x - world->x_fineparts[sv->llf.x];  
                if((d1 > R) && (d2 > R)) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }
   /* go (-Y and +X) */
      neighbor_index = sv->index  - (world->nz_parts-1) + (world->nz_parts-1)*(world->ny_parts-1);
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
          new_sv = &(world->subvol[neighbor_index]);
          new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
          new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
          new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
          new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
          new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
          new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
          if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){
            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                /* if molecule is farther than R from the shared SV edge */ 
                d1 = world->y_fineparts[sv->llf.y] - mp->pos.y; 
                d2 = mp->pos.x - world->x_fineparts[sv->llf.x];  
                if((d1 > R) && (d2 > R)) continue; 
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


          }
   
      }

   /* let's reach subvolumes located "corner-to-corner" from the top face
      of the current subvolume. */
   /* go (-X and -Y and +Z) */
      neighbor_index = sv->index - (world->nz_parts-1)*(world->ny_parts-1) - (world->nz_parts-1) + 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
        new_sv = &(world->subvol[neighbor_index]);
        new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
        new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
        new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
        new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
        new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
        new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
        if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){

                /* find coordinates of the shared corner */
                v0.x = world->x_fineparts[sv->llf.x];
                v0.y = world->y_fineparts[sv->llf.y];
                v0.z = world->z_fineparts[sv->urb.z];

            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                 
                /* check for the distance from the shared corner */
                vectorize(&(mp->pos), &v0,&v1);
                d = vect_length(&v1);
                if (d > R) continue;
 
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


        }
   
      }
   /* go (-X and +Y and +Z) */
      neighbor_index = sv->index - (world->nz_parts-1)*(world->ny_parts-1) + (world->nz_parts-1) + 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
        new_sv = &(world->subvol[neighbor_index]);
        new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
        new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
        new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
        new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
        new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
        new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
        if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){

                /* find coordinates of the shared corner */
                v0.x = world->x_fineparts[sv->llf.x];
                v0.y = world->y_fineparts[sv->urb.y];
                v0.z = world->z_fineparts[sv->urb.z];

            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                
                /* check for the distance from the shared corner */
                vectorize(&(mp->pos), &v0,&v1);
                d = vect_length(&v1);
                if (d > R) continue;
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


        }
   
      }
  /* go (+X and +Y and +Z) */
      neighbor_index = sv->index + (world->nz_parts-1)*(world->ny_parts-1) + (world->nz_parts-1) + 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
        new_sv = &(world->subvol[neighbor_index]);
        new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
        new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
        new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
        new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
        new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
        new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
        if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){

                /* find coordinates of the shared corner */
                v0.x = world->x_fineparts[sv->urb.x];
                v0.y = world->y_fineparts[sv->urb.y];
                v0.z = world->z_fineparts[sv->urb.z];

            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                  
                /* check for the distance from the shared corner */
                vectorize(&(mp->pos), &v0,&v1);
                d = vect_length(&v1);
                if (d > R) continue;

      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


        }
   
      }
  /* go (+X and -Y and +Z) */
      neighbor_index = sv->index + (world->nz_parts-1)*(world->ny_parts-1) - (world->nz_parts-1) + 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
        new_sv = &(world->subvol[neighbor_index]);
        new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
        new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
        new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
        new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
        new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
        new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
        if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){

                /* find coordinates of the shared corner */
                v0.x = world->x_fineparts[sv->urb.x];
                v0.y = world->y_fineparts[sv->llf.y];
                v0.z = world->z_fineparts[sv->urb.z];

            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                
                /* check for the distance from the shared corner */
                vectorize(&(mp->pos), &v0,&v1);
                d = vect_length(&v1);
                if (d > R) continue;
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


        }
   
      }
   /* go (-X and -Y and -Z) */
      neighbor_index = sv->index - (world->nz_parts-1)*(world->ny_parts-1) - (world->nz_parts-1) - 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
        new_sv = &(world->subvol[neighbor_index]);
        new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
        new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
        new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
        new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
        new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
        new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
        if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){

                /* find coordinates of the shared corner */
                v0.x = world->x_fineparts[sv->llf.x];
                v0.y = world->y_fineparts[sv->llf.y];
                v0.z = world->z_fineparts[sv->llf.z];

            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  

                /* check for the distance from the shared corner */
                vectorize(&(mp->pos), &v0,&v1);
                d = vect_length(&v1);
                if (d > R) continue;

      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


        }
   
      }
   /* go (-X and +Y and -Z) */
      neighbor_index = sv->index - (world->nz_parts-1)*(world->ny_parts-1) + (world->nz_parts-1) - 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
        new_sv = &(world->subvol[neighbor_index]);
        new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
        new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
        new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
        new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
        new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
        new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
        if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){

                /* find coordinates of the shared corner */
                v0.x = world->x_fineparts[sv->llf.x];
                v0.y = world->y_fineparts[sv->urb.y];
                v0.z = world->z_fineparts[sv->llf.z];

            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                  
                /* check for the distance from the shared corner */
                vectorize(&(mp->pos), &v0,&v1);
                d = vect_length(&v1);
                if (d > R) continue;

      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


        }
   
      }
  /* go (+X and +Y and -Z) */
      neighbor_index = sv->index + (world->nz_parts-1)*(world->ny_parts-1) + (world->nz_parts-1) - 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
        new_sv = &(world->subvol[neighbor_index]);
        new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
        new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
        new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
        new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
        new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
        new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
        if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){

                /* find coordinates of the shared corner */
                v0.x = world->x_fineparts[sv->urb.x];
                v0.y = world->y_fineparts[sv->urb.y];
                v0.z = world->z_fineparts[sv->llf.z];

            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                  
                /* check for the distance from the shared corner */
                vectorize(&(mp->pos), &v0,&v1);
                d = vect_length(&v1);
                if (d > R) continue;

      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


        }
   
      }
  /* go (+X and -Y and -Z) */
      neighbor_index = sv->index + (world->nz_parts-1)*(world->ny_parts-1) - (world->nz_parts-1) - 1;
      if((neighbor_index > 0) && (neighbor_index < world->n_subvols))
      {
        new_sv = &(world->subvol[neighbor_index]);
        new_sv_llf.x = world->x_fineparts[new_sv->llf.x];
        new_sv_llf.y = world->y_fineparts[new_sv->llf.y];
        new_sv_llf.z = world->z_fineparts[new_sv->llf.z];
        new_sv_urb.x = world->x_fineparts[new_sv->urb.x];
        new_sv_urb.y = world->y_fineparts[new_sv->urb.y];
        new_sv_urb.z = world->z_fineparts[new_sv->urb.z];
        if(test_bounding_boxes(&path_llf, &path_urb, &new_sv_llf, &new_sv_urb)){

                /* find coordinates of the shared corner */
                v0.x = world->x_fineparts[sv->urb.x];
                v0.y = world->y_fineparts[sv->llf.y];
                v0.z = world->z_fineparts[sv->llf.z];


            /* include molecules from this SV */
    	    for (mp = new_sv->mol_head ; mp != NULL ; mp = mp->next_v)
    	    {
      		if (mp->properties == NULL) continue;  
                
                /* check for the distance from the shared corner */
                vectorize(&(mp->pos), &v0,&v1);
                d = vect_length(&v1);
                if (d > R) continue;
                  
      		/* check for the possible reaction. */
      		rx = trigger_bimolecular(
               		m->properties->hashval,mp->properties->hashval,
               		(struct abstract_molecule*)m,
                        (struct abstract_molecule*)mp,0,0);
      
      		if (rx != NULL)
      		{
        		smash = mem_get(sv->local_storage->coll);
        		if (smash == NULL)
			{
	  			fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  			i = emergency_output();
	  			fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",m->properties->sym->name);
	  			exit( EXIT_FAILURE );
        		}
        		smash->target = (void*) mp;
        		smash->intermediate = rx;
        		smash->next = shead1;
        		smash->what = COLLIDE_MOL;
        		shead1 = smash;
      		}
    	     }


        }
   
      }

  return shead1;
}


/*************************************************************************
diffuse_3D:
  In: molecule that is moving
      maximum time we can spend diffusing
      are we inert (nonzero) or can we react (zero)?
  Out: Pointer to the molecule if it still exists (may have been
       reallocated), NULL otherwise.
       Position and time are updated, but molecule is not rescheduled.
*************************************************************************/

struct molecule* diffuse_3D(struct molecule *m,double max_time,int inert)
{
  /*const double TOL = 10.0*EPS_C;*//* Two walls are coincident if this close */
  struct vector3 displacement;           /* Molecule moves along this vector */
  struct collision *smash;     /* Thing we've hit that's under consideration */
  struct collision *shead = NULL;        /* Things we might hit (can interact with) */
  struct collision *shead2 = NULL;     /* Things that we will hit, given our motion */ 
  struct subvolume *sv;
  struct wall *w;
  struct wall *reflectee;        /* Bounced off this one, don't hit it again */
  struct rxn *rx;
  struct molecule *mp,*old_mp;
  struct grid_molecule *g;
  struct abstract_molecule *am;
  struct species *sm;
  double steps=1.0;
  double t_steps=1.0;
  double factor;           /* return value from 'exact_disk()' function */
  double scaling = 1.0;          /* scales reaction cumulative_probabilitities array */
  double rate_factor=1.0;
  double f;

  /* this flag is set to 1 only after reflection from a wall and only with expanded lists. */
  int redo_expand_collision_list_flag = 0; 

  int i,j,k,l;
    
  int calculate_displacement = 1;
  
  sm = m->properties;
  if (sm==NULL) {
	fprintf(world->err_file,"File '%s', Line %ld: This molecule should not diffuse!\n", __FILE__, (long)__LINE__);
	return NULL;
  }
  if (sm->space_step <= 0.0)
  {
    m->t += max_time;
    return m;
  }
  
  if (!world->surface_reversibility)
  {
    if (m->flags&ACT_CLAMPED) /* Pretend we were already moving */
    {
      m->birthday -= 5*sm->time_step; /* Pretend to be old */
    }
    else
    {
      /* Newly created particles that have long time steps gradually increase */
      /* their timestep to the full value */
      if (sm->time_step > 1.0)
      {
        f = 1.0 + 0.2*(m->t - m->birthday);
        if (f<1) printf("I don't think so.\n");
        if (max_time > f) max_time=f;
      }
    }
  }
  
/* Even if we can't react, let's do a little bit of clean up. */
/* We'll just clear out anyone after us who is defunct and */
/* not worry about the whole list.  (Tom's idea, impl. by Rex) */

  while (m->next_v != NULL && m->next_v->properties == NULL)
  {
    mp = m->next_v;
    m->next_v = mp->next_v;
    
    if ((mp->flags & IN_MASK)==IN_VOLUME)
    {
      mem_put(mp->birthplace,mp);
    }
    else if ((mp->flags & IN_VOLUME) != 0)
    {
      mp->flags -=IN_VOLUME;
      mp->next_v = NULL;
    }
  }
  
/* Done housekeeping, now let's do something fun! */

  if ((sm->flags&COUNT_ENCLOSED)!=0) m->flags |= COUNT_ME;
  else if ((m->flags&COUNT_ME)!=0) m->flags -= COUNT_ME;
  
pretend_to_call_diffuse_3D:   /* Label to allow fake recursion */

  sv = m->subvol;
  
  shead = NULL;
  old_mp = NULL;
  if ( (sm->flags & (CAN_MOLMOL | CANT_INITIATE)) == CAN_MOLMOL )
  {
    for (mp = sv->mol_head ; mp != NULL ; old_mp = mp , mp = mp->next_v)
    {
continue_special_diffuse_3D:   /* Jump here instead of looping if old_mp,mp already set */

      if (mp==m) continue;
      
      if (mp->properties == NULL)  /* Reclaim storage */
      {
        if (old_mp != NULL) old_mp->next_v = mp->next_v;
        else sv->mol_head = mp->next_v;
        
        if ((mp->flags & IN_MASK)==IN_VOLUME) mem_put(mp->birthplace,mp);
        else if ((mp->flags & IN_VOLUME) != 0) mp->flags -= IN_VOLUME;
        
        mp = mp->next_v;

        if (mp==NULL) break;
        else goto continue_special_diffuse_3D;  /*continue without incrementing pointer*/
      }
      
      rx = trigger_bimolecular(
               sm->hashval,mp->properties->hashval,
               (struct abstract_molecule*)m,(struct abstract_molecule*)mp,0,0
             );
      
      if (rx != NULL)
      {
        smash = mem_get(sv->local_storage->coll);
        if (smash == NULL)
	{
	  fprintf(world->err_file,"File '%s', Line %ld: out of memory.  Trying to save intermediate states.\n", __FILE__, (long)__LINE__);
	  i = emergency_output();
	  fprintf(world->err_file,"Out of memory while finding collisions for a molecule of type %s\n",sm->sym->name);
	  exit( EXIT_FAILURE );
        }
        smash->target = (void*) mp;
        smash->what = COLLIDE_MOL;
        smash->intermediate = rx;
        smash->next = shead;
        shead = smash;
      }
/*      else printf("Rx between %s and %s is NULL\n",sm->sym->name,mp->properties->sym->name); */
    }
  }
  
  if (calculate_displacement)
  {
    if (m->flags&ACT_CLAMPED)
    {
      steps = rng_dbl(world->rng);
      t_steps = sm->time_step;
      pick_clamped_displacement(&displacement,m);
      m->flags-=ACT_CLAMPED;
      m->previous_wall=NULL;
      m->index=-1;
      rate_factor=1.0;
    }
    else
    {
      if (max_time > MULTISTEP_WORTHWHILE) steps = safe_diffusion_step(m,shead);
      else steps = 1.0;
   
      t_steps = steps * sm->time_step;
      if (t_steps > max_time)
      {
        t_steps = max_time;
        steps = max_time / sm->time_step;
      }
      if (steps < EPS_C)
      {
        steps = EPS_C;
        t_steps = EPS_C*sm->time_step;
      }
      
      if (steps == 1.0)
      {
        pick_displacement(&displacement,sm->space_step);
        rate_factor = 1.0;
      }
      else
      {
        rate_factor = sqrt(steps);
        pick_displacement(&displacement,rate_factor*sm->space_step);
      }
    }

    world->diffusion_number += 1.0;
    world->diffusion_cumsteps += steps;
  }
  
  reflectee = NULL;
  
  if(world->use_expanded_list && (m->properties->flags & CAN_MOLMOL) != 0)
  { 
    shead = expand_collision_list(m, &displacement, sv, shead);
  }   

#define CLEAN_AND_RETURN(x) if (shead2!=NULL) mem_put_list(sv->local_storage->coll,shead2); if (shead!=NULL) mem_put_list(sv->local_storage->coll,shead); return (x)
#define ERROR_AND_QUIT fprintf(world->err_file,"File '%s', Line %ld: out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__); i=emergency_output(); fprintf(world->err_file,"Fatal error: out of memory during diffusion of a %s molecule\nAttempt to write intermediate results had %d errors\n",sm->sym->name,i); exit(EXIT_FAILURE)
  do
  {
    if(world->use_expanded_list && redo_expand_collision_list_flag)
    {
      if((m->properties->flags & CAN_MOLMOL) != 0){ 
    	  shead = expand_collision_list(m, &displacement, sv, shead);
      }  
    }
 
    shead2 = ray_trace(m,shead,sv,&displacement,reflectee);
    
    if (shead2==NULL) { ERROR_AND_QUIT; }

    if (shead2->next!=NULL)  /* Could be sped up/combined */
    {
      shead2 = (struct collision*)ae_list_sort((struct abstract_element*)shead2);
/*      shead2 = gather_walls_first(shead2,TOL); */
    }
    
    for (smash = shead2; smash != NULL; smash = smash->next)
    {
      
      if (smash->t >= 1.0 || smash->t < 0.0)
      {
	if ((smash->what&COLLIDE_MOL)!=0) printf("File '%s', Line %ld: Unexpected behavior. Iteration %lld, time of collision %.8e\n", __FILE__, (long)__LINE__,  world->it_time,smash->t);
        smash = NULL;
        break;
      }
      
      rx = smash->intermediate;

      if ( (smash->what & COLLIDE_MOL) != 0 && !inert )
      {
	if (smash->t < EPS_C) continue;
	
	world->mol_mol_colls++;

        am = (struct abstract_molecule*)smash->target;
        if ((am->flags & ACT_INERT) != 0)  /* FIXME */
        {
          /*if (smash->t < am->t + am->t2) continue;*/
        }  
	factor = exact_disk(
          &(smash->loc),&displacement,world->rx_radius_3d,m->subvol,m,
	  (struct molecule*)am
        );
      
	if (factor<0) /* Probably hit a wall, might have run out of memory */
	{
	  if (factor==EXD_OUT_OF_MEMORY) { ERROR_AND_QUIT; }
	  continue; /* Reaction blocked by a wall */
	}
        
        scaling = factor / rate_factor;
        if (rx->prob_t != NULL) check_probs(rx,m->t);

        i = test_bimolecular(rx,scaling);
        
        if (i < RX_LEAST_VALID_PATHWAY) continue;
	
        j = outcome_bimolecular(
                rx,i,(struct abstract_molecule*)m,
                am,0,0,m->t+t_steps*smash->t,&(smash->loc)
              );
	      
	if (j==RX_NO_MEM) { ERROR_AND_QUIT; }
	if (j!=RX_DESTROY) continue;
        else { CLEAN_AND_RETURN( NULL ); }
      }

      else if ( (smash->what & COLLIDE_WALL) != 0 )
      {
	w = (struct wall*) smash->target;
	
	world->ray_polygon_colls++;
	
	if ( (smash->what & COLLIDE_MASK) == COLLIDE_FRONT ) k = 1;
	else k = -1;
	
	if ( w->effectors != NULL && (sm->flags&CAN_MOLGRID) != 0)
	{
	  j = xyz2grid( &(smash->loc) , w->effectors );
	  if (w->effectors->mol[j] != NULL)
	  {
	    if (m->index != j || m->previous_wall != w )
	    {
	      g = w->effectors->mol[j];
	      rx = trigger_bimolecular(
		sm->hashval,g->properties->hashval,
		(struct abstract_molecule*)m,(struct abstract_molecule*)g,
		k,g->orient
	      );
	      if (rx!=NULL)
	      {
		if (rx->prob_t != NULL) check_probs(rx,m->t);
                scaling = 1.0 / (rate_factor * w->effectors->binding_factor);
		
		i = test_bimolecular(rx, scaling);
		
                if (i >= RX_LEAST_VALID_PATHWAY)
		{
		  l = outcome_bimolecular(
		    rx,i,(struct abstract_molecule*)m,
		    (struct abstract_molecule*)g,
		    k,g->orient,m->t+t_steps*smash->t,&(smash->loc)
		  );
		  
		  if (l==RX_NO_MEM) { ERROR_AND_QUIT; }
		  if (l==RX_FLIP)
		  {
		    if ( (sm->flags & w->flags & COUNT_HITS) )
		    {
		      update_collision_count(sm,w->counting_regions,k,1,rate_factor * w->effectors->binding_factor,&(smash->loc),smash->t);
		    }
		    if ((m->flags&COUNT_ME)!=0)
		    {
		      m->flags-=COUNT_ME;
		      count_me_by_region((struct abstract_molecule*)m,-1,NULL,smash->t);
		    }
		    
		    continue; /* pass through */
		  }
		  else if (l==RX_DESTROY)
		  {
		    if ( (sm->flags & w->flags & COUNT_HITS) )
		      update_collision_count(sm,w->counting_regions,k,0,rate_factor,&(smash->loc),smash->t);
		    
		    CLEAN_AND_RETURN(NULL);
		  }
		}
	      }
	    }
	    else m->index = -1; /* Avoided rebinding, but next time it's OK */
	  }
	}
	
	if ( (sm->flags&CAN_MOLWALL) != 0 )
	{
	  m->index = -1;
	  rx = trigger_intersect(
		  sm->hashval,(struct abstract_molecule*)m,k,w
		);
	  
	  if (rx != NULL)
	  {
	    if (rx->n_pathways == RX_TRANSP)
	    {
	      rx->n_occurred++;
	      if ( (sm->flags & COUNT_HITS) )
	      {
		update_collision_count(sm,w->counting_regions,k,1,rate_factor,&(smash->loc),smash->t);
	      }
	      if ((m->flags&COUNT_ME)!=0)
	      {
		m->flags-=COUNT_ME;
		count_me_by_region((struct abstract_molecule*)m,-1,NULL,smash->t);
	      }

	      continue; /* Ignore this wall and keep going */
	    }
	    else if (rx->n_pathways != RX_REFLEC)
	    {
	      if (rx->prob_t != NULL) check_probs(rx,m->t);
	      i = test_intersect(rx,1.0/rate_factor);
	      if (i > RX_NO_RX)
	      {
		j = outcome_intersect(
			rx,i,w,(struct abstract_molecule*)m,
			k,m->t + t_steps*smash->t,&(smash->loc)
		      );
		      
		if (j==RX_NO_MEM) { ERROR_AND_QUIT; } 
		if (j==RX_FLIP)
		{
		  if ( (sm->flags & COUNT_HITS) )
		  {
		    update_collision_count(sm,w->counting_regions,k,1,rate_factor,&(smash->loc),smash->t);
		  }
		  if ((m->flags&COUNT_ME)!=0)
		  {
		    m->flags-=COUNT_ME;
		    count_me_by_region((struct abstract_molecule*)m,-1,NULL,smash->t);
		  }
  
		  continue; /* pass through */
		}
		else if (j==RX_DESTROY)
		{
		  if ( (sm->flags & COUNT_HITS) )
		    update_collision_count(sm,w->counting_regions,k,0,rate_factor,&(smash->loc),smash->t);
  
		  CLEAN_AND_RETURN(NULL);
		}
	      }
	    }
	    /* RX_REFLEC will just fall through here and be reflected by default case below */
	  }
	}
	
        /* default is to reflect */
        
	if ( (sm->flags & COUNT_HITS)) update_collision_count(sm,w->counting_regions,k,0,rate_factor,&(smash->loc),smash->t);
	if (m->flags&COUNT_ME)
	{
	  m->flags -= COUNT_ME;
	  count_me_by_region((struct abstract_molecule*)m,-1,NULL,smash->t);
	}

	m->pos.x = smash->loc.x;
	m->pos.y = smash->loc.y;
	m->pos.z = smash->loc.z;
        m->t += t_steps*smash->t;
	reflectee = w;

        t_steps *= (1.0-smash->t);
        
        factor = -2.0 * (displacement.x*w->normal.x + displacement.y*w->normal.y + displacement.z*w->normal.z);
        displacement.x = (displacement.x + factor*w->normal.x) * (1.0-smash->t);
        displacement.y = (displacement.y + factor*w->normal.y) * (1.0-smash->t);
        displacement.z = (displacement.z + factor*w->normal.z) * (1.0-smash->t);

        redo_expand_collision_list_flag = 0; /* Only useful if we're using expanded lists, but easier to always set it */
        
        break;
      }
      else if ((smash->what & COLLIDE_SUBVOL) != 0)
      {
        struct subvolume *nsv;

        m->pos.x = smash->loc.x;
        m->pos.y = smash->loc.y;
        m->pos.z = smash->loc.z;
        
        displacement.x *= (1.0-smash->t);
        displacement.y *= (1.0-smash->t);
        displacement.z *= (1.0-smash->t);
        
        m->t += t_steps*smash->t;
        t_steps *= (1.0-smash->t);
        if (t_steps < EPS_C) t_steps = EPS_C;

        nsv = traverse_subvol(sv,&(m->pos),smash->what - COLLIDE_SV_NX - COLLIDE_SUBVOL);
        if (nsv==NULL)
        {
          fprintf(world->log_file,"Error: a %s molecule escaped the world at (%.2e,%.2e,%.2e)\n",
                  sm->sym->name,m->pos.x*world->length_unit,
                  m->pos.y*world->length_unit,m->pos.z*world->length_unit);
          if ((sm->flags&COUNT_CONTENTS)!=0 && (m->flags&COUNT_ME)!=0)
	    count_me_by_region((struct abstract_molecule*)m,-1,NULL,smash->t);
          sm->population--;
          m->properties = NULL;
	  
	  CLEAN_AND_RETURN(NULL);
        }
        else m = migrate_molecule(m,nsv);

        if (shead2 != NULL) mem_put_list(sv->local_storage->coll,shead2);
        if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);
        calculate_displacement = 0;
        
        if (m->properties==NULL) fprintf(world->err_file,"File '%s', Line %ld: This molecule should not be jumping.\n", __FILE__, (long)__LINE__);
        goto pretend_to_call_diffuse_3D;  /* Jump to beginning of function */        
      }
    }
    
    if (shead2 != NULL) mem_put_list(sv->local_storage->coll,shead2);

  }
  while (smash != NULL);

#undef ERROR_AND_QUIT
#undef CLEAN_AND_RETURN
  
  m->pos.x += displacement.x;
  m->pos.y += displacement.y;
  m->pos.z += displacement.z;
  m->t += t_steps;
  
  m->index = -1;
  m->previous_wall=NULL;
  if ((sm->flags&COUNT_ENCLOSED)!=0 && (m->flags&COUNT_ME)==0) count_me_by_region((struct abstract_molecule*)m,1,NULL,m->t);
  
  if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);

  return m;
}



/*************************************************************************
diffuse_2D:
  In: molecule that is moving
      maximum time we can spend diffusing
  Out: Pointer to the molecule, or NULL if there was an error (right now
       there is no reallocation)
       Position and time are updated, but molecule is not rescheduled,
       nor does it react
*************************************************************************/

struct grid_molecule* diffuse_2D(struct grid_molecule *g,double max_time)
{
  struct species *sg;
  struct vector2 displacement,new_loc;
  struct wall *new_wall;
  double f;
  double steps,t_steps;
  double space_factor;
  int find_new_position,new_idx;
  
  sg = g->properties;
  
  if (sg==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Error!  Surface molecule has no properties?  Ignoring!\n", __FILE__, (long)__LINE__);
    return NULL;
  }
  
  if (sg->space_step <= 0.0)
  {
    g->t += max_time;
    return g;
  }
  
  if (sg->time_step > 1.0)
  {
    f = 1.0 + 0.2*(g->t - g->birthday);
    if (f<1) fprintf(world->log_file, "File '%s', Line %ld: Unexpected behavior.\n", __FILE__, (long)__LINE__);
    if (max_time>f) max_time=f;
  }
  
  if ((sg->flags&COUNT_ENCLOSED) != 0) g->flags |= COUNT_ME;
  else if ((g->flags&COUNT_ME) != 0) g->flags -= COUNT_ME;
  
  /* Where are we going? */
  if (sg->time_step > max_time)
  {
    t_steps = max_time;
    steps = max_time / sg->time_step;
  }
  else
  {
    t_steps = sg->time_step;
    steps = 1.0;
  }
  if (steps < EPS_C)
  {
    steps = EPS_C;
    t_steps = EPS_C*sg->time_step;
  }
  
  if (steps==1.0) space_factor = sg->space_step;
  else space_factor = sg->space_step*sqrt(steps);
  
  world->diffusion_number += 1.0;
  world->diffusion_cumsteps += steps;
  
  for (find_new_position=(SURFACE_DIFFUSION_RETRIES+1) ; find_new_position > 0 ; find_new_position--)
  {
    pick_2d_displacement(&displacement,space_factor);
    
    new_wall = ray_trace_2d(g,&displacement,&new_loc);
    
    if (new_wall==NULL) continue;  /* Something went wrong--try again */
    
    if (new_wall == g->grid->surface)
    {
      new_idx = uv2grid(&new_loc,new_wall->effectors);
      if (new_idx < 0 || new_idx >= g->grid->n_tiles)
      {
	fprintf(world->log_file, "File '%s', Line %ld: Unexpected behaviour, iteration %d.\n", __FILE__, (long)__LINE__, (int)world->it_time);
      }
      if (new_idx != g->grid_index)
      {
	if (g->grid->mol[new_idx]!=NULL)
	{
	  continue; /* Pick again--full here */
	}
	
	count_me_by_region((struct abstract_molecule*)g,-1,NULL,g->t);

	g->grid->mol[g->grid_index]=NULL;
	g->grid->mol[new_idx] = g;
	g->grid_index = new_idx;
      }
      else count_me_by_region((struct abstract_molecule*)g,-1,NULL,g->t);
      
      g->s_pos.u = new_loc.u;
      g->s_pos.v = new_loc.v;
      count_me_by_region((struct abstract_molecule*)g,1,NULL,g->t);
      
      find_new_position = 0;
    }
    else 
    {
      if (new_wall->effectors==NULL)
      { 
	if (create_grid(new_wall,NULL))
	{
	  fprintf(world->err_file,"File '%s', Line %ld: Failed to create surface grid for diffusing molecule.\n", __FILE__, (long)__LINE__);
	  return NULL;
	}
      }

      /* Move to new tile */
      new_idx = uv2grid(&new_loc,new_wall->effectors);
      if (new_idx < 0 || new_idx >= new_wall->effectors->n_tiles) fprintf(world->log_file, "File '%s', Line %ld: Unexpected behaviour, iteration %d.\n", __FILE__, (long)__LINE__, (int)world->it_time);
      if (new_wall->effectors->mol[new_idx] != NULL) continue; /* Pick again */
      
      count_me_by_region((struct abstract_molecule*)g,-1,NULL,g->t);
      g->grid->mol[g->grid_index]=NULL;
      g->grid->n_occupied--;
      g->grid = new_wall->effectors;
      g->grid_index=new_idx;
      g->grid->mol[new_idx] = g;

      g->s_pos.u = new_loc.u;
      g->s_pos.v = new_loc.v;
      count_me_by_region((struct abstract_molecule*)g,1,NULL,g->t);
      
      find_new_position=0;
    }
  }
  
  g->t += t_steps;
  return g;
}
    

/*************************************************************************
react_2D:
  In: molecule that may react
      maximum duration we have to react
  Out: Pointer to the molecule if it still exists (may have been
       destroyed), NULL otherwise.
  Note: Time is not updated--assume that's already taken care of
        elsewhere.  Only nearest neighbors can react.
*************************************************************************/

struct grid_molecule* react_2D(struct grid_molecule *g,double t)
{
  struct surface_grid *sg[3];    /* Neighboring surface grids */
  int si[3];                     /* Indices on those grids of neighbor molecules */
  struct grid_molecule *gm[3];   /* Neighboring molecules */
  struct rxn *rx[3];             /* Reactions we can perform with those molecules */
  double cf[3];                  /* Correction factor for area for those molecules */
  int i,j,k,n;
  long long ii;
  
  grid_neighbors(g->grid,g->grid_index,sg,si);
  
  for (i=n=0 ; i<3 ; i++)
  {
    if (sg[i]!=NULL)
    {
      gm[n] = sg[i]->mol[ si[i] ];
      if (gm[n]!=NULL)
      {
	rx[n] = trigger_bimolecular(
	  g->properties->hashval,gm[n]->properties->hashval,
	  (struct abstract_molecule*)g,(struct abstract_molecule*)gm[n],
	  g->orient,gm[n]->orient
	);
	if (rx[n]!=NULL)
	{
	  cf[n] = sg[i]->binding_factor/t;
	  n++;
	}
      }
    }
  }
  
  if (n==0) return g;  /* Nobody to react with */
  else if (n==1)
  {
    i = test_bimolecular(rx[0],cf[0]);
    j = 0;
  }
  else
  {
    ii = test_many_bimolecular(rx,cf,n);
    if (ii>=RX_LEAST_VALID_PATHWAY)
    {
      i = (int)(ii & RX_GREATEST_VALID_PATHWAY);
      j = (int)(ii >> RX_PATHWAY_BITS);
    }
    else
    {
      i = (int)ii;
      j = 0;
    }
  }
  
  if (i<RX_LEAST_VALID_PATHWAY) return g;  /* No reaction */
  
  k = outcome_bimolecular(
    rx[j],i,
    (struct abstract_molecule*)g,(struct abstract_molecule*)gm[j],
    g->orient,gm[j]->orient,g->t,NULL
  );

  if (k==RX_NO_MEM)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory.  Trying to save intermediate results.\n", __FILE__, (long)__LINE__);
    k = emergency_output();
    fprintf(world->err_file,"Out of memory during bimolecular surface reaction %s...\n",rx[j]->sym->name);
    fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",k);
    exit( EXIT_FAILURE );
  }
  
  if (j==RX_DESTROY)
  {
    mem_put(g->birthplace,g);
    return NULL;
  }
  
  return g;
}


/*************************************************************************
run_timestep:
  In: local storage area to use
      time of the next release event
      time of the next checkpoint
  Out: No return value.  Every molecule in the subvolume is updated in
       position and rescheduled at least one timestep ahead.  The
       current_time of the subvolume is also incremented.
  Note: This also occasionally does garbage collection on the scheduling
        queue.
*************************************************************************/

void run_timestep(struct storage *local,double release_time,double checkpt_time)
{
  struct abstract_molecule *a;
  struct wall *w;
  struct rxn *r,*r2;
  double t,tt;
  double max_time;
  int i,j,err,special;
  
  /* Check for garbage collection first */
  if ( local->timer->defunct_count > MIN_DEFUNCT_FOR_GC &&
       MAX_DEFUNCT_FRAC*(local->timer->count) < local->timer->defunct_count )
  {
    struct abstract_molecule *temp;
    j=0;
    i=local->timer->defunct_count;
    a = (struct abstract_molecule*) schedule_cleanup(local->timer , *is_defunct_molecule);
    while (a!=NULL)
    {
      temp = a;
      a = a->next;
/*      if (temp->properties!=NULL) fprintf(world->err_file,"Removed a non-defunct molecule from scheduler!\n"); */
      if ((temp->flags&IN_MASK)==IN_SCHEDULE)
      {
	temp->next = NULL;
	mem_put(temp->birthplace,temp);
	j++;
      }
      else temp->flags -= IN_SCHEDULE;
    }
/*    fprintf(world->log_file,"Cleaning up memory: removed %d (actually only %d) unused molecules.\n",i,j); */
  }
  /* Now run the timestep */
  while ( (a = (struct abstract_molecule*)schedule_next(local->timer)) != NULL )
  {
    if (a->properties == NULL)  /* Defunct!  Remove molecule. */
    {
      if ((a->flags & IN_MASK) == IN_SCHEDULE)
      {
        a->next = NULL; 
        mem_put(a->birthplace,a);
      }
      else a->flags -= IN_SCHEDULE;
      if (local->timer->defunct_count>0) local->timer->defunct_count--;
      
      continue;
    }
    
    a->flags -= IN_SCHEDULE;
    
    /* Check for a unimolecular event */
    if (a->t2 < EPS_C || a->t2 < EPS_C*a->t)
    {
      if ((a->flags & (ACT_INERT+ACT_NEWBIE+ACT_CHANGE)) != 0)
      {
        a->flags -= (a->flags & (ACT_INERT + ACT_NEWBIE + ACT_CHANGE));
        if ((a->flags & ACT_REACT) != 0)
        {
          r = trigger_unimolecular(a->properties->hashval,a);
	  if (r!=NULL)
	  {
            if (r->prob_t != NULL) check_probs(r,(a->t + a->t2)*(1.0+EPS_C));
	  }
	  
	  tt=FOREVER; /* When will rates change? */
	  
	  r2=NULL;
	  if (a->properties->flags&CAN_GRIDWALL) r2=trigger_surface_unimol(a,NULL);
	  if ( r2!=NULL )
	  {
	    if (r2->prob_t != NULL) check_probs(r2,(a->t + a->t2)*(1.0+EPS_C));
	    a->t2 = (r==NULL) ? timeof_unimolecular(r2) : timeof_special_unimol(r,r2);
	    if (r!=NULL && r->prob_t!=NULL) tt = r->prob_t->time;
	    if (r2->prob_t!=NULL && tt > r2->prob_t->time) tt = r2->prob_t->time;
	  }
	  else if (r!=NULL)
	  {
	    a->t2 = timeof_unimolecular(r);
	    if (r->prob_t!=NULL) tt = r->prob_t->time;
	  }
	  else a->t2 = FOREVER;
	  
	  if (a->t + a->t2 > tt)
	  {
	    a->t2 = tt - a->t;
	    a->flags |= ACT_CHANGE;
	  }
        }
      }
      else if ((a->flags & ACT_REACT) != 0)
      {
	special = 0;
	r2 = NULL;
        r = trigger_unimolecular(a->properties->hashval,a);

	if (a->properties->flags&CAN_GRIDWALL)
	{
	  r2 = trigger_surface_unimol(a,NULL);
	  if (r2!=NULL)
	  {
	    special = 1;
	    if (r==NULL || is_surface_unimol(r,r2))
	    {
	      special = 2;
	      r = r2; /* Do surface-limited rx instead */
	    }
	  }
	}
	
	if (r!=NULL)
	{
	  i = which_unimolecular(r);
	  j = outcome_unimolecular(r,i,a,a->t);
	  if (j==RX_NO_MEM)
	  {
	    fprintf(world->err_file,"File '%s', Line %ld: Out of memory.  Trying to save intermediate results.\n", __FILE__, (long)__LINE__);
	    i = emergency_output();
	    fprintf(world->err_file,"Out of memory during unimolecular reaction %s...\n",r->sym->name);
	    fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",i);
	    exit( EXIT_FAILURE );
	  }
	}
	else j=RX_NO_RX;
	
        if (j!=RX_DESTROY) /* We still exist */
        {
	  tt = FOREVER;
	  if (special)
	  {
	    if (special==2)
	    {
	      r2 = r;
	      r = trigger_unimolecular(a->properties->hashval,a);
	    }
	    a->t2 = (r==NULL) ? timeof_unimolecular(r2) : timeof_special_unimol(r,r2);
	    if (r!=NULL && r->prob_t!=NULL) tt=r->prob_t->time;
	    if (r2!=NULL && r2->prob_t!=NULL && r2->prob_t->time < tt) tt = r2->prob_t->time;
	  }
	  else if (r!=NULL)
	  {
	    a->t2 = timeof_unimolecular(r);
	    if (r->prob_t != NULL) tt=r->prob_t->time;
	  }
	  else a->t2 = FOREVER;
	  
	  if (a->t + a->t2 > tt)
	  {
	    a->t2 = tt - a->t;
	    a->flags |= ACT_CHANGE;
	  }
        }
        else /* We don't exist.  Try to recover memory. */
	{
	  if (!(a->flags&IN_VOLUME)) mem_put(a->birthplace,a);
	  continue;
	}
      }
    }
    
    t = a->t;

    if ((a->flags & ACT_DIFFUSE) != 0)
    {
      max_time = checkpt_time - a->t;
      if (local->max_timestep < max_time) max_time = local->max_timestep;
      if (a->t2<max_time && (a->flags&(ACT_REACT|ACT_INERT))!=0) max_time = a->t2;

      if ((a->flags & TYPE_3D) != 0)
      {
        if (max_time > release_time - a->t) max_time = release_time - a->t;
        a = (struct abstract_molecule*)diffuse_3D((struct molecule*)a , max_time , a->flags & ACT_INERT);
        if (a!=NULL) /* We still exist */
        {
          a->t2 -= a->t - t;
        }
        else continue;
      }
      else
      {
	if (max_time > release_time - a->t) max_time = release_time - a->t;
	w = ((struct grid_molecule*)a)->grid->surface;
	a = (struct abstract_molecule*)diffuse_2D((struct grid_molecule*)a , max_time);
	
	if (a!=NULL)
	{
	  if ( (a->properties->flags&CAN_GRIDWALL)==0 ||
	       w==((struct grid_molecule*)a)->grid->surface )
	  {
	    a->t2 -= a->t - t;
	  }
	  else if (w->surf_class==((struct grid_molecule*)a)->grid->surface->surf_class)
	  {
	    a->t2 -= a->t - t;
	  }
	  else
	  {
	    a->t2 = 0;
	    a->flags |= ACT_CHANGE; /* Reschedule reaction time */
	  }
	}
	else continue;
      }
    }
    
    if ( (a->flags&TYPE_GRID)!=0 && (a->properties->flags&CAN_GRIDGRID) && !(a->flags&ACT_INERT))
    {
      if ((a->flags&ACT_DIFFUSE)==0) /* Didn't move, so we need to figure out how long to react for */
      {
	max_time = checkpt_time - a->t;
	if (a->t2<max_time && (a->flags &(ACT_REACT|ACT_INERT))!=0) max_time = a->t2;
	if (max_time > release_time - a->t) max_time = release_time - a->t;
	if (a->properties->time_step < max_time) max_time = a->properties->time_step;
      }
      else max_time = a->t - t;
       
      a = (struct abstract_molecule*)react_2D((struct grid_molecule*)a , max_time );
      
      if (a==NULL) continue;
      
      if ((a->flags&ACT_DIFFUSE)==0) /* Advance time if diffusion hasn't already done it */
      {
       a->t2 -= max_time;
       a->t += max_time;
      }
    }
    else if ((a->flags&ACT_DIFFUSE)==0)
    {
      if (a->t2==0) a->t += MAX_UNI_TIMESKIP;
      else 
      {
        a->t += a->t2;
        a->t2 = 0;
      }
    }

    a->flags += IN_SCHEDULE;
    
    /* If we're near an integer boundary, advance to the next integer */
    t = ceil(a->t)*(1.0+0.1*EPS_C);
    if (!distinguishable(t,a->t,EPS_C)) a->t=t;
    
    if (a->flags & TYPE_GRID) err = schedule_add(local->timer,a);
    else err = schedule_add(((struct molecule*)a)->subvol->local_storage->timer,a);
    
    if (err)
    {
      fprintf(world->err_file,"File '%s', Line %ld: Out of memory.  Trying to save intermediate results.\n", __FILE__, (long)__LINE__);
      i = emergency_output();
      fprintf(world->err_file,"Out of memory while scheduling molecule of type %s\n",a->properties->sym->name);
      fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",i);
      exit( EXIT_FAILURE );
    }
  }
  if (local->timer->error)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory.  Trying to save intermediate results.\n", __FILE__, (long)__LINE__);
    i = emergency_output();
    fprintf(world->err_file,"Out of memory while retrieving molecules to move.\n");
    fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",i);
    exit( EXIT_FAILURE );
  }
  
  local->current_time += 1.0;
}


/*************************************************************************
run_concentration_clamp:
  In: The current time.
  Out: No return value.  Molecules are released at concentration-clamped
       surfaces to maintain the desired concentation.
*************************************************************************/

void run_concentration_clamp(double t_now)
{
  struct ccn_clamp_data *ccd;
  struct ccn_clamp_data *ccdo;
  struct ccn_clamp_data *ccdm;
  double n_collisions;
  int n_emitted,idx,i;
  struct wall *w;
  struct vector3 v;
  double s1,s2,eps;
  struct molecule m;
  struct molecule *mp;
  
  int this_count = 0;
  static int total_count = 0;
  
  for (ccd=world->clamp_list ; ccd!=NULL ; ccd=ccd->next)
  {
    if (ccd->objp==NULL) continue;
    for (ccdo=ccd ; ccdo!=NULL ; ccdo=ccdo->next_obj)
    {
      for (ccdm=ccdo ; ccdm!=NULL ; ccdm=ccdm->next_mol)
      {
        n_collisions = ccdo->scaling_factor * ccdm->mol->space_step * 
                       ccdm->concentration / ccdm->mol->time_step;
        n_emitted = poisson_dist( n_collisions , rng_dbl(world->rng) );
        
        if (n_emitted==0) continue;
        
        m.t = t_now+0.5;
        m.t2 = 0;
        m.flags = IN_SCHEDULE | ACT_NEWBIE | TYPE_3D | IN_VOLUME | ACT_CLAMPED | ACT_DIFFUSE;
        m.properties = ccdm->mol;
        m.birthplace=NULL;
        m.birthday = t_now;
        m.subvol=NULL;
        m.previous_wall=NULL;
        m.index=0;
        mp = NULL;
        
        this_count+=n_emitted;
        while (n_emitted>0)
        {
          idx = bisect_high(ccdo->cum_area,ccdo->n_sides,rng_dbl(world->rng)*ccdo->cum_area[ccd->n_sides-1]);
          w = ccdo->objp->wall_p[ ccdo->side_idx[idx] ];
          
          s1 = sqrt(rng_dbl(world->rng));
          s2 = rng_dbl(world->rng)*s1;
          
          v.x = w->vert[0]->x + s1*(w->vert[1]->x - w->vert[0]->x) + s2*(w->vert[2]->x - w->vert[1]->x);
          v.y = w->vert[0]->y + s1*(w->vert[1]->y - w->vert[0]->y) + s2*(w->vert[2]->y - w->vert[1]->y);
          v.z = w->vert[0]->z + s1*(w->vert[1]->z - w->vert[0]->z) + s2*(w->vert[2]->z - w->vert[1]->z);
          
          if (ccdm->orient==1) m.index=1;
          else if (ccdm->orient==-1) m.index=-1;
          else if (rng_uint(world->rng)&1) m.index=-1;
          else m.index=1;
          
          eps = EPS_C*m.index;
          
          if (v.x>0) s1=v.x; else s1=-v.x;
          if (v.y>0) s2=v.y; else s2=-v.y;
          if (s1<s2) s1=s2;
          if (v.z>0) s2=v.z; else s2=-v.z;
          if (s1<s2) s1=s2;
          if (s1>1.0) eps *= s1;
          
          m.pos.x = v.x + w->normal.x*eps;
          m.pos.y = v.y + w->normal.y*eps;
          m.pos.z = v.z + w->normal.z*eps;
          m.previous_wall = w;
          
          if (mp==NULL)
          {
            mp = insert_molecule(&m,mp);
            if (mp==NULL)
            {
              i = emergency_output();
              fprintf(world->err_file,"File '%s', Line %ld: Out of memory while concentration clamping molecule of type %s\n", __FILE__, (long)__LINE__, m.properties->sym->name);
              fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",i);
              exit( EXIT_FAILURE );
            }
            if (trigger_unimolecular(ccdm->mol->hashval , (struct abstract_molecule*)mp) != NULL)
            {
              m.flags |= ACT_REACT;
              mp->flags |= ACT_REACT;
            }
          }
          else
          {
            mp=insert_molecule(&m,mp);
            if (mp==NULL)
            {
              i = emergency_output();
              fprintf(world->err_file,"File '%s', Line %ld: Out of memory while concentration clamping molecule of type %s\n",__FILE__, (long)__LINE__, m.properties->sym->name);
              fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",i);
              exit( EXIT_FAILURE );
            }
          }
          
          n_emitted--;
        }
      }
    }
  }
  
  total_count += this_count;
//  printf("Emitted %d\n",total_count);
}

