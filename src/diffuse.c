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

/*#define USE_EXPANDED_COLLISION_LIST */     

/* SOLVE_QF is a local #define in exact_disk (solves the quadratic formula) */
/* REGISTER_PARTITION is a local #define in exact_disk */
/* CHECK_PARTITON is a local #define in exact_disk */

/* CLEAN_AND_RETURN(x) is a local #define in diffuse_3D */
/* ERROR_AND_QUIT is a local #define in diffuse_3D */


extern struct volume *world;

/*************************************************************************
pick_displacement:
  In: vector3 to store the new displacement
      scale factor to apply to the displacement
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         molecule, scaled by the scaling factor.
*************************************************************************/

void pick_displacement(struct vector3 *v,double scale)
{
  double r;
  
#if 0
  if (world->fully_random)
  {
    double x,y,z;
    do
    {
      x = 2.0*rng_dbl(world->rng)-1.0;
      y = 2.0*rng_dbl(world->rng)-1.0;
      z = 2.0*rng_dbl(world->rng)-1.0;
    } while (x*x + y*y + z*z >= 1.0 || x*x + y*y + z*z < 0.001);
    r = scale * world->r_step[ rng_uint(world->rng) & (world->radial_subdivisions-1) ] / sqrt(x*x + y*y + z*z);
    v->x = r*x;
    v->y = r*y;
    v->z = r*z;
    return;
  }
#endif
    
  
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
  double dx,dy,dz,tx,ty,tz;
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

  if (v->x < 0.0)
  {
    dx = world->x_fineparts[ sv->llf.x ] - m->pos.x;
    i = 0;
  }
  else
  {
    dx = world->x_fineparts[ sv->urb.x ] - m->pos.x;
    i = 1;
  }

  if (v->y < 0.0)
  {
    dy = world->y_fineparts[ sv->llf.y ] - m->pos.y;
    j = 0;
  }
  else
  {
    dy = world->y_fineparts[ sv->urb.y ] - m->pos.y;
    j = 1;
  }

  if (v->z < 0.0)
  {
    dz = world->z_fineparts[ sv->llf.z ] - m->pos.z;
    k = 0;
  }
  else
  {
    dz = world->z_fineparts[ sv->urb.z ] - m->pos.z;
    k = 1;
  }
  
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



#if 0
/* This function will do what estimate_disk does but do an exact computation.*/
/* This function is currently under construction, and may not be used. */
double exact_disk(struct vector3 *loc,struct vector3 *mv,double R,struct subvolume *sv,
                  struct molecule *moving,struct molecule *target)
{
  int relevant_walls = 0;
  struct wall_list *wl;
  
  /* How many things we might intersect with? */
  
  for ( wl = sv->wall_head ; wl!=NULL ; wl = wl->next ) relevant_walls++;
  if ( world->x_fineparts[ sv->urb.x ] - loc->x < R ) relevant_walls++;
  if ( world->y_fineparts[ sv->urb.y ] - loc->y < R ) relevant_walls++;
  if ( world->z_fineparts[ sv->urb.z ] - loc->z < R ) relevant_walls++;
  if ( loc->x - world->x_fineparts[ sv->llf.x ] < R ) relevant_walls++;
  if ( loc->y - world->y_fineparts[ sv->llf.y ] < R ) relevant_walls++;
  if ( loc->z - world->z_fineparts[ sv->llf.z ] < R ) relevant_walls++;
  
  if (relevant_walls==0) return 1.0;
  else
  {
    struct vector3 n,u,v,hit;
    double a,b,c,d,f,g,h;
    struct vector2 startpoints[relevant_walls];
    struct vector2 endpoints[ relevant_walls ];
    struct vector2 fix_loc;
    int max_intersecting = 0;
    int idx = 0;
    int direct_hit = 0;
    int added_line;
    struct rxn *rx;
    struct wall *w;
    
    /* Set up an orthogonal coordinate system */
    /* first axis = direction of motion */
    /* second axis = direction from collision location to target particle */
    /* third axis = cross product of first two */
    
    a = 1.0/sqrt(mv->x*mv->x + mv->y*mv->y + mv->z*mv->z);
    n.x = a * mv->x;
    n.y = a * mv->y;
    n.z = a * mv->z;
    
    u.x = target->pos.x - loc->x;
    u.y = target->pos.y - loc->y;
    u.z = target->pos.z - loc->z;
    a = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
    if (a < 0.5*EPS_C)
    {
      direct_hit = 1;
      if (n.x * n.x < EPS_C) { u.x = 0; u.y = -n.z; u.z = n.y; }
      else if (n.y*n.y<EPS_C) { u.x = -n.z; u.y = 0; u.z = n.x; }
      else { u.x = -n.y; u.y = n.x; u.z = 0; }
      a = 1.0/sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
    }
    else a = 1.0/a;
    u.x *= a;
    u.y *= a;
    u.z *= a;
    
    v.x = n.y*u.z - n.z*u.y;
    v.y = n.z*u.x - n.x*u.z;
    v.z = n.x*u.y - n.y*u.x;
    
    d = loc->x*mv->x + loc->y*mv->y + loc->z*mv->z;
    
    /* First we check all the partitions */
    /* We touch partition z if R*u.z*cos(t) + R*v.z*sin(t) == distance to partition */
    /* Solve quadratically--define in a macro to save space */
    
    /*
    Macro does the equivalent of the following (except for setting the once-used b and c)
	  a = u.z*u.z + v.z*v.z;
	  b = f * v.z;              This is actually -b in the QF            
	  c = f*f - u.z*u.z;
	  g = b/a;                  Leading term for R*sin(t)
	  h = (g*g - c/a);          Square of plus/minus term                   
    */  
    #define SOLVE_QF( uw , vw ) a = (uw)*(uw) + (vw)*(vw); g = f*(vw)/a; h = g*g - (f*f - (uw)*(uw))
	
    /*
    And we want to add the intersection points of the  partition neatly, so we macro
    the equivalent of
	  f = world->z_fineparts[ sv->urb.z ] - loc->z;
	  if (f<R)         Might touch
	  {
	    SOLVE_QF( u.z , v.z );
	    if (h>0)       Root exists
	    {
	      h = sqrt(h);
	      if (g-h > -1.0 && g+h < 1.0)     In range
	      {
		startpoints[ idx ].v = g-h;
		startpoints[ idx ].u = sqrt( 1.0 - (g-h)*(g-h) );
		endpoints[ idx ].v = g+h;
		endpoints[ idx ].u = sqrt( 1.0 - (g+h)*(g+h) );
		idx++;
	      }
	    }
	  }
    */
    #define REGISTER_PARTITION startpoints[idx].v=g-h; startpoints[idx].u=sqrt(1-(g-h)*(g-h)); endpoints[idx].v=g+h; endpoints[idx].u=sqrt(1-(g+h)*(g+h)); idx++
    #define CHECK_PARTITION( uw , vw , fw , Rw ) f = (fw); if (f<(Rw)) { SOLVE_QF((uw),(vw)); if (h>0) { h=sqrt(h); if (g-h>-R&&g+h<R) { REGISTER_PARTITION; } } }

    CHECK_PARTITION( u.x , v.x , world->x_fineparts[ sv->urb.x ] - loc->x , R )
    CHECK_PARTITION( u.x , v.x , world->x_fineparts[ sv->llf.x ] - loc->x , -R)
    CHECK_PARTITION( u.y , v.y , world->y_fineparts[ sv->urb.y ] - loc->y , R )
    CHECK_PARTITION( u.y , v.y , world->y_fineparts[ sv->llf.y ] - loc->y , -R)
    CHECK_PARTITION( u.z , v.z , world->z_fineparts[ sv->urb.z ] - loc->z , R )
    CHECK_PARTITION( u.z , v.z , world->z_fineparts[ sv->llf.z ] - loc->z , -R)
    
    #undef CHECK_PARTITION
    #undef REGISTER_PARTITON
    #undef SOLVE_QF
    
    /* Blech, okay, I'm glad that's done with. */
    
    /* Now we check all the walls and return -1 if we can't reach our target */
    /* Ignoring very rare cases of parallel wall and vertices stuck in */
    /* the interaction disk--that makes the code much bulkier */
    
    if (sv->wall_head != NULL)
    {
      fix_loc.u = u.x*target->pos.x + u.y*target->pos.y + u.z*target->pos.z;
      fix_loc.v = v.x*target->pos.x + v.y*target->pos.y + v.z*target->pos.z;
    }
    
    for (wl = sv->wall_head ; wl!=NULL ; wl = wl->next)
    {
      w = wl->this_wall;
      added_line = 0;
      
      if ( (moving->properties->flags & CAN_MOLWALL) != 0 )
      {
	a = w->normal.x * loc->x + w->normal.y * loc->y + w->normal.z * loc->z;
	rx = trigger_intersect(moving->properties->hashval,(struct abstract_molecule*)moving,(a>w->d)?1:-1,w);
	if (rx != NULL && (rx->n_pathways==RX_WINDOW || rx->n_pathways==RX_GHOST) )
	{
	  continue; /* Ghost/window, we can see through it! */
	}
      }
      
      /* Point normal locations */
      a = w->vert[0]->x * n.x + w->vert[0]->y * n.y + w->vert[0]->z * n.z;
      b = w->vert[1]->x * n.x + w->vert[1]->y * n.y + w->vert[1]->z * n.z;
      c = w->vert[2]->x * n.x + w->vert[2]->y * n.y + w->vert[2]->z * n.z;
      
      /* Intersection time of edges */
      if ( a != b ) f = (d - a)/(b - a); else f = -1;
      if ( b != c ) g = (d - b)/(c - b); else g = -1;
      if ( c != a ) h = (d - c)/(a - c); else h = -1;
      
      if (f>0 && f<1)
      {
	if (g>0 && g<1)
	{
	  hit.x = (1-f) * w->vert[0]->x + f * w->vert[1]->x;
	  hit.y = (1-f) * w->vert[0]->y + f * w->vert[1]->y;
	  hit.z = (1-f) * w->vert[0]->z + f * w->vert[1]->z;
	  startpoints[idx].u = hit.x*u.x + hit.y*u.y + hit.z*u.z;
	  startpoints[idx].v = hit.x*v.x + hit.y*v.y + hit.z*v.z;
	  
	  hit.x = (1-g) * w->vert[0]->x + g * w->vert[1]->x;
	  hit.y = (1-g) * w->vert[0]->y + g * w->vert[1]->y;
	  hit.z = (1-g) * w->vert[0]->z + g * w->vert[1]->z;
	  endpoints[idx].u = hit.x*u.x + hit.y*u.y + hit.z*u.z;
	  endpoints[idx].v = hit.x*v.x + hit.y*v.y + hit.z*v.z;
	  
	  added_line = 1;
	}
	else if (h>0 && h<1)
	{
	  hit.x = (1-f) * w->vert[0]->x + f * w->vert[1]->x;
	  hit.y = (1-f) * w->vert[0]->y + f * w->vert[1]->y;
	  hit.z = (1-f) * w->vert[0]->z + f * w->vert[1]->z;
	  startpoints[idx].u = hit.x*u.x + hit.y*u.y + hit.z*u.z;
	  startpoints[idx].v = hit.x*v.x + hit.y*v.y + hit.z*v.z;
	  
	  hit.x = (1-h) * w->vert[0]->x + h * w->vert[1]->x;
	  hit.y = (1-h) * w->vert[0]->y + h * w->vert[1]->y;
	  hit.z = (1-h) * w->vert[0]->z + h * w->vert[1]->z;
	  endpoints[idx].u = hit.x*u.x + hit.y*u.y + hit.z*u.z;
	  endpoints[idx].v = hit.x*v.x + hit.y*v.y + hit.z*v.z;
	  
	  added_line = 1;
	}
      }
      else if (g>0 && g<1 && h>0 && h<1)
      {
	hit.x = (1-g) * w->vert[0]->x + g * w->vert[1]->x;
	hit.y = (1-g) * w->vert[0]->y + g * w->vert[1]->y;
	hit.z = (1-g) * w->vert[0]->z + g * w->vert[1]->z;
	startpoints[idx].u = hit.x*u.x + hit.y*u.y + hit.z*u.z;
	startpoints[idx].v = hit.x*v.x + hit.y*v.y + hit.z*v.z;
	
	hit.x = (1-h) * w->vert[0]->x + h * w->vert[1]->x;
	hit.y = (1-h) * w->vert[0]->y + h * w->vert[1]->y;
	hit.z = (1-h) * w->vert[0]->z + h * w->vert[1]->z;
	endpoints[idx].u = hit.x*u.x + hit.y*u.y + hit.z*u.z;
	endpoints[idx].v = hit.x*v.x + hit.y*v.y + hit.z*v.z;
	
	added_line = 1;
      }
	
      if (added_line && !direct_hit) /* Check if this blocked us from target */
      {
	if (fix_loc.u == 0)
	{
	  if ( startpoints[idx].u != endpoints[idx].u )
	  {
	    a = -startpoints[idx].u / (endpoints[idx].u - startpoints[idx].u);
	    if (a>0 && a<1)
	    {
	      b = ((1-a)*startpoints[idx].v + a*endpoints[idx].v) / fix_loc.v;
	      if (b>0 && b<1) return -1;
	    }
	  }
	}
	else if (fix_loc.v == 0)
	{
	  if ( startpoints[idx].v != endpoints[idx].v )
	  {
	    a = -startpoints[idx].v / (endpoints[idx].v - startpoints[idx].v);
	    if (a>0 && a<1)
	    {
	      b = ((1-a)*startpoints[idx].u + a*endpoints[idx].u) / fix_loc.u;
	      if (b>0 && b<1) return -1;
	    }
	  }
	}
	else
	{
	  c = endpoints[idx].v - startpoints[idx].v + startpoints[idx].u - endpoints[idx].u;
	  if (c != 0)
	  {
	    a = (fix_loc.v * startpoints[idx].u - fix_loc.u * startpoints[idx].v)/c;
	    if (a>0 && a<1)
	    {
	      b = ((1-b)*startpoints[idx].v + a*endpoints[idx].v) / fix_loc.v;
	      if (b>0 && b<1) return -1;
	    }
	  }
	}
      }
      
      idx += added_line;
    }
    
    /* Now we have to compute the area given all the line segments */
    max_intersecting = idx;
  }
 
  return 1.0;
}
#endif

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
       integration.
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
  
  idx = world->num_directions;
  while (idx >= world->num_directions)
  {
    bits = rng_uint(world->rng);
    idx = bits & world->directions_mask;
  }
  if (bits&0x80000000) u.x = world->d_step[idx]; else u.x = -world->d_step[idx];
  if (bits&0x40000000) u.y = world->d_step[idx+1]; else u.y = -world->d_step[idx+1];
  if (bits&0x20000000) u.z = world->d_step[idx+2]; else u.z = -world->d_step[idx+2];
  
  a = (u.x*mv->x + u.y*mv->y + u.z*mv->z);
  
  if (a*a*d2_mv_i < 0.9)  /* Vectors too closely aligned */
  {
    rpt--;
    continue;
  }
  
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
      if (t>0 && t<1) return -1;  /* This wall blocked us! */
    }
  }

#if !defined (USE_EXPANDED_COLLISION_LIST)
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
#endif

  if (upperU < 0 || upperU > 1.0 ||
      lowerU < 0 || lowerU > 1.0 ||
      upperV < 0 || upperV > 1.0 ||
      lowerV < 0 || lowerV > 1.0)
  {
    printf("Crazy!\n");
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
	in should have already been quicksorted, though, so it shouldn't
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
      subvolume that we start in
      displacement vector from current to new location
      linked list of potential collisions with molecules from the 
                starting subvolume
  Out: To the current collision list of molecules we may react with in our 
       subvolume added molecules from the 3 neighboring subvolumes in the 
       direction of movement. Returns new list of collisions.
	
****************************************************************************/
#if 0
struct collision* expand_collision_list(struct molecule *m, struct subvolume *sv, struct vector3 *v, struct collision *shead1)
{
  struct collision *smash;
  struct molecule *mp;
  struct subvolume *nsv1, *nsv2, *nsv3;
  struct rxn *rx;
  /* distance from the molecule to the closest subvolume wall. */
  double distance = 0;
  struct vector3 new_pos; /* new location of the molecule */
  int i;

  if((v == NULL) || (vect_length(v) == 0)) return shead1;
  new_pos.x = m->pos.x + v->x;
  new_pos.y = m->pos.y + v->y;
  new_pos.z = m->pos.z + v->z;
  

  /* find out in what subvolume direction the molecule is moving. */
  if (v->x < 0.0)
  {
    nsv1 = traverse_subvol(sv, &new_pos,X_NEG);
  }
  else
  {
    nsv1 = traverse_subvol(sv, &new_pos,X_POS);
  }

  if (v->y < 0.0)
  {
    nsv2 = traverse_subvol(sv, &new_pos,Y_NEG);
  }
  else
  {
    nsv2 = traverse_subvol(sv,&new_pos,Y_POS);
  }

  if (v->z < 0.0)
  {
    nsv3 = traverse_subvol(sv,&new_pos,Z_NEG);
  }
  else
  {
    nsv3 = traverse_subvol(sv,&new_pos,Z_POS);
  }
     
    /* find the molecules from the neighboring subvolumes that may react with
       the starting molecule. */
    if((nsv1 != NULL) && (nsv1->mol_head != NULL))
    {
    	for (mp = nsv1->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                if(v->x < 0.0) {
    		     distance = world->x_fineparts[ nsv1->urb.x ] - mp->pos.x;
	        }else if(v->x > 0.0){
    		     distance = mp->pos.x - world->x_fineparts[ nsv1->llf.x ];
                }

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    
    if((nsv2 != NULL) && (nsv2->mol_head != NULL))
    {
    	for (mp = nsv2->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                if(v->y < 0.0) {
    		     distance = world->y_fineparts[ nsv2->urb.y ] - mp->pos.y;
	        }else if(v->y > 0.0){
    		     distance = mp->pos.y - world->y_fineparts[ nsv2->llf.y ];
                }

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance >  world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    

    if((nsv3 != NULL) && (nsv3->mol_head != NULL))
    {
    	for (mp = nsv3->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                if(v->z < 0.0) {
    		     distance = world->z_fineparts[ nsv3->urb.z ] - mp->pos.z;
	        }else if(v->z > 0.0){
    		     distance = mp->pos.z - world->z_fineparts[ nsv3->llf.z ];
                }

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance >  world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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

    return shead1;

}
#endif
struct collision* expand_collision_list(struct molecule *m, struct subvolume *sv, struct vector3 *v, struct collision *shead1)
{
  struct collision *smash;
  struct molecule *mp;
  /* neighbors of the current subvolume "face-to-face" */
  struct subvolume *face1, *face2, *face3, *face4, *face5, *face6;
  /* neighbors of the current subvolume "edge-to-edge" */
  struct subvolume *edge1, *edge2, *edge3, *edge4, *edge5, *edge6;
  struct subvolume *edge7, *edge8, *edge9, *edge10, *edge11, *edge12;
  /* neighbors of the current subvolume "corner-to-corner" */
  /*struct subvolume *corner1, *corner2, *corner3, *corner4, *corner5, *corner6;
  struct subvolume *corner7, *corner8;*/
  struct rxn *rx;
  /* distance from the molecule to the closest subvolume wall. */
  double distance = 0;
  /* index of the current subvolume */
  int cur_index;
  /* index of the neighbor subvolume */
  int neighbor_index;
  struct vector3 v0, v1; /* verices of the edge */
  int i;

  if((v == NULL) || (vect_length(v) == 0)) return shead1;
  cur_index = sv->index;

    /*
    nsv1 = traverse_subvol(sv,&(m->pos),X_NEG);
    nsv2 = traverse_subvol(sv,&(m->pos),X_POS);

    nsv3 = traverse_subvol(sv,&(m->pos),Y_NEG);
    nsv4 = traverse_subvol(sv,&(m->pos),Y_POS);

    nsv5 = traverse_subvol(sv,&(m->pos),Z_NEG);
    nsv6 = traverse_subvol(sv,&(m->pos),Z_POS);
     */
    face1 = (struct subvolume *)(sv->neighbor[X_NEG]);
    face2 = (struct subvolume *)(sv->neighbor[X_POS]);
    face3 = (struct subvolume *)(sv->neighbor[Y_NEG]);
    face4 = (struct subvolume *)(sv->neighbor[Y_POS]);
    face5 = (struct subvolume *)(sv->neighbor[Z_NEG]);
    face6 = (struct subvolume *)(sv->neighbor[Z_POS]);

    /* find the molecules from the neighboring subvolumes that may react with
       the starting molecule. */
    if((face1 != NULL) && (face1->mol_head != NULL))
    {
    	for (mp = face1->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
    		distance = world->x_fineparts[ face1->urb.x ] - mp->pos.x;

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    
    if((face2 != NULL) && (face2->mol_head != NULL))
    {
    	for (mp = face2->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
    		distance = mp->pos.x - world->x_fineparts[ face2->llf.x ];

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    
    if((face3 != NULL) && (face3->mol_head != NULL))
    {
    	for (mp = face3->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
    		distance = world->y_fineparts[ face3->urb.y ] - mp->pos.y;

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance >  world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    
    if((face4 != NULL) && (face4->mol_head != NULL))
    {
    	for (mp = face4->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
    		distance = mp->pos.y - world->y_fineparts[ face4->llf.y ];

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance >  world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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

    if((face5 != NULL) && (face5->mol_head != NULL))
    {
    	for (mp = face5->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
    		distance = world->z_fineparts[ face5->urb.z ] - mp->pos.z;

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance >  world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((face6 != NULL) && (face6->mol_head != NULL))
    {
    	for (mp = face6->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume wall than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
    	        distance = mp->pos.z - world->z_fineparts[ face6->llf.z ];

      		if(distance < 0){
			distance = -distance;
      		}      

      		if(distance >  world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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

   /* let's reach subvolumes located "edge-to-edge" from the top face
      of the current subvolume. */
   /* go (-X and +Z) */
   neighbor_index = cur_index -(world->nz_parts-1)*(world->ny_parts-1) + 1;
   edge1 = &(world->subvol[neighbor_index]);
   /* go (+Y and +Z) */
   neighbor_index = cur_index  + (world->nz_parts-1) + 1;
   edge2 = &(world->subvol[neighbor_index]);
   /* go (+X and +Z) */
   neighbor_index = cur_index  + (world->nz_parts-1)*(world->ny_parts-1) + 1;
   edge3 = &(world->subvol[neighbor_index]);
   /* go (-Y and +Z) */
   neighbor_index = cur_index  - (world->nz_parts-1) + 1;
   edge4 = &(world->subvol[neighbor_index]);
   /* let's reach subvolumes located "edge-to-edge" from the bottom face
      of the current subvolume. */
   /* go (-X and -Z) */
   neighbor_index = cur_index - (world->nz_parts-1)*(world->ny_parts-1) - 1;
   edge5 = &(world->subvol[neighbor_index]);
   /* go (+Y and -Z) */
   neighbor_index = cur_index  + (world->nz_parts-1) - 1;
   edge6 = &(world->subvol[neighbor_index]);
   /* go (+X and -Z) */
   neighbor_index = cur_index  + (world->nz_parts-1)*(world->ny_parts-1) - 1;
   edge7 = &(world->subvol[neighbor_index]);
   /* go (-Y and -Z) */
   neighbor_index = cur_index  - (world->nz_parts-1) - 1;
   edge8 = &(world->subvol[neighbor_index]);
   /* let's reach subvolumes located "edge-to-edge" from the verical edges
      of the current subvolume. */
   /* go (-Y and -X) */
   neighbor_index = cur_index  - (world->nz_parts-1) - (world->nz_parts-1)*(world->ny_parts-1);;
   edge9 = &(world->subvol[neighbor_index]);
   /* go (+Y and -X) */
   neighbor_index = cur_index  + (world->nz_parts-1) - (world->nz_parts-1)*(world->ny_parts-1);;
   edge10 = &(world->subvol[neighbor_index]);
   /* go (+Y and +X) */
   neighbor_index = cur_index  + (world->nz_parts-1) + (world->nz_parts-1)*(world->ny_parts-1);;
   edge11 = &(world->subvol[neighbor_index]);
   /* go (-Y and +X) */
   neighbor_index = cur_index  - (world->nz_parts-1) + (world->nz_parts-1)*(world->ny_parts-1);;
   edge12 = &(world->subvol[neighbor_index]);

    if((edge1 != NULL) && (edge1->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->urb.z];
        v1.x = world->x_fineparts[sv->llf.x];
        v1.y = world->y_fineparts[sv->urb.y];
        v1.z = world->z_fineparts[sv->urb.z];
        
    	for (mp = edge1->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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

    if((edge2 != NULL) && (edge2->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->urb.y];
        v0.z = world->z_fineparts[sv->urb.z];
        v1.x = world->x_fineparts[sv->urb.x];
        v1.y = world->y_fineparts[sv->urb.y];
        v1.z = world->z_fineparts[sv->urb.z];
        
    	for (mp = edge2->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge3 != NULL) && (edge3->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->urb.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->urb.z];
        v1.x = world->x_fineparts[sv->urb.x];
        v1.y = world->y_fineparts[sv->urb.y];
        v1.z = world->z_fineparts[sv->urb.z];
        
    	for (mp = edge3->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge4 != NULL) && (edge4->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->urb.z];
        v1.x = world->x_fineparts[sv->urb.x];
        v1.y = world->y_fineparts[sv->llf.y];
        v1.z = world->z_fineparts[sv->urb.z];
        
    	for (mp = edge4->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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

    if((edge5 != NULL) && (edge5->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->llf.z];
        v1.x = world->x_fineparts[sv->llf.x];
        v1.y = world->y_fineparts[sv->urb.y];
        v1.z = world->z_fineparts[sv->llf.z];
        
    	for (mp = edge5->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge6 != NULL) && (edge6->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->urb.y];
        v0.z = world->z_fineparts[sv->llf.z];
        v1.x = world->x_fineparts[sv->urb.x];
        v1.y = world->y_fineparts[sv->urb.y];
        v1.z = world->z_fineparts[sv->llf.z];
        
    	for (mp = edge6->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge7 != NULL) && (edge7->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->urb.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->llf.z];
        v1.x = world->x_fineparts[sv->urb.x];
        v1.y = world->y_fineparts[sv->urb.y];
        v1.z = world->z_fineparts[sv->llf.z];
        
    	for (mp = edge7->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge8 != NULL) && (edge8->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->llf.z];
        v1.x = world->x_fineparts[sv->urb.x];
        v1.y = world->y_fineparts[sv->llf.y];
        v1.z = world->z_fineparts[sv->llf.z];
        
    	for (mp = edge8->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge9 != NULL) && (edge9->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->llf.z];
        v1.x = world->x_fineparts[sv->llf.x];
        v1.y = world->y_fineparts[sv->llf.y];
        v1.z = world->z_fineparts[sv->urb.z];
        
    	for (mp = edge9->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge10 != NULL) && (edge10->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->urb.y];
        v0.z = world->z_fineparts[sv->llf.z];
        v1.x = world->x_fineparts[sv->llf.x];
        v1.y = world->y_fineparts[sv->urb.y];
        v1.z = world->z_fineparts[sv->urb.z];
        
    	for (mp = edge10->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge11 != NULL) && (edge11->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->urb.x];
        v0.y = world->y_fineparts[sv->urb.y];
        v0.z = world->z_fineparts[sv->llf.z];
        v1.x = world->x_fineparts[sv->urb.x];
        v1.y = world->y_fineparts[sv->urb.y];
        v1.z = world->z_fineparts[sv->urb.z];
        
    	for (mp = edge11->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((edge12 != NULL) && (edge12->mol_head != NULL))
    {
        /* find coordinates of the shared edge */
        v0.x = world->x_fineparts[sv->urb.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->llf.z];
        v1.x = world->x_fineparts[sv->urb.x];
        v1.y = world->y_fineparts[sv->llf.y];
        v1.z = world->z_fineparts[sv->urb.z];
        
    	for (mp = edge12->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                distance = distance_point_line(&(mp->pos),&v0,&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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

#if 0
   /* let's reach subvolumes located "corner-to-corner" from the top face
      of the current subvolume. */
   /* go (-X and -Y and +Z) */
   neighbor_index = cur_index - (world->nz_parts-1)*(world->ny_parts-1) - (world->nz_parts-1) + 1;
   corner1 = &(world->subvol[neighbor_index]);
   /* go (-X and +Y and +Z) */
   neighbor_index = cur_index - (world->nz_parts-1)*(world->ny_parts-1) + (world->nz_parts-1) + 1;
   corner2 = &(world->subvol[neighbor_index]);
  /* go (+X and +Y and +Z) */
   neighbor_index = cur_index + (world->nz_parts-1)*(world->ny_parts-1) + (world->nz_parts-1) + 1;
   corner3 = &(world->subvol[neighbor_index]);
  /* go (+X and -Y and +Z) */
   neighbor_index = cur_index + (world->nz_parts-1)*(world->ny_parts-1) - (world->nz_parts-1) + 1;
   corner4 = &(world->subvol[neighbor_index]);

   /* go (-X and -Y and -Z) */
   neighbor_index = cur_index - (world->nz_parts-1)*(world->ny_parts-1) - (world->nz_parts-1) - 1;
   corner5 = &(world->subvol[neighbor_index]);
   /* go (-X and +Y and -Z) */
   neighbor_index = cur_index - (world->nz_parts-1)*(world->ny_parts-1) + (world->nz_parts-1) - 1;
   corner6 = &(world->subvol[neighbor_index]);
  /* go (+X and +Y and -Z) */
   neighbor_index = cur_index + (world->nz_parts-1)*(world->ny_parts-1) + (world->nz_parts-1) - 1;
   corner7 = &(world->subvol[neighbor_index]);
  /* go (+X and -Y and -Z) */
   neighbor_index = cur_index + (world->nz_parts-1)*(world->ny_parts-1) - (world->nz_parts-1) - 1;
   corner8 = &(world->subvol[neighbor_index]);

    if((corner1 != NULL) && (corner1->mol_head != NULL))
    {
        /* find coordinates of the shared corner */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->urb.z];
               
 
    	for (mp = corner1->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                vectorize(&(mp->pos),&v0,&v1);
                distance = vect_length(&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((corner2 != NULL) && (corner2->mol_head != NULL))
    {
        /* find coordinates of the shared corner */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->urb.y];
        v0.z = world->z_fineparts[sv->urb.z];
               
 
    	for (mp = corner2->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                vectorize(&(mp->pos),&v0,&v1);
                distance = vect_length(&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((corner3 != NULL) && (corner3->mol_head != NULL))
    {
        /* find coordinates of the shared corner */
        v0.x = world->x_fineparts[sv->urb.x];
        v0.y = world->y_fineparts[sv->urb.y];
        v0.z = world->z_fineparts[sv->urb.z];
               
 
    	for (mp = corner3->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                vectorize(&(mp->pos),&v0,&v1);
                distance = vect_length(&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((corner4 != NULL) && (corner4->mol_head != NULL))
    {
        /* find coordinates of the shared corner */
        v0.x = world->x_fineparts[sv->urb.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->urb.z];
               
 
    	for (mp = corner4->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                vectorize(&(mp->pos),&v0,&v1);
                distance = vect_length(&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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

    if((corner5 != NULL) && (corner5->mol_head != NULL))
    {
        /* find coordinates of the shared corner */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->llf.z];
               
 
    	for (mp = corner5->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                vectorize(&(mp->pos),&v0,&v1);
                distance = vect_length(&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((corner6 != NULL) && (corner6->mol_head != NULL))
    {
        /* find coordinates of the shared corner */
        v0.x = world->x_fineparts[sv->llf.x];
        v0.y = world->y_fineparts[sv->urb.y];
        v0.z = world->z_fineparts[sv->llf.z];
               
 
    	for (mp = corner6->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                vectorize(&(mp->pos),&v0,&v1);
                distance = vect_length(&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((corner7 != NULL) && (corner7->mol_head != NULL))
    {
        /* find coordinates of the shared corner */
        v0.x = world->x_fineparts[sv->urb.x];
        v0.y = world->y_fineparts[sv->urb.y];
        v0.z = world->z_fineparts[sv->llf.z];
               
 
    	for (mp = corner7->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                vectorize(&(mp->pos),&v0,&v1);
                distance = vect_length(&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
    if((corner8 != NULL) && (corner8->mol_head != NULL))
    {
        /* find coordinates of the shared corner */
        v0.x = world->x_fineparts[sv->urb.x];
        v0.y = world->y_fineparts[sv->llf.y];
        v0.z = world->z_fineparts[sv->llf.z];
               
 
    	for (mp = corner8->mol_head ; mp != NULL ; mp = mp->next_v)
    	{
      		if (mp->properties == NULL) continue;  
      
      		/* if the molecule is farther from the subvolume edge than 
         (world->rx_radius_3d) do not account for the interaction with it.	
      		*/
                vectorize(&(mp->pos),&v0,&v1);
                distance = vect_length(&v1); 

      		if(distance > world->rx_radius_3d){
			continue;
      		}

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
	  			fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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

#endif

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
  double factor;           /* return value from 'estimate_disk()' function */
  double scaling;     /* scales reaction cumulative_probabilitities array */
  double rate_factor=1.0;

  int i,j,k,l;
  
  int calculate_displacement = 1;

  sm = m->properties;
  if (sm==NULL) {
	fprintf(world->err_file,"BROKEN!!!!!\n");
	return NULL;
  }
  if (sm->space_step <= 0.0)
  {
    m->t += max_time;
    return m;
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
	  fprintf(world->err_file,"Out of memory.  Trying to save intermediate states.\n");
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
 
#ifdef USE_EXPANDED_COLLISION_LIST
    shead = expand_collision_list(m, sv, &displacement, shead); 
#endif
 
  if (calculate_displacement)
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


    world->diffusion_number += 1.0;
    world->diffusion_cumsteps += steps;
  }
  
   reflectee = NULL;
 
#define CLEAN_AND_RETURN(x) if (shead2!=NULL) mem_put_list(sv->local_storage->coll,shead2); if (shead!=NULL) mem_put_list(sv->local_storage->coll,shead); return (x)
#define ERROR_AND_QUIT fprintf(world->err_file,"Out of memory: trying to save intermediate results.\n"); i=emergency_output(); fprintf(world->err_file,"Fatal error: out of memory during diffusion of a %s molecule\nAttempt to write intermediate results had %d errors\n",sm->sym->name,i); exit(EXIT_FAILURE)
  do
  {
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
        smash = NULL;
        break;
      }
      
      rx = smash->intermediate;

      if ( (smash->what & COLLIDE_MOL) != 0 && !inert )
      {
	if (smash->t < EPS_C) continue;
	
        m->collisions++;
	world->mol_mol_colls++;

        am = (struct abstract_molecule*)smash->target;
        if ((am->flags & ACT_INERT) != 0)  /* FIXME */
        {
          /*if (smash->t < am->t + am->t2) continue;*/
        }
        
	factor = estimate_disk(
          &(smash->loc),&displacement,world->rx_radius_3d,m->subvol,m,
	  (struct molecule*)am
        );
	if (factor<0) continue; /* Reaction blocked by a wall */
	
#ifdef USE_EXPANDED_COLLISION_LIST
        scaling = 1.0; 
#else
        scaling = factor / rate_factor;
#endif
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
	
	if ( w->effectors != NULL && (sm->flags&CAN_MOLGRID) != 0 )
	{
	  j = xyz2grid( &(smash->loc) , w->effectors );
	  if (w->effectors->mol[j] != NULL)
	  {
	    if (m->index != j || m->previous_grid != w->effectors)
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
		      update_collision_count(sm,w->regions,k,1,rate_factor * w->effectors->binding_factor);
		    }
		    if ((m->flags&COUNT_ME)!=0)
		    {
		      m->flags-=COUNT_ME;
		      count_me_by_region((struct abstract_molecule*)m,-1,NULL);
		    }
		    
		    continue; /* pass through */
		  }
		  else if (l==RX_DESTROY)
		  {
		    if ( (sm->flags & w->flags & COUNT_HITS) )
		      update_collision_count(sm,w->regions,k,0,rate_factor);
		    
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
	      if ( (sm->flags & COUNT_HITS) )
	      {
		update_collision_count(sm,w->regions,k,1,rate_factor);
	      }
	      if ((m->flags&COUNT_ME)!=0)
	      {
		m->flags-=COUNT_ME;
		count_me_by_region((struct abstract_molecule*)m,-1,NULL);
	      }

	      continue; /* Ignore this wall and keep going */
	    }
	    if (rx->prob_t != NULL) check_probs(rx,m->t);
	    i = test_intersect(rx,rate_factor);
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
		  update_collision_count(sm,w->regions,k,1,rate_factor);
		}
		if ((m->flags&COUNT_ME)!=0)
		{
		  m->flags-=COUNT_ME;
		  count_me_by_region((struct abstract_molecule*)m,-1,NULL);
		}

		continue; /* pass through */
	      }
	      else if (j==RX_DESTROY)
	      {
		if ( (sm->flags & COUNT_HITS) )
		  update_collision_count(sm,w->regions,k,0,rate_factor);

		CLEAN_AND_RETURN(NULL);
	      }
	    }
	  }
	}
	
        /* default is to reflect */
        
	if ( (sm->flags & COUNT_HITS) ) update_collision_count(sm,w->regions,k,0,rate_factor);
	if (m->flags&COUNT_ME)
	{
	  m->flags -= COUNT_ME;
	  count_me_by_region((struct abstract_molecule*)m,-1,NULL);
	}

	m->pos.x = smash->loc.x;
	m->pos.y = smash->loc.y;
	m->pos.z = smash->loc.z;
        m->t += t_steps*smash->t;
	reflectee = w;
/*
        m->path_length += t_steps * smash->t *
	                  sqrt(( displacement.x * displacement.x +
                                 displacement.y * displacement.y +
                                 displacement.z * displacement.z ) );
*/
        t_steps *= (1.0-smash->t);
        
        factor = -2.0 * (displacement.x*w->normal.x + displacement.y*w->normal.y + displacement.z*w->normal.z);
        displacement.x = (displacement.x + factor*w->normal.x) * (1.0-smash->t);
        displacement.y = (displacement.y + factor*w->normal.y) * (1.0-smash->t);
        displacement.z = (displacement.z + factor*w->normal.z) * (1.0-smash->t);
        
        break;
      }
      else if ((smash->what & COLLIDE_SUBVOL) != 0)
      {
        struct subvolume *nsv;
/*        
        m->path_length += sqrt( (m->pos.x - smash->loc.x)*(m->pos.x - smash->loc.x)
                               +(m->pos.y - smash->loc.y)*(m->pos.y - smash->loc.y)
                               +(m->pos.z - smash->loc.z)*(m->pos.z - smash->loc.z));
*/
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
	    count_me_by_region((struct abstract_molecule*)m,-1,NULL);
          sm->population--;
          m->properties = NULL;
	  
	  CLEAN_AND_RETURN(NULL);
        }
        else m = migrate_molecule(m,nsv);

        if (shead2 != NULL) mem_put_list(sv->local_storage->coll,shead2);
        if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);
        
        calculate_displacement = 0;
        if (m->properties==NULL) fprintf(world->err_file,"This molecule shouldn't be jumping.\n");
        goto pretend_to_call_diffuse_3D;  /* Jump to beginning of function */        
      }
    }
    
    if (shead2 != NULL) mem_put_list(sv->local_storage->coll,shead2);

  }
  while (smash != NULL);
#undef ERROR_AND_QUIT
#undef CLEAN_AND_RETURN
  
  if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);

  m->pos.x += displacement.x;
  m->pos.y += displacement.y;
  m->pos.z += displacement.z;
  m->t += t_steps;
/*
  m->path_length += sqrt( displacement.x*displacement.x +
                          displacement.y*displacement.y +
                          displacement.z*displacement.z );
*/
  m->index = -1;
  if ((sm->flags&COUNT_ENCLOSED)!=0 && (m->flags&COUNT_ME)==0) count_me_by_region((struct abstract_molecule*)m,1,NULL);
  
  return m;
}
    


/*************************************************************************
run_timestep:
  In: local storage area to use
      time of the next release event
      time of the next checkpoint
  Out: No return value.  Every molecule in the subvolume is updated in
       position and rescheduled at least one timestep ahead.  The
       current_time of the subvolume is also incremented.
*************************************************************************/

void run_timestep(struct storage *local,double release_time,double checkpt_time)
{
  struct abstract_molecule *a;
  struct rxn *r;
  double t;
  double max_time;
  int i,j,err;
  
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
      
      continue;
    }

    a->flags -= IN_SCHEDULE;
    
    if (a->t2 < EPS_C || a->t2 < EPS_C*a->t)
    {
      if ((a->flags & (ACT_INERT+ACT_NEWBIE+ACT_CHANGE)) != 0)
      {
        a->flags -= (a->flags & (ACT_INERT + ACT_NEWBIE + ACT_CHANGE));
        if ((a->flags & ACT_REACT) != 0)
        {
          r = trigger_unimolecular(a->properties->hashval,a);
          if (r->prob_t != NULL) check_probs(r,(a->t + a->t2)*(1.0+EPS_C));
          a->t2 = timeof_unimolecular(r);
	  if (r->prob_t != NULL)
	  {
	    if (a->t + a->t2 > r->prob_t->time)
	    {
	      a->t2 = r->prob_t->time - a->t;
	      a->flags |= ACT_CHANGE;
	    }
	  }
        }
      }
      else if ((a->flags & ACT_REACT) != 0)
      {
        r = trigger_unimolecular(a->properties->hashval,a);
        i = which_unimolecular(r);
        j = outcome_unimolecular(r,i,a,a->t);
	if (j==RX_NO_MEM)
	{
	  fprintf(world->err_file,"Out of memory.  Trying to save intermediate results.\n");
	  i = emergency_output();
	  fprintf(world->err_file,"Out of memory during unimolecular reaction %s...\n",r->sym->name);
	  fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",i);
	  exit( EXIT_FAILURE );
	}
        if (j!=RX_DESTROY) /* We still exist */
        {
          a->t2 = timeof_unimolecular(r);
          if ( r->prob_t != NULL )
          {
            if (a->t + a->t2 > r->prob_t->time)
	    {
              a->t2 = r->prob_t->time - a->t;
              a->flags |= ACT_CHANGE;
	    }
          }
        }
        else /* We don't exist.  Try to recover memory. */
	{
	  if (!(a->flags&IN_VOLUME)) mem_put(a->birthplace,a);
	  continue;
	}
      }
    }

    if ((a->flags & ACT_DIFFUSE) != 0)
    {
      max_time = checkpt_time - a->t;
      if (local->max_timestep < max_time) max_time = local->max_timestep;
      if (a->t2<max_time && (a->flags&(ACT_REACT|ACT_INERT))!=0) max_time = a->t2;

      if ((a->flags & TYPE_3D) != 0)
      {
        if (max_time > release_time - a->t) max_time = release_time - a->t;
        t = a->t;
        a = (struct abstract_molecule*)diffuse_3D((struct molecule*)a , max_time , a->flags & ACT_INERT);
        if (a!=NULL) /* We still exist */
        {
          a->t2 -= a->t - t;
        }
        else continue;
      }
      else  /* Surface diffusion not yet implemented */
      {
        a->t += max_time;
        a->t2 -= max_time;
      }
    }
    else
    {
      if (a->t2==0)
      {
        if ((checkpt_time - a->t) < MAX_UNI_TIMESKIP)
          a->t = checkpt_time;
        else
          a->t += MAX_UNI_TIMESKIP;
      }
      else if (a->t2 + a->t + EPS_C < checkpt_time)
      {
        a->t += a->t2;
        a->t2 = 0;
      }
      else
      {
        a->t2 -= checkpt_time - a->t;
        a->t = checkpt_time;
      }
    }

    a->flags += IN_SCHEDULE;
    if (a->flags & TYPE_GRID) err = schedule_add(local->timer,a);
    else err = schedule_add(((struct molecule*)a)->subvol->local_storage->timer,a);
    
    if (err)
    {
      fprintf(world->err_file,"Out of memory.  Trying to save intermediate results.\n");
      i = emergency_output();
      fprintf(world->err_file,"Out of memory while scheduling molecule of type %s\n",a->properties->sym->name);
      fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",i);
      exit( EXIT_FAILURE );
    }
  }
  if (local->timer->error)
  {
    fprintf(world->err_file,"Out of memory.  Trying to save intermediate results.\n");
    i = emergency_output();
    fprintf(world->err_file,"Out of memory while retrieving molecules to move.\n");
    fprintf(world->err_file,"%d errors while trying to save intermediate results.\n",i);
    exit( EXIT_FAILURE );
  }
  
  local->current_time += 1.0;
}
