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



#define MULTISTEP_WORTHWHILE 2
#define MULTISTEP_PERCENTILE 0.99
#define MULTISTEP_FRACTION 0.9



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
      x = 2.0*rng_double(world->seed++)-1.0;
      y = 2.0*rng_double(world->seed++)-1.0;
      z = 2.0*rng_double(world->seed++)-1.0;
    } while (x*x + y*y + z*z >= 1.0 || x*x + y*y + z*z < 0.001);
    r = scale * world->r_step[ rng_uint(world->seed++) & (world->radial_subdivisions-1) ] / sqrt(x*x + y*y + z*z);
    v->x = r*x;
    v->y = r*y;
    v->z = r*z;
    return;
  }
#endif
    
  
  if (world->fully_random)
  {
    double h,r_sin_phi,theta;
    
    r = scale * world->r_step[ rng_uint(world->seed++) & (world->radial_subdivisions-1) ];
    h = 2.0*rng_double(world->seed++) - 1.0;
    theta = 2.0*MY_PI*rng_double(world->seed++);
    
    r_sin_phi = r * sqrt(1.0 - h*h);
    
    v->x = r_sin_phi * cos(theta);
    v->y = r_sin_phi * sin(theta);
    v->z = r*h;
  }
  else
  {
    u_int x_bit,y_bit,z_bit;
    uint r_bit,thetaphi_bit;
    uint bits;
    uint idx;
    
    bits = rng_uint(world->seed++);
    
    x_bit =        (bits & 0x80000000);
    y_bit =        (bits & 0x40000000);
    z_bit =        (bits & 0x20000000);
    thetaphi_bit = (bits & 0x1FFFF000) >> 12;
    r_bit =        (bits & 0x00000FFF);
    
    r = scale * world->r_step[ r_bit & (world->radial_subdivisions - 1) ];
    
    idx = (thetaphi_bit & world->directions_mask);
    while ( idx >= world->num_directions)
    {
      idx = ( rng_uint(world->seed++) & world->directions_mask);
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
  Out: collision list of walls and molecules we intersected along our ray
         (current subvolume only)
*************************************************************************/

struct collision* ray_trace(struct molecule *m, struct collision *c,
                            struct subvolume *sv, struct vector3 *v)
{
  struct collision *smash,*shead;
  struct wall_list *wlp;
  struct wall_list fake_wlp;
  double dx,dy,dz,tx,ty,tz,tmin;
  int i,j,k;
  
  shead = NULL;
  smash = (struct collision*) mem_get(sv->mem->coll);
  
  fake_wlp.next = sv->wall_head;
    
  for (wlp = sv->wall_head ; wlp != NULL; wlp = wlp->next)
  {
    i = collide_wall(&(m->pos),v,wlp->this_wall,&(smash->t),&(smash->loc));

    if (i==COLLIDE_REDO)
    {
      if (shead != NULL) mem_put_list(sv->mem->coll,shead);
      shead = NULL;
      wlp = &fake_wlp;
      continue;
    }
    else if (i!=COLLIDE_MISS)
    {
      if (smash->t < tmin) tmin = smash->t;
      smash->what = COLLIDE_WALL + i;
      smash->target = (void*) wlp->this_wall;
      smash->next = shead;
      shead = smash;
      smash = (struct collision*) mem_get(sv->mem->coll);
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
    i = collide_mol(&(m->pos),v,(struct abstract_molecule*)c->target,
                    &(c->t),&(c->loc));
    if (i != COLLIDE_MISS)
    {
      smash = (struct collision*) mem_get(sv->mem->coll);
      memcpy(smash,c,sizeof(struct collision));
      
      smash->what = i + COLLIDE_MOL;

      smash->next = shead;
      shead = smash;
    }
  }
  
  return shead;
}


void tell_loc(struct molecule *m,char *s)
{
  if (0 || s[0] == '\0')
  printf("%sMy name is %x and I live at %.3f,%.3f,%.3f\n",
         s,(int)m,m->pos.x*world->length_unit,m->pos.y*world->length_unit,m->pos.z*world->length_unit);
}


double estimate_disk(struct vector3 *loc,struct vector3 *mv,double R,struct subvolume *sv)
{
  struct vector3 u,v;
  double d2_mv_i,a,b,t;
  double upperU = 1.0;
  double upperV = 1.0;
  double lowerU = 1.0;
  double lowerV = 1.0;
  struct wall_list *wl;
  
  pick_displacement(&u,1.0);
  d2_mv_i = 1.0/(mv->x*mv->x + mv->y*mv->y + mv->z*mv->z);
  a = d2_mv_i*(u.x*mv->x + u.y*mv->y + u.z*mv->z);
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
  
/* FIXME--this is horribly inefficient; should create modified collide_wall */

  for (wl = sv->wall_head ; wl!=NULL ; wl = wl->next)
  {
    t = touch_wall(loc,&u,wl->this_wall);
    if (t>0.0 && t<upperU) upperU = t;
    if (t<0.0 && -t<lowerU) lowerU = -t;
    
    t = touch_wall(loc,&v,wl->this_wall);
    if (t>0.0 && t<upperV) upperV = t;
    if (t<0.0 && -t<lowerV) lowerV = -t;
  }


  if (u.x > EPS_C)
  {
    u.x = 1/u.x;
    t = (world->x_fineparts[sv->urb.x] - loc->x)*u.x;
    if (t < upperU) upperU = t;
    t = (loc->x - world->x_fineparts[sv->llf.x])*u.x;
    if (t < lowerU) lowerU = t;
  }
  else if (u.x < EPS_C)
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
  else if (u.y < EPS_C)
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
  else if (u.z < EPS_C)
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
  else if (v.x < EPS_C)
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
  else if (v.y < EPS_C)
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
  else if (v.z < EPS_C)
  {
    v.z = 1/v.z;
    t = (world->z_fineparts[sv->llf.z] - loc->z)*v.z;
    if (t < upperV) upperV = t;
    t = (loc->z - world->z_fineparts[sv->urb.z])*v.z;
    if (t < lowerV) lowerV = t;
  }

  a = 4.0/(upperU*upperU + lowerU*lowerU + upperV*upperV + lowerV*lowerV);
/*
  if (a > 1.1) printf("Correction factor %.2f\n",a);
  if (a < 1.0-EPS_C) printf("MUDDY BLURDER! a=%.2f R=%.2f u=[%.2f %.2f %.2f] %.2f %.2f %.2f %.2f\n",
                            a,R,u.x,u.y,u.z,upperU,lowerU,upperV,lowerV);
*/
  return a;
}


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


int search_memory_for_me(struct mem_helper *mh,struct abstract_list *al)
{
  int i;
  
  for (i=0;i<mh->buf_len;i++)
  {
    if (al == (struct abstract_list*)(mh->storage + mh->record_size*i)) return 1;
  }
  
  if (mh->next_helper == NULL) return 0;
  else return search_memory_for_me(mh->next_helper,al);
}

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


/* This isn't very efficient for lots of coincident objects; it uses
a bubble-sort-like algorithm (albeit on a nearly sorted list). */

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
#define TOL 10.0*EPS_C
  struct vector3 displacement;
  struct collision *smash,*shead,*shead2;
  struct subvolume *sv;
  struct wall_list *wl;
  struct wall *w;
  struct rxn *rx;
  struct molecule *mp,*old_mp;
  struct grid_molecule *g;
  struct abstract_molecule *am;
  struct species *sm;
  double d2;
  double d2_nearmax;
  double d2min = GIGANTIC;
  double steps=1.0;
  double factor;
  double rate_factor=1.0;
  double x,xx;
  
  int i,j,k,l;
  
  int calculate_displacement = 1;

  sm = m->properties;
  if (sm==NULL) printf("BROKEN!!!!!\n");
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
  
pretend_to_call_diffuse_3D:   /* Label to allow fake recursion */

  sv = m->subvol;
  
  shead = NULL;
  old_mp = NULL;
  if ( (m->properties->flags & CAN_MOLMOL) != 0 )
  {
    for (mp = sv->mol_head ; mp != NULL ; old_mp = mp , mp = mp->next_v)
    {
continue_special_diffuse_3D:   /* Jump here instead of looping if old_mp,mp already set */

      if (mp==m) continue;
      if (m->properties!=sm)
      {
        if (calculate_displacement) printf("WHAAA????\n");
        else printf("HWAAAAA????\n");
      }
      
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
               m->properties->hashval,mp->properties->hashval,
               (struct abstract_molecule*)m,(struct abstract_molecule*)mp,0,0
             );
      
      if (rx != NULL)
      {
        smash = mem_get(sv->mem->coll);
        smash->target = (void*) mp;
        smash->intermediate = rx;
        smash->next = shead;
        smash->what = COLLIDE_MOL;
        shead = smash;
      }
    }
  }
  
  if (calculate_displacement)
  {
    if (max_time > MULTISTEP_WORTHWHILE)
    {
      d2_nearmax = sm->space_step * world->r_step[ (int)(world->radial_subdivisions * MULTISTEP_PERCENTILE) ];
      d2_nearmax *= d2_nearmax;
    
      if ( (m->properties->flags & CAN_MOLMOL) != 0 )
      {
        for (smash = shead ; smash != NULL ; smash = smash->next)
        {
          mp = (struct molecule*)smash->target;
          d2 = (m->pos.x - mp->pos.x)*(m->pos.x - mp->pos.x) +
               (m->pos.y - mp->pos.y)*(m->pos.x - mp->pos.y) +
               (m->pos.z - mp->pos.z)*(m->pos.x - mp->pos.z);
          if (d2 < d2min) d2min = d2;
        }
      }
      for (wl = sv->wall_head ; wl!=NULL ; wl = wl->next)
      {
        w = wl->this_wall;
        d2 = w->normal.x*m->pos.x + w->normal.y*m->pos.y + w->normal.z*m->pos.z;
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
      else if ( d2_nearmax*max_time < d2min ) steps = (1.0+EPS_C)*max_time;
      else steps = d2min / d2_nearmax;
      
      if (steps < MULTISTEP_WORTHWHILE) steps = 1.0;
    }
    else if (max_time < MULTISTEP_FRACTION) steps = max_time;
    else steps = 1.0;
    
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
  
  d2 = displacement.x*displacement.x + displacement.y*displacement.y +
       displacement.z*displacement.z;
       
  do
  {
    shead2 = ray_trace(m,shead,sv,&displacement);
    
    if (shead2!=NULL && shead2->next!=NULL)  /* Could be sped up/combined */
    {
      shead2 = (struct collision*)ae_list_sort((struct abstract_element*)shead2);
      shead2 = gather_walls_first(shead2,TOL);
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
        m->collisions++;

        am = (struct abstract_molecule*)smash->target;
        if ((am->flags & ACT_INERT) != 0)  /* FIXME */
        {
          if (smash->t < am->t + am->t2) continue;
        }

        factor = rate_factor / estimate_disk(
          &(smash->loc),&displacement,world->rx_radius_3d,m->subvol
        );
        
        if (rx->rate_t != NULL) check_rates(rx,m->t);

        i = test_bimolecular(rx,factor);
        if (i<0) continue;
        
        j = outcome_bimolecular(
                rx,i,(struct abstract_molecule*)m,
                am,0,0,m->t+steps*smash->t,&(smash->loc)
              );
        
        if (j) continue;
        else
        {
          if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
          if (shead != NULL) mem_put_list(sv->mem->coll,shead);
          return NULL;
        }
      }

      else if ( (smash->what & COLLIDE_WALL) != 0 )

      {
  /* We may have to handle coincident walls, which is a pain. */
        if (smash->next!=NULL && (smash->next->what&COLLIDE_WALL)!=0 &&
            smash->next->t - smash->t <= TOL)
        {
          struct collision *smish,*g_head,*gp,*w_head,*wp,*c;
          int is_window = 0;
          int is_ghost = 0;
          int count_type = 0;
          
          g_head = gp = w_head = wp = NULL;
          
          for (smish=smash ; smish != NULL ; smish=smish->next)
          {
            w = (struct wall*) smish->target;
            world->ray_polygon_colls++;
            if ( (smish->what & COLLIDE_MASK) == COLLIDE_FRONT ) k = 1;
            else k = -1;
            
            if ( w->effectors != NULL && (sm->flags&CAN_MOLGRID) != 0 && !is_window )
            {
              j = xyz2grid( &(smish->loc) , w->effectors );
              if (w->effectors->mol[j] != NULL)
              {
                if (m->index != j || m->releaser != w->effectors)
                {
                  g = w->effectors->mol[j];
                  rx = trigger_bimolecular(
                    sm->hashval,g->properties->hashval,
                    (struct abstract_molecule*)m,(struct abstract_molecule*)g,
                    k,g->orient
                  );
                  if (rx!=NULL)
                  {
                    if (rx->rate_t != NULL) check_rates(rx,m->t);
                    c = (struct collision*)mem_get(sv->mem->coll);
                    c->intermediate = rx;
                    c->target = g;
                    c->t = rx->cum_rates[ rx->n_pathways - 1 ] * w->effectors->binding_factor * rate_factor;
                    c->next = NULL;
                    if (gp==NULL)
                    {
                      g_head = gp = c;
                      c->loc.x = c->t;
                    }
                    else
                    {
                      gp->next = c;
                      c->loc.x = gp->loc.x + c->t;
                      gp = c;
                    }
                  }
                }
              }
            }
            if ( !is_window )
            {
              m->index = -1;
              
              if ( (sm->flags&CAN_MOLWALL) == 0 ) rx = NULL;
              else rx = trigger_intersect(
                       m->properties->hashval,(struct abstract_molecule*)m,k,w
                     );
                     
              if (rx!=NULL && rx->n_pathways==RX_WINDOW) is_window++;
              else if (rx==NULL || rx->n_pathways!=RX_GHOST)
              {
                c = (struct collision*)mem_get(sv->mem->coll);
                c->intermediate = rx;
                c->target = w;
                if (rx!=NULL)
                {
                  c->t = rx->cum_rates[ rx->n_pathways - 1 ] * rate_factor;
                }
                else c->t = 0.0;
                
                c->next = NULL;
                if (wp==NULL)
                {
                  w_head = wp = c;
                  c->loc.x = c->t;
                }
                else
                {
                  wp->next = c;
                  c->loc.x = wp->loc.x + c->t;
                  wp = c;
                }
              }
            }
            if (smish->next == NULL) break;
            if (smish->next->t - smash->t > TOL) break;
          }
          
          if (is_window)
          {
            for (c=smash;c!=smish->next;c=c->next)
            {
              w = (struct wall*) c->target;
              if ( (c->what & COLLIDE_MASK) == COLLIDE_FRONT ) k = 1;
              else k = -1;
              
              if ( (sm->flags & w->flags & COUNT_HITS) )
                update_collision_count(m->properties,w->regions,k,1);
            }
                
            if (g_head!=NULL) mem_put_list(sv->mem->coll,g_head);
            if (w_head!=NULL) mem_put_list(sv->mem->coll,w_head);
            
            smash=smish;
            continue;
          }
          
          if (gp!=NULL) /* Possibility of reaction with grid molecule */
          {
            x = gp->loc.x;
            xx = rng_double(world->seed++);
            if (x>1.0) xx *= x;
            
            if (x >= xx)
            {
              for (c=g_head;c!=NULL;c=c->next)
              {
                if (c->loc.x >= x)
                {
                  rx = c->intermediate;
                  g = (struct grid_molecule*)c->target;
                  
                  i = test_bimolecular(rx,1/c->t); /* Always works */
                  if (i >= 0)
                  {
                    l = outcome_bimolecular(
                      rx,i,(struct abstract_molecule*)m,
                      (struct abstract_molecule*)g,
                      k,g->orient,m->t+steps*smash->t,&(smash->loc)
                    );
                    
                    if (l==1)
                    {
                      if (w_head!=NULL) mem_put_list(sv->mem->coll,w_head);
                      count_type = 1;  /* pass through */
                      break;
                    }
                    else if (l==0)
                    {
                      count_type = -1;  /* destroyed */
                      break;
                    }
                  }
                }
              }
            }
            mem_put_list(sv->mem->coll,g_head);
          }
          
          if (count_type==0) /* Possibility of a reaction with a wall */
          {
            if (wp->loc.x==0.0) count_type = 0; /* All reflective */
            else
            {
              x = wp->loc.x;
              xx = rng_double(world->seed++);
              if (x>1.0) xx *= x;
              
              if (x >= xx)
              {
                for (c=w_head;c!=NULL;c=c->next)
                {
                  if (c->loc.x >= x)
                  {
                    rx = c->intermediate;
                    if (rx==NULL) continue; /* Was reflective */
                    
                    w = (struct wall*)c->target;
                    
                    i = test_intersect(rx,rate_factor);
                    if (i>=0)
                    {
                      j = outcome_intersect(
                              rx,i,w,(struct abstract_molecule*)m,
                              k,m->t + steps*smash->t,&(smash->loc)
                            );
                      if (j==1)
                      {
                        count_type = 1; /* pass through */
                        break;
                      }
                      else if (j==0)
                      {
                        count_type = -1; /* destroyed */
                        break;
                      }
                    }
                  }
                }
              }
            }
            mem_put_list(sv->mem->coll,w_head);
          }
          
          for (c=smash;c!=smish->next;c=c->next)
          {
            w = (struct wall*) c->target;
            if ( (c->what & COLLIDE_MASK) == COLLIDE_FRONT ) k = 1;
            else k = -1;
            
            if ( (sm->flags & w->flags & COUNT_HITS) )
            {
              if (count_type==1) update_collision_count(sm,w->regions,k,1);
              else update_collision_count(sm,w->regions,k,0);
            }
          }
                
          smash = smish;
          if (count_type==-1)
          {
            if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
            if (shead != NULL) mem_put_list(sv->mem->coll,shead);
            return NULL;
          }
          else if (count_type==1 || is_ghost) continue;
        }
        else /* Simple case, no coincident walls */
        {
          w = (struct wall*) smash->target;
          
          world->ray_polygon_colls++;
          
          if ( (smash->what & COLLIDE_MASK) == COLLIDE_FRONT ) k = 1;
          else k = -1;
          
          if ( (sm->flags & CAN_MOLWALL) != 0 )
          {
            rx = trigger_intersect(sm->hashval,(struct abstract_molecule*)m,k,w);
            if (rx!=NULL && rx->n_pathways == RX_WINDOW)
            {
              if ( (sm->flags & w->flags & COUNT_HITS) )
                update_collision_count(sm,w->regions,k,1);

              continue; /* pass through */
            }
          }
                    
          if ( w->effectors != NULL && (sm->flags&CAN_MOLGRID) != 0 )
          {
            j = xyz2grid( &(smash->loc) , w->effectors );
            if (w->effectors->mol[j] != NULL)
            {
              if (m->index != j || m->releaser != w->effectors)
              {
                g = w->effectors->mol[j];
                rx = trigger_bimolecular(
                  sm->hashval,g->properties->hashval,
                  (struct abstract_molecule*)m,(struct abstract_molecule*)g,
                  k,g->orient
                );
                if (rx!=NULL)
                {
                  if (rx->rate_t != NULL) check_rates(rx,m->t);
                  i = test_bimolecular(rx,rate_factor * w->effectors->binding_factor);
                  if (i >= 0)
                  {
                    l = outcome_bimolecular(
                      rx,i,(struct abstract_molecule*)m,
                      (struct abstract_molecule*)g,
                      k,g->orient,m->t+steps*smash->t,&(smash->loc)
                    );
                    
                    if (l==1)
                    {
                      if ( (m->properties->flags & w->flags & COUNT_HITS) )
                        update_collision_count(m->properties,w->regions,k,1);
                      
                      continue; /* pass through */
                    }
                    else if (l==0)
                    {
                      if ( (sm->flags & w->flags & COUNT_HITS) )
                        update_collision_count(m->properties,w->regions,k,0);
                                        
                      if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
                      if (shead != NULL) mem_put_list(sv->mem->coll,shead);

                      return NULL;
                    }
                  }
                }
              }
              else
              {
  /*              printf("Nope!  I'm not rebinding.\n");  */
                m->index = -1;
              }
            }
          }
          
          if ( (m->properties->flags&CAN_MOLWALL) != 0 )
          {
            m->index = -1;
            rx = trigger_intersect(
                    m->properties->hashval,(struct abstract_molecule*)m,k,w
                  );
            
            if (rx != NULL)
            {
              if (rx->n_pathways == RX_GHOST)
              {
                if ( (sm->flags & w->flags & COUNT_HITS) )
                  update_collision_count(sm,w->regions,k,1);

                continue;
              }
              if (rx->rate_t != NULL) check_rates(rx,m->t);
              i = test_intersect(rx,rate_factor);
              if (i>=0)
              {
                j = outcome_intersect(
                        rx,i,w,(struct abstract_molecule*)m,
                        k,m->t + steps*smash->t,&(smash->loc)
                      );
                if (j==1)
                {
                  if ( (m->properties->flags & w->flags & COUNT_HITS) )
                    update_collision_count(m->properties,w->regions,k,1);

                  continue; /* pass through */
                }
                else if (j==0)
                {
                  if ( (sm->flags & w->flags & COUNT_HITS) )
                    update_collision_count(m->properties,w->regions,k,0);

                  if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
                  if (shead != NULL) mem_put_list(sv->mem->coll,shead);

                  return NULL;
                }
              }
            }
          }
          
          if ( (m->properties->flags & w->flags & COUNT_HITS) )
            update_collision_count(m->properties,w->regions,k,0);
        }
        /* default is to reflect */
        
        m->pos.x += displacement.x * smash->t;
        m->pos.y += displacement.y * smash->t;
        m->pos.z += displacement.z * smash->t;
        tell_loc(m,"Boing!!  ");
        m->t += steps*smash->t;
        m->path_length += sqrt( smash->t*smash->t*
                                ( displacement.x * displacement.x +
                                  displacement.y * displacement.y +
                                  displacement.z * displacement.z ) );
        steps *= (1.0-smash->t);
        
        factor = -2.0 * (displacement.x*w->normal.x + displacement.y*w->normal.y + displacement.z*w->normal.z);
        displacement.x = (displacement.x + factor*w->normal.x) * (1.0-smash->t);
        displacement.y = (displacement.y + factor*w->normal.y) * (1.0-smash->t);
        displacement.z = (displacement.z + factor*w->normal.z) * (1.0-smash->t);
        
        break;
      }
      else if ((smash->what & COLLIDE_SUBVOL) != 0)
      {
        struct subvolume *nsv;
        
        m->path_length += sqrt( (m->pos.x - smash->loc.x)*(m->pos.x - smash->loc.x)
                               +(m->pos.y - smash->loc.y)*(m->pos.y - smash->loc.y)
                               +(m->pos.z - smash->loc.z)*(m->pos.z - smash->loc.z));

        m->pos.x = smash->loc.x;
        m->pos.y = smash->loc.y;
        m->pos.z = smash->loc.z;
        displacement.x *= (1.0-smash->t);
        displacement.y *= (1.0-smash->t);
        displacement.z *= (1.0-smash->t);
        
        m->t += steps*smash->t;
        steps *= (1.0-smash->t);
        if (steps < EPS_C) steps = EPS_C;

        nsv = traverse_subvol(sv,&(m->pos),smash->what - COLLIDE_SV_NX - COLLIDE_SUBVOL);
        if (nsv==NULL)
        {
          fprintf(world->log_file,"Error: a %s molecule escaped the world at (%.2e,%.2e,%.2e)\n",
                  m->properties->sym->name,m->pos.x*world->length_unit,
                  m->pos.y*world->length_unit,m->pos.z*world->length_unit);
          m->properties->population--;
          m->properties = NULL;
          if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
          if (shead != NULL) mem_put_list(sv->mem->coll,shead);
          
          return NULL;
        }
        else m = migrate_molecule(m,nsv);

        if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
        if (shead != NULL) mem_put_list(sv->mem->coll,shead);
        
        calculate_displacement = 0;
        if (m->properties==NULL) printf("This molecule shouldn't be jumping.\n");
        goto pretend_to_call_diffuse_3D;  /* Jump to beginning of function */        
      }
    }
    
    if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
  }
  while (smash != NULL);
  
  if (shead != NULL) mem_put_list(sv->mem->coll,shead);
  m->pos.x += displacement.x;
  m->pos.y += displacement.y;
  m->pos.z += displacement.z;
  tell_loc(m,"Whew...  ");
  m->t += steps;
  m->path_length += sqrt( displacement.x*displacement.x +
                          displacement.y*displacement.y +
                          displacement.z*displacement.z );
  m->index = -1;
  return m;
#undef TOL
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
  double stop_time,max_time;
  int i,j;
  
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
    
    if (local->current_time + local->max_timestep < checkpt_time) stop_time = local->max_timestep;
    else stop_time = checkpt_time - local->current_time;
    
    if (a->t2 < EPS_C || a->t2 < EPS_C*a->t)
    {
      if ((a->flags & (ACT_INERT+ACT_NEWBIE+ACT_CHANGE)) != 0)
      {
        a->flags -= (a->flags & (ACT_INERT + ACT_NEWBIE + ACT_CHANGE));
        if ((a->flags & ACT_REACT) != 0)
        {
          r = trigger_unimolecular(a->properties->hashval,a);
          if (r->rate_t != NULL) check_rates(r,(a->t + a->t2)*(1.0+EPS_C));
          a->t2 = timeof_unimolecular(r);
        }
      }
      else if ((a->flags & ACT_REACT) != 0)
      {
        r = trigger_unimolecular(a->properties->hashval,a);
        i = which_unimolecular(r);
        j = outcome_unimolecular(r,i,a,a->t);
        if (j) /* We still exist */
        {
          a->t2 = timeof_unimolecular(r);
          if ( r->rate_t != NULL )
          {
            if (a->t + a->t2 > r->rate_t->time)
            a->t2 = r->rate_t->time - a->t;
            a->flags |= ACT_CHANGE;
          }
        }
        else continue;
      }
    }

    if ((a->flags & (ACT_INERT+ACT_REACT)) == 0) max_time = stop_time;
    else if (a->t2 < stop_time) max_time = a->t2;
    else max_time = stop_time;
    
    if ((a->flags & ACT_DIFFUSE) != 0)
    {
      if ((a->flags & TYPE_3D) != 0)
      {
        if (max_time > release_time - local->current_time) max_time = release_time - local->current_time;
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
      if (a->t2==0) a->t = checkpt_time;
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
    if (a->flags & TYPE_GRID) schedule_add(local->timer,a);
    else schedule_add(((struct molecule*)a)->subvol->mem->timer,a);
  }
  
  local->current_time += 1.0;
}
