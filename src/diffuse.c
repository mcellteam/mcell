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
#include "vol_util.h"
#include "wall_util.h"
#include "react.h"



#define MULTISTEP_WORTHWHILE 2.0
#define MULTISTEP_PERCENTILE 0.99
#define MULTISTEP_FRACTION 0.98
#define SET_ME_PROPERLY 1.0



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
    
    x_bit =        (bits & 0x8000000);
    y_bit =        (bits & 0x4000000);
    z_bit =        (bits & 0x2000000);
    thetaphi_bit = (bits & 0x1FFFF000) >> 12;
    r_bit =        (bits & 0x0000FFF);
    
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
    i = collide_wall(&(m->pos),v,wlp->this_wall,&(smash->t),&(smash->loc),m->collisions==-74);

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
  if (m->collisions == -74)
  printf("%sMy name is %x and I live at %.3f,%.3f,%.3f\n",
         s,(int)m,m->pos.x*world->length_unit,m->pos.y*world->length_unit,m->pos.z*world->length_unit);
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
  struct vector3 displacement;
  struct collision *smash,*shead,*shead2;
  struct subvolume *sv;
  struct wall *w;
  struct rxn *r;
  struct molecule *mp,*old_mp;
  struct species *sm;
  double d2;
  double d2_nearmax;
  double d2min = GIGANTIC;
  double steps;
  double factor;
  
  int i,j,k;
  
  int calculate_displacement = 1;
  
  int listed = 0;
  
  sm = m->properties;
  if (sm->space_step <= 0.0)
  {
    m->t += max_time;
    return m;
  }
  
pretend_to_call_diffuse_3D:   /* Label to allow fake recursion */

  sv = m->subvol;
  
  shead = NULL;
  old_mp = NULL;
  for (mp = sv->mol_head ; mp != NULL ; old_mp = mp , mp = mp->next_v)
  {
    if (mp==m) listed++;
    
    if (mp->properties == NULL)  /* Reclaim storage */
    {
      if (old_mp != NULL) old_mp->next_v = mp->next_v;
      else sv->mol_head = mp->next_v;
      
      if ((mp->flags & IN_MASK)==IN_VOLUME) mem_put(mp->birthplace,mp);
      else if ((mp->flags & IN_VOLUME) != NULL) mp->flags -= IN_VOLUME;
      
      continue;
    }
    
    r = trigger_bimolecular(
            m->properties->hashval,mp->properties->hashval,
            (struct abstract_molecule*)m,(struct abstract_molecule*)mp,0,0
          );
    
    if (r != NULL)
    {
      smash = mem_get(sv->mem->coll);
      smash->target = (void*) mp;
      smash->intermediate = r;
      smash->next = shead;
      smash->what = COLLIDE_MOL;
      shead = smash;
    }
  }
  
  if (calculate_displacement)
  {
    if (max_time > MULTISTEP_WORTHWHILE)
    {
      d2_nearmax = sm->space_step * world->r_step[ (int)(world->radial_subdivisions * MULTISTEP_PERCENTILE) ];
      d2_nearmax *= d2_nearmax;
    
      for (smash = shead ; smash != NULL ; smash = smash->next)
      {
        mp = (struct molecule*)smash->target;
        d2 = (m->pos.x - mp->pos.x)*(m->pos.x - mp->pos.x) +
             (m->pos.y - mp->pos.y)*(m->pos.x - mp->pos.y) +
             (m->pos.z - mp->pos.z)*(m->pos.x - mp->pos.z);
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
    
    if (steps == 1.0) pick_displacement(&displacement,sm->space_step);
    else pick_displacement(&displacement,sqrt(steps)*sm->space_step);
    
  }
  
  d2 = displacement.x*displacement.x + displacement.y*displacement.y +
       displacement.z*displacement.z;
       
#if 0
  if (sqrt(d2)*world->length_unit > 0.5)
  {
    if (calculate_displacement)
      printf("Yikes, molecule %x wants to travel %.3f on seed %d!\n",
             (int)m,sqrt(d2),world->seed-1);
    else
      printf("Yowzers, molecule %x wants to travel %.3f more!\n",
             (int)m,sqrt(d2));
  }
#endif
  
  do
  {
    shead2 = ray_trace(m,shead,sv,&displacement);
    
    shead2 = (struct collision*)ae_list_sort((struct abstract_element*)shead2);
    
    for (smash = shead2; smash != NULL; smash = smash->next)
    {
      if (smash->t >= 1.0 || smash->t < 0.0)
      {
        smash = NULL;
        break;
      }
      
      if (smash->next != NULL && smash->next->t - smash->t < 10*EPS_C &&
          (smash->what & COLLIDE_SUBVOL)!= 0)
      {
        struct collision *temp;
        temp = smash->next;
        smash->next = temp->next;
        temp->next = smash;
        smash = temp;
      }

      if ( (smash->what & COLLIDE_MOL) != 0 && !inert )
      {
/*        m->collisions++; */
        i = test_bimolecular(smash->intermediate,SET_ME_PROPERLY);
        if (i<0) continue;
        
        j = outcome_bimolecular(
                smash->intermediate,i,(struct abstract_molecule*)m,
                (struct abstract_molecule*)smash->target,
                0,0,m->t+steps*smash->t
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
        w = (struct wall*) smash->target;
        
        if ( (smash->what & COLLIDE_MASK) == COLLIDE_FRONT ) k = 1;
        else k = -1;
        r = trigger_intersect(
                m->properties->hashval,(struct abstract_molecule*)m,k,w
              );
        
        if (r != NULL)
        {
          i = test_intersect(r,SET_ME_PROPERLY);
          if (i < 0)
          {
            tell_loc(m,"(Pass)  ");
            continue; /* pass through--set counters here! */
          }
          else
          {
            j = outcome_intersect(
                    r,i,w,(struct abstract_molecule*)m,k,m->t + steps*smash->t
                  );
            if (j==1) continue; /* pass through -- set counters! */
            else if (j==0)
            {
              if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
              if (shead != NULL) mem_put_list(sv->mem->coll,shead);
              return NULL;
            }
          }
        }
        /* default is to reflect */
        smash->t *= (1.0 - EPS_C);
        
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
        tell_loc(m,"Whoosh!  ");

#if 0
        if ((fabs(m->pos.x)-10.0)>EPS_C ||
            (fabs(m->pos.y)-10.0)>EPS_C ||
            (fabs(m->pos.z)-10.0)>EPS_C)
        {
          struct wall_list *twlp;
          printf("  We shouldn't be whooshing at %.3f!\n",smash->t);
          printf("  We are: %.3f %.3f %.3f -> %.3f %.3f %.3f\n",
                 m->pos.x,m->pos.y,m->pos.z,
                 smash->loc.x,smash->loc.y,smash->loc.z);
          printf("  LLF corner: %.3f %.3f %.3f\n",
                 world->x_fineparts[sv->llf.x],
                 world->y_fineparts[sv->llf.y],
                 world->z_fineparts[sv->llf.z]);
          printf("  URB corner: %.3f %.3f %.3f\n",
                 world->x_fineparts[sv->urb.x],
                 world->y_fineparts[sv->urb.y],
                 world->z_fineparts[sv->urb.z]);
          for (twlp = sv->wall_head; twlp != NULL; twlp = twlp->next)
          {
            printf("    Wall %x: [%.2f %.2f %.2f] [%.2f %.2f %.2f] [%.2f %.2f %.2f]\n",
                   (int)twlp->this_wall,
                   twlp->this_wall->vert[0]->x,
                   twlp->this_wall->vert[0]->y,
                   twlp->this_wall->vert[0]->z,
                   twlp->this_wall->vert[1]->x,
                   twlp->this_wall->vert[1]->y,
                   twlp->this_wall->vert[1]->z,
                   twlp->this_wall->vert[2]->x,
                   twlp->this_wall->vert[2]->y,
                   twlp->this_wall->vert[2]->z);
          }
        }
#endif

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
          m->properties->population--;
          m->properties = NULL;
          if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
          if (shead != NULL) mem_put_list(sv->mem->coll,shead);
          
          return NULL;
        }
        else
        {
/*          printf("Moving molecule %x from subvolume %x to %x ",
                 (int)m,(int)sv,(int)nsv);*/
          m = migrate_molecule(m,nsv);
/*          printf("and identity is now %x.\n",(int)m);*/
        }

        if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
        if (shead != NULL) mem_put_list(sv->mem->coll,shead);
        
        calculate_displacement = 0;
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
  double stop_time,max_time;
  int i,j;
  
  while ( (a = (struct abstract_molecule*)schedule_next(local->timer)) != NULL )
  {

    if (a->properties == NULL)  /* Defunct!  Remove. */
    {
      if ((a->flags & IN_MASK) == IN_SCHEDULE) mem_put(a->birthplace,a);
      else a->flags -= IN_SCHEDULE;
      
      continue;
    }

    a->flags -=IN_SCHEDULE;
    
    if (local->current_time + local->max_timestep < checkpt_time) stop_time = local->max_timestep;
    else stop_time = checkpt_time - local->current_time;
    
    if (a->t2 < EPS_C)
    {
      if ((a->flags & (ACT_INERT+ACT_NEWBIE)) != 0)
      {
        a->flags -= (a->flags & (ACT_INERT + ACT_NEWBIE));
        if ((a->flags & ACT_REACT) != 0)
        {
          r = trigger_unimolecular(a->properties->hashval,a);
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
      }
    }
    else a->t += max_time;

    a->flags += IN_SCHEDULE;
    if (a->flags & TYPE_GRID) schedule_add(local->timer,a);
    else schedule_add(((struct molecule*)a)->subvol->mem->timer,a);
    
  }
  
  local->current_time += 1.0;
}
