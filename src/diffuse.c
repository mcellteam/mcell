/**************************************************************************\
** File: diffuse.c                                                        **
**                                                                        **
** Purpose: Moves molecules around the world with reactions and collisions**
**                                                                        **
** Testing status: compiles.                                              **
\**************************************************************************/



#include <math.h>
#include <stdio.h>

#include "rng.h"
#include "mem_util.h"
#include "sched_util.h"

#include "mcell_structs.h"
#include "vol_util.h"
#include "wall_util.h"
#include "react.h"



#define MULTISTEP_WORTHWHILE 1.414
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
  struct collision *smash,*shead,*shead2;
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
      smash->next = shead2;
      shead = smash;
      smash = (struct collision*) mem_get(sv->mem->coll);
    }
  }

  if (v->x < 0.0)
  {
    dx = m->pos.x - world->x_fineparts[ sv->llf.x ];
    i = 0;
  }
  else
  {
    dx = world->x_fineparts[ sv->urb.x ] - m->pos.x;
    i = 1;
  }

  if (v->y < 0.0)
  {
    dy = m->pos.y - world->y_fineparts[ sv->llf.y ];
    j = 0;
  }
  else
  {
    dy = world->y_fineparts[ sv->urb.y ] - m->pos.y;
    j = 1;
  }

  if (v->z < 0.0)
  {
    dz = m->pos.z - world->z_fineparts[ sv->llf.z ];
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
    if (ty < dz)
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
  

/*************************************************************************
diffuse_3D:
  In: molecule that is moving
      linked list of potential collisions with molecules (we could react)
      subvolume that we start in
      displacement vector from current to new location
  Out: collision list of walls and molecules we intersected along our ray
         (current subvolume only)
*************************************************************************/

int diffuse_3D(struct molecule *m,double target_time)
{
  struct vector3 displacement;
  struct collision *smash,*shead,*shead2;
  struct subvolume *sv,*oldsv;
  struct wall *w;
  struct rxn *r;
  struct molecule *mp;
  struct species *sm;
  double d2;
  double d2min = GIGANTIC;
  double steps;
  double factor;
  
  int i,j,k;
  
  sm = m->properties;
  if (sm->space_step <= 0)
  {
    m->t += m->subvol->mem->max_timestep;
    return 1;
  }
  
  oldsv = m->subvol;
  
  shead = NULL;
  for (mp = oldsv->mol_head ; mp != NULL ; mp = mp->next_v)
  {
    if (mp==m) continue;
    
    r = trigger_bimolecular(
            m->properties->hashval,mp->properties->hashval,
            (struct abstract_molecule*)m,(struct abstract_molecule*)mp,0,0
          );
    
    if (r != NULL)
    {
      smash = mem_get(oldsv->mem->coll);
      smash->target = (void*) mp;
      smash->intermediate = r;
      smash->next = shead;
      smash->what = COLLIDE_MOL;
      shead = smash;
    }
  }
  
  if (target_time - m->t > MULTISTEP_WORTHWHILE)
  {
    for (smash = shead ; smash != NULL ; smash = smash->next)
    {
      mp = (struct molecule*)smash->target;
      d2 = (m->pos.x - mp->pos.x)*(m->pos.x - mp->pos.x) +
           (m->pos.y - mp->pos.y)*(m->pos.x - mp->pos.y) +
           (m->pos.z - mp->pos.z)*(m->pos.x - mp->pos.z);
      if (d2 < d2min) d2min = d2;
    }
  
    d2 = (m->pos.x - world->x_fineparts[ oldsv->llf.x ]);
    d2 *= d2;
    if (d2 < d2min) d2min = d2;
    
    d2 = (m->pos.x - world->x_fineparts[ oldsv->urb.x ]);
    d2 *= d2;
    if (d2 < d2min) d2min = d2;

    d2 = (m->pos.y - world->y_fineparts[ oldsv->llf.y ]);
    d2 *= d2;
    if (d2 < d2min) d2min = d2;
    
    d2 = (m->pos.y - world->y_fineparts[ oldsv->urb.y ]);
    d2 *= d2;
    if (d2 < d2min) d2min = d2;

    d2 = (m->pos.z - world->z_fineparts[ oldsv->llf.z ]);
    d2 *= d2;
    if (d2 < d2min) d2min = d2;
    
    d2 = (m->pos.z - world->z_fineparts[ oldsv->urb.z ]);
    d2 *= d2;
    if (d2 < d2min) d2min = d2;
    
    if (d2 < sm->space_step * sm->space_step) steps = 1.0;
    else steps = d2 / (sm->space_step * sm->space_step);
    
    if (steps > target_time - m->t) steps = (1.0+EPS_C)*(target_time - m->t);
  }
  else steps = 1.0;
  
  if (steps > sqrt(10.0)) steps = sqrt(10.0);

  pick_displacement(&displacement,steps*m->properties->space_step);
  sv = m->subvol;
  
  steps *= steps;  /* Convert distance steps to time steps */
  
  do
  {
    shead2 = ray_trace(m,shead,sv,&displacement);
    
    shead2 = (struct collision*)ae_list_sort((struct abstract_element*)shead2);
    
    for (smash = shead2; smash != NULL; smash = smash->next)
    {
      if (smash->t >= 1.0)
      {
        smash = NULL;
        break;
      }

      if ( (smash->what & COLLIDE_MOL) != 0 )
      {
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
          if (shead != NULL) mem_put_list(oldsv->mem->coll,shead);
          return 0;
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
          if (i < 0) continue; /* pass through--set counters here! */
          else
          {
            j = outcome_intersect(
                    r,i,w,(struct abstract_molecule*)m,k,m->t + steps*smash->t
                  );
            if (j==1) continue; /* pass through -- set counters! */
            else if (j==0)
            {
              if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
              if (shead != NULL) mem_put_list(oldsv->mem->coll,shead);
              return 0;
            }
          }
        }
        /* default is to reflect */
        
        m->pos.x += displacement.x * smash->t;
        m->pos.y += displacement.y * smash->t;
        m->pos.z += displacement.z * smash->t;
        m->t += steps*smash->t;
        steps *= (1.0-smash->t);
        
        factor = -2.0 * (displacement.x*w->normal.x + displacement.y*w->normal.y + displacement.z*w->normal.z);
        displacement.x = (displacement.x + factor*w->normal.x) * (1.0-smash->t);
        displacement.y = (displacement.y + factor*w->normal.y) * (1.0-smash->t);
        displacement.z = (displacement.z + factor*w->normal.z) * (1.0-smash->t);
        
        break;
      }
      else /* subvolume */
      {
        /* TODO: add this case. */
      }
    }
    
    if (shead2 != NULL) mem_put_list(sv->mem->coll,shead2);
  }
  while (smash != NULL);
  
  if (shead != NULL) mem_put_list(oldsv->mem->coll,shead);
  m->pos.x += displacement.x;
  m->pos.y += displacement.y;
  m->pos.z += displacement.z;
  m->t += steps;
  return 1;
}
    


/*

Note: t_inert>0 means advance reaction times by t_inert and reschedule
      t_inert<0 means reschedule unimolecular from scratch
      t_inert==0 means a reaction happened or we're ready to move again
      t2>0 means a displacement will happen then
      t2<0 means a reaction will happen at -then
      t2==0 means that the particle only reacts and doesn't move

      Consult the following table for what is supposed to happen
                   t2>0            t2<0            t2==0
t_inert>0       rxn,resched.     inert-move       rxn,resched.
t_inert==0          rxn            move             rxn
t_inert<0         rxn+move        rxn+move        rxn+move
*/

void run_timestep(struct subvolume *sv)
{
  struct abstract_molecule *a;
  struct molecule *m;
  struct surface_molecule *s;
  struct rxn *r;
  double t;
  int i,j;
  
  while ( (a = (struct abstract_molecule*)schedule_next(sv->mem->timer)) != NULL )
  {
    if (a->properties == NULL)
    {
      a = a->next;
      continue;
    }
    
    if ((a->properties->flags & ON_GRID) != 0)  /* Grid mol, can't move */
    {
      if (a->t_inert > 0.0)
      {
        a->t += a->t_inert;
        a->t_inert = 0.0;
        schedule_add(sv->mem->timer,a);
      }
      else if (a->t_inert < 0.0)
      {
        r = trigger_unimolecular(a->properties->hashval,a);
        if (r != NULL)
        {
          a->t += timeof_unimolecular(r);
          a->t_inert = 0.0;
        }
        else a->t += sv->mem->max_timestep;
        schedule_add(sv->mem->timer,a);
      }
      else
      {
        r = trigger_unimolecular(a->properties->hashval,a);
        if (r != NULL)
        {
          i = which_unimolecular(r);
          j = outcome_unimolecular(r,i,a,a->t);
          if (j)  /* mol. still exists */
          {
            a->t += timeof_unimolecular(r);
            schedule_add(sv->mem->timer,a);
          }
        }
        else
        {
          a->t += sv->mem->max_timestep;
          schedule_add(sv->mem->timer,a);
        }
      }
    }
    
    else if ((a->properties->flags & ON_SURFACE) != 0)
    {
      s = (struct surface_molecule*) a;

      if (s->t_inert > 0.0)
      {
        if (s->t2 >= 0.0)
        {
          t += s->t + s->t_inert;
          if (s->t2 > 0.0 && s->t2 < t)
          {
            s->t = s->t2;
            s->t2 = - s->t * (1.0 + EPS_C);
            s->t_inert = t + s->t2;
          }
          else
          {
            s->t = t;
            s->t_inert = 0;
          }
          schedule_add(sv->mem->timer,a);
        }
        else
        {
          /* 2D moving not implemented. */
          s->t = -(s->t2 + s->t_inert);
          s->t2 = 0.0;
          schedule_add(sv->mem->timer,a);
        }
      }
      else if (s->t_inert < 0.0)
      {
        r = trigger_unimolecular(a->properties->hashval,a);
        if (r != NULL)
        {
          a->t += timeof_unimolecular(r);
          a->t_inert = 0.0;
        }
        else a->t += sv->mem->max_timestep;
        /* Add 2D movement code here! */
        
        schedule_add(sv->mem->timer,a);
      }
      else
      {
        if (s->t2 >= 0)
        {
          r = trigger_unimolecular(a->properties->hashval,a);
          if (r != NULL)
          {
            i = which_unimolecular(r);
            j = outcome_unimolecular(r,i,a,a->t);
            if (j) /* still exists */
            {
              a->t += timeof_unimolecular(r);
              schedule_add(sv->mem->timer,a);
            }
          }
          else
          {
            a->t += sv->mem->max_timestep;
            schedule_add(sv->mem->timer,a);
          }
        }
        else
        {
          /* 2D moving not implemented. */
          s->t = -(s->t2 + s->t_inert);
          s->t2 = 0.0;
          schedule_add(sv->mem->timer,a);
        }
      }
    }
    
    else /* freely diffusing molecule */
    {

      m = (struct molecule*) a;
      
      if (m->t_inert > 0.0) /* inert */
      {

        if (m->t2 >= 0.0)  /* react */
        {
          t = m->t + m->t_inert;
          if (m->t2 > 0.0 && m->t2 < t)
          {
            m->t = m->t2;
            m->t2 = - m->t * (1.0 + EPS_C);
            m->t_inert = t + m->t2;
          }
          else
          {
            m->t = t;
            m->t_inert = 0;
          }
          schedule_add(sv->mem->timer,a);
        }
        else /* diffuse */
        {
          t = m->t;
          j = diffuse_3D(m,m->t_inert);
          if (j) /* still exists */
          {
            if (m->t - t > m->t_inert)
            {
              m->t2 -= m->t_inert;
              m->t_inert = 0.0;
            }
            else
            {
              m->t2 -= m->t - t;
              m->t_inert -= m->t - t;
            }
            
            if (m->t < -m->t2) schedule_add(sv->mem->timer,a);
            else
            {
              t = m->t;
              m->t = -m->t2;
              m->t2 = t;
              schedule_add(sv->mem->timer,a);
            }
          }
        }
      
      }
      else if (m->t_inert < 0.0) /* schedule */
      {

        r = trigger_unimolecular(a->properties->hashval,a);
        if (r != NULL) t = m->t + timeof_unimolecular(r);
        else t = m->t + sv->mem->max_timestep;
        m->t_inert = 0.0;
        
        j = diffuse_3D(m,t);
        
        if (j) /* still exists */
        {
          if (m->t < t) m->t2 = -t;
          else
          {
            m->t2 = m->t;
            m->t = t;
          }
          
          schedule_add(sv->mem->timer,a);
        }
        
      }
      else /* diffuse */
      {
        j = diffuse_3D(m,-m->t2);
        if (j) /* still exists */
        {
          if (m->t > -m->t2)
          {
            t = m->t;
            m->t = -m->t2;
            m->t2 = t;
          }
          
          schedule_add(sv->mem->timer,a);
        }
        
      }

    }
    
    a = a->next;
  }
}
