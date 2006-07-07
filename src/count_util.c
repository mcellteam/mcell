/**************************************************************************\
** File: count_util.c                                                     **
**                                                                        **
** Purpose: Handles counting of interesting events                        **
**                                                                        **
** Testing status: untested.                                              **
\**************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rng.h"
#include "grid_util.h"
#include "mcell_structs.h"
#include "util.h"
#include "wall_util.h"
#include "vol_util.h"
#include "count_util.h"
#include "react_output.h"

extern struct volume *world;


/*************************************************************************
eps_equals:
   In: two doubles
   Out: 1 if they are equal to within some small tolerance, 0 otherwise
*************************************************************************/

int eps_equals(double x,double y)
{
  double mag;
  double diff;
  
  if (x<0) mag = -x;
  else mag = x;
  if (y<0) { if (-y > mag) mag = -y; }
  else { if (y > mag) mag = y; }
  
  if (x < y) diff = y-x;
  else diff = x-y;
  
  return diff < EPS_C * (mag + 1.0);
}


/*************************************************************************
update_collision_count:
   In: species of thing that hit
       region list for the wall we hit
       direction of impact relative to surface normal
       whether we crossed or not
       scaling factor for reaction probabilities (for estimating ccn)
   Out: No return value.  Appropriate counters are updated.
*************************************************************************/

void update_collision_count(struct species *sp,struct region_list *rl,int direction,int crossed, double factor,struct vector3 *loc,double t)
{
  int j;
  struct counter *hit_count;
  double hits_to_ccn=0;
  
  if (sp->flags&COUNT_HITS)
  {
    hits_to_ccn = sp->time_step * 2.9432976599069717358e-3 /  /* 1e6*sqrt(MY_PI)/(1e-15*N_AV) */ 
                  (sp->space_step*factor*world->length_unit*world->length_unit*world->length_unit);
  }
  
  hit_count = NULL;  
  for ( ; rl != NULL ; rl = rl->next)
  {
    if (rl->reg->flags & COUNT_SOME)
    {
      j = (rl->reg->hashval ^ sp->hashval)&world->count_hashmask;
      if (j==0) j = sp->hashval & world->count_hashmask;
      
      for (hit_count=world->count_hash[j] ; hit_count!=NULL ; hit_count=hit_count->next)
      {
        if (hit_count->reg_type == rl->reg && hit_count->target == sp)
        {
          if (rl->reg->flags & sp->flags & COUNT_HITS)
          {
            if (crossed)
            {
              if (direction==1)
              {
                if (hit_count->counter_type&TRIG_COUNTER)
                {
                  hit_count->data.trig.t_event=t;
                  fire_count_event(hit_count,1,NULL,loc,REPORT_FRONT_HITS|REPORT_TRIGGER);
                  fire_count_event(hit_count,1,NULL,loc,REPORT_FRONT_CROSSINGS|REPORT_TRIGGER);
                }
                else
                {
                  hit_count->data.move.front_hits++;
                  hit_count->data.move.front_to_back++;
                }
              }
              else
              {
                if (hit_count->counter_type&TRIG_COUNTER)
                {
                  hit_count->data.trig.t_event=t;
                  fire_count_event(hit_count,1,NULL,loc,REPORT_BACK_HITS|REPORT_TRIGGER);
                  fire_count_event(hit_count,1,NULL,loc,REPORT_BACK_CROSSINGS|REPORT_TRIGGER);
                }
                else
                {
                  hit_count->data.move.back_hits++;
                  hit_count->data.move.back_to_front++;
                }
              }
            }
            else
            {
              if (direction==1)
              {
                if (hit_count->counter_type&TRIG_COUNTER)
                {
                  hit_count->data.trig.t_event=t;
                  fire_count_event(hit_count,1,NULL,loc,REPORT_FRONT_HITS|REPORT_TRIGGER);
                }
                else
                {
                  hit_count->data.move.front_hits++;
                }
              }
              else
              {
                if (hit_count->counter_type&TRIG_COUNTER)
                {
                  hit_count->data.trig.t_event=t;
                  fire_count_event(hit_count,1,NULL,loc,REPORT_BACK_HITS|REPORT_TRIGGER);
                }
                else hit_count->data.move.back_hits++;
              }
            }
	    if (rl->reg->area != 0.0)
	    {
	      if ((hit_count->counter_type&TRIG_COUNTER)==0)
              {
                hit_count->data.move.scaled_hits += factor*hits_to_ccn/rl->reg->area;
              }
	    }
          }
        }
      }
    }
  }
}


/*************************************************************************
find_enclosing_regions:
   In: location we want to end up
       starting position
       list of regions we're inside at the starting position
       list of inside-out regions we're "outside" at the starting position
       memory handler to store lists of regions
   Out: 0 on success, 1 on memory allocation error.  The region and
        inside-out region lists are updated to be correct at the ending
	position.
*************************************************************************/

int find_enclosing_regions(struct vector3 *loc,struct vector3 *start,
                            struct region_list** rlp,struct region_list** arlp,
                            struct mem_helper *rmem)
{
  struct vector3 outside,delta,hit;
  struct subvolume *sv,*svt;
  struct wall_list *wl;
  struct region_list *rl,*arl;
  struct region_list *trl,*tarl,*xrl,*yrl,*nrl;
  double t,t_hit_sv;
  int traveling;
  int i;
  struct wall_list dummy;
  
  rl = *rlp;
  arl = *arlp;
  
  if (start==NULL || loc->x!=start->x || loc->y!=start->y || loc->z < start->z)
  {
    outside.x = loc->x;
    outside.y = loc->y;
    outside.z = (world->z_partitions[0] + world->z_partitions[1])/2;
  }
  else
  {
    outside.x = start->x;
    outside.y = start->y;
    outside.z = start->z;
  }
  
  delta.x = 0.0;
  delta.y = 0.0;
  delta.z = loc->z - outside.z;
  
  sv = find_subvolume(&outside,NULL);
  svt = find_subvolume(loc,NULL);
  traveling = 1;

  while (traveling)
  {
    tarl = trl = NULL;
    t_hit_sv = collide_sv_time(&outside,&delta,sv);
    
    for (wl = sv->wall_head ; wl != NULL ; wl = wl->next)
    {
      i = collide_wall(&outside , &delta , wl->this_wall , &t , &hit);
      if (i==COLLIDE_REDO)
      {
        while (trl != NULL)
        {
          xrl = trl->next;
          mem_put(rmem,trl);
          trl = xrl;
        }
        while (tarl != NULL)
        {
          xrl = tarl->next;
          mem_put(rmem,tarl);
          tarl = xrl;
        }
        dummy.next = sv->wall_head;
        wl = &dummy;
        continue;  /* Trick to restart for loop */
      }
      else if (i==COLLIDE_MISS || !(t >= 0 && t < 1.0) || t > t_hit_sv || (wl->this_wall->flags & COUNT_CONTENTS) == 0 ||
	       (hit.x-outside.x)*delta.x + (hit.y-outside.y)*delta.y + (hit.z-outside.z)*delta.z < 0) continue;
      else
      {
        for (xrl=wl->this_wall->counting_regions ; xrl != NULL ; xrl = xrl->next)
        {
          if ((xrl->reg->flags & COUNT_CONTENTS) != 0)
          {
            nrl = (struct region_list*) mem_get(rmem);
	    if (nrl==NULL)
	    {
	      fprintf(stderr, "File '%s', Line %ld:  Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
	      i = emergency_output();
	      fprintf(stderr, "Fatal error: out of memory while finding enclosing regions.\nAttempt to write intermediate results had %d errors\n", i);
	      exit(EXIT_FAILURE);
            }

            nrl->reg = xrl->reg;
            
            if (i==COLLIDE_BACK) { nrl->next = tarl; tarl = nrl; }
            else { nrl->next = trl; trl = nrl; }
          }
        }
      }
    }
    
    xrl = trl;
    while (trl != NULL)
    {
      nrl = NULL;
      yrl = arl;
      while (yrl != NULL)
      {
        if (xrl->reg == yrl->reg)
        {
          if (nrl==NULL)
          {
            arl = yrl->next;
            mem_put(rmem,yrl);
            yrl = arl;
          }
          else
          {
            nrl->next = yrl->next;
            mem_put(rmem,yrl);
            yrl = nrl;
          }
          trl = trl->next;
          mem_put(rmem,xrl);
          xrl = NULL;
          break;
        }
        else
        {
          nrl = yrl;
          yrl = yrl->next;
        }
      }
      if (xrl!=NULL)
      {
        trl = trl->next;
        xrl->next = rl;
        rl = xrl;
        xrl = trl;
      }
      else xrl = trl;
    }
    
    xrl = tarl;
    while (tarl != NULL)
    {
      nrl = NULL;
      yrl = rl;
      while (yrl != NULL)
      {
        if (xrl->reg == yrl->reg)
        {
          if (nrl==NULL)
          {
            rl = yrl->next;
            mem_put(rmem,yrl);
            yrl = rl;
          }
          else
          {
            nrl->next = yrl->next;
            mem_put(rmem,yrl);
            yrl = nrl;
          }
          tarl = tarl->next;
          mem_put(rmem,xrl);
          xrl = NULL;
          break;
        }
        else
        {
          nrl = yrl;
          yrl = yrl->next;
        }
      }
      if (xrl!=NULL)
      {
        tarl = tarl->next;
        xrl->next = arl;
        arl = xrl;
        xrl = tarl;
      }
      else xrl = tarl;
    }
    
    if (sv==svt) traveling = 0;
    else
    {
      sv = next_subvol( &outside , &delta , sv );
      delta.x = loc->x - outside.x;
      delta.y = loc->y - outside.y;
      delta.z = loc->z - outside.z;
      
      if (sv == NULL)
      {
	if ((delta.x*delta.x + delta.y*delta.y + delta.z*delta.z) < EPS_C*EPS_C)
	{
	  fprintf(world->log_file, "File '%s', Line %ld: Didn't quite reach waypoint target, fudging.\n", __FILE__, (long)__LINE__);
	  traveling = 0;
	}
	else
	{
	  fprintf(world->log_file, "File '%s', Line %ld: Couldn't reach waypoint target.\n", __FILE__, (long)__LINE__);
	  sv = find_subvolume(&outside , NULL);
	}
      }
    }
  }
  
  *rlp = rl;
  *arlp = arl;

  return 0;
}



/*************************************************************************
dup_region_list:
   In: a list of regions
       memory handler to use for duplicated regions
   Out: The duplicated list of regions, or NULL on a memory allocation
        error.
*************************************************************************/

struct region_list* dup_region_list(struct region_list *r,struct mem_helper *mh)
{
  struct region_list *nr,*rp,*r0;
  
  if (r==NULL) return NULL;
  
  r0 = rp = NULL;
  while (r!=NULL)
  {
    nr = (struct region_list*) mem_get(mh);
    if(nr == NULL) return NULL;

    nr->next = NULL;
    nr->reg = r->reg;
    if (rp==NULL) r0 = rp = nr;
    else { rp->next = nr; rp = nr; }
    
    r = r->next;
  }
  
  return r0;
}



/*************************************************************************
place_waypoints:
   In: No arguments.
   Out: Returns 1 if malloc fails, 0 otherwise.
        Allocates waypoints to SSVs, if any are needed.
   Note: you must have initialized SSVs before calling this routine!
*************************************************************************/

int place_waypoints()
{
  int g,h,i,j,k;
  int i_will_use_waypoints = 0;
  int waypoint_in_wall = 0;
  struct waypoint *wp;
  struct wall_list *wl;
  struct subvolume *sv;
  double d;
  
/* Being exactly in the center of a subdivision can be bad. */
/* Define "almost center" positions for X, Y, Z */
#define W_Xa (0.5 + 0.0005*MY_PI)
#define W_Ya (0.5 + 0.0002*MY_PI*MY_PI)
#define W_Za (0.5 - 0.00007*MY_PI*MY_PI*MY_PI)
#define W_Xb (1.0 - W_Xa)
#define W_Yb (1.0 - W_Ya)
#define W_Zb (1.0 - W_Za)

  /* Probably ought to check for whether you really need waypoints for releases */
  if (!world->releases_on_regions_flag)
  {
    for (i=0;i<world->n_species;i++)
    {
      if ((world->species_list[i]->flags & (NOT_FREE | COUNT_CONTENTS)) == COUNT_CONTENTS
          || (world->species_list[i]->flags & COUNT_ENCLOSED)!=0 )
        i_will_use_waypoints++;
    }
    
    if (i_will_use_waypoints==0)
    {
      world->n_waypoints = 0;
      world->waypoints = NULL;
      return 0;  /* Waypoints not needed. */
    }
  }
  
  world->n_waypoints = world->n_subvols;
  world->waypoints = (struct waypoint*)malloc(sizeof(struct waypoint)*world->n_waypoints);
  if (!world->waypoints) return 1;

  for (i=0;i<world->nx_parts-1;i++)
  {
    for (j=0;j<world->ny_parts-1;j++)
    {
      for (k=0;k<world->nz_parts-1;k++)
      {
        h = k + (world->nz_parts-1)*(j + (world->ny_parts-1)*i);
        wp = &(world->waypoints[h]);
        
        sv = &(world->subvol[h]);
        
        /* Place waypoint near center of subvolume (W_#a=W_#b=0.5 gives center) */
        wp->loc.x = W_Xa*world->x_fineparts[ sv->llf.x ] + W_Xb*world->x_fineparts[ sv->urb.x ];
        wp->loc.y = W_Ya*world->y_fineparts[ sv->llf.y ] + W_Yb*world->y_fineparts[ sv->urb.y ];
        wp->loc.z = W_Za*world->z_fineparts[ sv->llf.z ] + W_Zb*world->z_fineparts[ sv->urb.z ];
        
        do
        {
          waypoint_in_wall = 0;
          for (wl = sv->wall_head ; wl != NULL ; wl = wl->next)
          {
            d = dot_prod( &(wp->loc) , &(wl->this_wall->normal) ); 
            if ( eps_equals( d , wl->this_wall->d ) )
            { 
              waypoint_in_wall++;
              d = EPS_C * (double)((rng_uint(world->rng)&0xF) - 8);
              if (d==0) d = 8*EPS_C;
              wp->loc.x += d * wl->this_wall->normal.x;
              wp->loc.y += d * wl->this_wall->normal.y;
              wp->loc.z += d * wl->this_wall->normal.z;
              break;
            }
          }
        } while (waypoint_in_wall);
        
        if (k>0)
        {
	  if (world->waypoints[h-1].regions != NULL)
	  {
            wp->regions = dup_region_list(world->waypoints[h-1].regions,sv->local_storage->regl);
	    if (wp->regions == NULL) return 1;
	  }
	  else wp->regions = NULL;
	  
	  if (world->waypoints[h-1].antiregions != NULL)
	  {
            wp->antiregions = dup_region_list(world->waypoints[h-1].antiregions,sv->local_storage->regl);
	    if (wp->antiregions == NULL) return 1;
	  }
	  else wp->antiregions = NULL;
          
          g = find_enclosing_regions(&(wp->loc),&(world->waypoints[h-1].loc),
                                     &(wp->regions),&(wp->antiregions),sv->local_storage->regl);
        }
        else
        {
          wp->regions = NULL;
          wp->antiregions = NULL;
          g = find_enclosing_regions(&(wp->loc),NULL,&(wp->regions),
                                 &(wp->antiregions),sv->local_storage->regl);	  
        }
	if (g) return 1;
      }
    }
  }
  
  return 0;
#undef W_Zb
#undef W_Yb
#undef W_Xb
#undef W_Za
#undef W_Ya
#undef W_Xa  
}



/*************************************************************************
region_listed:
   In: list of regions
       one specific region we're interested in
   Out: 1 if the region is in the list.  0 if not.
*************************************************************************/

int region_listed(struct region_list *rl,struct region *r)
{
  while (rl!=NULL)
  {
    if (rl->reg==r) return 1;
    rl=rl->next;
  }
  return 0;
}



/*************************************************************************
count_me_by_region:
   In: abstract molecule we are supposed to count (or a representative one)
       number by which to update the counter (usually +1 or -1)
       named reaction pathway to count at location of abstract molecule
       time of event (used for triggers)
   Out: No return value.  Appropriate counters are updated.
   Note: This handles all types of molecules, from grid to free.  Grid
         and free are implemented, surface is not.  If the reaction
	 pathname is NULL, the molecule is counted.  Otherwise, the
	 reaction is counted at the location of the molecule.
*************************************************************************/

void count_me_by_region(struct abstract_molecule *me,int n,struct rxn_pathname *rxp,double t)
{
  int i,j,k,h;
  struct region_list *rl;
  struct species *sp = me->properties;
  struct counter *c;
  
  u_int desired_hash;
  u_int COUNT_flag;
  
  if (rxp!=NULL)
  {
    desired_hash = rxp->hashval;
    COUNT_flag = COUNT_RXNS;
  }
  else
  {
    desired_hash = sp->hashval;
    COUNT_flag = COUNT_CONTENTS;
  }
  
  //printf("Counting %s by region (up by %d) iter %d!\n",me->properties->sym->name,n,(int)world->it_time);
  
  if ((sp->flags & ON_GRID) != 0)
  {
    struct grid_molecule *g = (struct grid_molecule*)me;
    struct wall *w = g->grid->surface;

    if (w->flags & COUNT_flag)
    {
      for (rl=w->counting_regions ; rl!=NULL ; rl=rl->next)
      {
        i = (rl->reg->hashval ^ desired_hash) & world->count_hashmask;
        if (i==0) i = desired_hash & world->count_hashmask;
        
        for ( c = world->count_hash[i] ; c != NULL ; c = c->next )
        {
          if (c->counter_type&TRIG_COUNTER) c->data.trig.t_event=t;
          
          if (c->reg_type == rl->reg && (c->counter_type&ENCLOSING_COUNTER)==0)
	  {
            if (c->counter_type&TRIG_COUNTER) fire_count_event(c,n,g,NULL,REPORT_CONTENTS|REPORT_TRIGGER);
            else c->data.move.n_at += n;
          }
	  else if (c->target == rxp)
          {
            if (c->counter_type&TRIG_COUNTER) fire_count_event(c,n,g,NULL,REPORT_RXNS|REPORT_TRIGGER);            
            c->data.rx.n_rxn_at += n;
	  }
        }
      }
    } 
  }
  
  if ((sp->flags & NOT_FREE)==0 || (sp->flags&COUNT_ENCLOSED)!=0) /* Free molecule */
  {
    struct molecule *m;
    struct grid_molecule *g;
    struct subvolume *sv;
    struct vector3 here,delta,hit;
    struct waypoint *wp;
    struct wall_list *wl;
    struct region_list *rl;
    struct vector3 loc;
    double t;
    double t_sv_hit;
    struct wall *my_wall;
    
    if ((sp->flags&NOT_FREE)==0)
    {
      m = (struct molecule*)me;
      g = NULL;
      my_wall=NULL;
      
      loc.x = m->pos.x;
      loc.y = m->pos.y;
      loc.z = m->pos.z;
    }
    else /* Grid mol */
    {
      g = (struct grid_molecule*)me;
      m = NULL;
      my_wall=g->grid->surface;
      uv2xyz(&(g->s_pos),my_wall,&loc);
    }
      
    i = bisect(world->x_partitions,world->nx_parts,loc.x);
    j = bisect(world->y_partitions,world->ny_parts,loc.y);
    k = bisect(world->z_partitions,world->nz_parts,loc.z);
    
    h = k + (world->nz_parts-1)*( j + (world->ny_parts-1)*i );
    wp = &(world->waypoints[h]);
    for (rl=wp->regions ; rl!=NULL ; rl=rl->next)
    {
      if ( (rl->reg->flags & COUNT_flag) != 0 )
      {
        i = (rl->reg->hashval ^ desired_hash) & world->count_hashmask;
        if (i==0) i = desired_hash & world->count_hashmask;
        
        for ( c = world->count_hash[i] ; c != NULL ; c = c->next )
        {
          if (c->counter_type&TRIG_COUNTER) c->data.trig.t_event=t;
          
          if (c->reg_type==rl->reg && (my_wall==NULL || !region_listed(my_wall->counting_regions,rl->reg)))
	  {
	    if ( rxp==NULL && c->target==sp &&
	         (g==NULL || (c->counter_type&ENCLOSING_COUNTER)!=0) )
	    {
              if (c->counter_type&TRIG_COUNTER) fire_count_event(c,n,NULL,&loc,REPORT_CONTENTS|REPORT_ENCLOSED|REPORT_TRIGGER);
	      else c->data.move.n_enclosed += n;
	    }
	    else if ( c->target==rxp && (c->counter_type&ENCLOSING_COUNTER)!=0 )
	    {
              if (c->counter_type&TRIG_COUNTER) fire_count_event(c,n,NULL,&loc,REPORT_RXNS|REPORT_ENCLOSED|REPORT_TRIGGER);
	      else c->data.rx.n_rxn_enclosed += n;
	    }
	  }
        }
      }
    }
    for (rl=wp->antiregions ; rl!=NULL ; rl=rl->next)
    {
      if ( (rl->reg->flags & COUNT_flag) != 0 )
      {
        i = (rl->reg->hashval ^ desired_hash) & world->count_hashmask;
        if (i==0) i = desired_hash & world->count_hashmask;
        
        for ( c = world->count_hash[i] ; c != NULL ; c = c->next )
        {
          if (c->counter_type&TRIG_COUNTER) c->data.trig.t_event=t;
          
          if (c->reg_type==rl->reg && (my_wall==NULL || !region_listed(my_wall->counting_regions,rl->reg)))
	  {
	    if ( rxp==NULL && c->target==sp &&
	         (g==NULL || (c->counter_type&ENCLOSING_COUNTER)!=0) )
	    {
              if (c->counter_type&TRIG_COUNTER) fire_count_event(c,-n,NULL,&loc,REPORT_CONTENTS|REPORT_ENCLOSED|REPORT_TRIGGER);
	      else c->data.move.n_enclosed -= n;
	    }
	    else if ( c->target==rxp && (c->counter_type&ENCLOSING_COUNTER)!=0 )
	    {
              if (c->counter_type&TRIG_COUNTER) fire_count_event(c,-n,NULL,&loc,REPORT_RXNS|REPORT_ENCLOSED|REPORT_TRIGGER);
	      else c->data.rx.n_rxn_enclosed -= n;
	    }
	  }
        }
      }
    }
    
    here.x = wp->loc.x;
    here.y = wp->loc.y;
    here.z = wp->loc.z;
    
    for ( sv = &(world->subvol[h]) ; sv != NULL ; sv = next_subvol(&here,&delta,sv) )
    {
      delta.x = loc.x - here.x;
      delta.y = loc.y - here.y;
      delta.z = loc.z - here.z;
      
      t_sv_hit = collide_sv_time(&here,&delta,sv);
      if (t_sv_hit > 1.0) t_sv_hit = 1.0;

      for (wl = sv->wall_head ; wl != NULL ; wl = wl->next)
      {
	if (wl->this_wall==my_wall) continue;  /* If we're on a wall, skip it */
	
        if (wl->this_wall->flags & COUNT_flag)
        {
          j = collide_wall(&here,&delta,wl->this_wall,&t,&hit);
          
          if (j!=COLLIDE_MISS && t <= t_sv_hit &&
	    (hit.x-loc.x)*delta.x + (hit.y-loc.y)*delta.y + (hit.z-loc.z)*delta.z < 0)
          {
            for (rl=wl->this_wall->counting_regions ; rl!=NULL ; rl=rl->next)
            {
              if ( (rl->reg->flags & COUNT_flag) != 0 )
              {
                i = (rl->reg->hashval ^ desired_hash) & world->count_hashmask;
                if (i==0) i = desired_hash & world->count_hashmask;
                
                for ( c = world->count_hash[i] ; c != NULL ; c = c->next )
                {
                  if (c->counter_type&TRIG_COUNTER) c->data.trig.t_event=t;
                  
                  if (c->reg_type==rl->reg && (my_wall==NULL || !region_listed(my_wall->counting_regions,rl->reg)))
		  {
		    if ( rxp==NULL && c->target==sp &&
		         (g==NULL || (c->counter_type&ENCLOSING_COUNTER)!=0) )
		    {
		      if (j==COLLIDE_FRONT)
		      {
                        if (c->counter_type&TRIG_COUNTER) fire_count_event(c,n,NULL,&loc,REPORT_CONTENTS|REPORT_ENCLOSED|REPORT_TRIGGER);
			c->data.move.n_enclosed += n;
		      }
		      else if (j==COLLIDE_BACK)
		      {
                        if (c->counter_type&TRIG_COUNTER) fire_count_event(c,-n,NULL,&loc,REPORT_CONTENTS|REPORT_ENCLOSED|REPORT_TRIGGER);
			c->data.move.n_enclosed -= n;
		      }
		    }
		    else if (c->target==rxp && (c->counter_type&ENCLOSING_COUNTER)!=0)
		    {
		      if (j==COLLIDE_FRONT)
                      {
                        if (c->counter_type&TRIG_COUNTER) fire_count_event(c,n,NULL,&loc,REPORT_RXNS|REPORT_ENCLOSED|REPORT_TRIGGER);
                        c->data.rx.n_rxn_enclosed += n;
                      }
		      else if (j==COLLIDE_BACK)
                      {
                        if (c->counter_type&TRIG_COUNTER) fire_count_event(c,-n,NULL,&loc,REPORT_RXNS|REPORT_ENCLOSED|REPORT_TRIGGER);                        
                        c->data.rx.n_rxn_enclosed -=n;
                      }
		    }
		  }
                }
              }
            }            
          }
        }
      }
    }
  }
}


/******************************************************************
prepare_counters:
  In: No arguments.
  Out: 0 if counter statements are correct, 1 otherwise.
  Note: A statement is incorrect if a non-closed manifold region
        tries to count a freely diffusing molecule.  Fixes up all
        count requests to point at the data we care about.
********************************************************************/
int prepare_counters()
{
  struct output_request *request;
  struct counter *cp;
  struct region *rp;
  struct output_block *block;
  struct output_set *set;
  struct output_column *column;
  u_int i;
  
  /* First give everything a sensible name, if needed */
  for (block=world->output_block_head ; block!=NULL ; block=block->next)
  {
    for (set=block->data_set_head ; set!=NULL ; set=set->next)
    {
      if (set->header_comment==NULL) continue;
      for (column=set->column_head ; column!=NULL ; column=column->next)
      {
        if (column->data_type==TRIG_STRUCT) continue;
        if (column->expr->title==NULL) column->expr->title = oexpr_title(column->expr);
        if (column->expr->title==NULL)
        {
          fprintf(world->err_file,"Out of memory: file %s, line %d\n  Unable to create title for data output.",__FILE__,__LINE__);
          return 1;
        }
      }
    }
  }

  /* Then convert all requests to real counts */
  for (request=world->output_request_head ; request!=NULL ; request=request->next)
  {
    if (request->count_location!=NULL && request->count_location->sym_type==OBJ)
    {
      i = expand_object_output(request,(struct object*)request->count_location->value);
      if (i)
      {
        fprintf(world->err_file,"Error: unable to expand request to count on object");
        return 1;
      }
    }
    
    i = instantiate_request(request);
    if (i)
    {
      fprintf(world->err_file,"Error: unable to count as requested\n");
      return 1;
    }
  }
  /* Need to keep all the requests for now...could repackage them to save memory */
  
  /* Now check to make sure what we've created is geometrically sensible */
  for (i=0;i<world->count_hashmask+1;i++)
  {
    for (cp=world->count_hash[i];cp!=NULL;cp=cp->next)
    {
      if ( (cp->counter_type & ENCLOSING_COUNTER) != 0)
      {
        rp=cp->reg_type;
	
	if (rp->manifold_flag==MANIFOLD_UNCHECKED)
        {
	  if (is_manifold(rp)) rp->manifold_flag=IS_MANIFOLD;
	  else rp->manifold_flag=NOT_MANIFOLD;
	}
	
	if (rp->manifold_flag==NOT_MANIFOLD)
        {
	  fprintf(world->err_file,"File '%s', Line %ld: error, cannot count molecules or events inside non-manifold object region: %s\n", __FILE__, (long)__LINE__, rp->sym->name); 
	  return (1);
	}
      }
    }
  }

  return(0);
}


int expand_object_output(struct output_request *request,struct object *obj)
{
  struct output_request *new_request;
  struct output_expression *oe,*oel,*oer; /* Original expression and two children */
  struct object *child;
  struct region_list *rl;
  int n_expanded;
  
  switch (obj->object_type)
  {
    case META_OBJ:
      n_expanded=0;
      for (child=obj->first_child ; child!=NULL ; child=child->next)
      {
        if (!object_has_geometry(child)) continue;  /* NOTE -- for objects nested N deep, we check this N+(N-1)+...+2+1 times (slow) */
        if (n_expanded>0)
        {
          new_request = (struct output_request*)mem_get(world->outp_request_mem);
          oe = request->requester;
          oel = new_output_expr(world->oexpr_mem);
          oer = new_output_expr(world->oexpr_mem);
          if (new_request==NULL || oel==NULL || oer==NULL)
          {
            fprintf(world->err_file,"Out of memory while expanding count expression on object %s\n",obj->sym->name);
            return 1;
          }
          oel->column=oer->column=oe->column;
          oel->expr_flags=oer->expr_flags=oe->expr_flags;
          oel->up=oer->up=oe;
          oel->left=request;
          oer->left=new_request;
          oel->oper=oer->oper='#';
          oe->expr_flags=(oe->expr_flags&OEXPR_TYPE_MASK)|OEXPR_LEFT_OEXPR|OEXPR_RIGHT_OEXPR;
          oe->left=oel;
          oe->right=oer;
          oe->oper='+';
          
          new_request->report_type=request->report_type;
          new_request->count_target=request->count_target;
          new_request->requester=oer;
          request->requester=oel;
          new_request->next=request->next;
          request->next=new_request;
          request=new_request;
        }
        if (expand_object_output(request,child)) return 1;
      }
      if (n_expanded==0)
      {
        fprintf(world->err_file,"Error: trying to count on object %s but it has no geometry\n",obj->sym->name);
        return 1;
      }
      break;
    case BOX_OBJ:
    case POLY_OBJ:
      for (rl=obj->regions ; rl!=NULL ; rl=rl->next)
      {
        if (is_reverse_abbrev(",ALL",rl->reg->sym->name)) break;
      }
      if (rl==NULL)
      {
        fprintf(world->err_file,"All region missing on object %s?\n  File %s, line %d\n",obj->sym->name,__FILE__,__LINE__);
        return 1;
      }
      request->count_location = rl->reg->sym;
      break;
    case REL_SITE_OBJ:
      break;
    default:
      fprintf(world->err_file,"Bad object type in count on object expansion\n  File %s, line %d\n",__FILE__,__LINE__);
      return 1;
      break;
  }
  return 0;
}


int object_has_geometry(struct object *obj)
{
  struct object *child;
  switch (obj->object_type)
  {
    case BOX_OBJ:
    case POLY_OBJ:
      return 1;
      break;
    case META_OBJ:
      for (child=obj->first_child ; child!=NULL ; child=child->next)
      {
        if (object_has_geometry(child)) return 1;
      }
      break;
    default:
      return 0;
      break;
  }
  return 0;
}


int instantiate_request(struct output_request *request)
{
  int request_hash;
  struct rxn_pathname *rxpn_to_count;
  struct rxn *rx_to_count;
  struct species *mol_to_count;
  void *to_count;
  struct region *reg_of_count;
  struct counter *count;
  struct trigger_request *trig_req;
  u_int report_type_only;
  byte count_type;
  int is_enclosed;
  
  
  /* Set up and figure out hash value */
  to_count=request->count_target->value;
  switch (request->count_target->sym_type)
  {
    case MOL:
      rxpn_to_count=NULL;
      rx_to_count=NULL;
      mol_to_count=(struct species*)to_count;
      request_hash=mol_to_count->hashval;
      break;
    case RXPN:
      rxpn_to_count=(struct rxn_pathname*)to_count;
      rx_to_count=rxpn_to_count->rx;
      mol_to_count=NULL;
      request_hash=rxpn_to_count->hashval;
      break;
    default:
      fprintf(world->err_file,"Error at file %s line %d\n  Invalid object type in count request.\n",__FILE__,__LINE__);
      return 1;
      break;
  }
  
  if (request->count_location!=NULL)
  {
    if (request->count_location->sym_type!=REG)
    {
      fprintf(world->err_file,"Error at file %s line %d\n  Non-region location in count request.\n",__FILE__,__LINE__);
      return 1;
    }
    reg_of_count=(struct region*)request->count_location->value;
    if ((request_hash^reg_of_count->hashval) != 0) request_hash ^= reg_of_count->hashval;
  }
  else reg_of_count=NULL;
  
  request_hash&=world->count_hashmask;
  
  /* Now create count structs and set output expression to point to data */
  report_type_only=request->report_type&REPORT_TYPE_MASK;
  request->requester->expr_flags-=OEXPR_LEFT_REQUEST;  
  if ((request->report_type&REPORT_TRIGGER)==0 && request->count_location==NULL) /* World count is easy! */
  {
    switch (report_type_only)
    {
      case REPORT_CONTENTS:
        request->requester->expr_flags|=OEXPR_LEFT_INT;
        request->requester->left=(void*)&(mol_to_count->population);
        break;
      case REPORT_RXNS:
        request->requester->expr_flags|=OEXPR_LEFT_DBL;
        request->requester->left=(void*)&(rx_to_count->info[rxpn_to_count->path_num].count);
        break;
      default:
        fprintf(world->err_file,"Error at file %s line %d\n  Invalid report type 0x%x in count request.\n",__FILE__,__LINE__,report_type_only);
        return 1;
        break;
    }
  }
  else /* Triggered count or count on region */
  {
    /* Set count type flags */
    if (report_type_only==REPORT_RXNS) count_type=RXN_COUNTER;
    else count_type=MOL_COUNTER;
    if (request->report_type&REPORT_ENCLOSED)
    {
      reg_of_count->flags|=ENCLOSING_COUNTER;
      count_type|=ENCLOSING_COUNTER;
      if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_ENCLOSED;
    }
    if (request->report_type&REPORT_TRIGGER) count_type|=TRIG_COUNTER;
    
    /* Find or add counter */
    for (count=world->count_hash[request_hash] ; count!=NULL ; count=count->next)
    {
      if (count->reg_type==reg_of_count && count->target==to_count && count_type==count->counter_type) break;
    }
    if (count==NULL)
    {
      count=create_new_counter(reg_of_count,request->count_target->value,count_type);
      if (count==NULL)
      {
        fprintf(world->err_file,"Error at file %s line %d\n  Out of memory allocating count request\n",__FILE__,__LINE__);
        return 1;
      }
      
      count->next=world->count_hash[request_hash];
      world->count_hash[request_hash]=count;
    }
    
    /* Point appropriately */
    if (request->report_type&REPORT_TRIGGER)
    {
      trig_req = (struct trigger_request*)mem_get(world->trig_request_mem);
      if (trig_req==NULL)
      {
        fprintf(world->err_file,"Error at file %s line %d\n  Out of memory setting notifications for a trigger\n",__FILE__,__LINE__);
        return 1;
      }
      
      trig_req->next=count->data.trig.listeners;
      count->data.trig.listeners=trig_req;
      trig_req->ear=request;
      
      request->requester->expr_flags|=OEXPR_TYPE_TRIG;
    }
    else /* Not trigger--set up for regular count */
    {
      is_enclosed = (mol_to_count==NULL || !(mol_to_count->flags&ON_GRID) || (request->report_type&REPORT_ENCLOSED));
      request->requester->expr_flags|=OEXPR_LEFT_DBL;          /* Assume double */
      switch (report_type_only)
      {
        case REPORT_CONTENTS:
          request->requester->expr_flags-=OEXPR_LEFT_DBL;
          request->requester->expr_flags|=OEXPR_LEFT_INT;
          
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_CONTENTS;
          reg_of_count->flags|=COUNT_CONTENTS;
          if (!is_enclosed) request->requester->left=(void*)&(count->data.move.n_at);
          else request->requester->left=(void*)&(count->data.move.n_enclosed);
          break;
        case REPORT_RXNS:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_RXNS;
          reg_of_count->flags|=COUNT_RXNS;
          if (!is_enclosed) request->requester->left=(void*)&(count->data.rx.n_rxn_at);
          else request->requester->left=(void*)&(count->data.rx.n_rxn_enclosed);
          break;
        case REPORT_FRONT_HITS:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          reg_of_count->flags|=COUNT_HITS;
          request->requester->left=(void*)&(count->data.move.front_hits);
          break;
        case REPORT_BACK_HITS:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          reg_of_count->flags|=COUNT_HITS;
          request->requester->left=(void*)&(count->data.move.back_hits);
          break;
        case REPORT_FRONT_CROSSINGS:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          reg_of_count->flags|=COUNT_HITS;
          request->requester->left=(void*)&(count->data.move.front_to_back);
          break;
        case REPORT_BACK_CROSSINGS:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          reg_of_count->flags|=COUNT_HITS;
          request->requester->left=(void*)&(count->data.move.back_to_front);
          break;
        case REPORT_ALL_HITS:
          request->requester->expr_flags|=OEXPR_RIGHT_DBL;
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          reg_of_count->flags|=COUNT_HITS;
          request->requester->left=(void*)&(count->data.move.front_hits);
          request->requester->right=(void*)&(count->data.move.back_hits);
          break;
        case REPORT_ALL_CROSSINGS:
          request->requester->expr_flags|=OEXPR_RIGHT_DBL;
          reg_of_count->flags|=COUNT_HITS;
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          request->requester->left=(void*)&(count->data.move.front_to_back);
          request->requester->right=(void*)&(count->data.move.back_to_front);
          break;
        case REPORT_CONCENTRATION:
          request->requester->expr_flags|=OEXPR_RIGHT_DBL;
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          reg_of_count->flags|=COUNT_HITS;
          request->requester->left=(void*)&(count->data.move.scaled_hits);
          request->requester->right=(void*)&(world->elapsed_time);
          request->requester->oper='/';
          break;
        default:
          fprintf(world->err_file,"Error at file %s line %d\n  Bad report type %d when creating counts\n",__FILE__,__LINE__,report_type_only);
          return 1;
          break;
      }
    }
  }
  
  return 0;
}


struct counter* create_new_counter(struct region *where,void *who,byte what)
{
  struct counter *c;
  
  c = (struct counter*)mem_get(world->counter_mem);
  if (c==NULL) return NULL;
  
  c->next=NULL;
  c->reg_type=where;
  c->target=who;
  c->counter_type=what;
  if (what&TRIG_COUNTER)
  {
    c->data.trig.t_event=0.0;
    c->data.trig.loc.x = c->data.trig.loc.y = c->data.trig.loc.z = 0.0;
    c->data.trig.listeners=NULL;
  }
  else if (what&RXN_COUNTER)
  {
    c->data.rx.n_rxn_at = c->data.rx.n_rxn_enclosed = 0.0;
  }
  else if (what&MOL_COUNTER)
  {
    c->data.move.n_at = c->data.move.n_enclosed = 0;
    c->data.move.front_hits = c->data.move.back_hits = 0.0;
    c->data.move.front_to_back = c->data.move.back_to_front = 0.0;
    c->data.move.scaled_hits = 0.0;
  }
  return c;
}


void fire_count_event(struct counter *event,int n,struct grid_molecule *g,struct vector3 *where,byte what)
{
  struct trigger_request *tr;
  byte whatelse=what;
  
  if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_HITS) whatelse = (what-REPORT_FRONT_HITS)|REPORT_ALL_HITS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_BACK_HITS) whatelse = (what-REPORT_BACK_HITS)|REPORT_ALL_HITS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_CROSSINGS) whatelse = (what-REPORT_FRONT_CROSSINGS)|REPORT_ALL_CROSSINGS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_BACK_CROSSINGS) whatelse = (what-REPORT_BACK_CROSSINGS)|REPORT_ALL_CROSSINGS;
  
  for (tr=event->data.trig.listeners ; tr!=NULL ; tr=tr->next)
  {
    if (tr->ear->report_type==what || tr->ear->report_type==whatelse)
    {
      if (where!=NULL) memcpy(&(event->data.trig.loc),where,sizeof(struct vector3));
      else uv2xyz(&(g->s_pos),g->grid->surface,&(event->data.trig.loc)); 
      add_trigger_output(event,tr->ear,n);
    }
  }
}

