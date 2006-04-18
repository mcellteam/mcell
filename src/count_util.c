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

void update_collision_count(struct species *sp,struct region_list *rl,int direction,int crossed, double factor)
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
        if (hit_count->reg_type == rl->reg && hit_count->data.move.mol_type == sp)
        {
#if 0
          if (crossed)
          {
            hit_count->data.move.n_enclosed += direction;

/*
            printf("Counted %s (%x) on %s (%x); %x has n_inside = %.1f (up by %d).\n",
                   sp->sym->name,sp->hashval,rl->reg->sym->name,rl->reg->hashval,
                   (int)hit_count,hit_count->data.move.n_inside,direction);
*/

          }
#endif  
          if (rl->reg->flags & sp->flags & COUNT_HITS)
          {
            if (crossed)
            {
              if (direction==1)
              {
                hit_count->data.move.front_hits++;
                hit_count->data.move.front_to_back++;
              }
              else
              {
                hit_count->data.move.back_hits++;
                hit_count->data.move.back_to_front++;
              }
            }
            else
            {
              if (direction==1) hit_count->data.move.front_hits++;
              else hit_count->data.move.back_hits++;
            }
	    if (rl->reg->area != 0.0)
	    {
	      hit_count->data.move.scaled_hits += factor*hits_to_ccn/rl->reg->area;
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
	      fprintf(stderr, "Out of memory: trying to save intermediate results.\n");
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
	  printf("Didn't quite reach waypoint target, fudging.\n");
	  traveling = 0;
	}
	else
	{
	  printf("Couldn't reach waypoint target.  What's wrong?\n");
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
   Out: No return value.  Appropriate counters are updated.
   Note: This handles all types of molecules, from grid to free.  Grid
         and free are implemented, surface is not.  If the reaction
	 pathname is NULL, the molecule is counted.  Otherwise, the
	 reaction is counted at the location of the molecule.
*************************************************************************/

void count_me_by_region(struct abstract_molecule *me,int n,struct rxn_pathname *rxp)
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
  
  //printf("Counting %x by region (up by %d)!\n",(int)me,n);
  
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
          if (c->reg_type == rl->reg && (c->counter_type&ENCLOSING_COUNTER)==0)
	  {
	    if (rxp==NULL && c->data.move.mol_type == sp) c->data.move.n_at += n;
	    else if (c->data.rx.rxn_type == rxp) c->data.rx.n_rxn_at += n;
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
          if (c->reg_type==rl->reg && (my_wall==NULL || !region_listed(my_wall->counting_regions,rl->reg)))
	  {
	    if ( rxp==NULL && c->data.move.mol_type==sp &&
	         (g==NULL || (c->counter_type&ENCLOSING_COUNTER)!=0) )
	    {
	      c->data.move.n_enclosed += n;
	    }
	    else if ( c->data.rx.rxn_type==rxp && 
	              (c->counter_type&ENCLOSING_COUNTER)!=0 )
	    {
	      c->data.rx.n_rxn_enclosed += n;
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
          if (c->reg_type==rl->reg && (my_wall==NULL || !region_listed(my_wall->counting_regions,rl->reg)))
	  {
	    if ( rxp==NULL && c->data.move.mol_type==sp &&
	         (g==NULL || (c->counter_type&ENCLOSING_COUNTER)!=0) )
	    {
	      c->data.move.n_enclosed -= n;
	    }
	    else if ( c->data.rx.rxn_type==rxp && 
	              (c->counter_type&ENCLOSING_COUNTER)!=0 )
	    {
	      c->data.rx.n_rxn_enclosed -= n;
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
                  if (c->reg_type==rl->reg && (my_wall==NULL || !region_listed(my_wall->counting_regions,rl->reg)))
		  {
		    printf("Counting\n");
		    if ( rxp==NULL && c->data.move.mol_type==sp &&
		         (g==NULL || (c->counter_type&ENCLOSING_COUNTER)!=0) )
		    {
		      if (j==COLLIDE_FRONT)
		      {
			c->data.move.n_enclosed += n;
		      }
		      else if (j==COLLIDE_BACK)
		      {
			c->data.move.n_enclosed -= n;
		      }
		    }
		    else if (c->data.rx.rxn_type==rxp && (c->counter_type&ENCLOSING_COUNTER)!=0)
		    {
		      if (j==COLLIDE_FRONT) c->data.rx.n_rxn_enclosed += n;
		      else if (j==COLLIDE_BACK) c->data.rx.n_rxn_enclosed -=n;
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
check_region_counters:
  In: No arguments.
  Out: 0 if counter statements are correct, 1 otherwise.
  Note: A statement is incorrect if a non-closed manifold region
        tries to count a freely diffusing molecule.
********************************************************************/
int check_region_counters()
{
  FILE *log_file;
  struct counter *cp;
  struct region *rp;
  u_int i;

  log_file=world->log_file;  
  
  for (i=0;i<world->count_hashmask+1;i++)
  {
    for (cp=world->count_hash[i];cp!=NULL;cp=cp->next)
    {
      if ( (cp->counter_type & ENCLOSING_COUNTER) != 0)
      {
        rp=cp->reg_type;
	
	if (rp->manifold_flag==MANIFOLD_UNCHECKED) {
	  if (is_manifold(rp)) {
	    rp->manifold_flag=IS_MANIFOLD;
	  }
	  else {
	    rp->manifold_flag=NOT_MANIFOLD;
	  }
	}
	
	if (rp->manifold_flag==NOT_MANIFOLD) {
	  fprintf(log_file,"MCell: error, cannot count molecules or events inside non-manifold object region: %s\n",rp->sym->name); 
	  return (1);
	}
      }
    }
  }

  return(0);
}

