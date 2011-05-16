/**************************************************************************\
** File: count_util.c                                                     **
**                                                                        **
** Purpose: Handles counting of interesting events                        **
**                                                                        **
** Testing status: untested.                                              **
\**************************************************************************/


#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"
#include "logging.h"
#include "grid_util.h"
#include "mcell_structs.h"
#include "util.h"
#include "wall_util.h"
#include "vol_util.h"
#include "count_util.h"
#include "react_output.h"
#include "macromolecule.h"

extern struct volume *world;

/* Instantiate a request to track a particular quantity */
static int instantiate_request(struct output_request *request);

/* Create a new counter data structure */
static struct counter* create_new_counter(struct region *where,void *who,byte what);

/* Utility to resolve count requests for macromolecule states */
static int macro_convert_output_requests(void);

/* Pare down the region lists, annihilating any regions which appear in both
 * lists.  This code was moved out of count_region_from_scratch so that it can
 * be used in the macromolecules counting code as well.
 */
static void clean_region_lists(struct subvolume *my_sv,
                               struct region_list **p_all_regs,
                               struct region_list **p_all_antiregs);

/* Check if the object corresponding to a particular symbol has been referenced
 * directly or indirectly from one of the INSTANTIATE blocks in the model.
 */
static int is_object_instantiated(struct sym_table *entry);

/* Find the list of regions enclosing a particular point. given a particular
 * starting point and starting region list. */
static int find_enclosing_regions(struct vector3 *loc,struct vector3 *start,
                                  struct region_list** rlp,struct region_list** arlp,
                                  struct mem_helper *rmem);

/*************************************************************************
eps_equals:
   In: two doubles
   Out: 1 if they are equal to within some small tolerance, 0 otherwise
*************************************************************************/
static int eps_equals(double x,double y)
{
  double mag;
  double diff;
  
  mag = fabs(x);
  if (mag < fabs(y))
    mag = fabs(y);
  
  diff = fabs(x - y);
  return diff < EPS_C * (mag + 1.0);
}


/*************************************************************************
dup_region_list:
   In: a list of regions
       memory handler to use for duplicated regions
   Out: The duplicated list of regions, or NULL on a memory allocation
        error.
*************************************************************************/

static struct region_list* dup_region_list(struct region_list *r,struct mem_helper *mh)
{
  struct region_list *nr,*rp,*r0;
  
  if (r==NULL) return NULL;
  
  r0 = rp = NULL;
  while (r!=NULL)
  {
    nr = (struct region_list*) CHECKED_MEM_GET(mh, "region list entry");
    nr->next = NULL;
    nr->reg = r->reg;
    if (rp==NULL) r0 = rp = nr;
    else { rp->next = nr; rp = nr; }
    
    r = r->next;
  }
  
  return r0;
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
count_region_update:
   In: species of thing that hit
       region list for the wall we hit
       direction of impact relative to surface normal for volume molecule, or
         relative to the region border for surface molecule
         (inside out = 1, ouside in = 0)
       whether we crossed or not
       scaling factor for reaction probabilities (for estimating ccn)
       location of the hit (for triggers)
       time of the hit (for triggers)
   Out: Returns none.  
        Appropriate counters are updated, that is,
        hit counters are updated according to which side was hit,
	and crossings counters and counts within enclosed regions are
	updated if the surface was crossed.
*************************************************************************/
void count_region_update(struct species *sp,struct region_list *rl,int direction,int crossed, struct vector3 *loc, double t)
{
  struct counter *hit_count;
  double hits_to_ccn=0;
  int count_hits = 0;

  if ((sp->flags&COUNT_HITS) && ((sp->flags&NOT_FREE)==0))
  {
    count_hits = 1;
    hits_to_ccn = sp->time_step * 2.9432976599069717358e-3 /  /* 1e6*sqrt(MY_PI)/(1e-15*N_AV) */ 
                  (sp->space_step*world->length_unit*world->length_unit*world->length_unit);
  }
  
  hit_count = NULL;  
  for ( ; rl != NULL ; rl = rl->next)
  {
    if (rl->reg->flags & COUNT_SOME_MASK)
    {
      int hash_bin = (rl->reg->hashval + sp->hashval)&world->count_hashmask;
      
      for (hit_count=world->count_hash[hash_bin] ; hit_count!=NULL ; hit_count=hit_count->next)
      {
        if (hit_count->reg_type == rl->reg && hit_count->target == sp)
        {
          if (rl->reg->flags & sp->flags & (COUNT_HITS|COUNT_CONTENTS|COUNT_ENCLOSED))
          {
            if (crossed)
            {
              if (direction==1)
              {
                if (hit_count->counter_type&TRIG_COUNTER)
                {
                  hit_count->data.trig.t_event = (double)world->it_time + t; 
                  hit_count->data.trig.orient = 0; 
		  if (rl->reg->flags&sp->flags&COUNT_HITS)
		  {
                    fire_count_event(hit_count,1,loc,REPORT_FRONT_HITS|REPORT_TRIGGER);
                    fire_count_event(hit_count,1,loc,REPORT_FRONT_CROSSINGS|REPORT_TRIGGER);
		  }
		  if (rl->reg->flags&sp->flags&COUNT_CONTENTS)
		  {
                    fire_count_event(hit_count,1,loc,REPORT_ENCLOSED|REPORT_CONTENTS|REPORT_TRIGGER);
		  }
                }
                else
                {
		  if (rl->reg->flags&sp->flags&COUNT_HITS)
		  {
                    hit_count->data.move.front_hits++;
                    hit_count->data.move.front_to_back++;
		  }
		  if (rl->reg->flags&sp->flags&COUNT_CONTENTS)
		  {
		    hit_count->data.move.n_enclosed++;
		  }
                }
              }
              else
              {
                if (hit_count->counter_type&TRIG_COUNTER)
                {
                  hit_count->data.trig.t_event = (double)world->it_time + t; 
                  hit_count->data.trig.orient = 0; 
		  if (rl->reg->flags&sp->flags&COUNT_HITS)
		  {
                    fire_count_event(hit_count,1,loc,REPORT_BACK_HITS|REPORT_TRIGGER);
                    fire_count_event(hit_count,1,loc,REPORT_BACK_CROSSINGS|REPORT_TRIGGER);
		  }
		  if (rl->reg->flags&sp->flags&COUNT_CONTENTS)
		  {
		    fire_count_event(hit_count,-1,loc,REPORT_ENCLOSED|REPORT_CONTENTS|REPORT_TRIGGER);
		  }
                }
                else
                {
		  if (rl->reg->flags&sp->flags&COUNT_HITS)
		  {
                    hit_count->data.move.back_hits++;
                    hit_count->data.move.back_to_front++;
		  }
		  if (rl->reg->flags&sp->flags&COUNT_CONTENTS)
		  {
		    hit_count->data.move.n_enclosed--;
		  }
                }
              }
            }
            else if (rl->reg->flags & sp->flags & COUNT_HITS) /* Didn't cross, only hits might update */
            {
              if (direction==1)
              {
                if (hit_count->counter_type&TRIG_COUNTER)
                {
                  hit_count->data.trig.t_event = (double)world->it_time + t; 
                  hit_count->data.trig.orient = 0; 
                  fire_count_event(hit_count,1,loc,REPORT_FRONT_HITS|REPORT_TRIGGER);
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
                  hit_count->data.trig.t_event = (double)world->it_time + t; 
                  hit_count->data.trig.orient = 0; 
                  fire_count_event(hit_count,1,loc,REPORT_BACK_HITS|REPORT_TRIGGER);
                }
                else hit_count->data.move.back_hits++;
              }
            }
	    if ((count_hits  &&  rl->reg->area != 0.0) && ((sp->flags & NOT_FREE) == 0))
	    {
	      if ((hit_count->counter_type&TRIG_COUNTER)==0)
              {
                hit_count->data.move.scaled_hits += hits_to_ccn/rl->reg->area;
              }
	    }
          }
        }
      }
    }
  }
}


/**************************************************************************
count_region_border_update:
  In: species of thing that hit
      information about the hit (linked list of "hit_data")
  Out: Returns none.  
       Appropriate counters are updated, that is,
       hit counters are updated according to which side was hit,
       and crossings counters  are updated if the region border was crossed.
       We consider FRONT_HITS/FRONT_CROSSINGS when the molecule
       crosses region border "inside out" and BACK_HITS/BACK_CROSSINGS
       when the molecule hits region border "outside in".
**************************************************************************/
void count_region_border_update(struct species *sp,struct hit_data *hd_info)
{
  struct counter *hit_count;
  struct region_list *rl;
  struct hit_data * hd;
  int correct_orient; /* flag*/
 
  /* This function should be used for surface molecules only */
  if((sp->flags & NOT_FREE) == 0) mcell_internal_error("Function 'count_region_border_update()' is used with volume molecules.");
 
  hit_count = NULL;  
  
  for(hd = hd_info; hd != NULL; hd = hd->next)
  {
    for (rl = hd->count_regions ; rl != NULL ; rl = rl->next)
    {
      if (rl->reg->flags & COUNT_SOME_MASK)
      {
        int hash_bin = (rl->reg->hashval + sp->hashval)&world->count_hashmask;
      
        for (hit_count=world->count_hash[hash_bin] ; hit_count!=NULL ; hit_count=hit_count->next)
        {
          correct_orient = 0;
          if((hit_count->orientation == ORIENT_NOT_SET) || (hit_count->orientation == hd->orientation) || (hit_count->orientation == 0)) correct_orient = 1;

          if ((hit_count->reg_type == rl->reg) && (hit_count->target == sp) && correct_orient)
          {
            if (rl->reg->flags & sp->flags & COUNT_HITS)
            {
              if (hd->crossed)
              {
                if (hd->direction==1)
                {
                  if (hit_count->counter_type&TRIG_COUNTER)
                  {
                    hit_count->data.trig.t_event = hd->t; 
                    hit_count->data.trig.orient = 0; 
		    if (rl->reg->flags&sp->flags&COUNT_HITS)
		    {
                      fire_count_event(hit_count,1,&(hd->loc),REPORT_FRONT_HITS|REPORT_TRIGGER);
                      fire_count_event(hit_count,1,&(hd->loc),REPORT_FRONT_CROSSINGS|REPORT_TRIGGER);
		    }
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
                    hit_count->data.trig.t_event = hd->t; 
                    hit_count->data.trig.orient = 0; 
		    if (rl->reg->flags&sp->flags&COUNT_HITS)
		    {
                      fire_count_event(hit_count,1,&(hd->loc),REPORT_BACK_HITS|REPORT_TRIGGER);
                      fire_count_event(hit_count,1,&(hd->loc),REPORT_BACK_CROSSINGS|REPORT_TRIGGER);
		    }
                  }
                  else
                  {
                     hit_count->data.move.back_hits++;
                     hit_count->data.move.back_to_front++;
                  }
                }
              }
              else /* Didn't cross, only hits might update */
              {
                if (hd->direction==1)
                {
                  if (hit_count->counter_type&TRIG_COUNTER)
                  {
                    hit_count->data.trig.t_event = hd->t; 
                    hit_count->data.trig.orient = 0; 
                    fire_count_event(hit_count,1,&(hd->loc),REPORT_FRONT_HITS|REPORT_TRIGGER);
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
                    hit_count->data.trig.t_event = hd->t; 
                    hit_count->data.trig.orient = 0; 
                    fire_count_event(hit_count,1,&(hd->loc),REPORT_BACK_HITS|REPORT_TRIGGER);
                  }
                  else hit_count->data.move.back_hits++;
              
                }
              }
            }
          }
        }
      }
    }
  }  /* end for (hd...) */

}

/*************************************************************************
count_region_from_scratch:
   In: molecule to count, or NULL
       reaction pathname to count, or NULL
       number of these to count
       location at which to count them (may be NULL)
       wall at which this happened (may be NULL)
       time of the hit (for triggers)
   Out: Returns zero on success and 1 on failure.  
        Appropriate counters are updated and triggers are fired.
   Note: At least one of molecule or rxn pathname must be non-NULL; if
        other inputs are NULL, sensible values will be guessed (which
        may themselves be NULL).  This routine is not super-fast for
        volume counts (enclosed counts) since it has to dynamically create
        and test lists of enclosing regions.
*************************************************************************/

void count_region_from_scratch(struct abstract_molecule *am,struct rxn_pathname *rxpn,int n,struct vector3 *loc,struct wall *my_wall,double t)
{  
  struct region_list *rl,*arl,*nrl,*narl; /*a=anti p=previous n=new*/
  struct region_list *all_regs,*all_antiregs;
  struct wall_list *wl;
  struct waypoint *wp;
  struct subvolume *sv,*my_sv;
  struct counter *c;
  void *target;                   /* what we're counting: am->properties or rxpn */
  int hashval;                    /* Hash value of what we're counting */
  double t_hit,t_sv_hit;
  struct vector3 here,delta,hit;  /* For raytracing */
  struct vector3 xyz_loc;         /* Computed location of mol if loc==NULL */
  byte count_flags;
  int pos_or_neg;                 /* Sign of count (neg for antiregions) */
  int orient = SHRT_MIN;          /* orientation of the molecule 
                                     also serves as a flag for 
                                     triggering reactions  */ 
 
  /* Set up values and fill in things the calling function left out */
  if (rxpn!=NULL)
  {
    hashval=rxpn->hashval;
    target=rxpn;
    count_flags=REPORT_RXNS;
  }
  else
  {
    hashval=am->properties->hashval;
    target=am->properties;
    count_flags=REPORT_CONTENTS;
    if (loc==NULL)
    {
      if (am->properties->flags&ON_GRID)
      {
        uv2xyz(&(((struct grid_molecule*)am)->s_pos),((struct grid_molecule*)am)->grid->surface,&xyz_loc);
        loc=&xyz_loc;
      }
      else loc=&(((struct volume_molecule*)am)->pos);
    }
    
    if (my_wall==NULL && (am->properties->flags&ON_GRID)!=0)
    {
      my_wall=((struct grid_molecule*)am)->grid->surface;
    }

    if (am->properties->flags&ON_GRID)
    {
        orient = ((struct grid_molecule *)am)->orient;
    }else{
        orient = 0;
    } 
  }
  
  /* Count grid molecules and reactions on surfaces--easy */
  if (my_wall!=NULL && (my_wall->flags&COUNT_CONTENTS)!=0)
  {
    for (rl=my_wall->counting_regions ; rl!=NULL ; rl=rl->next)
    {
      int hash_bin = (hashval+rl->reg->hashval) & world->count_hashmask;
      for (c=world->count_hash[hash_bin] ; c!=NULL ; c=c->next)
      {
	if (c->target==target && c->reg_type==rl->reg && (c->counter_type&ENCLOSING_COUNTER)==0)
	{
	  if (c->counter_type&TRIG_COUNTER)
	  { 
	    c->data.trig.t_event=t;
	    c->data.trig.orient = orient;
	    fire_count_event(c,n,loc,count_flags|REPORT_TRIGGER);
	  }
	  else if (rxpn==NULL) 
          {
             if (am->properties->flags&ON_GRID)
             {
                if((c->orientation == ORIENT_NOT_SET) || (c->orientation == orient) || (c->orientation == 0))
                {  
                  c->data.move.n_at+=n;
                }
             }
             else
             {
               c->data.move.n_at+=n;
             }  
          }
	  else
            c->data.rx.n_rxn_at+=n;
	}
      }
    }
  }
  
  /* Waypoints must have been placed in order for the following code to work. */
  if (! world->place_waypoints_flag)
    return;

  /* Count volume molecules, vol reactions, and surface stuff that is enclosed--hard!!*/
  if (am==NULL || (am->properties->flags&COUNT_ENCLOSED)!=0 || (am->properties->flags&NOT_FREE)==0)
  {
    const int px = bisect(world->x_partitions,world->nx_parts,loc->x);
    const int py = bisect(world->y_partitions,world->ny_parts,loc->y);
    const int pz = bisect(world->z_partitions,world->nz_parts,loc->z);
    const int this_sv = pz + (world->nz_parts-1)*( py + (world->ny_parts-1)*px );
    wp = &(world->waypoints[this_sv]);
    my_sv = &(world->subvol[this_sv]);
    
    here.x = wp->loc.x;
    here.y = wp->loc.y;
    here.z = wp->loc.z;
    
    all_regs=NULL;
    all_antiregs=NULL;
    
    /* Copy all the potentially relevant regions from the nearest waypoint */
    for ( rl=wp->regions ; rl!=NULL ; rl=rl->next)
    {
      if (rl->reg == NULL) continue;
      int hash_bin = (hashval+rl->reg->hashval) & world->count_hashmask;
      if (world->count_hash[hash_bin]==NULL) continue; /* Won't count on this region so ignore it */
      
      nrl = (struct region_list*) CHECKED_MEM_GET(my_sv->local_storage->regl, "list of enclosing regions for count");
      nrl->reg=rl->reg;
      nrl->next=all_regs;
      all_regs=nrl;
    }
    
    /* And all the antiregions (regions crossed from inside to outside only) */
    for ( arl=wp->antiregions ; arl!=NULL ; arl=arl->next)
    {
      int hash_bin = (hashval+arl->reg->hashval) & world->count_hashmask;
      if (world->count_hash[hash_bin]==NULL) continue; /* Won't count on this region so ignore it */

      narl = (struct region_list*) CHECKED_MEM_GET(my_sv->local_storage->regl, "list of enclosing regions for count");
      narl->reg=arl->reg;
      narl->next=all_antiregs;
      all_antiregs=narl;    
    }
    
    /* Raytrace across any walls from waypoint to us and add to region lists */
    for ( sv = &(world->subvol[this_sv]) ; sv != NULL ; sv = next_subvol(&here,&delta,sv) )
    {
      delta.x = loc->x - here.x;
      delta.y = loc->y - here.y;
      delta.z = loc->z - here.z;
      
      t_sv_hit = collide_sv_time(&here,&delta,sv);
      if (t_sv_hit > 1.0) t_sv_hit = 1.0;
  
      for (wl = sv->wall_head ; wl != NULL ; wl = wl->next)
      {
        /* Skip wall that we are on unless we're a volume molecule */
        if (my_wall==wl->this_wall && (am==NULL || (am->properties->flags&NOT_FREE)))
        {
          continue;
        }
            
	if (wl->this_wall->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
	{
	  int hit_code = collide_wall(&here,&delta,wl->this_wall,&t_hit,&hit,0);
          if (hit_code != COLLIDE_MISS) world->ray_polygon_colls++;
 
	  if (hit_code != COLLIDE_MISS && t_hit <= t_sv_hit &&
	    (hit.x-loc->x)*delta.x + (hit.y-loc->y)*delta.y + (hit.z-loc->z)*delta.z < 0)
	  {
	    for (rl=wl->this_wall->counting_regions ; rl!=NULL ; rl=rl->next)
	    {
	      if ( (rl->reg->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0 )
	      {
		int hash_bin = (hashval+rl->reg->hashval) & world->count_hashmask;
		if (world->count_hash[hash_bin]==NULL) continue; /* Won't count on this region so ignore it */
		
                nrl = (struct region_list*) CHECKED_MEM_GET(my_sv->local_storage->regl, "list of enclosing regions for count");
		nrl->reg = rl->reg;
		if (hit_code == COLLIDE_FRONT)
		{
		  nrl->next=all_regs;
		  all_regs=nrl;
		}
		else if (hit_code == COLLIDE_BACK)
		{
		  nrl->next=all_antiregs;
		  all_antiregs=nrl;
		}
	      }
	    }
	  }
	}
      }
    }

    /* Clean up region lists */
    if (all_regs!=NULL && all_antiregs!=NULL)
      clean_region_lists(my_sv, &all_regs, &all_antiregs);

    /* Actually check the regions here */
    count_flags|=REPORT_ENCLOSED;
    
    for (nrl=all_regs ; nrl!=NULL ; nrl=(nrl==all_regs)?all_antiregs:NULL) /* Trick so we don't need to duplicate this code */
    {
      if (nrl==all_regs) pos_or_neg=1;
      else pos_or_neg=-1;
      for (rl=nrl ; rl!=NULL ; rl=rl->next)
      {
	int hash_bin = (hashval+rl->reg->hashval) & world->count_hashmask;
	for (c=world->count_hash[hash_bin] ; c!=NULL ; c=c->next)
	{
	  if ( c->target==target && c->reg_type==rl->reg &&
	       ((c->counter_type&ENCLOSING_COUNTER)!=0 || (am!=NULL && (am->properties->flags&ON_GRID)==0)) &&
	       (my_wall==NULL || 
	        (am!=NULL && (am->properties->flags&NOT_FREE)==0) ||
	        !region_listed(my_wall->counting_regions,rl->reg)) )
	  {
	    if (c->counter_type&TRIG_COUNTER)
	    {
	      c->data.trig.t_event=t;
	      c->data.trig.orient = orient;
	      fire_count_event(c,n*pos_or_neg,loc,count_flags|REPORT_TRIGGER);
	    }
	    else if (rxpn==NULL) {
               if (am->properties->flags&ON_GRID)
               {
                  if((c->orientation == ORIENT_NOT_SET) || (c->orientation == orient) || (c->orientation == 0)){  
                     c->data.move.n_enclosed += n*pos_or_neg;
                  }
               }else{
                     c->data.move.n_enclosed += n*pos_or_neg;
               }
            }
	    else c->data.rx.n_rxn_enclosed += n*pos_or_neg;
	  }
	}
      }
    }
    
    /* Free region memory */ 
    if (all_regs!=NULL) mem_put_list(my_sv->local_storage->regl,all_regs);
    if (all_antiregs!=NULL) mem_put_list(my_sv->local_storage->regl,all_antiregs);
  }
}


/*************************************************************************
count_moved_grid_mol:
   In: molecule to count
       new grid for molecule
       new location on that grid
   Out: Returns zero on success and 1 on failure.  
        Appropriate counters are updated and triggers are fired.
   Note: This routine is not super-fast for enclosed counts for
         surface molecules since it raytraces without using waypoints.
*************************************************************************/
void count_moved_grid_mol(struct grid_molecule *g, struct surface_grid *sg, struct vector2 *loc)
{
  struct region_list *rl,*prl,*nrl,*pos_regs,*neg_regs;
  struct storage *stor;
  struct counter *c;
  struct vector3 origin;
  struct vector3 target;
  struct vector3 *where = NULL;
  int delete_me;
  int n;
  int origin_loaded=0;
  int target_loaded=0;
  
  pos_regs = neg_regs = NULL;
  stor = g->grid->surface->birthplace;
 
  
  if (g->grid != sg) /* Different grids implies different walls, so we might have changed regions */
  {
    delete_me=0;
    if ((g->grid->surface->flags&COUNT_CONTENTS)!=0 && (sg->surface->flags&COUNT_CONTENTS)!=0)
    {
      delete_me=1;
      nrl = g->grid->surface->counting_regions;
      prl = sg->surface->counting_regions;
      while (prl!=NULL && nrl!=NULL)
      {
        if (prl->reg == nrl->reg) /* Skip identical regions */
        {
          prl=prl->next;
          nrl=nrl->next;
          continue;
        }
        while (prl!=NULL && prl->reg < nrl->reg) /* Entering these regions */
        {
          rl = (struct region_list*) CHECKED_MEM_GET(stor->regl, "region list entry");
          rl->next=pos_regs;
          rl->reg=prl->reg;
          pos_regs=rl;
          prl=prl->next;
        }
        while (nrl!=NULL && (prl==NULL || nrl->reg < prl->reg)) /* Leaving these regions */
        {
          rl = (struct region_list*) CHECKED_MEM_GET(stor->regl, "region list entry");
          rl->next=neg_regs;
          rl->reg=nrl->reg;
          neg_regs=rl;
          nrl=nrl->next;
        }
      }

      /* If we exhaust all negative regions before we've exhausted all
       * positive regions, the above loop will terminate, leaving some
       * regions uncounted. */
      while (prl != NULL)
      {
        rl = (struct region_list*) CHECKED_MEM_GET(stor->regl, "region list entry");
        rl->next=pos_regs;
        rl->reg=prl->reg;
        pos_regs=rl;
        prl=prl->next;
      }

      /* I don't think this can happen, but it could potentially happen
       * if, say, prl started off NULL (i.e. one of the grids belonged
       * to no counting regions at all). */
      while (nrl!=NULL)
      {
        rl = (struct region_list*) CHECKED_MEM_GET(stor->regl, "region list entry");
        rl->next=neg_regs;
        rl->reg=nrl->reg;
        neg_regs=rl;
        nrl=nrl->next;
      }
    }
    else if (g->grid->surface->flags&COUNT_CONTENTS) neg_regs = g->grid->surface->counting_regions;
    else if (sg->surface->flags&COUNT_CONTENTS) pos_regs = sg->surface->counting_regions;
    
    /* Sneaky way to go through both lists in one loop */
    n=1;
    if (pos_regs!=NULL)
    {
      uv2xyz(loc,sg->surface,&target);
      where = &target;
      target_loaded=1;
    }
    for (rl = (pos_regs != NULL) ? pos_regs : neg_regs;
         rl != NULL;
         rl = (rl->next == NULL && n>0) ? neg_regs : rl->next)
    {
      if (rl==neg_regs)
      {
        uv2xyz(&(g->s_pos),g->grid->surface,&origin);
        where = &origin;
        origin_loaded=1;
        n=-1;
      }

      int hash_bin = (g->properties->hashval+rl->reg->hashval) & world->count_hashmask;
      for (c=world->count_hash[hash_bin] ; c!=NULL ; c=c->next)
      {
        if (c->target==g->properties && c->reg_type==rl->reg && (c->counter_type&ENCLOSING_COUNTER)==0)
        {
          if (c->counter_type&TRIG_COUNTER)
          {
            c->data.trig.t_event = g->t;
            c->data.trig.orient = g->orient;
            fire_count_event(c,n,where,REPORT_CONTENTS|REPORT_TRIGGER);
          }
          else if((c->orientation == ORIENT_NOT_SET) || (c->orientation == g->orient) || (c->orientation == 0))
          {  
            c->data.move.n_at += n;
          }
        }
      }
    }
    
    
    if (delete_me)
    {
      if (pos_regs!=NULL) mem_put_list(stor->regl,pos_regs);
      if (neg_regs!=NULL) mem_put_list(stor->regl,neg_regs);
    }
  }

  if (g->properties->flags&COUNT_ENCLOSED) /* Have to raytrace */
  {
    struct vector3 delta;
    struct vector3 here;
    struct vector3 hit;
    struct subvolume *sv;
    struct wall_list *wl;
    double t_sv_hit,t;
    int j;
    
    pos_regs=neg_regs=NULL;
    
    if (!origin_loaded) uv2xyz(&(g->s_pos),g->grid->surface,&origin);
    if (!target_loaded) uv2xyz(loc,sg->surface,&target);
    delta.x = target.x-origin.x;
    delta.y = target.y-origin.y;
    delta.z = target.z-origin.z;
    
    here=origin;
    
    /* Collect all the relevant regions we pass through */
    for ( sv=find_subvolume(&origin,NULL) ; sv!=NULL ; sv = next_subvol(&here,&delta,sv) )
    {
      t_sv_hit = collide_sv_time(&here,&delta,sv);
      if (t_sv_hit>1.0) t_sv_hit=1.0;
      
      for (wl=sv->wall_head ; wl!=NULL ; wl=wl->next)
      {
        if (wl->this_wall==g->grid->surface || wl->this_wall==sg->surface) continue;  /* Don't count our own wall */
        
        j = collide_wall(&here,&delta,wl->this_wall,&t,&hit,0);
        
        if (j!=COLLIDE_MISS) world->ray_polygon_colls++;
        
        if (j!=COLLIDE_MISS && t<t_sv_hit && (hit.x-target.x)*delta.x + (hit.y-target.y)*delta.y + (hit.z-target.z)*delta.z < 0)
        {
          for (rl=wl->this_wall->counting_regions ; rl!=NULL ; rl=rl->next)
          {
            if ((rl->reg->flags&COUNT_ENCLOSED)==0) continue;  /* Only ENCLOSED counted here */
            
            if (j==COLLIDE_FRONT)
            {
              prl = (struct region_list*) CHECKED_MEM_GET(stor->regl, "region list entry");
              prl->reg = rl->reg;
              prl->next=pos_regs;
              pos_regs=prl;
            }
            else if (j==COLLIDE_BACK)
            {
              nrl = (struct region_list*) CHECKED_MEM_GET(stor->regl, "region list entry");
              nrl->reg = rl->reg;
              nrl->next=neg_regs;
              neg_regs=nrl;
            }
          }
        }
      }
    }
    
    if (pos_regs != NULL)
    {
      pos_regs = (struct region_list*)void_list_sort((struct void_list*)pos_regs);
    }
    if (neg_regs != NULL)
    {
      neg_regs = (struct region_list*)void_list_sort((struct void_list*)neg_regs);
    }
    
    prl = pos_regs;
    nrl = neg_regs;
    while (prl!=NULL || nrl!=NULL)
    {
      if (prl==NULL)
      {
        rl=nrl;
        nrl=nrl->next;
        n=-1;
        where=&origin;
      }
      else if (nrl==NULL)
      {
        rl=prl;
        prl=prl->next;
        n=1;
        where=&target;
      }
      else if (prl->reg < nrl->reg)
      {
        rl = prl;
        prl=prl->next;
        n=1;
        where=&origin;
      }
      else if (nrl->reg < prl->reg)
      {
        rl = nrl;
        nrl=nrl->next;
        n=-1;
        where=&target;
      }
      else
      {
        n = 1;              /* dummy init to silence compiler */
        rl = NULL;
        prl = prl->next;
        nrl = nrl->next;
      }
      
      if (rl!=NULL)
      {
        int hash_bin = (g->properties->hashval+rl->reg->hashval) & world->count_hashmask;
        for (c=world->count_hash[hash_bin] ; c!=NULL ; c=c->next)
        {
          if (c->target==g->properties && c->reg_type==rl->reg && (c->counter_type&ENCLOSING_COUNTER)!=0 &&
              !region_listed(g->grid->surface->counting_regions,rl->reg) && !region_listed(sg->surface->counting_regions,rl->reg))
          {
            if (c->counter_type&TRIG_COUNTER)
            {
              c->data.trig.t_event = g->t;
              c->data.trig.orient = g->orient;
              fire_count_event(c,n,where,REPORT_CONTENTS|REPORT_ENCLOSED|REPORT_TRIGGER);
            }
            else if((c->orientation == ORIENT_NOT_SET) || (c->orientation == g->orient) || (c->orientation == 0)){  
              c->data.move.n_enclosed += n;
            }
          }
        }
      }
    }
    
    if (pos_regs!=NULL) mem_put_list(stor->regl,pos_regs);
    if (neg_regs!=NULL) mem_put_list(stor->regl,neg_regs);
  }
}



/*************************************************************************
fire_count_event:
   In: counter of thing that just happened (trigger of some sort)
       number of times that thing happened (or hit direction for triggers)
       location where it happened
       what happened (Report Type Flags)   
   Out: None
*************************************************************************/

void fire_count_event(struct counter *event,int n,struct vector3 *where,byte what)
{
  struct trigger_request *tr;
  byte whatelse=what;
  short flags;

  if ((what&REPORT_TYPE_MASK)==REPORT_RXNS) flags = TRIG_IS_RXN;
  else if ((what&REPORT_TYPE_MASK)==REPORT_CONTENTS) flags = 0;
  else flags = TRIG_IS_HIT;
 
  if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_HITS) whatelse = (what-REPORT_FRONT_HITS)|REPORT_ALL_HITS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_BACK_HITS) whatelse = (what-REPORT_BACK_HITS)|REPORT_ALL_HITS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_CROSSINGS) whatelse = (what-REPORT_FRONT_CROSSINGS)|REPORT_ALL_CROSSINGS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_BACK_CROSSINGS) whatelse = (what-REPORT_BACK_CROSSINGS)|REPORT_ALL_CROSSINGS;
 
  for (tr=event->data.trig.listeners ; tr!=NULL ; tr=tr->next)
  {
    if (tr->ear->report_type==what)
    {
      memcpy(&(event->data.trig.loc),where,sizeof(struct vector3));
      if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_HITS || (what&REPORT_TYPE_MASK)==REPORT_FRONT_CROSSINGS){
          add_trigger_output(event,tr->ear,n,flags);
      }else if ((what&REPORT_TYPE_MASK)==REPORT_BACK_HITS || (what&REPORT_TYPE_MASK)==REPORT_BACK_CROSSINGS){
          add_trigger_output(event,tr->ear,-n,flags);
      }else{
          add_trigger_output(event,tr->ear,n,flags);
      }
      
    }
    else if (tr->ear->report_type==whatelse)
    {
      memcpy(&(event->data.trig.loc),where,sizeof(struct vector3));
      if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_HITS || (what&REPORT_TYPE_MASK)==REPORT_FRONT_CROSSINGS){
        add_trigger_output(event,tr->ear,n,flags);
      }else{
        add_trigger_output(event,tr->ear,-n,flags);
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
static int find_enclosing_regions(struct vector3 *loc,
                                  struct vector3 *start,
                                  struct region_list** rlp,
                                  struct region_list** arlp,
                                  struct mem_helper *rmem)
{
  struct vector3 outside,delta,hit;
  struct subvolume *sv,*svt;
  struct wall_list *wl;
  struct region_list *rl,*arl;
  struct region_list *trl,*tarl,*xrl,*yrl,*nrl;
  double t,t_hit_sv;
  int traveling;
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
      int hit_code = collide_wall(&outside , &delta , wl->this_wall , &t , &hit , 0);
      
      if((hit_code != COLLIDE_MISS) && (world->notify->final_summary == NOTIFY_FULL)){
          world->ray_polygon_colls++;
      }	 

      if (hit_code==COLLIDE_REDO)
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
      else if (hit_code==COLLIDE_MISS || !(t >= 0 && t < 1.0) || t > t_hit_sv || (wl->this_wall->flags & (COUNT_CONTENTS|COUNT_RXNS|COUNT_ENCLOSED)) == 0 ||
	       (hit.x-outside.x)*delta.x + (hit.y-outside.y)*delta.y + (hit.z-outside.z)*delta.z < 0) continue;
      else
      {
        for (xrl=wl->this_wall->counting_regions ; xrl != NULL ; xrl = xrl->next)
        {
          if ((xrl->reg->flags & (COUNT_CONTENTS|COUNT_RXNS|COUNT_ENCLOSED)) != 0)
          {
            nrl = (struct region_list*) CHECKED_MEM_GET(rmem, "region list entry");
            nrl->reg = xrl->reg;
            
            if (hit_code==COLLIDE_BACK) { nrl->next = tarl; tarl = nrl; }
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
	  mcell_log("Didn't quite reach waypoint target, fudging.");
	  traveling = 0;
	}
	else
	{
	  mcell_log("Couldn't reach waypoint target.");
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
place_waypoints:
   In: No arguments.
   Out: Returns 1 if malloc fails, 0 otherwise.
        Allocates waypoints to SSVs, if any are needed.
   Note: you must have initialized SSVs before calling this routine!
*************************************************************************/
int place_waypoints(void)
{
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

  /* Probably ought to check for whether you really need waypoints */
  
  if (world->waypoints != NULL) free(world->waypoints);
  world->n_waypoints = world->n_subvols;
  world->waypoints = CHECKED_MALLOC_ARRAY(struct waypoint,
                                          world->n_waypoints,
                                          "waypoints");
  memset(world->waypoints,0,world->n_waypoints*sizeof(struct waypoint));

  for (int px=0;px<world->nx_parts-1;px++)
  {
    for (int py=0;py<world->ny_parts-1;py++)
    {
      for (int pz=0;pz<world->nz_parts-1;pz++)
      {
        const int this_sv = pz + (world->nz_parts-1)*(py + (world->ny_parts-1)*px);
        wp = &(world->waypoints[this_sv]);
        
        sv = &(world->subvol[this_sv]);
        
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
        
        if (pz>0)
        {
	  if (world->waypoints[this_sv-1].regions != NULL)
	  {
            wp->regions = dup_region_list(world->waypoints[this_sv-1].regions,sv->local_storage->regl);
            if (wp->regions == NULL) return 1;
	  }
	  else wp->regions = NULL;
	  
	  if (world->waypoints[this_sv-1].antiregions != NULL)
	  {
            wp->antiregions = dup_region_list(world->waypoints[this_sv-1].antiregions,sv->local_storage->regl);
            if (wp->antiregions == NULL) return 1;
	  }
	  else wp->antiregions = NULL;
          
          if (find_enclosing_regions(&(wp->loc),&(world->waypoints[this_sv-1].loc),
                                     &(wp->regions),&(wp->antiregions),sv->local_storage->regl))
            return 1;
        }
        else
        {
          wp->regions = NULL;
          wp->antiregions = NULL;
          if (find_enclosing_regions(&(wp->loc),NULL,&(wp->regions),
                                     &(wp->antiregions),sv->local_storage->regl))
            return 1;
        }
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

/******************************************************************
prepare_counters:
  In: No arguments.
  Out: 0 if counter statements are correct, 1 otherwise.
  Note: A statement is incorrect if a non-closed manifold region
        tries to count a freely diffusing molecule.  Fixes up all
        count requests to point at the data we care about.
********************************************************************/
int prepare_counters(void)
{
  /* First give everything a sensible name, if needed */
  for (struct output_block *block=world->output_block_head;
       block != NULL;
       block = block->next)
  {
    for (struct output_set *set = block->data_set_head;
         set != NULL;
         set = set->next)
    {
      if (set->header_comment==NULL) continue;
      for (struct output_column *column = set->column_head;
           column != NULL;
           column = column->next)
      {
        if (column->expr->title==NULL) column->expr->title = oexpr_title(column->expr);
        if (column->expr->title==NULL)
          mcell_allocfailed("Unable to create title for reaction data output.");
      }
    }
  }

  /* Then convert all requests to real counts */
  for (struct output_request *request=world->output_request_head;
       request != NULL;
       request = request->next)
  {
     /* check whether the "count_location" refers to the instantiated
        object or region */ 
    if (request->count_location != NULL )
    {
      if (! is_object_instantiated(request->count_location))
        mcell_error("The object/region name '%s' in the COUNT/TRIGGER statement is not fully referenced.\n"
                    "  This occurs when a count is requested on an object which has not been referenced\n"
                    "  (directly or indirectly) from an INSTANTIATE block in the MDL file.",
                    request->count_location->name);
    }

    if (request->count_target->sym_type == MOL)
    {
      struct species *sp = (struct species *)(request->count_target->value);

      /* For volume molecules: */
      if ((sp->flags & ON_GRID) == 0)
      {
        /* Make sure orientation is not set */
        if (request->count_orientation != ORIENT_NOT_SET)
        {
          switch (world->notify->useless_vol_orient)
          {
            case WARN_WARN:
              mcell_warn("An orientation has been given for the molecule '%s', which is a volume molecule.\n"
                         "  Orientation is valid only for grid molecules, and will be ignored.",
                         request->count_target->name);
              /* Fall through */
            case WARN_COPE:
              request->count_orientation = ORIENT_NOT_SET;
              break;


            case WARN_ERROR:
              mcell_error("An orientation has been given for the molecule '%s', which is a volume molecule.\n"
                          "  Orientation is valid only for grid molecules.",
                          request->count_target->name);
              break;
          }
        }
      }

    }
 
    if (request->count_location!=NULL && request->count_location->sym_type==OBJ)
    {
      if (expand_object_output(request, (struct object*) request->count_location->value))
        mcell_error("Failed to expand request to count on object.");
    }
    
    if (instantiate_request(request))
      mcell_error("Failed to instantiate count request.");
  }

  /* Need to keep all the requests for now...could repackage them to save memory */
  macro_convert_output_requests();
  
  return 0;
}

/******************************************************************
is_object_instantiated:
  In: object
      symbol_table entry against which the object is tested 
  Out: 1 if the name of the object or one of its descendants matches the name 
       of the symbol passed, 0 otherwise.
  Note: Checking is performed for all instantiated objects
********************************************************************/
static int is_object_instantiated(struct sym_table *entry)
{
  struct object *obj = NULL;
  if (entry->sym_type == REG)
    obj = ((struct region *) (entry->value))->parent;
  else if (entry->sym_type == OBJ)
    obj = ((struct object *) (entry->value));
  else
    return 0;

  for (; obj != NULL; obj = obj->parent)
  {
    if (obj == world->root_instance)
      return 1;
  }

  return 0;
}

/*************************************************************************
check_counter_geometry:
   In: nothing  
   Out: 0 on success, 1 on failure.
        Checks all counters to make sure that if they are ENCLOSING,
        they count on closed regions.  If not, the function prints out
        the offending region name and returns 1.
*************************************************************************/
int check_counter_geometry(void)
{
  /* Check to make sure what we've created is geometrically sensible */
  for (int i = 0; i < world->count_hashmask+1; i++)
  {
    for (struct counter *cp = world->count_hash[i]; cp != NULL; cp = cp->next)
    {
      if ( (cp->counter_type & ENCLOSING_COUNTER) != 0)
      {
        struct region *rp = cp->reg_type;

	if (rp->manifold_flag==MANIFOLD_UNCHECKED)
        {
	  if (is_manifold(rp)) rp->manifold_flag=IS_MANIFOLD;
	  else rp->manifold_flag=NOT_MANIFOLD;
	}
	
	if (rp->manifold_flag==NOT_MANIFOLD)
	  mcell_error("Cannot count molecules or events inside non-manifold object region '%s'.  Please make sure that all objects/regions used to count 3D molecules are closed/watertight.",
                      rp->sym->name); 

	world->place_waypoints_flag=1;
      }
    }
  }

  return 0;
}

/*************************************************************************
expand_object_output:
   In: request for a count
       object upon which the request is made.
   Out: 0 on success, 1 on failure (memory allocation only?).
        Request is split into a separate request for each BOX and POLY
        object's ALL region that is a child of this object.  The result
        is then added up here.
   Note: This is probably broken for concentration.  It may also not be
         the most intuitive interpretation when used inside a large
         object with multiple layers of nesting--if one molecule is
         inside three sub-objects, it will be counted three times!
   PostNote: Checks that COUNT/TRIGGER statements are not allowed for
             metaobjects and release objects.
*************************************************************************/
int expand_object_output(struct output_request *request,struct object *obj)
{
#ifdef ALLOW_COUNTS_ON_METAOBJECT
  int n_expanded;
#endif

  switch (obj->object_type)
  {
    case REL_SITE_OBJ:
      mcell_error("COUNT and TRIGGER statements on release object '%s' are not allowed.\n",obj->sym->name);
      break;

    case META_OBJ:
      /* XXX: Should this really be disabled?  Some comments by Tom lead me to
       * believe that, despite the potential confusion for users, this should
       * not be disabled.
       */
#ifndef ALLOW_COUNTS_ON_METAOBJECT
      mcell_error("COUNT and TRIGGER statements on metaobject '%s' are not allowed.\n",obj->sym->name);
#else
#error "Support for counting in/on a metaobject doesn't work right now."
      n_expanded=0;
      for (struct object *child=obj->first_child ; child!=NULL ; child=child->next)
      {
        if (!object_has_geometry(child)) continue;  /* NOTE -- for objects nested N deep, we check this N+(N-1)+...+2+1 times (slow) */
        if (n_expanded>0)
        {
          struct output_request *new_request = (struct output_request*)mem_get(world->outp_request_mem);
          struct output_expression *oe = request->requester;
          struct output_expression *oel = new_output_expr(world->oexpr_mem);
          struct output_expression *oer = new_output_expr(world->oexpr_mem);
          if (new_request==NULL || oel==NULL || oer==NULL)
            mcell_allocfailed("Failed to expand count expression on object %s.",
                              obj->sym->name);
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
        ++ n_expanded;
      }
      if (n_expanded==0)
        mcell_error("Trying to count on object %s but it has no geometry.",
                    obj->sym->name);
#endif
      break;

    case BOX_OBJ:
    case POLY_OBJ:
      request->count_location = NULL;
      for (struct region_list *rl = obj->regions; rl != NULL; rl = rl->next)
      {
        if (is_reverse_abbrev(",ALL",rl->reg->sym->name))
        {
          request->count_location = rl->reg->sym;
          break;
        }
      }
      if (request->count_location == NULL)
        mcell_internal_error("ALL region missing on object %s", obj->sym->name);
      break;

    default:
      UNHANDLED_CASE(obj->object_type);
      return 1;
  }
  return 0;
}

/*************************************************************************
object_has_geometry:
   In: object (instantiated in world)  
   Out: 0 if there are no geometrical objects within that object (and it
        is not a geometrical object itself).  1 if there are such object.
*************************************************************************/

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

    case REL_SITE_OBJ:
    default:
      return 0;
      break;
  }
  return 0;
}



/*************************************************************************
instantiate_request:
   In: request for a count  
   Out: 0 on success, 1 on failure (memory allocation only?).
        Requesting output tree gets appropriate node pointed to the
        memory location where we will be collecting data.
*************************************************************************/
static int instantiate_request(struct output_request *request)
{
  int request_hash = 0;
  struct rxn_pathname *rxpn_to_count;
  struct rxn *rx_to_count = NULL;
  struct species *mol_to_count = NULL;
  void *to_count;
  struct region *reg_of_count;
  struct counter *count = NULL;
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
      if ((mol_to_count->flags&NOT_FREE)==0 && (request->report_type&REPORT_TYPE_MASK)==REPORT_CONTENTS)
      {
        request->report_type|=REPORT_ENCLOSED;
      }
      request_hash=mol_to_count->hashval;
      break;

    case RXPN:
      rxpn_to_count=(struct rxn_pathname*)to_count;
      rx_to_count=rxpn_to_count->rx;
      mol_to_count=NULL;
      if ((rx_to_count->players[0]->flags&NOT_FREE)==0 &&
	  (rx_to_count->n_reactants==1 || (rx_to_count->players[1]->flags&NOT_FREE)==0))
      {
	request->report_type|=REPORT_ENCLOSED;
      }
      request_hash=rxpn_to_count->hashval;
      break;

    default:
      UNHANDLED_CASE(request->count_target->sym_type);
      return 1;
  }
 
  if (request->count_location!=NULL)
  {
    if (request->count_location->sym_type != REG)
      mcell_internal_error("Non-region location symbol (type=%d) in count request.",
                     request->count_location->sym_type);
    reg_of_count = (struct region*) request->count_location->value;

    request_hash += reg_of_count->hashval;
   
  }
  else reg_of_count=NULL;
  request_hash&=world->count_hashmask;
  
  /* Now create count structs and set output expression to point to data */
  report_type_only=request->report_type&REPORT_TYPE_MASK;
  request->requester->expr_flags-=OEXPR_LEFT_REQUEST;  
  if ((request->report_type&REPORT_TRIGGER)==0 && request->count_location==NULL) /* World count is easy! */
  {
    request->report_type &= ~REPORT_ENCLOSED;
    switch (report_type_only)
    {
      case REPORT_CONTENTS:
        request->requester->expr_flags|=OEXPR_LEFT_INT;
        request->requester->left=(void*)&(mol_to_count->population);
        break;
      case REPORT_RXNS:
        assert(rx_to_count != NULL);
        request->requester->expr_flags|=OEXPR_LEFT_DBL;
        request->requester->left=(void*)&(rx_to_count->info[rxpn_to_count->path_num].count);
        break;
      default:
        mcell_internal_error("Invalid report type 0x%x in count request.", report_type_only);
        return 1;
    }
  }
  else /* Triggered count or count on region */
  {
    /* Set count type flags */
    if (report_type_only==REPORT_RXNS) count_type=RXN_COUNTER;
    else count_type=MOL_COUNTER;
    if (request->report_type&REPORT_ENCLOSED)
    {
      assert(reg_of_count != NULL);
      reg_of_count->flags|=COUNT_ENCLOSED;
      count_type|=ENCLOSING_COUNTER;
      if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_ENCLOSED;
    }
    if (request->report_type&REPORT_TRIGGER)
    {
      assert(reg_of_count != NULL);
      count_type|=TRIG_COUNTER;
      reg_of_count->flags|=COUNT_TRIGGER;
    }
    
    /* Find or add counter */
    for (count=world->count_hash[request_hash] ; count!=NULL ; count=count->next)
    {
      if (count->reg_type==reg_of_count && count->target==to_count && count_type==count->counter_type && count->orientation == request->count_orientation) break;
    }
    if (count==NULL)
    {
      count=create_new_counter(reg_of_count, request->count_target->value, count_type);
      if(request->count_orientation != ORIENT_NOT_SET)
      {
           count->orientation = request->count_orientation;
      }      
      
      count->next=world->count_hash[request_hash];
      world->count_hash[request_hash]=count;
    }
 
    is_enclosed = ((request->report_type&REPORT_ENCLOSED)!=0);
    
    /* Point appropriately */
    if (request->report_type&REPORT_TRIGGER)
    {
      trig_req = (struct trigger_request*) CHECKED_MEM_GET(world->trig_request_mem,
                                                           "trigger notification request");
      trig_req->next=count->data.trig.listeners;
      count->data.trig.listeners=trig_req;
      trig_req->ear=request;
      
      request->requester->expr_flags|=OEXPR_TYPE_TRIG;
      
      if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_TRIGGER;
      assert(reg_of_count != NULL);
      switch (report_type_only)
      {
        case REPORT_CONTENTS:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_CONTENTS;
          reg_of_count->flags|=COUNT_CONTENTS;
          break;
        case REPORT_RXNS:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_RXNS;
          reg_of_count->flags|=COUNT_RXNS;
          break;
        case REPORT_FRONT_HITS:
        case REPORT_BACK_HITS:
        case REPORT_FRONT_CROSSINGS:
        case REPORT_BACK_CROSSINGS:
        case REPORT_ALL_HITS:
        case REPORT_ALL_CROSSINGS:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          reg_of_count->flags|=COUNT_HITS;
          break;
        case REPORT_CONCENTRATION:
          if (mol_to_count!=NULL)
          {
             if(mol_to_count->flags & ON_GRID)
             {
               mcell_error("ESTIMATE_CONC counting on regions is implemented only for volume molecules, while %s is a surface molecule.", mol_to_count->sym->name);
             }else{
                mol_to_count->flags|=COUNT_HITS;
                reg_of_count->flags|=COUNT_HITS;
             }
          }
          break;
        default:
          UNHANDLED_CASE(report_type_only);
          return 1;
      }
    }
    else /* Not trigger--set up for regular count */
    {
      assert(reg_of_count != NULL);
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
          if (mol_to_count!=NULL) 
          {
             if(mol_to_count->flags & ON_GRID)
             {
               mcell_error("ESTIMATE_CONC counting on regions is implemented only for volume molecules, while %s is a surface molecule.", mol_to_count->sym->name);
             }else{
               mol_to_count->flags|=COUNT_HITS;
             }
          }
          reg_of_count->flags|=COUNT_HITS;
          request->requester->left=(void*)&(count->data.move.scaled_hits);
          request->requester->right=(void*)&(world->elapsed_time);
          request->requester->oper='/';
          break;

        default:
          UNHANDLED_CASE(report_type_only);
          return 1;
      }
    }
  }
  
  return 0;
}


/*************************************************************************
create_new_counter:
   In: region upon which to count
       target we're going to count (species or rxn pathname)
       what to count (*_COUNTER flags)   
   Out: Newly allocated counter initialized with the given region and
        target, or NULL if there is a memory allocation error.
   Note: memory is allocated from world->counter_mem using mem_get,
         not from the global heap using malloc.
*************************************************************************/
static struct counter* create_new_counter(struct region *where,void *who,byte what)
{
  struct counter *c;
  
  c = (struct counter*) CHECKED_MEM_GET(world->counter_mem, "counter");
  c->next=NULL;
  c->reg_type=where;
  c->target=who;
  c->orientation = ORIENT_NOT_SET;
  c->counter_type=what;
  if (what&TRIG_COUNTER)
  {
    c->data.trig.t_event=0.0;
    c->data.trig.loc.x = c->data.trig.loc.y = c->data.trig.loc.z = 0.0;
    c->data.trig.orient = SHRT_MIN;
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

/*************************************************************************
clean_region_lists:
   Cleans the region and antiregion lists, annihilating any items which appear
   on both lists.

   In:  struct subvolume *my_sv - subvolume containing waypoint
        struct region_list **p_all_regs - pointer to receive list of regions
        struct region_list **p_all_antiregs - pointer to receive list of antiregions
   Out: None
*************************************************************************/
static void clean_region_lists(struct subvolume *my_sv,
                               struct region_list **p_all_regs,
                               struct region_list **p_all_antiregs)
{
  if ((*p_all_regs)->next!=NULL || (*p_all_antiregs)->next!=NULL)
  {
    struct region_list pre_sentry,pre_antisentry;
    struct region_list *prl, *parl, *rl, *arl;

    /* Sort by memory address to make mutual annihilation faster */
    if ((*p_all_regs)->next!=NULL) *p_all_regs=(struct region_list*)void_list_sort((struct void_list*)*p_all_regs);
    if ((*p_all_antiregs)->next!=NULL) *p_all_antiregs=(struct region_list*)void_list_sort((struct void_list*)*p_all_antiregs);

    /* Need previous entry to fix up list, so we'll make an imaginary one for 1st list element */
    pre_sentry.next=*p_all_regs;
    pre_antisentry.next=*p_all_antiregs;
    prl=&pre_sentry;
    parl=&pre_antisentry;

    /* If we cross a region both ways, throw both out (once) */
    for (rl=*p_all_regs,arl=*p_all_antiregs ; rl!=NULL && arl!=NULL ; prl=rl,rl=rl->next,parl=arl,arl=arl->next)
    {
      if (rl->reg==arl->reg) /* Mutual annihilation */
      {
        prl->next=rl->next;
        parl->next=arl->next;
        mem_put(my_sv->local_storage->regl,rl);
        mem_put(my_sv->local_storage->regl,arl);
        rl=prl;
        arl=parl;
      }
    }
    *p_all_regs=pre_sentry.next;
    *p_all_antiregs=pre_antisentry.next;
  }
  else if ((*p_all_regs)->reg==(*p_all_antiregs)->reg)
  {
    /* Crossed one region both ways, toss them */
    mem_put(my_sv->local_storage->regl,*p_all_regs);
    mem_put(my_sv->local_storage->regl,*p_all_antiregs);
    *p_all_regs=NULL;
    *p_all_antiregs=NULL;
  }
}

/********************************************************************************/
/* Code for counting of complexes                                               */
/********************************************************************************/

/*************************************************************************
get_counting_regions_for_waypoint:
   Finds the regions and antiregions for a waypoint.  Note that this is only
   used for "macromol" counting -- presently, this is limited to counting the
   number of times a particular configuration of subunits occurs within a
   particular type of complex.

   In:  struct subvolume *my_sv - subvolume containing waypoint
        struct waypoint *wp - the waypoint
        struct region_list **p_all_regs - pointer to receive list of regions
        struct region_list **p_all_antiregs - pointer to receive list of antiregions
        struct pointer_hash *region_hash - hash whose keys are regions which
                                           are valid for counting
   Out: None
*************************************************************************/
static int get_counting_regions_for_waypoint(struct subvolume *my_sv,
                                             struct waypoint *wp,
                                             struct region_list **p_all_regs,
                                             struct region_list **p_all_antiregs,
                                             struct pointer_hash *region_hash)
{
  /* Copy all the potentially relevant regions from the nearest waypoint */
  struct region_list *rl;
  for ( rl=wp->regions ; rl!=NULL ; rl=rl->next)
  {
    if (pointer_hash_lookup(region_hash, rl->reg, rl->reg->hashval) == NULL)
      continue;

    struct region_list *nrl = (struct region_list*) CHECKED_MEM_GET(my_sv->local_storage->regl,
                                                                    "region list entry");
    nrl->reg=rl->reg;
    nrl->next=*p_all_regs;
    *p_all_regs=nrl;
  }

  /* And all the antiregions (regions crossed from inside to outside only) */
  for ( rl=wp->antiregions ; rl!=NULL ; rl=rl->next)
  {
    if (pointer_hash_lookup(region_hash, rl->reg, rl->reg->hashval) == NULL)
      continue;

    struct region_list *nrl = (struct region_list*) CHECKED_MEM_GET(my_sv->local_storage->regl,
                                                                    "region list entry");
    nrl->reg=rl->reg;
    nrl->next=*p_all_antiregs;
    *p_all_antiregs=nrl;    
  }

  return 0;
}

/*************************************************************************
get_counting_regions_for_point:
   Finds the regions and antiregions for a point.  Note that this is only used
   for "macromol" counting -- presently, this is limited to counting the number
   of times a particular configuration of subunits occurs within a particular
   type of complex.

   In:  struct subvolume *my_sv - subvolume containing point
        struct waypoint *wp - waypoint for subvolume containing point
        struct vector3 *loc - point for which to find counting regions
        struct region_list **p_all_regs - pointer to receive list of regions
        struct region_list **p_all_antiregs - pointer to receive list of antiregions
        struct pointer_hash *region_hash - hash whose keys are regions which
                                           are valid for counting
   Out: None
*************************************************************************/
static int get_counting_regions_for_point(struct subvolume *my_sv,
                                          struct waypoint *wp,
                                          struct vector3 *loc,
                                          struct region_list **p_all_regs,
                                          struct region_list **p_all_antiregs,
                                          struct pointer_hash *region_hash)
{
  struct region_list *all_regs=NULL, *all_antiregs=NULL;
  struct vector3 here;
  here.x = wp->loc.x;
  here.y = wp->loc.y;
  here.z = wp->loc.z;

  *p_all_regs=NULL;
  *p_all_antiregs=NULL;

  /* Get regions for waypoint */
  get_counting_regions_for_waypoint(my_sv, wp, &all_regs, &all_antiregs, region_hash);

  /* Raytrace across any walls from waypoint to us and add to region lists */
  struct vector3 delta;
  struct subvolume *sv;
  for ( sv = my_sv ; sv != NULL ; sv = next_subvol(&here,&delta,sv) )
  {
    delta.x = loc->x - here.x;
    delta.y = loc->y - here.y;
    delta.z = loc->z - here.z;

    /* When do we hit a subvolume boundary? */
    double t_sv_hit = collide_sv_time(&here,&delta,sv);
    if (t_sv_hit > 1.0) t_sv_hit = 1.0;

    /* Check for collision with each wall */
    struct wall_list *wl;
    for (wl = sv->wall_head ; wl != NULL ; wl = wl->next)
    {
      struct vector3 hit;
      double t_hit;

      /* Skip walls which do not participate in counting */
      if (! (wl->this_wall->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)))
        continue;

      /* Check for collision with wall, skip wall if we didn't hit it. */
      int j = collide_wall(&here,&delta,wl->this_wall,&t_hit,&hit,0);
      if (j == COLLIDE_MISS)
        continue;
      world->ray_polygon_colls++;

      /* Skip this collision if it's on the far side of the waypoint */
      if (t_hit > t_sv_hit)
        continue;

      /* This test might be superfluous, but I don't want to remove it until
       * we're sure.  Rex thinks it might be there to avoid problems with
       * roundoff error when the point is very near a wall, which is very
       * plausible.
       */
      if ((hit.x-loc->x)*delta.x + (hit.y-loc->y)*delta.y + (hit.z-loc->z)*delta.z >= 0)
        continue;

      /* Scan over all counting regions for the wall... */
      struct region_list *rl;
      for (rl=wl->this_wall->counting_regions ; rl!=NULL ; rl=rl->next)
      {
        /* Skip irrelevant regions */
        if (pointer_hash_lookup(region_hash, rl->reg, rl->reg->hashval) == NULL)
          continue;

        /* Add this region to the appropriate list */
        struct region_list *nrl = (struct region_list*) CHECKED_MEM_GET(my_sv->local_storage->regl,
                                                                        "region list entry");
        nrl->reg = rl->reg;
        if (j==COLLIDE_FRONT)
        {
          nrl->next=all_regs;
          all_regs=nrl;
        }
        else if (j==COLLIDE_BACK)
        {
          nrl->next=all_antiregs;
          all_antiregs=nrl;
        }
      }
    }
  }

  /* Clean up region lists */
  if (all_regs!=NULL && all_antiregs!=NULL)
    clean_region_lists(my_sv, &all_antiregs, &all_antiregs);

  *p_all_regs=all_regs;
  *p_all_antiregs=all_antiregs;
  return 0;
}

/*************************************************************************
scan_complex_update_table:
   Scans over the update table for a given counter in a given complex.  The
   update table we are scanning over is built for a particular subunit.  That
   is, once we get here, we're dealing only with rules regarding the relatives
   of the subunit in question, and do not need to be concerned with the state
   of the subunit itself.

   In:  struct species **relatives - the states of all related subunits (NULL
                                     if a subunit is empty)
        int num_relatives - the number of related subunits in the 'relatives'
                                     array
        struct complex_counter *counter - the counter containing the update
                                     table
        int rules_start - the start index within the table for our rules
        int rules_end - the end index within the table for our rules
        int n - the number to add to the count for any matching rules.
                                     generally, this is -1 for the "before"
                                     state, and +1 for the "after" state
   Out: None
*************************************************************************/
static void scan_complex_update_table(struct species **relatives,
                                      short *orients,
                                      int num_relatives,
                                      struct complex_counter *counter,
                                      int rules_start,
                                      int rules_end,
                                      int n)
{
  struct species **nptr = counter->neighbors + rules_start * num_relatives;
  int *iptr = counter->invert + rules_start * num_relatives;
  signed char *optr = NULL;

  if (orients)
    optr = counter->orientations + rules_start * num_relatives;

  /* For each rule, check whether we match, and update the count if we do */
  for (int rule_index = rules_start;
       rule_index < rules_end;
       ++ rule_index, nptr += num_relatives, iptr += num_relatives, optr += num_relatives)
  {

    /* For each clause of the rule, check whether we match */
    int neighbor_index;
    for (neighbor_index = 0; neighbor_index < num_relatives; ++ neighbor_index)
    {
      /* NULL means no clause in this rule mentions the neighbor, so the
       * species matches automatically.  Note that this clause also is used to
       * match the orientation of the reference subunit.
       */
      if (nptr[neighbor_index] == NULL)
      {
        /* XXX: invert is ignored if nptr is null.  This is ok right now
         * because no syntax can specify a null molecule with orientation, so
         * this is only useful for the subunit of reference.
         */

        /* If orients is NULL, this is a volume molecule.  Continue. */
        if (orients == NULL)
          continue;

        /* If optr[neighbor_index] is 0, we don't care about orientation.  Continue. */
        else if (optr[neighbor_index] == 0)
          continue;

        /* If optr[neighbor_index] matches orientation, continue. */
        else if (optr[neighbor_index] == orients[neighbor_index])
          continue;

        /* Orientation mismatch.  Next rule. */
        else
          break;
      }

      /* If a relative is NULL, that means there is no subunit there yet, so we
       * never match, unless nptr was NULL
       */
      if (relatives[neighbor_index] == NULL)
        break;

      /* Otherwise, check if we match this clause of the rule */
      if (iptr[neighbor_index])
      {
        if (nptr[neighbor_index] != relatives[neighbor_index])
          continue;

        if (orients == NULL  ||  optr[neighbor_index] * orients[neighbor_index] >= 0)
          break;
      }
      else
      {
        if (nptr[neighbor_index] != relatives[neighbor_index])
          break;

        if (orients != NULL  &&  optr[neighbor_index] * orients[neighbor_index] < 0)
          break;
      }
    }

    /* If we got through all of the clauses, then we matched.  Update the
     * count.
     */
    if (neighbor_index == num_relatives)
      counter->counts[rule_index] += n;
  }
}

/*************************************************************************
count_complex_for_single_region:
   Updates the counts of subunit states for a macromolecular complex within a
   particular region.  This is called from count_complex, which iterates over
   all regions which are doing counting on the complex.  This function, in
   turn, builds up a "before" and "after" state of all relevant related
   subunits, and then subtracts the counts from the old state and adds the
   counts from the new state.

   In:  struct complex_counter *c - the counter for the region in question
        struct complex_species *spec - the species of the complex we're counting
        short this_orient - the orientation of the complex as a whole
        struct species **before - states of all subunits before the update
        short *orient_before - orientations of all subunits before the update
        struct species **after - states of all subunits after the update
        short *orient_after - orientations of all subunits after the update
        int *update_subunit - an array of flags indicating whether each subunit
                              might have been affected by the update, in such a
                              way as to change the counts.  Presently, this
                              means the updated subunit, as well as any
                              subunits which have relationships leading to the
                              updated subunit.
        int amount - 1 for regions, -1 for antiregions
   Out: None
*************************************************************************/
static void count_complex_for_single_region(struct complex_counter *c,
                                            struct complex_species *spec,
                                            int this_orient,
                                            struct species **before,
                                            short *orient_before,
                                            struct species **after,
                                            short *orient_after,
                                            int *update_subunit,
                                            int amount)
{
  /* Now, add all relevant subunit updates to the counters */
  struct species *relatives_before[ spec->num_relations + 1], *relatives_after[ spec->num_relations + 1];
  short relatives_orient_before[ spec->num_relations + 1], relatives_orient_after[ spec->num_relations + 1];

  for (int subunit_index = 0; subunit_index < spec->num_subunits; ++ subunit_index)
  {
    /* Skip subunits that will not have changed */
    if (! update_subunit[subunit_index])
      continue;
    if (before[subunit_index] == NULL  &&  after[subunit_index] == NULL)
      continue;

    /* Set up tables of relations */
    int relation_idx;
    int offset = 0;
    if (orient_before)
    {
      relatives_before[0] = NULL;
      relatives_after[0]  = NULL;
      relatives_orient_before[0] = orient_before[subunit_index];
      relatives_orient_after[0] = orient_after[subunit_index];
      ++offset;
    }
    for (relation_idx = 0; relation_idx < spec->num_relations; ++ relation_idx)
    {
      int target_index = spec->relations[ relation_idx ].target[ subunit_index ];
      relatives_before[relation_idx + offset] = before[target_index];
      relatives_after[relation_idx + offset] = after[target_index];
      if (orient_before) relatives_orient_before[relation_idx + offset] = orient_before[target_index];
      if (orient_after)  relatives_orient_after[relation_idx + offset] = orient_after[target_index];
    }

    for (struct complex_counter *cCur = c; cCur != NULL; cCur = cCur->next)
    {
      if (this_orient != 0  &&  c->this_orient != 0  &&  c->this_orient != this_orient)
        continue;

      /* Remove "before" counts */
      if (before[subunit_index] != NULL)
      {
        int *before_indices = (int *) pointer_hash_lookup(&c->subunit_to_rules_range, before[subunit_index], before[subunit_index]->hashval);
        if (before_indices != NULL)
          scan_complex_update_table(relatives_before,
                                    orient_before ? relatives_orient_before : NULL,
                                    spec->num_relations + offset,
                                    c,
                                    before_indices[0],
                                    before_indices[1],
                                    -amount);
      }

      /* Add "after" counts */
      if (after[subunit_index] != NULL)
      {
        int *after_indices = (int *) pointer_hash_lookup(&c->subunit_to_rules_range, after[subunit_index], after[subunit_index]->hashval);
        if (after_indices != NULL)
          scan_complex_update_table(relatives_after,
                                    orient_before ? relatives_orient_after : NULL,
                                    spec->num_relations + offset,
                                    c,
                                    after_indices[0],
                                    after_indices[1],
                                    amount);
      }
    }
  }
}

/*************************************************************************
count_complex_new_for_single_region:
   Adds counts for a newly created macromolecule for counters in a
   single region.

   In:  struct complex_counter *c - the counter for the region in question
        struct complex_species *spec - the species of the complex we're counting
        short this_orient - the orientation of the complex as a whole
        struct species **specs - initial states of all subunits
        short *orients - initial orientations of all subunits
        int amount - 1 for regions, -1 for antiregions
   Out: None
*************************************************************************/
static void count_complex_new_for_single_region(struct complex_counter *c,
                                                struct complex_species *spec,
                                                int this_orient,
                                                struct species **specs,
                                                short *orients,
                                                int amount)
{
  /* Now, add all relevant subunit updates to the counters */
  struct species *relatives[ spec->num_relations + 1];
  short relatives_orient[ spec->num_relations + 1];

  for (int subunit_index = 0; subunit_index < spec->num_subunits; ++ subunit_index)
  {
    /* Really, this shouldn't happen, but for consistency with the
     * other counting rules, we'll do this. */
    if (specs[subunit_index] == NULL)
      continue;

    /* Set up tables of relations */
    int offset = 0;
    if (orients)
    {
      relatives[0] = NULL;
      relatives_orient[0] = orients[subunit_index];
      ++offset;
    }

    for (int relation_idx = 0; relation_idx < spec->num_relations; ++ relation_idx)
    {
      int target_index = spec->relations[ relation_idx ].target[ subunit_index ];
      relatives[relation_idx + offset] = specs[target_index];
      if (orients) relatives_orient[relation_idx + offset] = orients[target_index];
    }

    for (struct complex_counter *cCur = c; cCur != NULL; cCur = cCur->next)
    {
      int *indices;
      if (this_orient != 0  &&  c->this_orient != 0  &&  c->this_orient != this_orient)
        continue;

      /* Add counts */
      indices = (int *) pointer_hash_lookup(&c->subunit_to_rules_range, specs[subunit_index], specs[subunit_index]->hashval);
      if (indices != NULL)
        scan_complex_update_table(relatives,  relatives_orient, spec->num_relations + offset, c, indices[0],  indices[1],   amount);
    }
  }
}

/*************************************************************************
count_complex:
   Updates the counts of subunit states for a macromolecular complex.  Whenever
   one of the subunits of a complex changes state, this should be called.  Note
   that this should be called after the state of the complex has been updated.
   Currently, this includes only two events:
      - subunit is created (when a complex is created by a release event)
      - subunit changes to a different type of subunit (when a reaction occurs)

   In:  struct volume_molecule *cmplex - the molecule representing the complex
        struct volume_molecule *replaced_subunit - the subunit which has been
                        replaced (or NULL if there was no subunit before)
        int replaced_subunit_idx - the index within the subunits of the
                        updated subunit
   Out: 0 on success, 1 on failure
*************************************************************************/
int count_complex(struct volume_molecule *cmplex,
                  struct volume_molecule *replaced_subunit,
                  int replaced_subunit_idx)
{
  struct complex_species *spec = (struct complex_species *) cmplex->properties;
  if (spec->counters == NULL)
    return 0;

  const int px = bisect(world->x_partitions,world->nx_parts,cmplex->pos.x);
  const int py = bisect(world->y_partitions,world->ny_parts,cmplex->pos.y);
  const int pz = bisect(world->z_partitions,world->nz_parts,cmplex->pos.z);

  /* Find the waypoint and subvolume containing our point */
  const int this_sv = pz + (world->nz_parts-1)*( py + (world->ny_parts-1)*px );
  struct waypoint *wp = &(world->waypoints[this_sv]);
  struct subvolume *my_sv = &(world->subvol[this_sv]);

  /* Find out which regions contain this complex */
  struct region_list *all_regs;
  struct region_list *all_antiregs;
  get_counting_regions_for_point(my_sv,
                                 wp,
                                 &cmplex->pos,
                                 &all_regs,
                                 &all_antiregs,
                                 &spec->counters->region_to_counter);

  /* Figure out which subunits of this complex will need to be recounted */
  /* XXX: Restrict this to only relationships for which counting is done? */
  int update_subunit[ spec->num_subunits ];
  macro_count_inverse_related_subunits(spec, update_subunit, replaced_subunit_idx);
  update_subunit[replaced_subunit_idx] = 1;

  /* Build up array of before+after subunits */
  struct species *before[ spec->num_subunits ];
  struct species *after[ spec->num_subunits ];
  for (int subunit_index = 0; subunit_index < spec->num_subunits; ++ subunit_index)
  {
    struct volume_molecule *mol = cmplex->cmplx[subunit_index + 1];
    before[subunit_index] = after[subunit_index] = mol ? mol->properties : NULL;
  }
  before[replaced_subunit_idx] = replaced_subunit ? replaced_subunit->properties : NULL;

  /* Do any relevant counting for WORLD */
  count_complex_for_single_region(&spec->counters->in_world, spec, 0, before, NULL, after, NULL, update_subunit, 1);

  /* Now, for each region, do all relevant counting */
  for (struct region_list *rl = all_regs; rl != NULL; rl = rl->next)
  {
    /* Get the counters for the complex within this region */
    struct complex_counter *c = (struct complex_counter *) pointer_hash_lookup(&spec->counters->region_to_counter, rl->reg, rl->reg->hashval);
    if (c == NULL)
      continue;

    count_complex_for_single_region(c, spec, 0, before, NULL, after, NULL, update_subunit, 1);
  }
  for (struct region_list *rl = all_antiregs; rl != NULL; rl = rl->next)
  {
    /* Get the counters for the complex within this region */
    struct complex_counter *c = (struct complex_counter *) pointer_hash_lookup(&spec->counters->region_to_counter, rl->reg, rl->reg->hashval);
    if (c == NULL)
      continue;

    count_complex_for_single_region(c, spec, 0, before, NULL, after, NULL, update_subunit, -1);
  }

  /* Free region memory */ 
  if (all_regs!=NULL) mem_put_list(my_sv->local_storage->regl,all_regs);
  if (all_antiregs!=NULL) mem_put_list(my_sv->local_storage->regl,all_antiregs);
  return 0;
}

/*************************************************************************
count_complex_surface:
   Updates the counts of subunit states for a macromolecular complex.  Whenever
   one of the subunits of a complex changes state, this should be called.  Note
   that this should be called after the state of the complex has been updated.
   Currently, this includes only two events:
      - subunit is created (when a complex is created by a release event)
      - subunit changes to a different type of subunit (when a reaction occurs)

   In:  struct volume_molecule *cmplex - the molecule representing the complex
        struct volume_molecule *replaced_subunit - the subunit which has been
                        replaced (or NULL if there was no subunit before)
        int replaced_subunit_idx - the index within the subunits of the
                        updated subunit
   Out: 0 on success, 1 on failure
*************************************************************************/
int count_complex_surface(struct grid_molecule *cmplex,
                          struct grid_molecule *replaced_subunit,
                          int replaced_subunit_idx)
{
  struct complex_species *spec = (struct complex_species *) cmplex->properties;
  if (spec->counters == NULL)
    return 0;

  /* Figure out which subunits of this complex will need to be recounted */
  /* XXX: Restrict this to only relationships for which counting is done? */
  int update_subunit[ spec->num_subunits ];
  macro_count_inverse_related_subunits(spec, update_subunit, replaced_subunit_idx);
  update_subunit[replaced_subunit_idx] = 1;

  /* Build up array of before+after subunits */
  struct species *before[ spec->num_subunits ];
  struct species *after[ spec->num_subunits ];
  short orient_before[ spec->num_subunits ];
  short orient_after[ spec->num_subunits ];
  for (int subunit_index = 0; subunit_index < spec->num_subunits; ++ subunit_index)
  {
    struct grid_molecule *mol = cmplex->cmplx[subunit_index + 1];
    before[subunit_index] = after[subunit_index] = mol ? mol->properties : NULL;
    orient_before[subunit_index] = orient_after[subunit_index] = mol ? mol->orient : 0;
  }
  before[replaced_subunit_idx] = replaced_subunit ? replaced_subunit->properties : NULL;
  orient_before[replaced_subunit_idx] = replaced_subunit ? replaced_subunit->orient : 0;

  /* Do any relevant counting for WORLD */
  count_complex_for_single_region(&spec->counters->in_world, spec, cmplex->orient, before, orient_before, after, orient_after, update_subunit, 1);

  struct wall *my_wall = cmplex->grid->surface;
  if (my_wall!=NULL && (my_wall->flags&COUNT_CONTENTS)!=0)
  {
    for (struct region_list *rl=my_wall->counting_regions ; rl!=NULL ; rl=rl->next)
    {
      /* Get the counters for the complex within this region */
      struct complex_counter *c = (struct complex_counter *) pointer_hash_lookup(&spec->counters->region_to_counter, rl->reg, rl->reg->hashval);
      if (c == NULL)
        continue;

      count_complex_for_single_region(c, spec, cmplex->orient, before, orient_before, after, orient_after, update_subunit, 1);
    }
  }
  return 0;
}

/*************************************************************************
count_complex_surface_new:
   Adds a new macromolecular surface complex to our count.

   In:  struct volume_molecule *cmplex - the molecule representing the complex
   Out: 0 on success, 1 on failure
*************************************************************************/
int count_complex_surface_new(struct grid_molecule *cmplex)
{
  struct complex_species *spec = (struct complex_species *) cmplex->properties;
  if (spec->counters == NULL)
    return 0;

  /* Build up array of before+after subunits */
  struct species *specs[ spec->num_subunits ];
  short orients[ spec->num_subunits ];
  for (int subunit_index = 0; subunit_index < spec->num_subunits; ++ subunit_index)
  {
    struct grid_molecule *mol = cmplex->cmplx[subunit_index + 1];
    specs[subunit_index] = mol ? mol->properties : NULL;
    orients[subunit_index] = mol ? mol->orient : 0;
  }

  /* Do any relevant counting for WORLD */
  count_complex_new_for_single_region(&spec->counters->in_world, spec, cmplex->orient, specs, orients, 1);

  struct wall *my_wall = cmplex->grid->surface;
  if (my_wall!=NULL && (my_wall->flags&COUNT_CONTENTS)!=0)
  {
    for (struct region_list *rl = my_wall->counting_regions; rl!=NULL; rl=rl->next)
    {
      /* Get the counters for the complex within this region */
      struct complex_counter *c = (struct complex_counter *) pointer_hash_lookup(&spec->counters->region_to_counter, rl->reg, rl->reg->hashval);
      if (c == NULL)
        continue;

      count_complex_new_for_single_region(c, spec, cmplex->orient, specs, orients, 1);
    }
  }
  return 0;
}

/*************************************************************************
macro_collect_count_requests_by_subunit:
   Sorts the count requests out by subunit species type, storing them as lists
   in the pointer hash provided.

   In:  struct pointer_hash *h - the table to hold the sorted requests
        struct macro_count_request *requests - the requests to sort
   Out: 0 on success, 1 on failure
*************************************************************************/
static int macro_collect_count_requests_by_subunit(struct pointer_hash *h,
                                                   struct macro_count_request *requests)
{
  struct macro_count_request *mcr, *mcrnext;
  int total_entries = 0;
  for (mcr = requests; mcr != NULL; mcr = mcrnext)
  {
    mcrnext = mcr->next;
    mcr->next = (struct macro_count_request *) pointer_hash_lookup(h, mcr->subunit_state, mcr->subunit_state->hashval);
    if (pointer_hash_add(h, mcr->subunit_state, mcr->subunit_state->hashval, mcr))
      mcell_allocfailed("Failed to add complex counters to hash.");
    ++ total_entries;
  }

  return total_entries;
}

/*************************************************************************
macro_copy_count_requests_to_tables:
   Copies the rules from the count requests into tables.

   In:  struct pointer_hash *requests_by_subunit - the requests, sorted by subunit
        struct pointer_hash *subunit_to_rules_range - index giving range of rules in table by subunit
        struct species const **nptr - the neighbor species table
        int *iptr - the "invert" table
        int num_relations - number of relations (i.e. width of nptr/iptr table)
        int *su_rules_indices - array of indices within nptr and iptr of rules by subunit
        int *counts - array of ints which hold the actual counts during the simulation
   Out: 0 on success, 1 on failure
*************************************************************************/
static int macro_copy_count_requests_to_tables(struct pointer_hash *requests_by_subunit,
                                               struct pointer_hash *subunit_to_rules_range,
                                               struct species **nptr,
                                               int *iptr,
                                               signed char *optr,
                                               int num_relations,
                                               int *su_rules_indices,
                                               int *counts)
{
  int table_position = 0, start_pos;
  int su_index = 0;
  for (int bin_index = 0; bin_index < requests_by_subunit->table_size; ++ bin_index)
  {
    /* Skip empty bins */
    if (requests_by_subunit->keys[bin_index] == NULL  ||
        requests_by_subunit->values[bin_index] == NULL)
      continue;

    /* Mark the start of this entry */
    start_pos = table_position;

    /* Now, process all items for this subunit */
    struct macro_count_request *head = (struct macro_count_request *) requests_by_subunit->values[bin_index];
    struct macro_count_request *mcr, *mcrnext;
    int offset = optr ? 1 : 0;
    for (mcr = head; mcr != NULL; mcr = mcrnext)
    {
      mcrnext = mcr->next;

      /* Point the expression at the correct counter */
      mcr->paired_expression->left = &counts[table_position];
      mcr->paired_expression->expr_flags &= ~OEXPR_LEFT_MACROREQUEST;
      mcr->paired_expression->expr_flags |= OEXPR_LEFT_INT;

      /* Surface counters have an extra column for the reference subunit */
      if (optr)
      {
        nptr[0] = NULL;
        iptr[0] = 0;
        if (mcr->subunit_orientation > 0) optr[0] = 1;
        else if (mcr->subunit_orientation < 0) optr[0] = -1;
        else optr[0] = 0;
      }

      /* Copy into the table */
      struct macro_relation_state *msr, *msrnext;
      for (msr = mcr->relation_states; msr != NULL; msr = msrnext)
      {
        msrnext = msr->next;
        nptr[offset + msr->relation] = msr->mol;
        iptr[offset + msr->relation] = msr->invert;
        if (optr)
        {
          if (msr->orient > 0) optr[offset + msr->relation] = 1;
          else if (msr->orient < 0) optr[offset + msr->relation] = -1;
          else optr[offset + msr->relation] = 0;
        }
        free(msr);
      }

      /* Move to the next table slot */
      ++ table_position;
      nptr += num_relations;
      iptr += num_relations;
      if (optr)
        optr += num_relations;

      free(mcr);
    }

    /* Add to the index */
    struct species *key = (struct species *) requests_by_subunit->keys[bin_index];
    if (pointer_hash_add(subunit_to_rules_range, key, key->hashval, su_rules_indices + su_index))
      mcell_allocfailed("Failed to initialize macromolecule state count rules table.");
    su_rules_indices[su_index++] = start_pos;
    su_rules_indices[su_index++] = table_position;
  }

  return 0;
}

/*************************************************************************
macro_sort_output_requests_by_orientation:
    Sort out count requests by desired complex orientation.

   In:  struct macro_count_request *requests - the counts we've requested
        struct macro_count_request *(*by_orientation)[3] - the three lists to
                     contain the 'unoriented' requests, the 'positive' oriented
                     requests, and the 'negative' oriented requests.
   Out: None.  'by_orientation' array has been updated with requests separated
        by orientation.
*************************************************************************/
static void macro_sort_output_requests_by_orientation(struct macro_count_request *requests,
                                                      struct macro_count_request *(*by_orientation)[3])
{
  while (requests != NULL)
  {
    struct macro_count_request *next = requests->next;
    if (requests->master_orientation == 0)
    {
      requests->next = (*by_orientation)[0];
      (*by_orientation)[0] = requests;
    }
    else if (requests->master_orientation < 0)
    {
      requests->next = (*by_orientation)[2];
      (*by_orientation)[2] = requests;
    }
    else /* if (requests->master_orientation > 0) */
    {
      requests->next = (*by_orientation)[1];
      (*by_orientation)[1] = requests;
    }

    requests = next;
  }
}

/*************************************************************************
macro_initialize_counters_for_complex:
   Initializes the counters for a complex, copying the rules from the count
   requests into tables.

   In:  struct complex_species *spec - the species for which to initialize
        struct complex_counter *c - the counter to initialize
        struct macro_count_request *requests - the requests with which to initialize
   Out: 0 on success, 1 on failure
*************************************************************************/
static int macro_initialize_counters_for_complex(struct complex_species *spec,
                                                 struct complex_counter *c,
                                                 struct macro_count_request *requests)
{
  struct macro_count_request *by_orientation[3] = {NULL, NULL, NULL};    /* 0, 1, -1 */
  macro_sort_output_requests_by_orientation(requests, &by_orientation);

  /* Prepare a hash to sort our requests by subunit */
  struct pointer_hash requests_by_subunit;
  if (pointer_hash_init(&requests_by_subunit, 16))
    mcell_allocfailed("Failed to initialize subunit->request hash for complex counters.");

  struct complex_counter **cur = &c;

  /* Loop over orientation types (0, 1, -1). */
  for (int n_orients=0; n_orients<3; ++n_orients)
  {
    if (by_orientation[n_orients] == NULL)
      continue;

    if (*cur == NULL)
    {
      *cur = CHECKED_MALLOC_STRUCT(struct complex_counter,
                                   "macromolecular complex counter");
      memset(*cur, 0, sizeof(struct complex_counter));
      if (pointer_hash_init(&(*cur)->subunit_to_rules_range, 16))
        mcell_allocfailed("Failed to initialize subunit->rules hash for complex counters.");
    }
    c = *cur;
    c->this_orient = by_orientation[n_orients]->master_orientation;

    /* Sort our requests by subunit */
    int total_entries = macro_collect_count_requests_by_subunit(&requests_by_subunit, by_orientation[n_orients]);
    if (total_entries < 0)
    {
      pointer_hash_destroy(&requests_by_subunit);
      return 1;
    }

    /* Now, allocate space for tables */
    int is_surface = (spec->base.flags & ON_GRID) ? 1 : 0;
    int num_relations = is_surface ? (spec->num_relations + 1) : spec->num_relations;
    c->neighbors = CHECKED_MALLOC_ARRAY(struct species *,
                                        num_relations * total_entries,
                                        "macromolecular complex count table");
    c->invert = CHECKED_MALLOC_ARRAY(int,
                                     num_relations * total_entries,
                                     "macromolecular complex count table");
    c->counts = CHECKED_MALLOC_ARRAY(int,
                                     total_entries,
                                     "macromolecular complex count table");
    c->su_rules_indices = CHECKED_MALLOC_ARRAY(int,
                                               total_entries * 2,
                                               "macromolecular complex count table");
    memset(c->neighbors, 0, num_relations * total_entries * sizeof(struct species *));
    memset(c->invert, 0, num_relations * total_entries * sizeof(int));
    memset(c->counts, 0, total_entries * sizeof(int));
    memset(c->su_rules_indices, 0, total_entries * 2 * sizeof(int));

    /* If this is not a volume molecule, allocate space for orientations */
    if (is_surface)
    {
      c->orientations = CHECKED_MALLOC_ARRAY(signed char,
                                             num_relations * total_entries,
                                             "macromolecular complex count table");
      memset(c->orientations, 0, num_relations * total_entries * sizeof(signed char));
    }
    else
      c->orientations = NULL;

    /* Now, fill in the tables */
    if (macro_copy_count_requests_to_tables(&requests_by_subunit,
                                            &c->subunit_to_rules_range,
                                            c->neighbors,
                                            c->invert,
                                            c->orientations,
                                            num_relations,
                                            c->su_rules_indices,
                                            c->counts))
      goto failure;

    /* Advance to next orientation */
    cur = &(c->next);
  }

  pointer_hash_destroy(&requests_by_subunit);
  return 0;

failure:
  pointer_hash_destroy(&requests_by_subunit);
  return 1;
}

/*************************************************************************
macro_destroy_counters:
   Destroys the counters structure, freeing all memory.

   In:  struct complex_species *spec - the species to receive a counters struct
   Out: None
*************************************************************************/
static void macro_destroy_counters(struct complex_species *spec)
{
  pointer_hash_destroy(&spec->counters->region_to_counter);
  pointer_hash_destroy(&spec->counters->in_world.subunit_to_rules_range);
  if (spec->counters->in_world.su_rules_indices)
    free(spec->counters->in_world.su_rules_indices);
  if (spec->counters->in_world.neighbors)
    free(spec->counters->in_world.neighbors);
  if (spec->counters->in_world.invert)
    free((void *) spec->counters->in_world.invert);
  if (spec->counters->in_world.counts)
    free(spec->counters->in_world.counts);
  if (spec->counters->in_regions)
  {
    for (int region_idx = 0; region_idx < spec->counters->num_region_counters; ++ region_idx)
    {
      if (spec->counters->in_regions[region_idx].su_rules_indices)
        free(spec->counters->in_regions[region_idx].su_rules_indices);
      if (spec->counters->in_regions[region_idx].neighbors)
        free(spec->counters->in_regions[region_idx].neighbors);
      if (spec->counters->in_regions[region_idx].invert)
        free((void *) spec->counters->in_regions[region_idx].invert);
      if (spec->counters->in_regions[region_idx].counts)
        free(spec->counters->in_regions[region_idx].counts);
      pointer_hash_destroy(&spec->counters->in_regions[region_idx].subunit_to_rules_range);
    }
    free(spec->counters->in_regions);
  }
  free(spec->counters);
  spec->counters = NULL;
}

/*************************************************************************
macro_create_counters:
   Creates and initializes the counters structure for a species.  The structure
   is completely empty after this function finishes.  (i.e. no counters are
   defined).

   In:  struct complex_species *spec - the species to receive a counters struct
   Out: 0 on success, 1 on failure
*************************************************************************/
static int macro_create_counters(struct complex_species *spec, struct complex_counters **dest)
{
  /* Allocate the counters */
  *dest = CHECKED_MALLOC_STRUCT(struct complex_counters,
                                "macromolecular complex counters");
  memset(*dest, 0, sizeof(struct complex_counters));

  /* Initialize the hash tables in the counter */
  if (pointer_hash_init(&(*dest)->in_world.subunit_to_rules_range, 16)  ||
      pointer_hash_init(&(*dest)->region_to_counter, 16))
  {
    mcell_allocfailed("Failed to initialize top-level hashes for complex counters.");
    macro_destroy_counters(spec);
    return 1;
  }

  return 0;
}

/*************************************************************************
macro_collect_count_requests_by_location:
   Separate the count requests out based on their location.  Count requests for
   the entire world go into the "in_world" list, and other requests go into
   lists in the in_region hash table, keyed by region.

   In:  struct macro_count_request *requests - the requests to sort
        struct macro_count_request **in_world - list for WORLD requests
        struct pointer_hash *in_region - hash to contain per-region lists
   Out: 0 on success, 1 on failure
*************************************************************************/
static int macro_collect_count_requests_by_location(struct macro_count_request *requests,
                                                    struct macro_count_request **in_world,
                                                    struct pointer_hash *in_region)
{
  /* Scan over all requests */
  struct macro_count_request *mcr, *mcrnext;
  for (mcr = requests; mcr != NULL; mcr = mcrnext)
  {
    mcrnext = mcr->next;

    /* If location is NULL, request is for entire world */
    if (mcr->location == NULL)
    {
      mcr->next = *in_world;
      *in_world = mcr;
    }

    /* Request is for a specific region */
    else
    {
      struct region *r = (struct region *) mcr->location->value;
      mcr->next = (struct macro_count_request *) pointer_hash_lookup(in_region, r, r->hashval);
      if (pointer_hash_add(in_region, r, r->hashval, mcr))
      {
        mcell_allocfailed("Failed to add complex count request to region->request hash.");
        return 1;
      }
    }
  }

  return 0;
}

/*************************************************************************
macro_create_region_counters:
   For the specified set of complex counters, create structures to hold the
   counters for the requested number of regions.

   In:  struct complex_counters *c - the counters into which to add region counters
        int num_region_counters - the number of regions for which we're counting
   Out: 0 on success, 1 on failure
        All region-specific counters for the complex species are allocated.
*************************************************************************/
static int macro_create_region_counters(struct complex_counters *c,
                                        int num_region_counters)
{
  /* Allocate space for the counters */
  c->in_regions = CHECKED_MALLOC_ARRAY(struct complex_counter,
                                       num_region_counters,
                                       "macromolecular complex per-region counters");
  memset(c->in_regions, 0, num_region_counters * sizeof(struct complex_counter));
  c->num_region_counters = num_region_counters;

  /* Now initialize each counter */
  for (int counter_index = 0; counter_index < c->num_region_counters; ++ counter_index)
  {
    if (pointer_hash_init(&c->in_regions[counter_index].subunit_to_rules_range, 16))
      mcell_allocfailed("Failed to initialize subunit->rules counter for region.");
  }

  return 0;
}

/*************************************************************************
macro_convert_output_requests_for_complex:
   For the specified complex, convert all counter requests into actual
   counters attached to the complex species.

   In:  struct complex_species *spec - the species for which to fill in counters
        struct macro_count_request *requests - the counts we've requested
   Out: 0 on success, 1 on failure (memory allocation only?).
        All counters for the complex species are instantiated and populated.
*************************************************************************/
static int macro_convert_output_requests_for_complex(struct complex_species *spec,
                                                     struct macro_count_request *requests)
{
  /* Make sure we've got a counter for this guy */
  if (spec->counters == NULL  &&  macro_create_counters(spec, &spec->counters))
    return 1;

  struct macro_count_request *in_world = NULL;
  struct pointer_hash requests_by_region;
  if (pointer_hash_init(&requests_by_region, 16))
  {
    mcell_allocfailed("Failed to initialize region->requests hash.");
    goto failure;
  }

  /* Iterate over the requests, sorting them out by location */
  if (macro_collect_count_requests_by_location(requests, &in_world, &requests_by_region))
    goto failure;

  /* Fill in world counter, if appropriate */
  if (in_world  &&  macro_initialize_counters_for_complex(spec, &spec->counters->in_world, in_world))
    goto failure;

  /* Create counters for all regions */
  if (requests_by_region.num_items != 0)
  {
    if (macro_create_region_counters(spec->counters, requests_by_region.num_items))
      goto failure;

    /* Now, fill in each set of counters */
    int counter_index = 0;
    for (int bin_index = 0; bin_index < requests_by_region.table_size; ++ bin_index)
    {
      /* Skip empty bins */
      if (requests_by_region.keys[bin_index] == NULL  ||  requests_by_region.values[bin_index] == NULL)
        continue;

      struct complex_counter *my_counter = spec->counters->in_regions + counter_index++;
      struct macro_count_request *my_requests = (struct macro_count_request *) requests_by_region.values[bin_index];
      if (pointer_hash_add(&spec->counters->region_to_counter,
                           requests_by_region.keys[bin_index],
                           requests_by_region.hashes[bin_index],
                           my_counter))
        mcell_allocfailed("Failed to initialize macromolecule state count regions table.");
      if (macro_initialize_counters_for_complex(spec, my_counter, my_requests))
        goto failure;
    }
  }

  pointer_hash_destroy(&requests_by_region);
  return 0;

failure:
  pointer_hash_destroy(&requests_by_region);
  if (spec->counters)
    macro_destroy_counters(spec);
  return 1;
}

/*************************************************************************
macro_expand_object_output:
   Normalize the location of the given output request, converting all requests
   for counts on objects into requests for counts in the "ALL" region on the
   object.

   In:  request for a complex count
        object upon which the request is made.
   Out: 0 on success, 1 on failure (memory allocation only?).

        Request is split into a separate request for each BOX and POLY
        object's ALL region that is a child of this object.  The result
        is then added up here.
*************************************************************************/
static int macro_expand_object_output(struct macro_count_request *request,struct object *obj)
{
  switch (obj->object_type)
  {
    case META_OBJ:
#ifdef ALLOW_COUNTS_ON_METAOBJECT
#error "Support for counting macromolecules in/on a metaobject is unimplemented."
#endif
      mcell_error("COUNT and TRIGGER statements on metaobject '%s' are not allowed.", obj->sym->name);
      return 1;

    case REL_SITE_OBJ:
      mcell_error("COUNT and TRIGGER statements on release object '%s' are not allowed.", obj->sym->name);
      return 1;

    case BOX_OBJ:
    case POLY_OBJ:
      request->location = NULL;
      for (struct region_list *rl = obj->regions; rl != NULL; rl = rl->next)
      {
        if (is_reverse_abbrev(",ALL",rl->reg->sym->name))
        {
          rl->reg->flags|=COUNT_CONTENTS;
          request->location = rl->reg->sym;
          break;
        }
      }
      if (request->location == NULL)
        mcell_internal_error("ALL region missing on object %s", obj->sym->name);
      break;

    default:
      mcell_internal_error("Invalid object type (%d) in count on object expansion", obj->object_type);
      return 1;
  }

  return 0;
}

/*************************************************************************
macro_normalize_output_request_locations:
   Prepare all macromolecule output requests for the simulation.

   In:  None
   Out: 0 on success, 1 on failure
        Locations for macromolecule count requests are checked for validity,
        and normalized so that all locations are regions.
*************************************************************************/
static int macro_normalize_output_request_locations(void)
{
  /* Scan all requests, fixing up request locations */
  for (struct macro_count_request *mcr = world->macro_count_request_head; mcr != NULL; mcr = mcr->next)
  {
    /* If the location is "WORLD", we're done */
    if (mcr->location == NULL )
      continue;

    /* Now, make sure the object referenced is actually instantiated in the
     * world
     */
    if (! is_object_instantiated(mcr->location))
    {
      mcell_error("The object/region name '%s' in the COUNT/TRIGGER statement is not fully referenced.\n"
                  "  This occurs when a count is requested on an object which has not been referenced\n"
                  "  (directly or indirectly) from an INSTANTIATE block in the MDL file.",
                  mcr->location->name);
      return 1;
    }

    /* If the location is an object, convert it to a region */
    if (mcr->location->sym_type == OBJ)
    {
      if (macro_expand_object_output(mcr, (struct object*) mcr->location->value))
        mcell_error("Failed to expand request to count on object.");
    }
    else
    {
      struct region *reg = (struct region *) mcr->location->value;
      reg->flags |= COUNT_CONTENTS;
    }
  }
  return 0;
}

/*************************************************************************
macro_collect_count_requests_by_complex:
   Prepare all macromolecule output requests for the simulation.

   In:  struct pointer_hash *h - the hash table into which to collect requests
        struct macro_count_request *head - the requests to sort

   Out: 0 on success, 1 on failure
        The requests are sorted into the hash table keyed by the species of the
        complex associated with the request.  The value of each element in the
        hash table will be a macro_count_request list.
*************************************************************************/
static int macro_collect_count_requests_by_complex(struct pointer_hash *h,
                                                   struct macro_count_request *head)
{
  struct macro_count_request *mcr, *mcrnext;
  for (mcr = head; mcr != NULL; mcr = mcrnext)
  {
    mcrnext = mcr->next;
    struct complex_species *c = mcr->the_complex;
    mcr->next = (struct macro_count_request *) pointer_hash_lookup(h, c, c->base.hashval);
    if (pointer_hash_add(h, c, c->base.hashval, mcr))
      mcell_allocfailed("Failed to add count request to complex->counter hash.");
  }

  return 0;
}

/*************************************************************************
macro_convert_output_requests:
   Prepare all macromolecule output requests for the simulation.

   In: None
   Out: 0 on success, 1 on failure
        All macro_count_request objects are converted into counter structures
        attached to the associated complex, and the references in the output
        expressions are fixed to point to the appropriate counters.  Locations
        are normalized to refer to regions.
*************************************************************************/
static int macro_convert_output_requests(void)
{
  /* If we have no requests to process, skip all this */
  if (world->macro_count_request_head == NULL)
    return 0;

  /* Check that all locations are valid count locations */
  if (macro_normalize_output_request_locations())
    return 1;

  /* Scan over the requests, sorting them out by complex */
  struct pointer_hash complex_to_requests;
  if (pointer_hash_init(&complex_to_requests, 16))
    mcell_allocfailed("Failed to initialize complex->requests hash.");
  if (macro_collect_count_requests_by_complex(&complex_to_requests, world->macro_count_request_head))
    goto failure;
  world->macro_count_request_head = NULL;

  /* Now, handle the requests complex-by-complex */
  for (int n_bin = 0; n_bin < complex_to_requests.table_size; ++n_bin)
  {
    /* Skip empty bins */
    if (complex_to_requests.keys[n_bin] == NULL  ||  complex_to_requests.values[n_bin] == NULL)
      continue;

    if (macro_convert_output_requests_for_complex((struct complex_species *) complex_to_requests.keys[n_bin],
                                                  (struct macro_count_request *) complex_to_requests.values[n_bin]))
      goto failure;
  }

  pointer_hash_destroy(&complex_to_requests);
  return 0;

failure:
  pointer_hash_destroy(&complex_to_requests);
  return 1;
}
