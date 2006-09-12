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
       direction of impact relative to surface normal
       whether we crossed or not
       scaling factor for reaction probabilities (for estimating ccn)
       location of the hit (for triggers)
       time of the hit (for triggers)
   Out: No return value.  Appropriate counters are updated, that is,
        hit counters are updated according to which side was hit,
	and crossings counters and counts within enclosed regions are
	updated if the surface was crossed.
*************************************************************************/

int count_region_update(struct species *sp,struct region_list *rl,int direction,int crossed, double factor,struct vector3 *loc,double t)
{
  int i,j;
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
          if (rl->reg->flags & sp->flags & (COUNT_HITS|COUNT_CONTENTS|COUNT_ENCLOSED))
          {
            if (crossed)
            {
              if (direction==1)
              {
                if (hit_count->counter_type&TRIG_COUNTER)
                {
                  hit_count->data.trig.t_event=t;
		  if (rl->reg->flags&sp->flags&COUNT_HITS)
		  {
                    i=fire_count_event(hit_count,1,loc,REPORT_FRONT_HITS|REPORT_TRIGGER);
		    if (i) return 1;
                    i=fire_count_event(hit_count,1,loc,REPORT_FRONT_CROSSINGS|REPORT_TRIGGER);
		    if (i) return 1;
		  }
		  if (rl->reg->flags&sp->flags&COUNT_CONTENTS)
		  {
		    i=fire_count_event(hit_count,1,loc,REPORT_ENCLOSED|REPORT_CONTENTS|REPORT_TRIGGER);
		    if (i) return 1;
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
                  hit_count->data.trig.t_event=t;
		  if (rl->reg->flags&sp->flags&COUNT_HITS)
		  {
                    i=fire_count_event(hit_count,1,loc,REPORT_BACK_HITS|REPORT_TRIGGER);
		    if (i) return 1;
                    i=fire_count_event(hit_count,1,loc,REPORT_BACK_CROSSINGS|REPORT_TRIGGER);
		    if (i) return 1;
		  }
		  if (rl->reg->flags&sp->flags&COUNT_CONTENTS)
		  {
		    i=fire_count_event(hit_count,-1,loc,REPORT_ENCLOSED|REPORT_CONTENTS|REPORT_TRIGGER);
		    if (i) return 1;
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
                  hit_count->data.trig.t_event=t;
                  i=fire_count_event(hit_count,1,loc,REPORT_FRONT_HITS|REPORT_TRIGGER);
		  if (i) return 1;
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
                  i=fire_count_event(hit_count,1,loc,REPORT_BACK_HITS|REPORT_TRIGGER);
		  if (i) return 1;
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
  
  return 0;
}



/*************************************************************************
count_region_from_scratch:
   In: molecule to count, or NULL
       reaction pathname to count, or NULL
       number of these to count
       location at which to count them (may be NULL)
       wall at which this happened (may be NULL)
       time of the hit (for triggers)
   Out: No return value.  Appropriate counters are updated and triggers
        are fired.
   Note: At least one of molecule or rxn pathname must be non-NULL; if
        other inputs are NULL, sensible values will be guessed (which
        may themselves be NULL).  This routine is not super-fast for
        volume counts (enclosed counts) since it has to dynamically create
        and test lists of enclosing regions.
*************************************************************************/

int count_region_from_scratch(struct abstract_molecule *am,struct rxn_pathname *rxpn,int n,struct vector3 *loc,struct wall *my_wall,double t)
{  
  int i,j,k,h;
  struct region_list *rl,*arl,*prl,*parl,*nrl,*narl; /*a=anti p=previous n=new*/
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
  }
  
  /* Count grid molecules and reactions on surfaces--easy */
  if (my_wall!=NULL && (my_wall->flags&COUNT_CONTENTS)!=0)
  {
    for (rl=my_wall->counting_regions ; rl!=NULL ; rl=rl->next)
    {
      i=(hashval^rl->reg->hashval)&world->count_hashmask;
      if (i==0) i=hashval&world->count_hashmask;
      for (c=world->count_hash[i] ; c!=NULL ; c=c->next)
      {
	if (c->target==target && c->reg_type==rl->reg && (c->counter_type&ENCLOSING_COUNTER)==0)
	{
	  if (c->counter_type&TRIG_COUNTER)
	  {
	    c->data.trig.t_event=t;
	    k=fire_count_event(c,n,loc,count_flags|REPORT_TRIGGER);
	    if (k) return 1;
	  }
	  else if (rxpn==NULL) c->data.move.n_at+=n;
	  else c->data.rx.n_rxn_at+=n;
	}
      }
    }
  }
  
  /* Count volume molecules, vol reactions, and surface stuff that is enclosed--hard!!*/
  if (am==NULL || (am->properties->flags&COUNT_ENCLOSED)!=0 || (am->properties->flags&NOT_FREE)==0)
  {
    i = bisect(world->x_partitions,world->nx_parts,loc->x);
    j = bisect(world->y_partitions,world->ny_parts,loc->y);
    k = bisect(world->z_partitions,world->nz_parts,loc->z);
    
    h = k + (world->nz_parts-1)*( j + (world->ny_parts-1)*i );
    wp = &(world->waypoints[h]);
    my_sv = &(world->subvol[h]);
    
    here.x = wp->loc.x;
    here.y = wp->loc.y;
    here.z = wp->loc.z;
    
    all_regs=NULL;
    all_antiregs=NULL;
    
    /* Copy all the potentially relevant regions from the nearest waypoint */
    for ( rl=wp->regions ; rl!=NULL ; rl=rl->next)
    {
      i=(hashval^rl->reg->hashval)&world->count_hashmask;
      if (i==0) i=hashval&world->count_hashmask;
      if (world->count_hash[i]==NULL) continue; /* Won't count on this region so ignore it */
      
      nrl=(struct region_list*)mem_get(my_sv->local_storage->regl);
      if (nrl==NULL)
      {
	fprintf(world->err_file,"Error at file %s line %d\n  Out of memory making list of enclosing regions for count\n",__FILE__,__LINE__);
	return 1;
      }
      nrl->reg=rl->reg;
      nrl->next=all_regs;
      all_regs=nrl;
    }
    
    /* And all the antiregions (regions crossed from inside to outside only) */
    for ( arl=wp->antiregions ; arl!=NULL ; arl=arl->next)
    {
      i=(hashval^arl->reg->hashval)&world->count_hashmask;
      if (i==0) i=hashval&world->count_hashmask;
      if (world->count_hash[i]==NULL) continue; /* Won't count on this region so ignore it */

      narl=(struct region_list*)mem_get(my_sv->local_storage->regl);
      if (narl==NULL)
      {
	fprintf(world->err_file,"Error at file %s line %d\n  Out of memory making list of enclosing regions for count\n",__FILE__,__LINE__);
	return 1;
      }
      narl->reg=arl->reg;
      narl->next=all_antiregs;
      all_antiregs=narl;    
    }
    
    /* Raytrace across any walls from waypoint to us and add to region lists */
    for ( sv = &(world->subvol[h]) ; sv != NULL ; sv = next_subvol(&here,&delta,sv) )
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
	  j = collide_wall(&here,&delta,wl->this_wall,&t_hit,&hit);
          if((j != COLLIDE_MISS) && (world->notify->final_summary == NOTIFY_FULL)){
              world->ray_polygon_colls++;
          }	 

 
	  if (j!=COLLIDE_MISS && t_hit <= t_sv_hit &&
	    (hit.x-loc->x)*delta.x + (hit.y-loc->y)*delta.y + (hit.z-loc->z)*delta.z < 0)
	  {
	    for (rl=wl->this_wall->counting_regions ; rl!=NULL ; rl=rl->next)
	    {
	      if ( (rl->reg->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0 )
	      {
		i=(hashval^rl->reg->hashval)&world->count_hashmask;
		if (i==0) i=hashval&world->count_hashmask;
		if (world->count_hash[i]==NULL) continue; /* Won't count on this region so ignore it */
		
		nrl = (struct region_list*)mem_get(my_sv->local_storage->regl);
		if (nrl==NULL)
		{
		  fprintf(world->err_file,"Error at file %s line %d\n  Out of memory making list of enclosing regions for count\n",__FILE__,__LINE__);
		  return 1;
		}
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
	}
      }
    }
    
    /* Clean up region lists */
    if (all_regs!=NULL && all_antiregs!=NULL)
    {
      if (all_regs->next!=NULL || all_antiregs->next!=NULL)
      {
	struct region_list pre_sentry,pre_antisentry;
	
	/* Sort by memory address to make mutual annihilation faster */
	if (all_regs->next!=NULL) all_regs=(struct region_list*)void_list_sort((struct void_list*)all_regs);
	if (all_antiregs->next!=NULL) all_antiregs=(struct region_list*)void_list_sort((struct void_list*)all_antiregs);
	
	/* Need previous entry to fix up list, so we'll make an imaginary one for 1st list element */
	pre_sentry.next=all_regs;
	pre_antisentry.next=all_antiregs;
	prl=&pre_sentry;
	parl=&pre_antisentry;
	
	/* If we cross a region both ways, throw both out (once) */
	for (rl=all_regs,arl=all_antiregs ; rl!=NULL && arl!=NULL ; prl=rl,rl=rl->next,parl=arl,arl=arl->next)
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
	all_regs=pre_sentry.next;
	all_antiregs=pre_antisentry.next;
      }
      else if (all_regs->reg==all_antiregs->reg)
      {
	/* Crossed one region both ways, toss them */
	mem_put(my_sv->local_storage->regl,all_regs);
	mem_put(my_sv->local_storage->regl,all_antiregs);
	all_regs=NULL;
	all_antiregs=NULL;
      }
    }
    
    /* Actually check the regions here */
    count_flags|=REPORT_ENCLOSED;
    
    for (nrl=all_regs ; nrl!=NULL ; nrl=(nrl==all_regs)?all_antiregs:NULL) /* Trick so we don't need to duplicate this code */
    {
      if (nrl==all_regs) pos_or_neg=1;
      else pos_or_neg=-1;
      for (rl=nrl ; rl!=NULL ; rl=rl->next)
      {
	i = (hashval^rl->reg->hashval)&world->count_hashmask;
	if (i==0) i=hashval&world->count_hashmask;
	
	for (c=world->count_hash[i] ; c!=NULL ; c=c->next)
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
	      k = fire_count_event(c,n*pos_or_neg,loc,count_flags|REPORT_TRIGGER);
	      if (k) return 1;
	    }
	    else if (rxpn==NULL) c->data.move.n_enclosed += n*pos_or_neg;
	    else c->data.rx.n_rxn_enclosed += n*pos_or_neg;
	  }
	}
      }
    }
    
    /* Free region memory */ 
    if (all_regs!=NULL) mem_put_list(my_sv->local_storage->regl,all_regs);
    if (all_antiregs!=NULL) mem_put_list(my_sv->local_storage->regl,all_antiregs);
  }
  
  return 0;
}



/*************************************************************************
fire_count_event:
   In: counter of thing that just happened (trigger of some sort)
       number of times that thing happened
       location where it happened
       what happened (Report Type Flags)   
   Out: 0 on success, 1 on error (memory allocation or file I/O).
*************************************************************************/

int fire_count_event(struct counter *event,int n,struct vector3 *where,byte what)
{
  struct trigger_request *tr;
  byte whatelse=what;
  int i;
  
  if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_HITS) whatelse = (what-REPORT_FRONT_HITS)|REPORT_ALL_HITS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_BACK_HITS) whatelse = (what-REPORT_BACK_HITS)|REPORT_ALL_HITS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_CROSSINGS) whatelse = (what-REPORT_FRONT_CROSSINGS)|REPORT_ALL_CROSSINGS;
  else if ((what&REPORT_TYPE_MASK)==REPORT_BACK_CROSSINGS) whatelse = (what-REPORT_BACK_CROSSINGS)|REPORT_ALL_CROSSINGS;
  
  for (tr=event->data.trig.listeners ; tr!=NULL ; tr=tr->next)
  {
    if (tr->ear->report_type==what)
    {
      memcpy(&(event->data.trig.loc),where,sizeof(struct vector3));
      i=add_trigger_output(event,tr->ear,n);
      if (i) return 1;
    }
    else if (tr->ear->report_type==whatelse)
    {
      memcpy(&(event->data.trig.loc),where,sizeof(struct vector3));
      if ((what&REPORT_TYPE_MASK)==REPORT_FRONT_HITS || (what&REPORT_TYPE_MASK)==REPORT_FRONT_CROSSINGS)
      {
        i=add_trigger_output(event,tr->ear,n);
      }
      else
      {
        i=add_trigger_output(event,tr->ear,-n);
      }
      if (i) return 1;
    }
  }
  return 0;
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
      
      if((i != COLLIDE_MISS) && (world->notify->final_summary == NOTIFY_FULL)){
          world->ray_polygon_colls++;
      }	 

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
      else if (i==COLLIDE_MISS || !(t >= 0 && t < 1.0) || t > t_hit_sv || (wl->this_wall->flags & (COUNT_CONTENTS|COUNT_RXNS|COUNT_ENCLOSED)) == 0 ||
	       (hit.x-outside.x)*delta.x + (hit.y-outside.y)*delta.y + (hit.z-outside.z)*delta.z < 0) continue;
      else
      {
        for (xrl=wl->this_wall->counting_regions ; xrl != NULL ; xrl = xrl->next)
        {
          if ((xrl->reg->flags & (COUNT_CONTENTS|COUNT_RXNS|COUNT_ENCLOSED)) != 0)
          {
            nrl = (struct region_list*) mem_get(rmem);
	    if (nrl==NULL)
	    {
	      fprintf(world->err_file, "File '%s', Line %ld:  Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
	      i = emergency_output();
	      fprintf(world->err_file, "Fatal error: out of memory while finding enclosing regions.\nAttempt to write intermediate results had %d errors\n", i);
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
place_waypoints:
   In: No arguments.
   Out: Returns 1 if malloc fails, 0 otherwise.
        Allocates waypoints to SSVs, if any are needed.
   Note: you must have initialized SSVs before calling this routine!
*************************************************************************/

int place_waypoints()
{
  int g,h,i,j,k;
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
              if(world->notify->final_summary == NOTIFY_FULL){
                  world->random_number_use++;
              }
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
  struct output_block *block;
  struct output_set *set;
  struct output_column *column;
  int i;
  
  /* First give everything a sensible name, if needed */
  for (block=world->output_block_head ; block!=NULL ; block=block->next)
  {
    for (set=block->data_set_head ; set!=NULL ; set=set->next)
    {
      if (set->header_comment==NULL) continue;
      for (column=set->column_head ; column!=NULL ; column=column->next)
      {
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

int check_counter_geometry()
{
  int i;
  struct counter *cp;
  struct region *rp;
  
  /* Check to make sure what we've created is geometrically sensible */
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
	  fprintf(world->err_file,"Cannot count molecules or events inside non-manifold object region: %s\n", rp->sym->name); 
	  return (1);
	}
	
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
*************************************************************************/

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

int instantiate_request(struct output_request *request)
{
  int request_hash;
  struct rxn_pathname *rxpn_to_count;
  struct rxn *rx_to_count = NULL;
  struct species *mol_to_count = NULL;
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
    if (request->report_type&REPORT_ENCLOSED) request->report_type -= REPORT_ENCLOSED;
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
        fprintf(world->err_file,"Internal error at file %s line %d\n  Invalid report type 0x%x in count request.\n",__FILE__,__LINE__,report_type_only);
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
      reg_of_count->flags|=COUNT_ENCLOSED;
      count_type|=ENCLOSING_COUNTER;
      if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_ENCLOSED;
    }
    if (request->report_type&REPORT_TRIGGER)
    {
      count_type|=TRIG_COUNTER;
      reg_of_count->flags|=COUNT_TRIGGER;
    }
    
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
    
    is_enclosed = ((request->report_type&REPORT_ENCLOSED)!=0);
    
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
      
      if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_TRIGGER;
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
        case REPORT_CONCENTRATION:
          if (mol_to_count!=NULL) mol_to_count->flags|=COUNT_HITS;
          reg_of_count->flags|=COUNT_HITS;
          break;
        default:
          fprintf(world->err_file,"Error at file %s line %d\n  Bad report type %d when creating counts\n",__FILE__,__LINE__,report_type_only);
          return 1;
          break;
      }
    }
    else /* Not trigger--set up for regular count */
    {
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



