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
#include "mcell_structs.h"
#include "wall_util.h"
#include "vol_util.h"
#include "count_util.h"


extern struct volume *world;

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
   Out: No return value.  Appropriate counters are updated.
*************************************************************************/

void update_collision_count(struct species *sp,struct region_list *rl,int direction,int crossed)
{
  int j;
  struct counter *hit_count;
  
  hit_count = NULL;  
  for ( ; rl != NULL ; rl = rl->next)
  {
    if (rl->reg->flags & COUNT_SOME)
    {
      j = (rl->reg->hashval ^ sp->hashval)&world->count_hashmask;
      if (j==0) j = sp->hashval & world->count_hashmask;
      
      for (hit_count=world->count_hash[j] ; hit_count!=NULL ; hit_count=hit_count->next)
      {
        if (hit_count->reg_type == rl->reg && hit_count->mol_type == sp)
        {
          if (crossed)
          {
            hit_count->n_inside += direction;

/*            printf("Counted %s (%x) on %s (%x); %x has n_inside = %.1f (up by %d).\n",
                   sp->sym->name,sp->hashval,rl->reg->sym->name,rl->reg->hashval,
                   (int)hit_count,hit_count->n_inside,direction); */
          }
          if (rl->reg->flags & sp->flags & COUNT_HITS)
          {
            if (crossed)
            {
              if (direction==1)
              {
                hit_count->front_hits++;
                hit_count->front_to_back++;
              }
              else
              {
                hit_count->back_hits++;
                hit_count->back_to_front++;
              }
            }
            else
            {
              if (direction==1) hit_count->front_hits++;
              else hit_count->back_hits++;
            }
          }
        }
      }
    }
  }
}


void find_enclosing_regions(struct vector3 *loc,struct vector3 *start,
                            struct region_list** rlp,struct region_list** arlp,
                            struct mem_helper *rmem)
{
  struct vector3 outside,delta,hit;
  struct subvolume *sv,*svt;
  struct wall_list *wl;
  struct region_list *rl,*arl;
  struct region_list *trl,*tarl,*xrl,*yrl,*nrl;
  double t;
  int traveling;
  int i;
  
  rl = *rlp;
  arl = *arlp;
  
  if (start==NULL || loc->x!=start->x || loc->y!=start->y || loc->z < start->z)
  {
    outside.x = loc->x;
    outside.y = loc->y;
    if (world->bb_min.z < 0) outside.z = world->bb_min.z * (1 + EPS_C);
    else outside.z = world->bb_min.z * (1 - EPS_C);
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
          trl = xrl;
        }
      }
      else if (i==COLLIDE_MISS || t > 1.0 || (wl->this_wall->flags & COUNT_CONTENTS) == 0) continue;
      else
      {
        for (xrl=wl->this_wall->regions ; xrl != NULL ; xrl = xrl->next)
        {
          if ((xrl->reg->flags & COUNT_CONTENTS) != 0)
          {
            nrl = (struct region_list*) mem_get(rmem);
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
    }
  }
  
  *rlp = rl;
  *arlp = arl;
}



struct region_list* dup_region_list(struct region_list *r,struct mem_helper *mh)
{
  struct region_list *nr,*rp,*r0;
  
  if (r==NULL) return NULL;
  
  r0 = rp = NULL;
  while (r!=NULL)
  {
    nr = (struct region_list*) mem_get(mh);
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
  int h,i,j,k;
  int i_will_use_waypoints = 0;
  int waypoint_in_wall = 0;
  struct waypoint *wp;
  struct wall_list *wl;
  struct subvolume *sv;
  double d;
  
  for (i=0;i<world->n_species;i++)
  {
    if ((world->species_list[i]->flags & (NOT_FREE | COUNT_CONTENTS)) == COUNT_CONTENTS)
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
  
  if (!world->waypoints) return 1;  /* Malloc failure */
  
  for (i=0;i<world->nx_parts-1;i++)
  {
    for (j=0;j<world->ny_parts-1;j++)
    {
      for (k=0;k<world->nz_parts-1;k++)
      {
        h = k + (world->nz_parts-1)*(j + (world->ny_parts-1)*i);
        wp = &(world->waypoints[h]);
        
        sv = &(world->subvol[h]);
        wp->loc.x = 0.5*( world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ] );
        wp->loc.y = 0.5*( world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ] );
        wp->loc.z = 0.5*( world->z_fineparts[ sv->llf.z ] + world->z_fineparts[ sv->urb.z ] );
        
        do
        {
          waypoint_in_wall = 0;
          for (wl = sv->wall_head ; wl != NULL ; wl = wl->next)
          {
            d = dot_prod( &(wp->loc) , &(wl->this_wall->normal) ); 
            if ( eps_equals( d , wl->this_wall->d ) )
            { 
              waypoint_in_wall++;
              d = EPS_C * (double)((rng_uint(world->seed++)&0xF) - 8);
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
          wp->regions = dup_region_list(world->waypoints[h-1].regions,sv->mem->regl);
          wp->antiregions = dup_region_list(world->waypoints[h-1].antiregions,sv->mem->regl);
          
          find_enclosing_regions(&(wp->loc),&(world->waypoints[h-1].loc),
                                 &(wp->regions),&(wp->antiregions),sv->mem->regl);
        }
        else
        {
          wp->regions = NULL;
          wp->antiregions = NULL;
          find_enclosing_regions(&(wp->loc),NULL,&(wp->regions),
                                 &(wp->antiregions),sv->mem->regl);
        }
        
/*        if (wp->regions != NULL) printf("We have a region on waypoint %d\n",h);
        if (wp->antiregions != NULL) printf("We have an antiregion on waypoint %d\n",h); */
        
      }
    }
  }
  
  return 0;
}



/*************************************************************************
count_me_by_region:
   In: abstract molecule we are supposed to count (or a representative one)
       number by which to update the counter (usually +1 or -1)
   Out: No return value.  Appropriate counters are updated.
   Note: This handles all types of molecules, from grid to free.  Right
         now, only grid is implemented.
*************************************************************************/

void count_me_by_region(struct abstract_molecule *me,int n)
{
  int i,j,k,h;
  struct region_list *rl;
  struct species *sp = me->properties;
  struct counter *c;
  
  //printf("Counting %x by region (up by %d)!\n",(int)me,n);
  
  if ((sp->flags & ON_GRID) != 0)
  {
    struct grid_molecule *g = (struct grid_molecule*)me;
    struct wall *w = g->grid->surface;

    if (w->flags & COUNT_CONTENTS)
    {
      for (rl=w->regions ; rl!=NULL ; rl=rl->next)
      {
        i = (rl->reg->hashval ^ sp->hashval) & world->count_hashmask;
        if (i==0) i = sp->hashval & world->count_hashmask;
        
        for ( c = world->count_hash[i] ; c != NULL ; c = c->next )
        {
          if (c->reg_type == rl->reg && c->mol_type == sp) c->n_inside += n;
        }
      }
    } 
  }
  else if ((sp->flags & ON_SURFACE) != 0)
  {
    /* Not implemented */
  }
  else /* Free molecule */
  {
    struct molecule *m = (struct molecule*)me;
    struct subvolume *sv;
    struct vector3 here,delta,hit;
    struct waypoint *wp;
    struct wall_list *wl;
    struct region_list *rl;
    double t;
    
    i = bisect(world->x_partitions,world->nx_parts,m->pos.x);
    j = bisect(world->y_partitions,world->ny_parts,m->pos.y);
    k = bisect(world->z_partitions,world->nz_parts,m->pos.z);
    h = k + (world->nz_parts-1)*( j + (world->ny_parts-1)*i );
    
    wp = &(world->waypoints[h]);
    for (rl=wp->regions ; rl!=NULL ; rl=rl->next)
    {
      if ( (rl->reg->flags & COUNT_CONTENTS) != 0 )
      {
        i = (rl->reg->hashval ^ sp->hashval) & world->count_hashmask;
        if (i==0) i = sp->hashval & world->count_hashmask;
/*        printf("Trying to count %s (%x) on %s (%x) with hashval %x\n",
               sp->sym->name,sp->hashval,rl->reg->sym->name,rl->reg->hashval,i); */
        
        for ( c = world->count_hash[i] ; c != NULL ; c = c->next )
        {
          if (c->reg_type==rl->reg && c->mol_type==sp)
          {
            c->n_inside += n;
/*            printf("Actually counted; %x has n_inside = %.1f (up by %d).\n",(int)c,c->n_inside,n); */
          }
        }
      }
    }
    for (rl=wp->antiregions ; rl!=NULL ; rl=rl->next)
    {
      if ( (rl->reg->flags & COUNT_CONTENTS) != 0 )
      {
        i = (rl->reg->hashval ^ sp->hashval) & world->count_hashmask;
        if (i==0) i = sp->hashval & world->count_hashmask;
        
        for ( c = world->count_hash[i] ; c != NULL ; c = c->next )
        {
          if (c->reg_type==rl->reg && c->mol_type==sp) c->n_inside -= n;
        }
      }
    }
    
    here.x = wp->loc.x;
    here.y = wp->loc.y;
    here.z = wp->loc.z;
    
    for ( sv = &(world->subvol[h]) ; sv != NULL ; sv = next_subvol(&here,&delta,sv) )
    {
      delta.x = m->pos.x - here.x;
      delta.y = m->pos.y - here.y;
      delta.z = m->pos.z - here.z;

      for (wl = sv->wall_head ; wl != NULL ; wl = wl->next)
      {
        if (wl->this_wall->flags & COUNT_CONTENTS)
        {
          j = collide_wall(&here,&delta,wl->this_wall,&t,&hit);
          
          if (j!=COLLIDE_MISS)
          {
            for (rl=wl->this_wall->regions ; rl!=NULL ; rl=rl->next)
            {
              if ( (rl->reg->flags & m->flags & COUNT_CONTENTS) != 0 )
              {
                i = (rl->reg->hashval ^ sp->hashval) & world->count_hashmask;
                if (i==0) i = sp->hashval & world->count_hashmask;
                
                for ( c = world->count_hash[i] ; c != NULL ; c = c->next )
                {
                  if (c->reg_type==rl->reg && c->mol_type==sp)
                  {
                    if (j==COLLIDE_FRONT) c->n_inside += n;
                    else if (j==COLLIDE_BACK) c->n_inside -= n;
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


int check_region_counters()
{
  FILE *log_file;
  struct counter *cp;
  struct species *sp;
  struct region *rp;
  u_int i;


  log_file=world->log_file;  
  
  for (i=0;i<world->count_hashmask+1;i++) {
    for (cp=world->count_hash[i];cp!=NULL;cp=cp->next) {
      sp=cp->mol_type;
      rp=cp->reg_type;
      /* if species is freely diffusing
         make sure region is a closed manifold */
      if ((sp->flags & NOT_FREE)==0) {
        if (rp->manifold_flag==MANIFOLD_UNCHECKED) {
          if (is_manifold(rp)) {
            rp->manifold_flag=IS_MANIFOLD;
          }
          else {
            rp->manifold_flag=NOT_MANIFOLD;
          }
        }
        else {
          if (rp->manifold_flag==NOT_MANIFOLD) {
            fprintf(log_file,"MCell: error, cannot count diffusing molecules inside non-closed object region: %s\n",rp->sym->name); 
          }
        }
      }
    }
  }

  return(0);
}
