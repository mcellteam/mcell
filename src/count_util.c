/**************************************************************************\
** File: count_util.c                                                     **
**                                                                        **
** Purpose: Handles counting of interesting events                        **
**                                                                        **
** Testing status: untested.                                              **
\**************************************************************************/


#include <math.h>
#include <stdio.h>

#include "mcell_structs.h"
#include "wall_util.h"
#include "count_util.h"


extern struct volume *world;


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
      j = (rl->reg->hashval ^ sp->hashval)&world->collide_hashmask;
      if (j==0) j = sp->hashval & world->collide_hashmask;
      
      for (hit_count=world->collide_hash[j] ; hit_count!=NULL ; hit_count=hit_count->next)
      {
        if (hit_count->reg_type == rl->reg && hit_count->mol_type == sp)
        {
          if (crossed) hit_count->n_inside += direction;
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
  int i;
  struct region_list *rl;
  struct species *sp = me->properties;
  struct counter *c;
  
  if ((sp->flags & ON_GRID) != 0)
  {
    struct grid_molecule *g = (struct grid_molecule*)me;
    struct wall *w = g->grid->surface;
    if (w->flags & COUNT_CONTENTS)
    {
      for (rl=w->regions ; rl!=NULL ; rl=rl->next)
      {
        i = (rl->reg->hashval ^ sp->hashval) & world->collide_hashmask;
        if (i==0) i = sp->hashval & world->collide_hashmask;
        
        for ( c = world->collide_hash[i] ; c != NULL ; c = c->next )
        {
          if (c->reg_type == rl->reg && c->mol_type == sp) c->n_inside += n;
        }
      }
    } 
  }
}
