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
count_hit:
   In: molecule
       wall it hit
       direction of impact relative to normal
   Out: No return value.  Appropriate counters are updated.
*************************************************************************/

void count_hit(struct molecule *m,struct wall *w,int direction)
{
  int hashval;
  struct region_list *rl;
  struct counter *cp;
  
  world->ray_polygon_colls++;
  
  for (rl = w->regions ; rl != NULL ; rl = rl->next)
  {
    if (m->properties->hashval == rl->reg->hashval) hashval = m->properties->hashval;
    else hashval = m->properties->hashval ^ rl->reg->hashval;
    
    for (cp = world->collide_hash[hashval] ; cp != NULL ; cp = cp->next)
    {
      if (cp->mol_id == m->properties->hashval &&
          cp->wall_id == rl->reg->hashval)
      {
        if (direction==0) cp->impacts++;
        else cp->crossings += direction;
        break;
      }
    }
  }
}


/*************************************************************************
count_react:
   In: reaction
       path that reaction took
       the current timestep
   Out: No return value.  Appropriate counters are updated.
*************************************************************************/

void count_react(struct rxn *rx,int path,double timestep)
{
  rx->counter[path]++;
  
#if 0
  rx->rxn_count_cum[path]++;
  if ( floor(timestep) > rx->last_update ) rx->rxn_count_dt[path] = 1;
  else rx->rxn_count_dt[path]++;
#endif
}


/*************************************************************************
count_crossings:
   In: molecule
       subvolume to test
       displacement vector for the molecule's motion
   Out: No return value.  Appropriate counters are updated.
*************************************************************************/

void count_crossings(struct molecule *m,struct subvolume *sv,
                     struct vector3 *move)
{
  int i;
  struct wall_list *wl;
  double t;
  struct vector3 hitpt;
  
  for (wl = sv->wall_head ; wl != NULL ; wl = wl->next)
  {
    i = collide_wall(&(m->pos),move,wl->this_wall,&t,&hitpt);
    if (i==COLLIDE_REDO)
    {
      /* Ignore this for now.  Should fix. */
    }
    else if (i == COLLIDE_MISS) continue;
    else if (i == COLLIDE_FRONT) count_hit(m,wl->this_wall,1);
    else count_hit(m,wl->this_wall,-1);
  }
}
