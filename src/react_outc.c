/**************************************************************************\
** File: react_outc.c                                                     **
**                                                                        **
** Purpose: Implements specific reaction outcome pathways.                **
**                                                                        **
** Testing status: compiles.                                              **
**                                                                        **
** Note: Compartment entry/exit not implemented for flipping (will be)    **
** Note: Counting of reaction paths not implemented (simple, will be)     **
\**************************************************************************/

#include <string.h>

#include "rng.h"
#include "grid_util.h"
#include "mcell_structs.h"
#include "count_util.h"
#include "react.h"
#include "vol_util.h"

extern struct volume *world;


/*************************************************************************
outcome_products:
   In: relevant wall in the interaction
       first free molecule in the interaction, if any
       first surface molecule in the interaction, if any
       first grid molecule in the interaction, if any
       reaction that is occuring
       path that the reaction is taking
       local storage for creating new molecules
       orientations of the molecules in the reaction
       time that the reaction is occurring
   Out: Value depending on how orientations changed--
          -2 reaction blocked by full grid
          -1 moving molecule reflected
          0 everything went fine, nothing extra to do
          1 moving molecule needs to flip
        Products are created as necessary and scheduled.
*************************************************************************/

int outcome_products(struct wall *w,struct molecule *reac_m,
  struct surface_molecule *reac_s,struct grid_molecule *reac_g,
  struct rxn *rx,int path,struct storage *local,
  short orientA,short orientB,double t,struct vector3 *hitpt,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  struct abstract_molecule *moving)
{
  int blocked = 0;

  int i;
  struct molecule *m;
  struct surface_molecule *s;
  struct grid_molecule *g;
  struct species *p;
  struct surface_grid *sg;
  u_int bits;
  int j,k;
  double f;
  int i0 = rx->product_idx[path];
  int iN = rx->product_idx[path+1];
  
  struct abstract_molecule *plist[iN-i0];
  char ptype[iN-i0];
  short porient[iN-i0];
  
  if (reacA->properties == rx->players[1] && reacA->properties != rx->players[0])
  {
    plist[0] = reacA;
    reacA = reacB;
    reacB = plist[0];
    
    j = orientA;
    orientA = orientB;
    orientB = (short)j;
  }
  
  plist[0] = reacA;
  if ( (reacA->properties->flags&ON_GRID)!=0 ) ptype[0] = 'g';
  else if ( (reacA->properties->flags&ON_SURFACE)!=0 ) ptype[0] = 's';
  else if ( (reacA->properties->flags&IS_SURFACE)==0 ) ptype[0] = 'm';
  else ptype[0] = '!';
  if (rx->n_reactants > 1)
  {
    if (reacB == NULL)
    {
      ptype[1] = 'w';
      plist[1] = NULL;
    }
    else
    {
      plist[1] = reacB;
      if ( (reacB->properties->flags&ON_GRID)!=0 ) ptype[1] = 'g';
      else if ( (reacB->properties->flags&ON_SURFACE)!=0 ) ptype[1] = 's';
      else if ( (reacB->properties->flags&IS_SURFACE)==0 ) ptype[1] = 'm';
      else ptype[1] = '!';
    }
  }

  for (i=i0+rx->n_reactants;i<iN;i++)
  {
    p = rx->players[i];
  
    if ( (p->flags & ON_GRID) != 0 )
    {
      if (reac_g!=NULL || reac_s!=NULL || (reac_m!=NULL && w!=NULL))
      {
        j = -1;
        if (reac_g!=NULL) /* We will replace reac_g */
        {
          j = reac_g->grid_index;
          sg = reac_g->grid;
        }
        else if (reac_s!=NULL)
        {
          sg = reac_s->curr_wall->effectors;
          if (sg==NULL)
          {
            create_grid(reac_s->curr_wall,reac_s->subvol);
            sg = reac_s->curr_wall->effectors;
          }
          if (sg!=NULL)
          {
            j = uv2grid(&(reac_s->s_pos),sg);
            if (sg->mol[j]!=NULL)
            {
              j = -1; /* Slot full, no rxn */
              return -2;
            }
          }
        }
        else
        {
          sg = w->effectors;
          if (sg == NULL)
          {
            create_grid(w,reac_m->subvol);
            sg = w->effectors;
          }
          if (sg != NULL)
          {
            j = xyz2grid(hitpt,sg);
            if (sg->mol[j]!=NULL)
            {
              j = -1; /* Slot full, no rxn */
              return -2;
            }
          }
        }
        if (j>-1)
        {
          g = mem_get(local->gmol);
          g->birthplace = local->gmol;
          g->properties = p;
          p->population++;
          g->flags = TYPE_GRID + ACT_NEWBIE + IN_SCHEDULE;
          if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)g)!= NULL)
            g->flags += ACT_REACT;
          
          g->t = t;
          g->t2 = 0.0;
          g->grid_index = j;
          g->grid = sg;
          
          if (reac_g==NULL || sg->mol[j]!=reac_g) sg->n_occupied++;
          sg->mol[j] = g;
          
          plist[i-i0] = (struct abstract_molecule*)g;
          ptype[i-i0] = 'g';
          
          if (p->flags & COUNT_CONTENTS)
            count_me_by_region((struct abstract_molecule*)g,1);

          schedule_add(local->timer,g);
          
        }
        else
        {
          plist[i-i0] = NULL;
          ptype[i-i0] = 0;
        }
      }
      else
      {
        plist[i-i0] = NULL;
        ptype[i-i0] = 0;
        continue;
      }
    }
    else if ( (p->flags & ON_SURFACE) != 0 )
    {
      printf("Don't use me.\n");
      /* TODO: fix up so we use hitpt to set our location if avail */
      if ( reac_s != NULL || reac_g!=NULL || (reac_m!=NULL && w!=NULL))
      {
        s = mem_get(local->smol);
        s->birthplace = local->smol;
        s->properties = p;
        p->population++;
        s->flags = TYPE_SURF + ACT_NEWBIE + IN_SURFACE + IN_VOLUME + IN_SCHEDULE;
        if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)s) != NULL)
          s->flags += ACT_REACT;
        
        if (reac_s != NULL)
        {
          s->pos.x = reac_s->pos.x;
          s->pos.y = reac_s->pos.y;
          s->pos.z = reac_s->pos.z;
          s->s_pos.u = reac_s->s_pos.u;
          s->s_pos.v = reac_s->s_pos.v;
          s->curr_wall = reac_s->curr_wall;
          s->next_s = reac_s->next_s;
          reac_s->next_s = s;
          s->next_v = reac_s->next_v;
          reac_s->next_v = (struct molecule*)s;
        }
        else
        {
          if (reac_m != NULL)
          {
            s->pos.x = reac_m->pos.x;
            s->pos.y = reac_m->pos.y;
            s->pos.z = reac_m->pos.z;
            s->subvol = reac_m->subvol;
            s->next_v = s->subvol->mol_head;
            s->subvol->mol_head = (struct molecule*)s;
          }
          else
          {
            grid2xyz(reac_g->grid , reac_g->grid_index , &(s->pos));
            /* Could set uv too (faster)! */
            /* TODO: Look up subvolume. */
          }
          if (w != NULL) s->curr_wall = w;
          else s->curr_wall = reac_g->grid->surface;
          s->next_s = s->curr_wall->mol;
          s->curr_wall->mol = s;
          s->s_pos.u = s->pos.x * s->curr_wall->unit_u.x +
                       s->pos.y * s->curr_wall->unit_u.y +
                       s->pos.z * s->curr_wall->unit_u.z;
          s->s_pos.v = s->pos.x * s->curr_wall->unit_v.x +
                       s->pos.y * s->curr_wall->unit_v.y +
                       s->pos.z * s->curr_wall->unit_v.z;
        }
        s->curr_wall->mol_count++;
        plist[i-i0] = (struct abstract_molecule*)s;
        ptype[i-i0] = 's';
        s->t2 = 0.0;
        s->t = t;
        schedule_add( local->timer , s );
        
      }
      else
      {
        plist[i-i0] = NULL;
        ptype[i-i0] = 0;
        continue;
      }
    }
    else
    {
      m = mem_get(local->mol);
      m->birthplace = local->mol;
      m->properties = p;
      p->population++;
      m->collisions = 0;
      if (reac_g != NULL)
      {
        m->releaser = reac_g->grid;
        m->index = reac_g->grid_index;
      }
      else
      {
        m->releaser = 0;
        m->index = -1;
      }
      m->flags = TYPE_3D + ACT_NEWBIE + IN_VOLUME + IN_SCHEDULE;
      if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)m) != NULL)
        m->flags += ACT_REACT;
      if (p->space_step > 0.0) m->flags += ACT_DIFFUSE;
      
      if (hitpt != NULL)
      {
        m->pos.x = hitpt->x;
        m->pos.y = hitpt->y;
        m->pos.z = hitpt->z;
      }
      
      if (reac_m != NULL)
      {
        if (hitpt==NULL)
        {
          m->pos.x = reac_m->pos.x;
          m->pos.y = reac_m->pos.y;
          m->pos.z = reac_m->pos.z;
        }
        m->subvol = reac_m->subvol;
        
        if (w==NULL)
        {
          m->next_v = m->subvol->mol_head;
          m->subvol->mol_head = m;
          m->subvol->mol_count++;
        }
      }
      else if (reac_s != NULL)
      {
        if (hitpt==NULL)
        {
          m->pos.x = reac_s->pos.x;
          m->pos.y = reac_s->pos.y;
          m->pos.z = reac_s->pos.z;
        }
        m->subvol = reac_s->subvol;
        m->next_v = m->subvol->mol_head;
        m->subvol->mol_head = m;
      }
      else if (reac_g != NULL)
      {
        if (hitpt==NULL) grid2xyz(reac_g->grid , reac_g->grid_index , &(m->pos));
        m->subvol = find_subvolume(&(m->pos),reac_g->grid->subvol);
      }
      plist[i-i0] = (struct abstract_molecule*)m;
      ptype[i-i0] = 'm';
      m->t2 = 0.0;
      m->t = t;

      if (p->flags & COUNT_CONTENTS)
        count_me_by_region((struct abstract_molecule*)m,1);

      schedule_add( local->timer , m );
      
    }
  }
  
  bits = rng_uint( world->seed++ );
  for (i=i0;i<iN;i++,bits>>=1)
  {
    if (rx->players[i]==NULL) continue;
    
    if ( ptype[i-i0] != 0 && (ptype[i-i0]!='m' || w!=NULL) )
    {
      if (rx->geometries[i] == 0)
      {
        if ((bits&1)==0) porient[i-i0] = 1;
        else porient[i-i0] = -1;
      }
      else
      {
        if (rx->geometries[i] < 0)
        {
          j = -rx->geometries[i];
          k = -1;
        }
        else
        {
          j = rx->geometries[i];
          k = 1;
        }
        
        if (j > rx->n_reactants) porient[i-i0] = k*porient[j-(rx->n_reactants+1)];
        else if (j==1) porient[i-i0] = k*orientA;
        else if (j==2 && reacB!=NULL) porient[i-i0] = k*orientB;
        else porient[i-i0] = k;
        
      }
      
//      if (i-i0 < rx->n_reactants && porient[i-i0]==0) continue;

      if (ptype[i-i0]=='s')
      {
        ((struct surface_molecule*)plist[i-i0])->orient = porient[i-i0];
      }
      else if (ptype[i-i0]=='g')
      {
        ((struct grid_molecule*)plist[i-i0])->orient = porient[i-i0];
      }
      else if (moving == plist[i-i0])
      {
        if (moving==reacA)
        {
          if (orientA==porient[i-i0]) blocked = -1;
          else blocked = 1;
        }
        else
        {
          if (orientB==porient[i-i0]) blocked = -1;
          else blocked = 1;
        }
      }
      else if (ptype[i-i0]=='m')
      {
        m = (struct molecule*)plist[i-i0];
        if (porient[i-i0]>0) f = EPS_C;
        else f = -EPS_C;
/*        f *= m->properties->space_step*10.0/EPS_C; */
        	
        m->pos.x += f*w->normal.x;
        m->pos.y += f*w->normal.y;
        m->pos.z += f*w->normal.z;

#if 0
        if (strcmp(m->properties->sym->name,"MinX")==0)
        {
          printf("Just created a MinX at %.3e %.3e %.3e\n",m->pos.x * world->length_unit,m->pos.y * world->length_unit, m->pos.z*world->length_unit);
          printf("Hit by a particle at %.3e %.3e %.3e\n",reac_m->pos.x * world->length_unit,reac_m->pos.y*world->length_unit,reac_m->pos.z*world->length_unit);
        }
#endif
        
        m->subvol = find_subvolume(&(m->pos),m->subvol);
        m->next_v = m->subvol->mol_head;
        m->subvol->mol_head = m;
        m->subvol->mol_count++;
      }
    }
  }

  return blocked;
}


/*************************************************************************
outcome_unimolecular:
  In: the reaction that is occuring
      the path that the reaction is taking
      the molecule that is taking that path
      time that the reaction is occurring
  Out: 0 if molecule no longer exists. 1 if it does. 
       Products are created as needed.
*************************************************************************/

int outcome_unimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reac,double t)
{
  struct species *who_am_i;
  struct species *who_was_i = reac->properties;
  int blocked = 0;
  
  if ((reac->properties->flags & (ON_GRID | ON_SURFACE)) == 0)
  {
    struct molecule *m = (struct molecule*)reac;
    blocked = outcome_products(NULL,m,NULL,NULL,rx,path,m->subvol->mem,
                               0,0,t,NULL,reac,NULL,NULL);
  }
  else if ((reac->properties->flags & ON_GRID) != 0)
  {
    struct grid_molecule *g = (struct grid_molecule*) reac;
    blocked = outcome_products(g->grid->surface,NULL,NULL,g,rx,path,
                               g->grid->subvol->mem,
                               g->orient,0,t,NULL,reac,NULL,NULL);
  }
  else if ((reac->properties->flags & ON_SURFACE) != 0)
  {
    struct surface_molecule *s = (struct surface_molecule*)reac;
    blocked = outcome_products(s->curr_wall,NULL,s,NULL,rx,path,s->subvol->mem,
                               s->orient,0,t,NULL,reac,NULL,NULL);
  }
  
  if (blocked != -2) rx->counter[path]++;
  else printf("Huh?  Blocked in unimolecular?!\n");

  who_am_i = rx->players[rx->product_idx[path]];
  
  if (who_am_i == NULL)
  {
    if ((reac->properties->flags & ON_GRID) != 0)
    {
      struct grid_molecule *g = (struct grid_molecule*) reac;
      
      if (g->grid->mol[g->grid_index]==g) g->grid->mol[ g->grid_index ] = NULL;
      g->grid->n_occupied--;
    }
    else if ((reac->properties->flags & NOT_FREE)==0)
    {
      struct molecule *m = (struct molecule*)reac;
      m->subvol->mol_count--;
    }

    if ((reac->properties->flags & COUNT_CONTENTS) != 0)
      count_me_by_region(reac,-1);
    
    reac->properties->population--;
    reac->properties = NULL;
    return 0;
  }
  else if (who_am_i != who_was_i) return 0;
  else return 1;
}


/*************************************************************************
outcome_bimolecular:
  In: reaction that's occurring
      path the reaction's taking
      two molecules that are reacting (first is moving one)
      orientations of the two molecules
      time that the reaction is occurring
  Out: 0 of moving molecule no longer exists, 1 if it does,
       -1 if it's been transported across a membrane, and
       -2 if reaction failed due to lack of space for products.
       Products are created as needed.
*************************************************************************/

int outcome_bimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB,double t,struct vector3 *hitpt)
{
  struct grid_molecule *g = NULL;
  struct surface_molecule *s = NULL;
  struct molecule *m = NULL;
  struct wall *w = NULL;
  struct storage *x;
  int blocked;
  
  if ((reacA->properties->flags & (ON_GRID | ON_SURFACE)) == 0)
  {
    m = (struct molecule*) reacA;
    x = m->subvol->mem;
    if ((reacB->properties->flags & ON_SURFACE) != 0)
    {
      s = (struct surface_molecule*)reacB;
      w = s->curr_wall;
    }
    else if ((reacB->properties->flags & ON_GRID) != 0)
    {
      g = (struct grid_molecule*)reacB;
      w = g->grid->surface;
    }
  }
  else if ( (reacA->properties->flags & ON_GRID) == 0 )
  {
    s = (struct surface_molecule*)reacA;
    x = s->subvol->mem;
    w = s->curr_wall;
    
    if ((reacB->properties->flags & ON_GRID) != 0)
    {
      g = (struct grid_molecule*)reacB;
    }
    else if ((reacB->properties->flags & (ON_GRID | ON_SURFACE)) == 0)
    {
      m = (struct molecule*)reacB;
    }
  }
  else /* Grid molecule */
  {
    g = (struct grid_molecule*)reacA;
    x = g->grid->surface->birthplace;
    w = g->grid->surface;
    
    if ((reacB->properties->flags & ON_SURFACE) == 0)
    {
      m = (struct molecule*)reacB;
    }
    else /* Must be surface */
    {
      s = (struct surface_molecule*)reacA;
    }
  }
  
  blocked = outcome_products(w,m,s,g,rx,path,x,orientA,orientB,t,hitpt,reacA,reacB,reacA);
  
  if (blocked==-2) return blocked;

  rx->counter[path]++;
  
  if (rx->players[0]==reacA->properties)
  {
    if (rx->players[ rx->product_idx[path]+1 ] == NULL)
    {
      if ((reacB->properties->flags & ON_GRID) != 0)
      {
        g = (struct grid_molecule*)reacB;

        if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
        g->grid->n_occupied--;
      }
      else if ((reacB->properties->flags & NOT_FREE) == 0)
      {
        m = (struct molecule*)reacB;
        m->subvol->mol_count--;
      }

      if ((reacB->properties->flags & COUNT_CONTENTS) != 0)
        count_me_by_region(reacB,-1);
    
      reacB->properties->population--;
      reacB->properties = NULL;
    }
    if (rx->players[ rx->product_idx[path] ] == NULL)
    {
      if ((reacA->properties->flags & ON_GRID) != 0)
      {
        g = (struct grid_molecule*)reacA;

        if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
        g->grid->n_occupied--;
      }
      else if ((reacA->properties->flags & NOT_FREE) == 0)
      {
        m = (struct molecule*)reacA;
        m->subvol->mol_count--;
      }

      if ((reacA->properties->flags & COUNT_CONTENTS) != 0)
        count_me_by_region(reacA,-1);
    
      reacA->properties->population--;
      reacA->properties = NULL;
      return 0;
    }
  }
  else
  {
    if (rx->players[ rx->product_idx[path] ] == NULL)
    {
      if ((reacB->properties->flags & ON_GRID) != 0)
      {
        g = (struct grid_molecule*)reacB;

        if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
        g->grid->n_occupied--;
      }
      else if ((reacB->properties->flags & NOT_FREE) == 0)
      {
        m = (struct molecule*)reacB;
        m->subvol->mol_count--;
      }

      if ((reacB->properties->flags & COUNT_CONTENTS) != 0)
        count_me_by_region(reacB,-1);
    
      reacB->properties->population--;
      reacB->properties = NULL;
    }
    if (rx->players[ rx->product_idx[path]+1 ] == NULL)
    {
      if ((reacA->properties->flags & ON_GRID) != 0)
      {
        g = (struct grid_molecule*)reacA;

        if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
        g->grid->n_occupied--;
      }
      else if ((reacA->properties->flags & NOT_FREE) == 0)
      {
        m = (struct molecule*)reacA;
        m->subvol->mol_count--;
      }

      if ((reacA->properties->flags & COUNT_CONTENTS) != 0)
        count_me_by_region(reacA,-1);
    
      reacA->properties->population--;
      reacA->properties = NULL;
      return 0;
    }
  }
  
  return 1;
}


/*************************************************************************
outcome_intersect:
  In: reaction that's taking place
      path the reaction's taking
      wall that is being struck
      molecule that is hitting the wall
      orientation of the molecule
      time that the reaction is occurring
  Out: -1 if the molecule reflects
       1 if the molecule passes through
       0 if the molecule stops, is destroyed, changes identity, etc..
       Additionally, products are created as needed.
  Note: Can assume molecule is always first in the reaction.
*************************************************************************/

int outcome_intersect(struct rxn *rx, int path, struct wall *surface,
  struct abstract_molecule *reac,short orient,double t,struct vector3 *hitpt)
{
  int blocked,index;
  
  if (rx->n_pathways <= RX_SPECIAL) return 1;

  index = rx->product_idx[path];

  if ((reac->properties->flags & (ON_GRID | ON_SURFACE)) == 0)
  {
    struct molecule *m = (struct molecule*) reac;
    
    blocked = outcome_products(surface,m,NULL,NULL,rx,path,m->subvol->mem,orient,0,t,hitpt,reac,NULL,reac);

    if (blocked == -2) return -1;

    rx->counter[path]++;
    
    if (rx->players[ index ] == NULL)
    {
      m->subvol->mol_count--;
      if ( (reac->properties->flags & COUNT_CONTENTS) != 0 )
	count_me_by_region(reac,-1);
      reac->properties->population--;
      reac->properties = NULL;
      return 0;
    }
    else if (blocked==0) return -1;
    else return blocked;

  }
  else /* Grid can't intersect, so this must be a surface molecule */
  {
    /* TODO: fix this to look like 3D molecule case. */
    struct surface_molecule *s = (struct surface_molecule*) reac;

    if (index+2==rx->product_idx[path+1] &&
        rx->geometries[ rx->product_idx[path] ] == 1)
    {
      return -1;
    }
    else
    {
      blocked = outcome_products(surface,NULL,s,NULL,rx,path,s->subvol->mem,orient,0,t,hitpt,reac,NULL,reac);

      if (blocked == -2) return -1;

      rx->counter[path]++;

      if (rx->players[ rx->product_idx[path] ] == NULL)
      {
        reac->properties->population--;
        reac->properties = NULL;
        return 0;
      }
      else return blocked;
    }
  }
}

