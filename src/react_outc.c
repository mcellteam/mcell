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

#include "rng.h"
#include "grid_util.h"
#include "mcell_structs.h"
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
   Out: Zero if successful, one if a grid reaction was blocked by a full
        target grid tile.
        Creates products as necessary and schedules them.
*************************************************************************/

int outcome_products(struct wall *w,struct molecule *reac_m,
  struct surface_molecule *reac_s,struct grid_molecule *reac_g,
  struct rxn *rx,int path,struct storage *local,
  short orientA,short orientB,double t,struct vector3 *hitpt)
{
  int i,i0,iN;
  struct molecule *m;
  struct surface_molecule *s;
  struct grid_molecule *g;
  struct species *p;
  struct surface_grid *sg;
  int blocked = 1;

  i0 = rx->product_idx[path];
  iN = rx->product_idx[path+1];
  
  if (iN > i0)
  {
    u_int bits;
    int j;
    double f;
    struct abstract_molecule *plist[iN-i0];
    char ptype[iN-i0];
    short porient[iN-i0];

    for (i=i0;i<iN;i++)
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
              if (sg->mol[j]!=NULL) j = -1; /* Slot full, no rxn */
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
              printf("Hit @ [%.2f %.2f %.2f] and got %x:%d\n",
                     hitpt->x,hitpt->y,hitpt->z,(int)sg,j);
              if (sg->mol[j]!=NULL) j = -1; /* Slot full, no rxn */
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
            schedule_add(local->timer,g);
            
            blocked = 0;
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
        if ( reac_s != NULL || reac_g!=NULL || (reac_m!=NULL && w!=NULL))
        {
          s = mem_get(local->smol);
          s->birthplace = local->smol;
          s->properties = p;
          p->population++;
          s->flags = TYPE_SURF + ACT_NEWBIE + IN_SURFACE + IN_VOLUME + IN_SCHEDULE;
          if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)s) != NULL)
            s->flags += ACT_REACT;
          
          p->population++;
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
          
          blocked = 0;
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
        m->flags = TYPE_3D + ACT_NEWBIE + IN_VOLUME + IN_SCHEDULE;
        if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)m) != NULL)
          m->flags += ACT_REACT;
        if (p->space_step > 0.0) m->flags += ACT_DIFFUSE;
        
        if (reac_m != NULL)
        {
          m->pos.x = reac_m->pos.x;
          m->pos.y = reac_m->pos.y;
          m->pos.z = reac_m->pos.z;
          m->subvol = reac_m->subvol;
          m->next_v = reac_m->next_v;
          reac_m->next_v = m;
        }
        else if (reac_s != NULL)
        {
          m->pos.x = reac_s->pos.x;
          m->pos.y = reac_s->pos.y;
          m->pos.z = reac_s->pos.z;
          m->subvol = reac_s->subvol;
          m->next_v = m->subvol->mol_head;
          m->subvol->mol_head = m;
        }
        else if (reac_g != NULL)
        {
          grid2xyz(reac_g->grid , reac_g->grid_index , &(m->pos));
          m->subvol = find_subvolume(&(m->pos),reac_g->grid->subvol);
          m->next_v = m->subvol->mol_head;
          m->subvol->mol_head = m;
        }
        plist[i-i0] = (struct abstract_molecule*)m;
        ptype[i-i0] = 'm';
        m->t2 = 0.0;
        m->t = t;
        schedule_add( local->timer , m );
        
        blocked = 0;
      }
    }
    
    bits = rng_uint( world->seed++ );
    for (i=iN-1;i>=i0;i--,bits>>=1)
    {
      if ( ptype[i-i0] != 0 && (ptype[i-i0]!='m' || w!=NULL) )
      {
        if (rx->geometries[i] == 0)
        {
          if ((bits&1)==0) porient[i-i0] = 1;
          else porient[i-i0] = -1;
        }
        else
        {
          for (j=0;j<rx->n_reactants;j++)
          {
            if ( (rx->geometries[i]+rx->geometries[j])*
                 (rx->geometries[i]-rx->geometries[j]) == 0 ) break;
          }
          if (j>=rx->n_reactants)
          {
            for (j=iN-1;j>i;j--)
            {
              if ( (rx->geometries[i]+rx->geometries[j])*
                   (rx->geometries[i]-rx->geometries[j]) == 0 ) break;
            }
            if (j <= i)
            {
              if ((bits&1)==0) porient[i-i0] = 1;
              else porient[i-i0] = -1;
            }
            else
            {
              if (rx->geometries[j]+rx->geometries[i] == 0)
                porient[i-i0] = porient[j-i0];
              else porient[i-i0] = -porient[j-i0];
            }
          }
          else
          {
            if (rx->geometries[j]+rx->geometries[i] == 0)
            {
              if (j==0) porient[i-i0] = orientA;
              else if (j==2) porient[i-i0] = 1;
              else if (orientB!=0) porient[i-i0] = orientB;
              else porient[i-i0] = 1;
            }
            else
            {
              if (j==0) porient[i-i0] = -orientA;
              else if (j==2) porient[i-i0] = -1;
              else if (orientB!=0) porient[i-i0] = -orientB;
              else porient[i-i0] = -1;
            }
          }
        }
        
        if (ptype[i-i0]=='s')
        {
          ((struct surface_molecule*)plist[i-i0])->orient = porient[i-i0];
        }
        else if (ptype[i-i0]=='g')
        {
          ((struct grid_molecule*)plist[i-i0])->orient = porient[i-i0];
        }
        else
        {
          m = (struct molecule*)plist[i-i0];
          if (porient[i-i0]>0) f = -EPS_C;
          else f = EPS_C;
          
          m->pos.x += f*w->normal.x;
          m->pos.y += f*w->normal.y;
          m->pos.z += f*w->normal.z;
          
          printf("Creating a NEW MOLECULE at %.3e %.3e %.3e [%d]\n",m->pos.x,m->pos.y,m->pos.z,porient[i-i0]);
        }
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
  if ((reac->properties->flags & (ON_GRID | ON_SURFACE)) == 0)
  {
    struct molecule *m = (struct molecule*)reac;
    outcome_products(NULL,m,NULL,NULL,rx,path,m->subvol->mem,0,0,t,NULL);
  }
  else if ((reac->properties->flags & ON_GRID) != 0)
  {
    struct grid_molecule *g = (struct grid_molecule*) reac;
    outcome_products(g->grid->surface,NULL,NULL,g,rx,path,g->grid->subvol->mem,g->orient,0,t,NULL);
  }
  else if ((reac->properties->flags & ON_SURFACE) != 0)
  {
    struct surface_molecule *s = (struct surface_molecule*)reac;
    outcome_products(s->curr_wall,NULL,s,NULL,rx,path,s->subvol->mem,s->orient,0,t,NULL);
  }

  if ((rx->fates[path] & RX_DESTROY) != 0)
  {
    if ((reac->properties->flags & ON_GRID) != 0)
    {
      struct grid_molecule *g = (struct grid_molecule*) reac;
      if (g->grid->mol[g->grid_index]==g) g->grid->mol[ g->grid_index ] = NULL;
      g->grid->n_occupied--;
    }
    reac->properties->population--;
    reac->properties = NULL;
    return 0;
  }
  else if ((rx->fates[path] & RX_PROD) != 0) return 0;
  else return 1;
}


/*************************************************************************
outcome_bimolecular:
  In: reaction that's occurring
      path the reaction's taking
      two molecules that are reacting (first is moving one)
      orientations of the two molecules
      time that the reaction is occurring
  Out: 0 of moving molecule no longer exists, 1 if it does, -1 if reaction
       failed due to lack of space for products.
       Products are created as needed.
*************************************************************************/

int outcome_bimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB,double t)
{
  struct grid_molecule *g = NULL;
  struct surface_molecule *s = NULL;
  struct molecule *m = NULL;
  struct wall *w = NULL;
  struct storage *x;
  struct vector3 *hitpt;
  int blocked;
  
  if ((reacA->properties->flags & (ON_GRID | ON_SURFACE)) == 0)
  {
    m = (struct molecule*) reacA;
    hitpt = &(m->pos);
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
      hitpt = &(s->pos);
    }
    else if ((reacB->properties->flags & (ON_GRID | ON_SURFACE)) == 0)
    {
      m = (struct molecule*)reacB;
      hitpt = &(m->pos);
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
      hitpt = &(m->pos);
    }
    else /* Must be surface */
    {
      s = (struct surface_molecule*)reacA;
      hitpt = &(s->pos);
    }
  }
  
  blocked = outcome_products(w,m,s,g,rx,path,x,orientA,orientB,t,hitpt);
  
  if (blocked) return -1;
  
  if ((rx->fates[path] & RX_2DESTROY) != 0)
  {
    if ((reacB->properties->flags & ON_GRID) != 0)
    {
      g = (struct grid_molecule*)reacB;
      if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
      g->grid->n_occupied--;
    }
    reacB->properties->population--;
    reacB->properties = NULL;
  }
  else if ((rx->fates[path] & RX_2FLIP) != 0)
  {
    if ((reacB->properties->flags & ON_GRID) != 0)
      ((struct grid_molecule*)reacB)->orient *= -1;
    else if ((reacB->properties->flags & ON_SURFACE) != 0)
      ((struct surface_molecule*)reacB)->orient *= -1;
    else if (w != NULL)
    {
      double f;
      m = (struct molecule*)reacB;
      f = m->pos.x*w->normal.x + m->pos.y*w->normal.y + m->pos.z*w->normal.z;
      f *= 2.0;
      m->pos.x -= f*w->normal.x;
      m->pos.y -= f*w->normal.y;
      m->pos.z -= f*w->normal.z;      
    }
  }
  if ((rx->fates[path] & RX_DESTROY) != 0)
  {
    if ((reacA->properties->flags & ON_GRID) != 0)
    {
      g = (struct grid_molecule*)reacA;
      if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
      g->grid->n_occupied--;
    }
    reacA->properties->population--;
    reacA->properties = NULL;
    return 0;
  }
  else if ((rx->fates[path] & RX_FLIP) != 0)
  {
    if ((reacA->properties->flags & ON_GRID) != 0)
      ((struct grid_molecule*)reacA)->orient *= -1;
    else if ((reacA->properties->flags & ON_SURFACE) != 0)
      ((struct surface_molecule*)reacA)->orient *= -1;
    else if (w != NULL)
    {
      double f;
      m = (struct molecule*)reacA;
      f = m->pos.x*w->normal.x + m->pos.y*w->normal.y + m->pos.z*w->normal.z;
      f *= 2.0;
      m->pos.x -= f*w->normal.x;
      m->pos.y -= f*w->normal.y;
      m->pos.z -= f*w->normal.z;      
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
*************************************************************************/

int outcome_intersect(struct rxn *rx, int path, struct wall *surface,
  struct abstract_molecule *reac,short orient,double t,struct vector3 *hitpt)
{
  int blocked;
  
  if ((reac->properties->flags & (ON_GRID | ON_SURFACE)) == 0)
  {
    struct molecule *m = (struct molecule*) reac;

    if (rx->product_idx[path]==rx->product_idx[path+1] &&
        (rx->fates[path] == RX_REFL))
    {
      return -1;
    }
    else
    {
      blocked = outcome_products(surface,m,NULL,NULL,rx,path,m->subvol->mem,orient,0,t,hitpt);
      if ((rx->fates[path] & RX_REFL) != 0) return -1;
      else if ((rx->fates[path] & RX_FLIP) != 0) return 1;
      else if (blocked) return -1;
      else if ((rx->fates[path] & RX_DESTROY) != 0)
      {
        reac->properties->population--;
        reac->properties = NULL;
        return 0;
      }
      else return 0; 
    }
  }
  else /* Grid can't intersect, so this must be a surface molecule */
  {
    struct surface_molecule *s = (struct surface_molecule*) reac;

    if (rx->product_idx[path]==rx->product_idx[path+1] &&
        (rx->fates[path] == RX_REFL))
    {
      return -1;
    }
    else
    {
      blocked = outcome_products(surface,NULL,s,NULL,rx,path,s->subvol->mem,orient,0,t,hitpt);
      if ((rx->fates[path] & RX_REFL) != 0) return -1;
      else if (blocked || (rx->fates[path] & RX_DESTROY) == 0) return 1;
      else
      {
        reac->properties->population--;
        reac->properties = NULL;
        return 0;
      }
    }
  }
}

