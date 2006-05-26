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


int is_compatible_surface(void *req_species,struct wall *w)
{
  struct species *rs = (struct species*)req_species;
  
  if (rs==NULL) return 1;
  
  return (w->surf_class == rs);
}


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
          RX_BLOCKED reaction blocked by full grid
          RX_FLIP moving molecule passed through membrane
          RX_A_OK everything went fine, nothing extra to do
	  RX_NO_MEM out of memory error
        Products are created as necessary and scheduled.
*************************************************************************/

int outcome_products(struct wall *w,struct molecule *reac_m,
  struct grid_molecule *reac_g,struct rxn *rx,int path,struct storage *local,
  short orientA,short orientB,double t,struct vector3 *hitpt,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  struct abstract_molecule *moving)
{
  int bounce = RX_A_OK;

  int i;
  struct molecule *m;
  struct grid_molecule *g;
  struct species *p;
  struct surface_grid *sg;
  u_int bits;
  int j,k;
  double f;
  int i0 = rx->product_idx[path];
  int iN = rx->product_idx[path+1];
  int replace_p1 = 0;
  int replace_p2 = 0;
  struct vector2 uv_loc;
  
  struct abstract_molecule *plist[iN-i0];
  char ptype[iN-i0];
  short porient[iN-i0];
  struct surface_grid *glist[iN-i0];
  int xlist[iN-i0];
  byte flist[iN-i0];
  struct grid_molecule fake;
  int fake_idx = -1;
  
#define FLAG_NOT_SET 0
#define FLAG_USE_UV_LOC 1
#define FLAG_USE_REACA_UV 2
#define FLAG_USE_REACB_UV 3
#define FLAG_USE_RANDOM 4
  
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
  else if ( (reacA->properties->flags&NOT_FREE)==0 ) ptype[0] = 'm';
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
      else if ( (reacB->properties->flags&NOT_FREE)==0 ) ptype[1] = 'm';
      else ptype[1] = '!';
    }
  }
  
  /* Make sure there's space for the reaction to occur */
  /* FIXME--could speed this up with some pre-computation of reactions to at least see if we need to bother */
  k = -1;
  if (ptype[0]=='g' && rx->players[i0]==NULL) replace_p1=1;
  if (rx->n_reactants>1 && ptype[1]=='g' && rx->players[i0+1]==NULL) replace_p2=1;
  if (reac_g!=NULL || (reac_m!=NULL && w!=NULL))  /* Surface involved */
  {
    if (reac_g!=NULL) memcpy(&uv_loc , &(reac_g->s_pos) , sizeof(struct vector2));
    else xyz2uv(hitpt,w,&uv_loc);
    
    for (j=i0+rx->n_reactants;j<iN;j++)
    {
      if (rx->players[j]->flags&ON_GRID)
      {
	/* FIXME: decide on a policy for who is replaced first in grid-grid reactions */
	if (replace_p1)
	{
	  glist[j - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacA)->grid;
	  xlist[j - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacA)->grid_index;
	  flist[j - (i0+rx->n_reactants)] = FLAG_USE_REACA_UV;
	  replace_p1=0;
	  continue;
	}
	else if (replace_p2)
	{
	  glist[j - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacB)->grid;
	  xlist[j - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacB)->grid_index;
	  flist[j - (i0+rx->n_reactants)] = FLAG_USE_REACB_UV;
	  replace_p2=0;
	  continue;
	}
	else if (w->effectors==NULL)
	{
	  if (create_grid(w,reac_m->subvol)) return RX_NO_MEM;  /* No effectors means it must be a 3d mol/surface rx */
	  fake_idx = j - (i0+rx->n_reactants);
	  glist[fake_idx] = w->effectors;
	  xlist[fake_idx] = uv2grid(&uv_loc,w->effectors);
	  flist[fake_idx] = FLAG_USE_UV_LOC;
	  continue;
	}
	else
	{
	  struct wall *temp_w;
	  
	  if (fake_idx > -1) glist[fake_idx]->mol[ xlist[fake_idx] ] = &fake; /* Assumed empty! */

          fake_idx = j - (i0+rx->n_reactants);
	  if (k==-1)
	  {
	    k = uv2grid(&uv_loc,w->effectors);
	    if (w->effectors->mol[k]==NULL)
	    {
	      glist[fake_idx] = w->effectors;
	      xlist[fake_idx] = k;
	      flist[fake_idx] = FLAG_USE_UV_LOC;
	      continue;
	    }
	  }
	  
	  if (world->vacancy_search_dist2 > 0)
	  {
	    /* FIXME: only do this for surfaces compatible with the reaction! */
    	    temp_w = search_nbhd_for_free(w,&uv_loc,world->vacancy_search_dist2,&k,&is_compatible_surface,NULL);
            
	    if (temp_w != NULL)
	    {
	      glist[fake_idx] = temp_w->effectors;
	      xlist[fake_idx] = k;
	      flist[fake_idx] = FLAG_USE_RANDOM;
	      continue;
	    }
	  }
	  
	  /* Uh-oh--if we get to this point and we haven't found space, we're blocked */
	  for (k=0;k<fake_idx;k++)
	  {
	    if (glist[k]==NULL) continue;
	    if (glist[k]->mol[ xlist[k] ] == &fake) glist[k]->mol[xlist[k]]=NULL; /* Remove sentinels */
	  }
	  return RX_BLOCKED;
	}
      }
      else
      {
	glist[j - (i0+rx->n_reactants)]=NULL;
	xlist[j - (i0+rx->n_reactants)]=-1;
	flist[j - (i0+rx->n_reactants)]=FLAG_NOT_SET;
      }
    }
  }

  /* We know there's space, so now actually create everyone */
  for (i=i0+rx->n_reactants;i<iN;i++)
  {
    p = rx->players[i];
  
    if ( (p->flags & ON_GRID) != 0 )
    {
      if (reac_g!=NULL || (reac_m!=NULL && w!=NULL))
      {
	k = i-(i0+rx->n_reactants);
	
	g = mem_get(local->gmol);
	if (g==NULL) return RX_NO_MEM;
	g->birthplace = local->gmol;
	g->birthday = t;
	g->properties = p;
	p->population++;
	g->flags = TYPE_GRID | ACT_NEWBIE | IN_SCHEDULE;
	if (p->space_step>0) g->flags |= ACT_DIFFUSE;
	if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)g)!= NULL || (p->flags&CAN_GRIDWALL)!=0) g->flags |= ACT_REACT;
	
	g->t = t;
	g->t2 = 0.0;
	sg = g->grid = glist[k];
	j = g->grid_index = xlist[k];
	
	if (world->randomize_gmol_pos)
	{
	  switch (flist[k])
	  {
	    case FLAG_USE_REACA_UV:
	      memcpy(&(g->s_pos),&(((struct grid_molecule*)reacA)->s_pos),sizeof(struct vector2));
	      break;
	    case FLAG_USE_REACB_UV:
	      memcpy(&(g->s_pos),&(((struct grid_molecule*)reacB)->s_pos),sizeof(struct vector2));
	      break;
	    case FLAG_USE_UV_LOC:
	      memcpy(&(g->s_pos),&(uv_loc),sizeof(struct vector2));
	      break;
	    case FLAG_USE_RANDOM:
	      grid2uv_random(glist[k],xlist[k],&(g->s_pos));
	      break;
	    default:
	      fprintf(world->err_file,"Screwed up surface molecule placement badly!\n  Aborting execution (guessing out of memory?).\n");
	      return RX_NO_MEM;
	      break;
	  }
	}
	else grid2uv(sg,j,&(g->s_pos));

	if (reac_g==NULL || sg->mol[j]!=reac_g) sg->n_occupied++;
	sg->mol[j] = g;
	
	plist[i-i0] = (struct abstract_molecule*)g;
	ptype[i-i0] = 'g';
	
	if ( schedule_add(local->timer,g) ) return RX_NO_MEM;
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
      if (m==NULL) return RX_NO_MEM;
      m->birthplace = local->mol;
      m->birthday = t;
      m->properties = p;
      p->population++;
      if (reac_g != NULL)
      {
        m->previous_wall = reac_g->grid->surface;
        m->index = reac_g->grid_index;
      }
      else
      {
        m->previous_wall = NULL;
        m->index = -1;
      }
      m->flags = TYPE_3D + ACT_NEWBIE + IN_VOLUME + IN_SCHEDULE;
      if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)m) != NULL) m->flags |= ACT_REACT;
      if (p->space_step > 0.0) m->flags += ACT_DIFFUSE;
      
      if (hitpt != NULL)
      {
        m->pos.x = hitpt->x;
        m->pos.y = hitpt->y;
        m->pos.z = hitpt->z;
      }
      
      if (reac_m != NULL)
      {
        if (hitpt==NULL || (struct abstract_molecule*)reac_m != moving)
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
      else if (reac_g != NULL)
      {
        if (hitpt==NULL) uv2xyz(&(reac_g->s_pos) , reac_g->grid->surface , &(m->pos));
        m->subvol = find_subvolume(&(m->pos),reac_g->grid->subvol);
      }
      plist[i-i0] = (struct abstract_molecule*)m;
      ptype[i-i0] = 'm';
      m->t2 = 0.0;
      m->t = t;

      if ( schedule_add( local->timer , m ) ) return RX_NO_MEM;
      
    }
  }
  
/*
  if (rx->info[path].pathname != NULL)
  {
    count_rx_by_region(reacA,w,rx->info[path]->pathname,1);
  }
*/
  
  /* Finally, set orientations correctly */
  bits = rng_uint( world->rng );
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

      if (ptype[i-i0]=='g')
      {
        ((struct grid_molecule*)plist[i-i0])->orient = porient[i-i0];
      }
      else if (moving == plist[i-i0])
      {
        if (moving==reacA)
        {
          if (orientA==porient[i-i0]) bounce = RX_A_OK;
          else bounce = RX_FLIP;
        }
        else
        {
          if (orientB==porient[i-i0]) bounce = RX_A_OK;
          else bounce = RX_FLIP;
        }
      }
      else if (ptype[i-i0]=='m')
      {
        m = (struct molecule*)plist[i-i0];
        if (porient[i-i0]>0) f = EPS_C;
        else f = -EPS_C;
        	
        m->pos.x += f*w->normal.x;
        m->pos.y += f*w->normal.y;
        m->pos.z += f*w->normal.z;

        m->subvol = find_subvolume(&(m->pos),m->subvol);
        m->next_v = m->subvol->mol_head;
        m->subvol->mol_head = m;
        m->subvol->mol_count++;
      }
    }
    
    if (i >= i0+rx->n_reactants &&
        (plist[i-i0]->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
    {
      count_me_by_region(plist[i-i0],1,NULL);
    }
  }
  
  return bounce;
  
#undef FLAG_NOT_SET
#undef FLAG_USE_UV_LOC
#undef FLAG_USE_REACA_UV
#undef FLAG_USE_REACB_UV
#undef FLAG_USE_RANDOM
}


/*************************************************************************
outcome_unimolecular:
  In: the reaction that is occuring
      the path that the reaction is taking
      the molecule that is taking that path
      time that the reaction is occurring
  Out: Value based on outcome:
	 RX_DESTROY if molecule no longer exists.
	 RX_A_OK if it does.
	 RX_NO_MEM on an out-of-memory error.
       Products are created as needed.
*************************************************************************/

int outcome_unimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reac,double t)
{
  struct species *who_am_i;
  struct species *who_was_i = reac->properties;
  int result = RX_A_OK;
  
  if ((reac->properties->flags & NOT_FREE) == 0)
  {
    struct molecule *m = (struct molecule*)reac;
    result = outcome_products(NULL,m,NULL,rx,path,m->subvol->local_storage,
                              0,0,t,NULL,reac,NULL,NULL);
  }
  else
  {
    struct grid_molecule *g = (struct grid_molecule*) reac;
    result = outcome_products(g->grid->surface,NULL,g,rx,path,
                              g->grid->subvol->local_storage,
                              g->orient,0,t,NULL,reac,NULL,NULL);
  }
  
  if (result==RX_NO_MEM) return RX_NO_MEM;
  
  if (result != RX_BLOCKED) rx->info[path].count++;

  who_am_i = rx->players[rx->product_idx[path]];
  
  if (who_am_i == NULL)
  {
    if ((reac->properties->flags & NOT_FREE)==0)
    {
      struct molecule *m = (struct molecule*)reac;
      m->subvol->mol_count--;
      if (m->flags & IN_SCHEDULE)
      {
	m->subvol->local_storage->timer->defunct_count++;
      }
    }
    else
    {
      struct grid_molecule *g = (struct grid_molecule*) reac;
      
      if (g->grid->mol[g->grid_index]==g) g->grid->mol[ g->grid_index ] = NULL;
      g->grid->n_occupied--;
      if (g->flags & IN_SCHEDULE)
      {
	g->grid->subvol->local_storage->timer->defunct_count++;
      }
    }

    if ((reac->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
    {
      count_me_by_region(reac,-1,NULL);
    }
    
    reac->properties->n_deceased++;
    reac->properties->cum_lifetime += t - reac->birthday;
    reac->properties->population--;
    reac->properties = NULL;
    return RX_DESTROY;
  }
  else if (who_am_i != who_was_i) return RX_DESTROY;
  else return result;
}


/*************************************************************************
outcome_bimolecular:
  In: reaction that's occurring
      path the reaction's taking
      two molecules that are reacting (first is moving one)
      orientations of the two molecules
      time that the reaction is occurring
  Out: Value based on outcome:
	 RX_BLOCKED if there was no room to put products on grid
	 RX_FLIP if the molecule goes across the membrane
	 RX_DESTROY if the molecule no longer exists
	 RX_A_OK if everything proceeded smoothly
	 RX_NO_MEM on an out-of-memory error
       Products are created as needed.
*************************************************************************/

int outcome_bimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB,double t,struct vector3 *hitpt)
{
  struct grid_molecule *g = NULL;
  struct molecule *m = NULL;
  struct wall *w = NULL;
  struct storage *x;
  int result;
  
  if ((reacA->properties->flags & NOT_FREE) == 0)
  {
    m = (struct molecule*) reacA;
    x = m->subvol->local_storage;
    if ((reacB->properties->flags & ON_GRID) != 0)
    {
      g = (struct grid_molecule*)reacB;
      w = g->grid->surface;
    }
    else /* Prefer to use target */
    {
      m = (struct molecule*) reacB;
      x = m->subvol->local_storage;
    }
  }
  else /* Grid molecule */
  {
    g = (struct grid_molecule*)reacA;
    x = g->grid->surface->birthplace;
    w = g->grid->surface;
    
    if ((reacB->properties->flags & NOT_FREE) == 0)
    {
      m = (struct molecule*)reacB;
    }
    else /* Must be surface */
    {
      /* TODO: handle this case */
    }
  }
  
  result = outcome_products(w,m,g,rx,path,x,orientA,orientB,t,hitpt,reacA,reacB,reacA);
  
  if (result==RX_NO_MEM) return RX_NO_MEM;
  else if (result==RX_BLOCKED) return RX_BLOCKED;
  
  rx->info[path].count++;
  
  if (rx->players[0]==reacA->properties)
  {
    if (rx->players[ rx->product_idx[path]+1 ] == NULL)
    {
      if ((reacB->properties->flags & ON_GRID) != 0)
      {
        g = (struct grid_molecule*)reacB;

        if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
        g->grid->n_occupied--;
	if (g->flags&IN_SURFACE) g->flags -= IN_SURFACE;
	if (g->flags & IN_SCHEDULE)
	{
	  g->grid->subvol->local_storage->timer->defunct_count++;
	}
      }
      else if ((reacB->properties->flags & NOT_FREE) == 0)
      {
        m = (struct molecule*)reacB;
        m->subvol->mol_count--;
	if (m->flags & IN_SCHEDULE)
	{
	  m->subvol->local_storage->timer->defunct_count++;
	}
      }

      if ((reacB->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
        count_me_by_region(reacB,-1,NULL);
      
      reacB->properties->n_deceased++;
      reacB->properties->cum_lifetime += t - reacB->birthday;
      reacB->properties->population--;
      reacB->properties = NULL;
      if ((reacB->flags&IN_MASK)==0) mem_put(reacB->birthplace,reacB);
    }
    if (rx->players[ rx->product_idx[path] ] == NULL)
    {
      if ((reacA->properties->flags & ON_GRID) != 0)
      {
        g = (struct grid_molecule*)reacA;

        if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
        g->grid->n_occupied--;
	if (g->flags & IN_SCHEDULE)
	{
	  g->grid->subvol->local_storage->timer->defunct_count++;
	}
      }
      else if ((reacA->properties->flags & NOT_FREE) == 0)
      {
        m = (struct molecule*)reacA;
        m->subvol->mol_count--;
	if (m->flags & IN_SCHEDULE)
	{
	  m->subvol->local_storage->timer->defunct_count++;
	}
      }

      if ((reacA->properties->flags&(COUNT_CONTENTS|COUNT_ENCLOSED))!=0 && (reacA->flags&COUNT_ME)!=0)
	count_me_by_region(reacA,-1,NULL);
    
      reacA->properties->n_deceased++;
      reacA->properties->cum_lifetime += t - reacA->birthday;
      reacA->properties->population--;
      reacA->properties = NULL;
      return RX_DESTROY;
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
	if (g->flags&IN_SURFACE) g->flags -= IN_SURFACE;
	if (g->flags & IN_SCHEDULE)
	{
	  g->grid->subvol->local_storage->timer->defunct_count++;
	}
      }
      else if ((reacB->properties->flags & NOT_FREE) == 0)
      {
        m = (struct molecule*)reacB;
        m->subvol->mol_count--;
	if (m->flags & IN_SCHEDULE)
	{
	  m->subvol->local_storage->timer->defunct_count++;
	}
      }

      if ((reacB->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
        count_me_by_region(reacB,-1,NULL);
    
      reacB->properties->n_deceased++;
      reacB->properties->cum_lifetime += t - reacB->birthday;
      reacB->properties->population--;
      reacB->properties = NULL;
      if ((reacB->flags&IN_MASK)==0) mem_put(reacB->birthplace,reacB);
    }
    if (rx->players[ rx->product_idx[path]+1 ] == NULL)
    {
      if ((reacA->properties->flags & ON_GRID) != 0)
      {
        g = (struct grid_molecule*)reacA;

        if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
        g->grid->n_occupied--;
	if (g->flags & IN_SCHEDULE)
	{
	  g->grid->subvol->local_storage->timer->defunct_count++;
	}
      }
      else if ((reacA->properties->flags & NOT_FREE) == 0)
      {
        m = (struct molecule*)reacA;
        m->subvol->mol_count--;
	if (m->flags & IN_SCHEDULE)
	{
	  m->subvol->local_storage->timer->defunct_count++;
	}
      }

      if ((reacA->properties->flags&(COUNT_CONTENTS|COUNT_ENCLOSED))!=0 && (reacA->flags&COUNT_ME)!=0)
	count_me_by_region(reacA,-1,NULL);
    
      reacA->properties->n_deceased++;
      reacA->properties->cum_lifetime += t - reacA->birthday;
      reacA->properties->population--;
      reacA->properties = NULL;
      return RX_DESTROY;
    }
  }
  
  return result;
}


/*************************************************************************
outcome_intersect:
  In: reaction that's taking place
      path the reaction's taking
      wall that is being struck
      molecule that is hitting the wall
      orientation of the molecule
      time that the reaction is occurring
  Out: Value depending on outcome:
	 RX_A_OK if the molecule reflects
	 RX_FLIP if the molecule passes through
	 RX_DESTROY if the molecule stops, is destroyed, etc.
	 RX_NO_MEM on an out of memory error
       Additionally, products are created as needed.
  Note: Can assume molecule is always first in the reaction.
*************************************************************************/

int outcome_intersect(struct rxn *rx, int path, struct wall *surface,
  struct abstract_molecule *reac,short orient,double t,struct vector3 *hitpt)
{
  int result,index;
  
  if (rx->n_pathways <= RX_SPECIAL)
  {
    if (rx->n_pathways==RX_REFLEC) return RX_A_OK;
    else return RX_FLIP; /* Flip = transparent is default special case */
  }

  index = rx->product_idx[path];

  if ((reac->properties->flags & NOT_FREE) == 0)
  {
    struct molecule *m = (struct molecule*) reac;
    
    result = outcome_products(surface,m,NULL,rx,path,m->subvol->local_storage,orient,0,t,hitpt,reac,NULL,reac);

    if (result==RX_NO_MEM) return RX_NO_MEM;
    else if (result == RX_BLOCKED) return RX_A_OK;

    rx->info[path].count++;
    
    if (rx->players[ index ] == NULL)
    {
      m->subvol->mol_count--;
      if ((reac->properties->flags&(COUNT_CONTENTS|COUNT_ENCLOSED))!=0 && (reac->flags&COUNT_ME)!=0)
	count_me_by_region(reac,-1,NULL);
      reac->properties->n_deceased++;
      reac->properties->cum_lifetime += t - reac->birthday;
      reac->properties->population--;
      reac->properties = NULL;
      if (m->flags & IN_SCHEDULE)
      {
        m->subvol->local_storage->timer->defunct_count++;
      }
      
      return RX_DESTROY;
    }
    else return result; /* RX_A_OK or RX_FLIP */
  }
  else /* Grid can't intersect, so this must be a surface molecule */
  {
    /* TODO: write this for grid molecule case. */
    return RX_A_OK;
  }
}

