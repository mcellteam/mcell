/**************************************************************************\
** File: react_outc.c                                                     **
**                                                                        **
** Purpose: Implements specific reaction outcome pathways.                **
**                                                                        **
\**************************************************************************/

#include <string.h>
#include <math.h>

#include "rng.h"
#include "util.h"
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
   In: relevant wall in the interaction, if any
       first free molecule in the interaction, if any
       first surface molecule in the interaction, if any
       reaction that is occuring
       path that the reaction is taking
       local storage for creating new molecules
       orientations of the molecules in the reaction
       time that the reaction is occurring
       location of the reaction (may be NULL)
       the reactants
       molecule that is moving, if any
   Out: Value depending on how orientations changed--
          RX_BLOCKED reaction blocked by full grid
          RX_FLIP moving molecule passed through membrane
          RX_A_OK everything went fine, nothing extra to do
	  RX_NO_MEM out of memory error
        Products are created as necessary and scheduled.
*************************************************************************/

int outcome_products(struct wall *w,struct volume_molecule *reac_m,
  struct grid_molecule *reac_g,struct rxn *rx,int path,struct storage *local,
  short orientA,short orientB,double t,struct vector3 *hitpt,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  struct abstract_molecule *moving)
{
  int bounce = RX_A_OK;

  int i;
  struct volume_molecule *m;
  struct grid_molecule *g;
  struct species *p;
  struct surface_grid *sg;
  u_int bits;
  int j,k;
  double f;
  int i0 = rx->product_idx[path]; /*index of the first product for the pathway*/
  int iN = rx->product_idx[path+1];/*index of the last product for the pathway*/
  int replace_p1 = 0; /* flag for the product to replace position of reactant1 */
  int replace_p2 = 0; /* flag for the product to replace position of reactant2 */
  struct vector2 uv_loc;  /* where reaction happened */
  struct vector3 xyz_loc;
 
   
  struct abstract_molecule *plist[iN-i0]; /* array of products */
  /* array that decodes the type of each product */
  char ptype[iN-i0];
  /* array of orientations for each product */
  short porient[iN-i0];
  /* array of surface_grids for products (if they are grid_molecules) */
  struct surface_grid *glist[iN-i0];
  /* array of grid_indices for products (if they are grid_molecules) */
  int xlist[iN-i0];
  /* array of flags for products */
  byte flist[iN-i0];
  struct grid_molecule fake;
  int fake_idx = -1;
  
#define FLAG_NOT_SET 0
#define FLAG_USE_UV_LOC 1
#define FLAG_USE_REACA_UV 2
#define FLAG_USE_REACB_UV 3
#define FLAG_USE_RANDOM 4
 

  /* make sure that reacA correponds to rx->players[0], and
     reacB - to rx->players[1] */ 
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
      if(rx->n_reactants > 2){
         ptype[2] = 'w';
      }

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
	else if (w->grid==NULL)
	{
	  if (create_grid(w,reac_m->subvol)) return RX_NO_MEM;  /* No effectors means it must be a 3d mol/surface rx */
	  fake_idx = j - (i0+rx->n_reactants);
	  glist[fake_idx] = w->grid;
	  xlist[fake_idx] = uv2grid(&uv_loc,w->grid);
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
	    k = uv2grid(&uv_loc,w->grid);
	    if (w->grid->mol[k]==NULL)
	    {
	      glist[fake_idx] = w->grid;
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
	      glist[fake_idx] = temp_w->grid;
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
  if (hitpt!=NULL) memcpy(&xyz_loc,hitpt,sizeof(struct vector3));
  else if (reac_g!=NULL) uv2xyz(&(reac_g->s_pos),reac_g->grid->surface,&xyz_loc);
  else memcpy(&xyz_loc,&(reac_m->pos),sizeof(struct vector3));
  
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
               
        sg->n_occupied++;
                 
	sg->mol[j] = g;
	
	plist[i-i0] = (struct abstract_molecule*)g;
	ptype[i-i0] = 'g';
	
	if ( schedule_add(local->timer,g) ) return RX_NO_MEM;
      }
      else /* Should never happen, but it doesn't hurt to be safe */
      {
        plist[i-i0] = NULL;
        ptype[i-i0] = 0;
        continue;
      }
    }
    else /* volume molecule */
    {
      m = mem_get(local->mol);
      if (m==NULL) return RX_NO_MEM;
      m->birthplace = local->mol;
      m->birthday = t;
      m->properties = p;
      p->population++;

      m->flags = TYPE_3D | ACT_NEWBIE | IN_VOLUME | IN_SCHEDULE;
      if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)m) != NULL) m->flags |= ACT_REACT;
      if (p->space_step > 0.0) m->flags |= ACT_DIFFUSE;
      if (reac_g != NULL)
      {
        m->previous_wall = reac_g->grid->surface;
        m->index = reac_g->grid_index;  /* Overwrite this with orientation in CLAMPED case */
        if (world->surface_reversibility) m->flags |= ACT_CLAMPED;
      }
      else
      {
        m->previous_wall = NULL;
        m->index = -1;
      }
      
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
        
        if (w==NULL) /* place product of non-orientable reaction in volume */
        {
          m->next_v = m->subvol->mol_head;
          m->subvol->mol_head = m;
          m->subvol->mol_count++;
        }
        /* oriented case handled below after orientation is set */
      }
      else if (reac_g != NULL)
      {
        if (hitpt==NULL) uv2xyz(&(reac_g->s_pos) , reac_g->grid->surface , &(m->pos));
        m->subvol = find_subvolume(&(m->pos),reac_g->grid->subvol);
      }
      plist[i-i0] = (struct abstract_molecule*)m;
      ptype[i-i0] = 'm';
      m->t = t;
      m->t2 = 0.0;

      if ( schedule_add( local->timer , m ) ) return RX_NO_MEM;
      
    }
  }

  /* Finally, set orientations correctly */
  bits = 0;
  for (i=i0;i<iN;i++,bits>>=1)
  {
     if (rx->players[i]==NULL) continue; 
    
    /* generate 32 random bits every 32 times through this loop */
    if (((i-i0)&0x1F)==0) {
        bits = rng_uint( world->rng );
        if(world->notify->final_summary == NOTIFY_FULL){
           world->random_number_use++;
        }
    }

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
        m = (struct volume_molecule*)plist[i-i0];
        if (porient[i-i0]>0) f = EPS_C;
        else f = -EPS_C;
	
	if ((m->flags&ACT_CLAMPED) && world->surface_reversibility)
	{
	  m->index = (porient[i-i0]>0)?1:-1; /* Which direction do we move? */
	}
        	
        /* Note: no raytracing here so it is rarely possible to jump through closely spaced surfaces */
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
      j=count_region_from_scratch(plist[i-i0],NULL,1,NULL,w,t);
      if (j) return RX_NO_MEM;
    }
  }
  
  /* Handle events triggered off of named reactions */
  if (rx->info[path].pathname!=NULL)
  {
    /* No flags for reactions so we have to check regions if we have waypoints! Fix to be more efficient for WORLD-only counts? */
    if (world->place_waypoints_flag)
    {
      j=count_region_from_scratch(NULL,rx->info[path].pathname,1,&xyz_loc,w,t);
      if (j) return RX_NO_MEM;
    }
    
    /* Other magical stuff.  For now, can only trigger releases. */
    if (rx->info[path].pathname->magic!=NULL)
    {
      j=reaction_wizardry(rx->info[path].pathname->magic,w,&xyz_loc,t);
      if (j) return RX_NO_MEM;
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
  struct volume_molecule *m=NULL;
  struct grid_molecule *g=NULL;
  int i;
  
  if ((reac->properties->flags & NOT_FREE) == 0)
  {
    m = (struct volume_molecule*)reac;
    result = outcome_products(NULL,m,NULL,rx,path,m->subvol->local_storage,
                              0,0,t,NULL,reac,NULL,NULL);
  }
  else
  {
    g = (struct grid_molecule*) reac;
    result = outcome_products(g->grid->surface,NULL,g,rx,path,
                              g->grid->subvol->local_storage,
                              g->orient,0,t,NULL,reac,NULL,NULL);
  }
  
  if (result==RX_NO_MEM) return RX_NO_MEM;
  
  if (result != RX_BLOCKED) {
     rx->info[path].count++;
     rx->n_occurred++;
  }

  who_am_i = rx->players[rx->product_idx[path]];
  
  if (who_am_i == NULL)
  {
    if ((reac->properties->flags & NOT_FREE)==0)
    {
      m->subvol->mol_count--;
      if (m->flags & IN_SCHEDULE) m->subvol->local_storage->timer->defunct_count++;
      if (m->properties->flags&COUNT_SOME_MASK)
      {
	i = count_region_from_scratch((struct abstract_molecule*)m,NULL,-1,&(m->pos),NULL,m->t);
	if (i) return RX_NO_MEM;
      }
    }
    else
    {
      if (g->grid->mol[g->grid_index]==g) g->grid->mol[ g->grid_index ] = NULL;
      g->grid->n_occupied--;
      if (g->flags & IN_SCHEDULE)
      {
	g->grid->subvol->local_storage->timer->defunct_count++;
      }
      if (g->properties->flags&COUNT_SOME_MASK)
      {
	i = count_region_from_scratch((struct abstract_molecule*)g,NULL,-1,NULL,NULL,g->t);
	if (i) return RX_NO_MEM;
      }
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
      location of collision between molecules
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
  short orientA,short orientB,double t,struct vector3 *hitpt,
  struct vector3 *loc_okay)
{
  struct grid_molecule *g = NULL;
  struct volume_molecule *m = NULL;
  struct wall *w = NULL;
  struct storage *x;
  int result;
  int i;
  int reacB_was_free=0;
  int killA,killB;
  
  if ((reacA->properties->flags & NOT_FREE) == 0)
  {
    m = (struct volume_molecule*) reacA;
    x = m->subvol->local_storage;
    if ((reacB->properties->flags & ON_GRID) != 0)
    {
      g = (struct grid_molecule*)reacB;
      w = g->grid->surface;
    }
    else /* Prefer to use target */
    {
      m = (struct volume_molecule*) reacB;
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
      m = (struct volume_molecule*)reacB;
    }
  }

  result = outcome_products(w,m,g,rx,path,x,orientA,orientB,t,hitpt,reacA,reacB,reacA);
          
   
  if (result==RX_NO_MEM) return RX_NO_MEM;
  else if (result==RX_BLOCKED) return RX_BLOCKED;
  
  rx->n_occurred++;
  rx->info[path].count++;
  
  /* Figure out if either of the reactants was destroyed */
  if (rx->players[0]==reacA->properties)
  {
    killB = (rx->players[ rx->product_idx[path]+1 ] == NULL);
    killA = (rx->players[ rx->product_idx[path] ] == NULL);
  }
  else
  {
    killB = (rx->players[ rx->product_idx[path] ] == NULL);
    killA = (rx->players[ rx->product_idx[path]+1 ] == NULL);
  }
  
  if (killB)
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
      m = (struct volume_molecule*)reacB;
      m->subvol->mol_count--;
      if (m->flags & IN_SCHEDULE)
      {
	m->subvol->local_storage->timer->defunct_count++;
      }
      reacB_was_free=1;
    }

    if ((reacB->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
    {
      i = count_region_from_scratch(reacB,NULL,-1,NULL,NULL,t);
      if (i) return RX_NO_MEM;
    }
    
    reacB->properties->n_deceased++;
    reacB->properties->cum_lifetime += t - reacB->birthday;
    reacB->properties->population--;
    reacB->properties = NULL;
    if ((reacB->flags&IN_MASK)==0) mem_put(reacB->birthplace,reacB);
  }

  if (killA)
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
      m = (struct volume_molecule*)reacA;
      m->subvol->mol_count--;
      if (m->flags & IN_SCHEDULE)
      {
	m->subvol->local_storage->timer->defunct_count++;
      }
    }

    if ((reacA->properties->flags&ON_GRID)!=0)  /* Grid molecule is OK where it is, doesn't obey COUNT_ME */
    {
      if (reacA->properties->flags&COUNT_SOME_MASK)  /* If we're ever counted, try to count us now */
      {
        i=count_region_from_scratch(reacA,NULL,-1,NULL,NULL,t);	  
      }
    }
    else if (reacA->flags&COUNT_ME)
    {
      /* Subtlety: we made it up to hitpt, but our position is wherever we were before that! */
      if (hitpt==NULL || reacB_was_free || (reacB->properties!=NULL && (reacB->properties->flags&NOT_FREE)!=0))
      {
	/* Vol-vol rx should be counted at hitpt */
	i=count_region_from_scratch(reacA,NULL,-1,hitpt,NULL,t);
      }
      else /* Vol-surf but don't want to count exactly on a wall or we might count on the wrong side */
      {
	struct vector3 fake_hitpt;
	
	m = (struct volume_molecule*)reacA;
	
	/* Halfway in between where we were and where we react should be a safe away-from-wall place to remove us */
        if (loc_okay==NULL) loc_okay=&(m->pos);
	fake_hitpt.x = 0.5*hitpt->x + 0.5*loc_okay->x;
	fake_hitpt.y = 0.5*hitpt->y + 0.5*loc_okay->y;
	fake_hitpt.z = 0.5*hitpt->z + 0.5*loc_okay->z;
	
	i=count_region_from_scratch(reacA,NULL,-1,&fake_hitpt,NULL,t);
      }
      if (i) return RX_NO_MEM;
    }
  
    reacA->properties->n_deceased++;
    reacA->properties->cum_lifetime += t - reacA->birthday;
    reacA->properties->population--;
    reacA->properties = NULL;
    
    return RX_DESTROY;
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
      location of collision with wall
  Out: Value depending on outcome:
	 RX_A_OK if the molecule reflects
	 RX_FLIP if the molecule passes through
	 RX_DESTROY if the molecule stops, is destroyed, etc.
	 RX_NO_MEM on an out of memory error
       Additionally, products are created as needed.
  Note: Can assume molecule is always first in the reaction.
*************************************************************************/

int outcome_intersect(struct rxn *rx, int path, struct wall *surface,
  struct abstract_molecule *reac,short orient,double t,struct vector3 *hitpt,
  struct vector3 *loc_okay)
{
  int result,index,i;
  
  if (rx->n_pathways <= RX_SPECIAL)
  {
    rx->n_occurred++;
    if (rx->n_pathways==RX_REFLEC) return RX_A_OK;
    else return RX_FLIP; /* Flip = transparent is default special case */
  }

  index = rx->product_idx[path];

  if ((reac->properties->flags & NOT_FREE) == 0)
  {
    struct volume_molecule *m = (struct volume_molecule*) reac;
    
    result = outcome_products(surface,m,NULL,rx,path,m->subvol->local_storage,orient,0,t,hitpt,reac,NULL,reac);

    if (result==RX_NO_MEM) return RX_NO_MEM;
    else if (result == RX_BLOCKED) return RX_A_OK; /* reflect the molecule */

    rx->info[path].count++;
    rx->n_occurred++;
    
    if (rx->players[ index ] == NULL)
    {
      m->subvol->mol_count--;
      if (reac->flags&COUNT_ME)
      {
	if (hitpt==NULL) i=count_region_from_scratch(reac,NULL,-1,NULL,NULL,t);
	else
	{
	  struct vector3 fake_hitpt;
	  
	  /* Halfway in between where we were and where we react should be a safe away-from-wall place to remove us */
          if (loc_okay==NULL) loc_okay=&(m->pos);
	  fake_hitpt.x = 0.5*hitpt->x + 0.5*loc_okay->x;
	  fake_hitpt.y = 0.5*hitpt->y + 0.5*loc_okay->y;
	  fake_hitpt.z = 0.5*hitpt->z + 0.5*loc_okay->z;
	  
	  i=count_region_from_scratch(reac,NULL,-1,&fake_hitpt,NULL,t);
	}
	if (i) return RX_NO_MEM;
      }
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
  else
  {
    /* Should really be an error because we should never call outcome_intersect() on a grid molecule */
    return RX_A_OK;
  }
}


/*************************************************************************
reaction_wizardry:
  In: a list of releases to magically cause
      the wall associated with the release, if any
      the location of the release
      the time of the release
  Out: 0 if successful, 1 on failure (usually out of memory).
       Each release event in the list is triggered at a location that
       is relative to the location of the release and the surface normal
       of the wall.  The surface normal of the wall is considered to be
       the +Z direction.  Other coordinates are rotated in the "natural"
       way (i.e. the XYZ coordinate system has the Z-axis rotated directly
       to the direction of the normal and the other coordinates follow
       along naturally; if the normal is in the -Z direction, the rotation
       is about the X-axis.)
  Note: this function isn't all that cheap computationally because of
        all the math required to compute the right coordinate transform.
	If this gets really really heavily used, we should store the
	coordinate transform off of the wall data structure.
  Note: it would be more efficient to skip calculating the transform if
        the release type didn't use it (e.g. release by region).
  Note: if we wanted to be extra-super clever, we could actually schedule
        this event instead of running it and somehow have it start a
	time-shifted release pattern (so we could have delays and stuff).
*************************************************************************/

int reaction_wizardry(struct magic_list *incantation,struct wall *surface,struct vector3 *hitpt,double t)
{
  struct release_event_queue req; /* Create a release event on the fly */
  int i;
  
  /* Release event happens "now" */
  req.next=NULL;
  req.event_time=t;
  req.train_counter=0;
  req.train_high_time=t;
  
  /* Set up transform to place products at site of reaction */
  if (hitpt==NULL)
  {
    init_matrix(req.t_matrix);
  }
  else if (surface==NULL || !distinguishable(surface->normal.z,1.0,EPS_C)) /* Just need a translation */
  {
    init_matrix(req.t_matrix);
    req.t_matrix[3][0] = hitpt->x;
    req.t_matrix[3][1] = hitpt->y;
    req.t_matrix[3][2] = hitpt->z;
  }
  else /* Set up transform that will translate and then rotate Z axis to align with surface normal */
  {
    struct vector3 scale = {1.0,1.0,1.0};  /* No scaling */
    struct vector3 axis = {1.0,0.0,0.0};   /* X-axis is default */
    double cos_theta;
    double degrees;
    
    cos_theta = surface->normal.z;   /* (0,0,1) . surface->normal */
    if (!distinguishable(cos_theta,-1.0,EPS_C))
    {
      degrees=180.0;  /* Upside-down */
    }
    else
    {
      /* (0,0,1) x surface->normal */
      axis.x = -surface->normal.y;
      axis.y = surface->normal.x;
      axis.z = 0.0;
      
      degrees = acos(cos_theta)*180.0/MY_PI;
    }
    tform_matrix(&scale,hitpt,&axis,degrees,req.t_matrix);
  }
  
  /* Now we're ready to cast our spell! */
  for ( ; incantation!=NULL ; incantation=incantation->next )
  {
    if (incantation->type != magic_release) continue;  /* Only know how to magically release stuff */
    
    req.release_site = (struct release_site_obj*)incantation->data;
    
    i = release_molecules(&req);
    if (i) return i;
  }
  
  return 0;
}

