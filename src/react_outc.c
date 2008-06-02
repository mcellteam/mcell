/**************************************************************************\
** File: react_outc.c                                                     **
**                                                                        **
** Purpose: Implements specific reaction outcome pathways.                **
**                                                                        **
\**************************************************************************/

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "logging.h"
#include "rng.h"
#include "util.h"
#include "grid_util.h"
#include "mcell_structs.h"
#include "count_util.h"
#include "react.h"
#include "vol_util.h"
#include "macromolecule.h"

static int outcome_products(struct wall *w,struct volume_molecule *reac_m,
                            struct grid_molecule *reac_g,struct rxn *rx,int path,struct storage *local,
                            short orientA,short orientB,double t,struct vector3 *hitpt,
                            struct abstract_molecule *reacA,struct abstract_molecule *reacB,
                            struct abstract_molecule *moving);
static int outcome_products_trimol_reaction(struct wall *w,
                                            struct volume_molecule *reac_m, struct grid_molecule *reac_g,
                                            struct rxn *rx,int path,struct storage *local,
                                            short orientA, short orientB, short orientC,
                                            double t,struct vector3 *hitpt,
                                            struct abstract_molecule *reacA,struct abstract_molecule *reacB,
                                            struct abstract_molecule *reacC, struct abstract_molecule *moving);

extern struct volume *world;

static int reaction_wizardry(struct magic_list *incantation,
                             struct wall *surface,
                             struct vector3 *hitpt,
                             double t);

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
        Products are created as necessary and scheduled.
*************************************************************************/

static int outcome_products(struct wall *w,struct volume_molecule *reac_m,
  struct grid_molecule *reac_g,struct rxn *rx,int path,struct storage *local,
  short orientA,short orientB,double t,struct vector3 *hitpt,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  struct abstract_molecule *moving)
{
  int bounce = RX_A_OK;

  struct volume_molecule *m;
  struct grid_molecule *g;
  struct species *p;
  struct surface_grid *sg;
  u_int bits;
  int k;
  int i0 = rx->product_idx[path]; /*index of the first product for the pathway*/
  int iN = rx->product_idx[path+1];/*index of the first product for the next pathway*/
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
  int vol_rev_flag = 0;
  struct abstract_molecule *old_subunit = NULL;
  struct vector3 pos3d;
  struct subvolume *gsv = NULL;
  
  struct grid_molecule *surf_count_complex = NULL, *surf_count_subunit = NULL;
  int surf_count_idx = 0;

#define FLAG_NOT_SET 0
#define FLAG_USE_UV_LOC 1
#define FLAG_USE_REACA_UV 2
#define FLAG_USE_REACB_UV 3
#define FLAG_USE_RANDOM 4
     
  /* make sure that reacA corresponds to rx->players[0], and
     reacB - to rx->players[1] */ 
  if (reacA->properties == rx->players[1] && reacA->properties != rx->players[0])
  {
    plist[0] = reacA;
    reacA = reacB;
    reacB = plist[0];
    
    short tmp = orientA;
    orientA = orientB;
    orientB = tmp;
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

  if (rx->is_complex)
  {
    if (reacA->flags & COMPLEX_MEMBER)
      old_subunit = reacA;
    else if (reacB != NULL  &&  reacB->flags & COMPLEX_MEMBER)
      old_subunit = reacB;
    else
      mcell_internal_error("Macromolecular reaction [%s] occurred, but neither molecule is a subunit (%s and %s).",
                           rx->sym->name, reacA->properties->sym->name, reacB->properties->sym->name);
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

 
    for (int n_product=i0+rx->n_reactants; n_product<iN; n_product++)
    {
      if (rx->players[n_product]->flags&ON_GRID)
      {
        if(replace_p1 && replace_p2){
              glist[n_product - (i0+rx->n_reactants)] = reac_g->grid;
	      xlist[n_product - (i0+rx->n_reactants)] = reac_g->grid_index;
              if((struct abstract_molecule *)reac_g == reacA){
	           flist[n_product - (i0+rx->n_reactants)] = FLAG_USE_REACA_UV;
	           replace_p1=0;
              }else{
	           flist[n_product - (i0+rx->n_reactants)] = FLAG_USE_REACB_UV;
	           replace_p2=0;
              }
	      continue;
        }else if (replace_p1){
          glist[n_product - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacA)->grid;
	  xlist[n_product - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacA)->grid_index;
	  flist[n_product - (i0+rx->n_reactants)] = FLAG_USE_REACA_UV;
	  replace_p1=0;
	  continue;
	}
	else if (replace_p2)
	{
	  glist[n_product - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacB)->grid;
	  xlist[n_product - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacB)->grid_index;
	  flist[n_product - (i0+rx->n_reactants)] = FLAG_USE_REACB_UV;
	  replace_p2=0;
	  continue;
	}
	else if (w->grid==NULL)
	{
	  if (create_grid(w,reac_m->subvol))
            mcell_allocfailed("Failed to create a grid for a wall.");
	  fake_idx = n_product - (i0+rx->n_reactants);
	  glist[fake_idx] = w->grid;
	  xlist[fake_idx] = uv2grid(&uv_loc,w->grid);
	  flist[fake_idx] = FLAG_USE_UV_LOC;
	  continue;
	}
	else
	{

	  struct wall *temp_w = NULL;
	  
	  if (fake_idx > -1) glist[fake_idx]->mol[ xlist[fake_idx] ] = &fake; /* Assumed empty! */

          fake_idx = n_product - (i0+rx->n_reactants);
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
    	    temp_w = search_nbhd_for_free(w,&uv_loc,world->vacancy_search_dist2,&k,&is_compatible_surface,(void *)w->surf_class);
            
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
	glist[n_product - (i0+rx->n_reactants)]=NULL;
	xlist[n_product - (i0+rx->n_reactants)]=-1;
	flist[n_product - (i0+rx->n_reactants)]=FLAG_NOT_SET;
      }
    }
  }

  /* We know there's space, so now actually create everyone */
  if (hitpt!=NULL) memcpy(&xyz_loc,hitpt,sizeof(struct vector3));
  else if (reac_g!=NULL) uv2xyz(&(reac_g->s_pos),reac_g->grid->surface,&xyz_loc);
  else memcpy(&xyz_loc,&(reac_m->pos),sizeof(struct vector3));
  
  for (int n_product = i0+rx->n_reactants; n_product<iN; n_product++)
  {
    p = rx->players[n_product];
    
    if ( (p->flags & ON_GRID) != 0 )
    {
      if (reac_g!=NULL || (reac_m!=NULL && w!=NULL))
      {
        k = n_product-(i0+rx->n_reactants); 

        g = CHECKED_MEM_GET(local->gmol, "grid molecule");
        g->birthplace = local->gmol;
        g->birthday = t;
        g->properties = p;
        g->cmplx = NULL;
        p->population++;
        g->flags = TYPE_GRID | ACT_NEWBIE | IN_SCHEDULE;
        if (p->space_step>0) g->flags |= ACT_DIFFUSE;

        g->t = t;
        g->t2 = 0.0;
        sg = g->grid = glist[k];
        int grid_index = g->grid_index = xlist[k];

        if ((p->flags&COUNT_ENCLOSED) != 0) g->flags |= COUNT_ME;

        if (old_subunit  &&  rx->is_complex[n_product])
        {
          struct grid_molecule *old_g = (struct grid_molecule *) old_subunit;
          struct complex_species *cspec = (struct complex_species *)(old_subunit->cmplx[0]->properties);
          int num_subunits = cspec->num_subunits;
          int idx = macro_subunit_index(old_subunit);

          g->flags |= COMPLEX_MEMBER;
          g->flags &= ~ACT_DIFFUSE;

          /* Connect up new subunit to complex */
          old_g->cmplx[idx + 1] = g;
          g->cmplx = old_g->cmplx;

          /* Bind subunit to old molecule position */
          g->s_pos.u = old_g->s_pos.u;
          g->s_pos.v = old_g->s_pos.v;

          if (old_subunit->properties != g->properties  ||  old_g->orient != g->orient)
          {
            surf_count_complex = old_g->cmplx[0];
            surf_count_subunit = old_g;
            surf_count_idx = idx;

            int update_subunit[ num_subunits ];
            macro_count_inverse_related_subunits(cspec, update_subunit, idx);
            update_subunit[idx] = 0;

            int subunit_idx;
            for (subunit_idx = 1; subunit_idx <= num_subunits; ++ subunit_idx)
            {
              if (! update_subunit[subunit_idx - 1])
                continue;

              if (g->cmplx[subunit_idx]->flags & ACT_REACT)
              {
                g->cmplx[subunit_idx]->t2 = 0.0;
                g->cmplx[subunit_idx]->flags |= ACT_CHANGE;
                if (g->cmplx[subunit_idx]->flags & IN_SCHEDULE)
                {
                  uv2xyz(&g->cmplx[subunit_idx]->s_pos, g->cmplx[subunit_idx]->grid->surface, &pos3d);
                  gsv = find_subvolume(&pos3d, gsv);
                  schedule_reschedule(gsv->local_storage->timer, g->cmplx[subunit_idx], t);
                }
              }
            }
          }
          else
          {
            old_subunit->cmplx = NULL;
            old_subunit = NULL;
          }
        }
        else
        {
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
                UNHANDLED_CASE(flist[k]);
                break;
            }
          }
          else grid2uv(sg, grid_index, &(g->s_pos));
        }

        /* NOTE: This must be done after the macromolecular processing occurs --
         * otherwise, we may not be matching the right reactions... */
        if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)g)!= NULL || (p->flags&CAN_GRIDWALL)!=0) g->flags |= ACT_REACT;

        plist[n_product-i0] = (struct abstract_molecule*)g;
        ptype[n_product-i0] = 'g';
        sg->n_occupied++;
        sg->mol[grid_index] = g;

        uv2xyz(&g->s_pos, g->grid->surface, &pos3d);
        gsv = find_subvolume(&pos3d, gsv);
        if (schedule_add(gsv->local_storage->timer,g))
          mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                            g->properties->sym->name);
      }
      else /* Should never happen, but it doesn't hurt to be safe */
      {
        plist[n_product-i0] = NULL;
        ptype[n_product-i0] = 0;
        continue;
      }
    }
    else /* volume molecule */
    {
      m = CHECKED_MEM_GET(local->mol, "volume molecule");
      m->birthplace = local->mol;
      m->birthday = t;
      m->properties = p;
      m->cmplx = NULL;
      p->population++;
      m->prev_v = NULL;
      m->next_v = NULL;

      m->flags = TYPE_3D | ACT_NEWBIE | IN_VOLUME | IN_SCHEDULE;
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

      if ((p->flags&COUNT_SOME_MASK) != 0) m->flags |= COUNT_ME;
      
      if (old_subunit  &&  rx->is_complex[n_product])
      {
        struct volume_molecule *old_m = (struct volume_molecule *) old_subunit;
        struct complex_species *cspec = (struct complex_species *)(old_subunit->cmplx[0]->properties);
        int num_subunits = cspec->num_subunits;
        int idx = macro_subunit_index(old_subunit);

        /* Set appropriate flags - mol is a subunit, cannot diffuse */
        m->flags |= COMPLEX_MEMBER;
        m->flags &= ~ACT_DIFFUSE;

        /* Connect up new subunit to complex */
        old_m->cmplx[idx + 1] = m;
        m->cmplx = old_m->cmplx;

        m->pos.x = old_m->pos.x;
        m->pos.y = old_m->pos.y;
        m->pos.z = old_m->pos.z;
        m->subvol = old_m->subvol;
        ht_add_molecule_to_list(&m->subvol->mol_by_species, m);
        m->subvol->mol_count++;

        if (old_subunit->properties != m->properties)
        {
          if (count_complex(old_m->cmplx[0], old_m, idx))
            mcell_allocfailed("Failed to update counts for complex subunits after a reaction.");

          int update_subunit[ num_subunits ];
          macro_count_inverse_related_subunits(cspec, update_subunit, idx);
          update_subunit[idx] = 0;

          int subunit_idx;
          for (subunit_idx = 1; subunit_idx <= num_subunits; ++ subunit_idx)
          {
            if (! update_subunit[subunit_idx - 1])
              continue;

            if (m->cmplx[subunit_idx]->flags & ACT_REACT)
            {
              m->cmplx[subunit_idx]->t2 = 0.0;
              m->cmplx[subunit_idx]->flags |= ACT_CHANGE;
              if (m->cmplx[subunit_idx]->flags & IN_SCHEDULE)
                schedule_reschedule(m->cmplx[subunit_idx]->subvol->local_storage->timer, m->cmplx[subunit_idx], t);
            }
          }
        }

        old_subunit->cmplx = NULL;
        old_subunit = NULL;
      }
      else
      {
        m->cmplx = NULL;
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
            ht_add_molecule_to_list(&m->subvol->mol_by_species, m);
            m->subvol->mol_count++;
          }
          /* oriented case handled below after orientation is set */
        }
        else if (reac_g != NULL)
        {
          if (hitpt==NULL) uv2xyz(&(reac_g->s_pos) , reac_g->grid->surface , &(m->pos));
          m->subvol = find_subvolume(&(m->pos),reac_g->grid->subvol);
        }
      }

      /* NOTE: This must be done after the macromolecular processing occurs --
       * otherwise, we may not be matching the right reactions... */
      if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)m) != NULL) m->flags |= ACT_REACT;

      plist[n_product-i0] = (struct abstract_molecule*)m;
      ptype[n_product-i0] = 'm';
      m->t = t;
      m->t2 = 0.0;

      if (schedule_add(local->timer, m))
        mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                          m->properties->sym->name);
      
    }
  }

  /* Finally, set orientations correctly */
  bits = 0;
  for (int n_product=i0; n_product<iN; n_product++)
  {
    if (rx->players[n_product]==NULL) continue; 
   
    if ( ptype[n_product-i0] != 0 && (ptype[n_product-i0]!='m' || w!=NULL) )
    {
      if (rx->geometries[n_product] == 0)
      {
        porient[n_product-i0] = (rng_uint(world->rng) & 1) ? 1 : -1;
      }
      else
      {
        int geometry;
        if (rx->geometries[n_product] < 0)
        {
          geometry = -rx->geometries[n_product];
          k = -1;
        }
        else
        {
          geometry = rx->geometries[n_product];
          k = 1;
        }
        
        if (geometry > (int) rx->n_reactants) porient[n_product-i0] = k*porient[geometry-(rx->n_reactants+1)];
        else if (geometry==1) porient[n_product-i0] = k*orientA;
        else if (geometry==2 && reacB!=NULL) porient[n_product-i0] = k*orientB;
        else porient[n_product-i0] = k;
      }
      
      if (ptype[n_product-i0]=='g')
      {
        ((struct grid_molecule*)plist[n_product-i0])->orient = porient[n_product-i0];
      }
      else if (moving == plist[n_product-i0])
      {
        if (moving==reacA)
        {
          if (orientA==porient[n_product-i0]) bounce = RX_A_OK;
          else bounce = RX_FLIP;
        }
        else
        {
          if (orientB==porient[n_product-i0]) bounce = RX_A_OK;
          else bounce = RX_FLIP;
        }
      }
      else if (ptype[n_product-i0]=='m')
      {
        double bump;
        m = (struct volume_molecule*)plist[n_product-i0];
        if (porient[n_product-i0]>0) bump = EPS_C;
        else bump = -EPS_C;
	
	if ((m->flags&ACT_CLAMPED) && world->surface_reversibility)
	{
	  m->index = (porient[n_product-i0]>0)?1:-1; /* Which direction do we move? */
	}
        	
        /* Note: no raytracing here so it is rarely possible to jump through closely spaced surfaces */
        m->pos.x += bump*w->normal.x;
        m->pos.y += bump*w->normal.y;
        m->pos.z += bump*w->normal.z;

        m->subvol = find_subvolume(&(m->pos),m->subvol);
        ht_add_molecule_to_list(&m->subvol->mol_by_species, m);
        m->subvol->mol_count++;
      }
    }
    else if (world->volume_reversibility && reac_g==NULL && w==NULL && ptype[n_product-i0]=='m') /* Not orientable */
    {
      m = (struct volume_molecule*)plist[n_product-i0];
      m->index = world->dissociation_index;
      if (n_product-i0 >= (int) rx->n_reactants) m->flags |= ACT_CLAMPED;
      vol_rev_flag=1;
    }
    
    if (n_product >= i0 + (int) rx->n_reactants &&
        (plist[n_product-i0]->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
    {
      if (count_region_from_scratch(plist[n_product-i0],NULL,1,NULL,w,t))
        mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                          plist[n_product-i0]->properties->sym->name);
    }
  }

  if (surf_count_complex) {
    if (count_complex_surface(surf_count_complex, surf_count_subunit, surf_count_idx))
      mcell_allocfailed("Failed to update region counts for surface macromolecule subunit '%s/%s' after a reaction.",
                        surf_count_complex->properties->sym->name,
                        surf_count_subunit->properties->sym->name);
    old_subunit->cmplx = NULL;
    old_subunit = NULL;
  }
  
  if (vol_rev_flag)
  {
    world->dissociation_index--;
    if (world->dissociation_index < DISSOCIATION_MIN) world->dissociation_index=DISSOCIATION_MAX;
  }
  
  /* Handle events triggered off of named reactions */
  if (rx->info[path].pathname!=NULL)
  {
    /* No flags for reactions so we have to check regions if we have waypoints! Fix to be more efficient for WORLD-only counts? */
    if (world->place_waypoints_flag)
    {
      if (count_region_from_scratch(NULL,rx->info[path].pathname,1,&xyz_loc,w,t))
        mcell_allocfailed("Failed to update region counts for '%s' reactions.",
                          rx->info[path].pathname->sym->name);
    }
    
    /* Other magical stuff.  For now, can only trigger releases. */
    if (rx->info[path].pathname->magic!=NULL)
    {
      if (reaction_wizardry(rx->info[path].pathname->magic,w,&xyz_loc,t))
        mcell_allocfailed("Failed to complete reaction triggered release after a '%s' reaction.",
                          rx->info[path].pathname->sym->name);
    }
  }
  

#undef FLAG_NOT_SET
#undef FLAG_USE_UV_LOC
#undef FLAG_USE_REACA_UV
#undef FLAG_USE_REACB_UV
#undef FLAG_USE_RANDOM
  
  return bounce;
  
}

/*************************************************************************
outcome_products_trimol_reaction:
   In: relevant wall in the interaction, if any
       first free molecule in the interaction, if any
       first surface molecule in the interaction, if any
       reaction that is occuring
       path that the reaction is taking
       local storage for creating new molecules
       orientations of molecules in the reaction
       time that the reaction is occurring
       location of the reaction (may be NULL)
       the reactants (the last one is also the furthest
            one from the moving molecule)
   Out: Value depending on how orientations changed--
          RX_FLIP moving molecule passed through membrane
          RX_A_OK everything went fine, nothing extra to do
        Products are created as necessary and scheduled.

   Note: This function does not include macromolecules support yet, as
         macromolecules+trimol is not yet supported.
*************************************************************************/
static int outcome_products_trimol_reaction(struct wall *w,
  struct volume_molecule *reac_m, struct grid_molecule *reac_g,
  struct rxn *rx,int path, struct storage *local,
  short orientA, short orientB, short orientC,
  double t, struct vector3 *hitpt,
  struct abstract_molecule *reacA, struct abstract_molecule *reacB,
  struct abstract_molecule *reacC, struct abstract_molecule *moving)
{
  int bounce = RX_A_OK;
  struct volume_molecule *m;
  struct grid_molecule *g;
  struct species *p;
  struct surface_grid *sg;
  struct subvolume *gsv = NULL;
  struct vector3 pos3d;
  u_int bits;
  int k;
  int i0 = rx->product_idx[path]; /*index of the first product for the pathway*/
  int iN = rx->product_idx[path+1];/*index of the first product for the next pathway*/
  int replace_p1 = 0; /* flag for the product to replace position of reactant1 */
  int replace_p2 = 0; /* flag for the product to replace position of reactant2 */
  int replace_p3 = 0; /* flag for the product to replace position of reactant3 */
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
  int vol_rev_flag = 0;
  int trimol_reaction_flag = 0; /* checks whether the reaction is a
                                     trimol reaction */ 


#define FLAG_NOT_SET 0
#define FLAG_USE_UV_LOC 1
#define FLAG_USE_REACA_UV 2     /* use reacA position */
#define FLAG_USE_REACB_UV 3     /* use reacB position */
#define FLAG_USE_REACC_UV 4     /* use reacC position */
#define FLAG_USE_RANDOM 5
     
  if (rx->is_complex)
    mcell_internal_error("Macromolecular reaction [%s] occurred inside trimolecular reaction code.", rx->sym->name);

  /* make sure that reacA corresponds to rx->players[0], 
     reacB - to rx->players[1], and reacC - to rx->players[2] */
  if(reacA != NULL && reacB != NULL && reacC != NULL)
  {
   /* trimolecular reaction */
     if(reacA->properties == rx->players[0])
     {
       if(reacB->properties == rx->players[2] && reacB->properties != rx->players[1]){
          plist[0] = reacB;
          reacB = reacC;
          reacC = plist[0];
    
          short tmp = orientB;
          orientB = orientC;
          orientC = tmp;
       }
     }else if(reacA->properties == rx->players[1]){
    
       if(reacB->properties == rx->players[0] && reacB->properties != rx->players[1]){
          /* switch reacA and reacB */
          plist[0] = reacB;
          reacB = reacA;
          reacA = plist[0];
    
          short tmp = orientB;
          orientB = orientA;
          orientA = tmp;
       }else if(reacC->properties == rx->players[0]){
          /* switch reacA and reacC */
          plist[0] = reacA;
          reacA = reacC;
          reacC = plist[0];
    
          short tmp = orientA;
          orientA = orientC;
          orientC = tmp;

          /* now switch reacC and reacB  */ 
          plist[0] = reacB;
          reacB = reacC;
          reacC = plist[0];
    
          tmp = orientB;
          orientB = orientC;
          orientC = tmp;
       }
     }else if(reacA->properties == rx->players[2]){
        if(reacB->properties == rx->players[0])
        {
          /* switch reacA and reacB */
          plist[0] = reacB;
          reacB = reacA;
          reacA = plist[0];
    
          short tmp = orientB;
          orientB = orientA;
          orientA = tmp;
        
          /* switch reacB and reacC */
          plist[0] = reacB;
          reacB = reacC;
          reacC = plist[0];
    
          tmp = orientB;
          orientB = orientC;
          orientC = tmp;
    
        }else if ((reacC->properties == rx->players[0]) &&
           (reacC->properties != rx->players[2])){
          /* switch reacA and reacC */
          plist[0] = reacA;
          reacA = reacC;
          reacC = plist[0];
    
          short tmp = orientA;
          orientA = orientC;
          orientC = tmp;
       }
     }
  }else{
    /* bimolecular reaction */
    /* make sure that reacA corresponds to rx->players[0], and
     reacB - to rx->players[1] */
     if (reacA->properties == rx->players[1] && reacA->properties != rx->players[0])
     {
       plist[0] = reacA;
       reacA = reacB;
       reacB = plist[0];

       short tmp = orientA;
       orientA = orientB;
       orientB = tmp;
    }
  }

         
  if((reacA != NULL) && (reacB != NULL) && (reacC != NULL)){
     trimol_reaction_flag = 1;
  }

  

  plist[0] = reacA;
  
  if ( (reacA->properties->flags&ON_GRID)!=0 ) ptype[0] = 'g';
  else if ( (reacA->properties->flags&NOT_FREE)==0 ) ptype[0] = 'm';
  else ptype[0] = '!';

  if(rx->n_reactants > 1)
  {
     if(reacB == NULL)
     {
        ptype[1] = 'w';
        plist[1] = NULL;
     }else{
       plist[1] = reacB;
       if ( (reacB->properties->flags&ON_GRID)!=0 ) ptype[1] = 'g';
       else if ( (reacB->properties->flags&NOT_FREE)==0 ) ptype[1] = 'm';
       else ptype[1] = '!';
     }

     if(rx->n_reactants > 2)
     {
        if(reacC == NULL)
        {
           ptype[2] = 'w';
           plist[2] = NULL;
        }else{
           plist[2] = reacC;
           if ( (reacC->properties->flags&ON_GRID)!=0 ) ptype[2] = 'g';
           else if ( (reacC->properties->flags&NOT_FREE)==0 ) ptype[2] = 'm';
           else ptype[2] = '!';
        }
     }
  }
  
  /* Make sure there's space for the reaction to occur */
  /* FIXME--could speed this up with some pre-computation of reactions to at least see if we need to bother */
  k = -1;
  
  if (ptype[0]=='g' && rx->players[i0]==NULL) replace_p1=1;
  if (rx->n_reactants > 1 && ptype[1]=='g' && rx->players[i0+1]==NULL) replace_p2=1;
  if (rx->n_reactants > 2 && ptype[2]=='g' && rx->players[i0+2]==NULL) replace_p3=1;
   

  if (reac_g!=NULL  || (reac_m != NULL && w!=NULL))  /* Surface involved */
  {
    if(reac_g != NULL) memcpy(&uv_loc , &(reac_g->s_pos) , sizeof(struct vector2));
    else xyz2uv(hitpt,w,&uv_loc);
 
    for (int n_player=i0+rx->n_reactants; n_player<iN; n_player++)
    {
      if (rx->players[n_player]->flags&ON_GRID)
      {
  
        if(replace_p1 && replace_p2 && replace_p3){
          glist[n_player - (i0+rx->n_reactants)] = reac_g->grid;
	  xlist[n_player - (i0+rx->n_reactants)] = reac_g->grid_index;
          if((struct abstract_molecule *)reac_g == reacA){
              flist[n_player - (i0+rx->n_reactants)] = FLAG_USE_REACA_UV;
              replace_p1 = 0;
          }else if((struct abstract_molecule *)reac_g == reacB){
              flist[n_player - (i0+rx->n_reactants)] = FLAG_USE_REACB_UV;
              replace_p2 = 0;
          }else if((struct abstract_molecule *)reac_g == reacC){
              flist[n_player - (i0+rx->n_reactants)] = FLAG_USE_REACC_UV;
              replace_p3 = 0;
          }
          continue;
        }else if (replace_p1){
          glist[n_player - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacA)->grid;
	  xlist[n_player - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacA)->grid_index;
	  flist[n_player - (i0+rx->n_reactants)] = FLAG_USE_REACA_UV;
	  replace_p1=0;
	  continue;
	}
	else if (replace_p2)
	{
	  glist[n_player - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacB)->grid;
	  xlist[n_player - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacB)->grid_index;
	  flist[n_player - (i0+rx->n_reactants)] = FLAG_USE_REACB_UV;
	  replace_p2=0;
	  continue;
	}
	else if (replace_p3)
	{
	  glist[n_player - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacC)->grid;
	  xlist[n_player - (i0+rx->n_reactants)] = ((struct grid_molecule*)reacC)->grid_index;
	  flist[n_player - (i0+rx->n_reactants)] = FLAG_USE_REACC_UV;
	  replace_p3=0;
	  continue;
	}
        else if (w->grid==NULL)
	{
	  if (create_grid(w,reac_m->subvol))
            mcell_allocfailed("Failed to create a grid for a wall.");
	  fake_idx = n_player - (i0+rx->n_reactants);
	  glist[fake_idx] = w->grid;
	  xlist[fake_idx] = uv2grid(&uv_loc,w->grid);
	  flist[fake_idx] = FLAG_USE_UV_LOC;
	  continue;
	}
	else
	{

	  struct wall *temp_w = NULL;
	  
	  if (fake_idx > -1) glist[fake_idx]->mol[ xlist[fake_idx] ] = &fake; /* Assumed empty! */

          fake_idx = n_player - (i0+rx->n_reactants);
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
    	    temp_w = search_nbhd_for_free(w,&uv_loc,world->vacancy_search_dist2,&k,&is_compatible_surface,(void *)w->surf_class);
            
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
	glist[n_player - (i0+rx->n_reactants)]=NULL;
	xlist[n_player - (i0+rx->n_reactants)]=-1;
	flist[n_player - (i0+rx->n_reactants)]=FLAG_NOT_SET;
      }
    }
  }

  /* We know there's space, so now actually create everyone */
  if (hitpt!=NULL) memcpy(&xyz_loc,hitpt,sizeof(struct vector3));
  else if (reac_g!=NULL) uv2xyz(&(reac_g->s_pos),reac_g->grid->surface,&xyz_loc);
  else memcpy(&xyz_loc,&(reac_m->pos),sizeof(struct vector3));
  
  for (int n_player=i0+rx->n_reactants; n_player<iN; n_player++)
  {
    p = rx->players[n_player];
    
    if ( (p->flags & ON_GRID) != 0 )
    {
      if (reac_g != NULL || (reac_m != NULL && w != NULL))
      {
        k = n_player-(i0+rx->n_reactants); 
	
        g = CHECKED_MEM_GET(local->gmol, "grid molecule");
	g->birthplace = local->gmol;
	g->birthday = t;
	g->properties = p;
        g->cmplx = NULL;
	p->population++;
	g->flags = TYPE_GRID | ACT_NEWBIE | IN_SCHEDULE;
	if (p->space_step>0) g->flags |= ACT_DIFFUSE;
	if (trigger_unimolecular(p->hashval,(struct abstract_molecule*)g)!= NULL || (p->flags&CAN_GRIDWALL)!=0) g->flags |= ACT_REACT;
	
	g->t = t;
	g->t2 = 0.0;
	sg = g->grid = glist[k];
	int grid_index = g->grid_index = xlist[k];

        if ((p->flags&COUNT_ENCLOSED) != 0) g->flags |= COUNT_ME;
	
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

            case FLAG_USE_REACC_UV:
              memcpy(&(g->s_pos),&(((struct grid_molecule*)reacC)->s_pos),sizeof(struct vector2));
              break;

            case FLAG_USE_UV_LOC:
              memcpy(&(g->s_pos),&(uv_loc),sizeof(struct vector2));
              break;

            case FLAG_USE_RANDOM:
              grid2uv_random(glist[k],xlist[k],&(g->s_pos)); 
              break;

            default:
              UNHANDLED_CASE(flist[k]);
              break;
          }
	}
        else grid2uv(sg, grid_index, &(g->s_pos));

        sg->n_occupied++;
                 
        sg->mol[grid_index] = g;
	
        plist[n_player-i0] = (struct abstract_molecule*)g;
        ptype[n_player-i0] = 'g';
	
        uv2xyz(&g->s_pos, g->grid->surface, &pos3d);
        gsv = find_subvolume(&pos3d, gsv);
        if (schedule_add(gsv->local_storage->timer, g))
          mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                            g->properties->sym->name);
      }
      else /* Should never happen, but it doesn't hurt to be safe */
      {
        plist[n_player-i0] = NULL;
        ptype[n_player-i0] = 0;
        continue;
      }
    }
    else /* volume molecule */
    {
      m = CHECKED_MEM_GET(local->mol, "volume molecule");
      m->birthplace = local->mol;
      m->birthday = t;
      m->properties = p;
      m->cmplx = NULL;
      p->population++;
      m->prev_v = NULL;
      m->next_v = NULL;

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

      if ((p->flags&COUNT_SOME_MASK) != 0) m->flags |= COUNT_ME;
     
      if(hitpt != NULL)
      { 
         m->pos.x = hitpt->x;
         m->pos.y = hitpt->y;
         m->pos.z = hitpt->z;
         if(trimol_reaction_flag){
            m->subvol = find_coarse_subvol(hitpt);
         }
      }      
 
      if (reac_m != NULL)
      {
       
        if (hitpt==NULL || ((struct abstract_molecule*)reac_m != moving && !trimol_reaction_flag))
        {
          m->pos.x = reac_m->pos.x;
          m->pos.y = reac_m->pos.y;
          m->pos.z = reac_m->pos.z;
        }
        if(!trimol_reaction_flag){
           m->subvol = reac_m->subvol;
        }
 
        if (w==NULL) /* place product of non-orientable reaction in volume */
        {
          ht_add_molecule_to_list(&m->subvol->mol_by_species, m);
          m->subvol->mol_count++;
        }
        /* oriented case handled below after orientation is set */
      }
      else if (reac_g != NULL)
      {
        if (hitpt==NULL) uv2xyz(&(reac_g->s_pos) , reac_g->grid->surface , &(m->pos));
        if(!trimol_reaction_flag){
           m->subvol = find_subvolume(&(m->pos),reac_g->grid->subvol);
        }
        
      }
      plist[n_player-i0] = (struct abstract_molecule*)m;
      ptype[n_player-i0] = 'm';
      m->t = t;
      m->t2 = 0.0;

      if (schedule_add( local->timer, m))
        mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                          m->properties->sym->name);
      
    }
  }

  /* Finally, set orientations correctly */
  bits = 0;
  for (int n_player=i0; n_player<iN; n_player++)
  {
    if (rx->players[n_player]==NULL) continue; 
   
    if ( ptype[n_player-i0] != 0 && (ptype[n_player-i0]!='m' || w!=NULL) )
    {
      if (rx->geometries[n_player] == 0)
      {
        porient[n_player-i0] = (rng_uint(world->rng) & 1) ? 1 : -1;
      }
      else
      {
        int geometry;
        if (rx->geometries[n_player] < 0)
        {
          geometry = -rx->geometries[n_player];
          k = -1;
        }
        else
        {
          geometry = rx->geometries[n_player];
          k = 1;
        }
        
        if (geometry > (int) rx->n_reactants) porient[n_player-i0] = k*porient[geometry-(rx->n_reactants+1)];
        else if (geometry==1) porient[n_player-i0] = k*orientA;
        else if (geometry==2 && reacB!=NULL) porient[n_player-i0] = k*orientB;
        else if (geometry==3 && reacC!=NULL) porient[n_player-i0] = k*orientC;
        else porient[n_player-i0] = k;
        
      }
      
      if (ptype[n_player-i0]=='g')
      {
        ((struct grid_molecule*)plist[n_player-i0])->orient = porient[n_player-i0];
      }
      else if (moving == plist[n_player-i0])
      {
        if (moving==reacA)
        {
          if (orientA==porient[n_player-i0]) bounce = RX_A_OK;
          else bounce = RX_FLIP;
        }
        else if(moving == reacB)
        {
          if (orientB==porient[n_player-i0]) bounce = RX_A_OK;
          else bounce = RX_FLIP;
        }
        else
        {
          if (orientC==porient[n_player-i0]) bounce = RX_A_OK;
          else bounce = RX_FLIP;
        }
      }
      else if (ptype[n_player-i0]=='m')
      {
        double bump;
        m = (struct volume_molecule*)plist[n_player-i0];
        if (porient[n_player-i0]>0) bump = EPS_C;
        else bump = -EPS_C;
	
	if ((m->flags&ACT_CLAMPED) && world->surface_reversibility)
	{
          m->index = (porient[n_player-i0]>0)?1:-1; /* Which direction do we move? */
	}
        	
        /* Note: no raytracing here so it is rarely possible to jump through closely spaced surfaces */
        m->pos.x += bump*w->normal.x;
        m->pos.y += bump*w->normal.y;
        m->pos.z += bump*w->normal.z;

        m->subvol = find_subvolume(&(m->pos),m->subvol);
        ht_add_molecule_to_list(&m->subvol->mol_by_species, m);
        m->subvol->mol_count++;
      }
    }
    else if (world->volume_reversibility && reac_g==NULL && w==NULL && ptype[n_player-i0]=='m') /* Not orientable */
    {
      m = (struct volume_molecule*)plist[n_player-i0];
      m->index = world->dissociation_index;
      if (n_player-i0 >= (int) rx->n_reactants) m->flags |= ACT_CLAMPED;
      vol_rev_flag=1;
    }
    
    if (n_player >= i0 + (int) rx->n_reactants &&
        (plist[n_player-i0]->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
    {
      if (count_region_from_scratch(plist[n_player-i0],NULL,1,NULL,w,t))
        mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                          plist[n_player-i0]->properties->sym->name);
    }
  }

  if (vol_rev_flag)
  {
    world->dissociation_index--;
    if (world->dissociation_index < DISSOCIATION_MIN) world->dissociation_index=DISSOCIATION_MAX;
  }

  /* Handle events triggered off of named reactions */
  if (rx->info[path].pathname!=NULL)
  {
    /* No flags for reactions so we have to check regions if we have waypoints! Fix to be more efficient for WORLD-only counts? */
    if (world->place_waypoints_flag)
    {
      if (count_region_from_scratch(NULL, rx->info[path].pathname, 1, &xyz_loc, w, t))
        mcell_allocfailed("Failed to update region counts for '%s' reactions.",
                          rx->info[path].pathname->sym->name);
    }
    
    /* Other magical stuff.  For now, can only trigger releases. */
    if (rx->info[path].pathname->magic!=NULL)
    {
      if (reaction_wizardry(rx->info[path].pathname->magic, w, &xyz_loc, t))
        mcell_allocfailed("Failed to complete reaction triggered release after a '%s' reaction.",
                          rx->info[path].pathname->sym->name);
    }
  }

#undef FLAG_NOT_SET
#undef FLAG_USE_UV_LOC
#undef FLAG_USE_REACA_UV
#undef FLAG_USE_REACB_UV
#undef FLAG_USE_REACC_UV
#undef FLAG_USE_RANDOM
  
  return bounce;
  
}

/*************************************************************************
outcome_unimolecular:
  In: the reaction that is occuring
      the path that the reaction is taking
      the molecule that is taking that path
      time that the reaction is occurring
  Out: Value based on outcome:
	 RX_BLOCKED if there was no room to put products on grid
	 RX_DESTROY if molecule no longer exists.
	 RX_A_OK if it does.
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
  
  if (result==RX_BLOCKED) return RX_BLOCKED;
  
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
        if (count_region_from_scratch((struct abstract_molecule*)m, NULL, -1, &(m->pos), NULL, m->t))
          mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                            m->properties->sym->name);
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
        if (count_region_from_scratch((struct abstract_molecule*)g, NULL, -1, NULL, NULL, g->t))
          mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                            g->properties->sym->name);
      }
    }

    reac->properties->n_deceased++;
    reac->properties->cum_lifetime += t - reac->birthday;
    reac->properties->population--;
    if (m != NULL) collect_molecule(m);
    else
    {
      reac->properties = NULL;
      mem_put(reac->birthplace, reac);
    }
    return RX_DESTROY;
  }
  else if (who_am_i != who_was_i)
  {
    if (m != NULL) collect_molecule(m);
    else
      reac->properties = NULL;
    return RX_DESTROY;
  }
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
       Products are created as needed.
  Note: reacA is the triggering molecule (e.g. moving)
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
          
   
  if (result==RX_BLOCKED) return RX_BLOCKED;
  
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
    m = NULL;
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
      if (count_region_from_scratch(reacB, NULL, -1, NULL, NULL, t))
        mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                          reacB->properties->sym->name);
    }
    
    reacB->properties->n_deceased++;
    reacB->properties->cum_lifetime += t - reacB->birthday;
    reacB->properties->population--;
    if (m != NULL) collect_molecule(m);
    else reacB->properties = NULL;
  }

  if (killA)
  {
    m = NULL;
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
        if (count_region_from_scratch(reacA, NULL, -1, NULL, NULL, t))
          mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                            reacA->properties->sym->name);
      }
    }
    else if (reacA->flags&COUNT_ME)
    {
      /* Subtlety: we made it up to hitpt, but our position is wherever we were before that! */
      if (hitpt==NULL || reacB_was_free || (reacB->properties!=NULL && (reacB->properties->flags&NOT_FREE)==0))
      {
	/* Vol-vol rx should be counted at hitpt */
        if (count_region_from_scratch(reacA, NULL, -1, hitpt, NULL, t))
          mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                            reacA->properties->sym->name);
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
	
        if (count_region_from_scratch(reacA, NULL, -1, &fake_hitpt, NULL, t))
          mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                            reacA->properties->sym->name);
      }
    }
  
    reacA->properties->n_deceased++;
    reacA->properties->cum_lifetime += t - reacA->birthday;
    reacA->properties->population--;
    if (m != NULL) collect_molecule(m);
    else reacA->properties = NULL;
    
    return RX_DESTROY;
  }

  return result;
}

/*************************************************************************
outcome_trimolecular:
  In: reaction that's occurring
      path the reaction's taking
      three molecules that are reacting (first is moving one
          and the last one is the furthest from the moving molecule or
          the one that is hit the latest)
      orientations of the molecules
      time that the reaction is occurring
      location of collision between moving molecule and the furthest target
  Out: Value based on outcome:
	 RX_FLIP if the molecule goes across the membrane
	 RX_DESTROY if the molecule no longer exists
	 RX_A_OK if everything proceeded smoothly
       Products are created as needed.
  Note: reacA is the triggering molecule (e.g. moving)
        reacC is the target furthest from the reacA
*************************************************************************/
int outcome_trimolecular(struct rxn *rx,int path,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  struct abstract_molecule *reacC, short orientA, short orientB, short orientC, 
  double t, struct vector3 *hitpt, struct vector3 *loc_okay)
{
  struct wall *w = NULL;
  struct volume_molecule *m = NULL;
  struct grid_molecule *g = NULL;
  struct storage *x;
  int result;
  /* flags */
  int killA = 0,killB = 0, killC = 0;
  int reacB_is_free = 0;
  int reacC_is_free = 0;

   if ((reacA->properties->flags & NOT_FREE) == 0)
   {
       m = (struct volume_molecule*) reacA;
   } 
   if ((reacB->properties->flags & NOT_FREE) == 0) reacB_is_free = 1;
   if ((reacC->properties->flags & NOT_FREE) == 0) reacC_is_free = 1;

    /* we will use storage of the SV where the furthest target is located
       and products be placed  */
    if((reacC->properties->flags & ON_GRID) != 0){
       g = (struct grid_molecule *)reacC;
       x = g->grid->surface->birthplace;
       w = g->grid->surface;
    }else{
       x = ((struct volume_molecule *)reacC)->subvol->local_storage; 
       if((reacB->properties->flags & ON_GRID) != 0){
          g = (struct grid_molecule *)reacB;
          w = g->grid->surface;
       }
       
    }

     result = outcome_products_trimol_reaction(w,m,g,rx,path,x,orientA, orientB, orientC,t,hitpt,reacA,reacB,reacC, reacA);  
          
     
  if (result==RX_BLOCKED) return RX_BLOCKED;
             

  rx->n_occurred++;
  rx->info[path].count++;
  
  /* Figure out if either of the reactants was destroyed */

  if (rx->players[0]==reacA->properties)
  {
    if(rx->players[1] == reacB->properties)
    {
      killC = (rx->players[ rx->product_idx[path]+2 ] == NULL);
      killB = (rx->players[ rx->product_idx[path]+1 ] == NULL);
      killA = (rx->players[ rx->product_idx[path] ] == NULL);
    }else{
      killB = (rx->players[ rx->product_idx[path]+2 ] == NULL);
      killC = (rx->players[ rx->product_idx[path]+1 ] == NULL);
      killA = (rx->players[ rx->product_idx[path] ] == NULL);
    }
  }
  else if (rx->players[0]==reacB->properties)
  {
    if(rx->players[1] == reacA->properties)
    {
      killC = (rx->players[ rx->product_idx[path]+2 ] == NULL);
      killA = (rx->players[ rx->product_idx[path]+1 ] == NULL);
      killB = (rx->players[ rx->product_idx[path] ] == NULL);
    }else{
      killA = (rx->players[ rx->product_idx[path]+2 ] == NULL);
      killC = (rx->players[ rx->product_idx[path]+1 ] == NULL);
      killB = (rx->players[ rx->product_idx[path] ] == NULL);
    }
  }else if (rx->players[0]==reacC->properties)
  {
    if(rx->players[1] == reacA->properties)
    {
      killB = (rx->players[ rx->product_idx[path]+2 ] == NULL);
      killA = (rx->players[ rx->product_idx[path]+1 ] == NULL);
      killC = (rx->players[ rx->product_idx[path] ] == NULL);
    }else{
      killA = (rx->players[ rx->product_idx[path]+2 ] == NULL);
      killB = (rx->players[ rx->product_idx[path]+1 ] == NULL);
      killC = (rx->players[ rx->product_idx[path] ] == NULL);
    }
  }


  if (killC)
  {
    m = NULL;
    if((reacC->properties->flags & ON_GRID) != 0){
       g = (struct grid_molecule *)reacC;
       if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
       g->grid->n_occupied--;
       if (g->flags&IN_SURFACE) g->flags -= IN_SURFACE;

       if (g->flags & IN_SCHEDULE)
       {
          g->grid->subvol->local_storage->timer->defunct_count++;
       }
    }else{
       m = (struct volume_molecule*)reacC;
       m->subvol->mol_count--;
       if (m->flags & IN_SCHEDULE)
       {
          m->subvol->local_storage->timer->defunct_count++;
       }
    }

    if ((reacC->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
    {
      if (count_region_from_scratch(reacC, NULL, -1, NULL, NULL, t))
        mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                          reacC->properties->sym->name);
    }
    
    reacC->properties->n_deceased++;
    reacC->properties->cum_lifetime += t - reacC->birthday;
    reacC->properties->population--;
    if (m != NULL) collect_molecule(m);
    else
    {
      reacC->properties = NULL;
      if ((reacC->flags&IN_MASK)==0) mem_put(reacC->birthplace,reacC);
    }
  }
 
  if (killB)
  {
    m = NULL;
    if((reacB->properties->flags & ON_GRID) != 0){
       g = (struct grid_molecule *)reacB;
       if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
       g->grid->n_occupied--;
       if (g->flags&IN_SURFACE) g->flags -= IN_SURFACE;

       if (g->flags & IN_SCHEDULE)
       {
          g->grid->subvol->local_storage->timer->defunct_count++;
       }
    }else{
       m = (struct volume_molecule*)reacB;
       m->subvol->mol_count--;
       if (m->flags & IN_SCHEDULE)
       {
          m->subvol->local_storage->timer->defunct_count++;
       }
    }

    if ((reacB->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
    {
      if (count_region_from_scratch(reacB, NULL, -1, NULL, NULL, t))
        mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                          reacB->properties->sym->name);
    }
    
    reacB->properties->n_deceased++;
    reacB->properties->cum_lifetime += t - reacB->birthday;
    reacB->properties->population--;
    if (m != NULL) collect_molecule(m);
    else
    {
      reacB->properties = NULL;
      if ((reacB->flags&IN_MASK)==0) mem_put(reacB->birthplace,reacB);
    }
  }

  if (killA)
  {
    m = NULL;
    if((reacA->properties->flags & ON_GRID) != 0){
       g = (struct grid_molecule *)reacA;
       if (g->grid->mol[g->grid_index]==g) g->grid->mol[g->grid_index] = NULL;
       g->grid->n_occupied--;
       if (g->flags&IN_SURFACE) g->flags -= IN_SURFACE;

       if (g->flags & IN_SCHEDULE)
       {
          g->grid->subvol->local_storage->timer->defunct_count++;
       }
    }else{
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
        if (count_region_from_scratch(reacA, NULL, -1, NULL, NULL, t))
          mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                            reacA->properties->sym->name);
      }
    }
    else if ((reacA->flags&COUNT_ME) && world->place_waypoints_flag)
    {
      /* Subtlety: we made it up to hitpt, but our position is wherever we were before that! */
      if (hitpt==NULL || (reacB_is_free && reacC_is_free))
	   /* Vol-vol-vol rx should be counted at hitpt */
      {
        if (count_region_from_scratch(reacA, NULL, -1, hitpt, NULL, t))
          mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                            reacA->properties->sym->name);
      }
      else /* reaction involving surface or grid_molecule but we don't want to count exactly on a wall or we might count on the wrong side */
      {
        struct vector3 fake_hitpt;

        m = (struct volume_molecule*)reacA;

        /* Halfway in between where we were and where we react should be a safe away-from-wall place to remove us */
        if (loc_okay == NULL)
          loc_okay=&(m->pos);
        fake_hitpt.x = 0.5*hitpt->x + 0.5*loc_okay->x;
        fake_hitpt.y = 0.5*hitpt->y + 0.5*loc_okay->y;
        fake_hitpt.z = 0.5*hitpt->z + 0.5*loc_okay->z;

        if (count_region_from_scratch(reacA, NULL, -1, &fake_hitpt, NULL, t))
          mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                            reacA->properties->sym->name);
      }
    }
     reacA->properties->n_deceased++;
     reacA->properties->cum_lifetime += t - reacA->birthday;
     reacA->properties->population--;
    if (m != NULL) collect_molecule(m);
    else reacA->properties = NULL; 

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
       Additionally, products are created as needed.
  Note: Can assume molecule is always first in the reaction.
*************************************************************************/

int outcome_intersect(struct rxn *rx, int path, struct wall *surface,
  struct abstract_molecule *reac,short orient,double t,struct vector3 *hitpt,
  struct vector3 *loc_okay)
{
  int result, idx;
  
  if (rx->n_pathways <= RX_SPECIAL)
  {
    rx->n_occurred++;
    if (rx->n_pathways==RX_REFLEC) return RX_A_OK;
    else return RX_FLIP; /* Flip = transparent is default special case */
  }

  idx = rx->product_idx[path];

  if ((reac->properties->flags & NOT_FREE) == 0)
  {
    struct volume_molecule *m = (struct volume_molecule*) reac;
    
    result = outcome_products(surface,m,NULL,rx,path,m->subvol->local_storage,orient,0,t,hitpt,reac,NULL,reac);

    if (result == RX_BLOCKED) return RX_A_OK; /* reflect the molecule */

    rx->info[path].count++;
    rx->n_occurred++;
    
    if (rx->players[idx] == NULL)
    {
      m->subvol->mol_count--;
      if (reac->flags&COUNT_ME)
      {
        if (hitpt==NULL)
        {
          if (count_region_from_scratch(reac, NULL, -1, NULL, NULL, t))
            mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                              reac->properties->sym->name);
        }
	else
	{
	  struct vector3 fake_hitpt;
	  
	  /* Halfway in between where we were and where we react should be a safe away-from-wall place to remove us */
          if (loc_okay==NULL) loc_okay=&(m->pos);
	  fake_hitpt.x = 0.5*hitpt->x + 0.5*loc_okay->x;
	  fake_hitpt.y = 0.5*hitpt->y + 0.5*loc_okay->y;
	  fake_hitpt.z = 0.5*hitpt->z + 0.5*loc_okay->z;
	  
          if (count_region_from_scratch(reac, NULL, -1, &fake_hitpt, NULL, t))
            mcell_allocfailed("Failed to update region counts for '%s' molecules after a reaction.",
                              reac->properties->sym->name);
	}
      }
      reac->properties->n_deceased++;
      reac->properties->cum_lifetime += t - reac->birthday;
      reac->properties->population--;
      if (m->flags & IN_SCHEDULE)
      {
        m->subvol->local_storage->timer->defunct_count++;
      }
      collect_molecule(m);
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
static int reaction_wizardry(struct magic_list *incantation,struct wall *surface,struct vector3 *hitpt,double t)
{
  struct release_event_queue req; /* Create a release event on the fly */
  
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
    
    if (release_molecules(&req))
      return 1;
  }
  
  return 0;
}

