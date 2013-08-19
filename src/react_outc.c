/**************************************************************************\
** File: react_outc.c                                                     **
**                                                                        **
** Purpose: Implements specific reaction outcome pathways.                **
**                                                                        **
\**************************************************************************/

#include "config.h"

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "logging.h"
#include "rng.h"
#include "util.h"
#include "grid_util.h"
#include "mcell_structs.h"
#include "count_util.h"
#include "react.h"
#include "vol_util.h"
#include "macromolecule.h"
#include "wall_util.h"

static int outcome_products(struct wall *w,
                            struct vector3 *hitpt,
                            double t,
                            struct rxn *rx,
                            int path,
                            struct abstract_molecule *reacA,
                            struct abstract_molecule *reacB,
                            short orientA,
                            short orientB);
static int outcome_products_random(struct wall *w,
                            struct vector3 *hitpt,
                            double t,
                            struct rxn *rx,
                            int path,
                            struct abstract_molecule *reacA,
                            struct abstract_molecule *reacB,
                            short orientA,
                            short orientB);
static int outcome_products_trimol_reaction_random(struct wall *w,
                            struct vector3 *hitpt,
                            double t,
                            struct rxn *rx,
                            int path,
                            struct abstract_molecule *reacA,
                            struct abstract_molecule *reacB,
                            struct abstract_molecule *reacC,
                            short orientA,
                            short orientB,
                            short orientC);



extern struct volume *world;

static int reaction_wizardry(struct magic_list *incantation,
                             struct wall *surface,
                             struct vector3 *hitpt,
                             double t);

int is_compatible_surface(void *req_species,struct wall *w)
{
  struct surf_class_list *scl, *scl2;

  struct surf_class_list *rs_head = (struct surf_class_list *)req_species;

  if (rs_head == NULL) return 1;

  for(scl = w->surf_class_head; scl != NULL; scl = scl->next)
  {
    for(scl2 = rs_head; scl2 != NULL; scl2 = scl2->next)
    {
      if(scl->surf_class == scl2->surf_class) return 1;
    }
  }

  return 0;
}

enum {
  PLAYER_GRID_MOL = 'g',
  PLAYER_VOL_MOL  = 'm',
  PLAYER_WALL     = 'w',
  PLAYER_NONE     = '\0',
  PLAYER_INVALID  = '!'
};

enum {
  PRODUCT_FLAG_NOT_SET,
  PRODUCT_FLAG_USE_UV_LOC,
  PRODUCT_FLAG_USE_REACA_UV,
  PRODUCT_FLAG_USE_REACB_UV,
  PRODUCT_FLAG_USE_REACC_UV,
  PRODUCT_FLAG_USE_RANDOM
};

#define IS_GRID_MOL(g) ((g) != NULL  &&  ((g)->properties->flags & ON_GRID))

static void add_players_to_list(struct rxn *rx,
                                struct abstract_molecule *reacA,
                                struct abstract_molecule *reacB,
                                struct abstract_molecule *reacC,
                                struct abstract_molecule **player,
                                char                      *player_type)
{
  /* Add reacA to the list of players, saving the reactant's type. */
  player[0] = reacA;
  player_type[0] = IS_GRID_MOL(reacA) ? PLAYER_GRID_MOL : PLAYER_VOL_MOL;

  /* If we have a second reactant, add it to the list of players. */
  if (rx->n_reactants > 1)
  {

    /* If the second reactant is a wall... */
    if (reacB == NULL)
    {
       assert(rx->n_reactants == 2);
       player[1] = NULL;
       player_type[1] = PLAYER_WALL;
    }

    /* Else, the second reactant is a molecule */
    else
    {
       player[1] = reacB;
       player_type[1] = IS_GRID_MOL(reacB) ? PLAYER_GRID_MOL : PLAYER_VOL_MOL;
    }

    if (rx->n_reactants > 2)
    {
        if(reacC == NULL)
        {
          /* it's a wall. */
          player[2] = NULL;
          player_type[2] = PLAYER_WALL;
        }else{
          player[2] = reacC;
          player_type[2] = IS_GRID_MOL(reacC) ? PLAYER_GRID_MOL : PLAYER_VOL_MOL;
        }
    }
  }
}

static bool is_rxn_unimol(struct rxn *rx)
{
  if (rx->n_reactants == 1)
    return true;

  if (rx->n_reactants != 2)
    return false;

  if (! (rx->players[0]->flags & ON_GRID))
    return false;

  return (rx->players[1]->flags & IS_SURFACE) != 0;
}

static struct volume_molecule *place_volume_subunit(struct species *product_species,
                                                    struct volume_molecule *old_volume_mol,
                                                    double t)
{
  /* Make sure the new molecule is of a different species.  Otherwise, why bother? */
  assert(old_volume_mol->properties != product_species);

  /* Find who we're replacing */
  int const subunit_idx = macro_subunit_index((struct abstract_molecule *) old_volume_mol);

  /* Find the species of our complex */
  struct complex_species * c_species = (struct complex_species *) (old_volume_mol->cmplx[0]->properties);
  int const num_subunits_in_complex = c_species->num_subunits;

  /* Allocate and initialize the molecule. */
  struct storage *local = old_volume_mol->subvol->local_storage;
  struct volume_molecule *new_volume_mol;
  new_volume_mol = CHECKED_MEM_GET(local->mol, "volume molecule");
  new_volume_mol->birthplace = local->mol;
  new_volume_mol->birthday = t;
  new_volume_mol->id = world->current_mol_id++;
  new_volume_mol->t = t;
  new_volume_mol->t2 = 0.0;
  new_volume_mol->properties = product_species;
  new_volume_mol->cmplx = NULL;
  new_volume_mol->prev_v = NULL;
  new_volume_mol->next_v = NULL;
  new_volume_mol->pos = old_volume_mol->pos;
  new_volume_mol->subvol = old_volume_mol->subvol;
  new_volume_mol->flags = TYPE_3D | ACT_NEWBIE | IN_VOLUME | IN_SCHEDULE | COMPLEX_MEMBER;
  /* Do not set diffuse flag since subunits are stationary. */

  if ((product_species->flags & COUNT_SOME_MASK) != 0)
    new_volume_mol->flags |= COUNT_ME;

  /* As a macromol rxn, this cannot be the result of a grid reaction.  Either
   * handle volume reversibility, or clear reversibility state. */
  new_volume_mol->previous_wall = NULL;
  if (world->volume_reversibility)
  {
    new_volume_mol->index = world->dissociation_index;
    new_volume_mol->flags |= ACT_CLAMPED;
  }
  else
  {
    new_volume_mol->index = -1;
  }

  /* Connect up new subunit to complex */
  old_volume_mol->cmplx[subunit_idx + 1] = new_volume_mol;
  new_volume_mol->cmplx = old_volume_mol->cmplx;

  /* Bind subunit to old molecule position */
  new_volume_mol->pos = old_volume_mol->pos;
  new_volume_mol->subvol = old_volume_mol->subvol;

  /* Add new molecule to the appopriate subvolume */
  ht_add_molecule_to_list(&new_volume_mol->subvol->mol_by_species, new_volume_mol);
  new_volume_mol->subvol->mol_count++;

  /* Update counts for this complex. */
  if (count_complex(old_volume_mol->cmplx[0], old_volume_mol, subunit_idx))
    mcell_allocfailed("Failed to update counts for complex subunits after a reaction.");

  /* Decide which, if any, subunits need to have unimol rxn times recomputed. */
  int update_subunit[ num_subunits_in_complex ];
  macro_count_inverse_related_subunits(c_species, update_subunit, subunit_idx);
  update_subunit[subunit_idx] = 0;

  /* Iterate over subunits, rescheduling as needed. */
  for (int other_subunit_idx = 1;
       other_subunit_idx <= num_subunits_in_complex;
       ++ other_subunit_idx)
  {
    /* No reschedule required for unrelated subunits. */
    if (! update_subunit[other_subunit_idx - 1])
      continue;

    /* If the subunit doesn't do unimol. rxns, no reschedule required. */
    struct volume_molecule *related_subunit = new_volume_mol->cmplx[other_subunit_idx];
    if (! (related_subunit->flags & ACT_REACT))
      continue;

    /* Set unimol rxn time to now, flag mol for recomputation. */
    related_subunit->t2 = 0.0;
    related_subunit->flags |= ACT_CHANGE;

    /* If the molecule is in the schedule, we'll have to find and excise it. */
    if (related_subunit->flags & IN_SCHEDULE)
      schedule_reschedule(related_subunit->subvol->local_storage->timer,
                          related_subunit, t);
  }

  /* Detach the old subunit from the complex. */
  old_volume_mol->cmplx = NULL;

  /* Check whether the product can undergo unimolecular rxns; if so, mark it. */
  /* N.B. This must occur after we've been added to the complex or we might not
   * match properly. */
  if (trigger_unimolecular(product_species->hashval,
                           (struct abstract_molecule*) new_volume_mol) != NULL)
    new_volume_mol->flags |= ACT_REACT;

  /* Add to the schedule. */
  if (schedule_add(local->timer, new_volume_mol))
    mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                      product_species->sym->name);
  return new_volume_mol;
}


static struct volume_molecule *place_volume_product(struct species *product_species,
                                                    struct grid_molecule *grid_reactant,
                                                    struct wall *w,
                                                    struct subvolume *subvol,
                                                    struct vector3 *hitpt,
                                                    short orient,
                                                    double t)
{
  struct vector3 pos = *hitpt;

  /* For an orientable reaction, we need to bump products out from the surface
   * to ensure they end up on the correct side of the plane. */
  if (w)
  {
    /* Note: no raytracing here so it is rarely possible to jump through closely spaced surfaces */
    double bump = (orient > 0) ? EPS_C : -EPS_C;
    pos.x += bump*w->normal.x;
    pos.y += bump*w->normal.y;
    pos.z += bump*w->normal.z;
    subvol = find_subvolume(&pos, subvol);
  }

  /* Allocate and initialize the molecule. */
  struct volume_molecule *new_volume_mol;
  new_volume_mol = CHECKED_MEM_GET(subvol->local_storage->mol, "volume molecule");
  new_volume_mol->birthplace = subvol->local_storage->mol;
  new_volume_mol->birthday = t;
  new_volume_mol->id = world->current_mol_id++;
  new_volume_mol->t = t;
  new_volume_mol->t2 = 0.0;
  new_volume_mol->properties = product_species;
  new_volume_mol->cmplx = NULL;
  new_volume_mol->prev_v = NULL;
  new_volume_mol->next_v = NULL;
  new_volume_mol->pos = pos;
  new_volume_mol->subvol = subvol;
  new_volume_mol->index = 0;
  new_volume_mol->flags = TYPE_3D | ACT_NEWBIE | IN_VOLUME | IN_SCHEDULE;
  if (product_species->space_step > 0.0)
    new_volume_mol->flags |= ACT_DIFFUSE;
  if ((product_species->flags & COUNT_SOME_MASK) != 0)
    new_volume_mol->flags |= COUNT_ME;

  /* Check whether the product can undergo unimolecular rxns; if so, mark it. */
  if (trigger_unimolecular(product_species->hashval,
                           (struct abstract_molecule*) new_volume_mol) != NULL)
    new_volume_mol->flags |= ACT_REACT;

  /* If this product resulted from a surface rxn, store the previous wall position. */
  if (grid_reactant)
  {
    new_volume_mol->previous_wall = grid_reactant->grid->surface;

    /* This will be overwritten with orientation in the CLAMPED/surf.
     * reversibility case
     */
    new_volume_mol->index = grid_reactant->grid_index;
  }

  /* Else clear the previous wall position. */
  else
  {
    new_volume_mol->previous_wall = NULL;
    new_volume_mol->index = -1;
  }

  /* Set reversibility state for the new molecule. */
  if (w)
  {
    if (world->surface_reversibility)
    {
      /* Which direction did we move? */
      new_volume_mol->previous_wall = w;
      new_volume_mol->index = (orient > 0) ? 1 : -1;
      new_volume_mol->flags |= ACT_CLAMPED;
    }
  }
  else if (world->volume_reversibility)
  {
    new_volume_mol->index = world->dissociation_index;
    new_volume_mol->flags |= ACT_CLAMPED;
  }

  /* Add the molecule to the subvolume */
  ht_add_molecule_to_list(&new_volume_mol->subvol->mol_by_species, new_volume_mol);
  ++ new_volume_mol->subvol->mol_count;

  /* Add to the schedule. */
  if (schedule_add(subvol->local_storage->timer, new_volume_mol))
    mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                      product_species->sym->name);
  return new_volume_mol;
}

static struct grid_molecule *place_grid_subunit(struct species *product_species,
                                                struct grid_molecule *old_grid_mol,
                                                struct surface_grid *grid,
                                                int grid_index,
                                                short orient,
                                                double t)
{
  /* Make sure the new molecule is of a different species.  Otherwise, why bother? */
  assert(old_grid_mol->properties != product_species  ||  old_grid_mol->orient != orient);

  /* Find who we're replacing */
  int const subunit_idx = macro_subunit_index((struct abstract_molecule *) old_grid_mol);

  /* Find the species of our complex */
  struct complex_species * c_species = (struct complex_species *) (old_grid_mol->cmplx[0]->properties);
  int const num_subunits_in_complex = c_species->num_subunits;

  /* Allocate and initialize the molecule. */
  struct grid_molecule *new_grid_mol;
  new_grid_mol = CHECKED_MEM_GET(old_grid_mol->birthplace, "grid molecule");
  new_grid_mol->birthplace = old_grid_mol->birthplace;
  new_grid_mol->birthday = t;
  new_grid_mol->id = world->current_mol_id++;
  new_grid_mol->t = t;
  new_grid_mol->t2 = 0.0;
  new_grid_mol->properties = product_species;
  new_grid_mol->cmplx = NULL;
  new_grid_mol->flags = TYPE_GRID | ACT_NEWBIE | IN_SCHEDULE | COMPLEX_MEMBER;
  /* Do not set "diffuse" flag since subunits are, at present, stationary. */
  if (product_species->flags & COUNT_ENCLOSED)
    new_grid_mol->flags |= COUNT_ME;
  new_grid_mol->grid = grid;
  new_grid_mol->grid_index = grid_index;
  new_grid_mol->s_pos = old_grid_mol->s_pos;
  new_grid_mol->orient = orient;

  /* Connect up new subunit to complex */
  old_grid_mol->cmplx[subunit_idx + 1] = new_grid_mol;
  new_grid_mol->cmplx = old_grid_mol->cmplx;

  /* If the complex reactant changed state, update counts for this complex. */
  struct grid_molecule *old_complex = old_grid_mol->cmplx[0];
  if (count_complex_surface(old_complex, old_grid_mol, subunit_idx))
    mcell_allocfailed("Failed to update region counts for surface macromolecule subunit '%s/%s' after a reaction.",
                      old_complex->properties->sym->name,
                      product_species->sym->name);

  /* Find out which subunits may need to recompute unimolecular rxn times. */
  int update_subunit[ num_subunits_in_complex ];
  macro_count_inverse_related_subunits(c_species, update_subunit, subunit_idx);
  update_subunit[subunit_idx] = 0;

  /* If the complex reactant changed state, reschedule unimolecular rxns. */
  struct subvolume *last_subvol = NULL;
  struct vector3 pos3d;
  for (int this_subunit_idx = 1; this_subunit_idx <= num_subunits_in_complex; ++ this_subunit_idx)
  {
    if (! update_subunit[this_subunit_idx - 1])
      continue;

    struct grid_molecule *this_subunit = new_grid_mol->cmplx[this_subunit_idx];
    if (! (this_subunit->flags & ACT_REACT))
      continue;

    /* Flag this molecule for recomputation of unimolecular rate. */
    this_subunit->t2 = 0.0;
    this_subunit->flags |= ACT_CHANGE;
    if (this_subunit->flags & IN_SCHEDULE)
    {
      uv2xyz(&this_subunit->s_pos, this_subunit->grid->surface, &pos3d);
      last_subvol = find_subvolume(&pos3d, last_subvol);
      schedule_reschedule(last_subvol->local_storage->timer, this_subunit, t);
    }
  }

  /* Check whether the product can undergo unimolecular rxns; if so, mark it. */
  if (trigger_unimolecular(product_species->hashval, (struct abstract_molecule*) new_grid_mol) != NULL ||
      (product_species->flags&CAN_GRIDWALL) != 0)
    new_grid_mol->flags |= ACT_REACT;

  /* Add to the grid. */
  ++ grid->n_occupied;
  grid->mol[ grid_index ] = new_grid_mol;

  /* Add to the schedule. */
  uv2xyz(& new_grid_mol->s_pos, new_grid_mol->grid->surface, & pos3d);
  last_subvol = find_subvolume(&pos3d, last_subvol);
  if (schedule_add(last_subvol->local_storage->timer, new_grid_mol))
    mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                      new_grid_mol->properties->sym->name);

  return new_grid_mol;
}

static struct grid_molecule *place_grid_product(struct species *product_species,
                                                struct surface_grid *grid,
                                                int grid_index,
                                                struct vector2 *mol_uv_pos,
                                                short orient,
                                                double t)
{
  struct vector3 mol_xyz_pos;
  uv2xyz(mol_uv_pos, grid->surface, & mol_xyz_pos);
  struct subvolume *sv = find_subvolume(& mol_xyz_pos, grid->subvol);

  /* Allocate and initialize the molecule. */
  struct grid_molecule *new_grid_mol;
  new_grid_mol = CHECKED_MEM_GET(sv->local_storage->gmol, "grid molecule");
  new_grid_mol->birthplace = sv->local_storage->gmol;
  new_grid_mol->birthday = t;
  new_grid_mol->id = world->current_mol_id++;
  new_grid_mol->t = t;
  new_grid_mol->t2 = 0.0;
  new_grid_mol->properties = product_species;
  new_grid_mol->cmplx = NULL;
  new_grid_mol->flags = TYPE_GRID | ACT_NEWBIE | IN_SCHEDULE;
  if (product_species->space_step > 0)
    new_grid_mol->flags |= ACT_DIFFUSE;
  if (product_species->flags & COUNT_ENCLOSED)
    new_grid_mol->flags |= COUNT_ME;
  new_grid_mol->grid = grid;
  new_grid_mol->grid_index = grid_index;
  new_grid_mol->s_pos = *mol_uv_pos;
  new_grid_mol->orient = orient;

  /* Check whether the product can undergo unimolecular rxns; if so, mark it. */
  if (trigger_unimolecular(product_species->hashval, (struct abstract_molecule*) new_grid_mol) != NULL ||
      (product_species->flags&CAN_GRIDWALL) != 0)
    new_grid_mol->flags |= ACT_REACT;

  /* Add to the grid. */
  ++ grid->n_occupied;
  grid->mol[ grid_index ] = new_grid_mol;

  /* Add to the schedule. */
  if (schedule_add(sv->local_storage->timer, new_grid_mol))
    mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                      product_species->sym->name);


  return new_grid_mol;
}

/***************************************************************************
outcome_products:
   In: first wall in the reaction
       hit point (if any)
       time of the reaction
       reaction
       path of the reaction
       first reactant (moving molecule)
       second reactant
       orientation of the first reactant
       orientation of the second reactant
   Out: Returns RX_A_OK, RX_FLIP or RX_BLOCKED.
Note: This function replaces surface reactants (if needed) by the surface
       products picked in the deterministic order from the list of products.
       Also surface products are placed in the deterministic order
       in the surrounding empty tiles.
       After this function execution some walls that do not have surface
       molecules and therefore do not have a grid may get a grid as side
       effect of calling functions
       "grid_all_neigbors_across_walls_through_vertices()"
       and "grid_all_neighbors_across_walls_through_edges()".
****************************************************************************/
static int outcome_products(struct wall *w,
                            struct vector3 *hitpt,
                            double t,
                            struct rxn *rx,
                            int path,
                            struct abstract_molecule *reacA,
                            struct abstract_molecule *reacB,
                            short orientA,
                            short orientB)
{
  bool update_dissociation_index = false;         /* Do we need to advance the dissociation index? */
  bool cross_wall = false;                        /* Did the moving molecule cross the plane? */
  struct subvolume *last_subvol = NULL;           /* Last subvolume (guess used to speed sv finding) */

  int const i0 = rx->product_idx[path];           /* index of the first player for the pathway */
  int const iN = rx->product_idx[path+1];         /* index of the first player for the next pathway */
  assert(iN > i0);
  struct species **rx_players = rx->players + i0; /* Players array from the reaction. */

  int const n_players = iN - i0;                  /* number of reaction players */
  struct abstract_molecule *product[n_players];   /* array of products */
  char product_type[n_players];                   /* array that decodes the type of each product */
  short product_orient[n_players];                /* array of orientations for each product */
  struct surface_grid *product_grid[n_players];   /* array of surface_grids for products */
  int product_grid_idx[n_players];                /* array of grid indices for products */
  byte product_flag[n_players];                   /* array of placement flags for products */

  struct abstract_molecule *old_subunit = NULL;   /* Pointer to reactant which was a subunit, if any. */

  bool const is_unimol = is_rxn_unimol(rx);       /* Unimol rxn (not mol-mol, not mol-wall) */

  /* Clear the initial product info. */
  for (int i=0; i<n_players; ++i)
  {
    product[i] = NULL;
    product_type[i] = PLAYER_NONE;
    product_orient[i] = 0;
    product_grid[i] = NULL;
    product_grid_idx[i] = -1;
    product_flag[i] = PRODUCT_FLAG_NOT_SET;
  }

  /* Flag indicating that a surface is somehow involved with this reaction. */
  struct grid_molecule * const grid_1 = IS_GRID_MOL(reacA) ? (struct grid_molecule *) reacA : NULL;
  struct grid_molecule * const grid_2 = IS_GRID_MOL(reacB) ? (struct grid_molecule *) reacB : NULL;
  struct grid_molecule * const grid_reactant = grid_1 ? grid_1 : grid_2;
  bool const is_orientable = (w != NULL)  ||  (grid_reactant != NULL);

  /* reacA is the molecule which initiated the reaction. */
  struct abstract_molecule * const initiator = reacA;
  short const initiatorOrient = orientA;

  /* Ensure that reacA and reacB are sorted in the same order as the rxn players. */
  assert(reacA != NULL);
  if (reacA->properties != rx->players[0])
  {
    struct abstract_molecule *tmp_mol = reacA;
    reacA = reacB;
    reacB = tmp_mol;

    short tmp_orient = orientA;
    orientA = orientB;
    orientB = tmp_orient;
  }
  assert(reacA != NULL);

  /* Add the reactants (incl. any wall) to the list of players. */
  add_players_to_list(rx, reacA, reacB, NULL, product, product_type);

  /* If the reaction is complex, figure out which reactant is a subunit. */
  if (rx->is_complex)
  {
    if (reacA->flags & COMPLEX_MEMBER)
      old_subunit = reacA;
    else if (reacB != NULL  &&  reacB->flags & COMPLEX_MEMBER)
      old_subunit = reacB;
    else if (reacB != NULL)
      mcell_internal_error("Macromolecular reaction [%s] occurred, but neither molecule is a subunit (%s and %s).",
                           rx->sym->name, reacA->properties->sym->name, reacB->properties->sym->name);
    else
      mcell_internal_error("Macromolecular reaction [%s] occurred, but the molecule is not a subunit (%s).",
                           rx->sym->name, reacA->properties->sym->name);
  }

  /* If the reaction involves a surface, make sure there is room for each product. */
  struct vector2 rxn_uv_pos;
  if (is_orientable)
  {
    /* Determine whether any of the reactants can be replaced by a product. */
    int replace_p1 = (product_type[0] == PLAYER_GRID_MOL  &&  rx_players[0] == NULL);
    int replace_p2 = rx->n_reactants > 1  &&  (product_type[1] == PLAYER_GRID_MOL  &&  rx_players[1] == NULL);
    assert(! replace_p2  ||  reacB != NULL);

    /* Determine the point of reaction on the surface. */
    if (grid_reactant) rxn_uv_pos = grid_reactant->s_pos;
    else xyz2uv(hitpt, w, &rxn_uv_pos);

    /* For each product, find a position. */
    int last_placed = -1;
    struct grid_molecule sentinel;
    for (int n_product = 0; n_product < n_players; ++ n_product)
    {
      /* Skip NULL reactants. */
      if (rx_players[n_product] == NULL)
        continue;

      int this_geometry = rx->geometries[i0 + n_product];

      /* Geometry of 0 means "random orientation" */
      if (this_geometry == 0)
        product_orient[n_product] = (rng_uint(world->rng) & 1) ? 1 : -1;
      else
      {
        /* Geometry < 0 means inverted orientation */
        if (this_geometry < 0)
        {
          this_geometry = -this_geometry;
          if (this_geometry > (int) rx->n_reactants)
            product_orient[n_product] = - product_orient[this_geometry - rx->n_reactants - 1];
          else if (this_geometry == 1)
            product_orient[n_product] = - orientA;
          else if (this_geometry == 2 && reacB != NULL)
            product_orient[n_product] = - orientB;
          else
            product_orient[n_product] = -1;
        }

        /* Geometry > 0 means "positive" orientation. */
        else
        {
          if (this_geometry > (int) rx->n_reactants)
            product_orient[n_product] = product_orient[this_geometry - rx->n_reactants - 1];
          else if (this_geometry == 1)
            product_orient[n_product] = orientA;
          else if (this_geometry == 2 && reacB != NULL)
            product_orient[n_product] = orientB;
          else
            product_orient[n_product] = 1;
        }
      }

      /* If this is a reactant... */
      if (n_product < (int) rx->n_reactants)
      {
        /* If this is a grid molecule, we need to set its orientation. */
        if (rx_players[n_product]->flags & ON_GRID)
        {
          assert(IS_GRID_MOL(product[n_product]));
          struct grid_molecule *gm = (struct grid_molecule *) product[n_product];

          /* If the new orientation doesn't match the old, we've got some work to do. */
          if (gm->orient != product_orient[n_product])
          {
            int const subunit_idx = old_subunit ? macro_subunit_index((struct abstract_molecule *) gm) : -1;
            struct grid_molecule gm_old = *gm;

            /* We're about to update the molecule's orientation, so we will
             * first remove it from the counts in case we have any
             * orientation-sensitive counts.  Then, we will update the
             * orientation.  Finally, we will add the molecule back into the
             * counts in its new orientation.
             */

            /* Remove molecule from counts in old orientation, if mol is counted. */
            if (product[n_product]->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
              count_region_from_scratch(product[n_product],     /* molecule */
                                        NULL,                   /* rxn pathway */
                                        -1,                     /* remove count */
                                        NULL,                   /* Location at which to count */
                                        w,                      /* Wall on which this happened */
                                        t);                     /* Time of occurrence */

            /* Force check for the unimolecular reactions
               after changing orientation.
               There are two possible cases to be covered here:
               1) when (gm->t2) was previously set to FOREVER
               2) there may be two or more unimolecular
                  reactions involving surface class that have
                  different kinetics.
             */
            if(((gm->flags & ACT_REACT) != 0) && ((gm->properties->flags & CAN_GRIDWALL) != 0))  gm->t2 = 0;

            /* Set the molecule's orientation. */
            gm->orient = product_orient[n_product];

            /* Add molecule back to counts in new orientation, if mol is counted. */
            if (product[n_product]->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
              count_region_from_scratch(product[n_product],     /* molecule */
                                        NULL,                   /* rxn pathway */
                                        1,                      /* add count */
                                        NULL,                   /* Location at which to count */
                                        w,                      /* Wall on which this happened */
                                        t);                     /* Time of occurrence */

            /* Update macromolecular counts. */
            if (old_subunit  &&  count_complex_surface(gm->cmplx[0], & gm_old, subunit_idx))
              mcell_allocfailed("Failed to update region counts for surface macromolecule subunit '%s/%s' after a reaction.",
                                gm->cmplx[0]->properties->sym->name,
                                gm->properties->sym->name);
          }
        }

        /* Otherwise, check if we've crossed the plane. */
        else if (! is_unimol)
        {
          if (product[n_product] == initiator)
          {
            if (product_orient[n_product] != initiatorOrient)
              cross_wall = true;
          }
        }

        /* Skip placement for this molecule -- we're already placed.  */
        continue;
      }

      /* If the product is a volume product, no placement is required. */
      if (rx_players[n_product]->flags & ON_GRID)
      {
        /* If the first reactant should be replaced... */
        if (replace_p1  &&  (! replace_p2 || initiator == reacA))
        {
          product_grid[n_product]     = ((struct grid_molecule *) reacA)->grid;
          product_grid_idx[n_product] = ((struct grid_molecule *) reacA)->grid_index;
          product_flag[n_product]     = PRODUCT_FLAG_USE_REACA_UV;
          replace_p1 = 0;
        }

        /* Else if the second reactant (in rxn order) can be replaced, replace it. */
        else if (replace_p2)
        {
          product_grid[n_product]     = ((struct grid_molecule *) reacB)->grid;
          product_grid_idx[n_product] = ((struct grid_molecule *) reacB)->grid_index;
          product_flag[n_product]     = PRODUCT_FLAG_USE_REACB_UV;
          replace_p2 = 0;
        }

        /* Else, we'll need to place in a new location. */
        else
        {
          /* If the grid is nonexistent, create it and place the molecule. */
          assert(w != NULL);
          if (w->grid==NULL)
          {
            /* reacA must be a volume molecule, or this wall would have a grid already. */
            assert(! IS_GRID_MOL(reacA));

            if (create_grid(w, ((struct volume_molecule *) reacA)->subvol))
              mcell_allocfailed("Failed to create a grid for a wall.");

            /* This spot is empty because we just created the grid. */
            product_grid[n_product]     = w->grid;
            product_grid_idx[n_product] = uv2grid(& rxn_uv_pos, w->grid);
            product_flag[n_product]     = PRODUCT_FLAG_USE_UV_LOC;
            last_placed = n_product;
          }

          /* Else, search for a place to put the molecule. */
          else
          {
            /* Mark the last placed molecule slot as occupied. */
            if (last_placed >= 0)
              product_grid[last_placed]->mol[ product_grid_idx[last_placed] ] = & sentinel;

            /* If we've placed no molecule yet, and the desired spot is free, place product there. */
            int desired_pos;
            struct wall *desired_wall = NULL;
            if (last_placed < 0  &&  w->grid->mol[desired_pos = uv2grid(& rxn_uv_pos, w->grid)] == NULL)
            {
              product_grid[n_product]     = w->grid;
              product_grid_idx[n_product] = desired_pos;
              product_flag[n_product]     = PRODUCT_FLAG_USE_UV_LOC;
              last_placed = n_product;
            }

            /* Else if the vacancy search distance is non-zero, search nearby for a free spot. */
            else if (world->vacancy_search_dist2 > 0.0  &&
                     (desired_wall = search_nbhd_for_free(w,
                                                          & rxn_uv_pos,
                                                          world->vacancy_search_dist2,
                                                          & desired_pos,
                                                          &is_compatible_surface,
                                                          (void *) w->surf_class_head)) != NULL)
            {
              product_grid[n_product]     = desired_wall->grid;
              product_grid_idx[n_product] = desired_pos;
              product_flag[n_product]     = PRODUCT_FLAG_USE_RANDOM;
              last_placed = n_product;
            }

            /* If a spot isn't found, clean up and fail to react. */
            else
            {
              for (int n_placed = rx->n_reactants; n_placed < n_product; ++ n_placed)
              {
                if (product_grid[n_placed] == NULL)
                  continue;
                if (product_grid[n_placed]->mol[ product_grid_idx[n_placed] ] == & sentinel)
                  product_grid[n_placed]->mol[ product_grid_idx[n_placed] ] = NULL;
              }

              return RX_BLOCKED;
            }
          }
        }
      }
    }
  }

  /* Determine the location of the reaction for count purposes. */
  struct vector3 count_pos_xyz;
  if (hitpt != NULL) count_pos_xyz = *hitpt;
  else if (grid_reactant) uv2xyz(& grid_reactant->s_pos, grid_reactant->grid->surface, & count_pos_xyz);
  else count_pos_xyz = ((struct volume_molecule *) reacA)->pos;

  /* Create and place each product. */
  struct vector3 mol_pos_tmp;
  struct subvolume *product_subvol = NULL;
  for (int n_product = rx->n_reactants; n_product < n_players; ++ n_product)
  {
    struct abstract_molecule *this_product = NULL;
    struct species * const product_species = rx_players[n_product];

    bool const product_is_subunit = (old_subunit != NULL  &&  rx->is_complex[i0 + n_product]);

    /* If the product is a grid molecule, place it on the grid. */
    if (product_species->flags & ON_GRID)
    {
      if (! product_is_subunit)
      {
        struct vector2 prod_uv_pos;

        /* Pick an appropriate position for the new molecule. */
        if (world->randomize_gmol_pos)
        {
          switch (product_flag[n_product])
          {
            case PRODUCT_FLAG_USE_REACA_UV:
              assert(reacA != NULL);
              prod_uv_pos = ((struct grid_molecule*) reacA)->s_pos;
              break;

            case PRODUCT_FLAG_USE_REACB_UV:
              assert(reacB != NULL);
              prod_uv_pos = ((struct grid_molecule*) reacB)->s_pos;
              break;

            case PRODUCT_FLAG_USE_UV_LOC:
              prod_uv_pos = rxn_uv_pos;
              break;

            case PRODUCT_FLAG_USE_RANDOM:
              grid2uv_random(product_grid[n_product], product_grid_idx[n_product], &prod_uv_pos);
              break;

            default:
              UNHANDLED_CASE(product_flag[n_product]);
              break;
          }
        }
        else
          grid2uv(product_grid[n_product], product_grid_idx[n_product], & prod_uv_pos);

        this_product = (struct abstract_molecule *)
              place_grid_product(product_species,
                                 product_grid[n_product],
                                 product_grid_idx[n_product],
                                 & prod_uv_pos,
                                 product_orient[n_product],
                                 t);
      }
      else
        this_product = (struct abstract_molecule *)
              place_grid_subunit(product_species,
                                 (struct grid_molecule *) old_subunit,
                                 product_grid[n_product],
                                 product_grid_idx[n_product],
                                 product_orient[n_product],
                                 t);
    }

    /* else place the molecule in space. */
    else
    {
      if (! product_is_subunit)
      {
        /* Unless this is a unimolecular reaction, we will have a hitpt. */
        if (! hitpt)
        {
          /* If this is a unimolecular surface rxn... */
          if (reacA->properties->flags & ON_GRID)
          {
            uv2xyz(& ((struct grid_molecule *) reacA)->s_pos,
                   ((struct grid_molecule *) reacA)->grid->surface,
                   & mol_pos_tmp);
            product_subvol = find_subvolume(&mol_pos_tmp, last_subvol);
          }

          /* ... else a unimolecular volume rxn. */
          else
          {
            mol_pos_tmp = ((struct volume_molecule *) reacA)->pos;
            product_subvol = ((struct volume_molecule *) reacA)->subvol;
          }
          hitpt = & mol_pos_tmp;
        }
        else if (product_subvol == NULL)
          product_subvol = find_subvolume(hitpt, last_subvol);

        this_product = (struct abstract_molecule *) place_volume_product(product_species,
                                                                         grid_reactant,
                                                                         w,
                                                                         product_subvol,
                                                                         hitpt,
                                                                         product_orient[n_product],
                                                                         t);
      }
      else
        this_product = (struct abstract_molecule *) place_volume_subunit(product_species,
                                                                         (struct volume_molecule *) old_subunit,
                                                                         t);

      if (((struct volume_molecule *) this_product)->index < DISSOCIATION_MAX)
        update_dissociation_index = true;
    }

    /* Update molecule counts */
    ++ product_species->population;
    if (product_species->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
      count_region_from_scratch(this_product, NULL, 1, NULL, NULL, t);
  }

  /* If necessary, update the dissociation index. */
  if (update_dissociation_index)
  {
    if (-- world->dissociation_index < DISSOCIATION_MIN)
      world->dissociation_index = DISSOCIATION_MAX;
  }

  /* Handle events triggered off of named reactions */
  if (rx->info[path].pathname!=NULL)
  {
    /* No flags for reactions so we have to check regions if we have waypoints! Fix to be more efficient for WORLD-only counts? */
    if (world->place_waypoints_flag)
      count_region_from_scratch(NULL, rx->info[path].pathname, 1, &count_pos_xyz, w, t);

    /* Other magical stuff.  For now, can only trigger releases. */
    if (rx->info[path].pathname->magic!=NULL)
    {
      if (reaction_wizardry(rx->info[path].pathname->magic, w, &count_pos_xyz, t))
        mcell_allocfailed("Failed to complete reaction triggered release after a '%s' reaction.",
                          rx->info[path].pathname->sym->name);
    }
  }

  return cross_wall ? RX_FLIP : RX_A_OK;
}


/***************************************************************************
outcome_products_random:
   In: first wall in the reaction
       hit point (if any)
       time of the reaction
       reaction
       path of the reaction
       first reactant (moving molecule)
       second reactant
       orientation of the first reactant
       orientation of the second reactant
   Out: Returns RX_A_OK, RX_FLIP or RX_BLOCKED.
Note: This function replaces surface reactants (if needed) by the surface
       products picked in the random order from the list of products.
       Also surface products are placed in the random order
       in the surrounding empty tiles.
       After this function execution some walls that do not have surface
       molecules and therefore do not have a grid may get a grid as side
       effect of calling functions
       "grid_all_neigbors_across_walls_through_vertices()"
       and "grid_all_neighbors_across_walls_through_edges()".
Note: If both reactants are grid molecules, and they are both located within
      the same restricted region border (REFL/ABSORB),
      then reaction products - grid molecules
      for which this region border is restrictive will be placed inside
      this region. If any of the conditions above are not true, the reaction
      products - grid molecules can be placed on any tile from the list
      of vacant tiles.
Note: Policy on surface products placement is described in the document
      "policy_surf_products_placement.doc" (see "src/docs").

****************************************************************************/
static int outcome_products_random(struct wall *w,
                            struct vector3 *hitpt,
                            double t,
                            struct rxn *rx,
                            int path,
                            struct abstract_molecule *reacA,
                            struct abstract_molecule *reacB,
                            short orientA,
                            short orientB)
{
  assert(!rx->is_complex);

  bool update_dissociation_index = false;         /* Do we need to advance the dissociation index? */
  bool cross_wall = false;                        /* Did the moving molecule cross the plane? */
  struct subvolume *last_subvol = NULL;           /* Last subvolume (guess used to speed sv finding) */

  int const i0 = rx->product_idx[path];           /* index of the first player for the pathway */
  int const iN = rx->product_idx[path+1];         /* index of the first player for the next pathway */
  assert(iN > i0);
  struct species **rx_players = rx->players + i0; /* Players array from the reaction. */

  int const n_players = iN - i0;                  /* number of reaction players */
  struct abstract_molecule *product[n_players];   /* array of products */
  char product_type[n_players];                   /* array that decodes the type of each product */
  short product_orient[n_players];                /* array of orientations for each product */
  struct surface_grid *product_grid[n_players];   /* array of surface_grids for products */
  int product_grid_idx[n_players];                /* array of grid indices for products */
  byte product_flag[n_players];                   /* array of placement flags for products */

  bool const is_unimol = is_rxn_unimol(rx);       /* Unimol rxn (not mol-mol, not mol-wall) */

  struct tile_neighbor *tile_nbr_head = NULL;  /* list of neighbor tiles */
  struct tile_neighbor *tile_nbr;  /* iterator */
  /* head of the linked list of vacant neighbor tiles */
  struct tile_neighbor *tile_vacant_nbr_head = NULL;
  struct surface_grid *tile_grid;  /* surface grid the tile belongs to */
  int tile_idx;    /* index of the tile on the grid */
  unsigned int rnd_num;  /* random number */
  int num_vacant_tiles = 0; /* number of vacant tiles */

  int num_surface_products = 0; /* not counting reactants */
  int num_surface_static_products = 0; /* number of products with (D_2D == 0) */
  int num_surface_static_reactants = 0; /* number of reactants with (D_2D == 0) */
  int list_length = 0; /* length of the linked list tile_nbr_head */

  /* used for product placement for the reaction of type A->B+C[rate] */
  unsigned int reac_idx = UINT_MAX, mol_idx = UINT_MAX;
  struct surface_grid *reac_grid = NULL, *mol_grid = NULL;
  struct vector2 rxn_uv_pos; /* position of the reaction */
  int rxn_uv_idx = -1;  /* tile index of the reaction place */

  /* Clear the initial product info. */
  for (int i=0; i<n_players; ++i)
  {
    product[i] = NULL;
    product_type[i] = PLAYER_NONE;
    product_orient[i] = 0;
    product_grid[i] = NULL;
    product_grid_idx[i] = -1;
    product_flag[i] = PRODUCT_FLAG_NOT_SET;
  }

  /* Flag indicating that a surface is somehow involved with this reaction. */
  struct grid_molecule * const grid_1 = IS_GRID_MOL(reacA) ? (struct grid_molecule *) reacA : NULL;
  struct grid_molecule * const grid_2 = IS_GRID_MOL(reacB) ? (struct grid_molecule *) reacB : NULL;
  struct grid_molecule * const grid_reactant = grid_1 ? grid_1 : grid_2;
  bool const is_orientable = (w != NULL)  ||  (grid_reactant != NULL);

  /* list of the restricted regions for the reactants by wall */
  struct region_list *rlp_head_wall_1 = NULL, *rlp_head_wall_2 = NULL;
  /* list of the restricted regions for the reactants by object */
  struct region_list *rlp_head_obj_1 = NULL, *rlp_head_obj_2 = NULL;

  int grid_bitmask = \
    determine_molecule_region_topology(grid_1, grid_2,
                                       &rlp_head_wall_1,
                                       &rlp_head_wall_2,
                                       &rlp_head_obj_1,
                                       &rlp_head_obj_2,
                                       is_unimol);

  /* reacA is the molecule which initiated the reaction. */
  struct abstract_molecule * const initiator = reacA;
  short const initiatorOrient = orientA;

  /* Ensure that reacA and reacB are sorted in the same order as the rxn players. */
  assert(reacA != NULL);

  if (reacA->properties != rx->players[0])
  {
    struct abstract_molecule *tmp_mol = reacA;
    reacA = reacB;
    reacB = tmp_mol;

    short tmp_orient = orientA;
    orientA = orientB;
    orientB = tmp_orient;
  }

  assert(reacA != NULL);

  /* Add the reactants (incl. any wall) to the list of players. */
  add_players_to_list(rx, reacA, reacB, NULL, product, product_type);

  /* Determine whether any of the reactants can be replaced by a product. */
  int replace_p1 = (product_type[0] == PLAYER_GRID_MOL  &&  rx_players[0] == NULL);
  int replace_p2 = rx->n_reactants > 1  &&  (product_type[1] == PLAYER_GRID_MOL  &&  rx_players[1] == NULL);

  /* Determine the point of reaction on the surface. */
  if(is_orientable)
  {
     if (grid_reactant) rxn_uv_pos = grid_reactant->s_pos;
     else
     {
       xyz2uv(hitpt, w, &rxn_uv_pos);
     }

     assert(w != NULL);
     if(w->grid == NULL)
     {
           /* reacA must be a volume molecule, or this wall would have a grid already. */
        assert(! IS_GRID_MOL(reacA));

        if (create_grid(w, ((struct volume_molecule *) reacA)->subvol))
               mcell_allocfailed("Failed to create a grid for a wall.");
     }
        /* create list of neighbor tiles around rxn_uv_pos */
     rxn_uv_idx = uv2grid(&rxn_uv_pos, w->grid);

     /* find out number of static surface reactants */
     if((grid_1 != NULL) && (grid_1->properties->D == 0)) num_surface_static_reactants++;
     if((grid_2 != NULL) && (grid_2->properties->D == 0)) num_surface_static_reactants++;
  }

   /* find out number of surface products */
  for(int n_product = rx->n_reactants; n_product < n_players; ++ n_product)
  {
     if(rx_players[n_product] == NULL) continue;
     if(rx_players[n_product]->flags & ON_GRID) {
         num_surface_products++;
         if(rx_players[n_product]->D == 0) num_surface_static_products++;
     }
  }

  /* If the reaction involves a surface, make sure there is room for each product. */
  if (is_orientable)
  {

    if(num_surface_products > 0)
    {

       if(grid_reactant != NULL)
       {
         find_neighbor_tiles(grid_reactant, grid_reactant->grid, grid_reactant->grid_index, 1, 0, &tile_nbr_head, &list_length);
       }else{
         find_neighbor_tiles(grid_reactant, w->grid, rxn_uv_idx, 1, 0, &tile_nbr_head, &list_length);
       }

       /* Create list of vacant tiles */
       for(tile_nbr = tile_nbr_head; tile_nbr != NULL; tile_nbr = tile_nbr->next)
       {
          if(tile_nbr->grid->mol[tile_nbr->idx] == NULL)
          {
             num_vacant_tiles++;
             push_tile_neighbor_to_list(&tile_vacant_nbr_head, tile_nbr->grid, tile_nbr->idx);
          }
        }
     }

    /* Can this reaction happen at all? */
    if(replace_p1 && replace_p2)
    {
       if(num_surface_products > num_vacant_tiles + 2)
       {
          if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
          if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
          return RX_BLOCKED;
       }
    }else if(replace_p1 || replace_p2) {
       if(num_surface_products > num_vacant_tiles + 1)
       {
          if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
          if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
          return RX_BLOCKED;
       }
    }else{
       if(num_surface_products > num_vacant_tiles)
       {
          if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
          if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
          return RX_BLOCKED;
       }
    }

    /* set the orientations of the products. */
    for (int n_product = 0; n_product < n_players; ++ n_product)
    {
      /* Skip NULL reactants. */
      if (rx_players[n_product] == NULL)
      {
        continue;
      }

      int this_geometry = rx->geometries[i0 + n_product];
      int relative_orient = (this_geometry < 0) ? -1 : 1;
      this_geometry = abs(this_geometry);

      /* Geometry of 0 means "random orientation" */
      if (this_geometry == 0)
      {
        product_orient[n_product] = (rng_uint(world->rng) & 1) ? 1 : -1;
      }
      else
      {
        if (this_geometry > (int) rx->n_reactants)
        {
          product_orient[n_product] = relative_orient *
            product_orient[this_geometry - rx->n_reactants - 1];
        }
        else if (this_geometry == 1)
        {
          product_orient[n_product] = relative_orient * orientA;
        }
        else if (this_geometry == 2 && reacB != NULL)
        {
          product_orient[n_product] = relative_orient * orientB;
        }
        else
        {
          product_orient[n_product] = relative_orient * 1;
        }
      }

      /* If this is a reactant... */
      if (n_product < (int) rx->n_reactants)
      {
        /* If this is a grid molecule, we need to set its orientation. */
        if (rx_players[n_product]->flags & ON_GRID)
        {
          assert(IS_GRID_MOL(product[n_product]));
          struct grid_molecule *gm = (struct grid_molecule *) product[n_product];

          /* If the new orientation doesn't match the old, we've got some work to do. */
          if (gm->orient != product_orient[n_product])
          {

            /* We're about to update the molecule's orientation, so we will
             * first remove it from the counts in case we have any
             * orientation-sensitive counts.  Then, we will update the
             * orientation.  Finally, we will add the molecule back into the
             * counts in its new orientation.
             */

            /* Remove molecule from counts in old orientation, if mol is counted. */
            if (product[n_product]->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
            {
              count_region_from_scratch(product[n_product],     /* molecule */
                                        NULL,                   /* rxn pathway */
                                        -1,                     /* remove count */
                                        NULL,                   /* Location at which to count */
                                        w,                      /* Wall on which this happened */
                                        t);                     /* Time of occurrence */
            }

            /* Force check for the unimolecular reactions
               after changing orientation.
               There are two possible cases to be covered here:
               1) when (gm->t2) was previously set to FOREVER
               2) there may be two or more unimolecular
                  reactions involving surface class that have
                  different kinetics.
             */
            if(((gm->flags & ACT_REACT) != 0) && ((gm->properties->flags & CAN_GRIDWALL) != 0))
            {
              gm->t2 = 0;
            }

            /* Set the molecule's orientation. */
            gm->orient = product_orient[n_product];

            /* Add molecule back to counts in new orientation, if mol is counted. */
            if (product[n_product]->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
            {
              count_region_from_scratch(product[n_product],     /* molecule */
                                        NULL,                   /* rxn pathway */
                                        1,                      /* add count */
                                        NULL,                   /* Location at which to count */
                                        w,                      /* Wall on which this happened */
                                        t);                     /* Time of occurrence */
            }
          }
        }

        /* Otherwise, check if we've crossed the plane. */
        else if (!is_unimol)
        {
          if (product[n_product] == initiator)
          {
            if (product_orient[n_product] != initiatorOrient)
              cross_wall = true;
          }
        }
      }
    }

    /* find out where to place surface products */
    /* Some special cases are listed below. */

    if(num_surface_products == 1)
    {
      if (is_unimol)
      {
        if (replace_p1)
        {
          for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
          {
            if(rx_players[n_product] == NULL) continue;
            if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;

            if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
            {
               product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
               product_grid[n_product] = ((struct grid_molecule *)reacA)->grid;
               product_grid_idx[n_product] = ((struct grid_molecule *)reacA)->grid_index;
               replace_p1 = 0;
               break;
            }
          }
        }
      }
      else
      {  /* this is bimol reaction */

        if((num_surface_static_reactants == 1) &&
            (num_surface_static_products == 1) &&
            (replace_p1 || replace_p2))
        {
          /* the lonely static product always replaces lonely static reactant */
          for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
          {
            if(rx_players[n_product] == NULL) continue;
            if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;
            if(rx_players[n_product]->D != 0) continue;

            if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
            {
              if(replace_p1 && (reacA->properties->D == 0))
              {
                  product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
                  product_grid[n_product] = ((struct grid_molecule *)reacA)->grid;
                  product_grid_idx[n_product] = ((struct grid_molecule *)reacA)->grid_index;
                  replace_p1 = 0;
                  break;
              }
              else if(replace_p2 && (reacB->properties->D == 0))
              {
                  product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
                  product_grid[n_product] = ((struct grid_molecule *)reacB)->grid;
                  product_grid_idx[n_product] = ((struct grid_molecule *)reacB)->grid_index;
                  replace_p2 = 0;
                  break;
              }

            }
          }
        }
        else if (replace_p1 && replace_p2)
        {

          /* if both reactants should be  replaced and there is only one
              surface product here we make sure that initiator molecule is
              replaced */
          for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
          {
              if(rx_players[n_product] == NULL) continue;
              if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;

              if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
              {
                if(reacA == initiator)
                {
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacA)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacA)->grid_index;
                    replace_p1 = 0;
                }else{
                    mcell_log("are here\n");
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacB)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacB)->grid_index;
                    replace_p2 = 0;
                }
                break;
              }
          }
        }
        else if (replace_p1)
        {
          for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
          {
            if(rx_players[n_product] == NULL) continue;
            if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;
            if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
            {
              if(replace_p1)
              {
                  product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
                  product_grid[n_product] = ((struct grid_molecule *)reacA)->grid;
                  product_grid_idx[n_product] = ((struct grid_molecule *)reacA)->grid_index;
                  replace_p1 = 0;
                  break;
              }
            }
          }
        }
        else if (replace_p2)
        {
          for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
          {
            if(rx_players[n_product] == NULL) continue;
            if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;

            if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
            {
              if(replace_p2)
              {
                  product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
                  product_grid[n_product] = ((struct grid_molecule *)reacB)->grid;
                  product_grid_idx[n_product] = ((struct grid_molecule *)reacB)->grid_index;
                  replace_p2 = 0;
                  break;
              }
            }
          }
        }
      }
    }
    else if (num_surface_products > 1)
    {
      /* more than one surface products */
      int count;
      if(num_surface_static_reactants > 0)
      {
        bool replace_reacA = (reacA->properties->D == 0) && replace_p1;
        bool replace_reacB = (reacB == NULL) ? false : (reacB->properties->D == 0) && replace_p2;

        if(replace_reacA || replace_reacB)
        {
          int max_static_count = (num_surface_static_products < num_surface_static_reactants) ? num_surface_static_products : num_surface_static_reactants;

          count = 0;

          while(count < max_static_count)
          {
            rnd_num = rng_uint(world->rng) % n_players;
            /* pass reactants */
            if (rnd_num < rx->n_reactants) continue;

            if(rx_players[rnd_num] == NULL) continue;
            if((rx_players[rnd_num]->flags & NOT_FREE) == 0) continue;
            if(rx_players[rnd_num]->D != 0) continue;

            if(product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET)
            {
              if(replace_reacA)
              {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                product_grid[rnd_num] = ((struct grid_molecule *)reacA)->grid;
                product_grid_idx[rnd_num] = ((struct grid_molecule *)reacA)->grid_index;
                count++;
                replace_p1 = 0;
                continue;
              }
              if(replace_reacB)
              {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                product_grid[rnd_num] = ((struct grid_molecule *)reacB)->grid;
                product_grid_idx[rnd_num] = ((struct grid_molecule *)reacB)->grid_index;
                count++;
                replace_p2 = 0;
                continue;
              }
            }
          } /* end while */
        }
      }

        /* check whether there are any surface reactants left to replace
            with surface products since we are done with static
            reactants/products */

      if(replace_p1 || replace_p2)
      {
        /* are there any surface products left that have not yet replaced
            reactants? */

        int surf_prod_left = 0, surf_reactant_left = 0;

        for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
        {
          if(rx_players[n_product] == NULL) continue;
          if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;

          if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET) surf_prod_left++;
        }

        if(replace_p1) surf_reactant_left++;
        if(replace_p2) surf_reactant_left++;

        if(surf_prod_left > 0)
        {
          if(surf_prod_left >= surf_reactant_left)
          {
            count = 0;
            while(count < surf_reactant_left)
            {
              rnd_num = rng_uint(world->rng) % n_players;
              /* pass reactants */
              if (rnd_num < 2) continue;

              if(rx_players[rnd_num] == NULL) continue;
              if((rx_players[rnd_num]->flags & NOT_FREE) == 0) continue;

              if(product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET)
              {
                if(replace_p1)
                {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                  product_grid[rnd_num] = ((struct grid_molecule *)reacA)->grid;
                  product_grid_idx[rnd_num] = ((struct grid_molecule *)reacA)->grid_index;
                  count++;
                  replace_p1 = 0;
                  continue;
                }
                if(replace_p2)
                {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                  product_grid[rnd_num] = ((struct grid_molecule *)reacB)->grid;
                  product_grid_idx[rnd_num] = ((struct grid_molecule *)reacB)->grid_index;
                  replace_p2 = 0;
                  count++;
                  continue;
                }
              }
            } /* end while */

          }
          else   /* surf_prod_left < surf_reactant_left */
          {
            count = 0;
            while(count < surf_prod_left)
            {
              rnd_num = rng_uint(world->rng) % n_players;
              /* pass reactants */
              if (rnd_num < 2) continue;

              if(rx_players[rnd_num] == NULL) continue;
              if((rx_players[rnd_num]->flags & NOT_FREE) == 0) continue;

              if(product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET)
              {
                if(replace_p1)
                {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                  product_grid[rnd_num] = ((struct grid_molecule *)reacA)->grid;
                  product_grid_idx[rnd_num] = ((struct grid_molecule *)reacA)->grid_index;
                  replace_p1 = 0;
                  count++;
                  continue;
                }

                if(replace_p2)
                {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                  product_grid[rnd_num] = ((struct grid_molecule *)reacB)->grid;
                  product_grid_idx[rnd_num] = ((struct grid_molecule *)reacB)->grid_index;
                  replace_p2 = 0;
                  count++;
                  continue;
                }
              }
            } /* end while */
          }
        }
      }
    }

    /* here we will find placement for the case of the reaction
       of type "vol_mol + w -> surf_mol + ...[rate] " */
    if((grid_reactant == NULL) && (w != NULL) && (num_surface_products >= 1))
    {
       assert(! IS_GRID_MOL(reacA));
       assert(rxn_uv_idx != -1);

       while(true)
       {
          rnd_num = rng_uint(world->rng) % (n_players);
          /* skip the reactant in the players list */
          if(rnd_num == 0) continue;
          /* skip the wall in the players list */
          if(rnd_num == 1) continue;

          if(rx_players[rnd_num] == NULL) continue;
          if((rx_players[rnd_num]->flags & NOT_FREE) == 0) continue;

          if(product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET)
          {
              product_flag[rnd_num] = PRODUCT_FLAG_USE_UV_LOC;
              product_grid[rnd_num] = w->grid;
              product_grid_idx[rnd_num] = rxn_uv_idx;
              break;
          }
       }
    }

   /* we will implement special placement policy for reaction
        of  type of A->B+C[rate] */
    if(is_unimol && (grid_reactant != NULL) && (num_surface_products == 2))
    {
       reac_idx = grid_reactant->grid_index;
       reac_grid = grid_reactant->grid;
    }


    /* all other products are placed on one of the randomly chosen vacant
       tiles */
    int do_it_once = 0; /* flag */
    int num_attempts = 0;
    for (int n_product = rx->n_reactants; n_product < n_players; ++ n_product)
    {
      /* If the product is a volume product, no placement is required. */
      if (rx_players[n_product]->flags & ON_GRID)
      {
        if(product_flag[n_product] != PRODUCT_FLAG_NOT_SET)
        {
          continue;
        }

        /* can't place products - reaction blocked */
        if(num_vacant_tiles == 0)
        {
          if(tile_nbr_head != NULL)
          {
            delete_tile_neighbor_list(tile_nbr_head);
          }
          if(tile_vacant_nbr_head != NULL)
          {
            delete_tile_neighbor_list(tile_vacant_nbr_head);
          }
          return RX_BLOCKED;
        }

        num_attempts = 0;

        while(true)
        {
          if(num_attempts > SURFACE_DIFFUSION_RETRIES)
          {
            if(tile_nbr_head != NULL)
            {
              delete_tile_neighbor_list(tile_nbr_head);
            }
            if(tile_vacant_nbr_head != NULL)
            {
              delete_tile_neighbor_list(tile_vacant_nbr_head);
            }
            return RX_BLOCKED;
          }

          /* randomly pick a tile from the list */
          rnd_num = rng_uint(world->rng) % num_vacant_tiles;
          tile_idx = -1;
          tile_grid = NULL;

          if(get_tile_neighbor_from_list_of_vacant_neighbors(tile_vacant_nbr_head, rnd_num, &tile_grid, &tile_idx) == 0)
          {
            if(tile_nbr_head != NULL)
            {
              delete_tile_neighbor_list(tile_nbr_head);
            }
            if(tile_vacant_nbr_head != NULL)
            {
              delete_tile_neighbor_list(tile_vacant_nbr_head);
            }
            return RX_BLOCKED;
          }
          if(tile_idx < 0)
          {
            continue; /* this tile was probed already */
          }

          assert(tile_grid != NULL);

          if (!product_tile_can_be_reached(tile_grid->surface,
                                           rlp_head_wall_1,
                                           rlp_head_wall_2,
                                           rlp_head_obj_1,
                                           rlp_head_obj_2,
                                           grid_bitmask,
                                           is_unimol))
          {
            uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
            num_attempts++;
            continue;
          }

          product_grid[n_product]     = tile_grid;
          product_grid_idx[n_product] = tile_idx;
          product_flag[n_product]     = PRODUCT_FLAG_USE_RANDOM;
          if(!do_it_once && is_unimol && (grid_reactant != NULL) &&
              (num_surface_products == 2))
          {
            /*remember this tile (used for the reaction A->B+C[rate]) */
            mol_idx = tile_idx;
            mol_grid = tile_grid;
            do_it_once = 1;
          }
          break;
        } /* end while */
      }
    }

  } /* end if(is_orientable) */

  /* Determine the location of the reaction for count purposes. */
  struct vector3 count_pos_xyz;
  if (hitpt != NULL) count_pos_xyz = *hitpt;
  else if (grid_reactant) uv2xyz(& grid_reactant->s_pos, grid_reactant->grid->surface, & count_pos_xyz);
  else count_pos_xyz = ((struct volume_molecule *) reacA)->pos;

  /* Create and place each product. */
  struct vector3 mol_pos_tmp;
  struct subvolume *product_subvol = NULL;
  for (int n_product = rx->n_reactants; n_product < n_players; ++ n_product)
  {
    struct abstract_molecule *this_product = NULL;
    struct species * const product_species = rx_players[n_product];

    /* If the product is a grid molecule, place it on the grid. */
    if (product_species->flags & ON_GRID)
    {
        struct vector2 prod_uv_pos;

        /* Pick an appropriate position for the new molecule. */
        if (world->randomize_gmol_pos)
        {
          switch (product_flag[n_product])
          {
            case PRODUCT_FLAG_USE_REACA_UV:
              if(is_unimol && (num_surface_products == 2) && (grid_reactant != NULL))
              {
                 if(mol_grid == NULL) mcell_internal_error("Error in surface product placement for the unimolecular reaction.");
                 find_closest_position(product_grid[n_product], product_grid_idx[n_product], mol_grid, mol_idx, &prod_uv_pos);
              }else{
                 prod_uv_pos = ((struct grid_molecule*) reacA)->s_pos;
              }
              break;

            case PRODUCT_FLAG_USE_REACB_UV:
              assert(reacB != NULL);
              prod_uv_pos = ((struct grid_molecule*) reacB)->s_pos;
              break;

            case PRODUCT_FLAG_USE_UV_LOC:
              prod_uv_pos = rxn_uv_pos;
              break;

            case PRODUCT_FLAG_USE_RANDOM:
              if(is_unimol && replace_p1 && (num_surface_products == 2))
              {
                 find_closest_position(product_grid[n_product], product_grid_idx[n_product], reac_grid, reac_idx, &prod_uv_pos);
              }else{
                 grid2uv_random(product_grid[n_product], product_grid_idx[n_product], &prod_uv_pos);
              }
              break;

            default:
              UNHANDLED_CASE(product_flag[n_product]);
              break;
          }
        }
        else
          grid2uv(product_grid[n_product], product_grid_idx[n_product], & prod_uv_pos);

        this_product = (struct abstract_molecule *)
              place_grid_product(product_species,
                                 product_grid[n_product],
                                 product_grid_idx[n_product],
                                 & prod_uv_pos,
                                 product_orient[n_product],
                                 t);
    }

    /* else place the molecule in space. */
    else
    {
        /* Unless this is a unimolecular reaction, we will have a hitpt. */
        if (! hitpt)
        {
          /* If this is a unimolecular surface rxn... */
          if (reacA->properties->flags & ON_GRID)
          {
            uv2xyz(& ((struct grid_molecule *) reacA)->s_pos,
                   ((struct grid_molecule *) reacA)->grid->surface,
                   & mol_pos_tmp);
            product_subvol = find_subvolume(&mol_pos_tmp, last_subvol);
          }

          /* ... else a unimolecular volume rxn. */
          else
          {
            mol_pos_tmp = ((struct volume_molecule *) reacA)->pos;
            product_subvol = ((struct volume_molecule *) reacA)->subvol;
          }
          hitpt = & mol_pos_tmp;
        }
        else if (product_subvol == NULL)
          product_subvol = find_subvolume(hitpt, last_subvol);

        this_product = (struct abstract_molecule *) place_volume_product(product_species,
                                                                         grid_reactant,
                                                                         w,
                                                                         product_subvol,
                                                                         hitpt,
                                                                         product_orient[n_product],
                                                                         t);

      if (((struct volume_molecule *) this_product)->index < DISSOCIATION_MAX)
        update_dissociation_index = true;
    }

    /* Update molecule counts */
    ++ product_species->population;
    if (product_species->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
      count_region_from_scratch(this_product, NULL, 1, NULL, NULL, t);

    /* preserve molecule id if rxn is unimolecular with one product */
    if (is_unimol && (n_players == 2))
    {
      this_product->id = reacA->id;
      world->current_mol_id--;  /* give back id we used */
      continue;
    }
    /* preserve molecule id if rxn is surface rxn with one product */
    if ((n_players == 3) && product_type[1] == PLAYER_WALL)
    {
      this_product->id = reacA->id;
      world->current_mol_id--;  /* give back id we used */
      continue;
    }
    /* preserve molecule id if product is same species as reactant */
/* Disable this code to see if we don't really need it: */
#if 0
    if ((!matched_reacA) && (this_product->properties == reacA->properties))
    {
      matched_reacA = 1;
      this_product->id = reacA->id;
      world->current_mol_id--;  /* give back id we used */
      continue;
    }
    if ((reacB) && (!matched_reacB) && (this_product->properties == reacB->properties))
    {
      matched_reacB = 1;
      this_product->id = reacB->id;
      world->current_mol_id--;  /* give back id we used */
    }
#endif /* if 0 */

  }

  /* If necessary, update the dissociation index. */
  if (update_dissociation_index)
  {
    if (-- world->dissociation_index < DISSOCIATION_MIN)
      world->dissociation_index = DISSOCIATION_MAX;
  }

  /* Handle events triggered off of named reactions */
  if (rx->info[path].pathname!=NULL)
  {
    /* No flags for reactions so we have to check regions if we have waypoints! Fix to be more efficient for WORLD-only counts? */
    if (world->place_waypoints_flag)
      count_region_from_scratch(NULL, rx->info[path].pathname, 1, &count_pos_xyz, w, t);

    /* Other magical stuff.  For now, can only trigger releases. */
    if (rx->info[path].pathname->magic!=NULL)
    {
      if (reaction_wizardry(rx->info[path].pathname->magic, w, &count_pos_xyz, t))
        mcell_allocfailed("Failed to complete reaction triggered release after a '%s' reaction.",
                          rx->info[path].pathname->sym->name);
    }
  }


  /* recover memory */
  delete_tile_neighbor_list(tile_nbr_head);
  delete_tile_neighbor_list(tile_vacant_nbr_head);
  delete_region_list(rlp_head_wall_1);
  delete_region_list(rlp_head_wall_2);
  delete_region_list(rlp_head_obj_1);
  delete_region_list(rlp_head_obj_2);

  return cross_wall ? RX_FLIP : RX_A_OK;
}



/*************************************************************************
outcome_products_trimol_reaction_random:
   In: first wall in the reaction
       hit point (if any)
       time of the reaction
       reaction
       path of the reaction
       first reactant (moving molecule)
       second reactant
       third reactant
       orientation of the first reactant
       orientation of the second reactant
       orientation of the third reactant
Note: This function replaces surface reactants (if needed) by the surface
       products picked in the random order from the list of products.
       It also places surface products on the randomly selected tiles
       from the list of neighbors.
Note: Policy on surface products placement is described in the document
      "policy_surf_products_placement.doc" (see "src/docs").

************************************************************************/
static int outcome_products_trimol_reaction_random(struct wall *w,
                            struct vector3 *hitpt,
                            double t,
                            struct rxn *rx,
                            int path,
                            struct abstract_molecule *reacA,
                            struct abstract_molecule *reacB,
                            struct abstract_molecule *reacC,
                            short orientA,
                            short orientB,
                            short orientC)
{
  bool update_dissociation_index = false;         /* Do we need to advance the dissociation index? */
  bool cross_wall = false;                        /* Did the moving molecule cross the plane? */
  struct subvolume *last_subvol = NULL;           /* Last subvolume (guess used to speed sv finding) */

  int const i0 = rx->product_idx[path];           /* index of the first player for the pathway */
  int const iN = rx->product_idx[path+1];         /* index of the first player for the next pathway */
  assert(iN > i0);
  struct species **rx_players = rx->players + i0; /* Players array from the reaction. */

  int const n_players = iN - i0;                  /* number of reaction players */
  struct abstract_molecule *product[n_players];   /* array of products */
  char product_type[n_players];                   /* array that decodes the type of each product */
  short product_orient[n_players];                /* array of orientations for each product */
  struct surface_grid *product_grid[n_players];   /* array of surface_grids for products */
  int product_grid_idx[n_players];                /* array of grid indices for products */
  byte product_flag[n_players];                   /* array of placement flags for products */

  struct tile_neighbor *tile_nbr_head = NULL;  /* list of neighbor tiles */
  struct tile_neighbor *tile_nbr;  /* iterator */
  /* head of the linked list of vacant neighbor tiles */
  struct tile_neighbor *tile_vacant_nbr_head = NULL;
  struct surface_grid *tile_grid;  /* surface grid the tile belongs to */
  int tile_idx;    /* index of the tile on the grid */
  unsigned int rnd_num;  /* random number */
  int num_vacant_tiles = 0; /* number of vacant tiles */
  int num_surface_products = 0; /* not counting reactants */
  int num_surface_static_products = 0;/* # of products with (D_2D == 0) not counting reactants */
  int num_surface_static_reactants = 0;/* # of reactants with (D_2D == 0) */
  /* total number of surface reactants */
  int num_surface_reactants = 0;
  /* number of surface reactants that are not replaced in the reaction */
  int num_surface_reactants_to_stay = 0;
  int list_length; /* length of the linked list tile_nbr_head */
  /* flags */
  int replace_p1 = 0, replace_p2 = 0, replace_p3 = 0, only_one_to_replace = 0, two_to_replace = 0;
  int find_neighbor_tiles_flag = 0;
  struct wall *w_1, *w_2, *w_3;

  /* flag that indicates that all reactants lie inside their
     respective restrictive regions */
  int all_inside_restricted_boundary = 0;
  /* flag that indicates that all reactants lie outside
     their respective restrictive regions */
  int all_outside_restricted_boundary = 0;

  /* flag that indicates that  reactants "grid_1" and "grid_2"
     lies inside and reactant "grid_3" lies outside of their
     respective restrictive regions */
  int grid_1_inside_grid_2_inside_grid_3_outside = 0;
  /* flag that indicates that  reactants "grid_1" and "grid_3"
     lies inside and reactant "grid_2" lies outside of their
     respective restrictive regions */
  int grid_1_inside_grid_2_outside_grid_3_inside = 0;
  /* flag that indicates that  reactant "grid_1" lies outside
     and reactants "grid_2" and "grid_3" lie inside of their
     respective restrictive regions */
  int grid_1_outside_grid_2_inside_grid_3_inside = 0;
  /* flag that indicates that  reactant "grid_1" lies inside
     and reactants "grid_2" and "grid_3" lie outside of their
     respective restrictive regions */
  int grid_1_inside_grid_2_outside_grid_3_outside = 0;
  /* flag that indicates that  reactants "grid_1" and "grid_3" lie outside
     and reactant "grid_2" lies inside of their
     respective restrictive regions */
  int grid_1_outside_grid_2_inside_grid_3_outside = 0;
  /* flag that indicates that  reactants "grid_1" and "grid_2"
     lie outside and reactant "grid_3" lies inside of their
     respective restrictive regions */
  int grid_1_outside_grid_2_outside_grid_3_inside = 0;

  /* flag that indicates that  reactants "grid_1" and "grid_2"
     lie inside their respective restrictive regions
     and reactant "grid_3" has no restrictive region border properties */
  int only_grid_1_grid_2_inside = 0;
  /* flag that indicates that  reactant "grid_1" lies inside
     and "grid_2" lies outside their respective restrictive regions
     and reactant "grid_3" has no restrictive region border properties */
  int only_grid_1_inside_grid_2_outside = 0;
  /* flag that indicates that  reactant "grid_1" lies outside
     and "grid_2" lies inside their respective restrictive regions
     and reactant "grid_3" has no restrictive region border properties */
  int only_grid_1_outside_grid_2_inside = 0;
  /* flag that indicates that  reactants "grid_1" and "grid_2"
     lie outside their respective restrictive regions
     and reactant "grid_3" has no restrictive region border properties */
  int only_grid_1_grid_2_outside = 0;

  /* flag that indicates that  reactants "grid_1" and "grid_3"
     lie inside their respective restrictive regions
     and reactant "grid_2" has no restrictive region border properties */
  int only_grid_1_grid_3_inside = 0;
  /* flag that indicates that  reactant "grid_1" lies inside
     and "grid_3" lies outside their respective restrictive regions
     and reactant "grid_2" has no restrictive region border properties */
  int only_grid_1_inside_grid_3_outside = 0;
  /* flag that indicates that  reactant "grid_1" lies outside
     and "grid_3" lies inside their respective restrictive regions
     and reactant "grid_2" has no restrictive region border properties */
  int only_grid_1_outside_grid_3_inside = 0;
  /* flag that indicates that  reactants "grid_1" and "grid_3"
     lie outside their respective restrictive regions
     and reactant "grid_2" has no restrictive region border properties */
  int only_grid_1_grid_3_outside = 0;

  /* flag that indicates that  reactants "grid_2" and "grid_3"
     lie inside their respective restrictive regions
     and reactant "grid_1" has no restrictive region border properties */
  int only_grid_2_grid_3_inside = 0;
  /* flag that indicates that  reactant "grid_2" lies inside
     and "grid_3" lies outside their respective restrictive regions
     and reactant "grid_1" has no restrictive region border properties */
  int only_grid_2_inside_grid_3_outside = 0;
  /* flag that indicates that  reactant "grid_2" lies outside
     and "grid_3" lies inside their respective restrictive regions
     and reactant "grid_1" has no restrictive region border properties */
  int only_grid_2_outside_grid_3_inside = 0;
  /* flag that indicates that  reactants "grid_2" and "grid_3"
     lie outside their respective restrictive regions
     and reactant "grid_1" has no restrictive region border properties */
  int only_grid_2_grid_3_outside = 0;

  /* flag that indicates that only reactant "grid_1" has
     restrictive regions on the object and it lies
     inside it's restrictive region.  */
  int only_grid_1_inside = 0;
  /* flag that indicates that only reactant "grid_1" has
     restrictive regions on the object and it lies
     outside it's restrictive region.  */
  int only_grid_1_outside = 0;
  /* flag that indicates that only reactant "grid_2" has
     restrictive regions on the object and it lies
     inside it's restrictive region.  */
  int only_grid_2_inside = 0;
  /* flag that indicates that only reactant "grid_2" has
     restrictive regions on the object and it lies
     outside it's restrictive region.  */
  int only_grid_2_outside = 0;
  /* flag that indicates that only reactant "grid_3" has
     restrictive regions on the object and it lies
     inside it's restrictive region.  */
  int only_grid_3_inside = 0;
  /* flag that indicates that only reactant "grid_3" has
     restrictive regions on the object and it lies
     outside it's restrictive region.  */
  int only_grid_3_outside = 0;

  /* list of the restricted regions for the reactants by wall */
  struct region_list *rlp_head_wall_1 = NULL, *rlp_head_wall_2 = NULL,
  *rlp_head_wall_3 = NULL;
  /* list of the restricted regions for the reactants by object */
  struct region_list *rlp_head_obj_1 = NULL, *rlp_head_obj_2 = NULL,
  *rlp_head_obj_3 = NULL;

  struct vector2 rxn_uv_pos; /* position of the reaction */
  int rxn_uv_idx = -1;  /* tile index of the reaction place */

  struct abstract_molecule * tmp_mol;
  short tmp_orient;

  if (rx->is_complex)
  {
      mcell_internal_error("Function 'outcome_products_trimol_reaction_random()' does not support macromolecular reaction [%s].", rx->sym->name);
  }

  if((reacA == NULL) || (reacB == NULL) || (reacC == NULL))
  {
      mcell_internal_error("One of the reactants in 'outcome_products_trimol_reaction_random()' is NULL.");
  }

  /* Clear the initial product info. */
  for (int i=0; i<n_players; ++i)
  {
    product[i] = NULL;
    product_type[i] = PLAYER_NONE;
    product_orient[i] = 0;
    product_grid[i] = NULL;
    product_grid_idx[i] = -1;
    product_flag[i] = PRODUCT_FLAG_NOT_SET;
  }

  /* Flag indicating that a surface is somehow involved with this reaction. */
  struct grid_molecule * const grid_1 = IS_GRID_MOL(reacA) ? (struct grid_molecule *) reacA : NULL;
  struct grid_molecule * const grid_2 = IS_GRID_MOL(reacB) ? (struct grid_molecule *) reacB : NULL;
  struct grid_molecule * const grid_3 = IS_GRID_MOL(reacC) ? (struct grid_molecule *) reacC : NULL;
  struct grid_molecule * grid_reactant = NULL;
  if(grid_1 != NULL)
  {
    grid_reactant = grid_1;
  }else if(grid_2 != NULL){
    grid_reactant = grid_2;
  }else if(grid_3 != NULL){
    grid_reactant = grid_3;
  }

  bool const is_orientable = (w != NULL)  ||  (grid_reactant != NULL);

  /* Where are reactants relative to the restrictive region border? */
  if((grid_1 != NULL) && (grid_2 != NULL) && (grid_3 != NULL))
  {
     /* trimol_reaction */
     if((grid_1->properties->flags & CAN_REGION_BORDER) &&
        (grid_2->properties->flags & CAN_REGION_BORDER) &&
        (grid_3->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1) &&
        are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2) &&
        are_restricted_regions_for_species_on_object(grid_3->grid->surface->parent_object, grid_3))
     {
       w_1 = grid_1->grid->surface;
       w_2 = grid_2->grid->surface;
       w_3 = grid_3->grid->surface;
       rlp_head_wall_1 = find_restricted_regions_by_wall(w_1, grid_1);
       rlp_head_wall_2 = find_restricted_regions_by_wall(w_2, grid_2);
       rlp_head_wall_3 = find_restricted_regions_by_wall(w_3, grid_3);

       if((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL)  &&
           (rlp_head_wall_3 != NULL))
       {
          /* all reactants are inside their respective
             restricted regions */
             all_inside_restricted_boundary = 1;

       }else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL) &&
               (rlp_head_wall_3 == NULL)) {
          /* all reactants are outside their respective
             restricted regions */
             rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
             rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
             rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);

             all_outside_restricted_boundary = 1;
       }else if((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL)
             && (rlp_head_wall_3 == NULL))
       {
             rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);
             grid_1_inside_grid_2_inside_grid_3_outside = 1;
       }else if((rlp_head_wall_1 != NULL) && (rlp_head_wall_3 != NULL)
             && (rlp_head_wall_2 == NULL))
       {
             rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
             grid_1_inside_grid_2_outside_grid_3_inside = 1;
       }else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL)
             && (rlp_head_wall_3 == NULL))
       {
             rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
             rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);
             grid_1_outside_grid_2_inside_grid_3_outside = 1;
       }else if((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 == NULL)
             && (rlp_head_wall_3 == NULL))
       {
             rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
             rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);
             grid_1_inside_grid_2_outside_grid_3_outside = 1;
       }else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL)
             && (rlp_head_wall_3 != NULL))
       {
             rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
             grid_1_outside_grid_2_inside_grid_3_inside = 1;
       }else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL)
             && (rlp_head_wall_3 != NULL))
       {
             rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
             rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
             grid_1_outside_grid_2_outside_grid_3_inside = 1;
       }
     }
     else if((grid_1->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1) &&
        (grid_2->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2) &&
        (!(grid_3->properties->flags & CAN_REGION_BORDER) ||
        !are_restricted_regions_for_species_on_object(grid_3->grid->surface->parent_object, grid_3)))
     {
        /* only reactants "grid_1" and "grid_2" have restrictive
           region border property */
       w_1 = grid_1->grid->surface;
       rlp_head_wall_1 = find_restricted_regions_by_wall(w_1, grid_1);
       w_2 = grid_2->grid->surface;
       rlp_head_wall_2 = find_restricted_regions_by_wall(w_2, grid_2);
       if((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL)) {
           only_grid_1_grid_2_inside = 1;
       }else if((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 == NULL)){
           only_grid_1_inside_grid_2_outside = 1;
           rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
       }else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL)){
           only_grid_1_outside_grid_2_inside = 1;
           rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
       }else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL)){
           only_grid_1_grid_2_outside = 1;
           rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
           rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
       }
     }
     else if((grid_1->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1) &&
        (!(grid_2->properties->flags & CAN_REGION_BORDER) ||
        !are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2)) &&
        (grid_3->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(grid_3->grid->surface->parent_object, grid_3))
     {
        /* only reactants "grid_1" and "grid_3" have restrictive
           region border property */
       w_1 = grid_1->grid->surface;
       rlp_head_wall_1 = find_restricted_regions_by_wall(w_1, grid_1);
       w_3 = grid_3->grid->surface;
       rlp_head_wall_3 = find_restricted_regions_by_wall(w_3, grid_3);
       if((rlp_head_wall_1 != NULL) && (rlp_head_wall_3 != NULL)) {
           only_grid_1_grid_3_inside = 1;
       }else if((rlp_head_wall_1 != NULL) && (rlp_head_wall_3 == NULL)){
           only_grid_1_inside_grid_3_outside = 1;
           rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);
       }else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_3 != NULL)){
           only_grid_1_outside_grid_3_inside = 1;
           rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
       }else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_3 == NULL)){
           only_grid_1_grid_3_outside = 1;
           rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
           rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);
       }
     }
     else if((!(grid_1->properties->flags & CAN_REGION_BORDER) ||
        (!are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1))) &&
        (grid_2->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2) &&
        (grid_3->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(grid_3->grid->surface->parent_object, grid_3))
     {
        /* only reactants "grid_2" and "grid_3" have restrictive
           region border property */
       w_2 = grid_2->grid->surface;
       rlp_head_wall_2 = find_restricted_regions_by_wall(w_2, grid_2);
       w_3 = grid_3->grid->surface;
       rlp_head_wall_3 = find_restricted_regions_by_wall(w_3, grid_3);
       if((rlp_head_wall_2 != NULL) && (rlp_head_wall_3 != NULL)) {
           only_grid_2_grid_3_inside = 1;
       }else if((rlp_head_wall_2 != NULL) && (rlp_head_wall_3 == NULL)){
           only_grid_2_inside_grid_3_outside = 1;
           rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);
       }else if((rlp_head_wall_2 == NULL) && (rlp_head_wall_3 != NULL)){
           only_grid_2_outside_grid_3_inside = 1;
           rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
       }else if((rlp_head_wall_2 == NULL) && (rlp_head_wall_3 == NULL)){
           only_grid_2_grid_3_outside = 1;
           rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
           rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);
       }
     }
     else if((grid_1->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1) &&
        (!(grid_2->properties->flags & CAN_REGION_BORDER) ||
        !are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2)) &&
        (!(grid_3->properties->flags & CAN_REGION_BORDER) ||
        !are_restricted_regions_for_species_on_object(grid_3->grid->surface->parent_object, grid_3)))
     {
        /* only reactant "grid_1" has restrictive region border property */
       w_1 = grid_1->grid->surface;
       rlp_head_wall_1 = find_restricted_regions_by_wall(w_1, grid_1);
       if(rlp_head_wall_1 != NULL) {
           only_grid_1_inside = 1;
       }else{
           rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
           only_grid_1_outside = 1;
       }
     }
     else if ((grid_2->properties->flags & CAN_REGION_BORDER) &&
             are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2)  &&
            (!(grid_1->properties->flags & CAN_REGION_BORDER) ||
             !are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1))  &&
            (!(grid_3->properties->flags & CAN_REGION_BORDER) ||
             !are_restricted_regions_for_species_on_object(grid_3->grid->surface->parent_object, grid_3)))
     {
        /* only reactant "grid_2" has restrictive region border property */
       w_2 = grid_2->grid->surface;
       rlp_head_wall_2 = find_restricted_regions_by_wall(w_2, grid_2);
       if(rlp_head_wall_2 != NULL) only_grid_2_inside = 1;
       else{
          rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
          only_grid_2_outside = 1;
       }
     }
     else if ((grid_3->properties->flags & CAN_REGION_BORDER) &&
             are_restricted_regions_for_species_on_object(grid_3->grid->surface->parent_object, grid_3)  &&
            (!(grid_1->properties->flags & CAN_REGION_BORDER) ||
             !are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1))  &&
            (!(grid_2->properties->flags & CAN_REGION_BORDER) ||
             !are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2)))
     {
        /* only reactant "grid_3" has restrictive region border property */
       w_3 = grid_3->grid->surface;
       rlp_head_wall_3 = find_restricted_regions_by_wall(w_3, grid_3);
       if(rlp_head_wall_3 != NULL) only_grid_3_inside = 1;
       else{
          rlp_head_obj_3 = find_restricted_regions_by_object(w_3->parent_object, grid_3);
          only_grid_3_outside = 1;
       }
     }
  }


  /* reacA is the molecule which initiated the reaction. */
  struct abstract_molecule * const initiator = reacA;
  short const initiatorOrient = orientA;

  /* make sure that reacA corresponds to rx->players[0],
     reacB - to rx->players[1], and reacC - to rx->players[2]
     in trimolecular reaction */
  if(reacA->properties == rx->players[0])
  {
       if(reacB->properties == rx->players[2] && reacB->properties != rx->players[1]){
           /* switch B and C */
         tmp_mol = reacB;
         reacB = reacC;
         reacC = tmp_mol;

         tmp_orient = orientB;
         orientB = orientC;
         orientC = tmp_orient;
       }
  }else if(reacA->properties == rx->players[1]){

       if(reacB->properties == rx->players[0] && reacB->properties != rx->players[1]){
          /* switch reacA and reacB */
         tmp_mol = reacB;
         reacB = reacA;
         reacA = tmp_mol;

         tmp_orient = orientB;
         orientB = orientA;
         orientA = tmp_orient;
       }else if(reacC->properties == rx->players[0]){
          /* switch reacA and reacC */
         tmp_mol = reacA;
         reacA = reacC;
         reacC = tmp_mol;

         tmp_orient = orientA;
         orientA = orientC;
         orientC = tmp_orient;
          /* now switch reacC and reacB  */
         tmp_mol = reacB;
         reacB = reacC;
         reacC = tmp_mol;

         tmp_orient = orientB;
         orientB = orientC;
         orientC = tmp_orient;
       }
  }else if(reacA->properties == rx->players[2]){
        if(reacB->properties == rx->players[0])
        {
          /* switch reacA and reacB */
         tmp_mol = reacB;
         reacB = reacA;
         reacA = tmp_mol;

         tmp_orient = orientB;
         orientB = orientA;
         orientA = tmp_orient;

          /* switch reacB and reacC */
         tmp_mol = reacB;
         reacB = reacC;
         reacC = tmp_mol;

         tmp_orient = orientB;
         orientB = orientC;
         orientC = tmp_orient;

        }else if ((reacC->properties == rx->players[0]) &&
           (reacC->properties != rx->players[2])){
          /* switch reacA and reacC */
         tmp_mol = reacA;
         reacA = reacC;
         reacC = tmp_mol;

         tmp_orient = orientA;
         orientA = orientC;
         orientC = tmp_orient;

       }
  }

  /* Add the reactants (incl. any wall) to the list of players. */
  add_players_to_list(rx, reacA, reacB, reacC, product, product_type);

  /* Determine whether any of the reactants can be replaced by a product. */
  if (product_type[0] == PLAYER_GRID_MOL)
  {
    num_surface_reactants++;
    if(rx_players[0] == NULL) replace_p1 = 1;
    else num_surface_reactants_to_stay++;
  }
  if (product_type[1] == PLAYER_GRID_MOL)
  {
    num_surface_reactants++;
    if(rx_players[1] == NULL) replace_p2 = 1;
    else num_surface_reactants_to_stay++;
  }
  if (product_type[2] == PLAYER_GRID_MOL)
  {
    num_surface_reactants++;
    if(rx_players[2] == NULL) replace_p3 = 1;
    else num_surface_reactants_to_stay++;
  }


  if(replace_p1 && (!replace_p2) && (!replace_p3))
  {
      only_one_to_replace = 1;
  }else if((!replace_p1) && replace_p2 && (!replace_p3)){
      only_one_to_replace = 1;
  }else if((!replace_p1) && (!replace_p2) && replace_p3){
      only_one_to_replace = 1;
  }
  if(replace_p1 && (replace_p2) && (!replace_p3))
  {
      two_to_replace = 1;
  }else if(replace_p1 && (!replace_p2) && replace_p3){
      two_to_replace = 1;
  }else if((!replace_p1) && replace_p2 && replace_p3){
      two_to_replace = 1;
  }

   /* find out number of surface products */
  for(int n_product = rx->n_reactants; n_product < n_players; ++ n_product)
  {
     if(rx_players[n_product] == NULL) continue;
     if(rx_players[n_product]->flags & ON_GRID) {
        num_surface_products++;
        if(rx_players[n_product]->D == 0) num_surface_static_products++;
     }
  }

  if(num_surface_reactants >= 2) find_neighbor_tiles_flag = 1;
  if((num_surface_reactants == 1) && (num_surface_products > 1)) find_neighbor_tiles_flag = 1;


  /* Determine the point of reaction on the surface. */
  if(is_orientable)
  {
     if (grid_reactant) rxn_uv_pos = grid_reactant->s_pos;
     else
     {
       xyz2uv(hitpt, w, &rxn_uv_pos);
     }

     if((w == NULL) && (grid_reactant != NULL)) w = grid_reactant->grid->surface;

     assert(w != NULL);
     if(w->grid == NULL)
     {
           /* reacA must be a volume molecule, or this wall would have a grid already. */
        assert(! IS_GRID_MOL(reacA));

        if (create_grid(w, ((struct volume_molecule *) reacA)->subvol))
               mcell_allocfailed("Failed to create a grid for a wall.");
     }

     if(find_neighbor_tiles_flag)
     {
       /* create list of neighbor tiles around rxn_uv_pos */
       rxn_uv_idx = uv2grid(&rxn_uv_pos, w->grid);

       if(grid_reactant != NULL)
       {
         find_neighbor_tiles(grid_reactant, grid_reactant->grid, grid_reactant->grid_index, 1, 0, &tile_nbr_head, &list_length);
       }else{
         find_neighbor_tiles(grid_reactant, w->grid, rxn_uv_idx, 1, 0, &tile_nbr_head, &list_length);
       }


       /* Create list of vacant tiles */
       for(tile_nbr = tile_nbr_head; tile_nbr != NULL; tile_nbr = tile_nbr->next)
       {
          if(tile_nbr->grid->mol[tile_nbr->idx] == NULL)
          {
             num_vacant_tiles++;
             push_tile_neighbor_to_list(&tile_vacant_nbr_head, tile_nbr->grid, tile_nbr->idx);
          }
       }
     }

  }

  /* find out number of static surface reactants */
  if((grid_1 != NULL) && (grid_1->properties->D == 0)) num_surface_static_reactants++;
  if((grid_2 != NULL) && (grid_2->properties->D == 0)) num_surface_static_reactants++;
  if((grid_3 != NULL) && (grid_3->properties->D == 0)) num_surface_static_reactants++;



  /* If the reaction involves a surface, make sure there is room for each product. */
  if (is_orientable)
  {
    /* Can this reaction happen at all? */
    if(replace_p1 && replace_p2 && replace_p3)
    {
       if(num_surface_products > num_vacant_tiles + 3)
       {
          if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
          if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
          return RX_BLOCKED;
       }
    }else if(two_to_replace)
    {
       if(num_surface_products > num_vacant_tiles + 2)
       {
          if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
          if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
          return RX_BLOCKED;
       }
    }else if(only_one_to_replace) {
       if(num_surface_products > num_vacant_tiles + 1)
       {
          if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
          if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
          return RX_BLOCKED;
       }
    }else{
       /* none of the reactants is replaced */
       if(num_surface_products > num_vacant_tiles)
       {
          if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
          if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
          return RX_BLOCKED;
       }
    }

    /* set the orientations of the products. */
    for (int n_product = 0; n_product < n_players; ++ n_product)
    {
      /* Skip NULL reactants. */
      if (rx_players[n_product] == NULL)
        continue;

      int this_geometry = rx->geometries[i0 + n_product];

      /* Geometry of 0 means "random orientation" */
      if (this_geometry == 0)
        product_orient[n_product] = (rng_uint(world->rng) & 1) ? 1 : -1;
      else
      {
        /* Geometry < 0 means inverted orientation */
        if (this_geometry < 0)
        {
          this_geometry = -this_geometry;
          if (this_geometry > (int) rx->n_reactants)
            product_orient[n_product] = - product_orient[this_geometry - rx->n_reactants - 1];
          else if (this_geometry == 1)
            product_orient[n_product] = - orientA;
          else if ((this_geometry == 2) && (reacB != NULL))
            product_orient[n_product] = - orientB;
          else if ((this_geometry == 3) && (reacC != NULL))
            product_orient[n_product] = - orientC;
          else
            product_orient[n_product] = -1;
        }

        /* Geometry > 0 means "positive" orientation. */
        else
        {
          if (this_geometry > (int) rx->n_reactants)
            product_orient[n_product] = product_orient[this_geometry - rx->n_reactants - 1];
          else if (this_geometry == 1)
            product_orient[n_product] = orientA;
          else if ((this_geometry == 2) && (reacB != NULL))
            product_orient[n_product] = orientB;
          else if ((this_geometry == 3) && (reacC != NULL))
            product_orient[n_product] = orientC;
          else
            product_orient[n_product] = 1;
        }
      }

      /* If this is a reactant... */
      if (n_product < (int) rx->n_reactants)
      {
        /* If this is a grid molecule, we need to set its orientation. */
        if (rx_players[n_product]->flags & ON_GRID)
        {
          assert(IS_GRID_MOL(product[n_product]));
          struct grid_molecule *gm = (struct grid_molecule *) product[n_product];

          /* If the new orientation doesn't match the old, we've got some work to do. */
          if (gm->orient != product_orient[n_product])
          {

            /* We're about to update the molecule's orientation, so we will
             * first remove it from the counts in case we have any
             * orientation-sensitive counts.  Then, we will update the
             * orientation.  Finally, we will add the molecule back into the
             * counts in its new orientation.
             */

            /* Remove molecule from counts in old orientation, if mol is counted. */
            if (product[n_product]->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
              count_region_from_scratch(product[n_product],     /* molecule */
                                        NULL,                   /* rxn pathway */
                                        -1,                     /* remove count */
                                        NULL,                   /* Location at which to count */
                                        w,                      /* Wall on which this happened */
                                        t);                     /* Time of occurrence */

            /* Force check for the unimolecular reactions
               after changing orientation.
               There are two possible cases to be covered here:
               1) when (gm->t2) was previously set to FOREVER
               2) there may be two or more unimolecular
                  reactions involving surface class that have
                  different kinetics.
             */
            if(((gm->flags & ACT_REACT) != 0) && ((gm->properties->flags & CAN_GRIDWALL) != 0))  gm->t2 = 0;

            /* Set the molecule's orientation. */
            gm->orient = product_orient[n_product];

            /* Add molecule back to counts in new orientation, if mol is counted. */
            if (product[n_product]->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
              count_region_from_scratch(product[n_product],     /* molecule */
                                        NULL,                   /* rxn pathway */
                                        1,                      /* add count */
                                        NULL,                   /* Location at which to count */
                                        w,                      /* Wall on which this happened */
                                        t);                     /* Time of occurrence */

          }
        }

        /* Otherwise, check if we've crossed the plane. */
        else
        {
          if (product[n_product] == initiator)
          {
            if (product_orient[n_product] != initiatorOrient)
              cross_wall = true;
          }
        }

      }
    }

    /* find out where to place surface products */
    /* Some special cases are listed below. */
    if(num_surface_products == 1)
    {
       if((num_surface_static_reactants == 1) && (num_surface_static_products == 1) && (replace_p1 || replace_p2 || replace_p3))
       {
         /* the lonely static product always replaces the lonely static reactant */
         for(int n_product = rx->n_reactants; n_product < n_players; n_product++)
         {
             if(rx_players[n_product] == NULL) continue;
             if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;
             if(rx_players[n_product]->D != 0) continue;

             if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
             {
                if(replace_p1 && (reacA->properties->D == 0))
                {
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacA)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacA)->grid_index;
                    replace_p1 = 0;
                    break;
                }else if(replace_p2 && (reacB->properties->D == 0)){
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacB)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacB)->grid_index;
                    replace_p2 = 0;
                    break;
                }else if(replace_p3 && (reacC->properties->D == 0)){
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACC_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacC)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacC)->grid_index;
                    replace_p3 = 0;
                    break;
                }
             }
         }
       }
       else if(replace_p1 && replace_p2 && replace_p3)
       {
          /* if all reactants should be  replaced and there is only one
             surface product here we make sure that initiator molecule
             is replaced */
          for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
          {
             if(rx_players[n_product] == NULL) continue;
             if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;

             if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
             {
                 if(reacA == initiator)
                 {
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacA)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacA)->grid_index;
                    replace_p1 = 0;
                 }else if(reacB == initiator){
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacB)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacB)->grid_index;
                    replace_p2 = 0;
                 }else{
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACC_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacC)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacC)->grid_index;
                    replace_p3 = 0;
                 }
                 break;
             }
          }
       }else if(two_to_replace){
         /* replace one of the two reactants randomly */
         while(true)
         {
            rnd_num = rng_uint(world->rng) % (rx->n_reactants);

            if((rnd_num == 0) && replace_p1) break;

            if((rnd_num == 1) && replace_p2) break;

            if((rnd_num == 2) && replace_p3) break;
         }

         for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
         {
             if(rx_players[n_product] == NULL) continue;
             if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;

             if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
             {
                 if(rnd_num == 0)
                 {
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacA)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacA)->grid_index;
                    replace_p1 = 0;
                 }else if(rnd_num == 1){
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacB)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacB)->grid_index;
                    replace_p2 = 0;
                 }else if(rnd_num == 2){
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACC_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacC)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacC)->grid_index;
                    replace_p3 = 0;
                 }
                 break;
             }
         }

       }else if(only_one_to_replace){
            /* no need for a random number here */
          for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
          {
             if(rx_players[n_product] == NULL) continue;
             if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;

             if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
             {
                if(replace_p1)
                {
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacA)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacA)->grid_index;
                    replace_p1 = 0;
                    break;
                }else if(replace_p2){
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacB)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacB)->grid_index;
                    replace_p2 = 0;
                    break;
                }else{
                    product_flag[n_product] = PRODUCT_FLAG_USE_REACC_UV;
                    product_grid[n_product] = ((struct grid_molecule *)reacC)->grid;
                    product_grid_idx[n_product] = ((struct grid_molecule *)reacC)->grid_index;
                    replace_p3 = 0;
                    break;
                }
             }
          }
       }

    }else if(num_surface_products > 1){
       int count;
       if(num_surface_static_reactants > 0)
       {
         if(num_surface_static_products >= num_surface_static_reactants)
         {
           count = 0;
           while(count < num_surface_static_reactants)
           {
             rnd_num = rng_uint(world->rng) % n_players;
             /* pass reactants */
             if (rnd_num < 3) continue;

             if(rx_players[rnd_num] == NULL) continue;
             if((rx_players[rnd_num]->flags & NOT_FREE) == 0) continue;
             if(rx_players[rnd_num]->D != 0) continue;

             if(product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET)
             {
               if((reacA->properties->D == 0) && replace_p1)
               {
                 product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                 product_grid[rnd_num] = ((struct grid_molecule *)reacA)->grid;
                 product_grid_idx[rnd_num] = ((struct grid_molecule *)reacA)->grid_index;
                 replace_p1 = 0;
                 count++;
                 continue;
               }
               if((reacB->properties->D == 0) && replace_p2)
               {
                 product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                 product_grid[rnd_num] = ((struct grid_molecule *)reacB)->grid;
                 product_grid_idx[rnd_num] = ((struct grid_molecule *)reacB)->grid_index;
                 replace_p2 = 0;
                 count++;
                 continue;
               }
               if((reacC->properties->D == 0) && replace_p3)
               {
                 product_flag[rnd_num] = PRODUCT_FLAG_USE_REACC_UV;
                 product_grid[rnd_num] = ((struct grid_molecule *)reacC)->grid;
                 product_grid_idx[rnd_num] = ((struct grid_molecule *)reacC)->grid_index;
                 replace_p3 = 0;
                 count++;
                 continue;
               }
             }
           } /* end while */

         }else{ /*(num_surface_static_products<num_surface_static_reactants)*/
           count = 0;
           while(count < num_surface_static_products)
           {
             rnd_num = rng_uint(world->rng) % n_players;
             /* pass reactants */
             if (rnd_num < 3) continue;

             if(rx_players[rnd_num] == NULL) continue;
             if((rx_players[rnd_num]->flags & NOT_FREE) == 0) continue;
             if(rx_players[rnd_num]->D != 0) continue;

             if(product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET)
             {
               if((reacA->properties->D == 0) && replace_p1)
               {
                 product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                 product_grid[rnd_num] = ((struct grid_molecule *)reacA)->grid;
                 product_grid_idx[rnd_num] = ((struct grid_molecule *)reacA)->grid_index;
                 replace_p1 = 0;
                 count++;
                 continue;
               }
               if((reacB->properties->D == 0) && replace_p2)
               {
                 product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                 product_grid[rnd_num] = ((struct grid_molecule *)reacB)->grid;
                 product_grid_idx[rnd_num] = ((struct grid_molecule *)reacB)->grid_index;
                 replace_p2 = 0;
                 count++;
                 continue;
               }
               if((reacC->properties->D == 0) && replace_p3)
               {
                 product_flag[rnd_num] = PRODUCT_FLAG_USE_REACC_UV;
                 product_grid[rnd_num] = ((struct grid_molecule *)reacC)->grid;
                 product_grid_idx[rnd_num] = ((struct grid_molecule *)reacC)->grid_index;
                 replace_p3 = 0;
                 count++;
                 continue;
               }
             }

           } /* end while */
         }
       }

          /* check whether there are any surface reactants left to replace
             with surface products since we are done with static
             reactants/products */

       if(replace_p1 || replace_p2 || replace_p3)
       {
         /* are there any surface products left that have not replaced yet
             reactants? */

         int surf_prod_left = 0, surf_reactant_left = 0;

         for (int n_product = rx->n_reactants; n_product < n_players; n_product++)
         {
           if(rx_players[n_product] == NULL) continue;
           if((rx_players[n_product]->flags & NOT_FREE) == 0) continue;

           if(product_flag[n_product] == PRODUCT_FLAG_NOT_SET) surf_prod_left++;
         }

         if(replace_p1) surf_reactant_left++;
         if(replace_p2) surf_reactant_left++;
         if(replace_p3) surf_reactant_left++;

         if(surf_prod_left > 0)
         {

           if(surf_prod_left >= surf_reactant_left)
           {
             count = 0;
             while(count < surf_reactant_left)
             {
               rnd_num = rng_uint(world->rng) % n_players;
               /* pass reactants */
               if (rnd_num < 3) continue;

               if(rx_players[rnd_num] == NULL) continue;
               if((rx_players[rnd_num]->flags & NOT_FREE) == 0) continue;

               if(product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET)
               {
                 if(replace_p1)
                 {
                   product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                   product_grid[rnd_num] = ((struct grid_molecule *)reacA)->grid;
                   product_grid_idx[rnd_num] = ((struct grid_molecule *)reacA)->grid_index;
                   replace_p1 = 0;
                   count++;
                   continue;
                 }
                 if(replace_p2)
                 {
                   product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                   product_grid[rnd_num] = ((struct grid_molecule *)reacB)->grid;
                   product_grid_idx[rnd_num] = ((struct grid_molecule *)reacB)->grid_index;
                   replace_p2 = 0;
                   count++;
                   continue;
                 }
                 if(replace_p3)
                 {
                   product_flag[rnd_num] = PRODUCT_FLAG_USE_REACC_UV;
                   product_grid[rnd_num] = ((struct grid_molecule *)reacC)->grid;
                   product_grid_idx[rnd_num] = ((struct grid_molecule *)reacC)->grid_index;
                   replace_p3 = 0;
                   count++;
                   continue;
                 }
               }
             } /* end while */

           }else{   /* surf_prod_left < surf_reactant_left */
             count = 0;
             while(count < surf_prod_left)
             {
               rnd_num = rng_uint(world->rng) % n_players;
               /* pass reactants */
               if (rnd_num < 3) continue;

               if(rx_players[rnd_num] == NULL) continue;
               if((rx_players[rnd_num]->flags & NOT_FREE) == 0) continue;

               if(product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET)
               {
                 if(replace_p1)
                 {
                   product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                   product_grid[rnd_num] = ((struct grid_molecule *)reacA)->grid;
                   product_grid_idx[rnd_num] = ((struct grid_molecule *)reacA)->grid_index;
                   replace_p1 = 0;
                   count++;
                   continue;
                 }
                 if(replace_p2)
                 {
                   product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                   product_grid[rnd_num] = ((struct grid_molecule *)reacB)->grid;
                   product_grid_idx[rnd_num] = ((struct grid_molecule *)reacB)->grid_index;
                   replace_p2 = 0;
                   count++;
                   continue;
                 }
                 if(replace_p3)
                 {
                   product_flag[rnd_num] = PRODUCT_FLAG_USE_REACC_UV;
                   product_grid[rnd_num] = ((struct grid_molecule *)reacC)->grid;
                   product_grid_idx[rnd_num] = ((struct grid_molecule *)reacC)->grid_index;
                   replace_p3 = 0;
                   count++;
                   continue;
                 }
               }

             } /* end while */

           }
         }

       }

    }

    int num_attempts = 0;

    /* all other products are placed on one of the randomly chosen vacant
       tiles */
    for (int n_product = rx->n_reactants; n_product < n_players; ++ n_product)
    {
      /* If the product is a volume product, no placement is required. */
      if (rx_players[n_product]->flags & ON_GRID)
      {
          if(product_flag[n_product] != PRODUCT_FLAG_NOT_SET) continue;

          if(num_vacant_tiles == 0) {
             if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
             if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
             return RX_BLOCKED;
          }

          while(true)
          {
               if(num_attempts > SURFACE_DIFFUSION_RETRIES)
               {
                  if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
                  if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
                  return RX_BLOCKED;
               }

               /* randomly pick a tile from the list */
               rnd_num = rng_uint(world->rng) % num_vacant_tiles;
               tile_idx = -1;
               tile_grid = NULL;

               if(get_tile_neighbor_from_list_of_vacant_neighbors(tile_vacant_nbr_head, rnd_num, &tile_grid, &tile_idx) == 0) {
                   if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
                   if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);
                   return RX_BLOCKED;
               }
               if(tile_idx < 0) continue; /* this tile was checked out before */

               assert(tile_grid != NULL);

               if(all_inside_restricted_boundary)
               {
                 /* if this tile is not inside the restricted boundary
                    - try again */
                   int cond_1 = 0, cond_2 = 0, cond_3 = 0;
                   cond_1 = (!wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1));
                   cond_2 = (!wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2));
                   cond_3 = (!wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3));

                   if(cond_1 || cond_2 || cond_3)
                   {
                      uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                      num_attempts++;
                      continue;
                   }
               }
               else if(all_outside_restricted_boundary){
                 int cond_1 = 0, cond_2 = 0, cond_3 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1);
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2);
                 cond_3 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3);

                 if(cond_1 || cond_2 || cond_3)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(grid_1_inside_grid_2_inside_grid_3_outside){
                 int cond_1 = 0, cond_2 = 0, cond_3 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1));
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2));
                 cond_3 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3);

                 if(cond_1 || cond_2 || cond_3)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(grid_1_inside_grid_2_outside_grid_3_inside){
                 int cond_1 = 0, cond_2 = 0, cond_3 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1));
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2);
                 cond_3 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3));

                 if(cond_1 || cond_2 || cond_3)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(grid_1_inside_grid_2_outside_grid_3_outside){
                 int cond_1 = 0, cond_2 = 0, cond_3 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1));
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2);
                 cond_3 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3);

                 if(cond_1 || cond_2 || cond_3)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(grid_1_outside_grid_2_inside_grid_3_outside){
                 int cond_1 = 0, cond_2 = 0, cond_3 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1);
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2));
                 cond_3 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3);

                 if(cond_1 || cond_2 || cond_3)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(grid_1_outside_grid_2_inside_grid_3_inside){
                 int cond_1 = 0, cond_2 = 0, cond_3 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1);
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2));
                 cond_3 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3));

                 if(cond_1 || cond_2 || cond_3)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(grid_1_outside_grid_2_outside_grid_3_inside){
                 int cond_1 = 0, cond_2 = 0, cond_3 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1);
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2);
                 cond_3 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3));

                 if(cond_1 || cond_2 || cond_3)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_grid_2_inside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1));
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2));

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_inside_grid_2_outside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1));
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2);

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_outside_grid_2_inside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1);
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2));

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_grid_2_outside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1);
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2);

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_grid_3_inside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1));
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3));

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_inside_grid_3_outside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1));
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3);

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_outside_grid_3_inside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1);
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3));

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_grid_3_outside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1);
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3);

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_2_grid_3_inside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2));
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3));

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_2_inside_grid_3_outside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2));
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3);

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_2_outside_grid_3_inside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2);
                 cond_2 = !(wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3));

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_2_grid_3_outside){
                 int cond_1 = 0, cond_2 = 0;
                 cond_1 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2);
                 cond_2 = wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3);

                 if(cond_1 || cond_2)
                 {
                    uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                    num_attempts++;
                    continue;
                 }
               }
               else if(only_grid_1_inside){
                   if(!wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_1))
                   {
                     uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                     num_attempts++;
                     continue;
                   }
               }else if(only_grid_1_outside) {
                   if(wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_1))
                   {
                     uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                     num_attempts++;
                     continue;
                   }
               }else if(only_grid_2_inside) {
                   if(!wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_2))
                   {
                     uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                     num_attempts++;
                     continue;
                   }
               }else if(only_grid_2_outside) {
                   if(wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_2))
                   {
                     uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                     num_attempts++;
                     continue;
                   }
               }else if(only_grid_3_inside) {
                   if(!wall_belongs_to_all_regions_in_region_list(tile_grid->surface, rlp_head_wall_3))
                   {
                     uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                     num_attempts++;
                     continue;
                   }
               }else if(only_grid_3_outside) {
                   if(wall_belongs_to_any_region_in_region_list(tile_grid->surface, rlp_head_obj_3))
                   {
                     uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
                     num_attempts++;
                     continue;
                   }
               }

               product_grid[n_product]     = tile_grid;
               product_grid_idx[n_product] = tile_idx;
               product_flag[n_product]     = PRODUCT_FLAG_USE_RANDOM;

               break;
          }
      }
    }

  } /* end if(is_orientable) */


  /* Determine the location of the reaction for count purposes. */
  struct vector3 count_pos_xyz;
  if (hitpt != NULL) count_pos_xyz = *hitpt;
  else if (grid_reactant) uv2xyz(& grid_reactant->s_pos, grid_reactant->grid->surface, & count_pos_xyz);
  else count_pos_xyz = ((struct volume_molecule *) reacA)->pos;

  /* Create and place each product. */
  struct vector3 mol_pos_tmp;
  struct subvolume *product_subvol = NULL;
  for (int n_product = rx->n_reactants; n_product < n_players; ++ n_product)
  {
    struct abstract_molecule *this_product = NULL;
    struct species * const product_species = rx_players[n_product];

    /* If the product is a grid molecule, place it on the grid. */
    if (product_species->flags & ON_GRID)
    {
        struct vector2 prod_uv_pos;

        /* Pick an appropriate position for the new molecule. */
        if (world->randomize_gmol_pos)
        {
          switch (product_flag[n_product])
          {
            case PRODUCT_FLAG_USE_REACA_UV:
              prod_uv_pos = ((struct grid_molecule*) reacA)->s_pos;
              break;

            case PRODUCT_FLAG_USE_REACB_UV:
              prod_uv_pos = ((struct grid_molecule*) reacB)->s_pos;
              break;

            case PRODUCT_FLAG_USE_REACC_UV:
              prod_uv_pos = ((struct grid_molecule*) reacC)->s_pos;
              break;

            case PRODUCT_FLAG_USE_RANDOM:
              grid2uv_random(product_grid[n_product], product_grid_idx[n_product], &prod_uv_pos);
              break;

            default:
              UNHANDLED_CASE(product_flag[n_product]);
              break;
          }
        }
        else
          grid2uv(product_grid[n_product], product_grid_idx[n_product], & prod_uv_pos);

        this_product = (struct abstract_molecule *)
              place_grid_product(product_species,
                                 product_grid[n_product],
                                 product_grid_idx[n_product],
                                 & prod_uv_pos,
                                 product_orient[n_product],
                                 t);
    }

    /* else place the molecule in space. */
    else
    {
        /* Unless this is a unimolecular reaction, we will have a hitpt. */
        if (! hitpt)
        {
          /* If this is a unimolecular surface rxn... */
          if (reacA->properties->flags & ON_GRID)
          {
            uv2xyz(& ((struct grid_molecule *) reacA)->s_pos,
                   ((struct grid_molecule *) reacA)->grid->surface,
                   & mol_pos_tmp);
            product_subvol = find_subvolume(&mol_pos_tmp, last_subvol);
          }

          /* ... else a unimolecular volume rxn. */
          else
          {
            mol_pos_tmp = ((struct volume_molecule *) reacA)->pos;
            product_subvol = ((struct volume_molecule *) reacA)->subvol;
          }
          hitpt = & mol_pos_tmp;
        }
        else product_subvol = find_subvolume(hitpt, last_subvol);

        this_product = (struct abstract_molecule *) place_volume_product(product_species,
                                                                         grid_reactant,
                                                                         w,
                                                                         product_subvol,
                                                                         hitpt,
                                                                         product_orient[n_product],
                                                                         t);

      if (((struct volume_molecule *) this_product)->index < DISSOCIATION_MAX)
        update_dissociation_index = true;
    }

    /* Update molecule counts */
    ++ product_species->population;
    if (product_species->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
      count_region_from_scratch(this_product, NULL, 1, NULL, NULL, t);
  }

  /* If necessary, update the dissociation index. */
  if (update_dissociation_index)
  {
    if (-- world->dissociation_index < DISSOCIATION_MIN)
      world->dissociation_index = DISSOCIATION_MAX;
  }

  /* Handle events triggered off of named reactions */
  if (rx->info[path].pathname!=NULL)
  {
    /* No flags for reactions so we have to check regions if we have waypoints! Fix to be more efficient for WORLD-only counts? */
    if (world->place_waypoints_flag)
      count_region_from_scratch(NULL, rx->info[path].pathname, 1, &count_pos_xyz, w, t);

    /* Other magical stuff.  For now, can only trigger releases. */
    if (rx->info[path].pathname->magic!=NULL)
    {
      if (reaction_wizardry(rx->info[path].pathname->magic, w, &count_pos_xyz, t))
        mcell_allocfailed("Failed to complete reaction triggered release after a '%s' reaction.",
                          rx->info[path].pathname->sym->name);
    }
  }

  if(tile_nbr_head != NULL) delete_tile_neighbor_list(tile_nbr_head);
  if(tile_vacant_nbr_head != NULL) delete_tile_neighbor_list(tile_vacant_nbr_head);

  return cross_wall ? RX_FLIP : RX_A_OK;
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
    if(rx->is_complex)
    {
       result = outcome_products(NULL, NULL, t, rx, path, reac, NULL, 0, 0);
    }else{
       result = outcome_products_random(NULL, NULL, t, rx, path, reac, NULL, 0, 0);
    }
  }
  else
  {
    g = (struct grid_molecule*) reac;
    if(rx->is_complex)
    {
       result = outcome_products(g->grid->surface, NULL, t, rx, path, reac, NULL, g->orient, 0);

    }else{
       /* we will not create products if the reaction is with an ABSORPTIVE
          region border */

       if((strcmp(rx->players[0]->sym->name, "ALL_SURFACE_MOLECULES") == 0) ||
          (strcmp(rx->players[0]->sym->name, "ALL_MOLECULES") == 0))
       {
           /* do nothing */
       }else{
         result = outcome_products_random(g->grid->surface, NULL, t, rx, path, reac, NULL, g->orient, 0);
       }
    }
  }

  if (result==RX_BLOCKED) return RX_BLOCKED;

  if (result != RX_BLOCKED) {
     rx->info[path].count++;
     rx->n_occurred++;
  }

  who_am_i = rx->players[rx->product_idx[path]];

  if (who_am_i == NULL)
  {
    if (m != NULL)
    {
      m->subvol->mol_count--;
      if (m->flags & IN_SCHEDULE) m->subvol->local_storage->timer->defunct_count++;
      if (m->properties->flags&COUNT_SOME_MASK)
      {
        count_region_from_scratch((struct abstract_molecule*)m, NULL, -1, &(m->pos), NULL, m->t);
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
        count_region_from_scratch((struct abstract_molecule*)g, NULL, -1, NULL, NULL, g->t);
      }
    }

    who_was_i->n_deceased++;
    who_was_i->cum_lifetime += t - reac->birthday;

    who_was_i->population--;
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
  int killA, killB;

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

  if(rx->is_complex)
  {
     result = outcome_products(w, hitpt, t, rx, path, reacA, reacB, orientA, orientB);
  }else{
    result = outcome_products_random(w, hitpt, t, rx, path, reacA, reacB, orientA, orientB);
  }

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
      count_region_from_scratch(reacB, NULL, -1, NULL, NULL, t);
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
        count_region_from_scratch(reacA, NULL, -1, NULL, NULL, t);
      }
    }
    else if (reacA->flags&COUNT_ME)
    {
      /* Subtlety: we made it up to hitpt, but our position is wherever we were before that! */
      if (hitpt==NULL || reacB_was_free || (reacB->properties!=NULL && (reacB->properties->flags&NOT_FREE)==0))
      {
        /* Vol-vol rx should be counted at hitpt */
        count_region_from_scratch(reacA, NULL, -1, hitpt, NULL, t);
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

        count_region_from_scratch(reacA, NULL, -1, &fake_hitpt, NULL, t);
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
  int result;
  /* flags */
  int killA = 0,killB = 0, killC = 0;
  int reacA_is_free = 0;
  int reacB_is_free = 0;
  int reacC_is_free = 0;
  int num_surface_reactants = 0;

   if ((reacA->properties->flags & NOT_FREE) == 0)
   {
       reacA_is_free = 1;
   }else num_surface_reactants++;

   if ((reacB->properties->flags & NOT_FREE) == 0) reacB_is_free = 1;
   else num_surface_reactants++;

   if ((reacC->properties->flags & NOT_FREE) == 0) reacC_is_free = 1;
   else num_surface_reactants++;

    if(!reacA_is_free)
    {
       g = (struct grid_molecule *)reacA;
    }else if(!reacB_is_free){
       g = (struct grid_molecule *)reacB;
    }else if(!reacC_is_free){
       g = (struct grid_molecule *)reacC;
    }
    if(g != NULL) w = g->grid->surface;

    result = outcome_products_trimol_reaction_random(w,hitpt,t,rx,path,reacA, reacB, reacC, orientA, orientB, orientC);
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
      count_region_from_scratch(reacC, NULL, -1, NULL, NULL, t);
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
      count_region_from_scratch(reacB, NULL, -1, NULL, NULL, t);
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
        count_region_from_scratch(reacA, NULL, -1, NULL, NULL, t);
      }
    }
    else if ((reacA->flags&COUNT_ME) && world->place_waypoints_flag)
    {
      /* Subtlety: we made it up to hitpt, but our position is wherever we were before that! */
      if (hitpt==NULL || (reacB_is_free && reacC_is_free))
          /* Vol-vol-vol rx should be counted at hitpt */
      {
        count_region_from_scratch(reacA, NULL, -1, hitpt, NULL, t);
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

        count_region_from_scratch(reacA, NULL, -1, &fake_hitpt, NULL, t);
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

    /* If reaction object has ALL_MOLECULES or ALL_VOLUME_MOLECULES
       as the first reactant it means that reaction is of the type
       ABSORPTIVE = ALL_MOLECULES or ABSORPTIVE = ALL_VOLUME_MOLECULES
       since other cases (REFLECTIVE/TRANSPARENT) are taken care above.
       But there are no products for this reaction, so we do no need
       to go into "outcome_products()" function. */

           /*
    if(strcmp(rx->players[0]->sym->name, "ALL_MOLECULES") != 0)
    {
      if(rx->is_complex)
      {
         result = outcome_products(surface, hitpt, t, rx, path, reac, NULL, orient, 0);
      }else{
         result = outcome_products_random(surface, hitpt, t, rx, path, reac, NULL, orient, 0);
      }
    }else result = RX_DESTROY;
          */
    if((strcmp(rx->players[0]->sym->name, "ALL_MOLECULES") == 0) ||
       (strcmp(rx->players[0]->sym->name, "ALL_VOLUME_MOLECULES") == 0))
    {
       result = RX_DESTROY;
    }else{
      if(rx->is_complex)
      {
         result = outcome_products(surface, hitpt, t, rx, path, reac, NULL, orient, 0);
      }else{
         result = outcome_products_random(surface, hitpt, t, rx, path, reac, NULL, orient, 0);
      }
    }
    if (result == RX_BLOCKED) return RX_A_OK; /* reflect the molecule */

    rx->info[path].count++;
    rx->n_occurred++;

    if (rx->players[idx] == NULL)
    {
      /* The code below is also valid for the special reaction
         of the type ABSORPTIVE = ALL_MOLECULES (or ALL_VOLUME_MOLECULES) */
      m->subvol->mol_count--;
      if (world->place_waypoints_flag  &&  (reac->flags&COUNT_ME))
      {
        if (hitpt==NULL)
        {
          count_region_from_scratch(reac, NULL, -1, NULL, NULL, t);
        }
        else
        {
          struct vector3 fake_hitpt;

          /* Halfway in between where we were and where we react should be a safe away-from-wall place to remove us */
          if (loc_okay==NULL) loc_okay=&(m->pos);
          fake_hitpt.x = 0.5*hitpt->x + 0.5*loc_okay->x;
          fake_hitpt.y = 0.5*hitpt->y + 0.5*loc_okay->y;
          fake_hitpt.z = 0.5*hitpt->z + 0.5*loc_okay->z;

          count_region_from_scratch(reac, NULL, -1, &fake_hitpt, NULL, t);
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



/************************************************************************
 *
 * this function determines where reactants grid1 and grid2 are located
 * (inside/outside) with respect to their restrictive region borders if
 * they have any.
 *
 * in: grid molecule 1 (located on wall 1)
 *     grid molecule 2 (located on wall 2)
 *     pointer to array with restrictive regions which contain wall 1
 *     pointer to array with restrictive regions which contain wall 2
 *     pointer to array with restrictive regions which don't contain wall 1
 *     pointer to array with restrictive regions which don't contain wall 2
 *
 * out: the 4 arrays with pointers to restrictive regions will be filled
 *      and returned
 *
 ***********************************************************************/
int determine_molecule_region_topology(struct grid_molecule *grid_1,
                                       struct grid_molecule *grid_2,
                                       struct region_list **rlp_wall_1_ptr,
                                       struct region_list **rlp_wall_2_ptr,
                                       struct region_list **rlp_obj_1_ptr,
                                       struct region_list **rlp_obj_2_ptr,
                                       bool is_unimol)
{
  int grid_bitmask = 0;
  struct wall *w_1, *w_2;
  struct region_list *rlp_head_wall_1 = NULL;
  struct region_list *rlp_head_wall_2 = NULL;
  struct region_list *rlp_head_obj_1 = NULL;
  struct region_list *rlp_head_obj_2 = NULL;

  /* bimolecular reactions */
  if((grid_1 != NULL) && (grid_2 != NULL))
  {
    /* both reactants have restrictive region borders */
    if((grid_1->properties->flags & CAN_REGION_BORDER) &&
      (grid_2->properties->flags & CAN_REGION_BORDER) &&
      are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1) &&
      are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2))
    {
      w_1 = grid_1->grid->surface;
      w_2 = grid_2->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(w_1, grid_1);
      rlp_head_wall_2 = find_restricted_regions_by_wall(w_2, grid_2);

      /* both reactants are inside their respective restricted regions */
      if((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL))
      {
        grid_bitmask |= ALL_INSIDE;
      }
      /* both reactants are outside their respective restricted regions */
      else if((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL))
      {
        rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
        rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
        grid_bitmask |= ALL_OUTSIDE;
      }
      /* grid1 is inside and grid2 is outside of its respective
       * restrictive region */
      else if((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 == NULL))
      {
        rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
        grid_bitmask |= GRID1_IN_GRID2_OUT;
      }
      /* grid2 is inside and grid1 is outside of its respective
       * restrictive region */
      else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL))
      {
        rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
        grid_bitmask |= GRID1_OUT_GRID2_IN;
      }
    }

    /* only reactant grid_1 has restrictive region border property */
    else if((grid_1->properties->flags & CAN_REGION_BORDER) &&
            are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1) &&
            (!(grid_2->properties->flags & CAN_REGION_BORDER) ||
             !are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2)))
    {
      w_1 = grid_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(w_1, grid_1);
      if(rlp_head_wall_1 != NULL)
      {
        grid_bitmask |= GRID1_IN;
      }
      else
      {
        rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
        grid_bitmask |= GRID1_OUT;
      }
    }

    /* only reactant "grid_2" has restrictive region border property */
    else if ((grid_2->properties->flags & CAN_REGION_BORDER) &&
             are_restricted_regions_for_species_on_object(grid_2->grid->surface->parent_object, grid_2)  &&
            (!(grid_1->properties->flags & CAN_REGION_BORDER) ||
             !are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1)))
    {
      w_2 = grid_2->grid->surface;
      rlp_head_wall_2 = find_restricted_regions_by_wall(w_2, grid_2);
      if(rlp_head_wall_2 != NULL)
      {
        grid_bitmask |= GRID2_IN;
      }
      else
      {
        rlp_head_obj_2 = find_restricted_regions_by_object(w_2->parent_object, grid_2);
        grid_bitmask |= GRID2_OUT;
      }
    }
  }

  /* unimolecular reactions */
  else if((grid_1 != NULL) && is_unimol)
  {
    if ((grid_1->properties->flags & CAN_REGION_BORDER) &&
       are_restricted_regions_for_species_on_object(grid_1->grid->surface->parent_object, grid_1))
    {
      w_1 = grid_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(w_1, grid_1);
      if(rlp_head_wall_1 != NULL)
      {
        grid_bitmask |= ALL_INSIDE;
      }
      else
      {
        rlp_head_obj_1 = find_restricted_regions_by_object(w_1->parent_object, grid_1);
        grid_bitmask |= ALL_OUTSIDE;
      }
    }
  }

  *rlp_wall_1_ptr = rlp_head_wall_1;
  *rlp_wall_2_ptr = rlp_head_wall_2;
  *rlp_obj_1_ptr = rlp_head_obj_1;
  *rlp_obj_2_ptr = rlp_head_obj_2;

  return grid_bitmask;
}



/***********************************************************************
 *
 * this function tests if wall target can be reached for product placement
 * based on the previously stored reactant topology based on
 * grid_bitmask. Below, wall 1 is the wall containing reactant 1 and
 * wall 2 is the wall containing reactant 2.
 *
 * in: wall to test for product placement
 *     pointer to array with regions that contain wall 1
 *     pointer to array with regions that contain wall 2
 *     pointer to array with regions that do not contain wall 1
 *     pointer to array with regions that do not contain wall 2
 *
 * out: returns true or false depending if wall target can be
 *      used for product placement.
 *
 ***********************************************************************/
bool product_tile_can_be_reached(struct wall *target,
                                 struct region_list *rlp_head_wall_1,
                                 struct region_list *rlp_head_wall_2,
                                 struct region_list *rlp_head_obj_1,
                                 struct region_list *rlp_head_obj_2,
                                 int grid_bitmask,
                                 bool is_unimol)
{
  bool status = true;

  if(grid_bitmask & ALL_INSIDE)
  {
    if(is_unimol)
    {
      if(!wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_1))
      {
        status = false;
      }
    }
    else
    {
      /* bimol reaction */
      if(!wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_1) ||
         !wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_2))
      {
        status = false;
      }
    }
  }
  else if(grid_bitmask & ALL_OUTSIDE)
  {
    if(is_unimol)
    {
      if(wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_1))
      {
        status = false;
      }
    }
    else
    {
      if(wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_1) ||
         wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_2))
      {
        status = false;
      }
    }
  }
  else if(grid_bitmask & GRID1_IN_GRID2_OUT)
  {
    if(!wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_1) ||
       wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_2))
    {
      status = false;
    }
  }
  else if(grid_bitmask & GRID1_OUT_GRID2_IN)
  {
    if(wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_1) ||
      !wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_2))
    {
      status = false;
    }
  }
  else if(grid_bitmask & GRID1_IN)
  {
    if(!wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_1))
    {
      status = false;
    }
  }
  else if(grid_bitmask & GRID1_OUT)
  {
    if(wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_1))
    {
      status = false;
    }
  }
  else if(grid_bitmask & GRID2_IN)
  {
    if(!wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_2))
    {
      status = false;
    }
  }
  else if(grid_bitmask & GRID2_OUT)
  {
    if(wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_2))
    {
      status = false;
    }
  }

  return status;
}

