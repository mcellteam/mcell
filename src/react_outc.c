/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#include "config.h"

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "logging.h"
#include "rng.h"
#include "util.h"
#include "grid_util.h"
#include "count_util.h"
#include "react.h"
#include "vol_util.h"
#include "wall_util.h"
#include "diffuse.h"

static int outcome_products_random(struct volume *world, struct wall *w,
                                   struct vector3 *hitpt, double t,
                                   struct rxn *rx, int path,
                                   struct abstract_molecule *reacA,
                                   struct abstract_molecule *reacB,
                                   short orientA, short orientB);


static int cleanup_and_block_rx(struct tile_neighbor *tn1, struct tile_neighbor *tn2);

int is_compatible_surface(void *req_species, struct wall *w) {
  struct surf_class_list *scl, *scl2;

  struct surf_class_list *rs_head = (struct surf_class_list *)req_species;

  if (rs_head == NULL)
    return 1;

  for (scl = w->surf_class_head; scl != NULL; scl = scl->next) {
    for (scl2 = rs_head; scl2 != NULL; scl2 = scl2->next) {
      if (scl->surf_class == scl2->surf_class)
        return 1;
    }
  }

  return 0;
}

void add_reactants_to_product_list(struct rxn *rx, struct abstract_molecule *reacA,
  struct abstract_molecule *reacB, struct abstract_molecule *reacC,
  struct abstract_molecule **player, char *player_type) {

  /* Add reacA to the list of players, saving the reactant's type. */
  player[0] = reacA;
  player_type[0] = IS_SURF_MOL(reacA) ? PLAYER_SURF_MOL : PLAYER_VOL_MOL;

  /* If we have a second reactant, add it to the list of players. */
  if (rx->n_reactants > 1) {

    /* If the second reactant is a wall... */
    if (reacB == NULL) {
      assert(rx->n_reactants == 2);
      player[1] = NULL;
      player_type[1] = PLAYER_WALL;
    } else {  // else 2nd reactant is a wall
      player[1] = reacB;
      player_type[1] = IS_SURF_MOL(reacB) ? PLAYER_SURF_MOL : PLAYER_VOL_MOL;
    }

    if (rx->n_reactants > 2) {
      if (reacC == NULL) {
        /* it's a wall. */
        player[2] = NULL;
        player_type[2] = PLAYER_WALL;
      } else {
        player[2] = reacC;
        player_type[2] = IS_SURF_MOL(reacC) ? PLAYER_SURF_MOL : PLAYER_VOL_MOL;
      }
    }
  }
}

static bool is_rxn_unimol(struct rxn *rx) {
  if (rx->n_reactants == 1)
    return true;

  if (rx->n_reactants != 2)
    return false;

  if (!(rx->players[0]->flags & ON_GRID))
    return false;

  return (rx->players[1]->flags & IS_SURFACE) != 0;
}

void tiny_diffuse_3D(
    struct volume *world,
    struct subvolume *subvol,
    struct vector3 *displacement,
    struct vector3 *pos,
    struct wall *w) {

  struct vector3 temp_displacement = {
    .x = displacement->x,
    .y = displacement->y,
    .z = displacement->z
  };
  struct collision *shead = ray_trace(
      world, pos, NULL, subvol, &temp_displacement, w);
  if (shead->next != NULL) {
    shead = (struct collision *)ae_list_sort((struct abstract_element *)shead);
  }

  struct collision *smash = NULL;
  for (smash = shead; smash != NULL; smash = smash->next) {
    if ((smash->what & COLLIDE_WALL) != 0) {
      vectorize(pos, &(smash->loc), displacement);
      scalar_prod(displacement, 0.5, displacement);
      break;
    }
  }
  pos->x += displacement->x;
  pos->y += displacement->y;
  pos->z += displacement->z;
  subvol = find_subvolume(world, pos, subvol);
}

struct volume_molecule *
place_volume_product(struct volume *world, struct species *product_species,
                     struct surface_molecule *sm_reactant, struct wall *w,
                     struct subvolume *subvol, struct vector3 *hitpt,
                     short orient, double t, struct periodic_image *periodic_box) {
  struct vector3 pos = *hitpt;

  /* For an orientable reaction, we need to move products away from the surface
   * to ensure they end up on the correct side of the plane. */
  if (w) {
    double bump = (orient > 0) ? EPS_C : -EPS_C;
    struct vector3 displacement = { .x = 2 * bump * w->normal.x,
                                    .y = 2 * bump * w->normal.y,
                                    .z = 2 * bump * w->normal.z,
                                  };
    tiny_diffuse_3D(world, subvol, &displacement, &pos, w);
  }

  /* Allocate and initialize the molecule. */
  struct volume_molecule *new_volume_mol;
  new_volume_mol =
      CHECKED_MEM_GET(subvol->local_storage->mol, "volume molecule");
  new_volume_mol->birthplace = subvol->local_storage->mol;
  new_volume_mol->birthday = convert_iterations_to_seconds(
      world->start_iterations, world->time_unit,
      world->simulation_start_seconds, t);
  new_volume_mol->id = world->current_mol_id++;
  new_volume_mol->t = t;
  new_volume_mol->t2 = 0.0;

  new_volume_mol->periodic_box = CHECKED_MALLOC_STRUCT(struct periodic_image,
    "periodic image descriptor");
  new_volume_mol->periodic_box->x = periodic_box->x;
  new_volume_mol->periodic_box->y = periodic_box->y;
  new_volume_mol->periodic_box->z = periodic_box->z;

  new_volume_mol->properties = product_species;
  new_volume_mol->prev_v = NULL;
  new_volume_mol->next_v = NULL;
  new_volume_mol->pos = pos;
  new_volume_mol->subvol = subvol;
  new_volume_mol->index = 0;
  new_volume_mol->flags = TYPE_VOL | ACT_NEWBIE | IN_VOLUME | IN_SCHEDULE;
  if (product_species->space_step > 0.0)
    new_volume_mol->flags |= ACT_DIFFUSE;
  if ((product_species->flags & COUNT_SOME_MASK) != 0)
    new_volume_mol->flags |= COUNT_ME;

  /* Check whether the product can undergo unimolecular rxns; if so, mark it. */
  if (trigger_unimolecular(world->reaction_hash, world->rx_hashsize,
                           product_species->hashval,
                           (struct abstract_molecule *)new_volume_mol) != NULL)
    new_volume_mol->flags |= ACT_REACT;

  /* If this product resulted from a surface rxn, store the previous wall
   * position. */
  if (sm_reactant && distinguishable(product_species->D, 0, EPS_C)) {
    new_volume_mol->previous_wall = sm_reactant->grid->surface;

    /* This will be overwritten with orientation in the CLAMPED/surf.
     * reversibility case
     */
    new_volume_mol->index = sm_reactant->grid_index;
  }

  /* Else clear the previous wall position. */
  else {
    new_volume_mol->previous_wall = NULL;
    new_volume_mol->index = -1;
  }

  /* Set reversibility state for the new molecule. */
  if (w) {
    if (world->surface_reversibility) {
      /* Which direction did we move? */
      new_volume_mol->previous_wall = w;
      new_volume_mol->index = (orient > 0) ? 1 : -1;
      new_volume_mol->flags |= ACT_CLAMPED;
    }
  } else if (world->volume_reversibility) {
    new_volume_mol->index = world->dissociation_index;
    new_volume_mol->flags |= ACT_CLAMPED;
  }

  /* Add the molecule to the subvolume */
  ht_add_molecule_to_list(&new_volume_mol->subvol->mol_by_species,
                          new_volume_mol);
  ++new_volume_mol->subvol->mol_count;

  /* Add to the schedule. */
  if (schedule_add(subvol->local_storage->timer, new_volume_mol))
    mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                      product_species->sym->name);
  return new_volume_mol;
}

struct surface_molecule *
place_sm_product(struct volume *world, struct species *product_species,
                 struct surface_grid *grid, int grid_index,
                 struct vector2 *mol_uv_pos, short orient, double t,
                 struct periodic_image *periodic_box) {
  struct vector3 mol_xyz_pos;
  uv2xyz(mol_uv_pos, grid->surface, &mol_xyz_pos);
  struct subvolume *sv = find_subvolume(world, &mol_xyz_pos, grid->subvol);

  /* Allocate and initialize the molecule. */
  struct surface_molecule *new_surf_mol;
  new_surf_mol = CHECKED_MEM_GET(sv->local_storage->smol, "surface molecule");
  new_surf_mol->birthplace = sv->local_storage->smol;
  new_surf_mol->birthday = convert_iterations_to_seconds(
      world->start_iterations, world->time_unit,
      world->simulation_start_seconds, t);
  new_surf_mol->id = world->current_mol_id++;
  new_surf_mol->t = t;
  new_surf_mol->t2 = 0.0;
  new_surf_mol->properties = product_species;
  new_surf_mol->periodic_box = CHECKED_MALLOC_STRUCT(struct periodic_image,
    "periodic image descriptor");
  new_surf_mol->periodic_box->x = periodic_box->x;
  new_surf_mol->periodic_box->y = periodic_box->y;
  new_surf_mol->periodic_box->z = periodic_box->z;

  new_surf_mol->flags = TYPE_SURF | ACT_NEWBIE | IN_SCHEDULE;
  if (product_species->space_step > 0)
    new_surf_mol->flags |= ACT_DIFFUSE;
  if (product_species->flags & COUNT_ENCLOSED)
    new_surf_mol->flags |= COUNT_ME;
  new_surf_mol->grid = grid;
  new_surf_mol->grid_index = grid_index;
  new_surf_mol->s_pos = *mol_uv_pos;
  new_surf_mol->orient = orient;

  /* Check whether the product can undergo unimolecular rxns; if so, mark it. */
  if (trigger_unimolecular(world->reaction_hash, world->rx_hashsize,
                           product_species->hashval,
                           (struct abstract_molecule *)new_surf_mol) != NULL ||
      (product_species->flags & CAN_SURFWALL) != 0)
    new_surf_mol->flags |= ACT_REACT;

  /* Add to the grid. */
  ++grid->n_occupied;
  if (grid->sm_list[grid_index]) {
    remove_surfmol_from_list(
        &grid->sm_list[grid_index], grid->sm_list[grid_index]->sm);
  }
  grid->sm_list[grid_index] = add_surfmol_with_unique_pb_to_list(
    grid->sm_list[grid_index], new_surf_mol);

  /* Add to the schedule. */
  if (schedule_add(sv->local_storage->timer, new_surf_mol))
    mcell_allocfailed("Failed to add newly created %s molecule to scheduler.",
                      product_species->sym->name);

  return new_surf_mol;
}

/***************************************************************************
outcome_products_random:
   In: world: simulation state
       w: first wall in the reaction
       hitpt: hit point (if any)
       t: time of the reaction
       rx: reaction
       path: path of the reaction
       reacA: first reactant (moving molecule)
       reacB: second reactant
       orientA: orientation of the first reactant
       orientB: orientation of the second reactant
   Out: Returns RX_A_OK, RX_FLIP or RX_BLOCKED.
Note: This function replaces surface reactants (if needed) by the surface
      products picked in the random order from the list of products. Also
      surface products are placed in the random order in the surrounding empty
      tiles. After this function execution some walls that do not have surface
      molecules and therefore do not have a grid may get a grid as side effect
      of calling functions "grid_all_neigbors_across_walls_through_vertices()"
      and "grid_all_neighbors_across_walls_through_edges()".
Note: If both reactants are surface molecules, and they are both located within
      the same restricted region border (REFL/ABSORB), then reaction products -
      surface molecules for which this region border is restrictive will be
      placed inside this region. If any of the conditions above are not true,
      the reaction products - surface molecules can be placed on any tile from
      the list of vacant tiles.
Note: Policy on surface products placement is described in the document
      "policy_surf_products_placement.doc" (see "src/docs").

****************************************************************************/
static int outcome_products_random(struct volume *world, struct wall *w,
                                   struct vector3 *hitpt, double t,
                                   struct rxn *rx, int path,
                                   struct abstract_molecule *reacA,
                                   struct abstract_molecule *reacB,
                                   short orientA, short orientB) {

  /* Did the moving molecule cross the plane? */
  bool cross_wall = false; 

  /* index of the first player for the pathway */
  int const i0 = rx->product_idx[path]; 
  /* index of the first player for the next pathway */
  int const iN = rx->product_idx[path + 1]; 
  assert(iN > i0);
  /* Players array from the reaction. */
  struct species **rx_players = rx->players + i0; 

  int const n_players = iN - i0;                /* number of reaction players */
  struct abstract_molecule *product[n_players]; /* array of products */
  /* array that decodes the type of each product */
  char product_type[n_players]; 
  short product_orient[n_players]; /* array of orientations for each product */
  /* array of surface_grids for products */
  struct surface_grid *product_grid[n_players];
  int product_grid_idx[n_players]; /* array of grid indices for products */
  byte product_flag[n_players];    /* array of placement flags for products */

  /* Unimol rxn (not mol-mol, not mol-wall) */
  bool const is_unimol = is_rxn_unimol(rx);

  struct surface_grid *tile_grid; /* surface grid the tile belongs to */
  int num_vacant_tiles = 0;       /* number of vacant tiles */

  /* used for product placement for the reaction of type A->B+C[rate] */
  unsigned int reac_idx = UINT_MAX;
  struct surface_grid *reac_grid = NULL, *mol_grid = NULL;

  /* Clear the initial product info. */
  for (int i = 0; i < n_players; ++i) {
    product[i] = NULL;
    product_type[i] = PLAYER_NONE;
    product_orient[i] = 0;
    product_grid[i] = NULL;
    product_grid_idx[i] = -1;
    product_flag[i] = PRODUCT_FLAG_NOT_SET;
  }

  /* Flag indicating that a surface is somehow involved with this reaction. */
  struct surface_molecule *const sm_1 =
      IS_SURF_MOL(reacA) ? (struct surface_molecule *)reacA : NULL;
  struct surface_molecule *const sm_2 =
      IS_SURF_MOL(reacB) ? (struct surface_molecule *)reacB : NULL;
  struct surface_molecule *const sm_reactant = sm_1 ? sm_1 : sm_2;
  bool const is_orientable = (w != NULL) || (sm_reactant != NULL);

  /* list of the restricted regions for the reactants by wall */
  struct region_list *rlp_head_wall_1 = NULL, *rlp_head_wall_2 = NULL;

  /* list of the restricted regions for the reactants by object */
  struct region_list *rlp_head_obj_1 = NULL, *rlp_head_obj_2 = NULL;

  int sm_bitmask = determine_molecule_region_topology(
      world, sm_1, sm_2, &rlp_head_wall_1, &rlp_head_wall_2, &rlp_head_obj_1,
      &rlp_head_obj_2, is_unimol);

  /* reacA is the molecule which initiated the reaction. */
  struct abstract_molecule *const initiator = reacA;
  short const initiatorOrient = orientA;

  /* Ensure that reacA and reacB are sorted in the same order as the rxn players. */
  assert(reacA != NULL);

  if (reacA->properties != rx->players[0]) {
    struct abstract_molecule *tmp_mol = reacA;
    reacA = reacB;
    reacB = tmp_mol;
    short tmp_orient = orientA;
    orientA = orientB;
    orientB = tmp_orient;
  }
  assert(reacA != NULL);

  add_reactants_to_product_list(rx, reacA, reacB, NULL, product, product_type);

  /* Determine whether any of the reactants can be replaced by a product.
     This is only useful for surface molecules to make sure reactions of the
     type A + B -> A don't push A to neighboring tiles */
  int replace_p1 = (product_type[0] == PLAYER_SURF_MOL && rx_players[0] == NULL);
  int replace_p2 = rx->n_reactants > 1 && (product_type[1] == PLAYER_SURF_MOL &&
                                           rx_players[1] == NULL);

  /* Determine the point of reaction on the surface. */
  struct vector2 rxn_uv_pos; // position of reaction on wall
  int rxn_uv_idx = -1;       // tile index of where reaction occurred
  int num_surface_static_reactants = 0; // number of reactants with (D_2D == 0)
  if (is_orientable) {
    if (sm_reactant) {
      rxn_uv_pos = sm_reactant->s_pos;
    } else {
      xyz2uv(hitpt, w, &rxn_uv_pos);
    }

    assert(w != NULL);
    if (w->grid == NULL) {
      /* reacA must be a volume molecule, or this wall would have a grid already. */
      assert(!IS_SURF_MOL(reacA));
      if (create_grid(world, w, ((struct volume_molecule *)reacA)->subvol))
        mcell_allocfailed("Failed to create a grid for a wall.");
    }
    rxn_uv_idx = uv2grid(&rxn_uv_pos, w->grid);

    /* find out number of static surface reactants */
    if ((sm_1 != NULL) && (!distinguishable(sm_1->properties->D, 0, EPS_C))) {
      num_surface_static_reactants++;
    }
    if ((sm_2 != NULL) && (!distinguishable(sm_2->properties->D, 0, EPS_C))) {
      num_surface_static_reactants++;
    }
  }

  /* find out number of surface products */
  int num_surface_products = 0;
  int num_surface_static_products = 0; // number of products with (D_2D == 0)
  for (int n_product = rx->n_reactants; n_product < n_players; ++n_product) {
    if (rx_players[n_product] == NULL) {
      continue;
    }
    if (rx_players[n_product]->flags & ON_GRID) {
      num_surface_products++;
      if (!distinguishable(rx_players[n_product]->D, 0, EPS_C))
        num_surface_static_products++;
    }
  }

  int mol_idx = INT_MAX;
  /* If the reaction involves a surface, make sure there is room for each
   * product. */
  struct tile_neighbor *tile_nbr_head = NULL; // list of neighbor tiles
  int tile_nbr_list_length = 0;
  struct tile_neighbor *tile_vacant_nbr_head = NULL; // list of vacant neighbor tiles
  if (is_orientable) {
    if (num_surface_products > 0) {
      if (sm_reactant != NULL) {
        find_neighbor_tiles(world, sm_reactant, sm_reactant->grid,
                            sm_reactant->grid_index, 1, 0, &tile_nbr_head,
                            &tile_nbr_list_length);
      } else {
        find_neighbor_tiles(world, sm_reactant, w->grid, rxn_uv_idx, 1, 0,
                            &tile_nbr_head, &tile_nbr_list_length);
      }

      /* Create list of vacant tiles */
      for (struct tile_neighbor *tile_nbr = tile_nbr_head; tile_nbr != NULL;
           tile_nbr = tile_nbr->next) {
        struct surface_molecule_list *sm_list = tile_nbr->grid->sm_list[tile_nbr->idx]; 
        if (sm_list == NULL || sm_list->sm == NULL) {
          num_vacant_tiles++;
          push_tile_neighbor_to_list(&tile_vacant_nbr_head, tile_nbr->grid, tile_nbr->idx);
        }
      }
    }

    /* Can this reaction happen at all? */
    int num_recycled_tiles = 0;
    if (replace_p1 && replace_p2) {
      num_recycled_tiles = 2;
    } else if (replace_p1 || replace_p2) {
      num_recycled_tiles = 1;
    }
    if (num_surface_products > num_vacant_tiles + num_recycled_tiles) {
      return cleanup_and_block_rx(tile_nbr_head, tile_vacant_nbr_head);
    }

    /* set the orientations of the products. */
    for (int n_product = 0; n_product < n_players; ++n_product) {
      /* Skip NULL products */
      if (rx_players[n_product] == NULL) {
        continue;
      }

      int this_geometry = rx->geometries[i0 + n_product];
      int relative_orient = (this_geometry < 0) ? -1 : 1;
      this_geometry = abs(this_geometry);

      /* Geometry of 0 means "random orientation" */
      if (this_geometry == 0) {
        product_orient[n_product] = (rng_uint(world->rng) & 1) ? 1 : -1;
      } else {
        if (this_geometry > (int)rx->n_reactants) {
          product_orient[n_product] = relative_orient *
              product_orient[this_geometry - rx->n_reactants - 1];
        } else if (this_geometry == 1) {
          product_orient[n_product] = relative_orient * orientA;
        } else if (this_geometry == 2 && reacB != NULL) {
          product_orient[n_product] = relative_orient * orientB;
        } else {
          product_orient[n_product] = relative_orient * 1;
        }
      }

      /* If this is a reactant... */
      if (n_product < (int)rx->n_reactants) {
        /* If this is a surface molecule, we need to set its orientation. */
        if (rx_players[n_product]->flags & ON_GRID) {
          assert(IS_SURF_MOL(product[n_product]));
          struct surface_molecule *sm =
              (struct surface_molecule *)product[n_product];

          /* If the new orientation doesn't match the old, we've got some work
           * to do. */
          if (sm->orient != product_orient[n_product]) {

            /* We're about to update the molecule's orientation, so we will
             * first remove it from the counts in case we have any
             * orientation-sensitive counts. Then, we will update the
             * orientation. Finally, we will add the molecule back into the
             * counts in its new orientation. */

            /* Remove molecule from counts in old orientation, if mol is
             * counted. */
            if (product[n_product]->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) {
              count_region_from_scratch(world,
                                        product[n_product], /* molecule */
                                        NULL,               /* rxn pathway */
                                        -1,                 /* remove count */
                                        NULL, /* Location at which to count */
                                        w,    /* Wall on which this happened */
                                        t,    /* Time of occurrence */
                                        NULL);
            }

            /* Force check for the unimolecular reactions
               after changing orientation.
               There are two possible cases to be covered here:
               1) when (sm->t2) was previously set to FOREVER
               2) there may be two or more unimolecular reactions involving surface
                class that have different kinetics.
             */
            if (((sm->flags & ACT_REACT) != 0) &&
                ((sm->properties->flags & CAN_SURFWALL) != 0)) {
              sm->t2 = 0;
            }

            /* Set the molecule's orientation. */
            sm->orient = product_orient[n_product];

            /* Add molecule back to counts in new orientation, if mol is
             * counted. */
            if (product[n_product]->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) {
              count_region_from_scratch(world,
                                        product[n_product], /* molecule */
                                        NULL,               /* rxn pathway */
                                        1,                  /* add count */
                                        NULL, /* Location at which to count */
                                        w,    /* Wall on which this happened */
                                        t,    /* Time of occurrence */
                                        NULL);
            }
          }
        } else if (!is_unimol) { /* Otherwise, check if we've crossed the plane. */
          if (product[n_product] == initiator) {
            if (product_orient[n_product] != initiatorOrient)
              cross_wall = true;
          }
        }
      }
    }

    /* find out where to place surface products */
    /* Some special cases are listed below. */
    if (num_surface_products == 1) {
      if (is_unimol && replace_p1) {
        for (int n_product = rx->n_reactants; n_product < n_players; n_product++) {
          if (rx_players[n_product] == NULL ||
             (rx_players[n_product]->flags & NOT_FREE) == 0) {
            continue;
          }

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
            product_grid[n_product] = ((struct surface_molecule *)reacA)->grid;
            product_grid_idx[n_product] = ((struct surface_molecule *)reacA)->grid_index;
            replace_p1 = 0;
            break;
          }
        }
      } else if ((num_surface_static_reactants == 1) && (num_surface_static_products == 1)
          && (replace_p1 || replace_p2)) {
        /* the lonely static product always replaces lonely static reactant */
        for (int n_product = rx->n_reactants; n_product < n_players; n_product++) {
          if (rx_players[n_product] == NULL ||
              (rx_players[n_product]->flags & NOT_FREE) == 0 ||
              distinguishable(rx_players[n_product]->D, 0, EPS_C)) {
            continue;
          }

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            if (replace_p1 && (!distinguishable(reacA->properties->D, 0, EPS_C))) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
              product_grid[n_product] = ((struct surface_molecule *)reacA)->grid;
              product_grid_idx[n_product] = ((struct surface_molecule *)reacA)->grid_index;
              replace_p1 = 0;
              break;
            } else if (replace_p2 && (!distinguishable(reacB->properties->D, 0, EPS_C))) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
              product_grid[n_product] = ((struct surface_molecule *)reacB)->grid;
              product_grid_idx[n_product] = ((struct surface_molecule *)reacB)->grid_index;
              break;
            }
          }
        }
      } else if (replace_p1 && replace_p2) {
        /* if both reactants should be replaced and there is only one surface product
          here we make sure that the initiator molecule is replaced */
        for (int n_product = rx->n_reactants; n_product < n_players; n_product++) {
          if ((rx_players[n_product] == NULL) ||
              ((rx_players[n_product]->flags & NOT_FREE) == 0)) {
            continue;
          }

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            if (reacA == initiator) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
              product_grid[n_product] = ((struct surface_molecule *)reacA)->grid;
              product_grid_idx[n_product] = ((struct surface_molecule *)reacA)->grid_index;
              replace_p1 = 0;
            } else {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
              product_grid[n_product] = ((struct surface_molecule *)reacB)->grid;
              product_grid_idx[n_product] = ((struct surface_molecule *)reacB)->grid_index;
            }
            break;
          }
        }
      } else if (replace_p1) {
        for (int n_product = rx->n_reactants; n_product < n_players; n_product++) {
          if ((rx_players[n_product] == NULL) ||
              ((rx_players[n_product]->flags & NOT_FREE) == 0)) {
            continue;
          }

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
            product_grid[n_product] = ((struct surface_molecule *)reacA)->grid;
            product_grid_idx[n_product] = ((struct surface_molecule *)reacA)->grid_index;
            replace_p1 = 0;
            break;
          }
        }
      } else if (replace_p2) {
        for (int n_product = rx->n_reactants; n_product < n_players; n_product++) {
          if ((rx_players[n_product] == NULL) |
              ((rx_players[n_product]->flags & NOT_FREE) == 0)) {
            continue;
          }

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
            product_grid[n_product] = ((struct surface_molecule *)reacB)->grid;
            product_grid_idx[n_product] = ((struct surface_molecule *)reacB)->grid_index;
            break;
          }
        }
      }
    } else if (num_surface_products > 1) {
      /* more than one surface products */
      if (num_surface_static_reactants > 0) {
        bool replace_reacA = (!distinguishable(reacA->properties->D, 0, EPS_C)) && replace_p1;
        bool replace_reacB =
            (reacB == NULL) ? false : (!distinguishable(reacB->properties->D, 0, EPS_C)) && replace_p2;

        if (replace_reacA || replace_reacB) {
          int max_static_count = (num_surface_static_products < num_surface_static_reactants)
            ? num_surface_static_products : num_surface_static_reactants;

          int count = 0;
          while (count < max_static_count) {
            unsigned int rnd_num = rng_uint(world->rng) % n_players;
            /* pass reactants */
            if ((rnd_num < rx->n_reactants) || (rx_players[rnd_num] == NULL) ||
                ((rx_players[rnd_num]->flags & NOT_FREE) == 0) ||
                distinguishable(rx_players[rnd_num]->D, 0, EPS_C)) {
              continue;
            }
            if (product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET) {
              if (replace_reacA) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                product_grid[rnd_num] = ((struct surface_molecule *)reacA)->grid;
                product_grid_idx[rnd_num] = ((struct surface_molecule *)reacA)->grid_index;
                count++;
                replace_p1 = 0;
                replace_reacA = 0;
              } else if (replace_reacB) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                product_grid[rnd_num] = ((struct surface_molecule *)reacB)->grid;
                product_grid_idx[rnd_num] = ((struct surface_molecule *)reacB)->grid_index;
                count++;
                replace_p2 = 0;
                replace_reacB = 0;
              }
            }
          } /* end while */
        }
      }

      /* check whether there are any surface reactants left to replace with
       * surface products since we are done with static reactants/products */
      if (replace_p1 || replace_p2) {
        int surf_prod_left = 0, surf_reactant_left = 0;
        for (int n_product = rx->n_reactants; n_product < n_players; n_product++) {
          if ((rx_players[n_product] == NULL) ||
              ((rx_players[n_product]->flags & NOT_FREE) == 0)) {
            continue;
          }
          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            surf_prod_left++;
          }
        }
        if (replace_p1) {
          surf_reactant_left++;
        }
        if (replace_p2) {
          surf_reactant_left++;
        }

        if (surf_prod_left > 0) {
          int num_to_place = surf_prod_left;
          if (surf_prod_left >= surf_reactant_left) {
            num_to_place = surf_reactant_left;
          }
          int count = 0;
          while (count < num_to_place) {
            unsigned int rnd_num = rng_uint(world->rng) % n_players;
            if ((rnd_num < rx->n_reactants) || (rx_players[rnd_num] == NULL) ||
                (rx_players[rnd_num]->flags & NOT_FREE) == 0) {
             continue;
            }
            if (product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET) {
              if (replace_p1) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                product_grid[rnd_num] = ((struct surface_molecule *)reacA)->grid;
                product_grid_idx[rnd_num] = ((struct surface_molecule *)reacA)->grid_index;
                count++;
                replace_p1 = 0;
              } else if (replace_p2) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                product_grid[rnd_num] = ((struct surface_molecule *)reacB)->grid;
                product_grid_idx[rnd_num] = ((struct surface_molecule *)reacB)->grid_index;
                replace_p2 = 0;
                count++;
              }
            }
          } /* end while */
        }
      }
    }

    /* here we will find placement for the case of the reaction of type
     * "vol_mol + w -> surf_mol + ...[rate] " */
    if ((sm_reactant == NULL) && (w != NULL) && (num_surface_products >= 1)) {
      assert(!IS_SURF_MOL(reacA));
      assert(rxn_uv_idx != -1);

      while (true) {
        unsigned int rnd_num = rng_uint(world->rng) % (n_players);
        if (rnd_num <= 1 || (rx_players[rnd_num] == NULL) ||
            (rx_players[rnd_num]->flags & NOT_FREE) == 0) {
          continue;
        }

        if (product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET) {
          product_flag[rnd_num] = PRODUCT_FLAG_USE_UV_LOC;
          product_grid[rnd_num] = w->grid;
          product_grid_idx[rnd_num] = rxn_uv_idx;
          break;
        }
      }
    }

    /* we will implement special placement policy for reaction of type of *
     * A->B+C[rate] */
    if (is_unimol && (sm_reactant != NULL) && (num_surface_products == 2)) {
      reac_idx = sm_reactant->grid_index;
      reac_grid = sm_reactant->grid;
    }

    // all other products are placed on one of the randomly chosen vacant tiles 
    int do_it_once = 0; /* flag */
    int num_attempts = 0;
    for (int n_product = rx->n_reactants; n_product < n_players; ++n_product) {
      /* If the product is a volume product, no placement is required. */
      if (rx_players[n_product]->flags & ON_GRID) {
        if (product_flag[n_product] != PRODUCT_FLAG_NOT_SET) {
          continue;
        }

        /* can't place products - reaction blocked */
        if (num_vacant_tiles == 0) {
          return cleanup_and_block_rx(tile_nbr_head, tile_vacant_nbr_head);
        }

        num_attempts = 0;
        while (true) {
          if (num_attempts > SURFACE_DIFFUSION_RETRIES) {
            return cleanup_and_block_rx(tile_nbr_head, tile_vacant_nbr_head);
          }

          /* randomly pick a tile from the list */
          unsigned int rnd_num = rng_uint(world->rng) % num_vacant_tiles;
          int tile_idx = -1; /* index of the tile on the grid */
          tile_grid = NULL;
          if (get_tile_neighbor_from_list_of_vacant_neighbors(
                  tile_vacant_nbr_head, rnd_num, &tile_grid, &tile_idx) == 0) {
            return cleanup_and_block_rx(tile_nbr_head, tile_vacant_nbr_head);
          }
          if (tile_idx < 0) {
            continue; /* this tile was probed already */
          }
          assert(tile_grid != NULL);

          /* make sure we can get to the tile given the surface regions defined
           * in the model */
          if (!product_tile_can_be_reached(tile_grid->surface, rlp_head_wall_1,
            rlp_head_wall_2, rlp_head_obj_1, rlp_head_obj_2, sm_bitmask, is_unimol)) {
            uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
            num_attempts++;
            continue;
          }

          product_grid[n_product] = tile_grid;
          product_grid_idx[n_product] = tile_idx;
          product_flag[n_product] = PRODUCT_FLAG_USE_RANDOM;
          if (!do_it_once && is_unimol && (sm_reactant != NULL) && (num_surface_products == 2)) {
            /*remember this tile (used for the reaction A->B+C[rate]) */
            mol_idx = tile_idx;
            mol_grid = tile_grid;
            do_it_once = 1;
          }
          break;
        } /* end while */
      }
    }
  } /* end if (is_orientable) */

  /* Determine the location of the reaction for count purposes. */
  struct vector3 count_pos_xyz;
  struct periodic_image *periodic_box = ((struct volume_molecule *)reacA)->periodic_box;
  if (hitpt != NULL) {
    count_pos_xyz = *hitpt;
  } else if (sm_reactant) {
    uv2xyz(&sm_reactant->s_pos, sm_reactant->grid->surface, &count_pos_xyz);
  } else {
    count_pos_xyz = ((struct volume_molecule *)reacA)->pos;
  }

  /* Create and place each product. */
  struct vector3 mol_pos_tmp;
  struct subvolume *product_subvol = NULL;
  /* Do we need to advance the dissociation index? */
  bool update_dissociation_index = false; 
  for (int n_product = rx->n_reactants; n_product < n_players; ++n_product) {
    struct abstract_molecule *this_product = NULL;
    struct species *const product_species = rx_players[n_product];

    /* If the product is a surface molecule, place it on the grid. */
    if (product_species->flags & ON_GRID) {
      struct vector2 prod_uv_pos;

      /* Pick an appropriate position for the new molecule. */
      if (world->randomize_smol_pos) {
        switch (product_flag[n_product]) {
        case PRODUCT_FLAG_USE_REACA_UV:
          if (is_unimol && (num_surface_products == 2) && (sm_reactant != NULL)) {
            if (mol_grid == NULL) {
              mcell_internal_error("Error in surface product placement for the "
                                   "unimolecular reaction.");
            }
            find_closest_position(product_grid[n_product], product_grid_idx[n_product],
              mol_grid, mol_idx, &prod_uv_pos);
          } else {
            prod_uv_pos = ((struct surface_molecule *)reacA)->s_pos;
          }
          break;

        case PRODUCT_FLAG_USE_REACB_UV:
          assert(reacB != NULL);
          prod_uv_pos = ((struct surface_molecule *)reacB)->s_pos;
          break;

        case PRODUCT_FLAG_USE_UV_LOC:
          prod_uv_pos = rxn_uv_pos;
          break;

        case PRODUCT_FLAG_USE_RANDOM:
          if (is_unimol && replace_p1 && (num_surface_products == 2)) {
            find_closest_position(product_grid[n_product],
                                  product_grid_idx[n_product], reac_grid,
                                  reac_idx, &prod_uv_pos);
          } else {
            grid2uv_random(product_grid[n_product], product_grid_idx[n_product],
                           &prod_uv_pos, world->rng);
          }
          break;

        default:
          UNHANDLED_CASE(product_flag[n_product]);
          /*break;*/
        }
      } else {
        grid2uv(product_grid[n_product], product_grid_idx[n_product], &prod_uv_pos);
      }

      this_product = (struct abstract_molecule *)place_sm_product(
          world, product_species, product_grid[n_product],
          product_grid_idx[n_product], &prod_uv_pos, product_orient[n_product],
          t, reacA->periodic_box);
    } else { /* else place the molecule in space. */
      /* For either a unimolecular reaction, or a reaction between two surface
         molecules we don't have a hitpoint. */
      if (!hitpt) {
        if (reacA->properties->flags & ON_GRID) {
          /* Since we use reactA's position to compute the location of the reaction
             we also need to use its wall for picking the displacement later in
             place_volume_products */
          w = ((struct surface_molecule *)reacA)->grid->surface;
          uv2xyz(&((struct surface_molecule *)reacA)->s_pos,
                 w, &mol_pos_tmp);
          product_subvol = find_subvolume(world, &mol_pos_tmp, NULL);
        } else {
          mol_pos_tmp = ((struct volume_molecule *)reacA)->pos;
          product_subvol = ((struct volume_molecule *)reacA)->subvol;
        }
        hitpt = &mol_pos_tmp;
      } else if (product_subvol == NULL) {
        product_subvol = find_subvolume(world, hitpt, NULL);
      }

      this_product = (struct abstract_molecule *)place_volume_product(
          world, product_species, sm_reactant, w, product_subvol, hitpt,
          product_orient[n_product], t, reacA->periodic_box);

      if (((struct volume_molecule *)this_product)->index < DISSOCIATION_MAX)
        update_dissociation_index = true;
    }

    /* Update molecule counts */
    ++product_species->population;
    if (product_species->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
      count_region_from_scratch(world, this_product, NULL, 1, NULL, NULL, t, this_product->periodic_box);

    /* preserve molecule id if rxn is unimolecular with one product */
    if (is_unimol && (n_players == 1)) {
      this_product->id = reacA->id;
      world->current_mol_id--; /* give back id we used */
      continue;
    }
    /* preserve molecule id if rxn is surface rxn with one product */
    if ((n_players == 3) && product_type[1] == PLAYER_WALL) {
      this_product->id = reacA->id;
      world->current_mol_id--; /* give back id we used */
      continue;
    }
  }

  /* If necessary, update the dissociation index. */
  if (update_dissociation_index) {
    if (--world->dissociation_index < DISSOCIATION_MIN)
      world->dissociation_index = DISSOCIATION_MAX;
  }

  /* Handle events triggered off of named reactions */
  if (rx->info[path].pathname != NULL) {
    /* No flags for reactions so we have to check regions if we have waypoints!
     * Fix to be more efficient for WORLD-only counts? */
    if (world->place_waypoints_flag)
      count_region_from_scratch(world, NULL, rx->info[path].pathname, 1,
                                &count_pos_xyz, w, t, periodic_box);

    /* Other magical stuff.  For now, can only trigger releases. */
    if (rx->info[path].pathname->magic != NULL) {
      if (reaction_wizardry(world, rx->info[path].pathname->magic, w,
                            &count_pos_xyz, t))
        mcell_allocfailed("Failed to complete reaction triggered release after "
                          "a '%s' reaction.",
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
outcome_unimolecular:
  In: world: simulation state
      rx: the reaction that is occuring
      path: the path that the reaction is taking
      reac: the molecule that is taking that path
      t: time that the reaction is occurring
  Out: Value based on outcome:
       RX_BLOCKED if there was no room to put products on grid
       RX_DESTROY if molecule no longer exists.
       RX_A_OK if it does.
       Products are created as needed.
*************************************************************************/
int outcome_unimolecular(struct volume *world, struct rxn *rx, int path,
                         struct abstract_molecule *reac, double t) {
  struct species *who_was_i = reac->properties;
  int result = RX_A_OK;
  struct volume_molecule *vm = NULL;
  struct surface_molecule *sm = NULL;

  if ((reac->properties->flags & NOT_FREE) == 0) {
    vm = (struct volume_molecule *)reac;
    result = outcome_products_random(world, NULL, NULL, t, rx, path, reac,
                                     NULL, 0, 0);
  } else {
    sm = (struct surface_molecule *)reac;
      /* we will not create products if the reaction is with an ABSORPTIVE
         region border */

    if ((strcmp(rx->players[0]->sym->name, "ALL_SURFACE_MOLECULES") == 0) ||
        (strcmp(rx->players[0]->sym->name, "ALL_MOLECULES") == 0)) {
      /* do nothing */
    } else {
      result = outcome_products_random(world, sm->grid->surface, NULL, t, rx,
                                       path, reac, NULL, sm->orient, 0);
    }
  }

  if (result == RX_BLOCKED)
    return RX_BLOCKED;

  if (result != RX_BLOCKED) {
    rx->info[path].count++;
    rx->n_occurred++;
  }

  struct species *who_am_i = rx->players[rx->product_idx[path]];

  if (who_am_i == NULL) {
    if (vm != NULL) {
      vm->subvol->mol_count--;
      if (vm->flags & IN_SCHEDULE)
        vm->subvol->local_storage->timer->defunct_count++;
      if (vm->properties->flags & COUNT_SOME_MASK) {
        count_region_from_scratch(world, (struct abstract_molecule *)vm, NULL,
                                  -1, &(vm->pos), NULL, vm->t, vm->periodic_box);
      }
    } else {
      remove_surfmol_from_list(&sm->grid->sm_list[sm->grid_index], sm);
      sm->grid->n_occupied--;
      if (sm->flags & IN_SCHEDULE) {
        sm->grid->subvol->local_storage->timer->defunct_count++;
      }
      if (sm->properties->flags & COUNT_SOME_MASK) {
        count_region_from_scratch(world, (struct abstract_molecule *)sm, NULL,
                                  -1, NULL, NULL, sm->t, sm->periodic_box);
      }
    }

    free(reac->periodic_box);
    who_was_i->n_deceased++;
    double t_time = convert_iterations_to_seconds(
        world->start_iterations, world->time_unit,
        world->simulation_start_seconds, t);
    who_was_i->cum_lifetime_seconds += t_time - reac->birthday;

    who_was_i->population--;
    if (vm != NULL)
      collect_molecule(vm);
    else {
      reac->properties = NULL;
      mem_put(reac->birthplace, reac);
    }
    return RX_DESTROY;
  } else if (who_am_i != who_was_i) {
    if (vm != NULL)
      collect_molecule(vm);
    else
      reac->properties = NULL;
    return RX_DESTROY;
  } else
    return result;
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
int outcome_bimolecular(struct volume *world, struct rxn *rx, int path,
                        struct abstract_molecule *reacA,
                        struct abstract_molecule *reacB, short orientA,
                        short orientB, double t, struct vector3 *hitpt,
                        struct vector3 *loc_okay) {

  assert(periodic_boxes_are_identical(reacA->periodic_box, reacB->periodic_box));

  struct surface_molecule *sm = NULL;
  struct volume_molecule *vm = NULL;
  struct wall *w = NULL;
  /* struct storage *x; */
  int result;
  int reacB_was_free = 0;
  int killA, killB;

  if ((reacA->properties->flags & NOT_FREE) == 0) {
    if ((reacB->properties->flags & ON_GRID) != 0) {
      sm = (struct surface_molecule *)reacB;
      w = sm->grid->surface;
    }
  } else { /* Surface molecule */
    sm = (struct surface_molecule *)reacA;
    w = sm->grid->surface;
  }

  result = outcome_products_random(world, w, hitpt, t, rx, path, reacA, reacB,
                                   orientA, orientB);

  if (result == RX_BLOCKED)
    return RX_BLOCKED;

  rx->n_occurred++;
  rx->info[path].count++;

  /* Figure out if either of the reactants was destroyed */
  if (rx->players[0] == reacA->properties) {
    killB = (rx->players[rx->product_idx[path] + 1] == NULL);
    killA = (rx->players[rx->product_idx[path]] == NULL);
  } else {
    killB = (rx->players[rx->product_idx[path]] == NULL);
    killA = (rx->players[rx->product_idx[path] + 1] == NULL);
  }

  if (killB) {
    vm = NULL;
    if ((reacB->properties->flags & ON_GRID) != 0) {
      sm = (struct surface_molecule *)reacB;

      remove_surfmol_from_list(&sm->grid->sm_list[sm->grid_index], sm);
      sm->grid->n_occupied--;
      if (sm->flags & IN_SURFACE)
        sm->flags -= IN_SURFACE;
      if (sm->flags & IN_SCHEDULE) {
        sm->grid->subvol->local_storage->timer->defunct_count++;
      }
    } else if ((reacB->properties->flags & NOT_FREE) == 0) {
      vm = (struct volume_molecule *)reacB;
      vm->subvol->mol_count--;
      if (vm->flags & IN_SCHEDULE) {
        vm->subvol->local_storage->timer->defunct_count++;
      }
      reacB_was_free = 1;
    }

    if ((reacB->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) != 0) {
      count_region_from_scratch(world, reacB, NULL, -1, NULL, NULL, t, reacB->periodic_box);
    }

    free(reacB->periodic_box);
    reacB->properties->n_deceased++;
    double t_time = convert_iterations_to_seconds(
        world->start_iterations, world->time_unit,
        world->simulation_start_seconds, t);
    reacB->properties->cum_lifetime_seconds += t_time - reacB->birthday;
    reacB->properties->population--;

    if (vm != NULL)
      collect_molecule(vm);
    else
      reacB->properties = NULL;
  }

  if (killA) {
    vm = NULL;
    if ((reacA->properties->flags & ON_GRID) != 0) {
      sm = (struct surface_molecule *)reacA;

      remove_surfmol_from_list(&sm->grid->sm_list[sm->grid_index], sm);
      sm->grid->n_occupied--;
      if (sm->flags & IN_SCHEDULE) {
        sm->grid->subvol->local_storage->timer->defunct_count++;
      }
    } else if ((reacA->properties->flags & NOT_FREE) == 0) {
      vm = (struct volume_molecule *)reacA;
      vm->subvol->mol_count--;
      if (vm->flags & IN_SCHEDULE) {
        vm->subvol->local_storage->timer->defunct_count++;
      }
    }

    if ((reacA->properties->flags & ON_GRID) !=
        0) /* Surface molecule is OK where it is, doesn't obey COUNT_ME */
    {
      if (reacA->properties->flags &
          COUNT_SOME_MASK) /* If we're ever counted, try to count us now */
      {
        count_region_from_scratch(world, reacA, NULL, -1, NULL, NULL, t, reacA->periodic_box);
      }
    } else if (reacA->flags & COUNT_ME) {
      /* Subtlety: we made it up to hitpt, but our position is wherever we were
       * before that! */
      if (hitpt == NULL || reacB_was_free ||
          (reacB->properties != NULL &&
           (reacB->properties->flags & NOT_FREE) == 0)) {
        /* Vol-vol rx should be counted at hitpt */
        count_region_from_scratch(world, reacA, NULL, -1, hitpt, NULL, t, reacA->periodic_box);
      } else /* Vol-surf but don't want to count exactly on a wall or we might
                count on the wrong side */
      {
        struct vector3 fake_hitpt;

        vm = (struct volume_molecule *)reacA;

        /* Halfway in between where we were and where we react should be a safe
         * away-from-wall place to remove us */
        if (loc_okay == NULL)
          loc_okay = &(vm->pos);
        fake_hitpt.x = 0.5 * hitpt->x + 0.5 * loc_okay->x;
        fake_hitpt.y = 0.5 * hitpt->y + 0.5 * loc_okay->y;
        fake_hitpt.z = 0.5 * hitpt->z + 0.5 * loc_okay->z;

        count_region_from_scratch(world, reacA, NULL, -1, &fake_hitpt, NULL, t, reacA->periodic_box);
      }
    }

    free(reacA->periodic_box);
    reacA->properties->n_deceased++;
    double t_time = convert_iterations_to_seconds(
        world->start_iterations, world->time_unit,
        world->simulation_start_seconds, t);
    reacA->properties->cum_lifetime_seconds += t_time - reacA->birthday;
    reacA->properties->population--;

    if (vm != NULL)
      collect_molecule(vm);
    else
      reacA->properties = NULL;

    return RX_DESTROY;
  }

  return result;
}

/*************************************************************************
outcome_intersect:
  In: world: simulation state
      rx: reaction that's taking place
      path: path the reaction's taking
      surface: wall that is being struck
      reac: molecule that is hitting the wall
      orient: orientation of the molecule
      t: time that the reaction is occurring
      hitpt: location of collision with wall
      loc_okay:
  Out: Value depending on outcome:
       RX_A_OK if the molecule reflects
       RX_FLIP if the molecule passes through
       RX_DESTROY if the molecule stops, is destroyed, etc.
       Additionally, products are created as needed.
  Note: Can assume molecule is always first in the reaction.
*************************************************************************/
int outcome_intersect(struct volume *world, struct rxn *rx, int path,
                      struct wall *surface, struct abstract_molecule *reac,
                      short orient, double t, struct vector3 *hitpt,
                      struct vector3 *loc_okay) {

  if (rx->n_pathways <= RX_SPECIAL) {
    rx->n_occurred++;
    if (rx->n_pathways == RX_REFLEC)
      return RX_A_OK;
    else
      return RX_FLIP; /* Flip = transparent is default special case */
  }
  int idx = rx->product_idx[path];

  if ((reac->properties->flags & NOT_FREE) == 0) {
    struct volume_molecule *vm = (struct volume_molecule *)reac;

    /* If reaction object has ALL_MOLECULES or ALL_VOLUME_MOLECULES as the
     * first reactant it means that reaction is of the type ABSORPTIVE =
     * ALL_MOLECULES or ABSORPTIVE = ALL_VOLUME_MOLECULES since other cases
     * (REFLECTIVE/TRANSPARENT) are taken care above. But there are no products
     * for this reaction, so we do no need to go into "outcome_products()"
     * function. */

    int result;
    if ((strcmp(rx->players[0]->sym->name, "ALL_MOLECULES") == 0) ||
        (strcmp(rx->players[0]->sym->name, "ALL_VOLUME_MOLECULES") == 0)) {
      result = RX_DESTROY;
    } else {
      result = outcome_products_random(world, surface, hitpt, t, rx, path,
                                       reac, NULL, orient, 0);
    }
    if (result == RX_BLOCKED)
      return RX_A_OK; /* reflect the molecule */

    rx->info[path].count++;
    rx->n_occurred++;

    if (rx->players[idx] == NULL) {
      /* The code below is also valid for the special reaction of the type
       * ABSORPTIVE = ALL_MOLECULES (or ALL_VOLUME_MOLECULES) */
      vm->subvol->mol_count--;
      if (world->place_waypoints_flag && (reac->flags & COUNT_ME)) {
        if (hitpt == NULL) {
          count_region_from_scratch(
            world, reac, NULL, -1, NULL, NULL, t, reac->periodic_box);
        } else {
          struct vector3 fake_hitpt;

          /* Halfway in between where we were and where we react should be a
           * safe away-from-wall place to remove us */
          if (loc_okay == NULL)
            loc_okay = &(vm->pos);
          fake_hitpt.x = 0.5 * hitpt->x + 0.5 * loc_okay->x;
          fake_hitpt.y = 0.5 * hitpt->y + 0.5 * loc_okay->y;
          fake_hitpt.z = 0.5 * hitpt->z + 0.5 * loc_okay->z;

          count_region_from_scratch(world, reac, NULL, -1, &fake_hitpt, NULL,
                                    t, reac->periodic_box);
        }
      }
      free(reac->periodic_box);
      reac->properties->n_deceased++;
      double t_time = convert_iterations_to_seconds(
          world->start_iterations, world->time_unit,
          world->simulation_start_seconds, t);
      reac->properties->cum_lifetime_seconds += t_time - reac->birthday;
      reac->properties->population--;
      if (vm->flags & IN_SCHEDULE) {
        vm->subvol->local_storage->timer->defunct_count++;
      }
      collect_molecule(vm);
      return RX_DESTROY;
    } else
      return result; /* RX_A_OK or RX_FLIP */
  } else {
    /* Should really be an error because we should never call
     * outcome_intersect() on a surface molecule */
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
int reaction_wizardry(struct volume *world, struct magic_list *incantation,
                      struct wall *surface, struct vector3 *hitpt, double t) {
  struct release_event_queue req; /* Create a release event on the fly */

  /* Release event happens "now" */
  req.next = NULL;
  req.event_time = t;
  req.train_counter = 0;
  req.train_high_time = t;

  /* Set up transform to place products at site of reaction */
  if (hitpt == NULL) {
    init_matrix(req.t_matrix);
  } else if (surface == NULL ||
             !distinguishable(surface->normal.z, 1.0,
                              EPS_C)) /* Just need a translation */
  {
    init_matrix(req.t_matrix);
    req.t_matrix[3][0] = hitpt->x;
    req.t_matrix[3][1] = hitpt->y;
    req.t_matrix[3][2] = hitpt->z;
  } else /* Set up transform that will translate and then rotate Z axis to align
            with surface normal */
  {
    struct vector3 scale = { 1.0, 1.0, 1.0 }; /* No scaling */
    struct vector3 axis = { 1.0, 0.0, 0.0 };  /* X-axis is default */
    double cos_theta;
    double degrees;

    cos_theta = surface->normal.z; /* (0,0,1) . surface->normal */
    if (!distinguishable(cos_theta, -1.0, EPS_C)) {
      degrees = 180.0; /* Upside-down */
    } else {
      /* (0,0,1) x surface->normal */
      axis.x = -surface->normal.y;
      axis.y = surface->normal.x;
      axis.z = 0.0;

      degrees = acos(cos_theta) * 180.0 / MY_PI;
    }
    tform_matrix(&scale, hitpt, &axis, degrees, req.t_matrix);
  }

  /* Now we're ready to cast our spell! */
  for (; incantation != NULL; incantation = incantation->next) {
    if (incantation->type != magic_release)
      continue; /* Only know how to magically release stuff */

    req.release_site = (struct release_site_obj *)incantation->data;

    if (release_molecules(world, &req))
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
 * in: surface molecule 1 (located on wall 1)
 *     surface molecule 2 (located on wall 2)
 *     pointer to array with restrictive regions which contain wall 1
 *     pointer to array with restrictive regions which contain wall 2
 *     pointer to array with restrictive regions which don't contain wall 1
 *     pointer to array with restrictive regions which don't contain wall 2
 *
 * out: the 4 arrays with pointers to restrictive regions will be filled
 *      and returned
 *
 ***********************************************************************/
int determine_molecule_region_topology(
    struct volume *world, struct surface_molecule *sm_1,
    struct surface_molecule *sm_2, struct region_list **rlp_wall_1_ptr,
    struct region_list **rlp_wall_2_ptr, struct region_list **rlp_obj_1_ptr,
    struct region_list **rlp_obj_2_ptr, bool is_unimol) {
  int sm_bitmask = 0;
  struct wall *w_1, *w_2;
  struct region_list *rlp_head_wall_1 = NULL;
  struct region_list *rlp_head_wall_2 = NULL;
  struct region_list *rlp_head_obj_1 = NULL;
  struct region_list *rlp_head_obj_2 = NULL;

  /* bimolecular reactions */
  if ((sm_1 != NULL) && (sm_2 != NULL)) {
    /* both reactants have restrictive region borders */
    if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
        (sm_2->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(
            world, sm_1->grid->surface->parent_object, sm_1) &&
        are_restricted_regions_for_species_on_object(
            world, sm_2->grid->surface->parent_object, sm_2)) {
      w_1 = sm_1->grid->surface;
      w_2 = sm_2->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      rlp_head_wall_2 = find_restricted_regions_by_wall(world, w_2, sm_2);

      /* both reactants are inside their respective restricted regions */
      if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL)) {
        sm_bitmask |= ALL_INSIDE;
      }
      /* both reactants are outside their respective restricted regions */
      else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL)) {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        sm_bitmask |= ALL_OUTSIDE;
      }
      /* grid1 is inside and grid2 is outside of its respective
       * restrictive region */
      else if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 == NULL)) {
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        sm_bitmask |= SURF1_IN_SURF2_OUT;
      }
      /* grid2 is inside and grid1 is outside of its respective
       * restrictive region */
      else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL)) {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        sm_bitmask |= SURF1_OUT_SURF2_IN;
      }
    }

    /* only reactant sm_1 has restrictive region border property */
    else if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
             are_restricted_regions_for_species_on_object(
                 world, sm_1->grid->surface->parent_object, sm_1) &&
             (!(sm_2->properties->flags & CAN_REGION_BORDER) ||
              !are_restricted_regions_for_species_on_object(
                   world, sm_2->grid->surface->parent_object, sm_2))) {
      w_1 = sm_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      if (rlp_head_wall_1 != NULL) {
        sm_bitmask |= SURF1_IN;
      } else {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        sm_bitmask |= SURF1_OUT;
      }
    }

    /* only reactant "sm_2" has restrictive region border property */
    else if ((sm_2->properties->flags & CAN_REGION_BORDER) &&
             are_restricted_regions_for_species_on_object(
                 world, sm_2->grid->surface->parent_object, sm_2) &&
             (!(sm_1->properties->flags & CAN_REGION_BORDER) ||
              !are_restricted_regions_for_species_on_object(
                   world, sm_1->grid->surface->parent_object, sm_1))) {
      w_2 = sm_2->grid->surface;
      rlp_head_wall_2 = find_restricted_regions_by_wall(world, w_2, sm_2);
      if (rlp_head_wall_2 != NULL) {
        sm_bitmask |= SURF2_IN;
      } else {
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        sm_bitmask |= SURF2_OUT;
      }
    }
  }

  /* unimolecular reactions */
  else if ((sm_1 != NULL) && is_unimol) {
    if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(
            world, sm_1->grid->surface->parent_object, sm_1)) {
      w_1 = sm_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      if (rlp_head_wall_1 != NULL) {
        sm_bitmask |= ALL_INSIDE;
      } else {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        sm_bitmask |= ALL_OUTSIDE;
      }
    }
  }

  *rlp_wall_1_ptr = rlp_head_wall_1;
  *rlp_wall_2_ptr = rlp_head_wall_2;
  *rlp_obj_1_ptr = rlp_head_obj_1;
  *rlp_obj_2_ptr = rlp_head_obj_2;

  return sm_bitmask;
}

/***********************************************************************
 *
 * this function tests if wall target can be reached for product placement
 * based on the previously stored reactant topology based on
 * sm_bitmask. Below, wall 1 is the wall containing reactant 1 and
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
                                 int sm_bitmask, bool is_unimol) {
  bool status = true;

  if (sm_bitmask & ALL_INSIDE) {
    if (is_unimol) {
      if (!wall_belongs_to_all_regions_in_region_list(target,
                                                      rlp_head_wall_1)) {
        status = false;
      }
    } else {
      /* bimol reaction */
      if (!wall_belongs_to_all_regions_in_region_list(target,
                                                      rlp_head_wall_1) ||
          !wall_belongs_to_all_regions_in_region_list(target,
                                                      rlp_head_wall_2)) {
        status = false;
      }
    }
  } else if (sm_bitmask & ALL_OUTSIDE) {
    if (is_unimol) {
      if (wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_1)) {
        status = false;
      }
    } else {
      if (wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_1) ||
          wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_2)) {
        status = false;
      }
    }
  } else if (sm_bitmask & SURF1_IN_SURF2_OUT) {
    if (!wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_1) ||
        wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_2)) {
      status = false;
    }
  } else if (sm_bitmask & SURF1_OUT_SURF2_IN) {
    if (wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_1) ||
        !wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_2)) {
      status = false;
    }
  } else if (sm_bitmask & SURF1_IN) {
    if (!wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_1)) {
      status = false;
    }
  } else if (sm_bitmask & SURF1_OUT) {
    if (wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_1)) {
      status = false;
    }
  } else if (sm_bitmask & SURF2_IN) {
    if (!wall_belongs_to_all_regions_in_region_list(target, rlp_head_wall_2)) {
      status = false;
    }
  } else if (sm_bitmask & SURF2_OUT) {
    if (wall_belongs_to_any_region_in_region_list(target, rlp_head_obj_2)) {
      status = false;
    }
  }

  return status;
}

/*
 * cleanup_and_block_rx is a simple helper function which deletes the provided
 * linked lists of tile_neighbors and the return RX_BLOCKED
 */
int cleanup_and_block_rx(struct tile_neighbor *tn1, struct tile_neighbor *tn2) {
  if (tn1 != NULL) {
    delete_tile_neighbor_list(tn1);
  }
  if (tn2 != NULL) {
    delete_tile_neighbor_list(tn2);
  }
  return RX_BLOCKED;
}
