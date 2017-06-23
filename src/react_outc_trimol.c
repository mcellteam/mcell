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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "logging.h"
#include "rng.h"
#include "util.h"
#include "grid_util.h"
#include "mcell_structs.h"
#include "count_util.h"
#include "react.h"
#include "vol_util.h"
#include "wall_util.h"

static int outcome_products_trimol_reaction_random(
    struct volume *world, struct wall *w, struct vector3 *hitpt, double t,
    struct rxn *rx, int path, struct abstract_molecule *reacA,
    struct abstract_molecule *reacB, struct abstract_molecule *reacC,
    short orientA, short orientB, short orientC);

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
static int outcome_products_trimol_reaction_random(
    struct volume *world, struct wall *w, struct vector3 *hitpt, double t,
    struct rxn *rx, int path, struct abstract_molecule *reacA,
    struct abstract_molecule *reacB, struct abstract_molecule *reacC,
    short orientA, short orientB, short orientC) {

  if (reacA != NULL && reacB != NULL) {
    assert(periodic_boxes_are_identical(reacA->periodic_box, reacB->periodic_box));
  } else if (reacA != NULL && reacC != NULL) {
    assert(periodic_boxes_are_identical(reacA->periodic_box, reacC->periodic_box));
  } else if (reacB != NULL && reacC != NULL) {
    assert(periodic_boxes_are_identical(reacB->periodic_box, reacC->periodic_box));
  }

  bool update_dissociation_index =
      false;               /* Do we need to advance the dissociation index? */
  bool cross_wall = false; /* Did the moving molecule cross the plane? */
  struct subvolume *last_subvol =
      NULL; /* Last subvolume (guess used to speed sv finding) */

  int const i0 =
      rx->product_idx[path]; /* index of the first player for the pathway */
  int const iN =
      rx->product_idx[path + 1]; /* index of the first player for the next
                                    pathway */
  assert(iN > i0);
  struct species **rx_players =
      rx->players + i0; /* Players array from the reaction. */

  int const n_players = iN - i0;                /* number of reaction players */
  struct abstract_molecule *product[n_players]; /* array of products */
  char product_type[n_players]; /* array that decodes the type of each product
                                   */
  short product_orient[n_players]; /* array of orientations for each product */
  struct surface_grid *product_grid[n_players]; /* array of surface_grids for
                                                   products */
  int product_grid_idx[n_players]; /* array of grid indices for products */
  byte product_flag[n_players];    /* array of placement flags for products */

  struct tile_neighbor *tile_nbr_head = NULL; /* list of neighbor tiles */
  struct tile_neighbor *tile_nbr;             /* iterator */
  /* head of the linked list of vacant neighbor tiles */
  struct tile_neighbor *tile_vacant_nbr_head = NULL;
  struct surface_grid *tile_grid; /* surface grid the tile belongs to */
  int tile_idx;                   /* index of the tile on the grid */
  unsigned int rnd_num;           /* random number */
  int num_vacant_tiles = 0;       /* number of vacant tiles */
  int num_surface_products = 0;   /* not counting reactants */
  int num_surface_static_products =
      0; /* # of products with (D_2D == 0) not counting reactants */
  int num_surface_static_reactants = 0; /* # of reactants with (D_2D == 0) */
  /* total number of surface reactants */
  int num_surface_reactants = 0;
  /* number of surface reactants that are not replaced in the reaction */
  int num_surface_reactants_to_stay = 0;
  int list_length; /* length of the linked list tile_nbr_head */
  /* flags */
  int replace_p1 = 0, replace_p2 = 0, replace_p3 = 0, only_one_to_replace = 0,
      two_to_replace = 0;
  int find_neighbor_tiles_flag = 0;
  struct wall *w_1, *w_2, *w_3;

  /* flag that indicates that all reactants lie inside their
     respective restrictive regions */
  int all_inside_restricted_boundary = 0;
  /* flag that indicates that all reactants lie outside
     their respective restrictive regions */
  int all_outside_restricted_boundary = 0;

  /* flag that indicates that  reactants "sm_1" and "sm_2"
     lies inside and reactant "grid_3" lies outside of their
     respective restrictive regions */
  int sm_1_inside_sm_2_inside_grid_3_outside = 0;
  /* flag that indicates that  reactants "sm_1" and "grid_3"
     lies inside and reactant "sm_2" lies outside of their
     respective restrictive regions */
  int sm_1_inside_sm_2_outside_grid_3_inside = 0;
  /* flag that indicates that  reactant "sm_1" lies outside
     and reactants "sm_2" and "grid_3" lie inside of their
     respective restrictive regions */
  int sm_1_outside_sm_2_inside_grid_3_inside = 0;
  /* flag that indicates that  reactant "sm_1" lies inside
     and reactants "sm_2" and "grid_3" lie outside of their
     respective restrictive regions */
  int sm_1_inside_sm_2_outside_grid_3_outside = 0;
  /* flag that indicates that  reactants "sm_1" and "grid_3" lie outside
     and reactant "sm_2" lies inside of their
     respective restrictive regions */
  int sm_1_outside_sm_2_inside_grid_3_outside = 0;
  /* flag that indicates that  reactants "sm_1" and "sm_2"
     lie outside and reactant "grid_3" lies inside of their
     respective restrictive regions */
  int sm_1_outside_sm_2_outside_grid_3_inside = 0;

  /* flag that indicates that  reactants "sm_1" and "sm_2"
     lie inside their respective restrictive regions
     and reactant "grid_3" has no restrictive region border properties */
  int only_sm_1_sm_2_inside = 0;
  /* flag that indicates that  reactant "sm_1" lies inside
     and "sm_2" lies outside their respective restrictive regions
     and reactant "grid_3" has no restrictive region border properties */
  int only_sm_1_inside_sm_2_outside = 0;
  /* flag that indicates that  reactant "sm_1" lies outside
     and "sm_2" lies inside their respective restrictive regions
     and reactant "grid_3" has no restrictive region border properties */
  int only_sm_1_outside_sm_2_inside = 0;
  /* flag that indicates that  reactants "sm_1" and "sm_2"
     lie outside their respective restrictive regions
     and reactant "grid_3" has no restrictive region border properties */
  int only_sm_1_sm_2_outside = 0;

  /* flag that indicates that  reactants "sm_1" and "grid_3"
     lie inside their respective restrictive regions
     and reactant "sm_2" has no restrictive region border properties */
  int only_sm_1_grid_3_inside = 0;
  /* flag that indicates that  reactant "sm_1" lies inside
     and "grid_3" lies outside their respective restrictive regions
     and reactant "sm_2" has no restrictive region border properties */
  int only_sm_1_inside_grid_3_outside = 0;
  /* flag that indicates that  reactant "sm_1" lies outside
     and "grid_3" lies inside their respective restrictive regions
     and reactant "sm_2" has no restrictive region border properties */
  int only_sm_1_outside_grid_3_inside = 0;
  /* flag that indicates that  reactants "sm_1" and "grid_3"
     lie outside their respective restrictive regions
     and reactant "sm_2" has no restrictive region border properties */
  int only_sm_1_grid_3_outside = 0;

  /* flag that indicates that  reactants "sm_2" and "grid_3"
     lie inside their respective restrictive regions
     and reactant "sm_1" has no restrictive region border properties */
  int only_sm_2_grid_3_inside = 0;
  /* flag that indicates that  reactant "sm_2" lies inside
     and "grid_3" lies outside their respective restrictive regions
     and reactant "sm_1" has no restrictive region border properties */
  int only_sm_2_inside_grid_3_outside = 0;
  /* flag that indicates that  reactant "sm_2" lies outside
     and "grid_3" lies inside their respective restrictive regions
     and reactant "sm_1" has no restrictive region border properties */
  int only_sm_2_outside_grid_3_inside = 0;
  /* flag that indicates that  reactants "sm_2" and "grid_3"
     lie outside their respective restrictive regions
     and reactant "sm_1" has no restrictive region border properties */
  int only_sm_2_grid_3_outside = 0;

  /* flag that indicates that only reactant "sm_1" has
     restrictive regions on the object and it lies
     inside it's restrictive region.  */
  int only_sm_1_inside = 0;
  /* flag that indicates that only reactant "sm_1" has
     restrictive regions on the object and it lies
     outside it's restrictive region.  */
  int only_sm_1_outside = 0;
  /* flag that indicates that only reactant "sm_2" has
     restrictive regions on the object and it lies
     inside it's restrictive region.  */
  int only_sm_2_inside = 0;
  /* flag that indicates that only reactant "sm_2" has
     restrictive regions on the object and it lies
     outside it's restrictive region.  */
  int only_sm_2_outside = 0;
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
  int rxn_uv_idx = -1;       /* tile index of the reaction place */

  struct abstract_molecule *tmp_mol;
  short tmp_orient;

  if ((reacA == NULL) || (reacB == NULL) || (reacC == NULL)) {
    mcell_internal_error("One of the reactants in "
                         "'outcome_products_trimol_reaction_random()' is "
                         "NULL.");
  }

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
  struct surface_molecule *const grid_3 =
      IS_SURF_MOL(reacC) ? (struct surface_molecule *)reacC : NULL;
  struct surface_molecule *sm_reactant = NULL;
  if (sm_1 != NULL) {
    sm_reactant = sm_1;
  } else if (sm_2 != NULL) {
    sm_reactant = sm_2;
  } else if (grid_3 != NULL) {
    sm_reactant = grid_3;
  }

  bool const is_orientable = (w != NULL) || (sm_reactant != NULL);

  /* Where are reactants relative to the restrictive region border? */
  if ((sm_1 != NULL) && (sm_2 != NULL) && (grid_3 != NULL)) {
    /* trimol_reaction */
    if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
        (sm_2->properties->flags & CAN_REGION_BORDER) &&
        (grid_3->properties->flags & CAN_REGION_BORDER) &&
        are_restricted_regions_for_species_on_object(
            world, sm_1->grid->surface->parent_object, sm_1) &&
        are_restricted_regions_for_species_on_object(
            world, sm_2->grid->surface->parent_object, sm_2) &&
        are_restricted_regions_for_species_on_object(
            world, grid_3->grid->surface->parent_object, grid_3)) {
      w_1 = sm_1->grid->surface;
      w_2 = sm_2->grid->surface;
      w_3 = grid_3->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      rlp_head_wall_2 = find_restricted_regions_by_wall(world, w_2, sm_2);
      rlp_head_wall_3 = find_restricted_regions_by_wall(world, w_3, grid_3);

      if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL) &&
          (rlp_head_wall_3 != NULL)) {
        /* all reactants are inside their respective
           restricted regions */
        all_inside_restricted_boundary = 1;

      } else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL) &&
                 (rlp_head_wall_3 == NULL)) {
        /* all reactants are outside their respective
           restricted regions */
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);

        all_outside_restricted_boundary = 1;
      } else if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL) &&
                 (rlp_head_wall_3 == NULL)) {
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);
        sm_1_inside_sm_2_inside_grid_3_outside = 1;
      } else if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_3 != NULL) &&
                 (rlp_head_wall_2 == NULL)) {
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        sm_1_inside_sm_2_outside_grid_3_inside = 1;
      } else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL) &&
                 (rlp_head_wall_3 == NULL)) {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);
        sm_1_outside_sm_2_inside_grid_3_outside = 1;
      } else if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 == NULL) &&
                 (rlp_head_wall_3 == NULL)) {
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);
        sm_1_inside_sm_2_outside_grid_3_outside = 1;
      } else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL) &&
                 (rlp_head_wall_3 != NULL)) {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        sm_1_outside_sm_2_inside_grid_3_inside = 1;
      } else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL) &&
                 (rlp_head_wall_3 != NULL)) {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        sm_1_outside_sm_2_outside_grid_3_inside = 1;
      }
    } else if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, sm_1->grid->surface->parent_object, sm_1) &&
               (sm_2->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, sm_2->grid->surface->parent_object, sm_2) &&
               (!(grid_3->properties->flags & CAN_REGION_BORDER) ||
                !are_restricted_regions_for_species_on_object(
                     world, grid_3->grid->surface->parent_object, grid_3))) {
      /* only reactants "sm_1" and "sm_2" have restrictive
         region border property */
      w_1 = sm_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      w_2 = sm_2->grid->surface;
      rlp_head_wall_2 = find_restricted_regions_by_wall(world, w_2, sm_2);
      if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 != NULL)) {
        only_sm_1_sm_2_inside = 1;
      } else if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_2 == NULL)) {
        only_sm_1_inside_sm_2_outside = 1;
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
      } else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 != NULL)) {
        only_sm_1_outside_sm_2_inside = 1;
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
      } else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_2 == NULL)) {
        only_sm_1_sm_2_outside = 1;
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
      }
    } else if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, sm_1->grid->surface->parent_object, sm_1) &&
               (!(sm_2->properties->flags & CAN_REGION_BORDER) ||
                !are_restricted_regions_for_species_on_object(
                     world, sm_2->grid->surface->parent_object, sm_2)) &&
               (grid_3->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, grid_3->grid->surface->parent_object, grid_3)) {
      /* only reactants "sm_1" and "grid_3" have restrictive
         region border property */
      w_1 = sm_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      w_3 = grid_3->grid->surface;
      rlp_head_wall_3 = find_restricted_regions_by_wall(world, w_3, grid_3);
      if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_3 != NULL)) {
        only_sm_1_grid_3_inside = 1;
      } else if ((rlp_head_wall_1 != NULL) && (rlp_head_wall_3 == NULL)) {
        only_sm_1_inside_grid_3_outside = 1;
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);
      } else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_3 != NULL)) {
        only_sm_1_outside_grid_3_inside = 1;
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
      } else if ((rlp_head_wall_1 == NULL) && (rlp_head_wall_3 == NULL)) {
        only_sm_1_grid_3_outside = 1;
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);
      }
    } else if ((!(sm_1->properties->flags & CAN_REGION_BORDER) ||
                (!are_restricted_regions_for_species_on_object(
                      world, sm_1->grid->surface->parent_object, sm_1))) &&
               (sm_2->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, sm_2->grid->surface->parent_object, sm_2) &&
               (grid_3->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, grid_3->grid->surface->parent_object, grid_3)) {
      /* only reactants "sm_2" and "grid_3" have restrictive
         region border property */
      w_2 = sm_2->grid->surface;
      rlp_head_wall_2 = find_restricted_regions_by_wall(world, w_2, sm_2);
      w_3 = grid_3->grid->surface;
      rlp_head_wall_3 = find_restricted_regions_by_wall(world, w_3, grid_3);
      if ((rlp_head_wall_2 != NULL) && (rlp_head_wall_3 != NULL)) {
        only_sm_2_grid_3_inside = 1;
      } else if ((rlp_head_wall_2 != NULL) && (rlp_head_wall_3 == NULL)) {
        only_sm_2_inside_grid_3_outside = 1;
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);
      } else if ((rlp_head_wall_2 == NULL) && (rlp_head_wall_3 != NULL)) {
        only_sm_2_outside_grid_3_inside = 1;
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
      } else if ((rlp_head_wall_2 == NULL) && (rlp_head_wall_3 == NULL)) {
        only_sm_2_grid_3_outside = 1;
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);
      }
    } else if ((sm_1->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, sm_1->grid->surface->parent_object, sm_1) &&
               (!(sm_2->properties->flags & CAN_REGION_BORDER) ||
                !are_restricted_regions_for_species_on_object(
                     world, sm_2->grid->surface->parent_object, sm_2)) &&
               (!(grid_3->properties->flags & CAN_REGION_BORDER) ||
                !are_restricted_regions_for_species_on_object(
                     world, grid_3->grid->surface->parent_object, grid_3))) {
      /* only reactant "sm_1" has restrictive region border property */
      w_1 = sm_1->grid->surface;
      rlp_head_wall_1 = find_restricted_regions_by_wall(world, w_1, sm_1);
      if (rlp_head_wall_1 != NULL) {
        only_sm_1_inside = 1;
      } else {
        rlp_head_obj_1 =
            find_restricted_regions_by_object(world, w_1->parent_object, sm_1);
        only_sm_1_outside = 1;
      }
    } else if ((sm_2->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, sm_2->grid->surface->parent_object, sm_2) &&
               (!(sm_1->properties->flags & CAN_REGION_BORDER) ||
                !are_restricted_regions_for_species_on_object(
                     world, sm_1->grid->surface->parent_object, sm_1)) &&
               (!(grid_3->properties->flags & CAN_REGION_BORDER) ||
                !are_restricted_regions_for_species_on_object(
                     world, grid_3->grid->surface->parent_object, grid_3))) {
      /* only reactant "sm_2" has restrictive region border property */
      w_2 = sm_2->grid->surface;
      rlp_head_wall_2 = find_restricted_regions_by_wall(world, w_2, sm_2);
      if (rlp_head_wall_2 != NULL)
        only_sm_2_inside = 1;
      else {
        rlp_head_obj_2 =
            find_restricted_regions_by_object(world, w_2->parent_object, sm_2);
        only_sm_2_outside = 1;
      }
    } else if ((grid_3->properties->flags & CAN_REGION_BORDER) &&
               are_restricted_regions_for_species_on_object(
                   world, grid_3->grid->surface->parent_object, grid_3) &&
               (!(sm_1->properties->flags & CAN_REGION_BORDER) ||
                !are_restricted_regions_for_species_on_object(
                     world, sm_1->grid->surface->parent_object, sm_1)) &&
               (!(sm_2->properties->flags & CAN_REGION_BORDER) ||
                !are_restricted_regions_for_species_on_object(
                     world, sm_2->grid->surface->parent_object, sm_2))) {
      /* only reactant "grid_3" has restrictive region border property */
      w_3 = grid_3->grid->surface;
      rlp_head_wall_3 = find_restricted_regions_by_wall(world, w_3, grid_3);
      if (rlp_head_wall_3 != NULL)
        only_grid_3_inside = 1;
      else {
        rlp_head_obj_3 = find_restricted_regions_by_object(
            world, w_3->parent_object, grid_3);
        only_grid_3_outside = 1;
      }
    }
  }

  /* reacA is the molecule which initiated the reaction. */
  struct abstract_molecule *const initiator = reacA;
  short const initiatorOrient = orientA;

  /* make sure that reacA corresponds to rx->players[0],
     reacB - to rx->players[1], and reacC - to rx->players[2]
     in trimolecular reaction */
  if (reacA->properties == rx->players[0]) {
    if (reacB->properties == rx->players[2] &&
        reacB->properties != rx->players[1]) {
      /* switch B and C */
      tmp_mol = reacB;
      reacB = reacC;
      reacC = tmp_mol;

      tmp_orient = orientB;
      orientB = orientC;
      orientC = tmp_orient;
    }
  } else if (reacA->properties == rx->players[1]) {

    if (reacB->properties == rx->players[0] &&
        reacB->properties != rx->players[1]) {
      /* switch reacA and reacB */
      tmp_mol = reacB;
      reacB = reacA;
      reacA = tmp_mol;

      tmp_orient = orientB;
      orientB = orientA;
      orientA = tmp_orient;
    } else if (reacC->properties == rx->players[0]) {
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
  } else if (reacA->properties == rx->players[2]) {
    if (reacB->properties == rx->players[0]) {
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

    } else if ((reacC->properties == rx->players[0]) &&
               (reacC->properties != rx->players[2])) {
      /* switch reacA and reacC */
      tmp_mol = reacA;
      reacA = reacC;
      reacC = tmp_mol;

      tmp_orient = orientA;
      orientA = orientC;
      orientC = tmp_orient;
    }
  }

  add_reactants_to_product_list(rx, reacA, reacB, reacC, product, product_type);

  /* Determine whether any of the reactants can be replaced by a product. */
  if (product_type[0] == PLAYER_SURF_MOL) {
    num_surface_reactants++;
    if (rx_players[0] == NULL)
      replace_p1 = 1;
    else
      num_surface_reactants_to_stay++;
  }
  if (product_type[1] == PLAYER_SURF_MOL) {
    num_surface_reactants++;
    if (rx_players[1] == NULL)
      replace_p2 = 1;
    else
      num_surface_reactants_to_stay++;
  }
  if (product_type[2] == PLAYER_SURF_MOL) {
    num_surface_reactants++;
    if (rx_players[2] == NULL)
      replace_p3 = 1;
    else
      num_surface_reactants_to_stay++;
  }

  if (replace_p1 && (!replace_p2) && (!replace_p3)) {
    only_one_to_replace = 1;
  } else if ((!replace_p1) && replace_p2 && (!replace_p3)) {
    only_one_to_replace = 1;
  } else if ((!replace_p1) && (!replace_p2) && replace_p3) {
    only_one_to_replace = 1;
  }
  if (replace_p1 && (replace_p2) && (!replace_p3)) {
    two_to_replace = 1;
  } else if (replace_p1 && (!replace_p2) && replace_p3) {
    two_to_replace = 1;
  } else if ((!replace_p1) && replace_p2 && replace_p3) {
    two_to_replace = 1;
  }

  /* find out number of surface products */
  for (int n_product = rx->n_reactants; n_product < n_players; ++n_product) {
    if (rx_players[n_product] == NULL)
      continue;
    if (rx_players[n_product]->flags & ON_GRID) {
      num_surface_products++;
      if (!distinguishable(rx_players[n_product]->D, 0, EPS_C))
        num_surface_static_products++;
    }
  }

  if (num_surface_reactants >= 2)
    find_neighbor_tiles_flag = 1;
  if ((num_surface_reactants == 1) && (num_surface_products > 1))
    find_neighbor_tiles_flag = 1;

  /* Determine the point of reaction on the surface. */
  if (is_orientable) {
    if (sm_reactant)
      rxn_uv_pos = sm_reactant->s_pos;
    else {
      xyz2uv(hitpt, w, &rxn_uv_pos);
    }

    if ((w == NULL) && (sm_reactant != NULL))
      w = sm_reactant->grid->surface;

    assert(w != NULL);
    if (w->grid == NULL) {
      /* reacA must be a volume molecule, or this wall would have a grid
       * already. */
      assert(!IS_SURF_MOL(reacA));

      if (create_grid(world, w, ((struct volume_molecule *)reacA)->subvol))
        mcell_allocfailed("Failed to create a grid for a wall.");
    }

    if (find_neighbor_tiles_flag) {
      /* create list of neighbor tiles around rxn_uv_pos */
      rxn_uv_idx = uv2grid(&rxn_uv_pos, w->grid);

      if (sm_reactant != NULL) {
        find_neighbor_tiles(world, sm_reactant, sm_reactant->grid,
                            sm_reactant->grid_index, 1, 0, &tile_nbr_head,
                            &list_length);
      } else {
        find_neighbor_tiles(world, sm_reactant, w->grid, rxn_uv_idx, 1, 0,
                            &tile_nbr_head, &list_length);
      }

      /* Create list of vacant tiles */
      for (tile_nbr = tile_nbr_head; tile_nbr != NULL;
           tile_nbr = tile_nbr->next) {
        struct surface_molecule_list *sm_list = tile_nbr->grid->sm_list[tile_nbr->idx];
        if (sm_list == NULL || sm_list->sm == NULL) {
          num_vacant_tiles++;
          push_tile_neighbor_to_list(&tile_vacant_nbr_head, tile_nbr->grid,
                                     tile_nbr->idx);
        }
      }
    }
  }

  /* find out number of static surface reactants */
  if ((sm_1 != NULL) && !distinguishable(sm_1->properties->D, 0, EPS_C))
    num_surface_static_reactants++;
  if ((sm_2 != NULL) && !distinguishable(sm_2->properties->D, 0, EPS_C))
    num_surface_static_reactants++;
  if ((grid_3 != NULL) && !distinguishable(grid_3->properties->D, 0, EPS_C))
    num_surface_static_reactants++;

  /* If the reaction involves a surface, make sure there is room for each
   * product. */
  if (is_orientable) {
    /* Can this reaction happen at all? */
    if (replace_p1 && replace_p2 && replace_p3) {
      if (num_surface_products > num_vacant_tiles + 3) {
        if (tile_nbr_head != NULL)
          delete_tile_neighbor_list(tile_nbr_head);
        if (tile_vacant_nbr_head != NULL)
          delete_tile_neighbor_list(tile_vacant_nbr_head);
        return RX_BLOCKED;
      }
    } else if (two_to_replace) {
      if (num_surface_products > num_vacant_tiles + 2) {
        if (tile_nbr_head != NULL)
          delete_tile_neighbor_list(tile_nbr_head);
        if (tile_vacant_nbr_head != NULL)
          delete_tile_neighbor_list(tile_vacant_nbr_head);
        return RX_BLOCKED;
      }
    } else if (only_one_to_replace) {
      if (num_surface_products > num_vacant_tiles + 1) {
        if (tile_nbr_head != NULL)
          delete_tile_neighbor_list(tile_nbr_head);
        if (tile_vacant_nbr_head != NULL)
          delete_tile_neighbor_list(tile_vacant_nbr_head);
        return RX_BLOCKED;
      }
    } else {
      /* none of the reactants is replaced */
      if (num_surface_products > num_vacant_tiles) {
        if (tile_nbr_head != NULL)
          delete_tile_neighbor_list(tile_nbr_head);
        if (tile_vacant_nbr_head != NULL)
          delete_tile_neighbor_list(tile_vacant_nbr_head);
        return RX_BLOCKED;
      }
    }

    /* set the orientations of the products. */
    for (int n_product = 0; n_product < n_players; ++n_product) {
      /* Skip NULL reactants. */
      if (rx_players[n_product] == NULL)
        continue;

      int this_geometry = rx->geometries[i0 + n_product];

      /* Geometry of 0 means "random orientation" */
      if (this_geometry == 0)
        product_orient[n_product] = (rng_uint(world->rng) & 1) ? 1 : -1;
      else {
        /* Geometry < 0 means inverted orientation */
        if (this_geometry < 0) {
          this_geometry = -this_geometry;
          if (this_geometry > (int)rx->n_reactants)
            product_orient[n_product] =
                -product_orient[this_geometry - rx->n_reactants - 1];
          else if (this_geometry == 1)
            product_orient[n_product] = -orientA;
          else if ((this_geometry == 2) && (reacB != NULL))
            product_orient[n_product] = -orientB;
          else if ((this_geometry == 3) && (reacC != NULL))
            product_orient[n_product] = -orientC;
          else
            product_orient[n_product] = -1;
        }

        /* Geometry > 0 means "positive" orientation. */
        else {
          if (this_geometry > (int)rx->n_reactants)
            product_orient[n_product] =
                product_orient[this_geometry - rx->n_reactants - 1];
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
             * orientation-sensitive counts.  Then, we will update the
             * orientation.  Finally, we will add the molecule back into the
             * counts in its new orientation.
             */

            /* Remove molecule from counts in old orientation, if mol is
             * counted. */
            if (product[n_product]->properties->flags &
                (COUNT_CONTENTS | COUNT_ENCLOSED))
              count_region_from_scratch(world,
                                        product[n_product], /* molecule */
                                        NULL,               /* rxn pathway */
                                        -1,                 /* remove count */
                                        NULL, /* Location at which to count */
                                        w,    /* Wall on which this happened */
                                        t,    /* Time of occurrence */
                                        NULL);

            /* Force check for the unimolecular reactions
               after changing orientation.
               There are two possible cases to be covered here:
               1) when (sm->t2) was previously set to FOREVER
               2) there may be two or more unimolecular
                  reactions involving surface class that have
                  different kinetics.
             */
            if (((sm->flags & ACT_REACT) != 0) &&
                ((sm->properties->flags & CAN_SURFWALL) != 0))
              sm->t2 = 0;

            /* Set the molecule's orientation. */
            sm->orient = product_orient[n_product];

            /* Add molecule back to counts in new orientation, if mol is
             * counted. */
            if (product[n_product]->properties->flags &
                (COUNT_CONTENTS | COUNT_ENCLOSED))
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

        /* Otherwise, check if we've crossed the plane. */
        else {
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
      if ((num_surface_static_reactants == 1) &&
          (num_surface_static_products == 1) &&
          (replace_p1 || replace_p2 || replace_p3)) {
        /* the lonely static product always replaces the lonely static reactant
         */
        for (int n_product = rx->n_reactants; n_product < n_players;
             n_product++) {
          if (rx_players[n_product] == NULL)
            continue;
          if ((rx_players[n_product]->flags & NOT_FREE) == 0)
            continue;
          if (distinguishable(rx_players[n_product]->D, 0, EPS_C))
            continue;

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            if (replace_p1 && !distinguishable(reacA->properties->D, 0, EPS_C)) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacA)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacA)->grid_index;
              replace_p1 = 0;
              break;
            } else if (replace_p2 && !distinguishable(reacB->properties->D, 0, EPS_C)) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacB)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacB)->grid_index;
              replace_p2 = 0;
              break;
            } else if (replace_p3 && !distinguishable(reacC->properties->D, 0, EPS_C)) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACC_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacC)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacC)->grid_index;
              replace_p3 = 0;
              break;
            }
          }
        }
      } else if (replace_p1 && replace_p2 && replace_p3) {
        /* if all reactants should be  replaced and there is only one
           surface product here we make sure that initiator molecule
           is replaced */
        for (int n_product = rx->n_reactants; n_product < n_players;
             n_product++) {
          if (rx_players[n_product] == NULL)
            continue;
          if ((rx_players[n_product]->flags & NOT_FREE) == 0)
            continue;

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            if (reacA == initiator) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacA)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacA)->grid_index;
              replace_p1 = 0;
            } else if (reacB == initiator) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacB)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacB)->grid_index;
              replace_p2 = 0;
            } else {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACC_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacC)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacC)->grid_index;
              replace_p3 = 0;
            }
            break;
          }
        }
      } else if (two_to_replace) {
        /* replace one of the two reactants randomly */
        while (true) {
          rnd_num = rng_uint(world->rng) % (rx->n_reactants);

          if ((rnd_num == 0) && replace_p1)
            break;

          if ((rnd_num == 1) && replace_p2)
            break;

          if ((rnd_num == 2) && replace_p3)
            break;
        }

        for (int n_product = rx->n_reactants; n_product < n_players;
             n_product++) {
          if (rx_players[n_product] == NULL)
            continue;
          if ((rx_players[n_product]->flags & NOT_FREE) == 0)
            continue;

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            if (rnd_num == 0) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacA)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacA)->grid_index;
              replace_p1 = 0;
            } else if (rnd_num == 1) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacB)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacB)->grid_index;
              replace_p2 = 0;
            } else if (rnd_num == 2) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACC_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacC)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacC)->grid_index;
              replace_p3 = 0;
            }
            break;
          }
        }

      } else if (only_one_to_replace) {
        /* no need for a random number here */
        for (int n_product = rx->n_reactants; n_product < n_players;
             n_product++) {
          if (rx_players[n_product] == NULL)
            continue;
          if ((rx_players[n_product]->flags & NOT_FREE) == 0)
            continue;

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET) {
            if (replace_p1) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACA_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacA)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacA)->grid_index;
              replace_p1 = 0;
              break;
            } else if (replace_p2) {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACB_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacB)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacB)->grid_index;
              replace_p2 = 0;
              break;
            } else {
              product_flag[n_product] = PRODUCT_FLAG_USE_REACC_UV;
              product_grid[n_product] =
                  ((struct surface_molecule *)reacC)->grid;
              product_grid_idx[n_product] =
                  ((struct surface_molecule *)reacC)->grid_index;
              replace_p3 = 0;
              break;
            }
          }
        }
      }

    } else if (num_surface_products > 1) {
      int count;
      if (num_surface_static_reactants > 0) {
        if (num_surface_static_products >= num_surface_static_reactants) {
          count = 0;
          while (count < num_surface_static_reactants) {
            rnd_num = rng_uint(world->rng) % n_players;
            /* pass reactants */
            if (rnd_num < 3)
              continue;

            if (rx_players[rnd_num] == NULL)
              continue;
            if ((rx_players[rnd_num]->flags & NOT_FREE) == 0)
              continue;
            if (distinguishable(rx_players[rnd_num]->D, 0, EPS_C))
              continue;

            if (product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET) {
              if ((!distinguishable(reacA->properties->D, 0, EPS_C)) && replace_p1) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                product_grid[rnd_num] =
                    ((struct surface_molecule *)reacA)->grid;
                product_grid_idx[rnd_num] =
                    ((struct surface_molecule *)reacA)->grid_index;
                replace_p1 = 0;
                count++;
                continue;
              }
              if ((!distinguishable(reacB->properties->D, 0, EPS_C)) && replace_p2) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                product_grid[rnd_num] =
                    ((struct surface_molecule *)reacB)->grid;
                product_grid_idx[rnd_num] =
                    ((struct surface_molecule *)reacB)->grid_index;
                replace_p2 = 0;
                count++;
                continue;
              }
              if ((!distinguishable(reacC->properties->D, 0, EPS_C)) && replace_p3) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACC_UV;
                product_grid[rnd_num] =
                    ((struct surface_molecule *)reacC)->grid;
                product_grid_idx[rnd_num] =
                    ((struct surface_molecule *)reacC)->grid_index;
                replace_p3 = 0;
                count++;
                continue;
              }
            }
          } /* end while */

        } else { /*(num_surface_static_products<num_surface_static_reactants)*/
          count = 0;
          while (count < num_surface_static_products) {
            rnd_num = rng_uint(world->rng) % n_players;
            /* pass reactants */
            if (rnd_num < 3)
              continue;

            if (rx_players[rnd_num] == NULL)
              continue;
            if ((rx_players[rnd_num]->flags & NOT_FREE) == 0)
              continue;
            if (distinguishable(rx_players[rnd_num]->D, 0, EPS_C))
              continue;

            if (product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET) {
              if ((!distinguishable(reacA->properties->D, 0, EPS_C)) && replace_p1) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                product_grid[rnd_num] =
                    ((struct surface_molecule *)reacA)->grid;
                product_grid_idx[rnd_num] =
                    ((struct surface_molecule *)reacA)->grid_index;
                replace_p1 = 0;
                count++;
                continue;
              }
              if ((!distinguishable(reacB->properties->D, 0, EPS_C)) && replace_p2) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                product_grid[rnd_num] =
                    ((struct surface_molecule *)reacB)->grid;
                product_grid_idx[rnd_num] =
                    ((struct surface_molecule *)reacB)->grid_index;
                replace_p2 = 0;
                count++;
                continue;
              }
              if ((!distinguishable(reacC->properties->D, 0, EPS_C)) && replace_p3) {
                product_flag[rnd_num] = PRODUCT_FLAG_USE_REACC_UV;
                product_grid[rnd_num] =
                    ((struct surface_molecule *)reacC)->grid;
                product_grid_idx[rnd_num] =
                    ((struct surface_molecule *)reacC)->grid_index;
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

      if (replace_p1 || replace_p2 || replace_p3) {
        /* are there any surface products left that have not replaced yet
            reactants? */

        int surf_prod_left = 0, surf_reactant_left = 0;

        for (int n_product = rx->n_reactants; n_product < n_players;
             n_product++) {
          if (rx_players[n_product] == NULL)
            continue;
          if ((rx_players[n_product]->flags & NOT_FREE) == 0)
            continue;

          if (product_flag[n_product] == PRODUCT_FLAG_NOT_SET)
            surf_prod_left++;
        }

        if (replace_p1)
          surf_reactant_left++;
        if (replace_p2)
          surf_reactant_left++;
        if (replace_p3)
          surf_reactant_left++;

        if (surf_prod_left > 0) {

          if (surf_prod_left >= surf_reactant_left) {
            count = 0;
            while (count < surf_reactant_left) {
              rnd_num = rng_uint(world->rng) % n_players;
              /* pass reactants */
              if (rnd_num < 3)
                continue;

              if (rx_players[rnd_num] == NULL)
                continue;
              if ((rx_players[rnd_num]->flags & NOT_FREE) == 0)
                continue;

              if (product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET) {
                if (replace_p1) {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                  product_grid[rnd_num] =
                      ((struct surface_molecule *)reacA)->grid;
                  product_grid_idx[rnd_num] =
                      ((struct surface_molecule *)reacA)->grid_index;
                  replace_p1 = 0;
                  count++;
                  continue;
                }
                if (replace_p2) {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                  product_grid[rnd_num] =
                      ((struct surface_molecule *)reacB)->grid;
                  product_grid_idx[rnd_num] =
                      ((struct surface_molecule *)reacB)->grid_index;
                  replace_p2 = 0;
                  count++;
                  continue;
                }
                if (replace_p3) {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACC_UV;
                  product_grid[rnd_num] =
                      ((struct surface_molecule *)reacC)->grid;
                  product_grid_idx[rnd_num] =
                      ((struct surface_molecule *)reacC)->grid_index;
                  replace_p3 = 0;
                  count++;
                  continue;
                }
              }
            } /* end while */

          } else { /* surf_prod_left < surf_reactant_left */
            count = 0;
            while (count < surf_prod_left) {
              rnd_num = rng_uint(world->rng) % n_players;
              /* pass reactants */
              if (rnd_num < 3)
                continue;

              if (rx_players[rnd_num] == NULL)
                continue;
              if ((rx_players[rnd_num]->flags & NOT_FREE) == 0)
                continue;

              if (product_flag[rnd_num] == PRODUCT_FLAG_NOT_SET) {
                if (replace_p1) {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACA_UV;
                  product_grid[rnd_num] =
                      ((struct surface_molecule *)reacA)->grid;
                  product_grid_idx[rnd_num] =
                      ((struct surface_molecule *)reacA)->grid_index;
                  replace_p1 = 0;
                  count++;
                  continue;
                }
                if (replace_p2) {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACB_UV;
                  product_grid[rnd_num] =
                      ((struct surface_molecule *)reacB)->grid;
                  product_grid_idx[rnd_num] =
                      ((struct surface_molecule *)reacB)->grid_index;
                  replace_p2 = 0;
                  count++;
                  continue;
                }
                if (replace_p3) {
                  product_flag[rnd_num] = PRODUCT_FLAG_USE_REACC_UV;
                  product_grid[rnd_num] =
                      ((struct surface_molecule *)reacC)->grid;
                  product_grid_idx[rnd_num] =
                      ((struct surface_molecule *)reacC)->grid_index;
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
    for (int n_product = rx->n_reactants; n_product < n_players; ++n_product) {
      /* If the product is a volume product, no placement is required. */
      if (rx_players[n_product]->flags & ON_GRID) {
        if (product_flag[n_product] != PRODUCT_FLAG_NOT_SET)
          continue;

        if (num_vacant_tiles == 0) {
          if (tile_nbr_head != NULL)
            delete_tile_neighbor_list(tile_nbr_head);
          if (tile_vacant_nbr_head != NULL)
            delete_tile_neighbor_list(tile_vacant_nbr_head);
          return RX_BLOCKED;
        }

        while (true) {
          if (num_attempts > SURFACE_DIFFUSION_RETRIES) {
            if (tile_nbr_head != NULL)
              delete_tile_neighbor_list(tile_nbr_head);
            if (tile_vacant_nbr_head != NULL)
              delete_tile_neighbor_list(tile_vacant_nbr_head);
            return RX_BLOCKED;
          }

          /* randomly pick a tile from the list */
          rnd_num = rng_uint(world->rng) % num_vacant_tiles;
          tile_idx = -1;
          tile_grid = NULL;

          if (get_tile_neighbor_from_list_of_vacant_neighbors(
                  tile_vacant_nbr_head, rnd_num, &tile_grid, &tile_idx) == 0) {
            if (tile_nbr_head != NULL)
              delete_tile_neighbor_list(tile_nbr_head);
            if (tile_vacant_nbr_head != NULL)
              delete_tile_neighbor_list(tile_vacant_nbr_head);
            return RX_BLOCKED;
          }
          if (tile_idx < 0)
            continue; /* this tile was checked out before */

          assert(tile_grid != NULL);

          if (all_inside_restricted_boundary) {
            /* if this tile is not inside the restricted boundary
               - try again */
            int cond_1 = 0, cond_2 = 0, cond_3 = 0;
            cond_1 = (!wall_belongs_to_all_regions_in_region_list(
                           tile_grid->surface, rlp_head_wall_1));
            cond_2 = (!wall_belongs_to_all_regions_in_region_list(
                           tile_grid->surface, rlp_head_wall_2));
            cond_3 = (!wall_belongs_to_all_regions_in_region_list(
                           tile_grid->surface, rlp_head_wall_3));

            if (cond_1 || cond_2 || cond_3) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (all_outside_restricted_boundary) {
            int cond_1 = 0, cond_2 = 0, cond_3 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_1);
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_2);
            cond_3 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_3);

            if (cond_1 || cond_2 || cond_3) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (sm_1_inside_sm_2_inside_grid_3_outside) {
            int cond_1 = 0, cond_2 = 0, cond_3 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_1));
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_2));
            cond_3 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_3);

            if (cond_1 || cond_2 || cond_3) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (sm_1_inside_sm_2_outside_grid_3_inside) {
            int cond_1 = 0, cond_2 = 0, cond_3 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_1));
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_2);
            cond_3 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_3));

            if (cond_1 || cond_2 || cond_3) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (sm_1_inside_sm_2_outside_grid_3_outside) {
            int cond_1 = 0, cond_2 = 0, cond_3 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_1));
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_2);
            cond_3 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_3);

            if (cond_1 || cond_2 || cond_3) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (sm_1_outside_sm_2_inside_grid_3_outside) {
            int cond_1 = 0, cond_2 = 0, cond_3 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_1);
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_2));
            cond_3 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_3);

            if (cond_1 || cond_2 || cond_3) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (sm_1_outside_sm_2_inside_grid_3_inside) {
            int cond_1 = 0, cond_2 = 0, cond_3 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_1);
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_2));
            cond_3 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_3));

            if (cond_1 || cond_2 || cond_3) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (sm_1_outside_sm_2_outside_grid_3_inside) {
            int cond_1 = 0, cond_2 = 0, cond_3 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_1);
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_2);
            cond_3 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_3));

            if (cond_1 || cond_2 || cond_3) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_sm_2_inside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_1));
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_2));

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_inside_sm_2_outside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_1));
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_2);

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_outside_sm_2_inside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_1);
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_2));

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_sm_2_outside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_1);
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_2);

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_grid_3_inside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_1));
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_3));

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_inside_grid_3_outside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_1));
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_3);

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_outside_grid_3_inside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_1);
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_3));

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_grid_3_outside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_1);
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_3);

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_2_grid_3_inside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_2));
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_3));

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_2_inside_grid_3_outside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_2));
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_3);

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_2_outside_grid_3_inside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_2);
            cond_2 = !(wall_belongs_to_all_regions_in_region_list(
                          tile_grid->surface, rlp_head_wall_3));

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_2_grid_3_outside) {
            int cond_1 = 0, cond_2 = 0;
            cond_1 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_2);
            cond_2 = wall_belongs_to_any_region_in_region_list(
                tile_grid->surface, rlp_head_obj_3);

            if (cond_1 || cond_2) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_inside) {
            if (!wall_belongs_to_all_regions_in_region_list(tile_grid->surface,
                                                            rlp_head_wall_1)) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_1_outside) {
            if (wall_belongs_to_any_region_in_region_list(tile_grid->surface,
                                                          rlp_head_obj_1)) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_2_inside) {
            if (!wall_belongs_to_all_regions_in_region_list(tile_grid->surface,
                                                            rlp_head_wall_2)) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_sm_2_outside) {
            if (wall_belongs_to_any_region_in_region_list(tile_grid->surface,
                                                          rlp_head_obj_2)) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_grid_3_inside) {
            if (!wall_belongs_to_all_regions_in_region_list(tile_grid->surface,
                                                            rlp_head_wall_3)) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          } else if (only_grid_3_outside) {
            if (wall_belongs_to_any_region_in_region_list(tile_grid->surface,
                                                          rlp_head_obj_3)) {
              uncheck_vacant_tile(tile_vacant_nbr_head, rnd_num);
              num_attempts++;
              continue;
            }
          }

          product_grid[n_product] = tile_grid;
          product_grid_idx[n_product] = tile_idx;
          product_flag[n_product] = PRODUCT_FLAG_USE_RANDOM;

          break;
        }
      }
    }

  } /* end if (is_orientable) */

  /* Determine the location of the reaction for count purposes. */
  struct vector3 count_pos_xyz;
  if (hitpt != NULL)
    count_pos_xyz = *hitpt;
  else if (sm_reactant)
    uv2xyz(&sm_reactant->s_pos, sm_reactant->grid->surface, &count_pos_xyz);
  else
    count_pos_xyz = ((struct volume_molecule *)reacA)->pos;

  /* Create and place each product. */
  struct vector3 mol_pos_tmp;
  struct subvolume *product_subvol = NULL;
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
          prod_uv_pos = ((struct surface_molecule *)reacA)->s_pos;
          break;

        case PRODUCT_FLAG_USE_REACB_UV:
          prod_uv_pos = ((struct surface_molecule *)reacB)->s_pos;
          break;

        case PRODUCT_FLAG_USE_REACC_UV:
          prod_uv_pos = ((struct surface_molecule *)reacC)->s_pos;
          break;

        case PRODUCT_FLAG_USE_RANDOM:
          grid2uv_random(product_grid[n_product], product_grid_idx[n_product],
                         &prod_uv_pos, world->rng);
          break;

        default:
          UNHANDLED_CASE(product_flag[n_product]);
          /*break;*/
        }
      } else
        grid2uv(product_grid[n_product], product_grid_idx[n_product],
                &prod_uv_pos);

      this_product = (struct abstract_molecule *)place_sm_product(
          world, product_species, product_grid[n_product],
          product_grid_idx[n_product], &prod_uv_pos, product_orient[n_product],
          t, reacA->periodic_box);
    }

    /* else place the molecule in space. */
    else {
      /* Unless this is a unimolecular reaction, we will have a hitpt. */
      if (!hitpt) {
        /* If this is a unimolecular surface rxn... */
        if (reacA->properties->flags & ON_GRID) {
          uv2xyz(&((struct surface_molecule *)reacA)->s_pos,
                 ((struct surface_molecule *)reacA)->grid->surface,
                 &mol_pos_tmp);
          product_subvol = find_subvolume(world, &mol_pos_tmp, last_subvol);
        }

        /* ... else a unimolecular volume rxn. */
        else {
          mol_pos_tmp = ((struct volume_molecule *)reacA)->pos;
          product_subvol = ((struct volume_molecule *)reacA)->subvol;
        }
        hitpt = &mol_pos_tmp;
      } else
        product_subvol = find_subvolume(world, hitpt, last_subvol);

      this_product = (struct abstract_molecule *)place_volume_product(
          world, product_species, sm_reactant, w, product_subvol, hitpt,
          product_orient[n_product], t, reacA->periodic_box);

      if (((struct volume_molecule *)this_product)->index < DISSOCIATION_MAX)
        update_dissociation_index = true;
    }

    /* Update molecule counts */
    ++product_species->population;
    if (product_species->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
      count_region_from_scratch(world, this_product, NULL, 1, NULL, NULL, t, NULL);
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
                                &count_pos_xyz, w, t, NULL);

    /* Other magical stuff.  For now, can only trigger releases. */
    if (rx->info[path].pathname->magic != NULL) {
      if (reaction_wizardry(world, rx->info[path].pathname->magic, w,
                            &count_pos_xyz, t))
        mcell_allocfailed("Failed to complete reaction triggered release after "
                          "a '%s' reaction.",
                          rx->info[path].pathname->sym->name);
    }
  }

  if (tile_nbr_head != NULL)
    delete_tile_neighbor_list(tile_nbr_head);
  if (tile_vacant_nbr_head != NULL)
    delete_tile_neighbor_list(tile_vacant_nbr_head);

  return cross_wall ? RX_FLIP : RX_A_OK;
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
int outcome_trimolecular(struct volume *world, struct rxn *rx, int path,
                         struct abstract_molecule *reacA,
                         struct abstract_molecule *reacB,
                         struct abstract_molecule *reacC, short orientA,
                         short orientB, short orientC, double t,
                         struct vector3 *hitpt, struct vector3 *loc_okay) {
  struct wall *w = NULL;
  struct volume_molecule *vm = NULL;
  struct surface_molecule *sm = NULL;
  int result;
  /* flags */
  int killA = 0, killB = 0, killC = 0;
  int reacA_is_free = 0;
  int reacB_is_free = 0;
  int reacC_is_free = 0;
  int num_surface_reactants = 0;

  if ((reacA->properties->flags & NOT_FREE) == 0) {
    reacA_is_free = 1;
  } else
    num_surface_reactants++;

  if ((reacB->properties->flags & NOT_FREE) == 0)
    reacB_is_free = 1;
  else
    num_surface_reactants++;

  if ((reacC->properties->flags & NOT_FREE) == 0)
    reacC_is_free = 1;
  else
    num_surface_reactants++;

  if (!reacA_is_free) {
    sm = (struct surface_molecule *)reacA;
  } else if (!reacB_is_free) {
    sm = (struct surface_molecule *)reacB;
  } else if (!reacC_is_free) {
    sm = (struct surface_molecule *)reacC;
  }
  if (sm != NULL)
    w = sm->grid->surface;

  result = outcome_products_trimol_reaction_random(world, w, hitpt, t, rx, path,
                                                   reacA, reacB, reacC, orientA,
                                                   orientB, orientC);
  if (result == RX_BLOCKED)
    return RX_BLOCKED;

  rx->n_occurred++;
  rx->info[path].count++;

  /* Figure out if either of the reactants was destroyed */

  if (rx->players[0] == reacA->properties) {
    if (rx->players[1] == reacB->properties) {
      killC = (rx->players[rx->product_idx[path] + 2] == NULL);
      killB = (rx->players[rx->product_idx[path] + 1] == NULL);
      killA = (rx->players[rx->product_idx[path]] == NULL);
    } else {
      killB = (rx->players[rx->product_idx[path] + 2] == NULL);
      killC = (rx->players[rx->product_idx[path] + 1] == NULL);
      killA = (rx->players[rx->product_idx[path]] == NULL);
    }
  } else if (rx->players[0] == reacB->properties) {
    if (rx->players[1] == reacA->properties) {
      killC = (rx->players[rx->product_idx[path] + 2] == NULL);
      killA = (rx->players[rx->product_idx[path] + 1] == NULL);
      killB = (rx->players[rx->product_idx[path]] == NULL);
    } else {
      killA = (rx->players[rx->product_idx[path] + 2] == NULL);
      killC = (rx->players[rx->product_idx[path] + 1] == NULL);
      killB = (rx->players[rx->product_idx[path]] == NULL);
    }
  } else if (rx->players[0] == reacC->properties) {
    if (rx->players[1] == reacA->properties) {
      killB = (rx->players[rx->product_idx[path] + 2] == NULL);
      killA = (rx->players[rx->product_idx[path] + 1] == NULL);
      killC = (rx->players[rx->product_idx[path]] == NULL);
    } else {
      killA = (rx->players[rx->product_idx[path] + 2] == NULL);
      killB = (rx->players[rx->product_idx[path] + 1] == NULL);
      killC = (rx->players[rx->product_idx[path]] == NULL);
    }
  }

  if (killC) {
    vm = NULL;
    if ((reacC->properties->flags & ON_GRID) != 0) {
      sm = (struct surface_molecule *)reacC;
      remove_surfmol_from_list(&sm->grid->sm_list[sm->grid_index], sm);
      sm->grid->n_occupied--;
      if (sm->flags & IN_SURFACE)
        sm->flags -= IN_SURFACE;

      if (sm->flags & IN_SCHEDULE) {
        sm->grid->subvol->local_storage->timer->defunct_count++;
      }
    } else {
      vm = (struct volume_molecule *)reacC;
      vm->subvol->mol_count--;
      if (vm->flags & IN_SCHEDULE) {
        vm->subvol->local_storage->timer->defunct_count++;
      }
    }

    if ((reacC->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) != 0) {
      count_region_from_scratch(world, reacC, NULL, -1, NULL, NULL, t, NULL);
    }

    reacC->properties->n_deceased++;
    double t_time = convert_iterations_to_seconds(
        world->start_iterations, world->time_unit,
        world->simulation_start_seconds, t);
    reacC->properties->cum_lifetime_seconds += t_time - reacC->birthday;
    reacC->properties->population--;
    if (vm != NULL)
      collect_molecule(vm);
    else {
      reacC->properties = NULL;
      if ((reacC->flags & IN_MASK) == 0)
        mem_put(reacC->birthplace, reacC);
    }
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
    } else {
      vm = (struct volume_molecule *)reacB;
      vm->subvol->mol_count--;
      if (vm->flags & IN_SCHEDULE) {
        vm->subvol->local_storage->timer->defunct_count++;
      }
    }

    if ((reacB->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) != 0) {
      count_region_from_scratch(world, reacB, NULL, -1, NULL, NULL, t, NULL);
    }

    reacB->properties->n_deceased++;
    double t_time = convert_iterations_to_seconds(
        world->start_iterations, world->time_unit,
        world->simulation_start_seconds, t);
    reacB->properties->cum_lifetime_seconds += t_time - reacB->birthday;
    reacB->properties->population--;
    if (vm != NULL)
      collect_molecule(vm);
    else {
      reacB->properties = NULL;
      if ((reacB->flags & IN_MASK) == 0)
        mem_put(reacB->birthplace, reacB);
    }
  }

  if (killA) {
    vm = NULL;
    if ((reacA->properties->flags & ON_GRID) != 0) {
      sm = (struct surface_molecule *)reacA;
      remove_surfmol_from_list(&sm->grid->sm_list[sm->grid_index], sm);
      sm->grid->n_occupied--;
      if (sm->flags & IN_SURFACE)
        sm->flags -= IN_SURFACE;

      if (sm->flags & IN_SCHEDULE) {
        sm->grid->subvol->local_storage->timer->defunct_count++;
      }
    } else {
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
        count_region_from_scratch(world, reacA, NULL, -1, NULL, NULL, t, NULL);
      }
    } else if ((reacA->flags & COUNT_ME) && world->place_waypoints_flag) {
      /* Subtlety: we made it up to hitpt, but our position is wherever we were
       * before that! */
      if (hitpt == NULL || (reacB_is_free && reacC_is_free))
          /* Vol-vol-vol rx should be counted at hitpt */
      {
        count_region_from_scratch(world, reacA, NULL, -1, hitpt, NULL, t, NULL);
      } else /* reaction involving surface or surface_molecule but we don't want
                to
                count exactly on a wall or we might count on the wrong side */
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

        count_region_from_scratch(world, reacA, NULL, -1, &fake_hitpt, NULL, t, NULL);
      }
    }
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
