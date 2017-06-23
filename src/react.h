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

#pragma once

#include <stdbool.h>

#include "mcell_structs.h"

enum {
  PRODUCT_FLAG_NOT_SET,
  PRODUCT_FLAG_USE_UV_LOC,
  PRODUCT_FLAG_USE_REACA_UV,
  PRODUCT_FLAG_USE_REACB_UV,
  PRODUCT_FLAG_USE_REACC_UV,
  PRODUCT_FLAG_USE_RANDOM
};

enum {
  PLAYER_SURF_MOL = 'g',
  PLAYER_VOL_MOL = 'm',
  PLAYER_WALL = 'w',
  PLAYER_NONE = '\0',
};

#define IS_SURF_MOL(g) ((g) != NULL && ((g)->properties->flags & ON_GRID))

/* In react_trig.c */
struct rxn *trigger_unimolecular(struct rxn **reaction_hash, int hashsize,
                                 u_int hash, struct abstract_molecule *reac);

int trigger_surface_unimol(struct rxn **reaction_hash, int rx_hashsize,
                           struct species *all_mols,
                           struct species *all_volume_mols,
                           struct species *all_surface_mols,
                           struct abstract_molecule *reac, struct wall *w,
                           struct rxn **matching_rxns);

int trigger_bimolecular_preliminary(struct rxn **reaction_hash, int hashsize,
                                    u_int hashA, u_int hashB,
                                    struct species *reacA,
                                    struct species *reacB);

int trigger_bimolecular(struct rxn **reaction_hash, int rx_hashsize,
                        u_int hashA, u_int hashB,
                        struct abstract_molecule *reacA,
                        struct abstract_molecule *reacB, short orientA,
                        short orientB, struct rxn **matching_rxns);

int trigger_trimolecular(struct rxn **reaction_hash, int rx_hashsize,
                         u_int hashA, u_int hashB, u_int hashC,
                         struct species *reacA, struct species *reacB,
                         struct species *reacC, int orientA, int orientB,
                         int orientC, struct rxn **matching_rxns);

int trigger_intersect(struct rxn **reaction_hash, int rx_hashsize,
                      struct species *all_mols, struct species *all_volume_mols,
                      struct species *all_surface_mols, u_int hashA,
                      struct abstract_molecule *reacA, short orientA,
                      struct wall *w, struct rxn **matching_rxns,
                      int allow_rx_transp, int allow_rx_reflec,
                      int allow_rx_absorb_reg_border);

void compute_lifetime(struct volume *state,
                      struct rxn *r,
                      struct abstract_molecule *am);

int check_for_unimolecular_reaction(struct volume *state,
                                    struct abstract_molecule *am);

struct rxn *pick_unimolecular_reaction(struct volume *state,
                                       struct abstract_molecule *am);

int find_unimol_reactions_with_surf_classes(
    struct rxn **reaction_hash, int rx_hashsize,
    struct abstract_molecule *reacA, struct wall *w, u_int hashA, int orientA,
    int num_matching_rxns, int allow_rx_transp, int allow_rx_reflec,
    int allow_rx_absorb_reg_border, struct rxn **matching_rxns);

int find_surface_mol_reactions_with_surf_classes(
    struct rxn **reaction_hash, int rx_hashsize, struct species *all_mols,
    struct species *all_surface_mols, int orientA, struct species *scl,
    int num_matching_rxns, int allow_rx_transp, int allow_rx_reflec,
    int allow_rx_absorb_reg_border, struct rxn **matching_rxns);

int find_volume_mol_reactions_with_surf_classes(
    struct rxn **reaction_hash, int rx_hashsize, struct species *all_mols,
    struct species *all_volume_mols, int orientA, struct species *scl,
    int num_matching_rxns, int allow_rx_transp, int allow_rx_reflec,
    struct rxn **matching_rxns);

/* In react_cond.c */
double timeof_unimolecular(struct rxn *rx, struct abstract_molecule *a,
                           struct rng_state *rng);

int which_unimolecular(struct rxn *rx, struct abstract_molecule *a,
                       struct rng_state *rng);

int binary_search_double(double *A, double match, int max, double mult);

int test_bimolecular(struct rxn *rx, double scaling, double local_prob_factor,
                     struct abstract_molecule *a1, struct abstract_molecule *a2,
                     struct rng_state *rng);

int test_many_bimolecular(struct rxn **rx, double *scaling,
                          double local_prob_factor, int n, int *chosen_pathway,
                          struct rng_state *rng,
                          int all_neighbors_flag);

int test_many_reactions_all_neighbors(struct rxn **rx, double *scaling,
                                      double *local_prob_factor, int n,
                                      int *chosen_pathway,
                                      struct rng_state *rng);

int test_intersect(struct rxn *rx, double scaling, struct rng_state *rng);

int test_many_intersect(struct rxn **rx, double scaling, int n,
                        int *chosen_pathway, struct rng_state *rng);

struct rxn *test_many_unimol(struct rxn **rx, int n,
                             struct abstract_molecule *a,
                             struct rng_state *rng);

void update_probs(struct volume *world, struct rxn *rx, double t);

/* In react_outc.c */
int outcome_unimolecular(struct volume *world, struct rxn *rx, int path,
                         struct abstract_molecule *reac, double t);

int outcome_bimolecular(struct volume *world, struct rxn *rx, int path,
                        struct abstract_molecule *reacA,
                        struct abstract_molecule *reacB, short orientA,
                        short orientB, double t, struct vector3 *hitpt,
                        struct vector3 *loc_okay);

int outcome_trimolecular(struct volume *world, struct rxn *rx, int path,
                         struct abstract_molecule *reacA,
                         struct abstract_molecule *reacB,
                         struct abstract_molecule *reacC, short orientA,
                         short orientB, short orientC, double t,
                         struct vector3 *hitpt, struct vector3 *loc_okay);

int outcome_intersect(struct volume *world, struct rxn *rx, int path,
                      struct wall *surface, struct abstract_molecule *reac,
                      short orient, double t, struct vector3 *hitpt,
                      struct vector3 *loc_okay);

int is_compatible_surface(void *req_species, struct wall *w);

void add_reactants_to_product_list(struct rxn *rx, struct abstract_molecule *reacA,
  struct abstract_molecule *reacB, struct abstract_molecule *reacC,
  struct abstract_molecule **player, char *player_type);

struct surface_molecule *
place_sm_product(struct volume *world, struct species *product_species,
                 struct surface_grid *grid, int grid_index,
                 struct vector2 *mol_uv_pos, short orient, double t,
                 struct periodic_image *periodic_box);

int reaction_wizardry(struct volume *world, struct magic_list *incantation,
                      struct wall *surface, struct vector3 *hitpt, double t);

void tiny_diffuse_3D(
    struct volume *world,
    struct subvolume *subvol,
    struct vector3 *displacement,
    struct vector3 *pos,
    struct wall *w);

struct volume_molecule *
place_volume_product(struct volume *world, struct species *product_species,
                     struct surface_molecule *sm_reactant, struct wall *w,
                     struct subvolume *subvol, struct vector3 *hitpt,
                     short orient, double t, struct periodic_image *periodic_box);

/* ALL_INSIDE: flag that indicates that all reactants lie inside their
 *             respective restrictive regions
 * ALL_OUTSIDE: flag that indicates that all reactants lie outside
 *              their respective restrictive regions
 * SURF1_IN_SURF2_OUT: flag that indicates that  reactant "sm_1" lies
 *                     inside and reactant "sm_2" lies outside of their
 *                     respective restrictive regions
 * SURF1_OUT_SURF2_IN: flag that indicates that  reactant "sm_1" lies outside
 *                     and reactant "sm_2" lies inside of their
 *                     respective restrictive regions
 * SURF1_IN: flag that indicates that only reactant "sm_1" has
 *          restrictive regions on the object and it lies
 *          inside its restrictive region.
 * SURF1_OUT: flag that indicates that only reactant "sm_1" has
 *            restrictive regions on the object and it lies
 *            outside its restrictive region.
 * SURF2_IN: flag that indicates that only reactant "sm_2" has
 *           restrictive regions on the object and it lies
 *           inside its restrictive region.
 * SURF2_OUT: flag that indicates that only reactant "sm_2" has
 *            restrictive regions on the object and it lies
 *            outside its restrictive region.  */
#define ALL_INSIDE 0x01
#define ALL_OUTSIDE 0x02
#define SURF1_IN_SURF2_OUT 0x04
#define SURF1_OUT_SURF2_IN 0x08
#define SURF1_IN 0x10
#define SURF1_OUT 0x20
#define SURF2_IN 0x40
#define SURF2_OUT 0x80

int determine_molecule_region_topology(
    struct volume *world, struct surface_molecule *sm_1,
    struct surface_molecule *sm_2, struct region_list **rlp_wall_1_ptr,
    struct region_list **rlp_wall_2_ptr, struct region_list **rlp_obj_1_ptr,
    struct region_list **rlp_obj_2_ptr, bool is_unimol);

bool product_tile_can_be_reached(struct wall *target,
                                 struct region_list *rlp_head_wall_1,
                                 struct region_list *rlp_head_wall_2,
                                 struct region_list *rlp_head_obj_1,
                                 struct region_list *rlp_head_obj_2,
                                 int sm_bitmask, bool is_unimol);
