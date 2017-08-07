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

#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "diffuse_util.h"
#include "sym_table.h"
#include "logging.h"
#include "react_util.h"
#include "strfunc.h"
#include "react.h"
#include "mcell_reactions.h"
#include "diffuse.h"
#include "vol_util.h"
#include "mcell_structs.h"

/* static helper functions */
static char *concat_rx_name(char *name1, char *name2);

static MCELL_STATUS extract_reactants(struct pathway *path,
                                      struct mcell_species *reactants,
                                      int *num_reactants, int *num_vol_mols,
                                      int *num_surface_mols,
                                      int *all_3d,
                                      int *oriented_count);

static MCELL_STATUS extract_catalytic_arrow(struct pathway *path,
                                            struct reaction_arrow *react_arrow,
                                            int *num_reactants,
                                            int *num_vol_mols,
                                            int *num_surface_mols, int *all_3d,
                                            int *oriented_count);

static MCELL_STATUS extract_surface(struct pathway *path,
                                    struct mcell_species *surf_class,
                                    int *num_reactants,
                                    unsigned int *num_surfaces,
                                    int *oriented_count);

static MCELL_STATUS extract_products(struct notifications *notify,
                                     struct pathway *path,
                                     struct mcell_species *products,
                                     int *num_surf_products,
                                     int bidirectional,
                                     int all_3d);

static MCELL_STATUS check_surface_specs(struct notifications *notify,
                                        int num_reactants, int num_surfaces,
                                        int num_vol_mols, int all_3d,
                                        int oriented_count);

static MCELL_STATUS add_catalytic_species_to_products(struct pathway *path,
                                                      int catalytic,
                                                      int bidirectional,
                                                      int all_3d);

static MCELL_STATUS invert_current_reaction_pathway(
    struct sym_table_head *rxn_sym_table,
    double vacancy_search_dist2,
    struct pathway *pathp,
    struct reaction_rate *reverse_rate,
    const char *rate_filename);

static char *create_rx_name(struct pathway *p);

static char *create_prod_signature(struct product **product_head);

static void set_reaction_player_flags(struct rxn *rx);

static int build_reaction_hash_table(
  struct rxn ***reaction_hash, int *n_reactions,
  struct sym_table_head *rxn_sym_table, int *rx_hashsize, int num_rx);

static void check_reaction_for_duplicate_pathways(struct pathway **head);

static int load_rate_file(double time_unit, struct mem_helper *tv_rxn_mem,
                          struct rxn *rx, char *fname, int path, enum warn_level_t neg_reaction);

static void add_surface_reaction_flags(struct sym_table_head *mol_sym_table,
                                       struct species *all_mols,
                                       struct species *all_surface_mols,
                                       struct species *all_volume_mols);

static void alphabetize_pathway(struct pathway *path, struct rxn *reaction);

static struct rxn *split_reaction(struct rxn *rx);

static void check_duplicate_special_reactions(struct pathway *path);

static int set_product_geometries(struct pathway *path, struct rxn *rx,
                                  struct product *prod);

static int scale_probabilities(byte *reaction_prob_limit_flag,
                               struct notifications *notify,
                               struct pathway *path, struct rxn *rx,
                               double pb_factor);

static int sort_product_list_compare(struct product *list_item,
                                     struct product *new_item);

static struct product *sort_product_list(struct product *product_head);

/*************************************************************************
 *
 * mcell_modify_multiple_rate_constants - modifies the rate constant of multiple reactions with
 * names.
 *
 *************************************************************************/
MCELL_STATUS
mcell_modify_multiple_rate_constants(struct volume *world, char **names, double *rate_constants, int n_rxns) {

  // Store what each type of reaction is
  // diffusing
  struct rxn **reactions_ud = malloc(0 * sizeof(*reactions_ud)); // Empty
  // unimolecular non-diffusing, volume
  struct rxn **reactions_undv = malloc(0 * sizeof(*reactions_undv));
  // unimolecular non-diffusing, surface
  struct rxn **reactions_unds = malloc(0 * sizeof(*reactions_unds));
  // Keep track of the size
  int n_reactions_ud=0;
  int n_reactions_undv=0;
  int n_reactions_unds=0;

  // Go through all the reactions
  int i_rxn=0;
  while(i_rxn < n_rxns)
  {

    // Grab the reaction
    struct sym_table_head *rxpn_sym_table = world->rxpn_sym_table;
    struct sym_entry *sym = retrieve_sym(names[i_rxn], rxpn_sym_table);

    // If the reaction couldn't be found by name, return fail
    if (sym == NULL) 
    {
      return MCELL_FAIL;
    }

    // Found the reaction, now do the changing

    // What is the pathway that needs to be changed?
    struct rxn_pathname *rxpn = sym->value;  
    // The reaction that owns this pathway
    struct rxn *reaction = rxpn->rx;  
    // The index of the pathway in this reaction
    int j = rxpn->path_num;

    // Check what type of reaction this is; store this to update the scheduler later
    int can_diffuse = distinguishable(reaction->players[0]->D, 0, EPS_C);
    if (reaction->n_reactants == 1 && can_diffuse)
    {
      // unimolecular reactions w/ diffusable reactants
      reactions_ud = (struct rxn **) realloc(reactions_ud, (n_reactions_ud + 1)*sizeof(*reactions_ud));
      reactions_ud[n_reactions_ud] = reaction;
      n_reactions_ud++;
    }
    else if (((!can_diffuse) && (reaction->n_reactants == 1)) || 
    ((!can_diffuse) && (reaction->n_reactants == 2) && (reaction->players[1]->flags == IS_SURFACE))) 
    {
      // unimolecular reactions w/ non-diffusable reactants

      // Surface or volume
      if ((reaction->players[0]->flags & NOT_FREE) != 0) 
      {
        // Surface
        reactions_unds = (struct rxn **) realloc(reactions_unds, (n_reactions_unds + 1)*sizeof(*reactions_unds));
        reactions_unds[n_reactions_unds] = reaction;
        n_reactions_unds++;
      }
      else
      {
        // Volume
        reactions_undv = (struct rxn **) realloc(reactions_undv, (n_reactions_undv + 1)*sizeof(*reactions_undv));
        reactions_undv[n_reactions_undv] = reaction;
        n_reactions_undv++;
      }
    }
    
    // From the new rate constant, compute the NEW probability for this pathway
    double p = rate_constants[i_rxn] * reaction->pb_factor;

    // Find the delta_prob for this pathway
    double delta_prob;
    if (j == 0)
      delta_prob = p - reaction->cum_probs[0];
    else
      delta_prob = p - (reaction->cum_probs[j] - reaction->cum_probs[j - 1]);

    // Update the prob for this pathway, but ALSO all other pathways above it
    for (int k = j; k < reaction->n_pathways; k++) 
    {
      reaction->cum_probs[k] += delta_prob;
    }
    reaction->max_fixed_p += delta_prob;
    reaction->min_noreaction_p += delta_prob;

    // Go to the next reaction that needs to be changed
    i_rxn++;
  }

  // Now, reschedule all necessary reactions at once

  // Check: are there any reactions with diffusable reactants
  if (n_reactions_ud > 0) // There is at least one
  {
    for (struct storage_list *local = world->storage_head; local != NULL; local = local->next) 
    {
      struct abstract_element *head_molecule = local->store->timer->current;
      while (local->store->timer->current != NULL) 
      {
        struct abstract_molecule *am = (struct abstract_molecule *)schedule_peak(local->store->timer);

        // Go through all types of these reactions
        for (int i_reactions_ud=0; i_reactions_ud<n_reactions_ud; i_reactions_ud++)
        {
          // We only want to update molecules involved in this reaction.
          // Also, skip dead molecs (props=NULL). They'll be cleaned up later.
          if ((am->properties != NULL) && (am->properties->species_id == reactions_ud[i_reactions_ud]->players[0]->species_id)) 
          {
            // Setting t2=0 and ACT_CHANGE will cause the lifetime to be
            // recomputed during the next timestep
            am->t2 = 0.0;
            am->flags |= ACT_CHANGE;
          }
        }
      }
      // Reset current molecule in scheduler now that we're done "peaking"
      local->store->timer->current = head_molecule;
    }
  }

  // Check: are there any reactions with NON-diffusable reagents
  // Volume case
  if (n_reactions_undv > 0) // There is at least one
  {
    int n_subvols = world->n_subvols;
    for (int i = 0; i < n_subvols; i++) 
    {
      struct subvolume *sv = &(world->subvol[i]);

      for (struct per_species_list *psl = sv->species_head; psl != NULL; psl = psl->next) 
      {
        if (psl->properties == NULL) 
        {
          continue;
        }
        for (struct volume_molecule *vm = psl->head; vm != NULL; vm = vm->next_v) 
        {
          if ((vm->properties != NULL) && (vm->t > world->current_iterations)) 
          {
            // Go through all types of these reactions
            for (int i_reactions_undv=0; i_reactions_undv<n_reactions_undv; i_reactions_undv++)
            {
              // More efficient here would be a hash table from species_id to reaction
              if (vm->properties->species_id == reactions_undv[i_reactions_undv]->players[0]->species_id)
              { 
                for (struct storage_list *local = world->storage_head; local != NULL; local = local->next) 
                {
                  vm->flags |= ACT_CHANGE;
                  vm->t2 = 0.0;
                  schedule_reschedule(local->store->timer, vm, world->current_iterations);
                }
              }
            }
          }
        }
      }
    }
  }

  // Surface case
  if (n_reactions_unds > 0) // There is at least one
  {
    int n_subvols = world->n_subvols;
    for (int i = 0; i < n_subvols; i++) 
    {
      struct subvolume *sv = &(world->subvol[i]);

      for (struct wall_list *wl = sv->wall_head; wl != NULL; wl = wl->next) 
      {
        struct surface_grid *grid = wl->this_wall->grid;
        if (grid != NULL) 
        {
          for (u_int tile_idx = 0; tile_idx < grid->n_tiles; tile_idx++) 
          {
            if (grid->sm_list[tile_idx]) 
            {
              struct surface_molecule *sm = grid->sm_list[tile_idx]->sm;
              if ((sm->properties != NULL) && (sm->t > world->current_iterations)) 
              {
                // Go through all types of these reactions
                for (int i_reactions_unds=0; i_reactions_unds<n_reactions_unds; i_reactions_unds++)
                {
                  // More efficient here would be a hash table from species_id to reaction
                  if (sm->properties->species_id == reactions_unds[i_reactions_unds]->players[0]->species_id)
                  {
                    for (struct storage_list *local = world->storage_head; local != NULL; local = local->next) 
                    {
                      sm->flags |= ACT_CHANGE;
                      sm->t2 = 0.0;
                      schedule_reschedule(local->store->timer, sm, world->current_iterations);
                    }
                  }
                } 
              }
            }
          }
        }
      }
    }
  }

  free(reactions_ud);
  free(reactions_undv);
  free(reactions_unds);

  return MCELL_SUCCESS;
}

/*************************************************************************
 *
 * mcell_modify_rate_constant - modifies the rate constant of a reaction with a
 * specified name. For example, if you have this: 
 *
 * vm + vm -> NULL [1e7] : rxn
 *
 * then you can change the rate constant from 1e7 to 0 like this:
 *
 * mcell_modify_rate_constant(world, "rxn", 0)
 *
 * NOTE: This is inefficient and needs more extensive testing
 *
 *************************************************************************/
MCELL_STATUS
mcell_modify_rate_constant(struct volume *world, char *name, double rate_constant) {

  struct sym_table_head *rxpn_sym_table = world->rxpn_sym_table;
  struct sym_entry *sym = retrieve_sym(name, rxpn_sym_table);
  if (sym == NULL) 
  {
    return MCELL_FAIL;
  }
  else 
  {
    
    // What is the pathway that needs to be changed?
    struct rxn_pathname *rxpn = sym->value;  
    // The reaction that owns this pathway
    struct rxn *reaction = rxpn->rx;  
    // The index of the pathway in this reaction
    int j = rxpn->path_num;

    // From the new rate constant, compute the NEW probability for this pathway
    double p = rate_constant * reaction->pb_factor;

    // Find the delta_prob for this pathway
    double delta_prob;
    if (j == 0)
      delta_prob = p - reaction->cum_probs[0];
    else
      delta_prob = p - (reaction->cum_probs[j] - reaction->cum_probs[j - 1]);

    // Update the prob for this pathway, but ALSO all other pathways above it
    for (int k = j; k < reaction->n_pathways; k++) 
    {
      reaction->cum_probs[k] += delta_prob;
    }
    reaction->max_fixed_p += delta_prob;
    reaction->min_noreaction_p += delta_prob;
 
    // Print if the flags are set
    /*
    if (world->notify->time_varying_reactions == NOTIFY_FULL &&
        reaction->cum_probs[j] >= world->notify->reaction_prob_notify) {

      // Print the reaction probabilities
      double new_prob;
      if (j == 0)
      {
          new_prob = reaction->cum_probs[0];
      }
      else
      {
          new_prob = reaction->cum_probs[j] - reaction->cum_probs[j - 1];
      }
      // Print the new_prob
      if (reaction->n_reactants == 1) 
      {
        mcell_log_raw("Probability %.4e set for %s[%d] -> ", new_prob,
          reaction->players[0]->sym->name, reaction->geometries[0]);
      } 
      else if (reaction->n_reactants == 2) 
      {
        mcell_log_raw("Probability %.4e set for %s[%d] + %s[%d] -> ", new_prob,
          reaction->players[0]->sym->name, reaction->geometries[0],
          reaction->players[1]->sym->name, reaction->geometries[1]);
      } 
      else 
      {
        mcell_log_raw("Probability %.4e set for %s[%d] + %s[%d] + %s[%d] -> ",
          new_prob, reaction->players[0]->sym->name, reaction->geometries[0],
          reaction->players[1]->sym->name, reaction->geometries[1],
          reaction->players[2]->sym->name, reaction->geometries[2]);
      }
      for (unsigned int n_product = reaction->product_idx[j]; 
        n_product < reaction->product_idx[j + 1]; n_product++) 
      {
          if (reaction->players[n_product] != NULL)
          {
            mcell_log_raw("%s[%d] ", reaction->players[n_product]->sym->name,
                          reaction->geometries[n_product]);
          }
      }
      mcell_log_raw("\n");
    }
    */

    // This is the old code as of 01.27.2017. Fixed based on react_cond.c
    /*
    int num_path = reaction->n_pathways;
    double p = rate_constant * reaction->pb_factor;

    double delta_prob = 0;
    if (num_path > 1) {
      delta_prob = \
      p - (reaction->cum_probs[num_path-1] - reaction->cum_probs[num_path-2]);
    }
    else {
      delta_prob = p - (reaction->cum_probs[num_path-1]);
    }
    reaction->cum_probs[num_path-1] += delta_prob;
    reaction->max_fixed_p += delta_prob;
    reaction->min_noreaction_p += delta_prob;
    */

    // Now reschedule all the necessary reactions

    int can_diffuse = distinguishable(reaction->players[0]->D, 0, EPS_C);
    // Need to recompute lifetimes for unimolecular reactions w/ diffusable
    // reactants.
    if (reaction->n_reactants == 1 && can_diffuse) 
    {
      for (struct storage_list *local = world->storage_head; local != NULL; local = local->next) 
      {
        struct abstract_element *head_molecule = local->store->timer->current;
        while (local->store->timer->current != NULL) 
        {
          struct abstract_molecule *am = (struct abstract_molecule *)schedule_peak(local->store->timer);
          // We only want to update molecules involved in this reaction.
          // Also, skip dead molecs (props=NULL). They'll be cleaned up later.
          if ((am->properties != NULL) && (am->properties->species_id == reaction->players[0]->species_id)) 
          {
            // Setting t2=0 and ACT_CHANGE will cause the lifetime to be
            // recomputed during the next timestep
            am->t2 = 0.0;
          am->flags |= ACT_CHANGE;
          }
        }
        // Reset current molecule in scheduler now that we're done "peaking"
        local->store->timer->current = head_molecule;
      }
    }

    // Need to recompute lifetimes for non-diffusing molecules that are
    // unimolecular or where you have a surface molecule at a surface class
    // (e.g. sm@sc->whatever). These molecules won't come up next in the
    // scheduler, so we have to hunt them all down... :(
    if (((!can_diffuse) && (reaction->n_reactants == 1)) || 
      ((!can_diffuse) && (reaction->n_reactants == 2) && (reaction->players[1]->flags == IS_SURFACE))) 
    {
      for (struct storage_list *local = world->storage_head; local != NULL; local = local->next) 
      {
        int n_subvols = world->n_subvols;
        for (int i = 0; i < n_subvols; i++) 
        {
          struct subvolume *sv = &(world->subvol[i]);

          // Reschedule the surface molecules involved in the reaction
          if ((reaction->players[0]->flags & NOT_FREE) != 0) 
          {
            for (struct wall_list *wl = sv->wall_head; wl != NULL; wl = wl->next) 
            {
              struct surface_grid *grid = wl->this_wall->grid;
              if (grid != NULL) 
              {
                for (u_int tile_idx = 0; tile_idx < grid->n_tiles; tile_idx++) 
                {
                  if (grid->sm_list[tile_idx]) 
                  {
                    struct surface_molecule *sm = grid->sm_list[tile_idx]->sm;
                    if ((sm->properties != NULL) && 
                      (sm->properties->species_id == reaction->players[0]->species_id) &&
                      (sm->t > world->current_iterations)) 
                    {
                      sm->flags |= ACT_CHANGE;
                      sm->t2 = 0.0;
                      schedule_reschedule(local->store->timer, sm, world->current_iterations);
                    }
                  }
                } 
              }
            }   
          }
          // Reschedule the volume molecules involved in the reaction
          else 
          {
            for (struct per_species_list *psl = sv->species_head; psl != NULL; psl = psl->next) 
            {
              if (psl->properties == NULL) 
              {
                continue;
              }
              for (struct volume_molecule *vm = psl->head; vm != NULL; vm = vm->next_v) 
              {
                if ((vm->properties != NULL) && 
                  (vm->properties->species_id == reaction->players[0]->species_id)  &&
                  (vm->t > world->current_iterations)) 
                { 
                  vm->flags |= ACT_CHANGE;
                  vm->t2 = 0.0;
                  schedule_reschedule(local->store->timer, vm, world->current_iterations);
                }
              }
            }
          }
        }
      }
    }
  }

  return MCELL_SUCCESS;
}

MCELL_STATUS
mcell_add_reaction_simplified(
    struct volume *state, 
    struct mcell_species *reactants,
    struct reaction_arrow *arrow,
    struct mcell_species *surfs,
    struct mcell_species *products,
    struct reaction_rates *rates,
    struct sym_entry *pathname) {

  mcell_add_reaction(state->notify, &state->r_step_release,
                     state->rxn_sym_table, state->radial_subdivisions,
                     state->vacancy_search_dist2, reactants, arrow, surfs,
                     products, pathname, rates, NULL, NULL);

  return MCELL_SUCCESS;
}

/*************************************************************************
 *
 * mcell_add_reaction add a single reaction described by reaction_def to
 * the simulations.
 *
 *************************************************************************/
MCELL_STATUS
mcell_add_reaction(struct notifications *notify,
                   double **r_step_release,
                   struct sym_table_head *rxn_sym_table,
                   u_int radial_subdivisions,
                   double vacancy_search_dist2,
                   struct mcell_species *reactants,
                   struct reaction_arrow *react_arrow,
                   struct mcell_species *surf_class,
                   struct mcell_species *products, struct sym_entry *pathname,
                   struct reaction_rates *rates,
                   const char *forward_rate_filename,
                   const char *backward_rate_filename) {
  char *rx_name;
  struct sym_entry *symp;
  int bidirectional = 0;
  int num_surf_products = 0;
  struct rxn *rxnp;

  /* Create pathway */
  struct pathway *pathp = (struct pathway *)CHECKED_MALLOC_STRUCT(
      struct pathway, "reaction pathway");
  if (pathp == NULL) {
    return MCELL_FAIL;
  }
  memset(pathp, 0, sizeof(struct pathway));

  /* Scan reactants, copying into the new pathway */
  int num_vol_mols = 0;
  int num_surface_mols = 0;
  int all_3d = 1;
  int reactant_idx = 0;
  int oriented_count = 0;
  if (extract_reactants(pathp, reactants, &reactant_idx, &num_vol_mols,
                        &num_surface_mols, &all_3d,
                        &oriented_count) == MCELL_FAIL) {
    free(pathp);
    return MCELL_FAIL;
  }

  /* Grab info from the arrow */
  if (react_arrow->flags & ARROW_BIDIRECTIONAL) {
    bidirectional = 1;
  }

  int catalytic = -1;
  if (react_arrow->flags & ARROW_CATALYTIC) {
    if (extract_catalytic_arrow(pathp, react_arrow, &reactant_idx,
                                &num_vol_mols, &num_surface_mols, &all_3d,
                                &oriented_count) == MCELL_FAIL) {
      free(pathp);
      return MCELL_FAIL;
    }
    catalytic = reactant_idx - 1;
  }

  /* If a surface was specified, include it */
  int surface = -1;
  unsigned int num_surfaces = 0;
  if (surf_class->mol_type != NULL) {
    if (extract_surface(pathp, surf_class, &reactant_idx, &num_surfaces,
                        &oriented_count) == MCELL_FAIL) {
      free(pathp);
      return MCELL_FAIL;
    }
    surface = reactant_idx - 1;
    all_3d = 0;
  }

  /* Create a reaction name for the pathway we're creating */
  rx_name = create_rx_name(pathp);
  if (rx_name == NULL) {
    free(pathp);
    mcell_error("Out of memory while creating reaction.");
    return MCELL_FAIL;
  }

  /* If this reaction doesn't exist, create it */
  if ((symp = retrieve_sym(rx_name, rxn_sym_table)) != NULL) {
    /* do nothing */
  } else if ((symp = store_sym(rx_name, RX, rxn_sym_table, NULL)) ==
             NULL) {
    free(pathp);
    free(rx_name);
    mcell_error("Out of memory while creating reaction.");
    /*return MCELL_FAIL;*/
  }
  free(rx_name);

  rxnp = (struct rxn *)symp->value;
  rxnp->n_reactants = reactant_idx;
  ++rxnp->n_pathways;

  /* Check for invalid reaction specifications */
  if (check_surface_specs(notify, rxnp->n_reactants, num_surfaces,
                          num_vol_mols, all_3d, oriented_count) == MCELL_FAIL) {
    free(pathp);
    return MCELL_FAIL;
  }

  /* Add catalytic reagents to the product list.
   *    - For unidirectional catalytic reactions - copy catalyst to products
   *      only if catalyst is not a surface_clas.
   *    - For bidirectional catalytic reactions always copy catalyst to
   *      products and take care that surface_class will not appear in the
   *      products later after inverting the reaction
   */
  if (catalytic >= 0) {
    if (add_catalytic_species_to_products(pathp, catalytic, bidirectional,
                                          all_3d) == MCELL_FAIL) {
      free(pathp);
      return MCELL_FAIL;
    }
  }

  /* Add in all products */
  if (extract_products(notify, pathp, products, &num_surf_products,
                       bidirectional, all_3d) == MCELL_FAIL) {
    free(pathp);
    return MCELL_FAIL;
  }
  // mem_put_list(parse_state->mol_data_list_mem, products);

  /* Attach reaction pathway name, if we have one */
  if (pathname != NULL) {
    struct rxn_pathname *rxpnp = (struct rxn_pathname *)pathname->value;
    rxpnp->rx = rxnp;
    pathp->pathname = rxpnp;
  }

  if (pathp->product_head != NULL) {
    pathp->prod_signature = create_prod_signature(&pathp->product_head);
    if (pathp->prod_signature == NULL) {
      free(pathp);
      mcell_error(
          "Error creating 'prod_signature' field for the reaction pathway.");
      return MCELL_FAIL;
    }
  } else
    pathp->prod_signature = NULL;

  /* Copy in forward rate */
  switch (rates->forward_rate.rate_type) {
  case RATE_UNSET:
    mcell_error_raw("File %s, Line %d: Internal error: Rate is not set",
                    __FILE__, __LINE__);
    free(pathp);
    return MCELL_FAIL;

  case RATE_CONSTANT:
    pathp->km = rates->forward_rate.v.rate_constant;
    pathp->km_filename = NULL;
    break;

  case RATE_FILE:
    pathp->km = 0.0;
    pathp->km_filename = (char *)forward_rate_filename;
    free(rates->forward_rate.v.rate_file);
    rates->forward_rate.v.rate_file = NULL;
    break;

  default:
    UNHANDLED_CASE(rates->forward_rate.rate_type);
  }

  /* Add the pathway to the list for this reaction */
  if (rates->forward_rate.rate_type == RATE_FILE) {
    struct pathway *tpp;
    if (rxnp->pathway_head == NULL) {
      rxnp->pathway_head = pathp;
      pathp->next = NULL;
    } else /* Move varying reactions to the end of the list */
    {
      for (tpp = rxnp->pathway_head;
           tpp->next != NULL && tpp->next->km_filename == NULL;
           tpp = tpp->next) {
      }
      pathp->next = tpp->next;
      tpp->next = pathp;
    }
  } else {
    pathp->next = rxnp->pathway_head;
    rxnp->pathway_head = pathp;
  }

  /* If we're doing 3D releases, set up array so we can release reversibly */
  if (*r_step_release == NULL && all_3d && pathp->product_head != NULL) {
    *r_step_release = init_r_step_3d_release(radial_subdivisions);
    if (*r_step_release == NULL) {
      free(pathp);
      mcell_error("Out of memory building r_step array.");
      return MCELL_FAIL;
    }
  }

  /* If the vacancy search distance is zero and this reaction produces more
   * surface molecules than it comsumes, it can never succeed, except if it is
   * a volume molecule hitting the surface and producing a single surface
   * molecule.  Fail with an error message.
   */
  if ((!distinguishable(vacancy_search_dist2, 0, EPS_C)) &&
      (num_surf_products > num_surface_mols)) {
    /* The case with one volume molecule reacting with the surface and
     * producing one surface molecule is okay.
     */
    if (num_surface_mols == 0 && num_vol_mols == 1 && num_surf_products == 1) {
      /* do nothing */
    } else {
      free(pathp);
      mcell_error("number of surface products exceeds number of surface "
                  "reactants, but VACANCY_SEARCH_DISTANCE is not specified or "
                  "set to zero.");
      return MCELL_FAIL;
    }
  }

  /* A non-reversible reaction may not specify a reverse reaction rate */
  if (rates->backward_rate.rate_type != RATE_UNSET && !bidirectional) {
    free(pathp);
    mcell_error("reverse rate specified but the reaction isn't reversible.");
    return MCELL_FAIL;
  }

  /* Create reverse reaction if we need to */
  if (bidirectional) {
    /* A bidirectional reaction must specify a reverse rate */
    if (rates->backward_rate.rate_type == RATE_UNSET) {
      free(pathp);
      mcell_error("reversible reaction indicated but no reverse rate "
                  "supplied.");
      return MCELL_FAIL;
    }

    /* if "surface_class" is present on the reactant side of the reaction copy
     * it to the product side of the reaction.
     *
     * Reversible reaction of the type:
     *    A' @ surf' <---> C''[>r1,<r2]
     *
     * is equivalent now to the two reactions:
     *    A' @ surf' ---> C'' [r1]
     *    C'' @ surf' ----> A' [r2]
     *
     * Reversible reaction of the type:
     *    A' + B' @ surf' <---> C'' + D'' [>r1,<r2]
     *
     * is equivalent now to the two reactions:
     *    A' + B @ surf' ---> C'' + D'' [r1]
     *    C'' + D'' @ surf' ----> A' + B' [r2]
     */
    if (surface != -1 && surface != catalytic) {
      struct product *prodp;
      prodp = (struct product *)CHECKED_MALLOC_STRUCT(struct product,
                                                      "reaction product");
      if (prodp == NULL) {
        // mem_put(parse_state->prod_mem, prodp);
        free(pathp);
        return MCELL_FAIL;
      }

      switch (surface) {
      case 1:
        prodp->prod = pathp->reactant2;
        prodp->orientation = pathp->orientation2;
        break;

      case 2:
        prodp->prod = pathp->reactant3;
        prodp->orientation = pathp->orientation3;
        break;

      case 0:
      default:
        mcell_internal_error(
            "Surface appears in invalid reactant slot in reaction (%d).",
            surface);
        /*break;*/
      }
      prodp->next = pathp->product_head;
      pathp->product_head = prodp;
    }

    /* Invert the current reaction pathway */
    if (invert_current_reaction_pathway(
        rxn_sym_table, vacancy_search_dist2, pathp,
        &rates->backward_rate, backward_rate_filename)) {
      free(pathp);
      return MCELL_FAIL;
    }
  }

  return MCELL_SUCCESS;
}

/*************************************************************************
 *
 * mcell_add_surface_reaction adds a single surface reaction described
 * by reaction_def to the simulations.
 *
 *************************************************************************/
MCELL_STATUS
mcell_add_surface_reaction(struct sym_table_head *rxn_sym_table,
                           int reaction_type, struct species *surface_class,
                           struct sym_entry *reactant_sym, short orient) {
  struct species *reactant = (struct species *)reactant_sym->value;

  /* Make sure the other reactant isn't a surface */
  if (reactant->flags == IS_SURFACE) {
    // mdlerror_fmt(parse_state,
    //             "Illegal reaction between two surfaces in surface reaction:
    // %s -%s-> ...",
    //             reactant_sym->name,
    //             surface_class->sym->name);
    return MCELL_FAIL;
  }

  /* Build reaction name */
  char *rx_name =
      concat_rx_name(surface_class->sym->name, reactant_sym->name);
  if (rx_name == NULL) {
    // mdlerror_fmt(parse_state,
    //             "Out of memory while parsing surface reaction: %s -%s-> ...",
    //             surface_class->sym->name,
    //             reactant_sym->name);
    return MCELL_FAIL;
  }

  /* Find or create reaction */
  struct sym_entry *reaction_sym;
  if ((reaction_sym = retrieve_sym(rx_name, rxn_sym_table)) != NULL) {
    /* do nothing */
  } else if ((reaction_sym =
                  store_sym(rx_name, RX, rxn_sym_table, NULL)) == NULL) {
    free(rx_name);
    // mdlerror_fmt(parse_state,
    //             "Out of memory while creating surface reaction: %s -%s->
    // ...",
    //             reactant_sym->name,
    //             surface_class->sym->name);
    return MCELL_FAIL;
  }
  free(rx_name);

  /* Create pathway */
  struct pathway *pathp = (struct pathway *)CHECKED_MALLOC_STRUCT(
      struct pathway, "reaction pathway");

  if (pathp == NULL)
    return MCELL_FAIL;
  memset(pathp, 0, sizeof(struct pathway));

  struct rxn *rxnp = (struct rxn *)reaction_sym->value;
  rxnp->n_reactants = 2;
  ++rxnp->n_pathways;
  pathp->pathname = NULL;
  pathp->reactant1 = surface_class;
  pathp->reactant2 = (struct species *)reactant_sym->value;
  pathp->reactant3 = NULL;
  pathp->km = GIGANTIC;
  pathp->km_filename = NULL;
  pathp->prod_signature = NULL;
  pathp->flags = 0;

  pathp->orientation1 = 1;
  pathp->orientation3 = 0;
  if (orient == 0) {
    pathp->orientation2 = 0;
  } else {
    pathp->orientation2 = (orient < 0) ? -1 : 1;
  }

  struct name_orient *no;
  no = CHECKED_MALLOC_STRUCT(struct name_orient, "struct name_orient");
  no->name = CHECKED_STRDUP(reactant->sym->name, "reactant name");
  if (orient == 0) {
    no->orient = 0;
  } else {
    no->orient = (orient < 0) ? -1 : 1;
  }

  struct product *prodp;
  switch (reaction_type) {
  case RFLCT:
    prodp = (struct product *)CHECKED_MALLOC_STRUCT(struct product,
                                                    "reaction product");
    if (prodp == NULL) {
      free(no);
      free(pathp);
      return MCELL_FAIL;
    }

    pathp->flags |= PATHW_REFLEC;
    prodp->prod = pathp->reactant2;
    prodp->orientation = 1;
    prodp->next = NULL;
    pathp->product_head = prodp;
    if (pathp->product_head != NULL) {
      pathp->prod_signature = create_prod_signature(&pathp->product_head);
      if (pathp->prod_signature == NULL) {
        // mdlerror(parse_state, "Error creating 'prod_signature' field for the
        // reaction pathway.");
        free(no);
        free(pathp);
        return MCELL_FAIL;
      }
    }
    if (surface_class->refl_mols == NULL) {
      no->next = NULL;
      surface_class->refl_mols = no;
    } else {
      no->next = surface_class->refl_mols;
      surface_class->refl_mols = no;
    }

    break;
  case TRANSP:
    prodp = (struct product *)CHECKED_MALLOC_STRUCT(struct product,
                                                    "reaction product");
    if (prodp == NULL) {
      free(no);
      free(pathp);
      return MCELL_FAIL;
    }

    pathp->flags |= PATHW_TRANSP;
    prodp->prod = pathp->reactant2;
    prodp->orientation = -1;
    prodp->next = NULL;
    pathp->product_head = prodp;
    if (pathp->product_head != NULL) {
      pathp->prod_signature = create_prod_signature(&pathp->product_head);
      if (pathp->prod_signature == NULL) {
        free(no);
        free(pathp);
        return MCELL_FAIL;
      }
    }
    if (surface_class->transp_mols == NULL) {
      no->next = NULL;
      surface_class->transp_mols = no;
    } else {
      no->next = surface_class->transp_mols;
      surface_class->transp_mols = no;
    }
    break;
  case SINK:
    pathp->flags |= PATHW_ABSORP;
    pathp->product_head = NULL;
    if (surface_class->absorb_mols == NULL) {
      no->next = NULL;
      surface_class->absorb_mols = no;
    } else {
      no->next = surface_class->absorb_mols;
      surface_class->absorb_mols = no;
    }
    break;
  default:
    // mdlerror(parse_state, "Unknown special surface type.");
    free(no);
    free(pathp);
    return MCELL_FAIL;
  }

  pathp->next = rxnp->pathway_head;
  rxnp->pathway_head = pathp;

  return MCELL_SUCCESS;
}

/*************************************************************************
 *
 * mcell_add_concentration_clamp adds a surface clamp to the simulation
 *
 *************************************************************************/
MCELL_STATUS
mcell_add_concentration_clamp(struct sym_table_head *rxn_sym_table,
                              struct species *surface_class,
                              struct sym_entry *mol_sym, short orient,
                              double conc) {
  struct rxn *rxnp;
  struct pathway *pathp;
  struct sym_entry *stp3;
  struct species *specp = (struct species *)mol_sym->value;
  struct name_orient *no;

  if (specp->flags == IS_SURFACE) {
    //    mdlerror_fmt(parse_state,
    //                "Illegal reaction between two surfaces in surface
    // reaction: %s -%s-> ...",
    //               mol_sym->name, surface_class->sym->name);
    return MCELL_FAIL;
  }
  if (specp->flags & ON_GRID) {
    // mdlerror(parse_state, "Concentration clamp does not work on surface
    // molecules.");
    return MCELL_FAIL;
  }
  if (specp->flags & NOT_FREE || specp->D <= 0.0) {
    //    mdlerror(parse_state, "Concentration clamp must be applied to molecule
    // diffusing in 3D");
    return MCELL_FAIL;
  }
  if (conc < 0) {
    // mdlerror(parse_state, "Concentration can only be clamped to positive
    // values.");
    return MCELL_FAIL;
  }

  char *rx_name = concat_rx_name(surface_class->sym->name, mol_sym->name);
  if (rx_name == NULL) {
    //    mdlerror_fmt(parse_state,
    //                 "Memory allocation error: %s -%s-> ...",
    //                 surface_class->sym->name, mol_sym->name);
    return MCELL_FAIL;
  }
  if ((stp3 = retrieve_sym(rx_name, rxn_sym_table)) != NULL) {
    /* do nothing */
  } else if ((stp3 = store_sym(rx_name, RX, rxn_sym_table, NULL)) ==
             NULL) {
    free(rx_name);
    //    mdlerror_fmt(parse_state,
    //                 "Cannot store surface reaction: %s -%s-> ...",
    //                 mol_sym->name, surface_class->sym->name);
    return MCELL_FAIL;
  }
  free(rx_name);

  pathp = (struct pathway *)CHECKED_MALLOC_STRUCT(struct pathway,
                                                  "reaction pathway");
  if (pathp == NULL)
    return MCELL_FAIL;
  memset(pathp, 0, sizeof(struct pathway));

  rxnp = (struct rxn *)stp3->value;
  rxnp->n_reactants = 2;
  ++rxnp->n_pathways;
  pathp->pathname = NULL;
  pathp->reactant1 = surface_class;
  pathp->reactant2 = (struct species *)mol_sym->value;
  pathp->reactant3 = NULL;
  pathp->flags = 0;

  pathp->flags |= PATHW_CLAMP_CONC;

  pathp->km = conc;
  pathp->km_filename = NULL;

  pathp->orientation1 = 1;
  pathp->orientation3 = 0;
  if (orient == 0) {
    pathp->orientation2 = 0;
  } else {
    pathp->orientation2 = (orient < 0) ? -1 : 1;
  }

  pathp->product_head = NULL;
  pathp->prod_signature = NULL;

  pathp->next = rxnp->pathway_head;
  rxnp->pathway_head = pathp;

  no = CHECKED_MALLOC_STRUCT(struct name_orient, "struct name_orient");
  no->name = CHECKED_STRDUP(mol_sym->name, "molecule name");
  no->orient = pathp->orientation2;

  if (surface_class->clamp_conc_mols == NULL) {
    no->next = NULL;
    surface_class->clamp_conc_mols = no;
  } else {
    no->next = surface_class->clamp_conc_mols;
    surface_class->clamp_conc_mols = no;
  }

  return MCELL_SUCCESS;
}

/************************************************************************
 *
 * function for changing the reaction rate constant of a given named
 * reaction.
 *
 * The call expects:
 *
 * - MCELL_STATE
 * - reaction name: const char* containing the name of reaction
 * - new rate: a double with the new reaction rate constant
 *
 * NOTE: This function can be called anytime after the
 *       REACTION_DATA_OUTPUT has been either parsed or
 *       set up with API calls.
 *
 * Returns 1 on error and 0 on success
 *
 ************************************************************************/
MCELL_STATUS
mcell_change_reaction_rate(MCELL_STATE *state, const char *reaction_name,
                           double new_rate) {
  // sanity check
  if (new_rate < 0.0) {
    return MCELL_FAIL;
  }

  // retrive reaction corresponding to name if it exists
  struct rxn *rx = NULL;
  int path_id = 0;
  if (get_rxn_by_name(state->reaction_hash, state->rx_hashsize, reaction_name,
                      &rx, &path_id)) {
    return MCELL_FAIL;
  }

  // now change the rate
  if (change_reaction_probability(
      &state->reaction_prob_limit_flag, state->notify, rx, path_id,
      new_rate)) {
    return MCELL_FAIL;
  }

  return MCELL_SUCCESS;
}

/*******************************************************************************
 *
 * static helper functions
 *
 ******************************************************************************/

/*************************************************************************
 init_reactions:
    Postprocess the parsed reactions, moving them to the reaction hash table,
    and transferring information from the pathway structures to a more compact,
    runtime-optimized form.

 In: state: simulation state
 Out: Returns 1 on error, 0 on success.
      Reaction hash table is built and geometries are set properly.  Unlike in
      the parser, reactions with different reactant geometries are _different
      reactions_, and are stored as separate struct rxns.

 Note: The user inputs _geometric equivalence classes_, but here we convert
       from that to _output construction geometry_.  A geometry of 0 means to
       choose a random orientation.  A geometry of k means to adopt the
       geometry of the k'th species in the list (reactants start at #1,
       products are in order after reactants).  A geometry of -k means to adopt
       the opposite of the geometry of the k'th species.  The first n_reactants
       products determine the fate of the reactants (NULL = destroyed), and the
       rest are real products.
 PostNote: The reactants are used for triggering, and those have equivalence
       class geometry even in here.
*************************************************************************/
int init_reactions(MCELL_STATE *state) {
  struct pathway *path;
  struct product *prod = NULL;
  short geom;
  int num_rx = 0;

  state->vacancy_search_dist2 *= state->r_length_unit; /* Convert units */
  state->vacancy_search_dist2 *= state->vacancy_search_dist2; /* Take square */

  if (state->rx_radius_3d <= 0.0) {
    state->rx_radius_3d = 1.0 / sqrt(MY_PI * state->grid_density);
  }
  state->tv_rxn_mem = create_mem(sizeof(struct t_func), 100);
  if (state->tv_rxn_mem == NULL)
    return 1;

  for (int n_rxn_bin = 0; n_rxn_bin < state->rxn_sym_table->n_bins;
       n_rxn_bin++) {
    for (struct sym_entry *sym = state->rxn_sym_table->entries[n_rxn_bin];
         sym != NULL; sym = sym->next) {
      struct rxn *reaction = (struct rxn *)sym->value;
      reaction->next = NULL;

      for (path = reaction->pathway_head; path != NULL; path = path->next) {
        check_duplicate_special_reactions(path);

        /* if one of the reactants is a surface, move it to the last reactant.
         * Also arrange reactant1 and reactant2 in alphabetical order */
        if (reaction->n_reactants > 1) {
          struct species *temp_sp;
          /* Put surface last */
          if ((path->reactant1->flags & IS_SURFACE) != 0) {
            temp_sp = path->reactant1;
            path->reactant1 = path->reactant2;
            path->reactant2 = temp_sp;
            geom = path->orientation1;
            path->orientation1 = path->orientation2;
            path->orientation2 = geom;
          }
          if (reaction->n_reactants > 2) {
            if ((path->reactant2->flags & IS_SURFACE) != 0) {
              temp_sp = path->reactant3;
              path->reactant3 = path->reactant2;
              path->reactant2 = temp_sp;
              geom = path->orientation3;
              path->orientation3 = path->orientation2;
              path->orientation2 = geom;
            }
          }
          alphabetize_pathway(path, reaction);
        } /* end if (n_reactants > 1) */

      } /* end for (path = reaction->pathway_head; ...) */

      /* if reaction contains equivalent pathways, split this reaction into a
       * linked list of reactions each containing only equivalent pathways. */

      struct rxn *rx = split_reaction(reaction);

      /* set the symbol value to the head of the linked list of reactions */
      sym->value = (void *)rx;

      while (rx != NULL) {
        double pb_factor = 0.0;
        /* Check whether reaction contains pathways with equivalent product
         * lists.  Also sort pathways in alphabetical order according to the
         * "prod_signature" field. */
        check_reaction_for_duplicate_pathways(&rx->pathway_head);

        num_rx++;

        /* At this point we have reactions of the same geometry and can
         * collapse them and count how many non-reactant products are in each
         * pathway. */

        /* Search for reactants that appear as products */
        /* Any reactants that don't appear are set to be destroyed. */
        rx->product_idx = CHECKED_MALLOC_ARRAY(u_int, rx->n_pathways + 1,
                                               "reaction product index array");
        rx->cum_probs = CHECKED_MALLOC_ARRAY(
            double, rx->n_pathways, "reaction cumulative probabilities array");

        /* Note, that the last member of the array "rx->product_idx" contains
         * size of the array "rx->players" */

        if (rx->product_idx == NULL || rx->cum_probs == NULL)
          return 1;

        int n_prob_t_rxns = 0; /* # of pathways with time-varying rates */
        path = rx->pathway_head;

        for (int n_pathway = 0; path != NULL; n_pathway++, path = path->next) {

          rx->product_idx[n_pathway] = 0;

          /* Look for concentration clamp */
          if (path->reactant2 != NULL &&
              (path->reactant2->flags & IS_SURFACE) != 0 && path->km >= 0.0 &&
              path->product_head == NULL &&
              ((path->flags & PATHW_CLAMP_CONC) != 0)) {
            struct ccn_clamp_data *ccd;

            if (n_pathway != 0 || path->next != NULL)
              mcell_warn("Mixing surface modes with other surface reactions.  "
                         "Please don't.");

            if (path->km > 0) {
              ccd = CHECKED_MALLOC_STRUCT(struct ccn_clamp_data,
                                          "concentration clamp data");
              if (ccd == NULL)
                return 1;

              ccd->surf_class = path->reactant2;
              ccd->mol = path->reactant1;
              ccd->concentration = path->km;
              if (path->orientation1 * path->orientation2 == 0) {
                ccd->orient = 0;
              } else {
                ccd->orient =
                    (path->orientation1 == path->orientation2) ? 1 : -1;
              }
              ccd->sides = NULL;
              ccd->next_mol = NULL;
              ccd->next_obj = NULL;
              ccd->objp = NULL;
              ccd->n_sides = 0;
              ccd->side_idx = NULL;
              ccd->cum_area = NULL;
              ccd->scaling_factor = 0.0;
              ccd->next = state->clamp_list;
              state->clamp_list = ccd;
            }
            path->km = GIGANTIC;
          } else if ((path->flags & PATHW_TRANSP) != 0) {
            rx->n_pathways = RX_TRANSP;
            if (path->reactant2 != NULL &&
                (path->reactant2->flags & IS_SURFACE) &&
                (path->reactant1->flags & ON_GRID)) {
              path->reactant1->flags |= CAN_REGION_BORDER;
            }
          } else if ((path->flags & PATHW_REFLEC) != 0) {
            rx->n_pathways = RX_REFLEC;
            if (path->reactant2 != NULL &&
                (path->reactant2->flags & IS_SURFACE) &&
                (path->reactant1->flags & ON_GRID)) {
              path->reactant1->flags |= CAN_REGION_BORDER;
            }
          } else if (path->reactant2 != NULL &&
                     (path->reactant2->flags & IS_SURFACE) &&
                     (path->reactant1->flags & ON_GRID) &&
                     (path->product_head == NULL) &&
                     (path->flags & PATHW_ABSORP)) {
            rx->n_pathways = RX_ABSORB_REGION_BORDER;
            path->reactant1->flags |= CAN_REGION_BORDER;
          } else if ((strcmp(path->reactant1->sym->name,
                             "ALL_SURFACE_MOLECULES") == 0)) {
            if (path->reactant2 != NULL &&
                (path->reactant2->flags & IS_SURFACE) &&
                (path->product_head == NULL) && (path->flags & PATHW_ABSORP)) {
              rx->n_pathways = RX_ABSORB_REGION_BORDER;
              path->reactant1->flags |= CAN_REGION_BORDER;
            }
          }
          if (path->km_filename == NULL)
            rx->cum_probs[n_pathway] = path->km;
          else {
            rx->cum_probs[n_pathway] = 0;
            n_prob_t_rxns++;
          }

          /* flags that tell whether reactant_1 is also on the product list,
             same for reactant_2 and reactant_3 */
          int recycled1 = 0;
          int recycled2 = 0;
          int recycled3 = 0;

          for (prod = path->product_head; prod != NULL; prod = prod->next) {
            if (recycled1 == 0 && prod->prod == path->reactant1)
              recycled1 = 1;
            else if (recycled2 == 0 && prod->prod == path->reactant2)
              recycled2 = 1;
            else if (recycled3 == 0 && prod->prod == path->reactant3)
              recycled3 = 1;
            else
              rx->product_idx[n_pathway]++;
          }

        } /* end for (n_pathway=0,path=rx->pathway_head; ...) */

        /* Now that we know how many products there really are, set the index
         * array and malloc space for the products and geometries. */
        int num_players = rx->n_reactants;
        int kk = rx->n_pathways;
        if (kk <= RX_SPECIAL)
          kk = 1;
        for (int n_pathway = 0; n_pathway < kk; n_pathway++) {
          int k = rx->product_idx[n_pathway] + rx->n_reactants;
          rx->product_idx[n_pathway] = num_players;
          num_players += k;
        }
        rx->product_idx[kk] = num_players;

        rx->players = CHECKED_MALLOC_ARRAY(struct species *, num_players,
                                           "reaction players array");
        rx->geometries = CHECKED_MALLOC_ARRAY(short, num_players,
                                              "reaction geometries array");

        if (rx->players == NULL || rx->geometries == NULL)
          return 1;

        /* Load all the time-varying rates from disk (if any), merge them into
         * a single sorted list, and pull off any updates for time zero. */
        if (n_prob_t_rxns > 0) {
          path = rx->pathway_head;
          for (int n_pathway = 0; path != NULL;
               n_pathway++, path = path->next) {
            if (path->km_filename != NULL) {
              if (load_rate_file(state->time_unit, state->tv_rxn_mem, rx,
                                 path->km_filename, n_pathway, state->notify->neg_reaction))
                mcell_error("Failed to load rates from file '%s'.",
                            path->km_filename);
            }
            free(path->km_filename);
            path->km_filename = NULL;
          }
          rx->prob_t = (struct t_func *)ae_list_sort(
              (struct abstract_element *)rx->prob_t);

          while (rx->prob_t != NULL && rx->prob_t->time <= 0.0) {
            rx->cum_probs[rx->prob_t->path] = rx->prob_t->value;
            rx->prob_t = rx->prob_t->next;
          }
        } /* end if (n_prob_t_rxns > 0) */

        /* Set the geometry of the reactants.  These are used for triggering. */
        /* Since we use flags to control orientation changes, just tell everyone
         * to stay put. */
        path = rx->pathway_head;
        rx->players[0] = path->reactant1;
        rx->geometries[0] = path->orientation1;
        if (rx->n_reactants > 1) {
          rx->players[1] = path->reactant2;
          rx->geometries[1] = path->orientation2;
          if (rx->n_reactants > 2) {
            rx->players[2] = path->reactant3;
            rx->geometries[2] = path->orientation3;
          }
        }

        /* maximum number of surface products */
        path = rx->pathway_head;
        int max_num_surf_products = set_product_geometries(path, rx, prod);

        pb_factor = compute_pb_factor(
            state->time_unit, state->length_unit, state->grid_density,
            state->rx_radius_3d,
            &state->rxn_flags,
            &state->create_shared_walls_info_flag,
            rx, max_num_surf_products);
        rx->pb_factor = pb_factor;
        path = rx->pathway_head;

        if (scale_probabilities(&state->reaction_prob_limit_flag, state->notify,
                                path, rx, pb_factor))
          return 1;

        if (n_prob_t_rxns > 0) {
          for (struct t_func *tp = rx->prob_t; tp != NULL; tp = tp->next)
            tp->value *= pb_factor;
        }

        /* Move counts from list into array */
        if (rx->n_pathways > 0) {
          rx->info = CHECKED_MALLOC_ARRAY(struct pathway_info, rx->n_pathways,
                                          "reaction pathway info");
          if (rx->info == NULL)
            return 1;

          path = rx->pathway_head;
          for (int n_pathway = 0; path != NULL;
               n_pathway++, path = path->next) {
            rx->info[n_pathway].count = 0;
            rx->info[n_pathway].pathname =
                path->pathname; /* Keep track of named rxns */
            if (path->pathname != NULL) {
              rx->info[n_pathway].pathname->path_num = n_pathway;
              rx->info[n_pathway].pathname->rx = rx;
            }
          }
        } else /* Special reaction, only one exit pathway */
        {
          rx->info = CHECKED_MALLOC_STRUCT(struct pathway_info,
                                           "reaction pathway info");
          if (rx->info == NULL)
            return 1;
          rx->info[0].count = 0;
          rx->info[0].pathname = rx->pathway_head->pathname;
          if (rx->pathway_head->pathname != NULL) {
            rx->info[0].pathname->path_num = 0;
            rx->info[0].pathname->rx = rx;
          }
        }

        /* Compute cumulative properties */
        for (int n_pathway = 1; n_pathway < rx->n_pathways; ++n_pathway)
          rx->cum_probs[n_pathway] += rx->cum_probs[n_pathway - 1];
        if (rx->n_pathways > 0)
          rx->min_noreaction_p = rx->max_fixed_p =
              rx->cum_probs[rx->n_pathways - 1];
        else
          rx->min_noreaction_p = rx->max_fixed_p = 1.0;

        rx = rx->next;
      }
    }
  }

  if (state->rxn_flags.surf_surf_reaction_flag ||
      state->rxn_flags.surf_surf_surf_reaction_flag) {
    if (state->notify->reaction_probabilities == NOTIFY_FULL)
      mcell_log("For reaction between two (or three) surface molecules the "
                "upper probability limit is given. The effective reaction "
                "probability will be recalculated dynamically during "
                "simulation.");
  }

  if (build_reaction_hash_table(&state->reaction_hash, &state->n_reactions,
                                state->rxn_sym_table, &state->rx_hashsize,
                                num_rx))
    return 1;

  state->rx_radius_3d *= state->r_length_unit; /* Convert into length units */

  for (int n_rxn_bin = 0; n_rxn_bin < state->rx_hashsize; n_rxn_bin++) {
    for (struct rxn *this_rx = state->reaction_hash[n_rxn_bin]; this_rx != NULL;
         this_rx = this_rx->next) {
      /* Here we deallocate all memory used for creating pathways. */
      path = this_rx->pathway_head;
      struct pathway *next_path = path;
      while (next_path != NULL) {
        next_path = path->next;
        if (path->prod_signature != NULL) {
          free(path->prod_signature);
        }

        struct product *dead_prod = path->product_head;
        struct product *nxt = dead_prod;
        while (nxt != NULL) {
          nxt = dead_prod->next;
          free(dead_prod);
          dead_prod = nxt;
        }

        free(path);
        path = next_path;
      }

      set_reaction_player_flags(this_rx);
      this_rx->pathway_head = NULL;
    }
  }

  add_surface_reaction_flags(state->mol_sym_table, state->all_mols, state->all_surface_mols,
                             state->all_volume_mols);

  if (state->notify->reaction_probabilities == NOTIFY_FULL)
    mcell_log_raw("\n");

  return 0;
}

/*******************************************************************************
 *
 * static helper functions
 *
 ******************************************************************************/
/*************************************************************************
 *
 * extract_reactants extracts the reactant info into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS
extract_reactants(struct pathway *pathp, struct mcell_species *reactants,
                  int *num_reactants, int *num_vol_mols, int *num_surface_mols,
                  int *all_3d, int *oriented_count) {
  int reactant_idx = 0;
  struct mcell_species *current_reactant;
  for (current_reactant = reactants;
       reactant_idx < 3 && current_reactant != NULL;
       ++reactant_idx, current_reactant = current_reactant->next) {
    /* Extract orientation and species */
    short orient = current_reactant->orient_set ? current_reactant->orient : 0;
    struct species *reactant_species =
        (struct species *)current_reactant->mol_type->value;

    /* Count the type of this reactant */
    if (current_reactant->orient_set) {
      ++(*oriented_count);
    }

    if (reactant_species->flags & NOT_FREE) {
      *all_3d = 0;
      if (reactant_species->flags & ON_GRID) {
        ++(*num_surface_mols);
      }
    } else {
      ++(*num_vol_mols);
    }

    /* Sanity check this reactant */
    if (reactant_species->flags & IS_SURFACE) {
      mcell_error("surface class can be listed only as the last reactant on "
                  "the left-hand side of the reaction with the preceding '@' "
                  "sign.");
      return MCELL_FAIL;
    }

    /* Copy in reactant info */
    switch (reactant_idx) {
    case 0:
      pathp->reactant1 = reactant_species;
      pathp->orientation1 = orient;
      break;

    case 1:
      pathp->reactant2 = reactant_species;
      pathp->orientation2 = orient;
      break;

    case 2:
      pathp->reactant3 = reactant_species;
      pathp->orientation3 = orient;
      break;

    default:
      UNHANDLED_CASE(reactant_idx);
    }
  }
  *num_reactants = reactant_idx;

  /* we had more than 3 reactants */
  if (current_reactant != NULL) {
    return MCELL_FAIL;
  } else {
    return MCELL_SUCCESS;
  }
}

/*************************************************************************
 *
 * extract_catalytic_arrow extracts the info for a catalytic arrow
 * into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS
extract_catalytic_arrow(struct pathway *pathp,
                        struct reaction_arrow *react_arrow, int *reactant_idx,
                        int *num_vol_mols, int *num_surface_mols, int *all_3d,
                        int *oriented_count) {
  struct species *catalyst_species =
      (struct species *)react_arrow->catalyst.mol_type->value;
  short orient =
      react_arrow->catalyst.orient_set ? react_arrow->catalyst.orient : 0;

  /* XXX: Should surface class be allowed inside a catalytic arrow? */
  if (catalyst_species->flags & IS_SURFACE) {
     mcell_error("a surface class may not appear inside a catalytic arrow");
    return MCELL_FAIL;
  }

  /* Count the type of this reactant */
  if (react_arrow->catalyst.orient_set) {
    ++(*oriented_count);
  }

  if (catalyst_species->flags & NOT_FREE) {
    *all_3d = 0;
    if (catalyst_species->flags & ON_GRID) {
      ++(*num_surface_mols);
    }
  } else {
    ++(*num_vol_mols);
  }

  /* Copy in catalytic reactant */
  switch (*reactant_idx) {
  case 1:
    pathp->reactant2 = (struct species *)react_arrow->catalyst.mol_type->value;
    pathp->orientation2 = orient;
    break;

  case 2:
    pathp->reactant3 = (struct species *)react_arrow->catalyst.mol_type->value;
    pathp->orientation3 = orient;
    break;

  case 0:
  default:
    // mcell_internal_error("Catalytic reagent ended up in an invalid slot
    // (%d).", reactant_idx);
    return MCELL_FAIL;
  }
  ++(*reactant_idx);

  return MCELL_SUCCESS;
}

/*************************************************************************
 *
 * extract_surface extracts the info for a surface included in the
 * reaction specification
 *
 *************************************************************************/
MCELL_STATUS
extract_surface(struct pathway *path, struct mcell_species *surf_class,
                int *num_reactants, unsigned int *num_surfaces,
                int *oriented_count) {
  short orient = surf_class->orient_set ? surf_class->orient : 0;
  if (surf_class->orient_set) {
    (*oriented_count)++;
  }

  /* Copy reactant into next available slot */
  switch (*num_reactants) {
  case 0:
    // mdlerror(parse_state, "Before defining reaction surface class at least
    // one reactant should be defined.");
    return MCELL_FAIL;

  case 1:
    path->reactant2 = (struct species *)surf_class->mol_type->value;
    path->orientation2 = orient;
    break;

  case 2:
    path->reactant3 = (struct species *)surf_class->mol_type->value;
    path->orientation3 = orient;
    break;

  default:
    // mdlerror(parse_state, "Too many reactants--maximum number is two plus
    // reaction surface class.");
    return MCELL_FAIL;
  }

  (*num_reactants)++;
  (*num_surfaces)++;

  return MCELL_SUCCESS;
}

/*************************************************************************
 *
 * check_surface_specs performs a number of sanity checks to make sure
 * the surface specifications are sane
 *
 *************************************************************************/
MCELL_STATUS
check_surface_specs(struct notifications *notify, int num_reactants,
                    int num_surfaces, int num_vol_mols, int all_3d,
                    int oriented_count) {
  if (num_surfaces > 1) {
    /* Shouldn't happen */
    mcell_internal_error(
        "Too many surfaces--reactions can take place on at most one surface.");
    return MCELL_FAIL;
  }

  if (num_surfaces == num_reactants) {
    mcell_error("Reactants cannot consist entirely of surfaces.  Use a surface "
                "release site instead!");
    return MCELL_FAIL;
  }

  if ((num_vol_mols == 2) && (num_surfaces == 1)) {
    mcell_error(
        "Reaction between two volume molecules and a surface is not defined.");
    return MCELL_FAIL;
  }

  if (all_3d) {
    if (oriented_count != 0) {
      if (notify->useless_vol_orient == WARN_ERROR) {
        mcell_error("Orientation specified for molecule in reaction in volume");
        return MCELL_FAIL;
      } else if (notify->useless_vol_orient == WARN_WARN) {
        mcell_warn("Orientation specified for molecule in reaction in volume");
      }
    }
  } else {
    if (num_reactants != oriented_count) {
      if (notify->missed_surf_orient == WARN_ERROR) {
        mcell_error("Orientation not specified for molecule in reaction "
                    "at surface\n  (use ; or ', or ,' for random orientation)");
        return MCELL_FAIL;
      } else if (notify->missed_surf_orient == WARN_WARN) {
        mcell_warn("Orientation not specified for molecule in reaction at "
                   "surface\n  (use ; or ', or ,' for random orientation)");
      }
    }
  }

  return MCELL_SUCCESS;
}

/*************************************************************************
 *
 * add_catalytic_species_to_products adds all species that are part of a
 * catalytic reaction to the list of products.
 *
 *************************************************************************/
MCELL_STATUS
add_catalytic_species_to_products(struct pathway *path, int catalytic,
                                  int bidirectional, int all_3d) {
  struct species *catalyst;
  short catalyst_orient;
  switch (catalytic) {
  case 0:
    catalyst = path->reactant1;
    catalyst_orient = path->orientation1;
    break;
  case 1:
    catalyst = path->reactant2;
    catalyst_orient = path->orientation2;
    break;
  case 2:
    catalyst = path->reactant3;
    catalyst_orient = path->orientation3;
    break;
  default:
    mcell_internal_error("Catalytic reagent index is invalid.");
    return MCELL_FAIL;
  }

  if (bidirectional || !(catalyst->flags & IS_SURFACE)) {
    struct product *prodp = (struct product *)CHECKED_MALLOC_STRUCT(
        struct product, "reaction product");
    if (prodp == NULL) {
      return MCELL_FAIL;
    }

    prodp->prod = catalyst;
    if (all_3d) {
      prodp->orientation = 0;
    } else {
      prodp->orientation = catalyst_orient;
    }
    prodp->next = path->product_head;
    path->product_head = prodp;
  }

  return MCELL_SUCCESS;
}

/*************************************************************************
 *
 * extract_products extracts the product info into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS
extract_products(struct notifications *notify, struct pathway *pathp,
                 struct mcell_species *products, int *num_surf_products,
                 int bidirectional,
                 int all_3d) {
  struct mcell_species *current_product;
  for (current_product = products; current_product != NULL;
       current_product = current_product->next) {
    /* Nothing to do for NO_SPECIES */
    if (current_product->mol_type == NULL)
      continue;

    /* Create new product */
    struct product *prodp = (struct product *)CHECKED_MALLOC_STRUCT(
        struct product, "reaction product");

    if (prodp == NULL) {
      // mcell_error_raw("Out of memory while creating reaction: %s -> ... ",
      //                rxnp->sym->name);
      return MCELL_FAIL;
    }

    /* Set product species and orientation */
    prodp->prod = (struct species *)current_product->mol_type->value;
    if (all_3d) {
      prodp->orientation = 0;
    } else {
      prodp->orientation = current_product->orient;
    }

    /* Disallow surface as product unless reaction is bidirectional */
    if (!bidirectional) {
      if (prodp->prod->flags & IS_SURFACE) {
        mcell_error_raw("Surface_class '%s' is not allowed to be on the "
                        "product side of the reaction.",
                        prodp->prod->sym->name);
        free(prodp);
        return MCELL_FAIL;
      }
    }

    /* Append product to list */
    prodp->next = pathp->product_head;
    pathp->product_head = prodp;

    if (prodp->prod->flags & ON_GRID) {
      ++(*num_surf_products);
    }

    /* Add product if it isn't a surface */
    if (!(prodp->prod->flags & IS_SURFACE)) {
      if (all_3d == 0) {
        if (!current_product->orient_set) {
          if (notify->missed_surf_orient == WARN_ERROR) {
            mcell_error("Product orientation not specified for molecule in "
                        "reaction at surface\n  (use ; or ', or ,' for random "
                        "orientation)");
            return MCELL_FAIL;
          } else if (notify->missed_surf_orient == WARN_WARN) {
            mcell_warn("Product orientation not specified for molecule in "
                       "reaction at surface\n  (use ; or ', or ,' for random "
                       "orientation)");
          }
        }
      } else {
        if ((prodp->prod->flags & NOT_FREE) != 0) {
          mcell_error("Reaction has only volume reactants but is trying to "
                      "create a surface product");
          return MCELL_FAIL;
        }
        if (current_product->orient_set) {
          if (notify->useless_vol_orient == WARN_ERROR) {
            mcell_error("Orientation specified for molecule in reaction in "
                        "volume");
            return MCELL_FAIL;
          } else if (notify->useless_vol_orient == WARN_WARN) {
            mcell_warn("Orientation specified for molecule in reaction in "
                       "volume");
          }
        }
      }
    }
  }

  return MCELL_SUCCESS;
}

/*************************************************************************
 create_rx_name:
    Assemble reactants alphabetically into a reaction name string.

 In:  p: reaction pathway whose reaction name we are to create
 Out: a string to be used as a symbol name for the reaction
*************************************************************************/
char *create_rx_name(struct pathway *p) {

  struct species *reagents[3];
  int n_reagents = 0;

  /* Store reagents in an array. */
  reagents[0] = p->reactant1;
  reagents[1] = p->reactant2;
  reagents[2] = p->reactant3;

  /* Count non-null reagents. */
  for (n_reagents = 0; n_reagents < 3; ++n_reagents)
    if (reagents[n_reagents] == NULL)
      break;

  /* Sort reagents. */
  for (int i = 0; i < n_reagents; ++i) {
    for (int j = i + 1; j < n_reagents; ++j) {
      /* If 'j' precedes 'i', 'j' wins. */
      if (strcmp(reagents[j]->sym->name, reagents[i]->sym->name) < 0) {
        struct species *tmp = reagents[j];
        reagents[j] = reagents[i];
        reagents[i] = tmp;
      }
    }
  }

  /* Now, produce a name! */
  switch (n_reagents) {
  case 1:
    return alloc_sprintf("%s", reagents[0]->sym->name);
  case 2:
    return alloc_sprintf("%s+%s", reagents[0]->sym->name,
                         reagents[1]->sym->name);
  case 3:
    return alloc_sprintf("%s+%s+%s", reagents[0]->sym->name,
                         reagents[1]->sym->name, reagents[2]->sym->name);
  default:
    // mcell_internal_error("Invalid number of reagents in reaction pathway
    // (%d).", n_reagents);
    return NULL;
  }
}

/*************************************************************************
 concat_rx_name:
    Concatenates reactants onto a reaction name.

 In:  name1: name of first reactant (or first part of reaction name)
      name2: name of second reactant (or second part of reaction name)
 Out: reaction name as a string, or NULL if an error occurred
*************************************************************************/
static char *concat_rx_name(char *name1, char *name2) {
  char *rx_name;

  /* Sort them */
  if (strcmp(name2, name1) <= 0) {
    char *nametmp = name1;
    name1 = name2;
    name2 = nametmp;
  }

  /* Build the name */
  rx_name = CHECKED_SPRINTF("%s+%s", name1, name2);

  /* Die if we failed to allocate memory */
  if (rx_name == NULL)
    return NULL;

  return rx_name;
}

/***********************************************************************
 invert_current_reaction_pathway:
    Creates a new reversed pathway, where the reactants of new pathway are the
    products of the current pathway, and the products of new pathway are the
    reactants of the current pathway.

 In:  rxn_sym_table:
      vacancy_search_dist2:
      pathp: pathway to invert
      reverse_rate: the reverse reaction rate
      rate_filename:
 Out: Returns 1 on error and 0 - on success.  The new pathway is added to the
      linked list of the pathways for the current reaction.
***********************************************************************/
MCELL_STATUS invert_current_reaction_pathway(
    struct sym_table_head *rxn_sym_table,
    double vacancy_search_dist2,
    struct pathway *pathp,
    struct reaction_rate *reverse_rate,
    const char *rate_filename) {

  struct product *prodp;
  int num_surf_products = 0;
  int num_surface_mols = 0;
  int num_vol_mols = 0;

  /* flag that tells whether there is a surface_class
     among products in the direct reaction */
  int is_surf_class = 0;

  int all_3d = 1; // flag that tells whether all products are volume_molecules
  int nprods; /* number of products */
  for (nprods = 0, prodp = pathp->product_head; prodp != NULL;
       prodp = prodp->next) {
    nprods++;
    if ((prodp->prod->flags & NOT_FREE) != 0)
      all_3d = 0;
    if ((prodp->prod->flags & IS_SURFACE) != 0) {
      is_surf_class = 1;
    }
  }

  if (nprods == 0) {
    // mdlerror(parse_state, "Can't create a reverse reaction with no
    // products");
    return MCELL_FAIL;
  }
  if (nprods == 1 && (pathp->product_head->prod->flags & IS_SURFACE)) {
    // mdlerror(parse_state, "Can't create a reverse reaction starting from only
    // a surface");
    return MCELL_FAIL;
  }
  if (nprods > 3) {
    // mdlerror(parse_state, "Can't create a reverse reaction involving more
    // than three products. Please note that surface_class from the reaction
    // reactant side also counts as a product.");
    return MCELL_FAIL;
  }

  if (pathp->pathname != NULL) {
    // mdlerror(parse_state, "Can't name bidirectional reactions--write each
    // reaction and name them separately");
    return MCELL_FAIL;
  }
  if (all_3d) {
    if ((pathp->reactant1->flags & NOT_FREE) != 0)
      all_3d = 0;
    if (pathp->reactant2 != NULL && (pathp->reactant2->flags & NOT_FREE) != 0)
      all_3d = 0;
    if (pathp->reactant3 != NULL && (pathp->reactant3->flags & NOT_FREE) != 0)
      all_3d = 0;

    if (!all_3d) {
      // mdlerror(parse_state, "Cannot reverse orientable reaction with only
      // volume products");
      return MCELL_FAIL;
    }
  }

  prodp = pathp->product_head;
  char *inverse_name;
  if (nprods == 1) {
    inverse_name = strdup(prodp->prod->sym->name);

    if (inverse_name == NULL)
      return MCELL_FAIL;
  } else if (nprods == 2) {
    inverse_name =
        concat_rx_name(prodp->prod->sym->name, prodp->next->prod->sym->name);
  } else {
    char *tmp_inverse_name = concat_rx_name(
      prodp->prod->sym->name, prodp->next->prod->sym->name);
    if (tmp_inverse_name == NULL) {
      return MCELL_FAIL;
    }
    inverse_name =
        concat_rx_name(tmp_inverse_name, prodp->next->next->prod->sym->name);
    free(tmp_inverse_name);
  }
  if (inverse_name == NULL) {
    return MCELL_FAIL;
  }

  struct sym_entry *sym = retrieve_sym(inverse_name, rxn_sym_table);
  if (sym == NULL) {
    sym = store_sym(inverse_name, RX, rxn_sym_table, NULL);
    if (sym == NULL) {
      // mdlerror_fmt(parse_state, "File '%s', Line %ld: Out of memory while
      // storing reaction pathway.", __FILE__, (long)__LINE__);
      free(inverse_name);
      return MCELL_FAIL;
    }
  }
  free(inverse_name);
  struct rxn *rx = (struct rxn *)sym->value;
  rx->n_reactants = nprods;
  rx->n_pathways++;

  struct pathway *path = (struct pathway *)CHECKED_MALLOC_STRUCT(
    struct pathway, "reaction pathway");
  if (path == NULL) {
    return MCELL_FAIL;
  }
  path->pathname = NULL;
  path->flags = 0;
  path->reactant1 = prodp->prod;
  if ((path->reactant1->flags & NOT_FREE) == 0) {
    ++num_vol_mols;
  } else {
    if (path->reactant1->flags & ON_GRID) {
      ++num_surface_mols;
    }
  }
  path->orientation1 = prodp->orientation;
  path->reactant2 = NULL;
  path->reactant3 = NULL;
  path->prod_signature = NULL;
  if (nprods > 1) {
    path->reactant2 = prodp->next->prod;
    if ((path->reactant2->flags & NOT_FREE) == 0) {
      ++num_vol_mols;
    } else {
      if (path->reactant2->flags & ON_GRID) {
        ++num_surface_mols;
      }
    }
    path->orientation2 = prodp->next->orientation;
  }
  if (nprods > 2) {
    path->reactant3 = prodp->next->next->prod;
    if ((path->reactant3->flags & NOT_FREE) == 0) {
      ++num_vol_mols;
    } else {
      if (path->reactant3->flags & ON_GRID) {
        ++num_surface_mols;
      }
    }
    path->orientation3 = prodp->next->next->orientation;
  }

  switch (reverse_rate->rate_type) {
  case RATE_UNSET:
    // mdlerror_fmt(parse_state, "File %s, Line %d: Internal error: Reverse rate
    // is not set", __FILE__, __LINE__);
    free(path);
    return MCELL_FAIL;

  case RATE_CONSTANT:
    path->km = reverse_rate->v.rate_constant;
    path->km_filename = NULL;
    break;

  case RATE_FILE:
    path->km = 0.0;
    path->km_filename = (char *)rate_filename;
    free(reverse_rate->v.rate_file);
    reverse_rate->v.rate_file = NULL;
    break;

  default:
    UNHANDLED_CASE(reverse_rate->rate_type);
  }

  path->product_head = (struct product *)CHECKED_MALLOC_STRUCT(
      struct product, "reaction product");
  if (path->product_head == NULL) {
    free(path);
    return 1;
  }

  path->product_head->orientation = pathp->orientation1;
  path->product_head->prod = pathp->reactant1;
  path->product_head->next = NULL;
  if (path->product_head->prod->flags & ON_GRID)
    ++num_surf_products;

  if ((pathp->reactant2 != NULL) &&
      ((pathp->reactant2->flags & IS_SURFACE) == 0)) {
    path->product_head->next = (struct product *)CHECKED_MALLOC_STRUCT(
        struct product, "reaction product");
    if (path->product_head->next == NULL) {
      free(path);
      return 1;
    }
    path->product_head->next->orientation = pathp->orientation2;
    path->product_head->next->prod = pathp->reactant2;
    path->product_head->next->next = NULL;
    if (path->product_head->next->prod->flags & ON_GRID)
      ++num_surf_products;

    if ((pathp->reactant3 != NULL) &&
        ((pathp->reactant3->flags & IS_SURFACE) == 0)) {
      path->product_head->next->next = (struct product *)CHECKED_MALLOC_STRUCT(
          struct product, "reaction product");
      if (path->product_head->next->next == NULL) {
        free(path);
        return 1;
      }
      path->product_head->next->next->orientation = pathp->orientation3;
      path->product_head->next->next->prod = pathp->reactant3;
      path->product_head->next->next->next = NULL;
      if (path->product_head->next->next->prod->flags & ON_GRID)
        ++num_surf_products;
    }
  }

  path->prod_signature = create_prod_signature(&path->product_head);
  if (path->prod_signature == NULL) {
    // mdlerror(parse_state, "Error creating 'prod_signature' field for reaction
    // pathway.");
    free(path);
    return MCELL_FAIL;
  }

  if ((!distinguishable(vacancy_search_dist2, 0, EPS_C)) &&
      (num_surf_products > num_surface_mols)) {
    /* the case with one volume molecule reacting with the surface
       and producing one surface molecule is excluded */
    if (!((num_surface_mols == 0) && (num_vol_mols == 1))) {
      // mdlerror(parse_state, "Error: number of surface products exceeds number
      // of surface reactants, but VACANCY_SEARCH_DISTANCE is not specified or
      // set to zero.");
      free(path);
      return MCELL_FAIL;
    }
  }

  /* Now go back to the original reaction and if there is a "surface_class"
     among products - remove it.  We do not need it now on the product side
     of the reaction */
  if (is_surf_class) {
    prodp = pathp->product_head;
    if (prodp->prod->flags & IS_SURFACE) {
      pathp->product_head = prodp->next;
      prodp->next = NULL;
      // mem_put(parse_state->prod_mem, (void *)prodp);
    } else if (prodp->next->prod->flags & IS_SURFACE) {
      // struct product *temp = prodp->next;
      prodp->next = prodp->next->next;
      // mem_put(parse_state->prod_mem, temp);
    } else {
      // struct product *temp = prodp->next->next;
      prodp->next->next = prodp->next->next->next;
      // mem_put(parse_state->prod_mem, temp);
    }
  }

  path->next = rx->pathway_head;
  rx->pathway_head = path;
  return 0;
}

/************************************************************************
 * static helper functions
 ************************************************************************/

/*************************************************************************
 sort_product_list_compare:
    Comparison function for products to be sorted when generating the product
    signature.

 In:  list_item: first item to compare
      new_item:  second item to compare
 Out: -1 if list_item < new_item, 1 if list_item > new_item, 0 if they are
      equal

  XXX Currently this function also appears in mdlparse_util.c. It should
      eventually be removed from there and only appear in this file.
*************************************************************************/
static int sort_product_list_compare(struct product *list_item,
                                     struct product *new_item) {

  int cmp = strcmp(list_item->prod->sym->name, new_item->prod->sym->name);
  if (cmp == 0) {
    if (list_item->orientation > new_item->orientation)
      cmp = -1;
    else if (list_item->orientation < new_item->orientation)
      cmp = 1;
    else
      cmp = 0;
  }
  return cmp;
}

/*************************************************************************
 sort_product_list:
    Sorts product_head in alphabetical order, and descending orientation order.
    Current algorithm uses insertion sort.

 In:  product_head: list to sort
 Out: the new list head

  XXX Currently this function also appears in mdlparse_util.c. It should
      eventually be removed from there and only appear in this file.
*************************************************************************/
static struct product *sort_product_list(struct product *product_head) {
  struct product *next; /* Saved next item (next field in product is
                           overwritten) */
  struct product *iter;          /* List iterator */
  struct product *result = NULL; /* Sorted list */
  int cmp;

  /* Use insertion sort to sort the list of products */
  for (struct product *current = product_head; current != NULL;
       current = next) {
    next = current->next;

    /* First item added always goes at the head */
    if (result == NULL) {
      current->next = result;
      result = current;
      continue;
    }

    /* Check if the item belongs at the head */
    cmp = sort_product_list_compare(result, current);
    if (cmp >= 0) {
      current->next = result;
      result = current;
      continue;
    }

    /* Otherwise, if it goes after the current entry, scan forward to find the
       insert point */
    else {
      /* locate the node before the point of insertion */
      iter = result;
      while (iter->next != NULL && sort_product_list_compare(iter, current) < 0)
        iter = iter->next;

      current->next = iter->next;
      iter->next = current;
    }
  }

  return result;
}

/*************************************************************************
 create_prod_signature:
    Returns a string containing all products in the product_head list,
    separated by '+', and sorted in alphabetical order by name and descending
    orientation order.

 In:  product_head: list of products
 Out: product signature as a string.  *product_head list is sorted in
      alphabetical order by name, and descending order by orientation.  Returns
      NULL on failure.

  XXX Currently this function also appears in mdlparse_util.c. It should
      eventually be removed from there and only appear in this file.
*************************************************************************/
char *create_prod_signature(struct product **product_head) {
  /* points to the head of the sorted alphabetically list of products */
  char *prod_signature = NULL;

  *product_head = sort_product_list(*product_head);

  /* create prod_signature string */
  struct product *current = *product_head;
  if (current == NULL) {
    return NULL;
  }
  prod_signature = CHECKED_STRDUP(current->prod->sym->name, "product name");

  /* Concatenate to create product signature */
  char *temp_str = NULL;
  while (current->next != NULL) {
    temp_str = prod_signature;
    prod_signature = CHECKED_SPRINTF("%s+%s", prod_signature,
                                     current->next->prod->sym->name);

    if (prod_signature == NULL) {
      if (temp_str != NULL)
        free(temp_str);
      return NULL;
    }
    if (temp_str != NULL)
      free(temp_str);

    current = current->next;
  }

  return prod_signature;
}

/*************************************************************************
 * init_reactions and related machinery
 *************************************************************************/

/*************************************************************************
 check_duplicate_special_reactions:
   Check for duplicate special reaction pathways (e.g. TRANSPARENT = molecule).

 In: path: Parse-time structure for reaction pathways
 Out: Nothing.
 Note: I'm not sure if this code is ever actually called.
*************************************************************************/
void check_duplicate_special_reactions(struct pathway *path) {
  /* if it is a special reaction - check for the duplicates pathways */
  if (path->next != NULL) {
    if ((path->flags & PATHW_TRANSP) && (path->next->flags & PATHW_TRANSP)) {
      if ((path->orientation2 == path->next->orientation2) ||
          (path->orientation2 == 0) || (path->next->orientation2 == 0)) {
        mcell_error("Exact duplicates of special reaction TRANSPARENT = %s are "
                    "not allowed.  Please verify the contents of "
                    "DEFINE_SURFACE_CLASS statement.",
                    path->reactant2->sym->name);
      }
    }

    if ((path->flags & PATHW_REFLEC) && (path->next->flags & PATHW_REFLEC)) {
      if ((path->orientation2 == path->next->orientation2) ||
          (path->orientation2 == 0) || (path->next->orientation2 == 0)) {
        mcell_error("Exact duplicates of special reaction REFLECTIVE = %s are "
                    "not allowed.  Please verify the contents of "
                    "DEFINE_SURFACE_CLASS statement.",
                    path->reactant2->sym->name);
      }
    }
    if ((path->flags & PATHW_ABSORP) && (path->next->flags & PATHW_ABSORP)) {
      if ((path->orientation2 == path->next->orientation2) ||
          (path->orientation2 == 0) || (path->next->orientation2 == 0)) {
        mcell_error("Exact duplicates of special reaction ABSORPTIVE = %s are "
                    "not allowed.  Please verify the contents of "
                    "DEFINE_SURFACE_CLASS statement.",
                    path->reactant2->sym->name);
      }
    }
  }
}

/*************************************************************************
 set_product_geometries:

  Walk through the list, setting the geometries of each of the products. We do
  this by looking for an earlier geometric match and pointing there or we just
  point to 0 if there is no match.

 In: path: Parse-time structure for reaction pathways
     rx: Pathways leading away from a given intermediate
     prod: Parse-time structure for products of reaction pathways
 Out: max_num_surf_products: Maximum number of surface products
*************************************************************************/
int set_product_geometries(struct pathway *path, struct rxn *rx,
                           struct product *prod) {
  int recycled1, recycled2, recycled3;
  int k, kk, k2;
  short geom;
  struct product *prod2;
  int max_num_surf_products; /* maximum number of surface products */
  int num_surf_products_per_pathway;

  max_num_surf_products = 0;
  for (int n_pathway = 0; path != NULL; n_pathway++, path = path->next) {
    recycled1 = 0;
    recycled2 = 0;
    recycled3 = 0;
    k = rx->product_idx[n_pathway] + rx->n_reactants;
    num_surf_products_per_pathway = 0;
    for (prod = path->product_head; prod != NULL; prod = prod->next) {
      if (recycled1 == 0 && prod->prod == path->reactant1) {
        recycled1 = 1;
        kk = rx->product_idx[n_pathway] + 0;
      } else if (recycled2 == 0 && prod->prod == path->reactant2) {
        recycled2 = 1;
        kk = rx->product_idx[n_pathway] + 1;
      } else if (recycled3 == 0 && prod->prod == path->reactant3) {
        recycled3 = 1;
        kk = rx->product_idx[n_pathway] + 2;
      } else {
        kk = k;
        k++;
      }

      if (prod->prod->flags & ON_GRID)
        num_surf_products_per_pathway++;

      rx->players[kk] = prod->prod;

      if ((prod->orientation + path->orientation1) *
                  (prod->orientation - path->orientation1) ==
              0 &&
          prod->orientation * path->orientation1 != 0) {
        if (prod->orientation == path->orientation1)
          rx->geometries[kk] = 1;
        else
          rx->geometries[kk] = -1;
      } else if (rx->n_reactants > 1 &&
                 (prod->orientation + path->orientation2) *
                         (prod->orientation - path->orientation2) ==
                     0 &&
                 prod->orientation * path->orientation2 != 0) {
        if (prod->orientation == path->orientation2)
          rx->geometries[kk] = 2;
        else
          rx->geometries[kk] = -2;
      } else if (rx->n_reactants > 2 &&
                 (prod->orientation + path->orientation3) *
                         (prod->orientation - path->orientation3) ==
                     0 &&
                 prod->orientation * path->orientation3 != 0) {
        if (prod->orientation == path->orientation3)
          rx->geometries[kk] = 3;
        else
          rx->geometries[kk] = -3;
      } else {
        k2 = 2 * rx->n_reactants + 1; /* Geometry index of first non-reactant
                                         product, counting from 1. */
        geom = 0;
        for (prod2 = path->product_head;
             prod2 != prod && prod2 != NULL && geom == 0; prod2 = prod2->next) {
          if ((prod2->orientation + prod->orientation) *
                      (prod2->orientation - prod->orientation) ==
                  0 &&
              prod->orientation * prod2->orientation != 0) {
            if (prod2->orientation == prod->orientation)
              geom = 1;
            else
              geom = -1;
          } else
            geom = 0;

          if (recycled1 == 1) {
            if (prod2->prod == path->reactant1) {
              recycled1 = 2;
              geom *= rx->n_reactants + 1;
            }
          } else if (recycled2 == 1) {
            if (prod2->prod == path->reactant2) {
              recycled2 = 2;
              geom *= rx->n_reactants + 2;
            }
          } else if (recycled3 == 1) {
            if (prod2->prod == path->reactant3) {
              recycled3 = 2;
              geom *= rx->n_reactants + 3;
            }
          } else {
            geom *= k2;
            k2++;
          }
        }
        rx->geometries[kk] = geom;
      }
      if (num_surf_products_per_pathway > max_num_surf_products)
        max_num_surf_products = num_surf_products_per_pathway;
    }

    k = rx->product_idx[n_pathway];
    if (recycled1 == 0)
      rx->players[k] = NULL;
    if (recycled2 == 0 && rx->n_reactants > 1)
      rx->players[k + 1] = NULL;
    if (recycled3 == 0 && rx->n_reactants > 2)
      rx->players[k + 2] = NULL;
  } /* end for (n_pathway = 0, ...) */
  return max_num_surf_products;
}

/*************************************************************************
 alphabetize_pathway:
    The reaction pathway (path) is alphabetized.

 In: path: Parse-time structure for reaction pathways
     reaction: Reaction pathways leading away from a given intermediate
 Out: Nothing.
*************************************************************************/
void alphabetize_pathway(struct pathway *path, struct rxn *reaction) {
  short geom, geom2;
  struct species *temp_sp, *temp_sp2;

  /* Alphabetize if we have two molecules */
  if ((path->reactant2->flags & IS_SURFACE) == 0) {
    if (strcmp(path->reactant1->sym->name, path->reactant2->sym->name) > 0) {
      temp_sp = path->reactant1;
      path->reactant1 = path->reactant2;
      path->reactant2 = temp_sp;
      geom = path->orientation1;
      path->orientation1 = path->orientation2;
      path->orientation2 = geom;
    } else if (strcmp(path->reactant1->sym->name, path->reactant2->sym->name) ==
               0) {
      if (path->orientation1 < path->orientation2) {
        geom = path->orientation1;
        path->orientation1 = path->orientation2;
        path->orientation2 = geom;
      }
    }
  }

  /* Alphabetize if we have three molecules */
  if (reaction->n_reactants == 3) {
    if ((path->reactant3->flags & IS_SURFACE) == 0) {
      if (strcmp(path->reactant1->sym->name, path->reactant3->sym->name) > 0) {
        /* Put reactant3 at the beginning */
        temp_sp = path->reactant1;
        geom = path->orientation1;
        path->reactant1 = path->reactant3;
        path->orientation1 = path->orientation3;

        /* Put former reactant1 in place of reactant2 */
        temp_sp2 = path->reactant2;
        geom2 = path->orientation2;
        path->reactant2 = temp_sp;
        path->orientation2 = geom;

        /* Put former reactant2 in place of reactant3 */
        path->reactant3 = temp_sp2;
        path->orientation3 = geom2;

      } else if (strcmp(path->reactant2->sym->name,
                        path->reactant3->sym->name) > 0) {

        /* Put reactant3 after reactant1 */
        temp_sp = path->reactant2;
        path->reactant2 = path->reactant3;
        path->reactant3 = temp_sp;
        geom = path->orientation2;
        path->orientation2 = path->orientation3;
        path->orientation3 = geom;
      }
    } /*end */
  }
}

/*************************************************************************
 warn_about_high_rates:
    If HIGH_REACTION_PROBABILITY is set to WARNING or ERROR, and the reaction
    probability is high, give the user a warning or error respectively.

 In: notify:
     warn_file: The log/error file. Can be stdout/stderr
     rate_warn: If 1, warn the user about high reaction rates (or give error)
     print_once: If the warning has been printed once, don't repeat it
 Out: print_once. Also print out reaction probabilities (with warning/error)
*************************************************************************/
static int warn_about_high_rates(struct notifications *notify, FILE *warn_file,
                                 int rate_warn, int print_once) {
  if (rate_warn) {
    if (notify->high_reaction_prob == WARN_ERROR) {
      warn_file = mcell_get_error_file();
      if (!print_once) {
        fprintf(warn_file, "\n");
        fprintf(
            warn_file,
            "Reaction probabilities generated for the following reactions:\n");
        print_once = 1;
      }
      fprintf(warn_file, "\tError: High ");
    } else {
      if (!print_once) {
        fprintf(warn_file, "\n");
        fprintf(
            warn_file,
            "Reaction probabilities generated for the following reactions:\n");
        print_once = 1;
      }
      if (notify->high_reaction_prob == WARN_WARN)
        fprintf(warn_file, "\tWarning: High ");
      else
        fprintf(warn_file, "\t");
    }
  } else {
    if (!print_once) {
      fprintf(warn_file, "\n");
      fprintf(
          warn_file,
          "Reaction probabilities generated for the following reactions:\n");
      print_once = 1;
    }
    fprintf(warn_file, "\t");
  }
  return print_once;
}

/*************************************************************************
 add_surface_reaction_flags:

 In:  mol_sym_table:
      all_mols:
      all_surface_mols:
      all_volume_mols:
 Out: Nothing
*************************************************************************/
void add_surface_reaction_flags(struct sym_table_head *mol_sym_table,
                                struct species *all_mols,
                                struct species *all_surface_mols,
                                struct species *all_volume_mols) {
  struct species *temp_sp;

  /* Add flags for surface reactions with ALL_MOLECULES */
  if (all_mols->flags & (CAN_VOLWALL | CAN_SURFWALL)) {
    for (int n_mol_bin = 0; n_mol_bin < mol_sym_table->n_bins;
         n_mol_bin++) {
      for (struct sym_entry *symp = mol_sym_table->entries[n_mol_bin];
           symp != NULL; symp = symp->next) {
        temp_sp = (struct species *)symp->value;
        if (temp_sp == all_mols)
          continue;
        if (temp_sp == all_volume_mols)
          continue;
        if (temp_sp == all_surface_mols)
          continue;

        if (((temp_sp->flags & NOT_FREE) == 0) &&
            ((temp_sp->flags & CAN_VOLWALL) == 0)) {
          temp_sp->flags |= CAN_VOLWALL;
        } else if ((temp_sp->flags & ON_GRID) &&
                   ((temp_sp->flags & CAN_REGION_BORDER) == 0)) {
          temp_sp->flags |= CAN_REGION_BORDER;
        }
      }
    }
  }

  /* Add flags for surface reactions with ALL_VOLUME_MOLECULES */
  if (all_volume_mols->flags & CAN_VOLWALL) {
    for (int n_mol_bin = 0; n_mol_bin < mol_sym_table->n_bins;
         n_mol_bin++) {
      for (struct sym_entry *symp = mol_sym_table->entries[n_mol_bin];
           symp != NULL; symp = symp->next) {
        temp_sp = (struct species *)symp->value;
        if (temp_sp == all_mols)
          continue;
        if (temp_sp == all_volume_mols)
          continue;
        if (temp_sp == all_surface_mols)
          continue;
        if (((temp_sp->flags & NOT_FREE) == 0) &&
            ((temp_sp->flags & CAN_VOLWALL) == 0)) {
          temp_sp->flags |= CAN_VOLWALL;
        }
      }
    }
  }

  /* Add flags for surface reactions with ALL_SURFACE_MOLECULES */
  if (all_surface_mols->flags & CAN_SURFWALL) {
    for (int n_mol_bin = 0; n_mol_bin < mol_sym_table->n_bins;
         n_mol_bin++) {
      for (struct sym_entry *symp = mol_sym_table->entries[n_mol_bin];
           symp != NULL; symp = symp->next) {
        temp_sp = (struct species *)symp->value;
        if (temp_sp == all_mols)
          continue;
        if (temp_sp == all_volume_mols)
          continue;
        if (temp_sp == all_surface_mols)
          continue;
        if (((temp_sp->flags & ON_GRID) &&
             ((temp_sp->flags & CAN_REGION_BORDER) == 0))) {
          temp_sp->flags |= CAN_REGION_BORDER;
        }
      }
    }
  }
}

/*************************************************************************
 scale_probabilities:

  Scale probabilities, notifying and warning as appropriate.

 In: reaction_prob_limit_flag:
     path: Parse-time structure for reaction pathways
     rx: Pathways leading away from a given intermediate
     pb_factor:
 Out: Return 1 if rates are high and HIGH_REACTION_PROBABILITY is set to ERROR
 Note: This does not work properly right now. Even if rates are high and
       HIGH_REACTION_PROBABILITY is set to ERROR, the error is ignored
*************************************************************************/
int scale_probabilities(byte *reaction_prob_limit_flag,
                        struct notifications *notify,
                        struct pathway *path, struct rxn *rx,
                        double pb_factor) {
  int print_once = 0; /* flag */
  FILE *warn_file;
  int is_gigantic;
  double rate;

  for (int n_pathway = 0; path != NULL; n_pathway++, path = path->next) {
    int rate_notify = 0, rate_warn = 0;
    if (!distinguishable(rx->cum_probs[n_pathway], GIGANTIC, EPS_C))
      is_gigantic = 1;
    else
      is_gigantic = 0;

    /* automatic surface reactions will be printed out from 'init_sim()'. */
    if (is_gigantic)
      continue;

    rate = pb_factor * rx->cum_probs[n_pathway];
    rx->cum_probs[n_pathway] = rate;

    if ((notify->reaction_probabilities == NOTIFY_FULL &&
         ((rate >= notify->reaction_prob_notify) ||
          (notify->reaction_prob_notify == 0.0))))
      rate_notify = 1;
    if ((notify->high_reaction_prob != WARN_COPE &&
         ((rate >= notify->reaction_prob_warn) ||
          ((notify->reaction_prob_warn == 0.0)))))
      rate_warn = 1;

    if ((rate > 1.0) && (!*reaction_prob_limit_flag)) {
      *reaction_prob_limit_flag = 1;
    }

    if (rate_warn || rate_notify) {

      warn_file = mcell_get_log_file();

      print_once =
          warn_about_high_rates(notify, warn_file, rate_warn,
                                print_once);

      fprintf(warn_file, "Probability %.4e set for ", rate);

      if (rx->n_reactants == 1)
        fprintf(warn_file, "%s{%d} -> ", rx->players[0]->sym->name,
                rx->geometries[0]);
      else if (rx->n_reactants == 2) {
        if (rx->players[1]->flags & IS_SURFACE) {
          fprintf(warn_file, "%s{%d} @ %s{%d} -> ", rx->players[0]->sym->name,
                  rx->geometries[0], rx->players[1]->sym->name,
                  rx->geometries[1]);
        } else {
          fprintf(warn_file, "%s{%d} + %s{%d} -> ", rx->players[0]->sym->name,
                  rx->geometries[0], rx->players[1]->sym->name,
                  rx->geometries[1]);
        }
      } else {
        if (rx->players[2]->flags & IS_SURFACE) {
          fprintf(warn_file, "%s{%d} + %s{%d}  @ %s{%d} -> ",
                  rx->players[0]->sym->name, rx->geometries[0],
                  rx->players[1]->sym->name, rx->geometries[1],
                  rx->players[2]->sym->name, rx->geometries[2]);
        } else {
          fprintf(warn_file, "%s{%d} + %s{%d}  + %s{%d} -> ",
                  rx->players[0]->sym->name, rx->geometries[0],
                  rx->players[1]->sym->name, rx->geometries[1],
                  rx->players[2]->sym->name, rx->geometries[2]);
        }
      }
      if (path->product_head == NULL) {
        fprintf(warn_file, "NULL ");
      } else {
        for (struct product *prod = path->product_head; prod != NULL;
             prod = prod->next) {
          fprintf(warn_file, "%s{%d} ", prod->prod->sym->name,
                  prod->orientation);
        }
      }

      fprintf(warn_file, "\n");

      if (rate_warn && notify->high_reaction_prob == WARN_ERROR)
        return 1;
    }
  }
  return 0;
}

/*************************************************************************
 equivalent_geometry_for_two_reactants:

 In: o1a: orientation of the first reactant from first reaction
     o1b: orientation of the second reactant from first reaction
     o2a: orientation of the first reactant from second reaction
     o2b: orientation of the second reactant from second reaction
 Out: Returns 1 if the two pathways (defined by pairs o1a-o1b and o2a-o2b)
      have equivalent geometry, 0 otherwise.
*************************************************************************/
static int equivalent_geometry_for_two_reactants(int o1a, int o1b, int o2a,
                                                 int o2b) {

  /* both reactants for each pathway are in the same
     orientation class and parallel one another */
  if ((o1a == o1b) && (o2a == o2b)) {
    return 1;
    /* both reactants for each pathway are in the same
       orientation class and opposite one another */
  } else if ((o1a == -o1b) && (o2a == -o2b)) {
    return 1;
  }
  /* reactants are not in the same orientation class */
  if (abs(o1a) != abs(o1b)) {
    if ((abs(o2a) != abs(o2b)) || ((o2a == 0) && (o2b == 0))) {
      return 1;
    }
  }
  if (abs(o2a) != abs(o2b)) {
    if ((abs(o1a) != abs(o1b)) || ((o1a == 0) && (o1b == 0))) {
      return 1;
    }
  }

  return 0;
}

/*************************************************************************
 equivalent_geometry:

 In: p1, p2: pathways to compare
     n: The number of reactants for the pathways
 Out: Returns 1 if the two pathways are the same (i.e. have equivalent
      geometry), 0 otherwise.
*************************************************************************/
static int equivalent_geometry(struct pathway *p1, struct pathway *p2, int n) {

  short o11, o12, o13, o21, o22, o23; /* orientations of individual reactants */
  /* flags for 3-reactant reactions signaling whether molecules orientations
   * are parallel one another and molecule and surface orientaions are parallel
   * one another
   */
  int mols_parallel_1 = SHRT_MIN + 1;     /* for first pathway */
  int mols_parallel_2 = SHRT_MIN + 2;     /* for second pathway */
  int mol_surf_parallel_1 = SHRT_MIN + 3; /* for first pathway */
  int mol_surf_parallel_2 = SHRT_MIN + 4; /* for second pathway */

  if (n < 2) {
    /* one reactant case */
    /* RULE: all one_reactant pathway geometries are equivalent */

    return 1;

  } else if (n < 3) {
    /* two reactants case */

    /* RULE - Two pathways have equivalent geometry when:
       1) Both pathways have exactly the same number of reactants;
       2) There exists an identity mapping between reactants from Pathway 1 and
          Pathway 2 such that for each pair of reactants, r1a and r1b from
       Pathway
          1, and r2a, and r2b from Pathway 2:
         - r1a is the same species as r2a (likewise for r1b and r2b);
         - r1a and r1b have the same orientation in the same orientation class
           if and only if r2a and r2b do;
         - r1a and r1b have the opposite orientation in the same orientation
           class if and only if r2a and r2b do;
         - r1a and r1b are not in the same orientation class, either because
           they have different orientation classes or both are in the zero
           orientation class, if and only if r2a and r2b are not in the same
           orientation class or both are in the zero orientation class
     */

    o11 = p1->orientation1;
    o12 = p1->orientation2;
    o21 = p2->orientation1;
    o22 = p2->orientation2;

    return equivalent_geometry_for_two_reactants(o11, o12, o21, o22);

  } else if (n < 4) {
    /* three reactants case */

    o11 = p1->orientation1;
    o12 = p1->orientation2;
    o13 = p1->orientation3;
    o21 = p2->orientation1;
    o22 = p2->orientation2;
    o23 = p2->orientation3;

    /* special case: two identical reactants */
    if ((p1->reactant1 == p1->reactant2) && (p2->reactant1 == p2->reactant2)) {

      /* Case 1: two molecules and surface are in the same orientation class */
      if ((abs(o11) == abs(o12)) && (abs(o11) == abs(o13))) {
        if (o11 == o12)
          mols_parallel_1 = 1;
        else
          mols_parallel_1 = 0;

        if (mols_parallel_1) {
          if ((o11 == -o13) || (o12 == -o13)) {
            mol_surf_parallel_1 = 0;
          } else {
            mol_surf_parallel_1 = 1;
          }
        } else {
          mol_surf_parallel_1 = 0;
        }

        if ((abs(o21) == abs(o22)) && (abs(o21) == abs(o23))) {
          if (o21 == o22)
            mols_parallel_2 = 1;
          else
            mols_parallel_2 = 0;

          if (mols_parallel_2) {
            if ((o21 == -o23) || (o22 == -o23)) {
              mol_surf_parallel_2 = 0;
            } else {
              mol_surf_parallel_2 = 1;
            }
          } else {
            mol_surf_parallel_2 = 0;
          }
        }

        if ((mols_parallel_1 == mols_parallel_2) &&
            (mol_surf_parallel_1 == mol_surf_parallel_2)) {
          return 1;
        }

      } /* end case 1 */

      /* Case 2: one molecule and surface are in the same orientation class */
      else if ((abs(o11) == abs(o13)) || (abs(o12) == abs(o13))) {
        if ((o11 == o13) || (o12 == o13))
          mol_surf_parallel_1 = 1;
        else
          mol_surf_parallel_1 = 0;

        /* check that pathway2 is also in the case2 */

        if ((abs(o21) != abs(o23)) || (abs(o22) != abs(o23))) {
          if ((abs(o21) == abs(o23)) || (abs(o22) == abs(o23))) {
            if ((o21 == o23) || (o22 == o23))
              mol_surf_parallel_2 = 1;
            else
              mol_surf_parallel_2 = 0;
          }
        }
        if (mol_surf_parallel_1 == mol_surf_parallel_2) {
          return 1;
        }

      } /* end case 2 */

      /* Case 3: two molecules but not surface are in the same
                 orientation class */
      else if ((abs(o11) == abs(o12)) && (abs(o11) != abs(o13))) {
        if (o11 == o12)
          mols_parallel_1 = 1;
        else
          mols_parallel_1 = 0;

        if ((abs(o21) == abs(o22)) && (abs(o21) != abs(o23))) {
          if (o21 == o22)
            mols_parallel_2 = 1;
          else
            mols_parallel_2 = 0;
        }
        if (mols_parallel_1 == mols_parallel_2) {
          return 1;
        }

      }
      /* Case 4: all molecules and surface are in different orientation classes
         */
      else if ((abs(o11) != abs(o13)) && (abs(o12) != abs(o13)) &&
               (abs(o11) != abs(o12))) {

        if ((abs(o21) != abs(o23)) && (abs(o22) != abs(o23)) &&
            (abs(o21) != abs(o22))) {
          return 1;
        }
      } /* end all cases */

    } else { /* no identical reactants */

      if ((equivalent_geometry_for_two_reactants(o11, o12, o21, o22)) &&
          (equivalent_geometry_for_two_reactants(o12, o13, o22, o23)) &&
          (equivalent_geometry_for_two_reactants(o11, o13, o21, o23))) {
        return 1;
      }
    }

  } // end if (n < 4)

  return 0;
}

/*************************************************************************
 create_sibling_reaction:
    Create a sibling reaction to the given reaction -- a reaction into which
    some of the pathways may be split by split_reaction.

 In:  rx:   reaction for whom to create sibling
 Out: sibling reaction, or NULL on error
*************************************************************************/
static struct rxn *create_sibling_reaction(struct rxn *rx) {

  struct rxn *reaction = CHECKED_MALLOC_STRUCT(struct rxn, "reaction");
  if (reaction == NULL)
    return NULL;
  reaction->next = NULL;
  reaction->sym = rx->sym;
  reaction->n_reactants = rx->n_reactants;
  reaction->n_pathways = 0;
  reaction->cum_probs = NULL;
  reaction->product_idx = NULL;
  reaction->max_fixed_p = 0.0;
  reaction->min_noreaction_p = 0.0;
  reaction->pb_factor = 0.0;
  reaction->players = NULL;
  reaction->geometries = NULL;
  reaction->n_occurred = 0;
  reaction->n_skipped = 0.0;
  reaction->prob_t = NULL;
  reaction->pathway_head = NULL;
  reaction->info = NULL;
  return reaction;
}

/*************************************************************************
 split_reaction:
 In:  rx: reaction to split
 Out: Returns head of the linked list of reactions where each reaction
      contains only geometrically equivalent pathways
*************************************************************************/
struct rxn *split_reaction(struct rxn *rx) {
  struct rxn *curr_rxn_ptr = NULL, *head = NULL, *end = NULL;
  struct rxn *reaction;
  struct pathway *to_place, *temp;

  /* keep reference to the head of the future linked_list */
  head = end = rx;
  to_place = head->pathway_head->next;
  head->pathway_head->next = NULL;
  head->n_pathways = 1;
  while (to_place != NULL) {
    if (to_place->flags &
        (PATHW_TRANSP | PATHW_REFLEC | PATHW_ABSORP | PATHW_CLAMP_CONC)) {
      reaction = create_sibling_reaction(rx);
      if (reaction == NULL)
        return NULL;

      reaction->pathway_head = to_place;
      to_place = to_place->next;
      reaction->pathway_head->next = NULL;
      ++reaction->n_pathways;

      end->next = reaction;
      end = reaction;
    } else {
      for (curr_rxn_ptr = head; curr_rxn_ptr != NULL;
           curr_rxn_ptr = curr_rxn_ptr->next) {
        if (curr_rxn_ptr->pathway_head->flags &
            (PATHW_TRANSP | PATHW_REFLEC | PATHW_ABSORP))
          continue;
        if (equivalent_geometry(to_place, curr_rxn_ptr->pathway_head,
                                curr_rxn_ptr->n_reactants))
          break;
      }

      if (!curr_rxn_ptr) {
        reaction = create_sibling_reaction(rx);
        if (reaction == NULL)
          return NULL;

        end->next = reaction;
        end = reaction;

        curr_rxn_ptr = end;
      }

      temp = to_place;
      to_place = to_place->next;

      temp->next = curr_rxn_ptr->pathway_head;
      curr_rxn_ptr->pathway_head = temp;
      ++curr_rxn_ptr->n_pathways;
    }
  }

  return head;
}

/*************************************************************************
 check_reaction_for_duplicate_pathways:
 In:  head: head of linked list of pathways
 Out: Sorts linked list of pathways in alphabetical order according to the
      "prod_signature" field.  Checks for the duplicate pathways.  Prints error
      message and exits simulation if duplicates found.
 Note: This function is called after 'split_reaction()' function so all
       pathways have equivalent geometry from the reactant side.  Here we check
       whether relative orientation of all players (both reactants and
       products) is the same for the two seemingly identical pathways.
 RULE: Two reactions pathways are duplicates if and only if
        (a) they both have the same number and species of reactants;
        (b) they both have the same number and species of products;
        (c) there exists a bijective mapping between the reactants and products
            of the two pathways such that reactants map to reactants, products
            map to products, and the two pathways have equivalent geometry
            under mapping.
            Two pathways R1 and R2 have an equivalent geometry under a mapping
            M if and only if for every pair of players "i" and "j" in R1, the
            corresponding players M(i) and M(j) in R2 have the same orientation
            relation as do "i" and "j" in R1.
            Two players "i" and "j" in a reaction pathway have the following
            orientation:
              parallel - if both "i" and "j" are in the same nonzero orientation
              class with the same sign;
              antiparallel (opposite) - if they are both in the same nonzero
              orientation class but have opposite sign;
              independent - if they are in different orientation classes or both
              in the zero orientation class.

 PostNote: In this function we check only the validity of Rule (c) since
           conditions of Rule (a) and (b) are already satisfied when the
           function is called.
*************************************************************************/
void check_reaction_for_duplicate_pathways(struct pathway **head) {

  struct pathway *result = NULL;      /* build the sorted list here */
  struct pathway *null_result = NULL; /* put pathways with NULL
                                         prod_signature field here */
  struct pathway *current, *next, **pprev;
  struct product *iter1, *iter2;
  int pathways_equivalent; /* flag */
  int i, j;
  int num_reactants; /* number of reactants in the pathway */
  int num_products;  /* number of products in the pathway */
  int num_players;   /* total number of reactants and products in the pathway */
  int *orient_players_1,
      *orient_players_2; /* array of orientations of players */
  int o1a, o1b, o2a, o2b;

  /* extract  pathways with "prod_signature" field equal to NULL into
   * "null_result" list */
  current = *head;
  pprev = head;
  while (current != NULL) {
    if (current->prod_signature == NULL) {
      *pprev = current->next;
      current->next = null_result;
      null_result = current;
      current = *pprev;
    } else {
      pprev = &current->next;
      current = current->next;
    }
  }

  /* check for duplicate pathways in null_result */
  current = null_result;
  if ((current != NULL) && (current->next != NULL)) {
    /* From the previously called function "split_reaction()" we know that
     * reactant-reactant pairs in two pathways are equivalent. Because there
     * are no products the pathways are duplicates.
       RULE: There may be no more than one pathway with zero (--->NULL)
             products in the reaction->pathway_head
             after calling the function "split_reaction()"
    */
    if (current->reactant2 == NULL)
      mcell_error("Exact duplicates of reaction %s  ----> NULL are not "
                  "allowed.  Please verify that orientations of reactants are "
                  "not equivalent.",
                  current->reactant1->sym->name);
    else if (current->reactant3 == NULL)
      mcell_error("Exact duplicates of reaction %s + %s  ----> NULL are not "
                  "allowed.  Please verify that orientations of reactants are "
                  "not equivalent.",
                  current->reactant1->sym->name, current->reactant2->sym->name);
    else
      mcell_error("Exact duplicates of reaction %s + %s + %s  ----> NULL are "
                  "not allowed.  Please verify that orientations of reactants "
                  "are not equivalent.",
                  current->reactant1->sym->name, current->reactant2->sym->name,
                  current->reactant3->sym->name);
  }

  /* now sort the remaining pathway list by "prod_signature" field and check
   * for the duplicates */
  current = *head;

  while (current != NULL) {
    next = current->next;

    /* insert in sorted order into the "result" */
    if (result == NULL ||
        (strcmp(result->prod_signature, current->prod_signature) >= 0)) {
      current->next = result;
      result = current;
    } else {
      struct pathway *iter = result;
      while (iter->next != NULL && (strcmp(iter->next->prod_signature,
                                           current->prod_signature) < 0)) {
        iter = iter->next;
      }
      current->next = iter->next;
      iter->next = current;
    }

    /* move along the original list */
    current = next;
  }

  /* Now check for the duplicate pathways */
  /* Since the list is sorted we can proceed down the list and compare the
   * adjacent nodes */

  current = result;

  if (current != NULL) {
    while (current->next != NULL) {
      if (strcmp(current->prod_signature, current->next->prod_signature) == 0) {

        pathways_equivalent = 1;
        /* find total number of players in the pathways */
        num_reactants = 0;
        num_products = 0;
        if (current->reactant1 != NULL)
          num_reactants++;
        if (current->reactant2 != NULL)
          num_reactants++;
        if (current->reactant3 != NULL)
          num_reactants++;

        iter1 = current->product_head;
        while (iter1 != NULL) {
          num_products++;
          iter1 = iter1->next;
        }

        num_players = num_reactants + num_products;

        /* create arrays of players orientations */
        orient_players_1 = CHECKED_MALLOC_ARRAY(int, num_players,
                                                "reaction player orientations");
        if (orient_players_1 == NULL)
          mcell_die();
        orient_players_2 = CHECKED_MALLOC_ARRAY(int, num_players,
                                                "reaction player orientations");
        if (orient_players_2 == NULL)
          mcell_die();

        if (current->reactant1 != NULL)
          orient_players_1[0] = current->orientation1;
        if (current->reactant2 != NULL)
          orient_players_1[1] = current->orientation2;
        if (current->reactant3 != NULL)
          orient_players_1[2] = current->orientation3;
        if (current->next->reactant1 != NULL)
          orient_players_2[0] = current->next->orientation1;
        if (current->next->reactant2 != NULL)
          orient_players_2[1] = current->next->orientation2;
        if (current->next->reactant3 != NULL)
          orient_players_2[2] = current->next->orientation3;

        iter1 = current->product_head;
        iter2 = current->next->product_head;

        for (i = num_reactants; i < num_players; i++) {
          orient_players_1[i] = iter1->orientation;
          orient_players_2[i] = iter2->orientation;
          iter1 = iter1->next;
          iter2 = iter2->next;
        }

        /* below we will compare only reactant-product and product-product
         * combinations because reactant-reactant combinations were compared
         * previously in the function "equivalent_geometry()" */

        /* Initial assumption - pathways are equivalent. We check whether this
         * assumption is valid by comparing pairs as described above */

        i = 0;
        while ((i < num_players) && (pathways_equivalent)) {
          if (i < num_reactants) {
            j = num_reactants;
          } else {
            j = i + 1;
          }
          for (; j < num_players; j++) {
            o1a = orient_players_1[i];
            o1b = orient_players_1[j];
            o2a = orient_players_2[i];
            o2b = orient_players_2[j];
            if (!equivalent_geometry_for_two_reactants(o1a, o1b, o2a, o2b)) {
              pathways_equivalent = 0;
              break;
            }
          }
          i++;
        }

        if (pathways_equivalent) {
          if (current->reactant1 != NULL) {
            if (current->reactant2 == NULL)
              mcell_error("Exact duplicates of reaction %s  ----> %s are not "
                          "allowed.  Please verify that orientations of "
                          "reactants are not equivalent.",
                          current->reactant1->sym->name, current->prod_signature);
            else if (current->reactant3 == NULL)
              mcell_error("Exact duplicates of reaction %s + %s  ----> %s are "
                          "not allowed.  Please verify that orientations of "
                          "reactants are not equivalent.",
                          current->reactant1->sym->name,
                          current->reactant2->sym->name, current->prod_signature);
            else
              mcell_error("Exact duplicates of reaction %s + %s + %s  ----> %s "
                          "are not allowed.  Please verify that orientations of "
                          "reactants are not equivalent.",
                          current->reactant1->sym->name,
                          current->reactant2->sym->name,
                          current->reactant3->sym->name, current->prod_signature);
          }
        }
        free(orient_players_1);
        free(orient_players_2);
      }

      current = current->next;
    }
  }

  if (null_result == NULL) {
    *head = result;
  } else if (result == NULL) {
    *head = null_result;
  } else {
    current = result;
    while (current->next != NULL) {
      current = current->next;
    }
    current->next = null_result;
    null_result->next = NULL;

    *head = result;
  }
}

/*************************************************************************
 set_reaction_player_flags:
    Set the reaction player flags for all participating species in this
    reaction.

 In:  rx: the reaction
 Out: Nothing.  Species flags may be updated.
*************************************************************************/
void set_reaction_player_flags(struct rxn *rx) {
  switch (rx->n_reactants) {
  case 1:
    /* do nothing */
    return;

  case 2:
    if (strcmp(rx->players[0]->sym->name, "ALL_MOLECULES") == 0) {
      rx->players[0]->flags |= (CAN_VOLWALL | CAN_SURFWALL);
    } else if (strcmp(rx->players[0]->sym->name, "ALL_VOLUME_MOLECULES") == 0) {
      rx->players[0]->flags |= CAN_VOLWALL;
    } else if (strcmp(rx->players[0]->sym->name, "ALL_SURFACE_MOLECULES") ==
               0) {
      rx->players[0]->flags |= CAN_SURFWALL;
    } else if ((rx->players[0]->flags & NOT_FREE) == 0) {
      /* two volume molecules */
      if ((rx->players[1]->flags & NOT_FREE) == 0) {
        rx->players[0]->flags |= CAN_VOLVOL;
        rx->players[1]->flags |= CAN_VOLVOL;
      }
      /* one volume molecules and one wall */
      else if ((rx->players[1]->flags & IS_SURFACE) != 0) {
        rx->players[0]->flags |= CAN_VOLWALL;
      }
      /* one volume molecule and one surface molecule */
      else if ((rx->players[1]->flags & ON_GRID) != 0) {
        rx->players[0]->flags |= CAN_VOLSURF;
      }
    } else if ((rx->players[0]->flags & IS_SURFACE) != 0) {
      /* one volume molecule and one wall */
      if ((rx->players[1]->flags & NOT_FREE) == 0) {
        rx->players[1]->flags |= CAN_VOLWALL;
      }
      /* one surface molecule and one wall */
      else if ((rx->players[1]->flags & ON_GRID) != 0) {
        rx->players[1]->flags |= CAN_SURFWALL;
      }
    } else if ((rx->players[0]->flags & ON_GRID) != 0) {
      /* one volume molecule and one surface molecule */
      if ((rx->players[1]->flags & NOT_FREE) == 0) {
        rx->players[1]->flags |= CAN_VOLSURF;
      }
      /* two surface molecules */
      else if ((rx->players[1]->flags & ON_GRID) != 0) {
        rx->players[0]->flags |= CAN_SURFSURF;
        rx->players[1]->flags |= CAN_SURFSURF;
      }
      /* one surface molecule and one wall */
      else if ((rx->players[1]->flags & IS_SURFACE) != 0) {
        rx->players[0]->flags |= CAN_SURFWALL;
      }
    }
    break;

  case 3:
    if ((rx->players[2]->flags & IS_SURFACE) != 0) {
      /* two molecules and surface  */
      if ((rx->players[0]->flags & NOT_FREE) == 0) {
        /* one volume molecule, one surface molecule, one surface */
        if ((rx->players[1]->flags & ON_GRID) != 0) {
          rx->players[0]->flags |= CAN_VOLSURF;
        }
      } else if ((rx->players[0]->flags & ON_GRID) != 0) {
        /* one volume molecule, one surface molecule, one surface */
        if ((rx->players[1]->flags & NOT_FREE) == 0) {
          rx->players[1]->flags |= CAN_VOLSURF;
        }
        /* two surface molecules, one surface */
        else if ((rx->players[1]->flags & ON_GRID) != 0) {
          rx->players[0]->flags |= CAN_SURFSURF;
          rx->players[1]->flags |= CAN_SURFSURF;
        }
      }
    } else {
      if ((rx->players[0]->flags & NOT_FREE) == 0) {
        if ((rx->players[1]->flags & NOT_FREE) == 0) {
          /* three volume molecules */
          if ((rx->players[2]->flags & NOT_FREE) == 0) {
            rx->players[0]->flags |= CAN_VOLVOLVOL;
            rx->players[1]->flags |= CAN_VOLVOLVOL;
            rx->players[2]->flags |= CAN_VOLVOLVOL;
          }
          /* two volume molecules and one surface molecule */
          else if ((rx->players[2]->flags & ON_GRID) != 0) {
            rx->players[0]->flags |= CAN_VOLVOLSURF;
            rx->players[1]->flags |= CAN_VOLVOLSURF;
          }
        } else if ((rx->players[1]->flags & ON_GRID) != 0) {
          /* two volume molecules and one surface molecule */
          if ((rx->players[2]->flags & NOT_FREE) == 0) {
            rx->players[0]->flags |= CAN_VOLVOLSURF;
            rx->players[2]->flags |= CAN_VOLVOLSURF;
          }
          /* one volume molecules and two surface molecules */
          else if ((rx->players[2]->flags & ON_GRID) != 0) {
            rx->players[0]->flags |= CAN_VOLSURFSURF;
          }
        }
      } else if ((rx->players[0]->flags & ON_GRID) != 0) {
        if ((rx->players[1]->flags & NOT_FREE) == 0) {
          /* two volume molecules and one surface molecule */
          if ((rx->players[2]->flags & NOT_FREE) == 0) {
            rx->players[1]->flags |= CAN_VOLVOLSURF;
            rx->players[2]->flags |= CAN_VOLVOLSURF;
          }
          /* one volume molecule and two surface molecules */
          else if ((rx->players[2]->flags & ON_GRID) != 0) {
            rx->players[1]->flags |= CAN_VOLSURFSURF;
          }
        } else if ((rx->players[1]->flags & ON_GRID) != 0) {
          /* one volume molecule and two surface molecules */
          if ((rx->players[2]->flags & NOT_FREE) == 0) {
            rx->players[2]->flags |= CAN_VOLSURFSURF;
          }
          /* three surface molecules */
          else if ((rx->players[2]->flags & ON_GRID) != 0) {
            rx->players[0]->flags |= CAN_SURFSURFSURF;
            rx->players[1]->flags |= CAN_SURFSURFSURF;
            rx->players[2]->flags |= CAN_SURFSURFSURF;
          }
        }
      }
    }
    break;

  default:
    // assert(0);
    break;
  }
}

/*************************************************************************
 build_reaction_hash_table:
    Scan the symbol table, copying all reactions found into the reaction hash.

 In:  reaction_hash:
      n_reactions:
      rxn_sym_table:
      rx_hashsize:
      num_rx: num reactions expected
 Out: 0 on success, 1 if we fail to allocate the table
*************************************************************************/
int build_reaction_hash_table(
    struct rxn ***reaction_hash, int *n_reactions,
    struct sym_table_head *rxn_sym_table, int *rx_hashsize, int num_rx) {
  struct rxn **rx_tbl = NULL;
  int rx_hash;
  for (rx_hash = 2; rx_hash <= num_rx && rx_hash != 0; rx_hash <<= 1)
    ;
  rx_hash <<= 1;

  if (rx_hash == 0)
    rx_hash = MAX_RX_HASH;
  if (rx_hash > MAX_RX_HASH)
    rx_hash = MAX_RX_HASH;
#ifdef REPORT_RXN_HASH_STATS
  mcell_log("Num rxns: %d", num_rx);
  mcell_log("Size of hash: %d", rx_hash);
#endif

  /* Create the reaction hash table */
  *rx_hashsize = rx_hash;
  rx_hash -= 1;
  rx_tbl = CHECKED_MALLOC_ARRAY(struct rxn *, *rx_hashsize,
                                "reaction hash table");
  if (rx_tbl == NULL)
    return 1;
  *reaction_hash = rx_tbl;
  for (int i = 0; i <= rx_hash; i++)
    rx_tbl[i] = NULL;

#ifdef REPORT_RXN_HASH_STATS
  int numcoll = 0;
#endif
  for (int i = 0; i < rxn_sym_table->n_bins; i++) {
    for (struct sym_entry *sym = rxn_sym_table->entries[i]; sym != NULL;
         sym = sym->next) {

      struct rxn *rx = (struct rxn *)sym->value;
      int table_slot;
      if (rx->n_reactants == 1) {
        table_slot = rx->players[0]->hashval & rx_hash;
      } else {
        table_slot =
            (rx->players[0]->hashval + rx->players[1]->hashval) & rx_hash;
      }

#ifdef REPORT_RXN_HASH_STATS
      if (rx_tbl[table_slot] != NULL) {
        mcell_log("Collision: %s and %s", rx_tbl[table_slot]->sym->name,
                  sym->name);
        ++numcoll;
      }
#endif
      *n_reactions = *n_reactions + 1;
      while (rx->next != NULL)
        rx = rx->next;
      rx->next = rx_tbl[table_slot];
      rx_tbl[table_slot] = (struct rxn *)sym->value;
    }
  }
#ifdef REPORT_RXN_HASH_STATS
  mcell_log("Num collisions: %d", numcoll);
#endif

  return 0;
}

/*****************************************************************************
 *
 * mcell_create_reaction_rates list creates a struct reaction_rates used
 * for creating reactions from a forward and backward reaction rate.
 * The backward rate is only needed for catalytic arrow and should be
 * RATE_UNUSED otherwise.
 *
 *****************************************************************************/
struct reaction_rates mcell_create_reaction_rates(int forwardRateType,
                                                  double forwardRateConstant,
                                                  int backwardRateType,
                                                  double backwardRateConstant) {
  struct reaction_rate forwardRate;
  forwardRate.rate_type = forwardRateType;
  forwardRate.v.rate_constant = forwardRateConstant;

  struct reaction_rate backwardRate;
  backwardRate.rate_type = backwardRateType;
  backwardRate.v.rate_constant = backwardRateConstant;

  struct reaction_rates rates = { forwardRate, backwardRate };

  return rates;
}

/*************************************************************************
 load_rate_file:
    Read in a time-varying reaction rate constant file.

 In:  time_unit:
      tv_rxn_mem:
      rx:    Reaction structure that we'll load the rates into.
      fname: Filename to read the rates from.
      path:  Index of the pathway that these rates apply to.
      neg_reaction: warning or error policy for negative reactions.
 Out: Returns 1 on error, 0 on success.
      Rate constants are added to the prob_t linked list. If there is a rate
      constant given for time <= 0, then this rate constant is stuck into
      cum_probs and the (time <= 0) entries are not added to the list.  If no
      initial rate constnat is given in the file, it is assumed to be zero.
 Note: The file format is assumed to be two columns of numbers; the first
       column is time (in seconds) and the other is rate constant (in
       appropriate units) that starts at that time.  Lines that are not numbers
       are ignored.
*************************************************************************/
int load_rate_file(double time_unit, struct mem_helper *tv_rxn_mem,
                   struct rxn *rx, char *fname, int path,
                   enum warn_level_t neg_reaction) {

  const char *RATE_SEPARATORS = "\f\n\r\t\v ,;";
  const char *FIRST_DIGIT = "+-0123456789";
  int i;
  FILE *f = fopen(fname, "r");

  if (!f)
    return 1;
  else {
    struct t_func *tp, *tp2;
    double t, rate_constant;
    char buf[2048];
    char *cp;
    int linecount = 0;
#ifdef DEBUG
    int valid_linecount = 0;
#endif

    tp2 = NULL;
    while (fgets(buf, 2048, f)) {
      linecount++;
      for (i = 0; i < 2048; i++) {
        if (!strchr(RATE_SEPARATORS, buf[i]))
          break;
      }

      if (i < 2048 && strchr(FIRST_DIGIT, buf[i])) {
        t = strtod((buf + i), &cp);
        if (cp == (buf + i))
          continue; /* Conversion error. */

        for (i = cp - buf; i < 2048; i++) {
          if (!strchr(RATE_SEPARATORS, buf[i]))
            break;
        }
        // This is kind of a silly corner case, but let's check for it to keep
        // coverity happy.
        if (i == 2048) {
          mcell_error(
            "a time in the rate constant file consists of too many characters "
            "(it uses 2048 or more characters).");
          return(1);
        }
        rate_constant = strtod((buf + i), &cp);
        if (cp == (buf + i))
          continue; /* Conversion error */

        /* at this point we need to handle negative reaction rate constants */
        if (rate_constant < 0.0)
        {
          if (neg_reaction == WARN_ERROR)
          {
            mcell_error("reaction rate constants should be zero or positive.");
            return 1;
          }
          else if (neg_reaction == WARN_WARN) {
            mcell_warn("negative reaction rate constant %f; setting to zero "
                       "and continuing.", rate_constant);
            rate_constant = 0.0;
          }
        }

        tp = CHECKED_MEM_GET(tv_rxn_mem,
                             "time-varying reaction rate constants");
        if (tp == NULL) {
          fclose(f);
          return 1;
        }
        tp->next = NULL;
        tp->path = path;
        tp->time = t / time_unit;
        tp->value = rate_constant;
#ifdef DEBUG
        valid_linecount++;
#endif

        if (rx->prob_t == NULL) {
          rx->prob_t = tp;
          tp2 = tp;
        } else {
          if (tp2 == NULL) {
            tp2 = tp;
            tp->next = rx->prob_t;
            rx->prob_t = tp;
          } else {
            if (tp->time < tp2->time)
              mcell_warn(
                  "In rate constants file '%s', line %d is out of sequence. "
                  "Resorting.", fname, linecount);
            tp->next = tp2->next;
            tp2->next = tp;
            tp2 = tp;
          }
        }
      }
    }

#ifdef DEBUG
    mcell_log("Read %d rate constants from file %s.", valid_linecount, fname);
#endif

    fclose(f);
  }
  return 0;
}

struct sym_entry *mcell_new_rxn_pathname(struct volume *state, char *name) {
  if ((retrieve_sym(name, state->rxpn_sym_table)) != NULL) {
    mcell_log("Named reaction pathway already defined: %s", name);
    return NULL;
  } else if ((retrieve_sym(name, state->mol_sym_table)) != NULL) {
    mcell_log("Named reaction pathway already defined as a molecule: %s", name);
    return NULL;
  }

  struct sym_entry *symp = store_sym(name, RXPN, state->rxpn_sym_table, NULL);
  if (symp == NULL) {
    mcell_log("Out of memory while creating reaction name: %s", name);
    return NULL;
  }
  return symp;
}
