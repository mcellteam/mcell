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

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 *
 * There is an exception, methods that schedule actions into the event's
 * new_diffuse_or_unimol_react_actions action queue stay in the diffuse_react_event_t
 * class.
 */

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

using namespace std;

namespace MCell {

namespace rx_util {

// ---------------------------------- bimolecular reactions ----------------------------------


/*************************************************************************
trigger_bimolecular:
   In: hash values of the two colliding molecules
       pointers to the two colliding molecules
       orientations of the two colliding molecules
         both zero away from a surface
         both nonzero (+-1) at a surface
       A is the moving molecule and B is the target
       array of pointers to the possible reactions
   Out: number of possible reactions for molecules reacA and reacB
        Also the first 'number' slots in the 'matching_rxns'
        array are filled with pointers to the possible reactions objects.
   Note: The target molecule is already scheduled and can be destroyed
         but not rescheduled.  Assume we have or will check separately that
         the moving molecule is not inert!
*************************************************************************/
static void trigger_bimolecular(
    const BimolecularReactionsMap& reactions,
    const Molecule& reacA, const Molecule& reacB,
    orientation_t orientA, orientation_t orientB,
    ReactionsVector& matching_rxns // items are appended
) {

  // for all reactions applicable to reacA and reacB
  BimolecularReactionsMap::const_iterator reactions_reacA_it
    = reactions.find(reacA.species_id);
  if (reactions_reacA_it == reactions.end()) {
    // no reactions at all for reacA
    return;
  }

  SpeciesReactionMap::const_iterator reactions_reacA_and_reacB_it
    = reactions_reacA_it->second.find(reacB.species_id);
  if (reactions_reacA_and_reacB_it == reactions_reacA_it->second.end()) {
    return;
  }

  // for now we can have just one reaction, TODO_PATHWAYS: allow multiple reactions..
  const Reaction* inter = reactions_reacA_and_reacB_it->second;

  /*for (struct rxn *inter = reaction_hash[hash]; inter != NULL; inter = inter->next)*/ {

    int right_walls_surf_classes = 0;  /* flag to check whether SURFACE_CLASSES
                                          of the walls for one or both reactants
                                          match the SURFACE_CLASS of the reaction
                                          (if needed) */

    /* skip irrelevant reactions (i.e. non vol-surf reactions) */
    assert(inter->reactants.size() == 2 && "We already checked that there must be 2 reactants");

    /* Check to see if orientation classes are zero/different */


    /* Check to see if orientation classes are zero/different */
    int test_wall = 0;
    orientation_t geomA = inter->reactants[0].orientation;
    orientation_t geomB = inter->reactants[1].orientation;
    if (geomA == ORIENTATION_NONE || geomB == ORIENTATION_NONE || (geomA + geomB) * (geomA - geomB) != 0) {
      matching_rxns.push_back(inter);
      //continue;
    }
    else if (orientA != ORIENTATION_NONE && orientA * orientB * geomA * geomB > 0) {
      matching_rxns.push_back(inter);
      //continue;
    }
  }
}

static void trigger_bimolecular_orientation_from_mols(
    const BimolecularReactionsMap& reactions,
    const Molecule& reacA, const Molecule& reacB,
    ReactionsVector& matching_rxns // items are appended
) {
  trigger_bimolecular(
      reactions,
      reacA, reacB,
      reacA.s.orientation, reacB.s.orientation,
      matching_rxns
  );
}


/*************************************************************************
binary_search_double

  In: A: A pointer to an array of doubles
      match: The value to match in the array
      max_idx: Initially, the size of the array
      mult: A multiplier for the comparison to the match.
            Set to 1 if not needed.
  Out: Returns the index of the match in the array
  Note: This should possibly be moved to util.c
*************************************************************************/
static int binary_search_double(const small_vector<float_t>& A, float_t match, int max_idx, double mult) {
  int min_idx = 0;

  while (max_idx - min_idx > 1) {
    int mid_idx = (max_idx + min_idx) / 2;
    if (match > (A[mid_idx] * mult)) {
      min_idx = mid_idx;
    }
    else {
      max_idx = mid_idx;
    }
  }

  if (match > A[min_idx]) {
    return max_idx;
  }
  else {
    return min_idx;
  }
}


/*************************************************************************
test_bimolecular
  In: the reaction we're testing
      a scaling coefficient depending on how many timesteps we've
        moved at once (1.0 means one timestep) and/or missing interaction area
      local probability factor (positive only for the reaction between two
        surface molecules, otherwise equal to zero)
      reaction partners
  Out: RX_NO_RX if no reaction occurs
       int containing which reaction pathway to take if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
*************************************************************************/
/*, double scaling - 1, double local_prob_factor = 0,*/
static int test_bimolecular(
    const Reaction& rx,
    rng_state& rng,
    const Molecule& a1, // unused for now
    const Molecule& a2, // unused for now
    const float_t scaling,
    const float_t local_prob_factor
) {
  /* rescale probabilities for the case of the reaction
     between two surface molecules */
  float_t min_noreaction_p = rx.min_noreaction_p; // local_prob_factor == 0
  if (local_prob_factor != 0) {
    min_noreaction_p = rx.min_noreaction_p * local_prob_factor;
  }

  float_t p;
  if (min_noreaction_p < scaling) {
    /* Instead of scaling rx->cum_probs array we scale random probability */
    p = rng_dbl(&rng) * scaling;

    if (p >= min_noreaction_p) {
      return RX_NO_RX;
    }
    else {
      return 0; // we have just one pathwayy
    }
  }
  else {
    float max_p = rx.max_fixed_p; //rx.cum_probs[rx->n_pathways - 1]; // TODO_PATHWAYS
    if (local_prob_factor > 0) {
      max_p *= local_prob_factor;
    }

    if (max_p >= scaling) /* we cannot scale enough. add missed rxns */
    {
      /* How may reactions will we miss? - only statistics?? */
      /*if (scaling == 0.0)
        rx->n_skipped += GIGANTIC4;
      else
        rx->n_skipped += (max_p / scaling) - 1.0;*/

      /* Keep the proportions of outbound pathways the same. */
      p = rng_dbl(&rng) * max_p;
    }
    else /* we can scale enough */
    {
      /* Instead of scaling rx->cum_probs array we scale random probability */
      p = rng_dbl(&rng) * scaling;

      if (p >= max_p)
        return RX_NO_RX;
    }
  }

  // TODO_PATHWAYS: find the right reaction for multiple pathways
  return 0; // return index of the only reaction
}



/*************************************************************************
test_many_bimolecular:
  In: an array of reactions we're testing
      scaling coefficients depending on how many timesteps we've moved
        at once (1.0 means one timestep) and/or missing interaction areas
      local probability factor for the corresponding reactions
      the number of elements in the array of reactions
      placeholder for the chosen pathway in the reaction (works as return
          value)
      a flag to indicate if
  Out: RX_NO_RX if no reaction occurs
       index in the reaction array corresponding to which reaction occurs
          if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
  Note: this uses only one call to get a random double, so you can't
        effectively sample events that happen less than 10^-9 of the
        time (for 32 bit random number).
  NOTE: This function was merged with test_many_bimolecular_all_neighbors.
        These two functions were almost identical, and the behavior of the
        "all_neighbors" version is preserved with a flag that can be passed in.
        For reactions between two surface molecules, set this flag to 1. For
        such reactions local_prob_factor > 0.
*************************************************************************/
static int test_many_bimolecular(
    const ReactionsVector& rx,
    const small_vector<float_t>& scaling,
    float_t local_prob_factor,
    //TODO_PATHWAYS: int& chosen_pathway - always 0 for now?
    rng_state& rng,
    bool all_neighbors_flag
) {
  assert(rx.size() == scaling.size());
  uint n = rx.size();

  small_vector<float_t> rxp;  /* array of cumulative rxn probabilities */
  rxp.resize(2 * n, 0.0);

  const Reaction* my_rx;
  int m, M;
  float_t p, f;

  if (all_neighbors_flag && local_prob_factor <= 0) {
    mcell_internal_error(
        "Local probability factor = %g in the function 'test_many_bimolecular_all_neighbors().",
        local_prob_factor
    );
  }

  if (n == 1) {
    Molecule dummy;
    if (all_neighbors_flag) {
      return test_bimolecular(*rx[0], rng, dummy, dummy, scaling[0], local_prob_factor);
    }
    else {
      return test_bimolecular(*rx[0], rng, dummy, dummy, scaling[0], 0);
    }
  }

  /* Note: lots of division here, if we're CPU-bound,could invert the
     definition of scaling_coefficients */
  if (all_neighbors_flag && local_prob_factor > 0) {
    rxp[0] = (rx[0]->max_fixed_p) * local_prob_factor / scaling[0];
  }
  else {
    rxp[0] = rx[0]->max_fixed_p / scaling[0];
  }

  for (uint i = 1; i < n; i++) {
    if (all_neighbors_flag && local_prob_factor > 0) {
      rxp[i] = rxp[i - 1] + (rx[i]->max_fixed_p) * local_prob_factor / scaling[i];
    }
    else {
      rxp[i] = rxp[i - 1] + rx[i]->max_fixed_p / scaling[i];
    }
  }

  if (rxp[n - 1] > 1.0) {
    //f = rxp[n - 1] - 1.0;   /* Number of failed reactions */
#if 0
    for (i = 0; i < n; i++) /* Distribute failures */
    {
      if (all_neighbors_flag && local_prob_factor > 0) {
        rx[i]->n_skipped += f * ((rx[i]->cum_probs[rx[i]->n_pathways - 1]) *
                                 local_prob_factor) /
                            rxp[n - 1];
      } else {
        rx[i]->n_skipped +=
            f * (rx[i]->cum_probs[rx[i]->n_pathways - 1]) / rxp[n - 1];
      }
    }
#endif
    p = rng_dbl(&rng) * rxp[n - 1];
  }
  else {
    p = rng_dbl(&rng);
    if (p > rxp[n - 1])
      return RX_NO_RX;
  }

  /* Pick the reaction that happens */
  int rx_index = binary_search_double(rxp, p, rxp.size() - 1, 1);

#if 0 // TODO_PATHWAYS
  my_rx = rx[i];
  if (i > 0) {
    p = (p - rxp[i - 1]);
  }
  p = p * scaling[i];

  /* Now pick the pathway within that reaction */
  M = my_rx->n_pathways - 1;

  if (all_neighbors_flag && local_prob_factor > 0)
    m = binary_search_double(my_rx->cum_probs, p, M, local_prob_factor);
  else
    m = binary_search_double(my_rx->cum_probs, p, M, 1);

  *chosen_pathway = m;
#endif

  return rx_index;
}


// might return nullptr if there is no unimolecular reaction for this species
// based on pick_unimolecular_reaction
static const Reaction* pick_unimol_rx(
    const World* world,
    const species_id_t species_id
) {
  const UnimolecularReactionsMap* unimol_rxs = world->world_constants.unimolecular_reactions_map;
  auto it = unimol_rxs->find(species_id);
  if (it == unimol_rxs->end()) {
    return nullptr;
  }
  else {
    return it->second;
  }
}


// based on timeof_unimolecular
static float_t time_of_unimol(const Reaction* rx, rng_state& rng) {
  double k_tot = rx->max_fixed_p;
  double p = rng_dbl(&rng);

  if ((k_tot <= 0) || (!distinguishable_f(p, 0, EPS)))
    return TIME_FOREVER;
  return -log(p) / k_tot;
}


// based on compute_lifetime
static float_t compute_unimol_lifetime(
    World* world,
    const Molecule& vm,
    const Reaction* rx
) {
  assert(rx != nullptr);

  float_t res = time_of_unimol(rx, world->rng);

#ifdef DEBUG_REACTIONS
  DUMP_CONDITION4(
      // calling rng for unimolecular
      vm.dump(world, "Assigned unimolecular time (prev rng):", "", world->current_iteration, res);
  );
#endif

  return res;
}


// note: this might be precomputed
static uint get_num_surface_products(const World* world, const Reaction* rx) {
  uint res = 0;
  for (const SpeciesWithOrientation& product: rx->products) {
    const Species& species = world->get_species(product.species_id);
    if (species.is_surf()) {
      res++;
    }
  }
  return res;
}


} // namespace rx_util

} // namespace mcell
