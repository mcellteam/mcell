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
#ifndef SRC4_REACTION_UTILS_INC_
#define SRC4_REACTION_UTILS_INC_

#include "bng/bng.h"

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

using namespace std;

namespace MCell {

namespace RxUtil {

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
    BNG::BNGEngine& bng_engine,
    const Molecule& reacA, const Molecule& reacB,
    orientation_t orientA, orientation_t orientB,
    BNG::RxnClassesVector& matching_rxn_classes // items are appended
) {
  const BNG::RxnClass* rxn_class = bng_engine.all_rxns.get_bimol_rxn_class(reacA.species_id, reacB.species_id);
  if (rxn_class == nullptr) {
    // no reaction
    return;
  }
  assert(rxn_class->reactants.size() == 2 && "We already checked that there must be 2 reactants");
  
  matching_rxn_classes.push_back(rxn_class);
  

  /* Check to see if orientation classes are zero/different */
  int test_wall = 0;
  orientation_t geomA = rxn_class->get_reactant_orientation(0);
  orientation_t geomB = rxn_class->get_reactant_orientation(1);
  if (geomA == ORIENTATION_NONE || geomB == ORIENTATION_NONE || (geomA + geomB) * (geomA - geomB) != 0) {
    matching_rxn_classes.push_back(rxn_class);
  }
  else if (orientA != ORIENTATION_NONE && orientA * orientB * geomA * geomB > 0) {
    matching_rxn_classes.push_back(rxn_class);
  }
}

static void trigger_bimolecular_orientation_from_mols(
    BNG::BNGEngine& bng_engine,
    const Molecule& reacA, const Molecule& reacB,
    BNG::RxnClassesVector& matching_rxn_classes // items are appended
) {
  trigger_bimolecular(
      bng_engine,
      reacA, reacB,
      reacA.s.orientation, reacB.s.orientation,
      matching_rxn_classes
  );
}



/*************************************************************************
 *
 * find all surface reactions for any surface molecule with orientation
 * orientA on a surface class triggered via the ALL_MOLECULES and
 * ALL_SURFACE_MOLECULE keywords
 *
 * in: orientation of surface molecule
 *     surface class species to test
 *     number of matching reactions before the function call
 *     flag signaling the presence of a transparent region border
 *     flag signaling the presence of a reflective region border
 *     flag signaling the presence of a absorbing region border
 *     array holding matching reactions
 *
 * out: returns number of matching reactions
 *      adds matching reactions to matching_rxns array
 *
 *************************************************************************/

// to be called only from find_surface_mol_reactions_with_surf_classes
static void find_surface_mol_reactions_with_surf_classes_for_rxn_class_map(
    const Partition& p,
    const Molecule& reacA,
    const bool allow_rx_transp,
    const bool allow_rx_reflec,
    const bool allow_rx_absorb_reg_border,
    const Region& reg,
    const BNG::SpeciesRxnClassesMap& potential_reactions,
    BNG::RxnClassesVector& matching_rxn_classes
) {

  BNG::SpeciesRxnClassesMap::const_iterator reactions_reacA_and_surface_it = potential_reactions.find(reg.species_id);
  if (reactions_reacA_and_surface_it == potential_reactions.end()) {
    // no reactions for this type of region
    return;
  }

  const BNG::RxnClass* rxn_class = reactions_reacA_and_surface_it->second;
  assert(rxn_class->reactants.size() == 2);

  orientation_t orient0 = rxn_class->get_reactant_orientation(0);
  orientation_t orient1 = rxn_class->get_reactant_orientation(1);

  // TODO: can we move this condition to some shared function?
  if ( (orient0 == 0) ||
       (orient1 == 0 || (orient0 + orient1) * (orient0 - orient1) != 0) ||
       (reacA.s.orientation * orient0 * orient1 > 0) ) {

    bool reject_reaction =
        (rxn_class->is_transparent() && !allow_rx_transp) ||
        (rxn_class->is_reflect() && !allow_rx_reflec) ||
        (rxn_class->is_absorb() && !allow_rx_absorb_reg_border);

    if (!reject_reaction) {
      matching_rxn_classes.push_back(rxn_class);
    }
  }
}


static void find_surface_mol_reactions_with_surf_classes(
    const Partition& p,
    const Molecule& reacA,
    const Wall& w,
    const bool allow_rx_transp,
    const bool allow_rx_reflec,
    const bool allow_rx_absorb_reg_border,
    BNG::RxnClassesVector& matching_rxns
) {

  assert(reacA.is_surf());

  // for all reactions applicable to reacA and and wall
  const BNG::SpeciesRxnClassesMap& species_rxns = p.bng_engine.all_rxns.get_bimol_rxns_for_reactant(reacA.species_id);
  const BNG::SpeciesRxnClassesMap& all_molecules_rxns = p.bng_engine.all_rxns.get_bimol_rxns_for_reactant(p.bng_engine.all_rxns.all_molecules_species_id);
  const BNG::SpeciesRxnClassesMap& all_vol_surf_rxns = p.bng_engine.all_rxns.get_bimol_rxns_for_reactant(p.bng_engine.all_rxns.all_surface_molecules_species_id);

  if (species_rxns.empty() && all_molecules_rxns.empty() && all_vol_surf_rxns.empty()) {
    // no reactions at all for reacA
    return;
  }

  // for all reactive regions of a wall
  for (region_index_t region_index: w.regions) {
    const Region& reg = p.get_region(region_index);
    if (!reg.is_reactive()) {
      continue;
    }

    find_surface_mol_reactions_with_surf_classes_for_rxn_class_map(
        p, reacA, allow_rx_transp, allow_rx_reflec, allow_rx_absorb_reg_border, reg,
        species_rxns,
        matching_rxns
    );

    find_surface_mol_reactions_with_surf_classes_for_rxn_class_map(
        p, reacA, allow_rx_transp, allow_rx_reflec, allow_rx_absorb_reg_border, reg,
        all_molecules_rxns,
        matching_rxns
    );

    find_surface_mol_reactions_with_surf_classes_for_rxn_class_map(
        p, reacA, allow_rx_transp, allow_rx_reflec, allow_rx_absorb_reg_border, reg,
        all_vol_surf_rxns,
        matching_rxns
    );
}
}


/*************************************************************************
trigger_intersect:
   In: hash value of molecule's species
       pointer to a molecule
       orientation of that molecule
       pointer to a wall
       array of matching reactions (placeholder for output)
       flags that tells whether we should include special reactions
          (REFL/TRANSP/ABSORB_REGION_BORDER) in the output array
   Out: number of matching reactions for this
        molecule/wall intersection, or for this mol/generic wall,
        or this wall/generic mol.  All matching reactions are placed in
        the array "matching_rxns" in the first "number" slots.
   Note: Moving molecule may be inert.

*************************************************************************/
static void trigger_intersect(
    const Partition& p,
    const Molecule& reacA,
    const Wall& w,
    BNG::RxnClassesVector& matching_rxns,
    const bool allow_rx_transp = true,
    const bool allow_rx_reflec = true,
    const bool allow_rx_absorb_reg_border = true
) {
  /*
   * ??? how are unimol reactions related to molecule - surface reactions?
  if (w.in_reactive_region()) {
    num_matching_rxns = find_unimol_reactions_with_surf_classes(
        reaction_hash, rx_hashsize, reacA, w, hashA, orientA, num_matching_rxns,
        allow_rx_transp, allow_rx_reflec, allow_rx_absorb_reg_border,
        matching_rxns);
  }
  */

  for (region_index_t region_index: w.regions) {

    if (reacA.is_vol()) {
      /*num_matching_rxns = find_volume_mol_reactions_with_surf_classes(
          reaction_hash, rx_hashsize, all_mols, all_volume_mols, orientA,
          scl->surf_class, num_matching_rxns, allow_rx_transp, allow_rx_reflec,
          matching_rxns);*/
      assert(false && "TODO");
    }
    else if (reacA.is_surf()) {

      find_surface_mol_reactions_with_surf_classes(
          p, reacA, w,  allow_rx_transp, allow_rx_reflec,
          allow_rx_absorb_reg_border, matching_rxns);
    }
  }
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
static int binary_search_double(const std::vector<float_t>& A, float_t match, int max_idx, double mult) {
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

  if (match > A[min_idx] * mult) {
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
static int test_bimolecular(
    const BNG::RxnClass* rxn_class,
    rng_state& rng,
    const Molecule& a1, // unused for now
    const Molecule& a2, // unused for now
    const float_t scaling,
    const float_t local_prob_factor
) {
  assert(rxn_class != nullptr);
  /* rescale probabilities for the case of the reaction
     between two surface molecules */
  float_t min_noreaction_p = rxn_class->min_noreaction_p; // local_prob_factor == 0
  if (local_prob_factor != 0) {
    min_noreaction_p = rxn_class->min_noreaction_p * local_prob_factor;
  }

  float_t p;
  if (min_noreaction_p < scaling) {
    /* Instead of scaling rx->cum_probs array we scale random probability */
    p = rng_dbl(&rng) * scaling;

    if (p >= min_noreaction_p) {
      return RX_NO_RX;
    }
    // continue below
  }
  else {
    float max_p = rxn_class->max_fixed_p; //rx.cum_probs[rx->n_pathways - 1]; // TODO_PATHWAYS
    if (local_prob_factor > 0) {
      max_p *= local_prob_factor;
    }

    if (max_p >= scaling) /* we cannot scale enough. add missed rxns */
    {
      // TODO: this statistics is important
      // where to store this info? - RXN is constant and should stay like this
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

#ifdef DEBUG_REACTION_PROBABILITIES
  mcell_log(
      "test_bimolecular: p = %.8f, scaling = %.8f, min_noreaction_p = %.8f, local_prob_factor = %.8f",
      p, scaling, min_noreaction_p, local_prob_factor
  );
#endif

  int M = rxn_class->get_num_reactions() - 1;
  if (local_prob_factor > 0) {
    return binary_search_double(rxn_class->cum_probs, p, M, local_prob_factor);
  }
  else {
    return binary_search_double(rxn_class->cum_probs, p, M, 1);
  }
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
    const BNG::RxnClassesVector& rxn_classes,
    const small_vector<float_t>& scaling,
    float_t local_prob_factor,
    rng_state& rng,
    bool all_neighbors_flag,
    reaction_index_t& chosen_reaction_index
) {
  assert(rxn_classes.size() == scaling.size());
  uint n = rxn_classes.size();

  //TODO: use small vector
  /* array of cumulative rxn probabilities */
  std::vector<float_t> cum_rxn_class_probs; // rxn in mcell3
  cum_rxn_class_probs.resize(2 * n, 0.0);

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
      return test_bimolecular(rxn_classes[0], rng, dummy, dummy, scaling[0], local_prob_factor);
    }
    else {
      return test_bimolecular(rxn_classes[0], rng, dummy, dummy, scaling[0], 0);
    }
  }

  /* Note: lots of division here, if we're CPU-bound,could invert the
     definition of scaling_coefficients */
  if (all_neighbors_flag && local_prob_factor > 0) {
    cum_rxn_class_probs[0] = (rxn_classes[0]->max_fixed_p) * local_prob_factor / scaling[0];
  }
  else {
    cum_rxn_class_probs[0] = rxn_classes[0]->max_fixed_p / scaling[0];
  }

  for (uint i = 1; i < n; i++) {
    if (all_neighbors_flag && local_prob_factor > 0) {
      cum_rxn_class_probs[i] = cum_rxn_class_probs[i - 1] + (rxn_classes[i]->max_fixed_p) * local_prob_factor / scaling[i];
    }
    else {
      cum_rxn_class_probs[i] = cum_rxn_class_probs[i - 1] + rxn_classes[i]->max_fixed_p / scaling[i];
    }
  }

  if (cum_rxn_class_probs[n - 1] > 1.0) {
    //f = rxp[n - 1] - 1.0;   /* Number of failed reactions */
    // TODO: important statistics
#if 0
    for (i = 0; i < n; i++) /* Distribute failures */
    {
      if (all_neighbors_flag && local_prob_factor > 0) {
        rxn_classes[i]->n_skipped += f * ((rxn_classes[i]->cum_probs[rxn_classes[i]->n_pathways - 1]) *
                                 local_prob_factor) /
                            cum_rxn_class_probs[n - 1];
      } else {
        rxn_classes[i]->n_skipped +=
            f * (rxn_classes[i]->cum_probs[rxn_classes[i]->n_pathways - 1]) / cum_rxn_class_probs[n - 1];
      }
    }
#endif
    p = rng_dbl(&rng) * cum_rxn_class_probs[n - 1];
  }
  else {
    p = rng_dbl(&rng);
    if (p > cum_rxn_class_probs[n - 1])
      return RX_NO_RX;
  }

  /* Pick the reaction class that happens */
  int rx_index = binary_search_double(cum_rxn_class_probs, p, cum_rxn_class_probs.size() - 1, 1);
  assert(rx_index >= 0);

  const BNG::RxnClass* my_rxn = rxn_classes[rx_index];
  if (rx_index > 0) {
    p = (p - cum_rxn_class_probs[rx_index - 1]);
  }
  p = p * scaling[rx_index];

  /* Now pick the pathway within that reaction */
  // NOTE: might optimize if there is just one rxn
  M = my_rxn->get_num_reactions() - 1;
  if (all_neighbors_flag && local_prob_factor > 0) {
    m = binary_search_double(my_rxn->cum_probs, p, M, local_prob_factor);
  }
  else {
    m = binary_search_double(my_rxn->cum_probs, p, M, 1);
  }

  chosen_reaction_index = m;

  return rx_index;
}


// might return nullptr if there is no unimolecular reaction for this species
// based on pick_unimolecular_reaction
static const BNG::RxnClass* pick_unimol_rx(
    World* world,
    const species_id_t species_id
) {
  return world->bng_engine.all_rxns.get_unimol_rxn_class(species_id);
}


// based on timeof_unimolecular
static float_t time_of_unimol(const BNG::RxnClass* rx, rng_state& rng) {
  double k_tot = rx->max_fixed_p;
  double p = rng_dbl(&rng);

  if ((k_tot <= 0) || (!distinguishable_f(p, 0, EPS)))
    return TIME_FOREVER;
  return -log(p) / k_tot;
}


// based on compute_lifetime
static float_t compute_unimol_lifetime(
    const Partition& p,
    rng_state& rng,
    const Molecule& vm,
    const BNG::RxnClass* rx
) {
  assert(rx != nullptr);

  float_t res = time_of_unimol(rx, rng);

#ifdef DEBUG_REACTIONS
  SimulationStats* world = &p.stats;
  DUMP_CONDITION4(
      // calling rng for unimolecular
      vm.dump(p, "Assigned unimolecular time (prev rng):", "", p.stats.get_current_iteration(), res);
  );
#endif

  return res;
}

/*************************************************************************
which_unimolecular:
  In: the reaction we're testing
  Out: int containing which unimolecular reaction occurs (one must occur)
*************************************************************************/
static reaction_index_t which_unimolecular(const BNG::RxnClass *rxn_class, rng_state& rng) {
  assert(rxn_class != nullptr);
  if (rxn_class->get_num_reactions() == 1) {
    return 0;
  }

  int max = rxn_class->get_num_reactions() - 1;
  float_t match = rng_dbl(&rng);
  match = match * rxn_class->cum_probs[max];
  return binary_search_double(rxn_class->cum_probs, match, max, 1);
}

} // namespace RxUtil

} // namespace MCell

#endif // SRC4_REACTION_UTILS_INC_
