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
#ifndef SRC4_RXN_UTILS_INC_
#define SRC4_RXN_UTILS_INC_

#include "bng/bng.h"

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

using namespace std;

namespace MCell {
namespace RxnUtils {

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
  BNG::RxnClass* rxn_class =
      bng_engine.get_all_rxns().get_bimol_rxn_class(reacA.species_id, reacB.species_id);
  if (rxn_class == nullptr) {
    // no reaction
    return;
  }
  assert(rxn_class->is_bimol() && "We already checked that there must be 2 reactants");

  /* Check to see if orientation classes are zero/different */
  int test_wall = 0;
  orientation_t geomA = rxn_class->get_reactant_orientation(0);
  orientation_t geomB = rxn_class->get_reactant_orientation(1);
  if (geomA == ORIENTATION_NONE || geomA == ORIENTATION_DEPENDS_ON_SURF_COMP ||
      geomB == ORIENTATION_NONE || geomB == ORIENTATION_DEPENDS_ON_SURF_COMP ||
      (geomA + geomB) * (geomA - geomB) != 0) {
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
static void find_reactions_with_surf_classes_for_rxn_class_map(
    const Partition& p,
    const Molecule& reacA,
	  const orientation_t reacA_orient,
    const Region& reg,
    const BNG::SpeciesRxnClassesMap& potential_reactions,
    const bool allow_rx_transp_reflec_absorb_reg_border,
    BNG::RxnClassesVector& matching_rxn_classes
) {

  auto reactions_reacA_and_surface_it = potential_reactions.find(reg.species_id);
  if (reactions_reacA_and_surface_it == potential_reactions.end()) {
    // no reactions for this type of region
    return;
  }

  // NOTE: do we need to handle compartments here?
  BNG::RxnClass* rxn_class = reactions_reacA_and_surface_it->second;

  if (!allow_rx_transp_reflec_absorb_reg_border &&
      rxn_class->is_reflect_transparent_or_absorb_region_border()
  ) {
    // do not allow this type
    return;
  }

  assert(rxn_class->is_bimol());
  if (rxn_class->get_num_reactions() == 0) {
    return;
  }

  orientation_t orient0 = rxn_class->get_reactant_orientation(0);
  orientation_t orient1 = rxn_class->get_reactant_orientation(1);

  // TODO: can we move this condition to some shared function?
  if ( (orient0 == 0) ||
       (orient1 == 0 || (orient0 + orient1) * (orient0 - orient1) != 0) ||
       (reacA_orient * orient0 * orient1 > 0) ) {

    matching_rxn_classes.push_back(rxn_class);
  }
}


/*************************************************************************
 *
 * find all volume reactions for any volume molecule with orientation
 * orientA with a surface class triggered via the ALL_MOLECULES and
 * ALL_VOLUME_MOLECULE keywords
 *
 * in: orientation of surface molecule
 *     surface class species to test
 *     number of matching reactions before the function call
 *     flag signalling the presence of transparent region border
 *     flag signalling the presence of a reflective region border
 *     flag signalling the presence of a absorbing region border
 *     array holding matching reactions
 *
 * out: returns number of matching reactions
 *      adds matching reactions to matching_rxns array
 *
 *************************************************************************/
// TODO LATER: practically the same as find_surface_mol_reactions_with_surf_classes, but let's wait for full vol/surface rxn support
static void find_volume_mol_reactions_with_surf_classes(
    Partition& p,
    const Molecule& reacA,
	  const orientation_t reacA_orient,
    const Wall& w,
    const bool allow_rx_transp_reflec_absorb_reg_border,
    BNG::RxnClassesVector& matching_rxns
) {

  assert(reacA.is_vol());

  // for all reactions applicable to reacA and and wall
  const BNG::SpeciesRxnClassesMap* species_rxns =
      p.get_all_rxns().get_bimol_rxns_for_reactant(reacA.species_id);
  const BNG::SpeciesRxnClassesMap* all_molecules_rxns =
      p.get_all_rxns().get_bimol_rxns_for_reactant(p.get_all_species().get_all_molecules_species_id());
  const BNG::SpeciesRxnClassesMap* all_vol_rxns =
      p.get_all_rxns().get_bimol_rxns_for_reactant(p.get_all_species().get_all_volume_molecules_species_id());

  if (species_rxns == nullptr && all_molecules_rxns == nullptr && all_vol_rxns == nullptr) {
    // no reactions at all for reacA
    return;
  }

  // for all reactive regions of a wall
  for (region_index_t region_index: w.regions) {
    const Region& reg = p.get_region(region_index);
    if (!reg.is_reactive()) {
      continue;
    }

    if (species_rxns != nullptr) {
      find_reactions_with_surf_classes_for_rxn_class_map(
          p, reacA, reacA_orient, reg,
          *species_rxns,
          allow_rx_transp_reflec_absorb_reg_border,
          matching_rxns
      );
    }

    if (all_molecules_rxns != nullptr) {
      find_reactions_with_surf_classes_for_rxn_class_map(
          p, reacA, reacA_orient, reg,
          *all_molecules_rxns,
          allow_rx_transp_reflec_absorb_reg_border,
          matching_rxns
      );
    }

    if (all_vol_rxns != nullptr) {
      find_reactions_with_surf_classes_for_rxn_class_map(
          p, reacA, reacA_orient, reg,
          *all_vol_rxns,
          allow_rx_transp_reflec_absorb_reg_border,
          matching_rxns
      );
    }
  }
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
 *     flag signalling the presence of transparent region border
 *     flag signalling the presence of a reflective region border
 *     flag signalling the presence of a absorbing region border
 *     array holding matching reactions
 *
 * out: returns number of matching reactions
 *      adds matching reactions to matching_rxns array
 *
 *************************************************************************/
template<typename WallOrObj>
void find_surface_mol_reactions_with_surf_classes(
    Partition& p,
    const Molecule& reacA,
    const WallOrObj& regions_list,
    const bool allow_rx_transp_reflec_absorb_reg_border,
    BNG::RxnClassesVector& matching_rxns
) {

  assert(reacA.is_surf());

  const BNG::SpeciesRxnClassesMap* species_rxns =
      p.get_all_rxns().get_bimol_rxns_for_reactant(reacA.species_id);
  const BNG::SpeciesRxnClassesMap* all_molecules_rxns =
      p.get_all_rxns().get_bimol_rxns_for_reactant(p.get_all_species().get_all_molecules_species_id());
  const BNG::SpeciesRxnClassesMap* all_surf_rxns =
      p.get_all_rxns().get_bimol_rxns_for_reactant(p.get_all_species().get_all_surface_molecules_species_id());

  if (species_rxns == nullptr && all_molecules_rxns == nullptr && all_surf_rxns == nullptr) {
    // no reactions at all for reacA
    return;
  }

  // for all reactive regions of a wall
  for (region_index_t region_index: regions_list.regions) {
    const Region& reg = p.get_region(region_index);
    if (!reg.is_reactive()) {
      continue;
    }

    if (species_rxns != nullptr) {
      find_reactions_with_surf_classes_for_rxn_class_map(
          p, reacA, reacA.s.orientation, reg,
          *species_rxns,
          allow_rx_transp_reflec_absorb_reg_border,
          matching_rxns
      );
    }

    if (all_molecules_rxns != nullptr) {
      find_reactions_with_surf_classes_for_rxn_class_map(
          p, reacA, reacA.s.orientation, reg,
          *all_molecules_rxns,
          allow_rx_transp_reflec_absorb_reg_border,
          matching_rxns
      );
    }

    if (all_surf_rxns != nullptr) {
      find_reactions_with_surf_classes_for_rxn_class_map(
          p, reacA, reacA.s.orientation, reg,
          *all_surf_rxns,
          allow_rx_transp_reflec_absorb_reg_border,
          matching_rxns
      );
    }
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
    Partition& p,
    const Molecule& reacA,
	  const orientation_t reacA_orient,
    const Wall& w,
    const bool allow_rx_transp_reflec_absorb_reg_border,
    BNG::RxnClassesVector& matching_rxns
) {
  /*if (w.in_reactive_region()) {
    find_surface_mol_reactions_with_surf_classes(
        p, reacA, w, matching_rxns, false
        reaction_hash, rx_hashsize, reacA, w, hashA, orientA, num_matching_rxns,
        allow_rx_transp, allow_rx_reflec, allow_rx_absorb_reg_border,
        matching_rxns);
  }*/


  if (reacA.is_vol()) {
    find_volume_mol_reactions_with_surf_classes(
        p, reacA, reacA_orient, w, allow_rx_transp_reflec_absorb_reg_border,
        matching_rxns
    );
  }
  else if (reacA.is_surf()) {
	  assert(reacA.s.orientation == reacA_orient);
    find_surface_mol_reactions_with_surf_classes(
        p, reacA, w, allow_rx_transp_reflec_absorb_reg_border,
        matching_rxns
    );
  }
  else {
    assert(false);
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
static int binary_search_double(const std::vector<double>& A, double match, int max_idx, double mult) {
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
  Out: PATHWAY_INDEX_NO_RXN if no reaction occurs
       int containing which reaction pathway to take if one does occur
  Note: If this reaction does not return PATHWAY_INDEX_NO_RXN, then we update
        counters appropriately assuming that the reaction does take place.
*************************************************************************/
static int test_bimolecular(
    BNG::RxnClass* rxn_class,
    rng_state& rng,
    const Molecule& a1, // unused for now
    const Molecule& a2, // unused for now
    const double scaling,
    const double local_prob_factor,
    const double current_time
) {
  assert(rxn_class != nullptr);
  assert(rxn_class->get_num_reactions() != 0);

  rxn_class->update_rxn_rates_if_needed(current_time);

  /* rescale probabilities for the case of the reaction
     between two surface molecules */
  double max_fixed_p = rxn_class->get_max_fixed_p(); // local_prob_factor == 0
  if (local_prob_factor != 0) {
    max_fixed_p = rxn_class->get_max_fixed_p() * local_prob_factor;
  }

  double p;
  if (max_fixed_p < scaling) {
    /* Instead of scaling rx->cum_probs array we scale random probability */
    p = rng_dbl(&rng) * scaling;

    if (p >= max_fixed_p) {
      return BNG::PATHWAY_INDEX_NO_RXN;
    }
    // continue below
  }
  else {
    float max_p = rxn_class->get_max_fixed_p(); //rx.cum_probs[rx->n_pathways - 1]; // TODO_PATHWAYS
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
        return BNG::PATHWAY_INDEX_NO_RXN;
    }
  }

#ifdef DEBUG_REACTION_PROBABILITIES
  mcell_log(
      "test_bimolecular: p = %.8f, scaling = %.8f, max_fixed_p = %.8f, local_prob_factor = %.8f",
      p, scaling, max_fixed_p, local_prob_factor
  );
#endif

  if (local_prob_factor > 0) {
    return rxn_class->get_pathway_index_for_probability(p, local_prob_factor);
  }
  else {
    return rxn_class->get_pathway_index_for_probability(p, 1);
  }
}


/*************************************************************************
test_many_unimol:
  In: an array of reactions we're testing
      rng state
  Out: index of the selected rxn class
*************************************************************************/
static uint test_many_unimol(
    const BNG::RxnClassesVector& rxn_classes,
    rng_state& rng) {

  assert(!rxn_classes.empty());

  size_t n = rxn_classes.size();

  if (n == 1) {
    return 0;
  }

  std::vector<double> cum_rxn_class_probs(n); /* array of cumulative rxn probabilities */
  cum_rxn_class_probs[0] = rxn_classes[0]->get_max_fixed_p();

  for (size_t i = 1; i < n; i++) {
    cum_rxn_class_probs[i] = cum_rxn_class_probs[i - 1] + rxn_classes[i]->get_max_fixed_p();
  }

  double p = rng_dbl(&rng) * cum_rxn_class_probs[n - 1];

  /* Pick the reaction that happens */
  uint res = binary_search_double(cum_rxn_class_probs, p, cum_rxn_class_probs.size() - 1, 1);

  return res;
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
  Out: PATHWAY_INDEX_NO_RXN if no reaction occurs
       index in the reaction array corresponding to which reaction occurs
          if one does occur
  Note: If this reaction does not return PATHWAY_INDEX_NO_RXN, then we update
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
    BNG::RxnClassesVector& rxn_classes,
    const small_vector<double>& scaling,
    const double local_prob_factor,
    rng_state& rng,
    const bool all_neighbors_flag,
    const double current_time,
    BNG::rxn_class_pathway_index_t& chosen_pathway_index
) {
  assert(rxn_classes.size() == scaling.size());
  uint n = rxn_classes.size();

  for (BNG::RxnClass* cls: rxn_classes) {
    cls->update_rxn_rates_if_needed(current_time);
  }

  /* array of cumulative rxn probabilities */
  std::vector<double> cum_rxn_class_probs; // rxn in mcell3
  cum_rxn_class_probs.resize(2 * n, 0.0);

  int m;
  double p;

  if (all_neighbors_flag && local_prob_factor <= 0) {
    mcell_internal_error(
        "Local probability factor = %g in the function 'test_many_bimolecular_all_neighbors().",
        local_prob_factor
    );
  }

  if (n == 1) {
    Molecule dummy;
    if (all_neighbors_flag) {
      return test_bimolecular(rxn_classes[0], rng, dummy, dummy, scaling[0], local_prob_factor, current_time);
    }
    else {
      return test_bimolecular(rxn_classes[0], rng, dummy, dummy, scaling[0], 0, current_time);
    }
  }

  /* Note: lots of division here, if we're CPU-bound,could invert the
     definition of scaling_coefficients */
  if (all_neighbors_flag && local_prob_factor > 0) {
    cum_rxn_class_probs[0] = (rxn_classes[0]->get_max_fixed_p()) * local_prob_factor / scaling[0];
  }
  else {
    cum_rxn_class_probs[0] = rxn_classes[0]->get_max_fixed_p() / scaling[0];
  }

  for (uint i = 1; i < n; i++) {
    if (all_neighbors_flag && local_prob_factor > 0) {
      cum_rxn_class_probs[i] = cum_rxn_class_probs[i - 1] + (rxn_classes[i]->get_max_fixed_p()) * local_prob_factor / scaling[i];
    }
    else {
      cum_rxn_class_probs[i] = cum_rxn_class_probs[i - 1] + rxn_classes[i]->get_max_fixed_p() / scaling[i];
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
    if (p > cum_rxn_class_probs[n - 1]) {
      return BNG::PATHWAY_INDEX_NO_RXN;
    }
  }

  /* Pick the reaction class that happens */
  int rx_index = binary_search_double(cum_rxn_class_probs, p, cum_rxn_class_probs.size() - 1, 1);
  assert(rx_index >= 0);

  BNG::RxnClass* selected_rxn_class = rxn_classes[rx_index];
  if (rx_index > 0) {
    p = (p - cum_rxn_class_probs[rx_index - 1]);
  }
  p = p * scaling[rx_index];

  /* Now pick the pathway within that reaction */
  // NOTE: might optimize if there is just one rxn
  if (all_neighbors_flag && local_prob_factor > 0) {
    m = selected_rxn_class->get_pathway_index_for_probability(p, local_prob_factor);
  }
  else {
    m = selected_rxn_class->get_pathway_index_for_probability(p, 1);
  }

  chosen_pathway_index = m;

  return rx_index;
}


/*************************************************************************
test_intersect
  In: the reaction we're testing
      a probability multiplier depending on how many timesteps we've
        moved at once (1.0 means one timestep)
  Out: PATHWAY_INDEX_NO_RXN if no reaction occurs (assume reflection)
       int containing which reaction occurs if one does occur
  Note: If not PATHWAY_INDEX_NO_RXN, and not the trasparency shortcut, then we
        update counters assuming the reaction will take place.
*************************************************************************/
// TODO LATER: see if it can be merged with test_bimolecular
static BNG::rxn_class_pathway_index_t test_intersect(
    BNG::RxnClass* rxn_class,
    const double scaling,
    const double current_time,
    rng_state& rng) {
  double p;

  assert(rxn_class->type == BNG::RxnType::Standard &&
      "Reflect and Transparent should be handled elsewhere, AbsorbRegionBorder is not applicable here");

  rxn_class->update_rxn_rates_if_needed(current_time);

  double max_prob = rxn_class->get_max_fixed_p();

  if (max_prob > scaling) {
    /*if (scaling <= 0.0)
      rxn_class->n_skipped += GIGANTIC;
    else
      rxn_class->n_skipped += rxn_class->cum_probs[rxn_class->n_pathways - 1] / scaling - 1.0;*/
    p = rng_dbl(&rng) * max_prob;
  }
  else {
    p = rng_dbl(&rng) * scaling;

    if (p > max_prob) {
      return BNG::PATHWAY_INDEX_NO_RXN;
    }
  }

  if (p > rxn_class->get_max_fixed_p()) {
    return BNG::PATHWAY_INDEX_NO_RXN;
  }

  double match = rng_dbl(&rng);
  match = match * rxn_class->get_max_fixed_p();

  return rxn_class->get_pathway_index_for_probability(match, 1);
}

/*************************************************************************
test_many_intersect:
  In: an array of reactions we're testing
      a probability multiplier depending on how many timesteps we've
        moved at once (1.0 means one timestep)
      the number of elements in the array of reactions
      placeholder for the chosen pathway in the reaction (return value)
  Out: PATHWAY_INDEX_NO_RXN if no reaction occurs (assume reflection)
       index in the reaction array if reaction does occur
  Note: If not PATHWAY_INDEX_NO_RXN, and not the trasparency shortcut, then we
        update counters assuming the reaction will take place.
*************************************************************************/
// TODO LATER: see if it can be merged with test_many bimolecular
static BNG::rxn_class_pathway_index_t test_many_intersect(
    BNG::RxnClassesVector& rxn_classes,
    const double scaling,
    const double current_time,
    BNG::rxn_class_index_t& selected_rxn_class_index,
    rng_state& rng) {

  for (BNG::RxnClass* cls: rxn_classes) {
    cls->update_rxn_rates_if_needed(current_time);
  }

  uint num_classes = rxn_classes.size();

  if (num_classes == 1) {
    selected_rxn_class_index = 0;
    return test_intersect(rxn_classes[0], scaling, current_time, rng);
  }

  // array of cumulative rxn probabilities
  std::vector<double> rxn_probs;
  rxn_probs.resize(num_classes);

  rxn_probs[0] = rxn_classes[0]->get_max_fixed_p() / scaling;
  for (uint i = 1; i < num_classes; i++) {
    rxn_probs[i] = rxn_probs[i - 1] + rxn_classes[i]->get_max_fixed_p() / scaling;
  }

  double p;
  if (rxn_probs[num_classes - 1] > 1.0) {
    /*double f = rxp[n - 1] - 1.0; // Number of failed reactions
    for (i = 0; i < n; i++)      // Distribute failures
    {
      rxn_classes[i]->n_skipped +=
          f * (rxn_classes[i]->cum_probs[rxn_classes[i]->n_pathways - 1]) / rxp[n - 1];
    }
    */
    p = rng_dbl(&rng) * rxn_probs[num_classes - 1];
  } else {
    p = rng_dbl(&rng);
    if (p > rxn_probs[num_classes - 1]) {
      selected_rxn_class_index = BNG::RNX_CLASS_INDEX_INVALID;
      return BNG::PATHWAY_INDEX_NO_RXN;
    }
  }

  /* Pick the reaction that happens */
  selected_rxn_class_index = binary_search_double(rxn_probs, p, num_classes - 1, 1);

  BNG::RxnClass *selected_rxn_class = rxn_classes[selected_rxn_class_index];

  if (selected_rxn_class_index > 0) {
    p = (p - rxn_probs[selected_rxn_class_index - 1]);
  }
  p = p * scaling;

  /* Now pick the pathway within that reaction */
  return selected_rxn_class->get_pathway_index_for_probability(p, 1);
}


// might return nullptr if there is no unimolecular reaction for this species
// based on pick_unimolecular_reaction
static void pick_unimol_rxn_classes(
    Partition& p,
    const Molecule& m,
    const double current_time,
    BNG::RxnClassesVector& matching_rxn_classes
) {
  // MCell3 returns mol+surf class rxn(s) as the first one(s), then the true unimol rxns
  // maintaining order only for compatibility
  const BNG::Species& species = p.get_species(m.species_id);
  if (species.has_flag(BNG::SPECIES_FLAG_CAN_SURFWALL)) {
    assert(m.is_surf());
    const Wall& w = p.get_wall(m.s.wall_index);
    trigger_intersect(p, m, m.s.orientation, w, false, matching_rxn_classes);
  }

  // get unimol rxn class if exists
  BNG::RxnClass* rxn_class = p.get_all_rxns().get_unimol_rxn_class(m.species_id);
  if (rxn_class != nullptr) {
    rxn_class->update_rxn_rates_if_needed(current_time);
    matching_rxn_classes.push_back(rxn_class);
  }
}


// based on timeof_unimolecular
static double time_of_unimol(BNG::RxnClass* rxn_class, rng_state& rng) {
  double k_tot = rxn_class->get_max_fixed_p();
#ifdef MCELL4_NO_RNG_FOR_UNIMOL_RXN_P_0
  if (k_tot <= 0) {
    // not calling random generator when p is 0 - for bompatibility with MCell3R
    return TIME_FOREVER;
  }
#endif  
  double p = rng_dbl(&rng);

  if ((k_tot <= 0) || (!distinguishable_f(p, 0, EPS))) {
    return TIME_FOREVER;
  }
  return -log(p) / k_tot;
}


// based on compute_lifetime
static double compute_unimol_lifetime(
    const Partition& p,
    rng_state& rng,
    BNG::RxnClass* rx,
    const double current_time,
    Molecule& m
) {
  assert(rx != nullptr);

  double unimol_time_from_now = time_of_unimol(rx, rng);

#ifdef DEBUG_RXNS
  SimulationStats* world = &p.stats;
  DUMP_CONDITION4(
      // calling rng for unimolecular
      m.dump(p, "Assigned unimolecular time (prev rng):", "", p.stats.get_current_iteration(), unimol_time_from_now);
  );
#endif

  // ignore if the is the next change if rxn probability closer than the scheduled time
  double update_time = rx->get_next_time_of_rxn_rate_update();
  if (current_time + unimol_time_from_now > update_time) {
    m.set_flag(MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE);
    unimol_time_from_now = update_time - current_time;
  }

  return unimol_time_from_now;
}

/*************************************************************************
which_unimolecular:
  In: the reaction we're testing
  Out: int containing which unimolecular reaction occurs (one must occur)
*************************************************************************/
static BNG::rxn_class_pathway_index_t which_unimolecular(const Molecule& m, BNG::RxnClass *rxn_class, rng_state& rng) {
  assert(rxn_class != nullptr);
  if (rxn_class->get_num_reactions() == 1) {
    return 0;
  }

  // TODO_COMPARTMENTS check compartments & pass to which_unimolecular only
  // those that match compartments

  double match = rng_dbl(&rng);
  match = match * rxn_class->get_max_fixed_p();
  return rxn_class->get_pathway_index_for_probability(match, 1);
}

} // namespace RxUtil

} // namespace MCell

#endif // SRC4_RXN_UTILS_INC_
