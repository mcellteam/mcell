/*
 * RxnClass.cpp
 *
 *  Created on: Mar 25, 2020
 *      Author: ahusar
 */

#include <iostream>
#include <cmath>

#include "bng/rxn_class.h"
#include "bng/species_container.h"
#include "debug_config.h"

using namespace std;

namespace BNG {

// might need to be different for NFsim
// not sure if this belongs here
float_t RxnClass::get_reactant_space_step(const uint reactant_index) const {
  assert(reactant_index < specific_reactants.size());

  const Species& s = all_species.get(specific_reactants[reactant_index]);
  return s.space_step;
}


float_t RxnClass::get_reactant_time_step(const uint reactant_index) const {
  assert(reactant_index < specific_reactants.size());

  const Species& s = all_species.get(specific_reactants[reactant_index]);
  return s.time_step;
}


float_t RxnClass::get_reactant_diffusion(const uint reactant_index) const {
  assert(reactant_index < specific_reactants.size());

  const Species& s = all_species.get(specific_reactants[reactant_index]);
  return s.D;
}


// function for computing the probability factor (pb_factor) used to
// convert reaction rate constants into probabilities
float_t RxnClass::compute_pb_factor() const {

#ifndef NDEBUG
  assert(get_num_reactions() >= 1);
  // checking that all reactions in the same rxn class have the same orientation,
  // this is used later
  if (is_bimol()) {
    orientation_t orient0 = rxn_rules[0]->reactants[0].get_orientation();
    orientation_t orient1 = rxn_rules[0]->reactants[1].get_orientation();
    for (uint i = 0; i < get_num_reactions(); i++) {
      assert(orient0 == rxn_rules[i]->reactants[0].get_orientation());
      assert(orient1 == rxn_rules[i]->reactants[1].get_orientation());
    }
  }
  else {
    // orientation does not make much sense for unimol rxns, but let's check it as well
    orientation_t orient0 = rxn_rules[0]->reactants[0].get_orientation();
    for (uint i = 0; i < get_num_reactions(); i++) {
      assert(orient0 == rxn_rules[i]->reactants[0].get_orientation());
    }
  }
#endif

  /* determine the number of volume and surface reactants as well
   * as the number of surfaces */
  uint num_vol_reactants = 0;
  uint num_surf_reactants = 0;
  uint num_surfaces = 0;

  // mcell3 divides the radius after reaction initialization,
  // we already have the value that was divided
  float_t rx_radius_3d_mul_length_unit = bng_config.rx_radius_3d * bng_config.length_unit;

  small_vector<const Species*> reactant_species;
  for (uint n_reactant = 0; n_reactant < specific_reactants.size(); n_reactant++) {
    const Species& s = all_species.get(specific_reactants[n_reactant]);
    reactant_species.push_back(&s);
    if (s.is_surf()) {
      num_surf_reactants++;
    }
    else if (s.is_vol()) {
      num_vol_reactants++;
    }
    else if (s.is_reactive_surface()) {
      num_surfaces++;
    }
  }

  /* probability for this reaction */
  float_t pb_factor = 0.0;

  /* determine reaction probability by proper conversion of the reaction rate constant */
  assert(specific_reactants.size() == 1 || specific_reactants.size() == 2);
  if (specific_reactants.size() == 1) {
    // unimolecular
    pb_factor = bng_config.time_unit;
  }
  else if (num_surf_reactants >= 1 || num_surfaces == 1) {

    if (num_surf_reactants == 2) {
      /* this is a reaction between two surface molecules */
      if (reactant_species[0]->cant_initiate() && reactant_species[1]->cant_initiate()) {
        warns() <<
            "There is a surface reaction between " << reactant_species[0]->name << " and " <<
            reactant_species[1]->name << ", but neither of them can initiate this reaction " <<
            "(they are marked as target_only)\n";
      }

      if (reactant_species[0]->cant_initiate() || reactant_species[1]->cant_initiate()) {
        pb_factor = bng_config.time_unit * bng_config.grid_density / 3; /* 3 neighbors */
      }
      else {
        pb_factor = bng_config.time_unit * bng_config.grid_density / 6; /* 2 molecules, 3 neighbors each */
      }
    }
    else if (
        (reactant_species[0]->is_reactive_surface() && reactant_species[1]->is_surf()) ||
        (reactant_species[1]->is_reactive_surface() && reactant_species[0]->is_surf())
    ) {

      /* This is actually a unimolecular reaction in disguise! */
      pb_factor = bng_config.time_unit;
    }
    else if (
        (num_vol_reactants == 1 && num_surfaces == 1) ||
        (num_vol_reactants == 1 && num_surf_reactants == 1)
    ) {
      /* this is a reaction between "vol_mol" and "surf_mol" */
      /* or reaction between "vol_mol" and SURFACE           */
      // NOTE: we do not care about cant_initiate here, is this correct?

      float_t D_tot = 0.0;
      float_t t_step = 0.0;
      if (reactant_species[0]->is_vol()) {
        D_tot = get_reactant_diffusion(0);
        t_step = get_reactant_time_step(0) * bng_config.time_unit;
      }
      else if (reactant_species[1]->is_vol()) {
        D_tot = get_reactant_diffusion(1);
        t_step = get_reactant_time_step(1) * bng_config.time_unit;
      } else {
        /* Should never happen. */
        assert(false);
        D_tot = 1.0;
        t_step = 1.0;
      }

      if (D_tot <= 0.0) {
        pb_factor = 0; /* Reaction can't happen (not true for vol--surface rxns) */
      }
      else {
        pb_factor = 1.0e11 * bng_config.grid_density / (2.0 * BNG_N_AV) * sqrt(BNG_PI * t_step / D_tot);
      }

      // reactant_species are general, the do not have orientation
      // assuming that all reactions in the same rxn class have the same orientation
      orientation_t orient0 = rxn_rules[0]->reactants[0].get_orientation();
      orientation_t orient1 = rxn_rules[0]->reactants[1].get_orientation();

      // double pb factor if both use orientation (both orientations are not zero)
      // TODO: search elsewhere in the code for explanation, the first condition
      // seems superfluous because already the multiplications tells us that neither of them is 0
      // (assuming the allowed values are -1, 0, 1)
      assert(orient0 == ORIENTATION_UP || orient0 == ORIENTATION_NONE || orient0 == ORIENTATION_DOWN);
      assert(orient1 == ORIENTATION_UP || orient1 == ORIENTATION_NONE || orient1 == ORIENTATION_DOWN);

      if ( ((orient0 + orient1) * (orient0 - orient1) == 0) && (orient0 * orient1 != 0) ) {
        pb_factor *= 2.0;
      }
    }
  }
  else if (num_vol_reactants == 2) {
    /* This is the reaction between two "vol_mols" */

    float_t eff_vel_a = get_reactant_space_step(0) / get_reactant_time_step(0);
    float_t eff_vel_b = get_reactant_space_step(1) / get_reactant_time_step(1);
    float_t eff_vel;

    if (reactant_species[0]->cant_initiate() && reactant_species[1]->cant_initiate()) {
      warns() <<
          "There is a volume reaction between " << reactant_species[0]->name << " and " <<
          reactant_species[1]->name << ", but neither of them can initiate this reaction " <<
          "(they are marked as target_only)\n";
    }
    else if (reactant_species[0]->cant_initiate()) {
      eff_vel_a = 0;
    }
    else if (reactant_species[1]->cant_initiate()) {
      eff_vel_b = 0;
    }

    if (eff_vel_a + eff_vel_b > 0) {
      eff_vel = (eff_vel_a + eff_vel_b) * bng_config.length_unit / bng_config.time_unit; /* Units=um/sec */
      pb_factor = 1.0 / (2.0 * sqrt(BNG_PI) * rx_radius_3d_mul_length_unit * rx_radius_3d_mul_length_unit * eff_vel);
      pb_factor *= 1.0e15 / BNG_N_AV; /* Convert L/mol.s to um^3/number.s */
    }
    else {
      pb_factor = 0.0; /* No rxn possible */
    }
  }
  else {
    assert(false);
  }

  return pb_factor;
}


// based on mcell3's implementation init_reactions
void RxnClass::compute_initial_rxn_rates() {

#ifdef ORDER_RXNS_IN_RXN_CLASS_BY_NAME
  sort(rxn_rules.begin(), rxn_rules.end(),
      [](const RxnRule* a, const RxnRule* b) -> bool {
          return a->name < b->name;
      }
  );
#endif

  // alphabetize? - order of reactants

  // TODO LATER: check_reaction_for_duplicate_pathways

  cum_probs.resize(rxn_rules.size());

  // initialize rates
  for (uint i = 0; i < rxn_rules.size(); i++) {
    cum_probs[i] = rxn_rules[i]->get_rate_constant();
  }

  float_t pb_factor = compute_pb_factor();

  // scale_rxn_probabilities
  // TODO LATER: info and warning printouts
  for (uint i = 0; i < rxn_rules.size(); i++) {
    if (fabs(cum_probs[i] - FLT_GIGANTIC) < EPS) {
      // special surface reactions are not scaled because their pb_factor is 0
      continue;
    }
    float_t rate = pb_factor * cum_probs[i];
    cum_probs[i] = rate;
  }

  // init_reactions - compute cumulative properties
  for (uint i = 1; i < rxn_rules.size(); i++) {
    cum_probs[i] += cum_probs[i - 1];
  }

  // NOTE: when can be these probabilities different?
  if (!rxn_rules.empty()) {
    max_fixed_p = cum_probs.back();
    min_noreaction_p = max_fixed_p;
  }
  else {
    max_fixed_p = 1.0;
    min_noreaction_p = 1.0;
  }

  // set class' rxn type
  type = RxnType::Invalid;
  for (uint i = 0; i < rxn_rules.size(); i++) {
    assert(rxn_rules[i]->type != RxnType::Invalid && "Type for individual rxns must be set");

    if (type == RxnType::Invalid) {
      type = rxn_rules[i]->type;
    }
    else {
      // type must be the same as before
      assert(type == rxn_rules[i]->type);
    }
  }
}


void RxnClass::update_variable_rxn_rates(const float_t current_time) {
  bool any_changed = false;
  vector<bool> specific_rxn_changed;

  for (RxnRule* rxn: rxn_rules) {
    bool current_changed = rxn->update_variable_rxn_rate(current_time, this);
    specific_rxn_changed.push_back(current_changed);
    any_changed |= current_changed;
  }
  if (any_changed) {
    compute_initial_rxn_rates();
  }

  // report
  for (size_t i = 0; i < specific_rxn_changed.size(); i++) {
    if (specific_rxn_changed[i]) {
      float_t prob = (i == 0) ? cum_probs[0] : cum_probs[i] - cum_probs[i - 1];
      notifys() <<
          "Probability " << prob << " set for " << rxn_rules[i]->to_str() <<
          " at time " << current_time << ".\n";
    }
  }
}


void RxnClass::dump_array(const vector<RxnClass>& vec) {
  cout << "Reaction class array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump("  ");
  }
}


void RxnClass::dump(const std::string ind) const {
  assert(specific_reactants.size() == 1 || specific_reactants.size() == 2);
  cout << ind << "species_matching_reactants: " <<
      all_species.get(specific_reactants[0]).name << " (" << specific_reactants[0] << ")";

  if (specific_reactants.size() == 2) {
    cout << " + " << all_species.get(specific_reactants[1]).name << " (" << specific_reactants[1] << ")\n";
  }
  else {
    cout << "\n";
  }

  if (rxn_rules.empty()) {
    return;
  }

  cout << ind << "max_fixed_p: \t\t" << max_fixed_p << " [float_t] \t\t\n";
  cout << ind << "min_noreaction_p: \t\t" << min_noreaction_p << " [float_t] \t\t\n";
  cout << ind << "cum_probs: ";
  for (float_t p: cum_probs) {
    cout << p << ", ";
  }
  cout << "\n";

  for (const RxnRule* rxn: rxn_rules) {
    rxn->dump(false, ind);
    cout << "\n";
  }
}

} /* namespace BNG */
