/*
 * RxnClass.cpp
 *
 *  Created on: Mar 25, 2020
 *      Author: ahusar
 */

#include <iostream>
#include <cmath>

#include "bng/rxn_class.h"
#include "bng/rxn_container.h"
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


// this assert cannot be in the header file
void RxnClass::debug_check_bimol_vol_rxn_flag() const {
  assert(!rxn_rule_ids.empty() &&
      all_rxns.get(rxn_rule_ids[0])->is_bimol_vol_rxn() == bimol_vol_rxn_flag);
}


bool RxnClass::is_simple() const {
  for (rxn_rule_id_t id: rxn_rule_ids) {
    const RxnRule* rxn = all_rxns.get(id);
    if (!rxn->is_simple()) {
      return false;
    }
  }
  return true;
}


RxnRule* RxnClass::get_rxn_for_pathway(const rxn_class_pathway_index_t pathway_index) {
  assert(pathway_index >= 0 && pathway_index < (int)pathways.size());
  return all_rxns.get(pathways[pathway_index].rxn_rule_id);
}


orientation_t RxnClass::get_reactant_orientation(uint reactant_index) const {
  assert(!rxn_rule_ids.empty());
  // all reactants for this class are the same
  const RxnRule* first_rxn = all_rxns.get(rxn_rule_ids[0]);
  assert(reactant_index < first_rxn->reactants.size());
  return first_rxn->reactants[reactant_index].get_orientation();
}


void RxnClass::update_rxn_rates_if_needed(const float_t current_time) {
  // check if any of the reactions needs update
  for (rxn_rule_id_t id: rxn_rule_ids) {
    RxnRule* rxn = all_rxns.get(id);
    if (rxn->may_update_rxn_rate()) {
      update_variable_rxn_rates(current_time);
      break;
    }
  }
}


// this function expects that update_rxn_rates_if_needed was called
// already for the current time
float_t RxnClass::get_next_time_of_rxn_rate_update() const {
  float_t min = TIME_FOREVER;
  for (rxn_rule_id_t id: rxn_rule_ids) {
    const RxnRule* rxn = all_rxns.get(id);

    float_t t = rxn->get_next_time_of_rxn_rate_update();
    if (t < min) {
      min = t;
    }
  }
  return min;
}

// based on MCell3's binary_search_double
rxn_class_pathway_index_t RxnClass::get_pathway_index_for_probability(
    const float_t prob, const float_t local_prob_factor) const {
  assert(!pathways.empty());
  int min_idx = 0;
  int max_idx = pathways.size() - 1;

  while (max_idx - min_idx > 1) {
    int mid_idx = (max_idx + min_idx) / 2;
    if (prob > (pathways[min_idx].cum_prob * local_prob_factor)) {
      min_idx = mid_idx;
    }
    else {
      max_idx = mid_idx;
    }
  }

  if (prob > pathways[min_idx].cum_prob * local_prob_factor) {
    return max_idx;
  }
  else {
    return min_idx;
  }
}


// function for computing the probability factor (pb_factor) used to
// convert reaction rate constants into probabilities
float_t RxnClass::compute_pb_factor() const {

#ifndef NDEBUG
  assert(get_num_reactions() >= 1);
  // checking that all reactions in the same rxn class have the same orientation,
  // this is used later
  if (is_bimol()) {
    orientation_t orient0 = all_rxns.get(rxn_rule_ids[0])->reactants[0].get_orientation();
    orientation_t orient1 = all_rxns.get(rxn_rule_ids[0])->reactants[1].get_orientation();
    for (uint i = 0; i < get_num_reactions(); i++) {
      assert(orient0 == all_rxns.get(rxn_rule_ids[i])->reactants[0].get_orientation());
      assert(orient1 == all_rxns.get(rxn_rule_ids[i])->reactants[1].get_orientation());
    }
  }
  else {
    // orientation does not make much sense for unimol rxns, but let's check it as well
    orientation_t orient0 = all_rxns.get(rxn_rule_ids[0])->reactants[0].get_orientation();
    for (uint i = 0; i < get_num_reactions(); i++) {
      assert(orient0 == all_rxns.get(rxn_rule_ids[0])->reactants[0].get_orientation());
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
      orientation_t orient0 = get_reactant_orientation(0);
      orientation_t orient1 = get_reactant_orientation(1);

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
// but added support for cases where one reaction rule can have multiple sets of products
void RxnClass::update_rxn_pathways() {

  assert(!specific_reactants.empty());

#ifdef ORDER_RXNS_IN_RXN_CLASS_BY_NAME
  sort(rxn_rules.begin(), rxn_rules.end(),
      [](const RxnRule* a, const RxnRule* b) -> bool {
          return a->name < b->name;
      }
  );
#endif

  // TODO LATER: check_reaction_for_duplicate_pathways?

  pathways.clear();

  // 1) compute binding probability factor
  float_t pb_factor = compute_pb_factor();

  // 2) define pathways
  for (rxn_rule_id_t id: rxn_rule_ids) {
    const RxnRule* rxn = all_rxns.get(id);
    rxn->define_rxn_pathways_for_specific_reactants(
        all_species,
        bng_config,
        specific_reactants[0],
        (is_bimol() ? specific_reactants[1] : SPECIES_ID_INVALID),
        pb_factor,
        pathways
    );
  }
  assert(!pathways.empty());

  // 3) compute cumulative properties
  pathways[0].cum_prob = pathways[0].pathway_prob;
  for (uint i = 1; i < pathways.size(); i++) {
    pathways[i].cum_prob = pathways[i].pathway_prob + pathways[i-1].cum_prob;
  }

  // 4) set remaining values for this rxn class
  // NOTE: when can be the max_fixed_p and min_noreaction_p probabilities different?
  if (!pathways.empty()) {
    max_fixed_p = pathways.back().cum_prob;
    min_noreaction_p = max_fixed_p;
  }
  else {
    max_fixed_p = 1.0;
    min_noreaction_p = 1.0;
  }

  // set class' rxn type
  type = RxnType::Invalid;
  for (rxn_rule_id_t id: rxn_rule_ids) {
    const RxnRule* rxn = all_rxns.get(id);
    assert(rxn->type != RxnType::Invalid && "Type for individual rxns must be set");

    if (type == RxnType::Invalid) {
      type = rxn->type;
    }
    else {
      // type must be the same as before
      assert(type == rxn->type);
    }
  }
}


void RxnClass::update_variable_rxn_rates(const float_t current_time) {
  bool any_changed = false;
  vector<rxn_rule_id_t> changed_rxn_rules;

  for (rxn_rule_id_t id: rxn_rule_ids) {
    RxnRule* rxn = all_rxns.get(id);
    bool current_changed = rxn->update_variable_rxn_rate(current_time, this);
    if (current_changed) {
      changed_rxn_rules.push_back(id);
      any_changed = true;
    }
  }
  if (any_changed) {
    update_rxn_pathways();
  }

  // report
  for (rxn_rule_id_t changed_id: changed_rxn_rules) {
    // find the first corresponding pathway
    rxn_class_pathway_index_t first_pw_changed = PATHWAY_INDEX_INVALID;
    for (size_t pwi = 0; pwi < pathways.size(); pwi++) {
      if (pathways[pwi].rxn_rule_id == changed_id) {
        first_pw_changed = pwi;
        break;
      }
    }
    release_assert(first_pw_changed != PATHWAY_INDEX_INVALID);

    float_t prob = pathways[first_pw_changed].pathway_prob;
    notifys() <<
        "Probability " << prob << " set for " << all_rxns.get(changed_id)->to_str() <<
        " at time " << current_time << ".\n";
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

  if (rxn_rule_ids.empty()) {
    return;
  }

  cout << ind << "max_fixed_p: \t\t" << max_fixed_p << " [float_t] \t\t\n";
  cout << ind << "min_noreaction_p: \t\t" << min_noreaction_p << " [float_t] \t\t\n";
  cout << ind << "cum_probs: ";
  for (const RxnClassPathway& pw: pathways) {
    cout << pw.cum_prob << ", ";
  }
  cout << "\n";

  for (const RxnClassPathway& pw: pathways) {
    RxnRule* rxn = all_rxns.get(pw.rxn_rule_id);
    rxn->dump(false, ind);
    cout << ", product ids: ";
    for (species_id_t sid: pw.product_species) {
      cout << sid << ", ";
    }
    cout << "\n";
  }
}

} /* namespace BNG */
