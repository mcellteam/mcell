/*
 * RxnClass.cpp
 *
 *  Created on: Mar 25, 2020
 *      Author: ahusar
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "bng/rxn_class.h"
#include "bng/rxn_container.h"
#include "bng/species_container.h"
#include "debug_config.h"

using namespace std;

namespace BNG {

RxnClass::~RxnClass() {
  for (rxn_rule_id_t id: rxn_rule_ids) {
    RxnRule* rxn = all_rxns.get(id);
    rxn->remove_rxn_class_where_used(this);
  }
}


// might need to be different for NFsim
// not sure if this belongs here
float_t RxnClass::get_reactant_space_step(const uint reactant_index) const {
  assert(reactant_index < reactant_ids.size());

  const Species& s = all_species.get(reactant_ids[reactant_index]);
  return s.space_step;
}


float_t RxnClass::get_reactant_time_step(const uint reactant_index) const {
  assert(reactant_index < reactant_ids.size());

  const Species& s = all_species.get(reactant_ids[reactant_index]);
  return s.time_step;
}


float_t RxnClass::get_reactant_diffusion(const uint reactant_index) const {
  assert(reactant_index < reactant_ids.size());

  const Species& s = all_species.get(reactant_ids[reactant_index]);
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
  if (!pathways_and_rates_initialized) {
    // when a rxn class has only 1 rxn rule, the total prob may not be queried before
    init_rxn_pathways_and_rates();
  }
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


species_id_t RxnClass::get_reactive_surface_reactant_species_id() const {
  assert(is_bimol() && "Reactive surface cannot be unimol");
  if (all_species.get(reactant_ids[0]).is_reactive_surface()) {
    return reactant_ids[0];
  }
  else if (all_species.get(reactant_ids[1]).is_reactive_surface()) {
    return reactant_ids[1];
  }
  else {
    assert(false);
    return SPECIES_ID_INVALID;
  }
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


void RxnClass::define_rxn_pathway_using_mapping(const rxn_class_pathway_index_t pathway_index) {
  release_assert(pathways_and_rates_initialized);

  vector<species_id_t> reactant_species;
  for (auto& r: reactant_ids) {
    reactant_species.push_back(r);
  }

  const RxnRule* rxn = all_rxns.get(pathways[pathway_index].rxn_rule_id);
  rxn->define_rxn_pathway_using_mapping(all_species, bng_config, reactant_species, pathways[pathway_index]);
}


// based on MCell3's binary_search_double
rxn_class_pathway_index_t RxnClass::get_pathway_index_for_probability(
    const float_t prob, const float_t local_prob_factor) {
  if (!pathways_and_rates_initialized) {
    // when a rxn class has only 1 rxn rule, the total prob may not be queried before
    init_rxn_pathways_and_rates();
  }

  assert(!pathways.empty());
  int min_idx = 0;
  int max_idx = pathways.size() - 1;

  while (max_idx - min_idx > 1) {
    int mid_idx = (max_idx + min_idx) / 2;
    if (prob > (pathways[mid_idx].cum_prob * local_prob_factor)) {
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
  for (uint n_reactant = 0; n_reactant < reactant_ids.size(); n_reactant++) {
    const Species& s = all_species.get(reactant_ids[n_reactant]);
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
  assert(reactant_ids.size() == 1 || reactant_ids.size() == 2);
  if (reactant_ids.size() == 1 || is_intermembrane_surf_surf_rxn_class()) {
    // unimolecular or
    // experimental surf-surf on different objects
    pb_factor = bng_config.time_unit;
  }
  else if (num_surf_reactants >= 1 || num_surfaces == 1) {

    if (num_surf_reactants == 2) {
      /* this is a reaction between two surface molecules */
      if (reactant_species[0]->is_target_only() && reactant_species[1]->is_target_only()) {
        warns() <<
            "There is a surface reaction between " << reactant_species[0]->name << " and " <<
            reactant_species[1]->name << ", but neither of them can initiate this reaction " <<
            "(they are marked as target_only)\n";
      }

      if (reactant_species[0]->is_target_only() || reactant_species[1]->is_target_only()) {
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

      /* The value of pb_factor above is calculated for the case
          when surface_molecule can be hit from either side Otherwise the
          reaction_rate should be doubled. So we check whether both of the
          volume_molecules are in the same orientation class as
          surface_molecule.
      */
      assert(orient0 == ORIENTATION_UP || orient0 == ORIENTATION_NONE ||
          orient0 == ORIENTATION_DOWN || orient0 == ORIENTATION_DEPENDS_ON_SURF_COMP);
      assert(orient1 == ORIENTATION_UP || orient1 == ORIENTATION_NONE ||
          orient1 == ORIENTATION_DOWN || orient1 == ORIENTATION_DEPENDS_ON_SURF_COMP);

      // original condition: ((orient0 + orient1) * (orient0 - orient1) == 0) && (orient0 * orient1 != 0)
      // the first condition is not required
      if (orient0 * orient1 != 0) {
        pb_factor *= 2.0;
      }
    }
  }
  else if (num_vol_reactants == 2) {
    /* This is the reaction between two "vol_mols" */

    float_t eff_vel_a = get_reactant_space_step(0) / get_reactant_time_step(0);
    float_t eff_vel_b = get_reactant_space_step(1) / get_reactant_time_step(1);
    float_t eff_vel;

    if (reactant_species[0]->is_target_only() && reactant_species[1]->is_target_only()) {
      warns() <<
          "There is a volume reaction between " << reactant_species[0]->name << " and " <<
          reactant_species[1]->name << ", but neither of them can initiate this reaction " <<
          "(they are marked as target_only)\n";
    }
    else if (reactant_species[0]->is_target_only()) {
      eff_vel_a = 0;
    }
    else if (reactant_species[1]->is_target_only()) {
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


// does not do pathways update
void RxnClass::add_rxn_rule_no_update(RxnRule* r) {

  if (rxn_rule_ids.empty()) {
    bimol_vol_rxn_flag = r->is_bimol_vol_rxn();
  }
  else {
    // this should not normally happen
    assert(bimol_vol_rxn_flag == r->is_bimol_vol_rxn());
  }

  if (rxn_rule_ids.empty()) {
    intermembrane_surf_surf_rxn_flag = r->is_intermembrane_surf_rxn();
  }
  else {
    release_assert(intermembrane_surf_surf_rxn_flag == r->is_intermembrane_surf_rxn() &&
        "Intermembrane and other reactions cannot be mixed in a single rxn class");
  }

  // check that the rule was not added already,
  // for now simple pointer comparison
  for (rxn_rule_id_t id: rxn_rule_ids) {
    if (r->id == id) {
      // reaction is already present
      return;
    }
  }

  rxn_rule_ids.push_back(r->id);

  // remember bidirectional mapping for rxn rate updates
  r->add_rxn_class_where_used(this);

  // and also set rxn class type
  const RxnRule* rxn = all_rxns.get(r->id);
  if (type == RxnType::Invalid) {
    type = rxn->type;
  }
  else {
    assert(type == rxn->type);
  }
  assert(type != RxnType::Invalid);
}


// based on mcell3's implementation init_reactions
// but added support for cases where one reaction rule can have multiple sets of products
void RxnClass::init_rxn_pathways_and_rates(const bool force_update) {
  if (!force_update && pathways_and_rates_initialized) {
    return;
  }

  assert(!reactant_ids.empty());

#ifdef ORDER_RXNS_IN_RXN_CLASS_BY_NAME
  sort(rxn_rules.begin(), rxn_rules.end(),
      [](const RxnRule* a, const RxnRule* b) -> bool {
          return a->name < b->name;
      }
  );
#endif

  pathways.clear();

  // 1) compute binding probability factor
  float_t pb_factor = compute_pb_factor();

  // 2) define pathways
  // sort rules by ID to make sure we get identical results all the time
  sort(rxn_rule_ids.begin(), rxn_rule_ids.end());
  for (rxn_rule_id_t id: rxn_rule_ids) {

    RxnRule* rxn = all_rxns.get(id);
    rxn->define_rxn_pathways_for_specific_reactants(
        all_species,
        bng_config,
        reactant_ids[0],
        (is_bimol() ? reactant_ids[1] : SPECIES_ID_INVALID),
        pb_factor,
        pathways
    );
  }
  assert(!pathways.empty());

/*#ifndef MCELL4_DO_NOT_SORT_PATHWAYS
  // and sort them, due to rxn class removals in cleanup events, in different runs can the
  // pathways be sorted differently, we must maintain the order
  // may provide different result than MCell3
  RxnPathwayComparator pathway_cmp(all_species);
  sort(pathways.begin(), pathways.end(), pathway_cmp);
#endif*/

#ifdef MCELL4_REVERSED_RXNS_IN_RXN_CLASS
  // reverse pathways
  RxnClassPathwayVector rev_pathways;
  for (int i = pathways.size() - 1; i >= 0; i--) {
    rev_pathways.push_back(pathways[i]);
  }
  pathways = rev_pathways;
#endif

  // 3) compute cumulative properties
  pathways[0].cum_prob = pathways[0].pathway_prob;
  for (uint i = 1; i < pathways.size(); i++) {
    pathways[i].cum_prob = pathways[i].pathway_prob + pathways[i-1].cum_prob;
  }

  // 4) set remaining values for this rxn class
  // NOTE: when can be the max_fixed_p and min_noreaction_p probabilities different?
  if (!pathways.empty()) {
    max_fixed_p = pathways.back().cum_prob;
  }
  else {
    max_fixed_p = 1.0;
  }

  if (bng_config.rxn_and_species_report) {
    append_to_report(bng_config.get_rxn_report_file_name(), to_str() + "\n\n");
  }

  // reactive surfaces have always maximum probability
  if (!rxn_rule_ids.empty() &&
      !all_rxns.get(rxn_rule_ids[0])->is_unimol() &&
      !all_rxns.get(rxn_rule_ids[0])->is_reactive_surface_rxn()) {

    if (max_fixed_p > 1.0) {
      stringstream ss;
      ss << "Warning: total probability of reaction is > 1 (" << max_fixed_p << ")";
      cout << ss.str() << ", for reactant(s) " << reactants_to_str() << ".\n";
      append_to_report(bng_config.get_warnings_report_file_name(), ss.str() + "\n" + to_str());
      bng_config.warnings.bimol_rxn_probability_over_1 = true;
    }
    else if (max_fixed_p > 0.5) {
      // print final report after simulation ended
      bng_config.warnings.bimol_rxn_probability_over_05_less_1 = true;
    }
  }

  if (is_unimol() && max_fixed_p > MAX_UNIMOL_RXN_PROBABILITY) {
    errs() << "Unimolecular reaction class probability is " + f_to_str(max_fixed_p) +
        " which is higher than the maximum allowed value " + f_to_str(MAX_UNIMOL_RXN_PROBABILITY) << ". " <<
        "Terminating simulation because this would cause simulation to slow down or stop completely and is probably not what was expected. " <<
        "Please check your reaction rates in reaction class:\n" << to_str() << "\n";
    exit(1);
  }

  pathways_and_rates_initialized = true;
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
  if (!pathways_and_rates_initialized || any_changed) {
    init_rxn_pathways_and_rates(true);
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
    if (bng_config.notifications.rxn_probability_changed) {
      notifys() <<
          "Probability " << prob << " set for " << all_rxns.get(changed_id)->to_str() <<
          " at time " << current_time << ".\n";
    }
  }
}


std::string RxnClass::reactants_to_str() const {
  stringstream ss;
  ss << all_species.get(reactant_ids[0]).name << " (" << reactant_ids[0] << ")";
  if (reactant_ids.size() == 2) {
    ss << " + " << all_species.get(reactant_ids[1]).name << " (" << reactant_ids[1] << ")";
  }
  return ss.str();
}


std::string RxnClass::to_str(const std::string ind) const {
  stringstream out;
  assert(reactant_ids.size() == 1 || reactant_ids.size() == 2);
  out << ind << "rxn class for reactants: \n    " << reactants_to_str() << "\n";

  if (!pathways_and_rates_initialized) {
    out << "  pathways were not initialized\n";
  }

  if (rxn_rule_ids.empty()) {
    out << "  no pathways\n";
    return out.str();
  }

  for (size_t i = 0; i < pathways.size(); i++) {
    out << i << ": ";
    const RxnClassPathway& pw = pathways[i];
    RxnRule* rxn = all_rxns.get(pw.rxn_rule_id);

    out << "products based on rule " << rxn->to_str(true, false) << "\n    ";
    for (size_t k = 0; k < pw.product_species_w_indices.size(); k++) {
      species_id_t sid = pw.product_species_w_indices[k].product_species_id;
      out << all_species.get(sid).to_str() << " (" << sid << ") ";
      if (k != pw.product_species_w_indices.size() - 1) {
        out << " + ";
      }
    }
    out << "\n";
  }

  out << ind << "cum_probs: ";
  for (const RxnClassPathway& pw: pathways) {
      out << pw.cum_prob << ", ";
  }

  out << ind << "max_fixed_p: " << max_fixed_p << "\n";

  return out.str();
}


void RxnClass::dump_array(const vector<RxnClass>& vec) {
  cout << "Reaction class array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump("  ");
  }
}


void RxnClass::dump(const std::string ind) const {
  cout << to_str(ind);
}


} /* namespace BNG */
