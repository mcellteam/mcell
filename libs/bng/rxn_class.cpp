/*
 * RxnClass.cpp
 *
 *  Created on: Mar 25, 2020
 *      Author: ahusar
 */

#include <iostream>
#include <cmath>

#include "rxn_class.h"

#include "species_container.h"

using namespace std;

namespace BNG {


// might need to be different for NFsim
// not sure if this belongs here
float_t RxnClass::get_reactant_space_step(const uint reactant_index) const {
  assert(reactant_index < reactants.size());

  const Species& s = all_species.get(reactants[reactant_index]);
  return s.space_step;
}

float_t RxnClass::get_reactant_time_step(const uint reactant_index) const {
  assert(reactant_index < reactants.size());

  const Species& s = all_species.get(reactants[reactant_index]);
  return s.time_step;
}

float_t RxnClass::get_reactant_diffusion(const uint reactant_index) const {
  assert(reactant_index < reactants.size());

  const Species& s = all_species.get(reactants[reactant_index]);
  return s.D;
}

/*************************************************************************
 *
 * function for computing the probability factor (pb_factor) used to
 * convert reaction rate constants into probabilities
 *
 * in: mcell state
 *     rx to compute pb_factor for
 *     maximum number of expected surface products for this reaction
 *
 * out: pb_factor
 *
 ************************************************************************/
float_t RxnClass::compute_pb_factor(const BNGConfig& bng_config) const {

  /* determine the number of volume and surface reactants as well
   * as the number of surfaces */
  uint num_vol_reactants = 0;
  uint num_surf_reactants = 0;
  uint num_surfaces = 0;

  // mcell3 divides the radius after reaction initialization,
  // we already have the value that was divided
  float_t rx_radius_3d_mul_length_unit = bng_config.rx_radius_3d * bng_config.length_unit;

  small_vector<const Species*> reactant_species;
  for (uint n_reactant = 0; n_reactant < reactants.size(); n_reactant++) {
    const Species& s = all_species.get(reactants[n_reactant]);
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
  assert(reactants.size() == 1 || reactants.size() == 2);
  if (reactants.size() == 1) {
    // unimolecular
    pb_factor = bng_config.time_unit;
  }
  else if (num_surf_reactants >= 1 || num_surfaces == 1) {

    if (num_surf_reactants == 2) {
      /* this is a reaction between two surface molecules */
      pb_factor = bng_config.time_unit * bng_config.grid_density / 6; /* 2 molecules, 3 neighbors each */
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

      float_t D_tot = 0.0;
      float_t t_step = 0.0;
      if (reactant_species[0]->is_surf()) {
        D_tot = get_reactant_diffusion(0);
        t_step = get_reactant_time_step(0) * bng_config.time_unit;
      }
      else if (reactant_species[1]->is_surf()) {
        D_tot = get_reactant_diffusion(1);
        t_step = get_reactant_time_step(1) * bng_config.time_unit;
      } else {
        /* Should never happen. */
        assert(false);
        D_tot = 1.0;
        t_step = 1.0;
      }

      if (D_tot <= 0.0) {
        pb_factor = 0; /* Reaction can't happen! */
      }
      else {
        pb_factor = 1.0e11 * bng_config.grid_density / (2.0 * BNG_N_AV) * sqrt(BNG_PI * t_step / D_tot);
      }

      /*
      // not sure what to do with this code from the orig implemetation
      if ((rx->geometries[0] + rx->geometries[1]) * (rx->geometries[0] - rx->geometries[1]) == 0 &&
          rx->geometries[0] * rx->geometries[1] != 0) {
        pb_factor *= 2.0;
      }
      */
    }
  }
  else if (num_vol_reactants == 2) {
    /* This is the reaction between two "vol_mols" */

    float_t eff_vel_a = get_reactant_space_step(0) / get_reactant_time_step(0);
    float_t eff_vel_b = get_reactant_space_step(1) / get_reactant_time_step(1);
    float_t eff_vel;

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
void RxnClass::update(const BNGConfig& bng_config) {

  // alphabetize?
  // also, we might need to sort the reactions somehow (later, when input is from Python)

  // TODO: check_reaction_for_duplicate_pathways

  cum_probs.resize(reactions.size());

  // initialize rates
  for (uint i = 0; i < reactions.size(); i++) {
    cum_probs[i] = reactions[i]->rate_constant;
  }

  float_t pb_factor = compute_pb_factor(bng_config);

  // scale_rxn_probabilities
  // TODO: info and warning printouts
  for (uint i = 0; i < reactions.size(); i++) {
    float_t rate = pb_factor * cum_probs[i];
    cum_probs[i] = rate;
  }

  // init_reactions - compute cumulative properties
  for (uint i = 1; i < reactions.size(); i++) {
    cum_probs[i] += cum_probs[i - 1];
  }

  // NOTE: when can be these probabilities different?
  if (!reactions.empty()) {
    max_fixed_p = cum_probs.back();
    min_noreaction_p = max_fixed_p;
  }
  else {
    max_fixed_p = 1.0;
    min_noreaction_p = 1.0;
  }

  // set class' rxn type
  type = RxnType::Invalid;
  for (uint i = 0; i < reactions.size(); i++) {
    assert(reactions[i]->type != RxnType::Invalid && "Type for individual rxns must be set");

    if (type == RxnType::Invalid) {
      type = reactions[i]->type;
    }
    else {
      // type must be the same as before
      assert(type == reactions[i]->type);
    }
  }
}


void RxnClass::dump_array(const BNGData& bng_data, const vector<RxnClass>& vec) {
  cout << "Reaction class array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump(bng_data, "  ");
  }
}


void RxnClass::dump(const BNGData& bng_data, const std::string ind) const {
  assert(reactants.size() == 1 || reactants.size() == 2);
  cout << ind <<
      all_species.get(reactants[0]).name << " (" << reactants[0] << ")";
  if (reactants.size() == 2) {
    cout << " + " << all_species.get(reactants[1]).name << " (" << reactants[1] << ")\n";
  }
  else {
    cout << "\n";
  }

  if (reactions.empty()) {
    return;
  }

  cout << ind << "max_fixed_p: \t\t" << max_fixed_p << " [float_t] \t\t\n";
  cout << ind << "min_noreaction_p: \t\t" << min_noreaction_p << " [float_t] \t\t\n";
  cout << ind << "cum_probs: ";
  for (float_t p: cum_probs) {
    cout << p << ", ";
  }
  cout << "\n";

  for (const RxnRule* rxn: reactions) {
    rxn->dump(bng_data, ind);
    cout << "\n";
  }
}

} /* namespace BNG */
