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

  for (uint n_reactant = 0; n_reactant < reactants.size(); n_reactant++) {
    const Species& s = all_species.get(reactants[n_reactant]);
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
    pb_factor = bng_config.time_unit;
  }
  else if (num_surf_reactants >= 1 || num_surfaces == 1) {

    if ((num_surf_reactants == 2) && (num_vol_reactants == 0) &&
        (num_surfaces < 2)) {
      /* this is a reaction between two surface molecules */
      /* with an optional SURFACE                         */

      assert(false && "TODO");
    }
#if 0
      if (rx->players[0]->flags & rx->players[1]->flags & CANT_INITIATE)
        mcell_error("Reaction between %s and %s listed, but both are marked "
                    "TARGET_ONLY.",
                    rx->players[0]->sym->name, rx->players[1]->sym->name);
      else if ((rx->players[0]->flags | rx->players[1]->flags) &
               CANT_INITIATE) {
        pb_factor =
            time_unit * grid_density / 3; /* 3 neighbors */
      } else {
        pb_factor = time_unit * grid_density /
                    6; /* 2 molecules, 3 neighbors each */
      }
    } else if ((((rx->players[0]->flags & IS_SURFACE) != 0 &&
                 (rx->players[1]->flags & ON_GRID) != 0) ||
                ((rx->players[1]->flags & IS_SURFACE) != 0 &&
                 (rx->players[0]->flags & ON_GRID) != 0)) &&
               (rx->n_reactants == 2)) {
      /* This is actually a unimolecular reaction in disguise! */
      pb_factor = time_unit;
      if (max_num_surf_products > 0)
        *create_shared_walls_info_flag = 1;
    } else if (((rx->n_reactants == 2) && (num_vol_reactants == 1) &&
                (num_surfaces == 1)) ||
               ((rx->n_reactants == 2) && (num_vol_reactants == 1) &&
                (num_surf_reactants == 1)) ||
               ((rx->n_reactants == 3) && (num_vol_reactants == 1) &&
                (num_surf_reactants == 1) && (num_surfaces == 1))) {
      /* this is a reaction between "vol_mol" and "surf_mol" */
      /* with an optional SURFACE                            */
      /* or reaction between "vol_mol" and SURFACE           */
      if (max_num_surf_products > 0)
        *create_shared_walls_info_flag = 1;
      if (((rx->n_reactants == 2) && (num_vol_reactants == 1) &&
           (num_surfaces == 1))) {
        /* do not take into acccount SPECIAL reactions */
        if (rx->n_pathways > RX_SPECIAL) {
          rxn_flags->vol_wall_reaction_flag = 1;
        }
      } else {
        rxn_flags->vol_surf_reaction_flag = 1;
      }

      float_t D_tot = 0.0;
      float_t t_step = 0.0;
      if ((rx->players[0]->flags & NOT_FREE) == 0) {
        D_tot = rx->get_reactant_diffusion(rx, 0);
        t_step = rx->get_reactant_time_step(rx,0) * time_unit;
      } else if ((rx->players[1]->flags & NOT_FREE) == 0) {
        D_tot = rx->get_reactant_diffusion(rx, 1);
        t_step = rx->get_reactant_time_step(rx,1) * time_unit;
      } else {
        /* Should never happen. */
        D_tot = 1.0;
        t_step = 1.0;
      }

      if (D_tot <= 0.0)
        pb_factor = 0; /* Reaction can't happen! */
      else
        pb_factor = 1.0e11 * grid_density / (2.0 * N_AV) *
                    sqrt(MY_PI * t_step / D_tot);

      if ((rx->geometries[0] + rx->geometries[1]) *
                  (rx->geometries[0] - rx->geometries[1]) ==
              0 &&
          rx->geometries[0] * rx->geometries[1] != 0) {
        pb_factor *= 2.0;
      }
    } /* end else */
#endif
  }
  else if (num_vol_reactants == 2) {
    /* This is the reaction between two "vol_mols" */

    const Species& s1 = all_species.get(reactants[0]);
    const Species& s2 = all_species.get(reactants[1]);

    float_t eff_vel_a = get_reactant_space_step(0) / get_reactant_time_step(0);
    float_t eff_vel_b = get_reactant_space_step(1) / get_reactant_time_step(1);
    float_t eff_vel;

    if (eff_vel_a + eff_vel_b > 0) {
      eff_vel = (eff_vel_a + eff_vel_b) * bng_config.length_unit / bng_config.time_unit; /* Units=um/sec */
      pb_factor = 1.0 / (2.0 * sqrt(BNG_PI) * bng_config.rx_radius_3d * bng_config.rx_radius_3d * eff_vel);
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
    float_t rate = pb_factor * cum_probs[1];
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
}


/*
void RxnClass::add_rxn_rule(RxnRule* r) {
  reactions.push_back(r);
  update();
}
*/


void RxnClass::dump_array(const BNGData& bng_data, const vector<RxnClass>& vec) {
  cout << "Reaction class array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  for (size_t i = 0; i < vec.size(); i++) {
    cout << i << ":\n";
    vec[i].dump(bng_data, "  ");
  }
}

void RxnClass::dump(const BNGData& bng_data, const std::string ind) const {
  cout << ind << "max_fixed_p: \t\t" << max_fixed_p << " [float_t] \t\t\n";
  cout << ind << "min_noreaction_p: \t\t" << min_noreaction_p << " [float_t] \t\t\n";

  for (const RxnRule* rxn: reactions) {
    rxn->dump(bng_data, ind);
  }
}

} /* namespace BNG */
