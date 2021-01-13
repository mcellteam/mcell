/*
 * elementary_molecule.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <bng/elem_mol_type.h>
#include <iostream>
#include <sstream>

#include "bng/bng_engine.h"

using namespace std;

namespace BNG {

// ------------- ComponentType -------------
std::string ComponentType::to_str(const BNGData& bng_data) const {
  stringstream out;
  out << name;

  for (state_id_t state_id: allowed_state_ids) {
    out << "~" << bng_data.get_state_name(state_id);
  }
  return out.str();
}


void ComponentType::dump(const BNGData& bng_data) const {
  cout << to_str(bng_data);
}


float_t get_default_space_step(const BNGConfig& config, const float_t D) {
  return sqrt_f(4.0 * 1.0e8 * D * config.time_unit) * config.rcp_length_unit;
}

// ------------- MoleculeType -------------
// based on assemble_mol_species
/**************************************************************************
assemble_mol_species:
   Helper function to assemble a molecule species from its component pieces.

   NOTE: A couple of comments regarding the unit conversions below:
   Internally, mcell works with with the per species length
   normalization factor

      new_spec->space_step = sqrt(4*D*t), D = diffusion constant (1)

   If the user supplies a CUSTOM_SPACE_STEP or SPACE_STEP then
   it is assumed to correspond to the average diffusion step and
   is hence equivalent to lr_bar in 2 or 3 dimensions for surface and
   volume molecules, respectively:

   lr_bar_2D = sqrt(pi*D*t)       (2)
   lr_bar_3D = 2*sqrt(4*D*t/pi)   (3)

   Hence, given a CUSTOM_SPACE_STEP/SPACE_STEP we need to
   solve eqs (2) and (3) for t and obtain new_spec->space_step
   via equation (1)

   2D:
    lr_bar_2D = sqrt(pi*D*t) => t = (lr_bar_2D^2)/(pi*D)

   3D:
    lr_bar_3D = 2*sqrt(4*D*t/pi) => t = pi*(lr_bar_3D^2)/(16*D)

   The remaining coefficients are:

    - 1.0e8 : needed to convert D from cm^2/s to um^2/s
    - global_time_unit, length_unit, r_length_unit: mcell
      internal time/length conversions.

In: state: the simulation state
    sym_ptr:   symbol for the species
    D:     diffusion constant
    is_2d: 1 if the species is a 2D molecule, 0 if 3D
    custom_time_step: time_step for the molecule (<0.0 for a custom space
                      step, >0.0 for custom timestep, 0.0 for default
                      timestep)
    target_only: 1 if the molecule cannot initiate reactions
Out: the species, or NULL if an error occurred
**************************************************************************/

void get_space_and_time_step(
    const BNGConfig& config,
    const bool is_surf, const float_t D,
    const float_t custom_time_step, const float_t custom_space_step,
    float_t& time_step, float_t& space_step) {

  assert(D != FLT_INVALID);
   // Immobile (boring)
   if (!distinguishable_f(D, 0, EPS)) {
     space_step = 0.0;
     time_step = 1.0;
   }
   // Custom timestep or spacestep
   else if (custom_space_step > 0) {
     float_t lr_bar = custom_space_step;
     if (is_surf) {
       time_step =
           lr_bar * lr_bar / (BNG_PI * 1.0e8 * D * config.time_unit);
     }
     else {
       time_step =
           lr_bar * lr_bar * BNG_PI / (16.0 * 1.0e8 * D * config.time_unit);
     }
     space_step =
         sqrt_f(4.0 * 1.0e8 * D * time_step * config.time_unit) * config.rcp_length_unit;
   }
   else if (custom_time_step > 0) {
     space_step =
         sqrt_f(4.0 * 1.0e8 * D * custom_time_step) * config.rcp_length_unit;
     time_step = custom_time_step / config.time_unit;
   }
   // Global timestep (this is the typical case)
   else /*if (!distinguishable(state->space_step, 0, EPS_C))*/ {
     space_step = get_default_space_step(config, D);
     time_step = DEFAULT_TIME_STEP;
   }
   /*// Global spacestep - not supported yet
   else {
     double space_step = state->space_step * state->length_unit;
     if (species->is_2d) {
       new_spec->time_step =
           space_step * space_step /
           (MY_PI * 1.0e8 * new_spec->D * global_time_unit);
     }
     else {
       new_spec->time_step =
           space_step * space_step * MY_PI /
           (16.0 * 1.0e8 * new_spec->D * global_time_unit);
     }
     new_spec->space_step = sqrt(4.0 * 1.0e8 * new_spec->D *
                                 new_spec->time_step * global_time_unit) *
                                 state->r_length_unit;
   }*/
}


// time_unit and length_unit are coming from the simulation configuration
void ElemMolType::compute_space_and_time_step(const BNGConfig& config) {
  get_space_and_time_step(
      config, 
      has_flag_no_finalized_check(SPECIES_CPLX_MOL_FLAG_SURF), D,
      custom_time_step, custom_space_step,
      time_step, space_step
  );
}


std::string ElemMolType::to_str(const BNGData& bng_data) const {
  stringstream out;

  out << name << "(";

  for (size_t i = 0; i < component_type_ids.size(); i++) {

    const ComponentType& ct = bng_data.get_component_type(component_type_ids[i]);
    out << ct.to_str(bng_data);

    if (i != component_type_ids.size() - 1) {
      out << ",";
    }
  }

  out << ")";
  return out.str();
}


void ElemMolType::dump(const BNGData& bng_data) const {

  cout << to_str(bng_data);
  cout <<
      " D=" << D <<
      ", custom_time_step=" << custom_time_step <<
      ", custom_space_step=" << custom_space_step;
}

} /* namespace BNG */
