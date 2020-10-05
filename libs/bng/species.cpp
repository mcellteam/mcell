/******************************************************************************
 *
 * Copyright (C) 2019 by
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

#include <iostream>
#include <sstream>

#include "bng/species.h"
#include "bng/bng_defines.h"
#include "bng/bng_data.h"
#include "bng/species_container.h"
#include "bng/rxn_container.h"

#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

using namespace std;

namespace BNG {

// sets SPECIES_FLAG_CAN_VOLVOL, SPECIES_FLAG_CAN_VOLSURF, SPECIES_FLAG_CAN_VOLWALL,
// SPECIES_FLAG_CAN_SURFSURF, and/or SPECIES_FLAG_CAN_REGION_BORDER
// flags according to reactions in the system,
// also sets SPECIES_MOL_FLAG_CANT_INITIATE
void Species::update_rxn_and_custom_flags(
    const SpeciesContainer& all_species, RxnContainer& all_rxns,
    const BaseCustomFlagsAnalyzer* flags_analyzer) {

  assert(id != SPECIES_ID_INVALID);
#ifndef NDEBUG
  all_species.get(id); // must not fail
#endif

  bool all_vol_mols_can_react_with_surface = all_rxns.all_vol_mols_can_react_with_surface;
  bool all_surf_mols_can_react_with_surface = all_rxns.all_surf_mols_can_react_with_surface;

  set_flag(SPECIES_FLAG_CAN_VOLVOL, false);
  set_flag(SPECIES_FLAG_CAN_VOLSURF, false);
  set_flag(SPECIES_FLAG_CAN_VOLWALL, false);
  set_flag(SPECIES_FLAG_CAN_SURFSURF, false);
  set_flag(SPECIES_FLAG_CAN_REGION_BORDER, false);

  // set any custom flags (only SPECIES_FLAG_NEEDS_COUNTED_VOLUME is set currently)
  if (flags_analyzer != nullptr) {
    uint mask_to_clear = flags_analyzer->get_custom_species_flags_mask();
    clear_flags(mask_to_clear);
    uint mask_to_set = flags_analyzer->get_custom_species_flags_to_set(*this);
    add_flags(mask_to_set);
  }

  if (is_vol() && all_vol_mols_can_react_with_surface) {
    set_flag(SPECIES_FLAG_CAN_VOLWALL);
  }

  if (is_surf() && all_surf_mols_can_react_with_surface) {
    set_flag(SPECIES_FLAG_CAN_REGION_BORDER);
  }

  // go through all applicable reactions
  // it does not make sense to create reaction classes now because
  // we might not have all the reactant species and when reactant species do not exist,
  // there are no bimol reaction classes, so we must really go through all the reactions
  for (RxnRule* rxn: all_rxns.get_rxn_rules_vector()) {
    if (!rxn->species_can_be_reactant(id, all_species)) {
      continue;
    }

    if (rxn->is_unimol()) {
      set_flag(SPECIES_FLAG_HAS_UNIMOL_RXN);

      if (rxn->is_counted_in_volume_regions()) {
        set_flag(SPECIES_FLAG_NEEDS_COUNTED_VOLUME, true);
      }
      continue;
    }

    assert(rxn->is_bimol());

    if (rxn->is_counted_in_volume_regions()) {
      set_flag(SPECIES_FLAG_NEEDS_COUNTED_VOLUME);
    }

    if (rxn->is_bimol_vol_rxn()) {
      set_flag(SPECIES_FLAG_HAS_BIMOL_VOL_RXN);
    }

    // second reactant - we must not define new species for reactant here
    // because adding it to the species array might invalidate *this,
    // we must analyze the cplx instance
    int my_index = rxn->get_reactant_index(*this, all_species);
    assert(my_index == 0 || my_index == 1);
    int second_index = (my_index + 1) % 2;

    const CplxInstance& second_reactant_inst = rxn->reactants[second_index];

    // we can use is_vol/is_surf for ALL_VOLUME_MOLECULES and ALL_SURFACE_MOLECULES
    if (is_vol()) {
      if (second_reactant_inst.is_vol()) {
        set_flag(SPECIES_FLAG_CAN_VOLVOL);
      }
      if (second_reactant_inst.is_surf()) {
        set_flag(SPECIES_FLAG_CAN_VOLSURF);
      }
      if (second_reactant_inst.is_reactive_surface()) {
        set_flag(SPECIES_FLAG_CAN_VOLWALL);
      }
    }

    if (is_surf()) {
      if (second_reactant_inst.is_surf()) {
        set_flag(SPECIES_FLAG_CAN_SURFSURF);
      }

      if (second_reactant_inst.is_reactive_surface()) {
        set_flag(SPECIES_FLAG_CAN_REGION_BORDER);
      }
    }
  }

  // SPECIES_MOL_FLAG_CANT_INITIATE == target_only
  const BNGData& bng_data = all_species.get_bng_data();
  bool all_mols_are_cant_initiate = true;
  for (const MolInstance& mi: mol_instances) {
    all_mols_are_cant_initiate = all_mols_are_cant_initiate &&
        bng_data.get_molecule_type(mi.mol_type_id).has_flag(SPECIES_MOL_FLAG_CANT_INITIATE);
  }
  set_flag(SPECIES_MOL_FLAG_CANT_INITIATE, all_mols_are_cant_initiate);

  rxn_flags_were_updated = true;
}


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
// time_unit and length_unit are coming from the simulation configuration
void Species::update_space_and_time_step(const BNGConfig& config) {
  // Immobile (boring)
  if (!distinguishable_f(D, 0, EPS)) {
    space_step = 0.0;
    time_step = 1.0;
  }
  // Custom timestep or spacestep
  else if (custom_space_step > 0) {
    float_t lr_bar = custom_space_step;
    assert(!has_flag_no_finalized_check(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE));
    if (has_flag_no_finalized_check(SPECIES_CPLX_MOL_FLAG_SURF)) {
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
    space_step = sqrt_f(4.0 * 1.0e8 * D * config.time_unit) * config.rcp_length_unit;
    time_step = 1.0;
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


void Species::update_diffusion_constant(const BNGData& data, const BNGConfig& config) {
  if (mol_instances.size() == 1) {
    // nothing to compute if we have just one molecule instance
    const MolType& mt = data.get_molecule_type(mol_instances[0].mol_type_id);
    D = mt.D;
    assert(D != FLT_INVALID);
    custom_space_step = mt.custom_space_step;
    custom_time_step = mt.custom_time_step;
  }
  else {
    // check for custom steps that are not supported yet for complexes
    for (const MolInstance& mi: mol_instances) {
      const MolType& mt = data.get_molecule_type(mi.mol_type_id);
      release_assert(mt.custom_space_step == 0 && mt.custom_time_step == 0 &&
          "Custom time and space steps are not yet supported for complexes having more than one elementary molecule"
      );
    }

    // based on NFSim DerivedDiffusion::getDiffusionValue
    // For large complex aggregates (1 um) diffusing on 2D surface:
    //D_2D = KB*T*LOG((eta_PM*h/(Rc*(eta_EC+eta_CP)/2))-gamma)/(4*PI*eta_PM*h) (Saffman Delbruck)
    //  combining rule for this case is not elaborated here.

    // For small complexes diffusing on 2D surface:
    //D_2D = ~1/Rc  (ref?)
    //  combining rule for this case is the squareroot of the sum of the squares

    // For complexes diffusing in 3D:
    //D_3D = KB*T/(6*PI*eta_EC*Rs) (Einstein Stokes)
    //  combining rule for this case is the cuberoot of the sum of the cubes

    if (is_surf()) {
      // if complex contains any 2D subunits then the whole complex is considered to be a surface complex.
      // in this case combine only the 2D subunits to derive the 2D diffusion constant
      bool is_immobile = false;
      float_t acc = 0;
      for (const MolInstance& mi: mol_instances) {
        if (mi.is_surf()) {
          float_t mol_type_D = data.get_molecule_type(mi.mol_type_id).D;

          // if diffusion constant of any 2D member is zero (i.e. is immobile) then whole complex should be immobile
          if (mol_type_D == 0) {
            is_immobile = true;
            break;
          }

          acc += pow_f(mol_type_D, -2.0);
        }
      }
      if (is_immobile) {
        D = 0;
      }
      else {
        D = pow(acc, -0.5);
      }
    }
    else {
      // Only if there are no 2D subunits should the complex be considered to be a volume complex.
      // in this case combine all the 3D subunits to derive the 3D diffusion constant

      // 3D combining rule: compute cuberoot of the sum of the cubes
      bool is_immobile = false;
      float_t acc = 0;
      for (const MolInstance& mi: mol_instances) {
        assert(!mi.is_surf());
        float_t mol_type_D = data.get_molecule_type(mi.mol_type_id).D;

        //  if diffusion constant of any 3D member is zero (i.e. is immobile) then whole complex should be immobile
        if (mol_type_D == 0) {
          is_immobile = true;
          break;
        }

        acc += pow_f(mol_type_D, -3.0);
      }

      if (is_immobile) {
        D = 0;
      }
      else {
        D = pow_f(acc, -1.0/3.0); // NFSim uses constant -0.3333333333333
      }
    }
  }

  update_space_and_time_step(config);
}

// temporary, to make the diff const same in dump as in data model
// TODO: move this to a shared library, the same function is in datamodel_defines.h
static inline std::string f_to_string(const float_t val, const int n = 17)
{
  std::stringstream out;
  if (val == 0.0) {
    return "0";
  }
  else if (val <= 0.01 || val >= 100000) {
    out << std::scientific;
  }
  out.precision(n);
  out << val;
  return out.str();
}


void Species::dump(const BNGData& bng_data, const string ind) const {
  cout << ind << "species_id: \t\t" << id << " [uint16_t] \t\t/* Unique ID for this species */\n";
  cout << ind << "name: *\t\t" << name << " [string] \t\t/* Symbol table entry (name) */\n";
  cout << ind << "D: \t\t" << f_to_string(D) << " [float_t] \t\t/* Diffusion constant */\n";
  cout << ind << "space_step: \t\t" << space_step << " [float_t] \t\t/* Characteristic step length */\n";
  cout << ind << "time_step: \t\t" << time_step << " [float_t] \t\t/* Minimum (maximum?) sensible timestep */\n";
  cout << ind << "flags: \t\t" << BaseSpeciesCplxMolFlag::to_str() << "\n";
  cout << ind << "CplxInstance:\n";
  CplxInstance::dump(true, ind + "  ");
  cout << "\n";
}


template <typename T>
vector<size_t> sort_indexes(const T &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1].name < v[i2].name;});

  return idx;
}


void Species::dump_array(const BNGData& bng_data, const SpeciesVector& vec, const bool sorted) {
  cout << "Species array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  if (sorted) {
    // dump sorted by name
    vector<size_t> sorted_indices = sort_indexes(vec);
    for (auto i: sorted_indices) {
      cout << i << ":\n";
      vec[i].dump(bng_data, "  ");
    }
  }
  else {
    for (size_t i = 0; i < vec.size(); i++) {
      cout << i << ":\n";
      vec[i].dump(bng_data, "  ");
    }
  }
}

} // namespace mcell
