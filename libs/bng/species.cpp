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
// flags according to reactions in the system
void Species::update_rxn_flags(const SpeciesContainer& all_species, RxnContainer& all_rxns) {
  if (rxn_flags_were_updated) {
    return;
  }

  assert(id != SPECIES_ID_INVALID);
#ifndef NDEBUG
  all_species.get(id); // must not fail
#endif
  species_id_t all_molecules_species_id = all_species.get_all_molecules_species_id();
  species_id_t all_volume_molecules_species_id = all_species.get_all_volume_molecules_species_id();
  species_id_t all_surface_molecules_species_id = all_species.get_all_surface_molecules_species_id();

  bool all_vol_mols_can_react_with_surface = all_rxns.all_vol_mols_can_react_with_surface;
  bool all_surf_mols_can_react_with_surface = all_rxns.all_surf_mols_can_react_with_surface;

  set_flag(BNG::SPECIES_FLAG_CAN_VOLVOL, false);
  set_flag(BNG::SPECIES_FLAG_CAN_VOLSURF, false);
  set_flag(BNG::SPECIES_FLAG_CAN_VOLWALL, false);
  set_flag(BNG::SPECIES_FLAG_CAN_SURFSURF, false);
  set_flag(BNG::SPECIES_FLAG_CAN_REGION_BORDER, false);

  if (is_vol() && all_vol_mols_can_react_with_surface) {
    set_flag(BNG::SPECIES_FLAG_CAN_VOLWALL);
  }

  if (is_surf() && all_surf_mols_can_react_with_surface) {
    set_flag(BNG::SPECIES_FLAG_CAN_REGION_BORDER);
  }

  // get reactions, this also creates all reaction classes for the species that we currently have
  BNG::SpeciesRxnClassesMap* rxn_classes =
      all_rxns.get_bimol_rxns_for_reactant(id);
  if (rxn_classes == nullptr) {
    rxn_flags_were_updated = true;
    return;
  }

  // go through all applicable reactants
  for (auto it: *rxn_classes) {
    const BNG::RxnClass* rxn_class = it.second;
    assert(rxn_class->is_bimol());

    species_id_t second_species_id;
    if (rxn_class->specific_reactants[0] != id) {
      second_species_id = rxn_class->specific_reactants[0];
    }
    else {
      second_species_id = rxn_class->specific_reactants[1];
    }
    const BNG::Species& sp2 = all_species.get(second_species_id);

    // we can use is_vol/is_surf for ALL_VOLUME_MOLECULES and ALL_SURFACE_MOLECULES
    if (is_vol()) {
      if (sp2.is_vol() || sp2.id == all_molecules_species_id) {
        set_flag(BNG::SPECIES_FLAG_CAN_VOLVOL);
      }
      if (sp2.is_surf() || sp2.id == all_molecules_species_id) {
        set_flag(BNG::SPECIES_FLAG_CAN_VOLSURF);
      }
      if (sp2.is_reactive_surface()) {
        set_flag(BNG::SPECIES_FLAG_CAN_VOLWALL);
      }
    }

    if ((is_surf() || id == all_molecules_species_id)) {
      if (sp2.is_surf() || sp2.id == all_molecules_species_id) {
        set_flag(BNG::SPECIES_FLAG_CAN_SURFSURF);
      }

      if (sp2.is_reactive_surface()) {
        set_flag(BNG::SPECIES_FLAG_CAN_REGION_BORDER);
      }
    }
  }

  rxn_flags_were_updated = true;
}


// based on assemble_mol_species
// time_unit and length_unit are coming from the simulation configuration
void Species::update_space_and_time_step(const BNGConfig& config) {
  // Immobile (boring)
  if (!distinguishable_f(D, 0, EPS)) {
    space_step = 0.0;
    time_step = 1.0;
  }
  // Custom timestep or spacestep
  // not supported yet
  /*else if (new_spec->time_step != 0.0) {
    // Hack--negative value means custom space step
    if (new_spec->time_step < 0) {
      double lr_bar = -new_spec->time_step;
      if (species->is_2d) {
        new_spec->time_step =
            lr_bar * lr_bar / (MY_PI * 1.0e8 * new_spec->D * global_time_unit);
      } else {
        new_spec->time_step =
            lr_bar * lr_bar * MY_PI /
            (16.0 * 1.0e8 * new_spec->D * global_time_unit);
      }
      new_spec->space_step =
          sqrt(4.0 * 1.0e8 * new_spec->D * new_spec->time_step *
               global_time_unit) *
          state->r_length_unit;
    }
    else {
      new_spec->space_step =
          sqrt(4.0 * 1.0e8 * new_spec->D * new_spec->time_step) *
          state->r_length_unit;
      new_spec->time_step /= global_time_unit;
    }
  }*/
  // Global timestep (this is the typical case)
  else /*if (!distinguishable(state->space_step, 0, EPS_C))*/ {
    space_step = sqrt_f(4.0 * 1.0e8 * D * config.time_unit) / config.length_unit;
    time_step = 1.0;
  }
  /*// Global spacestep
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
    D = data.get_molecule_type(mol_instances[0].mol_type_id).D;
    assert(D != FLT_INVALID);
  }
  else {
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
