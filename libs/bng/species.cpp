
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
    const SpeciesContainer& all_species, RxnContainer& all_rxns) {

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

      continue;
    }

    assert(rxn->is_bimol());

    if (rxn->is_bimol_vol_rxn()) {
      set_flag(SPECIES_FLAG_HAS_BIMOL_VOL_RXN);
    }

    if (rxn->is_intermembrane_surf_rxn()) {
      set_flag(SPECIES_FLAG_CAN_INTERMEMBRANE_SURFSURF);
    }

    // second reactant - we must not define new species for reactant here
    // because adding it to the species array might invalidate *this,
    // we must analyze the cplx
    std::vector<uint> indices;
    rxn->get_reactant_indices(id, all_species, indices);
    assert(!indices.empty() && (indices[0] == 0 || indices[0] == 1));
    int second_index = (indices[0] + 1) % 2;

    const Cplx& second_reactant = rxn->reactants[second_index];

    // we can use is_vol/is_surf for ALL_VOLUME_MOLECULES and ALL_SURFACE_MOLECULES
    if (is_vol()) {
      if (second_reactant.is_vol()) {
        set_flag(SPECIES_FLAG_CAN_VOLVOL);
      }
      if (second_reactant.is_surf()) {
        set_flag(SPECIES_FLAG_CAN_VOLSURF);
      }
      if (second_reactant.is_reactive_surface()) {
        set_flag(SPECIES_FLAG_CAN_VOLWALL);
      }
    }

    if (is_surf()) {
      if (second_reactant.is_surf()) {
        set_flag(SPECIES_FLAG_CAN_SURFSURF);
      }
      if (second_reactant.is_reactive_surface()) {
        set_flag(SPECIES_FLAG_CAN_REGION_BORDER);
        set_flag(SPECIES_FLAG_CAN_SURFWALL);
      }
    }
  }

  // SPECIES_MOL_FLAG_CANT_INITIATE == target_only
  const BNGData& bng_data = all_species.get_bng_data();
  bool all_mols_are_cant_initiate = true;
  for (const ElemMol& mi: elem_mols) {
    all_mols_are_cant_initiate = all_mols_are_cant_initiate &&
        bng_data.get_elem_mol_type(mi.elem_mol_type_id).has_flag(SPECIES_MOL_FLAG_TARGET_ONLY);
  }
  set_flag(SPECIES_MOL_FLAG_TARGET_ONLY, all_mols_are_cant_initiate);

  rxn_flags_were_updated = true;
}


void Species::compute_space_and_time_step(const BNGConfig& config) {

  // check if any of the used elementary molecules has custom time or space step
  // 0 means that it was not set
  float_t min_custom_time_step = FLT_GIGANTIC;
  float_t min_custom_space_step = FLT_GIGANTIC;
  bool has_mol_wo_custom_time_step = false;
  bool has_mol_wo_custom_space_step = false;

  for (uint i = 0; i < elem_mols.size(); i++) {
    const ElemMolType& emt = bng_data->get_elem_mol_type(elem_mols[i].elem_mol_type_id);

    if (emt.custom_time_step != 0) {
      if (emt.custom_time_step < min_custom_time_step) {
        min_custom_time_step = emt.custom_time_step;
      }
    }
    else {
      has_mol_wo_custom_time_step = true;
    }

    if (emt.custom_space_step != 0) {
      if (emt.custom_space_step < min_custom_space_step) {
        min_custom_space_step = emt.custom_space_step;
      }
    }
    else {
      has_mol_wo_custom_space_step = true;
    }
  }
  // if custom time step was either not set or it is > 1.0 but at least one of the
  // molecules does not have it set, so we must keep the default
  if (min_custom_time_step == FLT_GIGANTIC ||
      (has_mol_wo_custom_time_step && min_custom_time_step > DEFAULT_TIME_STEP)) {
    min_custom_time_step = 0;
  }

  float_t default_space_step = get_default_space_step(config, D);
  if (min_custom_space_step == FLT_GIGANTIC ||
      (has_mol_wo_custom_space_step && min_custom_space_step > default_space_step)) {
    min_custom_space_step = 0;
  }

  if (min_custom_time_step != 0 || min_custom_space_step != 0) {
    // custom time or space step is used either to:
    // 1) speedup simulation if the default probability of reactions is low enough and
    //    we don't care that the molecule diffuses further than usual or
    // 2) lower the probability of a reaction per collision so that it is lower than 1
    //
    // the goal when computing space and time step for a complex with a custom time or space step
    // is to lower the resulting probability factor so that when reaction of A + X has a certain
    // probability (pb_factor), then the probability of A.A + X is not higher
    // this is used only when the user specified custom time or space step, when not given,
    // we don't care
    //
    // pb_factor is proportional to sqrt(time_step)
    //
    // simplified equations (based on RxnClass::compute_pb_factor):
    //
    // vol:
    // pb_factor = 1/(C1 * (C2 + space_step/time_step))
    // space_step/time_step = C3/sqrt(custom_time_step)
    // -> pb_factor = C4 * sqrt(custom_time_step)  (setting C2 to 0)
    //
    // surf:
    // pb_factor = C5 * sqrt(time_step)
    //
    // therefore, when pb_factor would double from probability for A to probability for A.A:
    // pb_factor = sqrt(time_step_A)
    // pb_factor = 2*sqrt(time_step_AA) ->
    //
    // to keep the same probability factor we would need to set:
    //
    // time_step_AA = time_step_A / (2^2)
    //
    // This seems overly conservative and probably won't be needed. E.g. for
    // a molecule with 100 As -> the time step would have to be lowered by 100^2.
    //
    // We can come back to it when issues will arise, but for now,
    // let's use the minimum custom time or space step value for the space and time step settings

    if (min_custom_space_step != 0 && min_custom_time_step != 0) {
      // both are used in our complex, use the one giving smaller result
      float_t ts_time, ss_time;
      get_space_and_time_step(
          config,
          has_flag_no_finalized_check(SPECIES_CPLX_MOL_FLAG_SURF), D,
          min_custom_time_step, 0,
          ts_time, ss_time
      );

      float_t ts_space, ss_space;
      get_space_and_time_step(
          config,
          has_flag_no_finalized_check(SPECIES_CPLX_MOL_FLAG_SURF), D,
          0, min_custom_space_step,
          ts_space, ss_space
      );

      if (ts_time <= ts_space) {
        assert(ss_time <= ss_space);
        time_step = ts_time;
        space_step = ss_time;
      }
      else {
        assert(ss_space < ss_time);
        time_step = ts_space;
        space_step = ss_space;
      }
    }
    else {
      // one of them is 0
      get_space_and_time_step(
          config,
          has_flag_no_finalized_check(SPECIES_CPLX_MOL_FLAG_SURF), D,
          min_custom_time_step, min_custom_space_step,
          time_step, space_step
      );
    }
  }
  else {
    // use default computation based on diffusion constant of the whole complex
    space_step = get_default_space_step(config, D);
    time_step = DEFAULT_TIME_STEP;
  }
}


void Species::compute_diffusion_constant_and_space_time_step(const BNGConfig& config) {
  if (elem_mols.size() == 1) {
    // nothing to compute if we have just one molecule instance
    const ElemMolType& mt = bng_data->get_elem_mol_type(elem_mols[0].elem_mol_type_id);
    D = mt.D;
    assert(is_reactive_surface() || D != FLT_INVALID);
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
      for (const ElemMol& mi: elem_mols) {
        if (mi.is_surf()) {
          float_t mol_type_D = bng_data->get_elem_mol_type(mi.elem_mol_type_id).D;

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
      for (const ElemMol& mi: elem_mols) {
        assert(!mi.is_surf());
        float_t mol_type_D = bng_data->get_elem_mol_type(mi.elem_mol_type_id).D;

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

  compute_space_and_time_step(config);
}


void Species::dump(const string ind) const {
  cout << ind << "species_id: \t\t" << id << " [uint16_t] \t\t/* Unique ID for this species */\n";
  cout << ind << "name: *\t\t" << name << " [string] \t\t/* Symbol table entry (name) */\n";
  cout << ind << "D: \t\t" << f_to_str(D) << " [float_t] \t\t/* Diffusion constant */\n";
  cout << ind << "space_step: \t\t" << space_step << " [float_t] \t\t/* Characteristic step length */\n";
  cout << ind << "time_step: \t\t" << time_step << " [float_t] \t\t/* Minimum (maximum?) sensible timestep */\n";
  cout << ind << "flags: \t\t" << BaseSpeciesCplxMolFlag::to_str() << "\n";
  cout << ind << "CplxInstance:\n";
  Cplx::dump(true, ind + "  ");
  cout << "\n";
}


// TODO: do not use a template here
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
       [&v](size_t i1, size_t i2) {return v[i1]->name < v[i2]->name;});

  return idx;
}


void Species::dump_array(const SpeciesVector& vec, const bool sorted) {
  cout << "Species array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  if (sorted) {
    // dump sorted by name
    vector<size_t> sorted_indices = sort_indexes(vec);
    for (auto i: sorted_indices) {
      cout << i << ":\n";
      vec[i]->dump("  ");
    }
  }
  else {
    for (size_t i = 0; i < vec.size(); i++) {
      cout << i << ":\n";
      vec[i]->dump("  ");
    }
  }
}

} // namespace mcell
