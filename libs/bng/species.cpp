
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

  // check if any of the used elementaty molecules has custom time or space step
  bool has_custom_step = false;
  float_t min_space_step = FLT_GIGANTIC;
  elem_mol_type_id_t min_emt_id = UINT_INVALID;
#ifndef NDEBUG
  float_t min_time_step = FLT_GIGANTIC;
#endif

  for (uint i = 0; i < elem_mols.size(); i++) {
    const ElemMolType& emt = bng_data->get_elem_mol_type(elem_mols[i].elem_mol_type_id);
    assert(emt.time_step != FLT_INVALID && emt.space_step != FLT_INVALID);

    has_custom_step = has_custom_step || emt.has_custom_time_or_space_step();

    if (emt.space_step < min_space_step) {
      min_emt_id = elem_mols[i].elem_mol_type_id;

#ifndef NDEBUG
      min_time_step = emt.time_step;
#endif
    }
  }

  if (has_custom_step) {
    // custom time or space step is used either to:
    // 1) speedup simulation if the default probability of reactions is low enough and
    //    we don't care that the molecule diffuses further than usual or
    // 2) lower the probability of a reaction per collision so that it is lower than 1
    //
    // TODO
    release_assert(elem_mols.size() == 1);

    const ElemMolType& min_emt = bng_data->get_elem_mol_type(min_emt_id);

    assert(
        min_emt.time_step == min_time_step &&
        "The found elem. mol type with minimal space step must be also the one with minimal time step");

    space_step = min_emt.space_step;
    time_step = min_emt.time_step;
  }
  else {
    // use default computation based on diffusion constant
    space_step = sqrt_f(4.0 * 1.0e8 * D * config.time_unit) * config.rcp_length_unit;
    time_step = 1.0;
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


void Species::dump(const BNGData& bng_data, const string ind) const {
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


void Species::dump_array(const BNGData& bng_data, const SpeciesVector& vec, const bool sorted) {
  cout << "Species array: " << (vec.empty() ? "EMPTY" : "") << "\n";

  if (sorted) {
    // dump sorted by name
    vector<size_t> sorted_indices = sort_indexes(vec);
    for (auto i: sorted_indices) {
      cout << i << ":\n";
      vec[i]->dump(bng_data, "  ");
    }
  }
  else {
    for (size_t i = 0; i < vec.size(); i++) {
      cout << i << ":\n";
      vec[i]->dump(bng_data, "  ");
    }
  }
}

} // namespace mcell
