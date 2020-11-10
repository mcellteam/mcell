/*
 * bng_engine.cpp
 *
 *  Created on: Mar 26, 2020
 *      Author: ahusar
 */
#include <iostream>
#include <sstream>

#include "bng/bng_engine.h"
#include "bng/bngl_names.h"

using namespace std;

namespace BNG {

string BNGEngine::get_stats_report() const {
  stringstream res;

  uint num_active_species = 0;
  for (const Species& s: all_species.get_species_vector()) {
    if (s.was_instantiated()) {
      num_active_species++;
    }
  }



  res << "[" <<
      "active species " << num_active_species <<
      ", total species " << all_species.get_species_vector().size() <<
      ", total rxn classes " << all_rxns.get_num_rxn_classes() <<
      "]";
  return res.str();
}


Cplx BNGEngine::create_cplx_from_species(
    const species_id_t id, const orientation_t o, const compartment_id_t compartment_id) const {
  const Cplx& ref = all_species.get(id);
  Cplx copy = ref;
  copy.set_orientation(o);
  copy.set_compartment_id(compartment_id);
  return copy;
}



std::string BNGEngine::export_as_bngl(
    std::ostream& out_parameters,
    std::ostream& out_molecule_types,
    std::ostream& out_reaction_rules,
    const float_t volume_um3) const {

  out_parameters << IND << PARAM_V << " " << volume_um3 << " * 1e-15 # volume in um^3\n";

  export_molecule_types_as_bngl(out_parameters, out_molecule_types);
  export_reaction_rules_as_bngl(out_parameters, out_reaction_rules);

  return "";
}


void BNGEngine::export_molecule_types_as_bngl(std::ostream& out_parameters, std::ostream& out_molecule_types) const {
  out_molecule_types << BEGIN_MOLECULE_TYPES << "\n";

  for (const MolType& mt: data.get_molecule_types()) {
    if (mt.is_reactive_surface() || is_species_superclass(mt.name)) {
      continue;
    }

    // define as mol type
    out_molecule_types << IND << mt.to_str(data) << "\n";

    // and also set its diffusion constant as parameter
    if (mt.is_vol()) {
      out_parameters << IND << MCELL_DIFFUSION_CONSTANT_3D_PREFIX;
    }
    else if (mt.is_surf()){
      out_parameters << IND << MCELL_DIFFUSION_CONSTANT_2D_PREFIX;
    }
    out_parameters << mt.name << " " << f_to_str(mt.D) << "\n";
  }

  out_molecule_types << END_MOLECULE_TYPES << "\n";
}


void BNGEngine::export_reaction_rules_as_bngl(
    std::ostream& out_parameters,
    std::ostream& out_reaction_rules) const {
  out_reaction_rules << BEGIN_REACTION_RULES << "\n";

  // parameters to control rates MCell/BNG
  out_parameters << IND << PARAM_VOL_RXN << " 1\n";
  out_parameters << IND << PARAM_NA_V << " " << NA_VALUE_STR << " * " << PARAM_V << "\n";
  out_parameters << IND << MCELL_REDEFINE_PREFIX << PARAM_VOL_RXN << " " << PARAM_NA_V << "\n";

  for (size_t i = 0; i < get_all_rxns().get_rxn_rules_vector().size(); i++) {
    const RxnRule* rr = get_all_rxns().get_rxn_rules_vector()[i];

    string rate_param = "k" + to_string(i);
    out_parameters << IND << rate_param << " " << f_to_str(rr->base_rate_constant) << " / NA_V";
    if (rr->is_bimol_vol_rxn()) {
      out_parameters << " * " << PARAM_VOL_RXN << "\n";
    }
    // TODO: surf rxns

    out_reaction_rules << IND << rr->to_str(false, false, false);
    out_reaction_rules << " " << rate_param << "\n";
  }

  out_reaction_rules << END_REACTION_RULES << "\n";
}


} // namespace BNG
