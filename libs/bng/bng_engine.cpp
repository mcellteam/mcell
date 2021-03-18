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

#define BOOST_POOL_NO_MT
#include <boost/pool/pool_alloc.hpp>
#undef BOOST_POOL_NO_MT

using namespace std;

namespace BNG {

string BNGEngine::get_stats_report() const {
  stringstream res;

  std::set<reactant_class_id_t> active_reactant_classes;
  uint num_active_species = 0;
  for (const Species* s: all_species.get_species_vector()) {
    release_assert(s != nullptr);
    if (s->was_instantiated()) {
      num_active_species++;
      if (s->has_valid_reactant_class_id()) {
        active_reactant_classes.insert(s->get_reactant_class_id());
      }
    }
  }

  res << "[" <<
      "active/total species " << num_active_species << "/" << all_species.get_species_vector().size() <<
      ", rxn classes " << all_rxns.get_num_rxn_classes() <<
      ", active/total reactant classes " << active_reactant_classes.size() << "/" << all_rxns.get_num_existing_reactant_classes() <<
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



std::string BNGEngine::export_to_bngl(
    std::ostream& out_parameters,
    std::ostream& out_molecule_types,
    std::ostream& out_reaction_rules,
    const float_t volume_um3) const {

  out_parameters << IND << PARAM_V << " " << f_to_str(volume_um3) << " * 1e-15 # volume in litres\n";

  export_molecule_types_as_bngl(out_parameters, out_molecule_types);
  string err_msg = export_reaction_rules_as_bngl(out_parameters, out_reaction_rules);

  return err_msg;
}


void BNGEngine::export_molecule_types_as_bngl(std::ostream& out_parameters, std::ostream& out_molecule_types) const {
  out_molecule_types << BEGIN_MOLECULE_TYPES << "\n";

  out_parameters << "\n" << BNG::IND << "# diffusion constants\n";
  for (const ElemMolType& mt: data.get_elem_mol_types()) {
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


std::string BNGEngine::export_reaction_rules_as_bngl(
    std::ostream& out_parameters,
    std::ostream& out_reaction_rules) const {
  out_reaction_rules << BEGIN_REACTION_RULES << "\n";

  out_parameters << "\n" << BNG::IND << "# parameters to control rates in MCell and BioNetGen\n";
  out_parameters << IND << PARAM_NA_V << " " << NA_VALUE_STR << " * " << PARAM_V << "\n";
  out_parameters << IND << PARAM_VOL_RXN << " 1\n";
  out_parameters << IND << MCELL_REDEFINE_PREFIX << PARAM_VOL_RXN << " " << PARAM_NA_V << "\n";
  out_parameters << "\n" << BNG::IND << "# reaction rates\n";

  for (size_t i = 0; i < get_all_rxns().get_rxn_rules_vector().size(); i++) {
    const RxnRule* rr = get_all_rxns().get_rxn_rules_vector()[i];

    string rxn_as_bngl = rr->to_str(false, false, false);

    string rate_param = "k" + to_string(i);
    out_parameters << IND << rate_param << " " << f_to_str(rr->base_rate_constant);
    if (rr->is_bimol_vol_rxn()) {
      out_parameters  << " / NA_V * " << PARAM_VOL_RXN;
    }
    else if (rr->is_surf_rxn() || rr->is_reactive_surface_rxn()) {
      return "Export of surface reactions to BNGL is not supported yet, error for " + rxn_as_bngl + ".";
    }
    out_parameters << "\n";


    out_reaction_rules << IND << rxn_as_bngl;
    out_reaction_rules << " " << rate_param << "\n";
  }

  out_reaction_rules << END_REACTION_RULES << "\n";
  return "";
}


} // namespace BNG
