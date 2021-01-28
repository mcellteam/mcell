/*
 * bng_converter.h
 *
 *  Created on: Jan 28, 2021
 *      Author: ahusar
 */

#ifndef LIBMCELL_API_BNG_CONVERTER_H_
#define LIBMCELL_API_BNG_CONVERTER_H_

#include <string>
#include "bng/bng_defines.h"

namespace BNG {
class BNGData;
class BNGConfig;
class Component;
class ElemMol;
class Cplx;
}

namespace MCell {
namespace API {

class ComponentType;
class Component;
class ElementaryMoleculeType;
class ElementaryMolecule;
class Complex;


class BNGConverter {
public:
  BNGConverter(BNG::BNGData& bng_data_, const BNG::BNGConfig& bng_config_) :
    bng_data(bng_data_), bng_config(bng_config_) {
  }


  BNG::component_type_id_t convert_component_type(
      const std::string& elem_mol_type_name, API::ComponentType& api_ct);
  BNG::Component convert_component_instance(
      const std::string& elem_mol_type_name, API::Component& api_ci);

  BNG::elem_mol_type_id_t convert_elementary_molecule_type(
      API::ElementaryMoleculeType& mt, const bool in_rxn_or_observables = false);
  BNG::ElemMol convert_molecule_instance(
      API::ElementaryMolecule& mi, const bool in_rxn_or_observables = false);

  BNG::Cplx convert_complex(
      API::Complex& inst, const bool in_observables = false, const bool in_rxn = false);


private:

  BNG::BNGData& bng_data;
  const BNG::BNGConfig& bng_config;
};

} // namespace API
} // namespace MCell

#endif /* LIBMCELL_API_BNG_CONVERTER_H_ */
