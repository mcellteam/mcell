/*
 * cplx_species.h
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_CPLX_SPECIES_H_
#define LIBS_BNG_CPLX_SPECIES_H_

#include "cplx_instance.h"

namespace BNG {

class CplxSpecies : public CplxInstance {
public:
  CplxSpecies()
    : species_id(SPECIES_ID_INVALID), D(FLT_INVALID) {
  }

  species_id_t species_id;

  std::string name; // string representation of the complex instance

  float_t D; // diffusion constant

  bool equal_except_for_id(const CplxSpecies& s2) {
    return
        CplxInstance::equal(s2) &&
        name == s2.name &&
        D == s2.D;;
  }

};

} // namespace BNG

#endif /* LIBS_BNG_CPLX_SPECIES_H_ */
