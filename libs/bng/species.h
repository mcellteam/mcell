/*
 * cplx_species.h
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_SPECIES_H_
#define LIBS_BNG_SPECIES_H_

#include "bng/cplx_instance.h"
#include "bng/defines_shared.h"

namespace BNG {

class Species;
typedef small_vector<Species> SpeciesVector;

class Species: public CplxInstance {
public:
  Species()
    : species_id(SPECIES_ID_INVALID), D(FLT_INVALID),
      // MCell-specific
      space_step(FLT_INVALID), time_step(TIME_INVALID)
    {
  }

  species_id_t species_id;

  std::string name; // string representation of the complex instance

  float_t D; // diffusion constant

  bool equal_except_for_id_base(const Species& s2) const {
    return
        CplxInstance::equal(s2) &&  // FIXME: do we really want to distinguish by orientation at this stage?
        name == s2.name &&
        D == s2.D;;
  }

  // ----------- MCell-specific -----------
  float_t space_step;
  float_t time_step; // in standard time

  bool has_count_contents_flag() const {
    return has_flag(SPECIES_FLAG_COUNT_CONTENTS);
  }

  bool has_count_enclosed_flag() const {
    return has_flag(SPECIES_FLAG_COUNT_ENCLOSED);
  }

  // true if can interact with edge of an border
  bool can_interact_with_border() const {
    return has_flag(SPECIES_FLAG_CAN_REGION_BORDER);
  }

  bool can_diffuse() const {
    return D != 0;
  }

  float_t get_time_step() const {
    return time_step;
  }

  float_t get_space_step() const {
    return space_step;
  }

  void dump(const BNGData& bng_data, const std::string ind = "") const;
  static void dump_array(const BNGData& bng_data, const SpeciesVector& vec);

  // not virtual
  bool equal_except_for_id(const Species& s2) const {
    return
        equal_except_for_id_base(s2) &&
        space_step == s2.space_step &&
        time_step == s2.time_step &&
        get_flags() == s2.get_flags();
  }

  // ^^^^^^^^^^ MCell-specific ^^^^^^^^^^
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_H_ */
