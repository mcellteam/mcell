/*
 * cplx_species.h
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_SPECIES_H_
#define LIBS_BNG_SPECIES_H_

#include "common_defines.h"

#include "cplx_instance.h"

#ifdef BNG_ONLY_MCELL
//#include "mcell_structs.h"
#endif

#ifdef BNG_ONLY_MCELL

namespace MCell {
// same as in mcell_structs but renamed to make sure it is used correctly
enum species_flag_t {

  SPECIES_FLAG_ON_GRID = 0x01,
  SPECIES_FLAG_IS_SURFACE = 0x02,
  SPECIES_FLAG_CAN_VOLVOL = 0x10, // can vol vol react?
  SPECIES_FLAG_CAN_VOLSURF = 0x20,
  SPECIES_FLAG_CAN_SURFSURF = 0x80,
  SPECIES_FLAG_CANT_INITIATE = 0x400, // must not be set, not sure what to do with this yet (at least for some cases)
  SPECIES_FLAG_CAN_SURFSURFSURF = 0x20000, // 0x20000 - not supported
  SPECIES_FLAG_SET_MAX_STEP_LENGTH = 0x80000,
  SPECIES_FLAG_CAN_REGION_BORDER = 0x100000, // CAN_REGION_BORDER, // 0x100000
  SPECIES_FLAG_EXTERNAL_SPECIES = 0x400000 // 0x400000 - not supported
};
}
#endif // BNG_ONLY_MCELL

namespace BNG {

class Species;
typedef small_vector<Species> SpeciesVector;

class Species: public CplxInstance {
public:
  Species()
    : species_id(SPECIES_ID_INVALID), D(FLT_INVALID)
#ifdef BNG_ONLY_MCELL
       , space_step(FLT_INVALID), time_step(TIME_INVALID), flags(0)
#endif // BNG_ONLY_MCELL
    {
  }

  species_id_t species_id;

  std::string name; // string representation of the complex instance

  float_t D; // diffusion constant

  bool equal_except_for_id_base(const Species& s2) const {
    return
        CplxInstance::equal(s2) &&
        name == s2.name &&
        D == s2.D;;
  }

#ifdef BNG_ONLY_MCELL
  float_t space_step;
  float_t time_step; // in standard time

  uint flags;

  bool has_flag(uint flag) const {
    assert(__builtin_popcount(flag) == 1);
    return (flags & flag) != 0;
  }

  bool is_surf() const {
    return has_flag(MCell::SPECIES_FLAG_ON_GRID);
  }

  bool is_vol() const {
    return !has_flag(MCell::SPECIES_FLAG_ON_GRID);
  }

  bool is_reactive_surface() const {
    return has_flag(MCell::SPECIES_FLAG_IS_SURFACE);
  }

  // true if can interact with edge of an border
  bool can_interact_with_border() const {
    return has_flag(MCell::SPECIES_FLAG_CAN_REGION_BORDER);
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
        flags == s2.flags;
  }
#endif // BNG_ONLY_MCELL

};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_H_ */
