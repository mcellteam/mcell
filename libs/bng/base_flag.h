#ifndef LIBS_BNG_BASE_FLAG_H_
#define LIBS_BNG_BASE_FLAG_H_

#include "bng/bng_defines.h"

namespace BNG {

// this single enum defines flags for species, complex instances and molecule instances
// TODO: rename to somethign shorter
enum species_cplx_mol_rxn_flag_t {

  // maintaining the same values as in MCell

  // if both of the following flags are unset, the a volume molecule/species
  // the type of molecule will be bits when bitfield will be used
  SPECIES_CPLX_MOL_FLAG_SURF = 0x01, // surface or vol mol
  SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE = 0x02, // reactive surface

  // simple species
  SPECIES_CPLX_FLAG_ONE_MOL_NO_COMPONENTS = 0x04,

  SPECIES_FLAG_CAN_VOLVOL = 0x10, // can vol vol react? (unused for now)
  SPECIES_FLAG_CAN_VOLSURF = 0x20, // can vol-surf react
  SPECIES_FLAG_CAN_VOLWALL = 0x40, // can vol-surface react
  SPECIES_FLAG_CAN_SURFSURF = 0x80, // can surf-surf react
  //SPECIES_FLAG_CAN_SURFWALL = 0x100, // can surf-surface react
  SPECIES_FLAG_CANT_INITIATE = 0x400, // must not be set, not sure what to do with this yet (at least for some cases)
  SPECIES_FLAG_COUNT_ENCLOSED = 0x8000, // this species is marked to be counted inside of a volume
  SPECIES_FLAG_CAN_SURFSURFSURF = 0x20000, // 0x20000 - not supported - TODO LATER: remove
  SPECIES_FLAG_SET_MAX_STEP_LENGTH = 0x80000,
  SPECIES_FLAG_CAN_REGION_BORDER = 0x100000,
  SPECIES_FLAG_EXTERNAL_SPECIES = 0x400000, // 0x400000 - not supported - TODO LATER: remove

  RXN_FLAG_COUNTED = 0x1000000,
};


// use bitfield?
class BaseFlag {
private:
  bool finalized;

  union {
    uint flags;
    struct {

    };
  };

public:
  BaseFlag()
    : finalized(false), flags(0) {
  }

  bool has_flag(uint flag) const {
    assert(finalized);
    assert(__builtin_popcount(flag) == 1);
    return (flags & flag) != 0;
  }

  void set_flag(uint flag, bool value = true) {
    assert(__builtin_popcount(flag) == 1);
    if (value) {
      flags = flags | flag;
    }
    else {
      clear_flag(flag);
    }
  }

  void clear_flag(uint flag) {
    assert(__builtin_popcount(flag) == 1);
    flags = flags & ~flag;
  }

  void set_finalized() {
    finalized = true;
  }

  bool is_finalized() const {
    return finalized;
  }

  uint get_flags() const {
    return flags;
  }

  void set_flags(uint value) {
    flags = value;
  }
};


class BaseSpeciesCplxMolFlag: public BaseFlag {
public:
  void set_is_vol() {
    clear_flag(SPECIES_CPLX_MOL_FLAG_SURF);
    clear_flag(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE);
  }

  void set_is_surf() {
    set_flag(SPECIES_CPLX_MOL_FLAG_SURF);
    clear_flag(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE);
  }

  void set_is_reactive_surface() {
    clear_flag(SPECIES_CPLX_MOL_FLAG_SURF);
    set_flag(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE);
  }

  bool is_surf() const {
    return has_flag(SPECIES_CPLX_MOL_FLAG_SURF);
  }

  bool is_vol() const {
    return !has_flag(SPECIES_CPLX_MOL_FLAG_SURF) && !has_flag(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE);
  }

  bool is_reactive_surface() const {
    return has_flag(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE);
  }
};

} // namespace BNG

#endif // LIBS_BNG_BASE_FLAG_H_
