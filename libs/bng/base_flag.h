#ifndef LIBS_BNG_BASE_FLAG_H_
#define LIBS_BNG_BASE_FLAG_H_

#include "bng/bng_defines.h"

namespace BNG {

// this single enum defines flags for species, complexes and molecule instances
enum species_cplx_mol_rxn_flag_t {

  // maintaining the same values as in MCell

  // if both of the following flags are unset, the a volume molecule/species
  // the type of molecule will be bits when bitfield will be used
  SPECIES_CPLX_MOL_FLAG_SURF = 0x01, // surface or vol mol
  SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE = 0x02, // reactive surface

  // simple species
  SPECIES_CPLX_FLAG_ONE_MOL_NO_COMPONENTS = 0x04,

  // species or complex has canonical representation
  SPECIES_CPLX_FLAG_IS_CANONICAL = 0x08,

  SPECIES_FLAG_CAN_VOLVOL = 0x10, // can vol vol react? (unused for now)
  SPECIES_FLAG_CAN_VOLSURF = 0x20, // can vol-surf react
  SPECIES_FLAG_CAN_VOLWALL = 0x40, // can vol-surface react
  SPECIES_FLAG_CAN_SURFSURF = 0x80, // can surf-surf react
  SPECIES_FLAG_CAN_SURFWALL = 0x100, // can surf-surface react
  SPECIES_FLAG_CAN_INTERMEMBRANE_SURFSURF = 0x200, // can surf-surf react across intermembrane space
  SPECIES_MOL_FLAG_TARGET_ONLY = 0x400, // this molecule may not trigger a reaction with another molecule
  SPECIES_FLAG_CAN_DIFFUSE = 0x800, // value used as COUNT_TRIGGER in MCell3
  
  SPECIES_FLAG_IS_DEFUNCT = 0x1000, // these species were removed

  // these species were created on runtime and therefore nothing references them and if so, then
  // the code that references them must be able to handle their deletion
  SPECIES_FLAG_IS_REMOVABLE = 0x2000,

  SPECIES_FLAG_CAN_SURFSURFSURF = 0x20000, // not supported - TODO LATER: remove
  SPECIES_FLAG_SET_MAX_STEP_LENGTH = 0x80000, // not supported

  SPECIES_FLAG_CAN_REGION_BORDER = 0x100000,
  SPECIES_FLAG_HAS_UNIMOL_RXN = 0x200000, // collides with mcell3's REGION_PRESENT, but this flag is not copied
  SPECIES_FLAG_HAS_BIMOL_VOL_RXN = 0x400000, // collides with mcell3's SPECIES_FLAG_EXTERNAL_SPECIES, but this flag is not supported at alkl

  // set when the count of instances of this species gets higher than 0
  // may be reset in RxnClassCleanup event
  SPECIES_FLAG_WAS_INSTANTIATED = 0x800000,

  RXN_FLAG_COUNTED_IN_WORLD = 0x1000000,
  RXN_FLAG_COUNTED_IN_VOLUME_REGIONS = 0x2000000,
  RXN_FLAG_COUNTED_ON_SURFACE_REGIONS = 0x4000000,
  RXN_FLAG_SIMPLE = 0x8000000, // reactants and products are only simple complexes
  RXN_FLAG_MAY_PRODUCE_MUTLIPLE_IDENTICAL_PRODUCTS = 0x10000000,
  RXN_FLAG_CREATED_FOR_CONCENTRATION_CLAMP = 0x20000000,
  RXN_FLAG_CREATED_FOR_FLUX_CLAMP = 0x40000000,
  RXN_FLAG_INTERMEMBRANE = 0x80000000
};


class BaseFlag {
private:
  bool finalized;

  uint flags;

public:
  BaseFlag()
    : finalized(false), flags(0) {
  }

  virtual ~BaseFlag() {
  }

  // Species has an extra check that when asked for rxn flags,
  // it checks whether the rxn flags were initialized,
  // the rxn flags must be initialized in runtime because
  // they check reactions
  virtual bool has_flag(uint flag) const {
    assert(finalized);
    assert(__builtin_popcount(flag) == 1);
    return (flags & flag) != 0;
  }

  // flag retrieval for circumstances when the flags were not finalized
  // but we know that it is ok to ask, should be used only when necessary
  virtual bool has_flag_no_finalized_check(uint flag) const {
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

  void add_flags(uint flags_to_set) {
    flags = flags | flags_to_set;
  }

  void clear_flag(uint flag) {
    assert(__builtin_popcount(flag) == 1);
    flags = flags & ~flag;
  }

  // flags is a mask of flags to be cleared
  void clear_flags(uint flags_to_clear) {
    flags = flags & ~flags_to_clear;
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

  bool is_vol_no_finalized_check() const {
    return !has_flag_no_finalized_check(SPECIES_CPLX_MOL_FLAG_SURF) &&
        !has_flag_no_finalized_check(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE);
  }

  bool is_reactive_surface() const {
    return has_flag(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE);
  }

  std::string to_str() const;
};


class Species;

} // namespace BNG

#endif // LIBS_BNG_BASE_FLAG_H_
