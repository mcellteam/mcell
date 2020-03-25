#ifndef LIBS_BNG_BASE_FLAG_H_
#define LIBS_BNG_BASE_FLAG_H_

#include "bng_defines.h"


enum cplx_mol_flag_t {
  CPLX_MOL_FLAG_SURF = 1 << 0,
  CPLX_MOL_FLAG_VOL = 1 << 1,

  CPLX_FLAG_HAS_SINGLE_ORIENTATION = 1 << 2,
  CPLX_FLAG_SINGLE_ORIENTATION_IS_UP = 1 << 3,
  CPLX_FLAG_ONE_MOL_NO_COMPONENTS = 1 << 4,
};


// NOTE: maybe use bitfield instead?
class BaseFlag {
private:
  bool finalized;
  uint flags;

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
    assert(finalized);
    assert(__builtin_popcount(flag) == 1);
    if (value) {
      flags = flags | flag;
    }
    else {
      clear_flag(flag);
    }
  }

  void clear_flag(uint flag) {
    assert(finalized);
    assert(__builtin_popcount(flag) == 1);
    flags = flags & ~flag;
  }

  void set_finalized() {
    finalized = true;
  }
};


#endif // LIBS_BNG_BASE_FLAG_H_
