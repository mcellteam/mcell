/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#ifndef SRC4_MOLECULE_H_
#define SRC4_MOLECULE_H_

#include "bng/bng.h"

#include "defines.h"


namespace MCell {

class Partition;

#if 0
// from mcell3, copied for reference
#define TYPE_SURF 0x001
#define TYPE_VOL 0x002
#define TYPE_MASK 0x003

#define ACT_DIFFUSE 0x008
#define ACT_REACT 0x020

// IN_VOLUME
// IN_SURFACE

// #define ACT_NEWBIE 0x040  // does not have unimolecular time specified
// #define ACT_CHANGE 0x080
// -> these two above were merged into: ACT_RESCHEDULE_UNIMOL_RX

#define ACT_CLAMPED 0x1000

/* Flags telling us which linked lists the molecule appears in. */
#define IN_SCHEDULE 0x100
#define IN_SURFACE 0x200
#define IN_VOLUME 0x400
/* And a mask to pick off all three IN_ flags */
#define IN_MASK 0x700

/* Flags telling us what our counting status is */
#define COUNT_ME 0x800

/* Flag indicating that a molecule is old enough to take the maximum timestep */
#define MATURE_MOLECULE 0x2000
#endif



enum molecule_flag_t {
  // volume/surface information is only cached from BNG CplxInstance
  MOLECULE_FLAG_SURF = 1 << 0, // originally TYPE_SURF
  MOLECULE_FLAG_VOL = 1 << 1, // originally TYPE_VOL

  MOLECULE_FLAG_SCHEDULE_UNIMOL_RXN = 1 << 16,
  MOLECULE_FLAG_RESCHEDULE_UNIMOL_RXN_ON_NEXT_RXN_RATE_UPDATE = 1 << 17,

  MOLECULE_FLAG_DEFUNCT = 1 << 31,
};

/**
 * Base class for all molecules.
 */
// BNG: do we need more information from BNG or species suffices?
// TODO: remove orientation -> get this info from species
class Molecule {
public:
  Molecule()
    : id(MOLECULE_ID_INVALID), species_id(SPECIES_ID_INVALID), flags(0),
      diffusion_time(TIME_INVALID), unimol_rx_time(TIME_FOREVER) {
  }

  Molecule(const Molecule& m) {
    *this = m;
  }

  // volume molecule
  Molecule(
      const molecule_id_t id_, const species_id_t species_id_,
      const Vec3& pos_
    )
    : id(id_), species_id(species_id_), flags(MOLECULE_FLAG_VOL),
      diffusion_time(TIME_INVALID), unimol_rx_time(TIME_INVALID) {
    v.pos = pos_;
    v.subpart_index = SUBPART_INDEX_INVALID;
    v.reactant_subpart_index = SUBPART_INDEX_INVALID;
    v.counted_volume_index = COUNTED_VOLUME_INDEX_INVALID;
  }

  // surface molecule
  Molecule(
      const molecule_id_t id_, const species_id_t species_id_,
      const Vec2& pos2d
    )
    : id(id_), species_id(species_id_), flags(MOLECULE_FLAG_SURF),
      diffusion_time(TIME_INVALID), unimol_rx_time(TIME_INVALID) {
    s.pos = pos2d;
    s.orientation = ORIENTATION_NONE;
    s.wall_index = WALL_INDEX_INVALID;
    s.grid_tile_index = TILE_INDEX_INVALID;
  }

  // assuming that this function has no virtual methods and has only POD types
  // gcc9 reports a warning but, the memcpy here is safe
  void operator = (const Molecule& m) {
    memcpy(this, &m, sizeof(Molecule));
  }

  // data is ordered to avoid alignment holes (for 64-bit floats)
  molecule_id_t id; // unique molecule id (for now it is unique per partition but should be world-wide unique)
  species_id_t species_id;
  uint flags;

  // time for which it was scheduled, based on this value Partition creates 'ready list'
  // for DiffuseAndReactEvent
  float_t diffusion_time;

  // time assigned for unimol rxn, TIME_INVALID if time has not been set or mo.ecule has not unimol rxn,
  // TIME_FOREVER if the probability of an existing unimol rxn is 0
  float_t unimol_rx_time;

  // update assignment operator when modifying this
  union {
    // volume molecule data
    struct {
      Vec3 pos;
      subpart_index_t subpart_index;
      // during diffusion the molecules' subpart index might change but the reactant_subpart_index
      // stays the same until its moved in the Partition's volume_molecule_reactants_per_subpart[] array
      subpart_index_t reactant_subpart_index;
      geometry_object_id_t counted_volume_index;
    } v;

    // surface molecule data
    struct {
      Vec2 pos;
      // we probably do not want subpart index, wall index serves this purpose
      //subpart_index_t subpart_index;
      orientation_t orientation;
      wall_index_t wall_index;
      tile_index_t grid_tile_index;
    } s;
  };

  bool has_flag(uint flag) const {
    assert(__builtin_popcount(flag) == 1);
    return (flags & flag) != 0;
  }

  void set_flag(uint flag) {
    assert(__builtin_popcount(flag) == 1);
    flags = flags | flag;
  }

  void clear_flag(uint flag) {
    assert(__builtin_popcount(flag) == 1);
    flags = flags & ~flag;
  }

  bool is_vol() const {
    bool res = has_flag(MOLECULE_FLAG_VOL);
    if (res) {
      assert(!is_surf());
    }
    return res;
  }

  bool is_surf() const {
    bool res = has_flag(MOLECULE_FLAG_SURF);
    if (res) {
      assert(!is_vol());
    }
    return res;
  }

  bool is_defunct() const {
    // TODO LATER: molecules with release_delay > 0 were actually not released yet,
    // but is seems that mcell3 does not care, so let's keep it like this for now
    return has_flag(MOLECULE_FLAG_DEFUNCT) != 0;
  }

  void set_is_defunct() {
    assert(!is_defunct() && "We really should not be defuncting one molecule multiple times");
    flags |= MOLECULE_FLAG_DEFUNCT;
  }

  orientation_t get_orientation() const {
    if (is_surf()) {
      return s.orientation;
    }
    else {
      return ORIENTATION_NONE; // or not set?
    }
  }

  void dump(const std::string ind) const;
  void dump(
      const Partition& p,
      const std::string extra_comment,
      const std::string ind = "  ",
      const uint64_t iteration = 0,
      const float_t time = 0,
      const bool print_position = true
  ) const;
  std::string to_string() const;
  static void dump_array(const std::vector<Molecule>& vec);
};

} // namespace mcell

#endif // SRC4_MOLECULE_H_
