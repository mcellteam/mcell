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
    : id(SPECIES_ID_INVALID), D(FLT_INVALID),
      // MCell-specific
      space_step(FLT_INVALID), time_step(TIME_INVALID),
      color_set(false), color_r(1), color_g(0), color_b(0), scale(1)
    {
  }

  // create species from a complex instance
  // id is not set and name is determined automatically
  Species(const CplxInstance& cplx_inst, const BNGData& data, const BNGConfig& config)
    : id(SPECIES_ID_INVALID), D(FLT_INVALID),
      // MCell-specific
      space_step(FLT_INVALID), time_step(TIME_INVALID),
      color_set(false), color_r(1), color_g(0), color_b(0), scale(1)  {

    mol_instances = cplx_inst.mol_instances;
    update_diffusion_constant(data, config);
    name = cplx_inst.to_str(data);
    finalize();
  }

  species_id_t id;

  std::string name; // string representation of the complex instance

  float_t D; // diffusion constant, entered by user in MCell3 mode, computed in MCell4 BNG model

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

  // true if a molecule of this species cannot initiate a reaction
  bool cant_initiate() const {
    return has_flag(SPECIES_FLAG_CANT_INITIATE);
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
  static void dump_array(const BNGData& bng_data, const SpeciesVector& vec, const bool sorted = false);

  // not virtual
  bool equal_except_for_id_and_flags(const Species& s2) const {
    return
        equal_except_for_id_base(s2) &&
        space_step == s2.space_step &&
        time_step == s2.time_step;
  }

  // for initialization
  void update_space_and_time_step(const BNGConfig& config);

  // use information from contained molecule types to compute diffusion constant
  // calls update_space_and_time_step
  void update_diffusion_constant(const BNGData& data, const BNGConfig& config);


  // ^^^^^^^^^^ MCell-specific ^^^^^^^^^^

  // visualization
  void set_color(float_t r, float_t g, float_t b) {
    color_r = r;
    color_g = g;
    color_b = b;
    color_set = true;
  }
  void set_scale(float_t s) {
    scale = s;
  }

  bool color_set;
  float_t color_r, color_g, color_b ;  // mol color default is red
  float_t scale; // scale = 1 by default
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_H_ */
