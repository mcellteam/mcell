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
class SpeciesContainer;
class RxnContainer;

typedef small_vector<Species> SpeciesVector;

class Species: public CplxInstance {
public:
  Species(const BNGData& data)
    : CplxInstance(&data),
      id(SPECIES_ID_INVALID), D(FLT_INVALID),
      // MCell-specific
      space_step(FLT_INVALID), time_step(TIME_INVALID),
      color_set(false), color_r(1), color_g(0), color_b(0), scale(1),
      rxn_flags_were_updated(false)
    {
  }

  // create species from a complex instance
  // id is not set and name is determined automatically
  Species(const CplxInstance& cplx_inst, const BNGData& data, const BNGConfig& config)
    : CplxInstance(&data),
      id(SPECIES_ID_INVALID), D(FLT_INVALID),
      // MCell-specific
      space_step(FLT_INVALID), time_step(TIME_INVALID),
      color_set(false), color_r(1), color_g(0), color_b(0), scale(1),
      rxn_flags_were_updated(false)
  {

    mol_instances = cplx_inst.mol_instances;
    // the only finalize method, but showing that we are finalizing
    // just the CplxInstance part of the Species
    CplxInstance::finalize();
    update_diffusion_constant(data, config);
    name = cplx_inst.to_str(data);

  }


  void finalize() {
    // set flag on whether D is 0
    set_flag(BNG::SPECIES_FLAG_CAN_DIFFUSE, D != 0);
    CplxInstance::finalize();
  }

  species_id_t id;

  std::string name; // string representation of the complex instance

  float_t D; // diffusion constant, entered by user in MCell3 mode, computed in MCell4 BNG model

  // ----------- MCell-specific -----------
  float_t space_step;
  float_t time_step; // in standard time

  // sets SPECIES_FLAG_CAN_VOLVOL, SPECIES_FLAG_CAN_VOLSURF, SPECIES_FLAG_CAN_VOLWALL,
  // SPECIES_FLAG_CAN_SURFSURF, and/or SPECIES_FLAG_CAN_REGION_BORDER
  // flags according to reactions in the system
  void update_flags_based_on_rxns(const SpeciesContainer& all_species, RxnContainer& all_rxns, const bool force_update = false);

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
    return has_flag(SPECIES_FLAG_CAN_DIFFUSE);
  }

  float_t get_time_step() const {
    return time_step;
  }

  float_t get_space_step() const {
    return space_step;
  }

  bool needs_counted_colume() const {
    return has_flag(SPECIES_FLAG_NEEDS_COUNTED_VOLUME);
  }

  bool has_flag(uint flag) const override {
    // check that rxn flags are up-to-date
    if (flag == SPECIES_FLAG_CAN_VOLVOL ||
        flag == SPECIES_FLAG_CAN_VOLSURF ||
        flag == SPECIES_FLAG_CAN_VOLWALL ||
        flag == SPECIES_FLAG_CAN_SURFSURF
    ) {
      assert(rxn_flags_were_updated);
    }
    return BaseFlag::has_flag(flag);
  }

  void dump(const BNGData& bng_data, const std::string ind = "") const;
  static void dump_array(const BNGData& bng_data, const SpeciesVector& vec, const bool sorted = false);

  // not virtual
  bool matches_fully_ignore_name_id_and_flags(const Species& s2) const {
    // we do not want to compare by name because name is defined
    // by the complex and depends on ordering of molecules and components
    return
        CplxInstance::matches_fully(s2) &&
        D == s2.D &&
        space_step == s2.space_step &&
        time_step == s2.time_step;
  }

  bool cplx_matches_fully_ignore_orientation_and_flags(const CplxInstance& cplx_inst) const {
    return CplxInstance::matches_fully(cplx_inst, true);
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

private:
  // rxn flags are updated when a molecule of this species is added to world
  bool rxn_flags_were_updated;
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_H_ */
