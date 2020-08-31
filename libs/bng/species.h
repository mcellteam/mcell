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
      custom_time_step(0), custom_space_step(0),
      space_step(FLT_INVALID), time_step(TIME_INVALID),
      color_set(false), color_r(1), color_g(0), color_b(0), scale(1),
      rxn_flags_were_updated(false) {
  }

  // create species from a complex instance
  // id is not set and name is determined automatically
  Species(
      const CplxInstance& cplx_inst, const BNGData& data, const BNGConfig& config,
      const bool do_update_diffusion_constant = true)
    : CplxInstance(&data),
      id(SPECIES_ID_INVALID), D(FLT_INVALID),
      custom_time_step(0), custom_space_step(0),
      space_step(FLT_INVALID), time_step(TIME_INVALID),
      color_set(false), color_r(1), color_g(0), color_b(0), scale(1),
      rxn_flags_were_updated(false) {

    mol_instances = cplx_inst.mol_instances;
    // the only finalize method, but showing that we are finalizing
    // just the CplxInstance part of the Species
    CplxInstance::finalize();
    if (do_update_diffusion_constant) {
      update_diffusion_constant(data, config);
    }
    set_flag(BNG::SPECIES_FLAG_CAN_DIFFUSE, D != 0);
    finalize();
    name = cplx_inst.to_str(data);
  }

  // we need explicit copy ctor to call CplxInstance's copy ctor
  Species(const Species& other)
    : CplxInstance(other),
      id(other.id), name(other.name), D(other.D),
      custom_time_step(other.custom_time_step), custom_space_step(other.custom_space_step),
      space_step(other.space_step), time_step(other.time_step),
      color_set(other.color_set), color_r(other.color_r), color_g(other.color_g), color_b(other.color_b), scale(other.scale),
      rxn_flags_were_updated(other.rxn_flags_were_updated) {
  }

  // TODO: maybe an assignment operator is needed, e.g. in the CplxInstance case, the copy ctor was not
  // called all the time (not sure why...) and assignment operator was needed to fix an issue

  // TODO: why is this called from the Species ctor, can we remove it?
  void finalize() {
    CplxInstance::finalize();
    set_flag(BNG::SPECIES_FLAG_CAN_DIFFUSE, D != 0);
  }

  void canonicalize(const BNGData& bng_data) {
    CplxInstance::canonicalize(); // calls also CplxInstance::finalize
    name = to_str(bng_data);
  }

  species_id_t id;

  std::string name; // string representation of the complex instance

  float_t D; // diffusion constant, entered by user in MCell3 mode, computed in MCell4 BNG model

  // when the user supplied a custom step, the attribute
  // is set to non-zero value, max one of them can be set to a non-zero value
  float_t custom_time_step;
  float_t custom_space_step;

  // ----------- MCell-specific -----------

  float_t space_step;
  float_t time_step; // in standard time

  // sets SPECIES_FLAG_CAN_VOLVOL, SPECIES_FLAG_CAN_VOLSURF, SPECIES_FLAG_CAN_VOLWALL,
  // SPECIES_FLAG_CAN_SURFSURF, and/or SPECIES_FLAG_CAN_REGION_BORDER
  // flags according to reactions in the system
  bool are_rxn_and_custom_flags_uptodate() const {
    return rxn_flags_were_updated;
  }

  void update_rxn_and_custom_flags(
      const SpeciesContainer& all_species, RxnContainer& all_rxns,
      const BaseCustomFlagsAnalyzer* flags_analyzer = nullptr
  );

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

  bool has_unimol_rxn() const {
    return has_flag(SPECIES_FLAG_HAS_UNIMOL_RXN);
  }

  bool has_bimol_vol_rxn() const {
    return has_flag(SPECIES_FLAG_HAS_BIMOL_VOL_RXN);
  }

  bool can_vol_react() const {
    if (cant_initiate()) {
      return false;
    }
    else {
      return has_bimol_vol_rxn();
    }
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
