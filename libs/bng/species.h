/*
 * cplx_species.h
 *
 *  Created on: Mar 24, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_SPECIES_H_
#define LIBS_BNG_SPECIES_H_

#include "bng/elem_mol_type.h"
#include "bng/cplx.h"
#include "bng/shared_defines.h"

namespace BNG {

class Species;
class SpeciesContainer;
class RxnContainer;

typedef std::vector<Species*> SpeciesVector;

class Species: public Cplx, public ElemMolTypeSpeciesCommonData {
public:
  species_id_t id;
  std::string name; // string representation of the complex

  // ----------- MCell-specific -----------
  float_t space_step;
  float_t time_step; // in standard time

  // must call finalize afterwards
  Species(const BNGData& data)
    : Cplx(&data),
      id(SPECIES_ID_INVALID),
      space_step(FLT_INVALID), time_step(TIME_INVALID),
      rxn_flags_were_updated(false), num_instantiations(0),
      reactant_class_id(REACTANT_CLASS_ID_INVALID) {
  }

  void finalize_species(const BNGConfig& config, const bool update_diffusion_constant = true) {
    // species must not use IN/OUT, remove it automatically when defining species
    for (auto& em: elem_mols) {
      if (is_in_out_compartment_id(em.compartment_id)) {
        em.compartment_id = COMPARTMENT_ID_NONE;
      }
    }

    canonicalize(); // sets name as well, calls also Cplx::finalize
    set_flag(SPECIES_FLAG_CAN_DIFFUSE, D != 0);
    if (is_reactive_surface()) {
      // surfaces are always assumed to be instantiated
      set_flag(SPECIES_FLAG_WAS_INSTANTIATED);
    }
    if (update_diffusion_constant) {
      compute_diffusion_constant_and_space_time_step(config);
    }
    set_flag(BNG::SPECIES_FLAG_CAN_DIFFUSE, D != 0); // TODO: can this be removed when we set it in finalize?
  }

  // create species from a complex instance
  // id is not set and name is determined automatically
  Species(
      const Cplx& cplx_inst, const BNGData& data, const BNGConfig& config,
      const bool update_diffusion_constant = true)
    : Cplx(&data),
      id(SPECIES_ID_INVALID),
      space_step(FLT_INVALID), time_step(TIME_INVALID),
      rxn_flags_were_updated(false), num_instantiations(0), reactant_class_id(REACTANT_CLASS_ID_INVALID) {

    elem_mols = cplx_inst.elem_mols;
    finalize_species(config, update_diffusion_constant);
  }

  // we need explicit copy ctor to call CplxInstance's copy ctor
  Species(const Species& other)
    : Cplx(other), ElemMolTypeSpeciesCommonData(other),
      id(other.id), name(other.name),
      space_step(other.space_step), time_step(other.time_step),
      rxn_flags_were_updated(other.rxn_flags_were_updated), num_instantiations(other.num_instantiations),
      reactant_class_id(other.reactant_class_id) {
  }

  // used when these species are added as new to the species container
  void reset_num_instantiations() {
    num_instantiations = 0;
  }

  // TODO: maybe an assignment operator is needed, e.g. in the CplxInstance case, the copy ctor was not
  // called all the time (not sure why...) and assignment operator was needed to fix an issue

  // default sorting of components is according to molecule types
  void canonicalize(const bool sort_components_by_name_do_not_finalize = false) {
    Cplx::canonicalize(sort_components_by_name_do_not_finalize); // calls also CplxInstance::finalize
    name = "";
    to_str(name);
  }

  // sets SPECIES_FLAG_CAN_VOLVOL, SPECIES_FLAG_CAN_VOLSURF, SPECIES_FLAG_CAN_VOLWALL,
  // SPECIES_FLAG_CAN_SURFSURF, and/or SPECIES_FLAG_CAN_REGION_BORDER
  // flags according to reactions in the system
  bool are_rxn_and_custom_flags_uptodate() const {
    return rxn_flags_were_updated;
  }

  void update_rxn_and_custom_flags(
      const SpeciesContainer& all_species, RxnContainer& all_rxns
  );

  uint get_num_instantiations() const {
    if (is_reactive_surface()) {
      // we are not counting reactive surfaces
      return 1;
    }
    else {
      // superspecies have their num_instantiations set to 1
      return num_instantiations;
    }
  }

  void inc_num_instantiations() {
    assert(!is_reactive_surface());
    num_instantiations++;
    set_was_instantiated(true);
  }

  void dec_num_instantiations() {
    assert(!is_reactive_surface());
    assert(num_instantiations > 0);
    num_instantiations--;
    // does not reset flag SPECIES_FLAG_WAS_INSTANTIATED because this flag also means that
    // species with which this species can react in a bimol rxn
    // have classes with these species in their rxn
    // TODO: improve explanation
  }

  // do not call this directly, should be called only from
  // SpeciesContainer::remove
  void set_is_defunct() {
    set_flag(SPECIES_FLAG_IS_DEFUNCT);
  }

  bool is_defunct() const {
    return has_flag(SPECIES_FLAG_IS_DEFUNCT);
  }

  void set_is_removable() {
    set_flag(SPECIES_FLAG_IS_REMOVABLE);
  }

  bool is_removable() const {
    return has_flag(SPECIES_FLAG_IS_REMOVABLE);
  }

  void set_was_instantiated(const bool value) {
    assert(!is_reactive_surface() && "Reactive surfaces are automatically considered to be instantiated");
    if (value) {
      set_flag(SPECIES_FLAG_WAS_INSTANTIATED);
    }
    else {
      clear_flag(SPECIES_FLAG_WAS_INSTANTIATED);
    }
  }

  bool was_instantiated() const {
    if (is_reactive_surface()) {
      return true;
    }
    else {
      return has_flag(SPECIES_FLAG_WAS_INSTANTIATED);
    }
  }

  // true if can interact with edge of an border
  bool can_interact_with_border() const {
    return has_flag(SPECIES_FLAG_CAN_REGION_BORDER);
  }

  // true if a molecule of this species cannot initiate a reaction
  bool is_target_only() const {
    return has_flag(SPECIES_MOL_FLAG_TARGET_ONLY);
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

  bool has_unimol_rxn() const {
    return has_flag(SPECIES_FLAG_HAS_UNIMOL_RXN);
  }

  bool has_bimol_vol_rxn() const {
    return has_flag(SPECIES_FLAG_HAS_BIMOL_VOL_RXN);
  }

  bool can_vol_react() const {
    if (is_target_only()) {
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

  bool has_valid_reactant_class_id() const {
    return reactant_class_id != REACTANT_CLASS_ID_INVALID;
  }

  reactant_class_id_t get_reactant_class_id() const {
    assert(has_valid_reactant_class_id());
    return reactant_class_id;
  }

  void set_reactant_class_id(const reactant_class_id_t id) {
    assert(id != REACTANT_CLASS_ID_INVALID);
    reactant_class_id = id;
  }

  void dump(const std::string ind = "") const;
  static void dump_array(const SpeciesVector& vec, const bool sorted = false);

  // not virtual
  bool matches_fully_ignore_name_id_and_flags(const Species& s2) const {
    // we do not want to compare by name because name is defined
    // by the complex and depends on ordering of molecules and components
    return
        Cplx::matches_fully(s2) &&
        cmp_eq(D, s2.D) &&
        cmp_eq(space_step, s2.space_step) &&
        cmp_eq(time_step, s2.time_step);
  }

  bool cplx_matches_fully_ignore_orientation_and_flags(const Cplx& cplx) const {
    return Cplx::matches_fully(cplx, true);
  }

  // use information from contained molecule types to compute diffusion constant
  // calls update_space_and_time_step
  void compute_diffusion_constant_and_space_time_step(const BNGConfig& config);

private:
  // used in initialization
  void compute_space_and_time_step(const BNGConfig& config);

  // rxn flags are updated when a molecule of this species is added to world
  bool rxn_flags_were_updated;

  uint num_instantiations;

  reactant_class_id_t reactant_class_id;
};

} // namespace BNG

#endif /* LIBS_BNG_SPECIES_H_ */
