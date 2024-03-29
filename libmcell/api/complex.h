/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef API_COMPLEX_H
#define API_COMPLEX_H

#include "generated/gen_complex.h"
#include "api/api_common.h"
#include "api/api_utils.h"
#include "api/compartment_utils.h"
#include "bng/bngl_names.h"

namespace BNG {
class BNGData;
class Cplx;
}

namespace MCell {
namespace API {

class Species;

// WARNING: do not set compartment_name through attribute, use set_compartment_name
class Complex: public GenComplex {
public:
  COMPLEX_CTOR()

  static std::shared_ptr<API::Complex> construct_empty() {
    // to avoid Complex object semantic check, we need to insert a dummy name
    // when creating the object
    auto res_cplx_inst = std::make_shared<API::Complex>(STR_UNSET);
    res_cplx_inst->name = STR_UNSET;
    return res_cplx_inst;
  }

  static std::shared_ptr<API::Complex> construct_from_bng_cplx(
      const BNG::BNGData& bng_data,
      const BNG::Cplx& bng_inst);

  static std::shared_ptr<API::Complex> construct_from_bng_cplx_w_orientation(
      const BNG::BNGData& bng_data,
      const BNG::Cplx& bng_inst,
      const Orientation orientation);

  void set_name(const std::string& name_) override {
    BaseDataClass::set_name(name_);
    // rerun initialization because the name is parsed as a BNGL string
    elementary_molecules.clear();
    postprocess_in_ctor();
  }

  void postprocess_in_ctor() override;

  void check_semantics() const override {
    if (is_species_object()) {
      // all semantic checks will be done in Species
      return;
    }

    if (compartment_name == BNG::DEFAULT_COMPARTMENT_NAME) {
      throw ValueError("Compartment name '" + compartment_name + "' is reserved, please use a different name.");
    }

    GenComplex::check_semantics();

  }

  // using shorter printout when all_details is false
  std::string to_str(const bool all_details=false, const std::string ind="") const override;

  bool __eq__(const Complex& other) const override;

  std::string to_bngl_str() const override {
    return to_bngl_str_w_custom_orientation();
  }

  std::shared_ptr<Species> as_species() override;

  // do not export elementary_molecules (they are represented by name)
  bool skip_vectors_export() const override {
    return true;
  }

  // export into a single line
  bool export_as_string_without_newlines() const override {
    return true;
  }

  std::string export_to_python(std::ostream& out, PythonExportContext& ctx) {
    // we must set name for export if it was not set
    if (!is_set(name)) {
      name = to_bngl_str_w_custom_orientation();
    }
    return GenComplex::export_to_python(out, ctx);
  }

  // complexes can be only either surf or vol, there is no other option
  // WARNING: information on whether this is a surface or volume complex
  // may not be available all the time, e.g. when complex is constructed with m.Complex('A')
  bool is_vol() const {
    return !is_surf();
  }
  // WARNING: same as above
  bool is_surf() const;

  std::string to_bngl_str_w_custom_orientation(const bool include_mcell_orientation = false) const;

  // not really const, sets mutable members that serve as cache
  const std::string& get_canonical_name() const;

  // WARNING: information on whether this is a surface or volume complex - used when determining the
  // compartment may not be available all the time, e.g. when complex is constructed with m.Complex('A')
  const std::string& get_primary_compartment_name() const;


  void set_compartment_name(const std::string& new_compartment_name) {
    compartment_name = new_compartment_name;
    set_unset_compartments_of_elementary_molecules();
  }

private:
  void set_unset_compartments_of_elementary_molecules();
  bool is_species_object() const;

  // set when __eq__ is called, valid if cached_data_are_uptodate is true
  mutable std::string canonical_name;
};

} // namespace API
} // namespace MCell

#endif // API_COMPLEX_H
