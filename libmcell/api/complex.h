/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#ifndef API_COMPLEX_H
#define API_COMPLEX_H

#include "generated/gen_complex.h"
#include "api/common.h"
#include "api/api_utils.h"
#include "api/compartment_utils.h"

namespace MCell {
namespace API {

class Species;

class Complex: public GenComplex {
public:
  COMPLEX_CTOR()

  static std::shared_ptr<API::Complex> make_shared_empty() {
    // to avoid ComplexInstance semantic check, we need to insert a dummy name
    // when creating the object
    auto res_cplx_inst = std::make_shared<API::Complex>("TMP_NAME");
    res_cplx_inst->name = STR_UNSET;
    return res_cplx_inst;
  }

  void postprocess_in_ctor() override {
    // set compartment_name and check that there is only one
    if (is_set(name)) {
      std::vector<std::string> compartments;
      // we do not want to use the BNG parser at thsi point
      // or do we?
      get_compartment_names(name, compartments);
      if (!compartments.empty()) {
        if (is_set(compartment_name)) {
          throw ValueError("Complex " + name + " is defined with both compartment name in " +
              NAME_NAME + " and in " + NAME_COMPARTMENT_NAME + ", only one is allowed.");
        }

        std::string single_compartment_name = compartments[0];
        for (size_t i = 1; i < compartments.size(); i++) {
          if (single_compartment_name != compartments[i]) {
            throw ValueError("Complex cannot be in multiple compartments, error for " + name + ".");
          }
        }
        compartment_name = single_compartment_name;
      }
    }
  }

  void check_semantics() const override {
    if (is_species_object()) {
      // all semantic checks will be done in Species
      return;
    }
    GenComplex::check_semantics();
    if (get_num_set(name, elementary_molecule_instances) != 1) {
      throw ValueError(
          S("Exactly one of ") + NAME_NAME + " or " + NAME_ELEMENTARY_MOLECULE_INSTANCES +
          " must be set for " + NAME_CLASS_COMPLEX + ".");
    }

    if (is_set(name)) {
      if (is_simple_species(name) && name.find('.') != std::string::npos) {
        throw ValueError("Simple species name must not contain '.', this is incompatible with BNGL definition"
            ", error for " + name + ".");
      }
    }
  }

  bool __eq__(const Complex& other) const override;

  std::string to_bngl_str() const override;

  std::shared_ptr<Species> as_species() override;

  // complex instances can be only either surf or vol, there is no other option
  bool is_vol() const {
    return !is_surf();
  }
  bool is_surf() const;

private:
  bool is_species_object() const;

  // not really const, sets mutable members
  void set_canonical_name_if_needed() const;

  // set when __eq__ is called, valid if cached_data_are_uptodate is true
  mutable std::string canonical_name;
};

} // namespace API
} // namespace MCell

#endif // API_COMPLEX_H
