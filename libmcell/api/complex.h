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

    if (orientation != Orientation::NONE && is_set(compartment_name)) {
      throw ValueError(S(NAME_COMPARTMENT_NAME) + " must not be set when " + NAME_ORIENTATION +
          " is set to a value other than " + NAME_ENUM_ORIENTATION + "." + NAME_EV_NONE + ".");
    }

    // TODO: how can we check that the used molecule types were defined?
  }

  std::string to_bngl_str() override;

  std::shared_ptr<Species> as_species() override;

  // complex instances can be only either surf or vol, there is no other option
  bool is_vol() const {
    return !is_surf();
  }
  bool is_surf() const;

private:
  bool is_species_object() const;
};

} // namespace API
} // namespace MCell

#endif // API_COMPLEX_H
