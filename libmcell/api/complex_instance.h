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

#ifndef API_COMPLEX_INSTANCE_H
#define API_COMPLEX_INSTANCE_H

#include "../generated/gen_complex_instance.h"
#include "../api/common.h"

namespace MCell {
namespace API {

class Species;

class ComplexInstance: public GenComplexInstance {
public:
  COMPLEX_INSTANCE_CTOR()

  static std::shared_ptr<API::ComplexInstance> make_shared_empty() {
    // to avoid ComplexInstance semantic check, we need to insert a dummy name
    // when creating the object
    auto res_cplx_inst = std::make_shared<API::ComplexInstance>("TMP_NAME");
    res_cplx_inst->name = STR_UNSET;
    return res_cplx_inst;
  }

  void check_semantics() const override {
    GenComplexInstance::check_semantics();
    if (get_num_set(name, elementary_molecule_instances) != 1) {
      throw ValueError(
          S("Exactly one of ") + NAME_NAME + " or " + NAME_ELEMENTARY_MOLECULE_INSTANCES +
          " must be set for " + NAME_CLASS_COMPLEX_INSTANCE + ".");
    }
  }

  std::string to_bngl_str() override;

  std::shared_ptr<Species> as_species() override;

  // complex instances can be only either surf or vol, there is no other option
  bool is_vol() const {
    return !is_surf();
  }
  bool is_surf() const;
};

} // namespace API
} // namespace MCell

#endif // API_COMPLEX_INSTANCE_H
