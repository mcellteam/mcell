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

#ifndef API_ELEMENTARY_MOLECULE_TYPE_H
#define API_ELEMENTARY_MOLECULE_TYPE_H

#include "bng/bng_defines.h"

#include "generated/gen_elementary_molecule_type.h"
#include "api/api_common.h"
#include "api/elementary_molecule.h"

namespace BNG {
class BNGData;
class ElemMolType;
}

namespace MCell {
namespace API {

class ElementaryMoleculeType:
    public GenElementaryMoleculeType, public std::enable_shared_from_this<ElementaryMoleculeType> {
public:
  ELEMENTARY_MOLECULE_TYPE_CTOR()

  static std::shared_ptr<API::ElementaryMoleculeType> construct_from_bng_elem_mol_type(
    const BNG::BNGData& bng_data,
    const BNG::ElemMolType& bng_mt);

  void postprocess_in_ctor() override {
    mol_type_id = BNG::MOL_TYPE_ID_INVALID;
  }

  std::shared_ptr<ElementaryMolecule> inst(
      const std::vector<std::shared_ptr<Component>> components,
      const std::string& compartment_name = STR_UNSET) override {

    return std::make_shared<ElementaryMolecule>( shared_from_this(), components, compartment_name);
  }

  bool __eq__(const ElementaryMoleculeType& other) const override;

  std::string to_bngl_str() const override;

  bool skip_python_export() const override;

  // added methods
  std::string get_canonical_name() const;

  bool all_numerical_attributes_are_unset() const {
    // NOTE: this must be updated when attributtes are added
    return
      !is_set(diffusion_constant_2d) &&
      !is_set(diffusion_constant_3d) &&
      !is_set(custom_time_step) &&
      !is_set(custom_space_step) &&
      !target_only
    ;
  }

  bool warn_if_adding_identical_object() const { return false; }

  // mapping to MCell4
  BNG::elem_mol_type_id_t mol_type_id;
};

} // namespace API
} // namespace MCell

#endif // API_ELEMENTARY_MOLECULE_TYPE_H
