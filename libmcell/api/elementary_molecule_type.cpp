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

#include "api/elementary_molecule_type.h"
#include "api/component_type.h"
#include "bng/bng.h"

using namespace std;

namespace MCell {
namespace API {

std::shared_ptr<API::ElementaryMoleculeType> ElementaryMoleculeType::construct_from_bng_elem_mol_type(
    const BNG::BNGData& bng_data, const BNG::ElemMolType& bng_mt) {

  const string& name = bng_mt.name;
  auto res_mt = make_shared<API::ElementaryMoleculeType>(name);

  // using the MCELL_* parameters for now, see if the diffusion rate was
  // specified in the model
  double D2 = FLT_INVALID; // init to silence compiler warning
  bool found2 = bng_data.get_parameter_value(BNG::MCELL_DIFFUSION_CONSTANT_2D_PREFIX + name, D2);
  double D3 = FLT_INVALID;
  bool found3 = bng_data.get_parameter_value(BNG::MCELL_DIFFUSION_CONSTANT_3D_PREFIX + name, D3);

  if (found2 && found3) {
    throw RuntimeError("Molecule type '" + name + "' has both 2d and 3d diffusion constants specified.");
  }
  // no need to check whether the diffusion constant was set, it may be set later and
  // model conversion checks that it was set

  if (found2) {
    res_mt->diffusion_constant_2d = D2;
  }
  else if (found3) {
    res_mt->diffusion_constant_3d = D3;
  }

  // process components
  for (BNG::component_type_id_t ct_id: bng_mt.component_type_ids) {
    const BNG::ComponentType bng_ct = bng_data.get_component_type(ct_id);

    auto ct = make_shared<API::ComponentType>(bng_ct.name);

    // and allowed states
    for (BNG::state_id_t state_id: bng_ct.allowed_state_ids) {
      ct->states.push_back(bng_data.get_state_name(state_id));
    }

    res_mt->components.push_back(ct);
  }

  return res_mt;
}


static std::string get_components_str(
    const std::vector<std::shared_ptr<ComponentType>>& components,
    const bool canonical = false
) {
  string res;
  if (!components.empty()) {
    res += "(";
    for (size_t i = 0; i < components.size(); i++) {
      if (!canonical) {
        res += components[i]->to_bngl_str();
      }
      else {
        res += components[i]->get_canonical_name();
      }
      if (i + 1 != components.size()) {
        res += ",";
      }
    }
    res += ")";
  }
  return res;
}


std::string ElementaryMoleculeType::get_canonical_name() const {
  std::vector<std::shared_ptr<ComponentType>> sorted;
  sorted = components;
  std::sort(sorted.begin(), sorted.end(),
      [](const std::shared_ptr<ComponentType>& a, const std::shared_ptr<ComponentType>& b) -> bool {
          return *a < *b;
      });
  return name + get_components_str(sorted);
}


bool ElementaryMoleculeType::__eq__(const ElementaryMoleculeType& other) const {

  if (!eq_nonarray_attributes(other)) {
    return false;
  }

  return get_canonical_name() == other.get_canonical_name();
}


std::string ElementaryMoleculeType::to_bngl_str() const {
  return name + get_components_str(components);
}


bool ElementaryMoleculeType::skip_python_export() const {
  return BNG::is_species_superclass(name);
}


} // namespace API
} // namespace MCell
