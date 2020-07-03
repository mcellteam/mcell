/*
 * complex_species.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <iostream>
#include <sstream>

#include "bng/ast.h"
#include "bng/bng_engine.h"
#include "bng/cplx_instance.h"
#include "bng/mol_type.h"

using namespace std;

namespace BNG {

// ------------- ComponentInstance -------------
std::string ComponentInstance::to_str(const BNGData& bng_data) const {
  stringstream ss;

  const ComponentType& ct = bng_data.get_component_type(component_type_id);
  ss << ct.name;

  assert(state_id != STATE_ID_INVALID);
  if (state_id != STATE_ID_DONT_CARE) {
    ss << "~" << bng_data.get_state_name(state_id);
  }

  assert(state_id != BOND_VALUE_INVALID);
  if (bond_value == BOND_VALUE_ANY) {
    ss << "!" + BOND_STR_ANY;
  }
  else if (bond_value != BOND_VALUE_NO_BOND) {
    ss << "!"  << bond_value;
  }
  return ss.str();
}

void ComponentInstance::dump(const BNGData& bng_data, const string& ind) const {
  cout << ind << to_str(bng_data) << "\n";
}


// ------------- MoleculeInstance -------------
/*
void MolInstance::initialize_components_types(const MolType& mt) {
  for (component_type_id_t component_type_id: mt.component_type_ids) {
    // state is don't care, no bond
    component_instances.push_back(ComponentInstance(component_type_id));
  }
}
*/

/*
// searches for component with name
uint MolInstance::get_corresponding_component_index(
    const BNGData& bng_data,
    const MolType& mt,
    const std::string& name,
    const uint starting_index
) const {
  for (uint i = starting_index; i < mt.component_type_ids.size(); i++) {
    const ComponentType& ct = bng_data.get_component_type(mt.component_type_ids[i]);

    // we are looking for the first component with the same name
    if (ct.name == name) {
      return i;
    }
  }

  return INDEX_INVALID;
}
*/

void MolInstance::finalize_flags_and_sort_components(const BNGData& bng_data) {

  const MolType& mt = bng_data.get_molecule_type(mol_type_id);

  // copy flags from molecule type
  set_flags(mt.get_flags());

  set_finalized();
  // flag about molecule type must be set
  assert(is_vol() || is_surf() || is_reactive_surface());

  // sort components according to the order in molecule instance so that
  // all molecule instances are either identical or when used as patterns, they
  // can be easily compared
  if (component_instances.size() > 1) {

    // last item is skipped because it will be in the right position already
    size_t template_pos = 0;
    for (size_t pos_to_assign = 0; pos_to_assign < component_instances.size() - 1; pos_to_assign++) {

      assert(template_pos < mt.component_type_ids.size());

      // use template to find which item to swap
      size_t pos_to_swap = pos_to_assign;
      bool found = false;
      do {
        while (
            pos_to_swap < component_instances.size() &&
            mt.component_type_ids[template_pos] != component_instances[pos_to_swap].component_type_id
        ) {
          pos_to_swap++;
        }
        found = pos_to_swap != component_instances.size();
        template_pos++;
      } while (!found);

      if (pos_to_assign != pos_to_swap) {
        ComponentInstance tmp = component_instances[pos_to_assign];
        component_instances[pos_to_assign] = component_instances[pos_to_swap];
        component_instances[pos_to_swap] = tmp;
      }
    }
  }
}


std::string MolInstance::to_str(const BNGData& bng_data, const bool only_explicit) const {
  stringstream ss;
  const MolType& mt = bng_data.get_molecule_type(mol_type_id);

  ss << mt.name;
  if (!component_instances.empty()) {
    ss << "(";
  }

  bool first_component = true;
  for (size_t i = 0; i < component_instances.size(); i++) {

    if (!only_explicit || component_instances[i].explicitly_listed_in_pattern) {
      if (!first_component) {
        ss << ",";
      }

      ss << component_instances[i].to_str(bng_data);

      first_component = false;
    }
  }
  if (!component_instances.empty()) {
    ss << ")";
  }
  return ss.str();
}


void MolInstance::dump(const BNGData& bng_data, const bool for_diff, const bool only_explicit, const std::string ind) const {
  if (!for_diff) {
    cout << to_str(bng_data, only_explicit);
  }
  else {
    const MolType& mt = bng_data.get_molecule_type(mol_type_id);
    cout << ind << "mol_type_id: " << mol_type_id << " (" << mt.name << ")\n";
    cout << ind << "flags: " << BaseSpeciesCplxMolFlag::to_str() << "\n";
    cout << ind << "components: \n";
    for (const ComponentInstance& ci: component_instances) {
      ci.dump(bng_data, ind + "  ");
    }
  }
}


} /* namespace BNG */
