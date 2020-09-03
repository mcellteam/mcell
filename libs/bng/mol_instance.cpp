/*
 * complex_species.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <iostream>
#include <sstream>
#include <algorithm>

#include "bng/ast.h"
#include "bng/bng_engine.h"
#include "bng/cplx_instance.h"
#include "bng/mol_type.h"

using namespace std;

namespace BNG {

class CanonicalComponentComparator {
public:
  // need to pass bng data for state names
  // and also molecule type because we do not know to which mol type the components belong
  CanonicalComponentComparator(const BNGData& bng_data_, const MolType& mt_)
    : bng_data(bng_data_), mt(mt_) {
  }

  bool operator()(const ComponentInstance& ci1, const ComponentInstance& ci2) {
    // we must maintain the order of components if they have the same name
    if (ci1.component_type_id != ci2.component_type_id) {
      // different components - the ordering is given by the order in the molecule
      uint ci1_index = INDEX_INVALID;
      uint ci2_index = INDEX_INVALID;

      for (size_t i = 0; i < mt.component_type_ids.size(); i++) {
        if (ci1_index == INDEX_INVALID && ci1.component_type_id == mt.component_type_ids[i]) {
          ci1_index = 0;
        }
        if (ci2_index == INDEX_INVALID && ci2.component_type_id == mt.component_type_ids[i]) {
          ci2_index = 0;
        }
      }
      assert(ci1_index != INDEX_INVALID && ci2_index != INDEX_INVALID);

      return ci1_index < ci2_index;
    }

    // 2) bond
    if (ci1.bond_value != ci2.bond_value) {
      return ci1.bond_value < ci2.bond_value;
    }

    // 3) state name
    if (ci1.state_is_set() && ci2.state_is_set()) {
      const string& sn1 = bng_data.get_state_name(ci1.state_id);
      const string& sn2 = bng_data.get_state_name(ci2.state_id);
      return sn1 < sn2;
    }
    else {
      return ci1.state_id < ci2.state_id;
    }
  }

private:
  const BNGData& bng_data;
  const MolType& mt;
};


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
  if (bond_value == BOND_VALUE_BOUND) {
    ss << "!" + BOND_STR_BOUND;
  }
  else if (bond_value == BOND_VALUE_ANY) {
    ss << "!" + BOND_STR_ANY;
  }
  else if (bond_value == BOND_VALUE_UNBOUND) {
    // nothing to print
  }
  else {
    ss << "!"  << bond_value;
  }
  return ss.str();
}

void ComponentInstance::dump(const BNGData& bng_data, const string& ind) const {
  cout << ind << to_str(bng_data) << "\n";
}


// ------------- MoleculeInstance -------------

void MolInstance::canonicalize(const BNGData& bng_data) {
  // we need to sort components first
  CanonicalComponentComparator comp_cmp(bng_data, bng_data.get_molecule_type(mol_type_id));
  sort(component_instances.begin(), component_instances.end(), comp_cmp);
}


// TODO: component sorting overlaps with MolType::canonicalize
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

    vector<bool> used_component_instances;
    used_component_instances.resize(component_instances.size(), false);
    small_vector<ComponentInstance> new_component_instances;

    for (size_t template_index = 0; template_index < mt.component_type_ids.size(); template_index++) {
      // try to find component that matches the current template position
      size_t instance_index;
      for (instance_index = 0; instance_index < component_instances.size(); instance_index++) {
        if (!used_component_instances[instance_index] &&
            component_instances[instance_index].component_type_id == mt.component_type_ids[template_index]) {
          break;
        }
      }
      if (instance_index < component_instances.size()) {
        // found
        used_component_instances[instance_index] = true;
        new_component_instances.push_back(component_instances[instance_index]);
      }
    }

    // update the component instance array
    component_instances = new_component_instances;
  }
}


// this is a method used when computing rxn rate multiplier,
// we transform this mol instance pattern in a similar one only
// each added component has any state and any bond specifier
void MolInstance::insert_missing_components_as_any_state_pattern(const BNGData& bng_data) {
  assert(is_finalized() && "Components must be sorted");

  const MolType& mt = bng_data.get_molecule_type(mol_type_id);

  size_t instance_index = 0;
  for (size_t template_index = 0; template_index < mt.component_type_ids.size(); template_index++) {
    component_type_id_t template_component_id = mt.component_type_ids[template_index];
    if (instance_index < component_instances.size() &&
        component_instances[instance_index].component_type_id != template_component_id) {
      // missing - need to insert
      ComponentInstance ci = ComponentInstance(template_component_id);
      ci.state_id = STATE_ID_DONT_CARE;
      ci.bond_value = BOND_VALUE_ANY;
      component_instances.push_back(ci);
    }
    else {
      // ok, present
      instance_index++;
    }
  }

  // sort components
  finalize_flags_and_sort_components(bng_data);
}


std::string MolInstance::to_str(const BNGData& bng_data) const {
  stringstream ss;
  const MolType& mt = bng_data.get_molecule_type(mol_type_id);

  ss << mt.name;
  if (!component_instances.empty()) {
    ss << "(";
  }

  bool first_component = true;
  for (size_t i = 0; i < component_instances.size(); i++) {

    if (!first_component) {
      ss << ",";
    }

    ss << component_instances[i].to_str(bng_data);

    first_component = false;
  }
  if (!component_instances.empty()) {
    ss << ")";
  }
  return ss.str();
}


void MolInstance::dump(const BNGData& bng_data, const bool for_diff, const std::string ind) const {
  if (!for_diff) {
    cout << to_str(bng_data);
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
