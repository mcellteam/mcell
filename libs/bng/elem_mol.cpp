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
#include "bng/cplx.h"
#include "bng/bngl_names.h"
#include "bng/elem_mol_type.h"

using namespace std;

namespace BNG {

static bool state_bond_less(
    const BNGData& bng_data, const Component& ci1, const Component& ci2) {

  // bond
  if (ci1.bond_value != ci2.bond_value) {
    return ci1.bond_value < ci2.bond_value;
  }

  // state name
  if (ci1.state_is_set() && ci2.state_is_set()) {
    const string& sn1 = bng_data.get_state_name(ci1.state_id);
    const string& sn2 = bng_data.get_state_name(ci2.state_id);
    return sn1 < sn2;
  }
  else {
    return ci1.state_id < ci2.state_id;
  }
}

class CanonicalComponentComparator {
public:
  // need to pass bng data for state names
  // and also molecule type because we do not know to which mol type the components belong
  CanonicalComponentComparator(const BNGData& bng_data_, const ElemMolType& mt_)
    : bng_data(bng_data_), mt(mt_) {
  }

  bool operator()(const Component& ci1, const Component& ci2) {
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

    return state_bond_less(bng_data, ci1, ci2);
  }

private:
  const BNGData& bng_data;
  const ElemMolType& mt;
};


class ComponentNameComparator {
public:
  // need to pass bng data for component and state names
  ComponentNameComparator(const BNGData& bng_data_)
    : bng_data(bng_data_) {
  }

  bool operator()(const Component& ci1, const Component& ci2) {
    // component name
    if (ci1.component_type_id != ci2.component_type_id) {
      return bng_data.get_component_type(ci1.component_type_id).name <
          bng_data.get_component_type(ci2.component_type_id).name;
    }

    return state_bond_less(bng_data, ci1, ci2);
  }

private:
  const BNGData& bng_data;
};

// ------------- Component -------------

std::string Component::to_str(const BNGData& bng_data) const {
  std::string res;
  to_str(bng_data, res);
  return res;
}


void Component::to_str(const BNGData& bng_data, std::string& res) const {
  const ComponentType& ct = bng_data.get_component_type(component_type_id);
  res += ct.name;

  assert(state_id != STATE_ID_INVALID);
  if (state_id != STATE_ID_DONT_CARE) {
    res += "~" + bng_data.get_state_name(state_id);
  }

  assert(state_id != BOND_VALUE_INVALID);
  if (bond_value == BOND_VALUE_BOUND) {
    res += "!" + BOND_STR_BOUND;
  }
  else if (bond_value == BOND_VALUE_ANY) {
    res += "!" + BOND_STR_ANY;
  }
  else if (bond_value == BOND_VALUE_UNBOUND) {
    // nothing to print
  }
  else {
    res += "!"  + to_string(bond_value);
  }
}

void Component::dump(const BNGData& bng_data, const string& ind) const {
  cout << ind << to_str(bng_data) << "\n";
}


// ------------- MoleculeInstance -------------

bool ElemMol::is_fully_qualified(const BNGData& bng_data) const {
  multiset<component_type_id_t> this_inst_components;
  multiset<component_type_id_t> mol_type_components;

  for (const Component& ci: components) {
    component_type_id_t ct_id = ci.component_type_id;
    this_inst_components.insert(ct_id);
    // if a component has a state, it must be set
    if (!bng_data.get_component_type(ct_id).allowed_state_ids.empty() && !ci.state_is_set()) {
      return false;
    }
  }

  for (component_type_id_t ct_id: bng_data.get_elem_mol_type(elem_mol_type_id).component_type_ids) {
    mol_type_components.insert(ct_id);
  }

  // do we have the same components?
  return this_inst_components == mol_type_components;
}


void ElemMol::canonicalize(const BNGData& bng_data) {
  // sort components first
  CanonicalComponentComparator comp_cmp(bng_data, bng_data.get_elem_mol_type(elem_mol_type_id));
  sort(components.begin(), components.end(), comp_cmp);
}


void ElemMol::sort_components_by_name(const BNGData& bng_data) {
  ComponentNameComparator name_cmp(bng_data);
  sort(components.begin(), components.end(), name_cmp);
}

// TODO: component sorting overlaps with MolType::canonicalize
void ElemMol::finalize_flags_and_sort_components(const BNGData& bng_data) {

  const ElemMolType& mt = bng_data.get_elem_mol_type(elem_mol_type_id);

  // copy flags from molecule type
  set_flags(mt.get_flags());

  set_finalized();
  // flag about molecule type must be set
  assert(is_vol() || is_surf() || is_reactive_surface());

  // sort components according to the order in molecule instance so that
  // all molecule instances are either identical or when used as patterns, they
  // can be easily compared
  if (components.size() > 1) {

    vector<bool> used_component_instances;
    used_component_instances.resize(components.size(), false);
    small_vector<Component> new_component_instances;

    for (size_t template_index = 0; template_index < mt.component_type_ids.size(); template_index++) {
      // try to find component that matches the current template position
      size_t instance_index;
      for (instance_index = 0; instance_index < components.size(); instance_index++) {
        if (!used_component_instances[instance_index] &&
            components[instance_index].component_type_id == mt.component_type_ids[template_index]) {
          break;
        }
      }
      if (instance_index < components.size()) {
        // found
        used_component_instances[instance_index] = true;
        new_component_instances.push_back(components[instance_index]);
      }
    }

    // update the component instance array
    components = new_component_instances;
  }
}


// this is a method used when computing rxn rate multiplier,
// we transform this mol instance pattern in a similar one only
// each added component has any state and any bond specifier
void ElemMol::insert_missing_components_as_any_state_pattern(const BNGData& bng_data) {
  assert(is_finalized() && "Components must be sorted");

  const ElemMolType& mt = bng_data.get_elem_mol_type(elem_mol_type_id);

  size_t instance_index = 0;
  for (size_t template_index = 0; template_index < mt.component_type_ids.size(); template_index++) {
    component_type_id_t template_component_id = mt.component_type_ids[template_index];
    if (instance_index < components.size() &&
        components[instance_index].component_type_id != template_component_id) {
      // missing - need to insert
      Component ci = Component(template_component_id);
      ci.state_id = STATE_ID_DONT_CARE;
      ci.bond_value = BOND_VALUE_ANY;
      components.push_back(ci);
    }
    else {
      // ok, present
      instance_index++;
    }
  }

  // sort components
  finalize_flags_and_sort_components(bng_data);
}


bool ElemMol::has_bond() const {
  for (const auto& comp: components) {
    if (comp.bond_has_numeric_value()) {
      return true;
    }
  }
  return false;
}


std::string ElemMol::to_str(const BNGData& bng_data) const {
  std::string res;
  to_str(bng_data, res);
  return res;
}


void ElemMol::to_str(const BNGData& bng_data, std::string& res, const bool include_compartment) const {
  const ElemMolType& mt = bng_data.get_elem_mol_type(elem_mol_type_id);

  res += mt.name;
  if (!components.empty()) {
    res += "(";
  }

  bool first_component = true;
  for (size_t i = 0; i < components.size(); i++) {

    if (!first_component) {
      res += ",";
    }

    components[i].to_str(bng_data, res);

    first_component = false;
  }
  if (!components.empty()) {
    res += ")";
  }

  if (include_compartment) {
    if (is_in_out_compartment_id(compartment_id)) {
      res += "@" + compartment_id_to_str(compartment_id);
    }
    else if (compartment_id != COMPARTMENT_ID_NONE) {
      const string& compartment_name = bng_data.get_compartment(compartment_id).name;
      if (compartment_name != DEFAULT_COMPARTMENT_NAME) {
        res += "@" + bng_data.get_compartment(compartment_id).name;
      }
    }
  }
}


void ElemMol::dump(const BNGData& bng_data, const bool for_diff, const std::string ind) const {
  if (!for_diff) {
    cout << to_str(bng_data);
  }
  else {
    const ElemMolType& mt = bng_data.get_elem_mol_type(elem_mol_type_id);
    cout << ind << "mol_type_id: " << elem_mol_type_id << " (" << mt.name << ")\n";
    cout << ind << "flags: " << BaseSpeciesCplxMolFlag::to_str() << "\n";
    cout << ind << "compartment: " << compartment_id_to_str(compartment_id) << "\n";
    cout << ind << "components: \n";
    for (const Component& ci: components) {
      ci.dump(bng_data, ind + "  ");
    }
  }
}


} /* namespace BNG */
