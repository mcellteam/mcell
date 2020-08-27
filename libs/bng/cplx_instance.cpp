/*
 * complex_species.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <iostream>
#include <sstream>
#include <map>

#include "bng/ast.h"
#include "bng/bng_engine.h"
#include "bng/cplx_instance.h"
#include "bng/mol_type.h"
#include "bng/mol_instance.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "debug_config.h"

using namespace boost;

using namespace std;

namespace BNG {

// ------------------------------------ CplxInstance -------------------------

void CplxInstance::finalize() {
  if (mol_instances.empty()) {
    return; // empty complex, ignoring finalization
  }

  // volume or surface type
  bool surf_type = false;
  bool reactive_surf_type = false;
  for (MolInstance& mp: mol_instances) {
    // need to finalize flags - copy them from molecule type
    mp.finalize_flags_and_sort_components(*bng_data);
    // if at least one is a surface molecule then the whole cplx is surface molecule
    if (mp.is_surf()) {
      surf_type = true;
    }

    // if at least one is a reactive surface then the whole cplx is reactive surface
    if (mp.is_reactive_surface()) {
      reactive_surf_type = true;
    }
  }
  if (surf_type) {
    release_assert(!reactive_surf_type && "Species cannot be both reactive surface and surface molecule.");
    set_flag(SPECIES_CPLX_MOL_FLAG_SURF);
  }
  else if (reactive_surf_type) {
    set_flag(SPECIES_CPLX_MOL_FLAG_REACTIVE_SURFACE);
  }

  // CPLX_FLAG_SINGLE_MOL_NO_COMPONENTS
  bool is_simple = true;
  if (mol_instances.size() > 1) {
    is_simple = false;
  }
  if (is_simple) {
    for (MolInstance& mp: mol_instances) {
      if (!mp.component_instances.empty()) {
        is_simple = false;
        break;
      }
    }
  }
  set_flag(SPECIES_CPLX_FLAG_ONE_MOL_NO_COMPONENTS, is_simple);

  // we need graphs even for simple complexes because they can be used in reaction patterns
  graph.clear();
  create_graph();

  set_finalized();
}


void CplxInstance::create_graph() {
  graph.clear();

  // convert molecule instances and their bonds into the boost graph representation

  map<bond_value_t, vector<Graph::vertex_descriptor>> bonds_to_vertices_map;

  // add all molecules with their components and remember how they should be bound
  for (MolInstance& mi: mol_instances) {
    Graph::vertex_descriptor mol_desc = boost::add_vertex(MtVertexProperty(Node(&mi)), graph);

    for (ComponentInstance& ci: mi.component_instances) {
      // for patterns, only components that were explicitly listed are in component instances

      Graph::vertex_descriptor comp_desc = boost::add_vertex(MtVertexProperty(Node(&ci)), graph);

      // connect the component to its molecule
      boost::add_edge(mol_desc, comp_desc, graph);

      // and remember its bond
      if (ci.bond_has_numeric_value()) {
        bonds_to_vertices_map[ci.bond_value].push_back(comp_desc);
      }
    }
  }

  // connect components
  for (auto it: bonds_to_vertices_map) {
    release_assert(it.second.size() == 2 && "There must be exactly pairs of components connected with numbered bonds.");
    boost::add_edge(it.second[0], it.second[1], graph);
  }
}


bool CplxInstance::matches_complex_pattern_ignore_orientation(const CplxInstance& pattern) const {

#ifdef DEBUG_CPLX_MATCHING
  cout << "** matches_complex_pattern_ignore_orientation:\n";
  cout << "pattern:\n";
  pattern.dump(false); cout << "\n";
  pattern.dump(true); cout << "\n";
  dump_graph(pattern.graph, bng_data);
  cout << "this instance:\n";
  dump(false); cout << "\n";
  dump(true); cout << "\n";
  dump_graph(graph, bng_data);
#endif
  // this result cannot be cached because it might not be and applicable for other equivalent complexes,
  // we might need to impose some ordering on elementary molecules and then we can reuse the result when
  // creating products
  VertexMappingVector mappings;
  get_subgraph_isomorphism_mappings(pattern.graph, graph, true, mappings);
  assert((mappings.size() == 0 || mappings.size() == 1) && "We are searching only for the first match");

#ifdef DEBUG_CPLX_MATCHING
  cout << "** result: " << !mappings.empty() << "\n";
#endif
  // we need at least one match
  return !mappings.empty();
}


uint CplxInstance::get_pattern_num_matches(const CplxInstance& pattern) const {
  assert(is_finalized() && pattern.is_finalized());
  VertexMappingVector mappings;
  get_subgraph_isomorphism_mappings(pattern.graph, graph, false, mappings);
  return mappings.size();
}


bool CplxInstance::matches_complex_fully_ignore_orientation(const CplxInstance& pattern) const {
  if (graph.m_vertices.size() != pattern.graph.m_vertices.size()) {
    // we need full match
    return false;
  }

  VertexMappingVector mappings;
  get_subgraph_isomorphism_mappings(pattern.graph, graph, true, mappings);
  assert((mappings.size() == 0 || mappings.size()) == 1 && "We are searching only for the first match");

  return mappings.size() == 1 && mappings[0].size() == graph.m_vertices.size();
}


class CanonicalComponentComparator {
public:
  CanonicalComponentComparator(const MolType& mt_)
    : mt(mt_) {
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

    // 3) state
    return ci1.state_id < ci2.state_id;
  }

private:
  const MolType& mt;
};


static bool canonical_mol_instance_less(const MolInstance& mi1, const MolInstance& mi2) {
  // 1) id
  if (mi1.mol_type_id != mi2.mol_type_id) {
    return mi1.mol_type_id < mi2.mol_type_id;
  }

  // no let's go component by component, the canonicalization is used currently only
  // for species that have full set of components
  vector<uint> bonds1;
  vector<uint> bonds2;
  assert(mi1.component_instances.size() == mi2.component_instances.size());
  for (size_t i = 0; i < mi1.component_instances.size(); i++) {
    const ComponentInstance& ci1 = mi1.component_instances[i];
    const ComponentInstance& ci2 = mi2.component_instances[i];
    assert(ci1.component_type_id == ci2.component_type_id);

    // 2) state
    if (ci1.state_id != ci2.state_id) {
      return ci1.state_id < ci2.state_id;
    }

    if (ci1.bond_has_numeric_value()) {
      bonds1.push_back(ci1.bond_value);
    }
    if (ci2.bond_has_numeric_value()) {
      bonds2.push_back(ci2.bond_value);
    }
  }

  // NOTE: the bonds comparison does not probably give fully canonical variant,
  // should be improved but it is ok like this for now
  if (bonds1.size() != bonds2.size()) {
    return bonds1.size() < bonds2.size();
  }
  else {
    for (size_t i = 0; i < bonds1.size(); i++) {
      if (bonds1[i] != bonds2[i]) {
        return bonds1[i] < bonds2[i];
      }
    }
    return false; // same for our purposes
  }
}


void CplxInstance::sort_components_and_mols() {
  // we need to sort components first
  for (MolInstance& mi: mol_instances) {
    CanonicalComponentComparator cmp(bng_data->get_molecule_type(mi.mol_type_id));
    sort(mi.component_instances.begin(), mi.component_instances.end(), cmp);
  }

  // then molecules
  sort(mol_instances.begin(), mol_instances.end(), canonical_mol_instance_less);
}


void CplxInstance::renumber_bonds() {
  map<bond_value_t, bond_value_t> new_bond_values;

  for (MolInstance& mi: mol_instances) {
    for (ComponentInstance& ci: mi.component_instances) {
      if (ci.bond_has_numeric_value()) {
        auto it = new_bond_values.find(ci.bond_value);
        if (it == new_bond_values.end()) {
          // new bond value, numbering from 1
          bond_value_t new_value = new_bond_values.size() + 1;
          new_bond_values[ci.bond_value] = new_value;
          ci.bond_value = new_value;
        }
        else {
          // we have already seen this bond value, use the new value
          ci.bond_value = it->second;
        }
      }
    }
  }
}


void CplxInstance::canonicalize() {

  // first sorting to put molecules to their places
  sort_components_and_mols();

  // make sure that bonds indices increase from the left
  renumber_bonds();

  // sort again after renumbering
  sort_components_and_mols();

  // need to rebuild graph
  finalize();
}


std::string CplxInstance::to_str(const BNGData& bng_data, bool in_surf_reaction) const {
  stringstream ss;
  for (size_t i = 0; i < mol_instances.size(); i++) {
    ss << mol_instances[i].to_str(bng_data);

    if (i != mol_instances.size() - 1) {
      ss << ".";
    }
  }
  if (orientation == ORIENTATION_UP) {
    ss << "'";
  }
  else if (orientation == ORIENTATION_DOWN) {
    ss << ",";
  }
  else if (in_surf_reaction && orientation == ORIENTATION_NONE) {
    ss << ";";
  }
  return ss.str();
}


void CplxInstance::dump(const bool for_diff, const std::string ind) const {
  if (!for_diff) {
    cout << ind << to_str(*bng_data);
  }
  else {
    cout << ind << "orientation: " << orientation << "\n";
    cout << ind << "mol_instances:\n";
    for (size_t i = 0; i < mol_instances.size(); i++) {
      cout << ind << i << ":\n";
      mol_instances[i].dump(*bng_data, true, ind + "  ");
    }
  }
}

} /* namespace BNG */
