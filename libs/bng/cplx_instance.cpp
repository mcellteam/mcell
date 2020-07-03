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
      mol_instances[i].dump(*bng_data, true, false, ind + "  ");
    }
  }
}

} /* namespace BNG */
