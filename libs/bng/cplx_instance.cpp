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

using namespace boost;

using namespace std;

namespace BNG {

// ------------------------------------ CplxInstance -------------------------

void CplxInstance::finalize() {
  assert(!mol_instances.empty() && "There must be at least one molecule type");

  // finalize mol instances first
  for (MolInstance& mp: mol_instances) {
    mp.finalize_flags();
  }

  // volume or surface type
  bool vol_type = true;
  for (MolInstance& mp: mol_instances) {
    mp.finalize_flags();
    if (mp.is_surf()) {
      vol_type = false;
    }
  }
  if (!vol_type) {
    set_flag(SPECIES_CPLX_MOL_FLAG_SURF);
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
  create_graph();

  set_finalized();
}


void CplxInstance::create_graph() {
  // convert molecule instances and their bonds into the boost graph representation

  map<bond_value_t, vector<Graph::vertex_descriptor>> bonds_to_vertices_map;

  // add all molecules with their components and remomber how they should be bound
  for (const MolInstance& mi: mol_instances) {
    Graph::vertex_descriptor mol_desc = boost::add_vertex(MtVertexProperty(Node(&mi)), graph);

    for (const ComponentInstance& ci: mi.component_instances) {
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
    release_assert(it.second.size() != 2 && "There must be exactly pairs of components connected with numbered bonds.");
    boost::add_edge(it.second[0], it.second[1], graph);
  }
}



bool CplxInstance::matches_complex_pattern_ignore_orientation(const CplxInstance& pattern) const {
  assert(false && "Support for BNG style matching is not implemented yet");
  return false;
}


bool CplxInstance::matches_complex_fully_ignore_orientation(const CplxInstance& pattern) const {
  if (graph.m_vertices.size() != pattern.graph.m_vertices.size()) {
    // we need full match
    return false;
  }

  MultipleMappingsVector res;
  get_subgraph_isomorphism_mappings(pattern.graph, graph, res);

  // at least one mapping found and all nodes are covered
  return res.size() > 1 && res[0].size() == graph.m_vertices.size();
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

void CplxInstance::dump(const BNGData& bng_data, const bool for_diff, const std::string ind) const {
  if (!for_diff) {
    cout << ind << to_str(bng_data);
  }
  else {
    cout << ind << "orientation: " << orientation << "\n";
    cout << ind << "mol_instances:\n";
    for (size_t i = 0; i < mol_instances.size(); i++) {
      cout << ind << i << ":\n";
      mol_instances[i].dump(bng_data, true, false, ind + "  ");
    }
  }
}

} /* namespace BNG */
