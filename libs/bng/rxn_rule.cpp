/*
 * rule.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#include <iostream>
#include <sstream>


#include "bng/graph.h" // must be included before boost

#include <boost/graph/copy.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>

#include "bng/rxn_rule.h"
#include "bng/rxn_class.h"

#include "bng/species.h"
#include "bng/species_container.h"

#include "debug_config.h"


using namespace std;

namespace BNG {


// used in set to represent edges
class UnorderedPair {
public:
  UnorderedPair(const vertex_descriptor_t a, const vertex_descriptor_t b) : first(std::min(a,b)), second(std::max(a,b)) {
  }
  bool operator == (const UnorderedPair& b) const {
    return first == b.first && second == b.second;
  }
  bool operator < (const UnorderedPair& b) const {
    return first < b.first && second < b.second;
  }

  vertex_descriptor_t first;
  vertex_descriptor_t second;
};



void RxnRule::finalize() {
  assert(id != RXN_RULE_ID_INVALID);

  bool simple = true;

  // finalize all reactants and products
  for (CplxInstance& ci: reactants) {
    ci.finalize();
    simple &= ci.is_simple();
  }

  num_surf_products = 0;
  for (CplxInstance& ci: products) {
    ci.finalize();
    simple &= ci.is_simple();
    if (ci.is_surf()) {
      num_surf_products++;
    }
  }

  if (simple) {
    set_flag(RXN_FLAG_SIMPLE);
  }

  compute_reactants_products_mapping();

  // for MCell3 compatibility
  move_products_that_are_also_reactants_to_be_the_first_products();

  set_finalized();
}


static void merge_graphs(Graph& srcdst, const Graph& src) {
  boost::copy_graph(src, srcdst);
}

// TODO: use BGL_FORALL_VERTICES_T

static int get_rxn_component_instance_matching_score(
    const Node& prod_comp,
    const Node& reac_comp
) {
  assert(!prod_comp.is_mol);
  assert(!reac_comp.is_mol);

  const ComponentInstance& ci1 = *prod_comp.component;
  const ComponentInstance& ci2 = *reac_comp.component;

  if (ci1.component_type_id != ci2.component_type_id) {
    return -1; // not a match
  }

  int res = 0;
  if (ci1.state_id == ci2.state_id) {
    res++;
  }

  if (ci1.bond_value == ci2.bond_value) {
    res++;
  }

  // note: maybe we will need to make the scoring more fine grained
  return res;
}


static int get_rxn_mol_instance_matching_score(
    Graph& products_graph,
    const vertex_descriptor_t prod_mi_desc,
    Graph& reactants_graph,
    const vertex_descriptor_t reac_mi_desc,
    const VertexNameMap& index
) {
  // similar code as in find_best_component_mapping
  // FIXME: find a way how to merge it

  int res = 0;

  set<vertex_descriptor_t> mapped_components;

  boost::graph_traits<Graph>::out_edge_iterator prod_ei, prod_edge_end;
  for (boost::tie(prod_ei, prod_edge_end) = boost::out_edges(prod_mi_desc, products_graph); prod_ei != prod_edge_end; ++prod_ei) {
    vertex_descriptor_t prod_comp_desc = boost::target(*prod_ei, products_graph);
    const Node& prod_comp = index[prod_comp_desc];

    int best_comp_score = -1; // not found
    vertex_descriptor_t best_comp_reac_desc;

    // compare against components of the best match
    boost::graph_traits<Graph>::out_edge_iterator reac_ei, reac_edge_end;
    for (boost::tie(reac_ei, reac_edge_end) = boost::out_edges(reac_mi_desc, reactants_graph); reac_ei != reac_edge_end; ++reac_ei) {
      vertex_descriptor_t reac_comp_desc = boost::target(*reac_ei, products_graph);
      const Node& reac_comp = index[reac_comp_desc];

      if (!reac_comp.is_mol &&
          reac_comp.component->component_type_id == prod_comp.component->component_type_id &&
          mapped_components.count(reac_comp_desc) == 0 // not mapped yet
      ) {

        int current_score = get_rxn_component_instance_matching_score(
            prod_comp,
            reac_comp
        );
        if (current_score > best_comp_score) {
          best_comp_reac_desc = reac_comp_desc;
          best_comp_score = current_score;
        }
      }
    }

    if (best_comp_score >= 0) {
      // remember score of this mapping
      res += best_comp_score;

      // we used this reactant's component
      assert(mapped_components.count(best_comp_reac_desc) == 0);
      mapped_components.insert(best_comp_reac_desc);
    }
  }

  return res;
}


static void find_best_component_mapping(
    Graph& products_graph,
    vertex_descriptor_t prod_mi_desc,
    Graph& reactants_graph,
    vertex_descriptor_t reac_mi_desc,
    const VertexNameMap& index,
    VertexMapping& prod_reac_mapping
) {
  set<vertex_descriptor_t> mapped_components;

  // now also define mapping for the components of this molecule instance
  boost::graph_traits<Graph>::out_edge_iterator prod_ei, prod_edge_end;
  for (boost::tie(prod_ei, prod_edge_end) = boost::out_edges(prod_mi_desc, products_graph); prod_ei != prod_edge_end; ++prod_ei) {
    vertex_descriptor_t prod_comp_desc = boost::target(*prod_ei, products_graph);
    const Node& prod_comp = index[prod_comp_desc];

    int best_comp_score = -1; // not found
    vertex_descriptor_t best_comp_reac_desc;

    // compare against components of the best match
    boost::graph_traits<Graph>::out_edge_iterator reac_ei, reac_edge_end;
    for (boost::tie(reac_ei, reac_edge_end) = boost::out_edges(reac_mi_desc, reactants_graph); reac_ei != reac_edge_end; ++reac_ei) {
      vertex_descriptor_t reac_comp_desc = boost::target(*reac_ei, products_graph);
      const Node& reac_comp = index[reac_comp_desc];

      if (!reac_comp.is_mol &&
          reac_comp.component->component_type_id == prod_comp.component->component_type_id &&
          mapped_components.count(reac_comp_desc) == 0 // not mapped yet
      ) {

        int current_score = get_rxn_component_instance_matching_score(
            prod_comp,
            reac_comp
        );
        if (current_score > best_comp_score) {
          best_comp_reac_desc = reac_comp_desc;
          best_comp_score = current_score;
        }
      }
    }

    if (best_comp_score >= 0) {
      // remember this mapping
      assert(prod_reac_mapping.count(prod_comp_desc) == 0);
      prod_reac_mapping[prod_comp_desc] = best_comp_reac_desc;

      // we mapped this component reactant
      assert(mapped_components.count(best_comp_reac_desc) == 0);
      mapped_components.insert(best_comp_reac_desc);
    }
  }
}

static void find_best_product_to_pattern_mapping(
    Graph& products_graph,
    Graph& reactants_graph,
    VertexMapping& prod_reac_mapping) {

  prod_reac_mapping.clear();

  set<vertex_descriptor_t> mapped_reactants;

  // get the property map for vertex indices
  const VertexNameMap& index = boost::get(boost::vertex_name, reactants_graph);

  // for each molecule instance in pattern_graph
  typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> prod_it;
  for (prod_it = boost::vertices(products_graph); prod_it.first != prod_it.second; ++prod_it.first) {
    vertex_descriptor_t prod_desc = *prod_it.first;
    const Node& prod_mi = index[prod_desc];
    if (prod_mi.is_mol) {

      // find the best match for this molecule instance
      // in products graph
      int best_mol_score = -1; // not found
      vertex_descriptor_t best_mol_reac_desc;
      std::pair<vertex_iter, vertex_iter> reac_it;
      for (reac_it = boost::vertices(reactants_graph); reac_it.first != reac_it.second; ++reac_it.first) {
        const Node& reac_mi = index[*reac_it.first];
        if (reac_mi.is_mol &&
            reac_mi.mol->mol_type_id == prod_mi.mol->mol_type_id &&
            mapped_reactants.count(*reac_it.first) == 0 // not mapped yet
        ) {

          int current_score = get_rxn_mol_instance_matching_score(
              products_graph,
              *reac_it.first,
              reactants_graph,
              *prod_it.first,
              index
          );

          if (current_score > best_mol_score) {
            best_mol_reac_desc = *reac_it.first;
            best_mol_score = current_score;
          }
        }
      }

      if (best_mol_score >= 0) {
        // remember this mapping
        assert(prod_reac_mapping.count(prod_desc) == 0);
        prod_reac_mapping[prod_desc] = best_mol_reac_desc;

        // we mapped this reactant
        assert(mapped_reactants.count(best_mol_reac_desc) == 0);
        mapped_reactants.insert(best_mol_reac_desc);

        find_best_component_mapping(
            products_graph,
            prod_desc,
            reactants_graph,
            best_mol_reac_desc,
            index,
            prod_reac_mapping
        );
      } // ^^^ handling of found reactants
    } // is mol
  }
}


// vertex pointed to by reac_desc
// must be a component and must have a single bond
// to another component and also a bond to its molecule instance
// returns -1 when the component is not connected to another component,
// asserts in debug mode when must_exist is set to true and target was not found
const vertex_descriptor_t TARGET_NOT_FOUND = -1;
static vertex_descriptor_t get_bond_target(
    const Graph& graph,
    const VertexNameMap& index,
    const vertex_descriptor_t desc,
    const bool must_exist = true
) {

  bool comp_found = false;
  bool mol_found = false;
  vertex_descriptor_t res = TARGET_NOT_FOUND;

  boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
  for (boost::tie(ei,edge_end) = boost::out_edges(desc, graph); ei != edge_end; ++ei) {
    Graph::edge_descriptor e_desc = *ei;
    Graph::vertex_descriptor n_desc = boost::target(e_desc, graph);
    const Node& n = index[n_desc];

    if (n.is_mol) {
      assert(!mol_found);
      mol_found = true; // just for debug
    }
    else {
      assert(!comp_found);
      comp_found = true;
      res = n_desc;
    }
  }

  assert(mol_found && "Component must be connected to its molecule");
  if (must_exist) {
    assert(comp_found);
  }

  return res;
}


vertex_descriptor_t get_new_bond_target(
    const Graph& reactants_graph,
    const VertexMapping& pattern_reactant_mapping,
    const VertexMapping& prod_pattern_mapping,
    const Graph& products_graph,
    const VertexNameMap& products_index,
    const vertex_descriptor_t prod_desc
) {

  // to which node we should connect?
  vertex_descriptor_t prog_target_desc = get_bond_target(products_graph, products_index, prod_desc);

  // to which reactant we will point?
  auto target_prod_pat_it = prod_pattern_mapping.find(prog_target_desc);
  assert(target_prod_pat_it != prod_pattern_mapping.end() && "Mapping must exist");
  auto target_pat_reac_it = pattern_reactant_mapping.find(target_prod_pat_it->second);
  assert(target_pat_reac_it != pattern_reactant_mapping.end() && "Mapping must exist");

  return target_pat_reac_it->second;
}


static void apply_rxn_on_reactants_graph(
    Graph& reactants_graph,
    const VertexMapping& pattern_reactant_mapping,
    const Graph& pattern_graph,
    const VertexMapping& prod_pattern_mapping,
    Graph& products_graph
) {
  set<vertex_descriptor_t> mol_instances_to_keep;

  VertexNameMap reactants_index = boost::get(boost::vertex_name, reactants_graph);
  VertexNameMap products_index = boost::get(boost::vertex_name, products_graph);

  // TODO: needed when molecules are removed
  // remove all disconnected graphs whose molecule instance is not mapped to
  // by products


  set<UnorderedPair> bonds_to_remove;
  set<UnorderedPair> bonds_to_add;

  //
  // for each component in product graph:
  //   find corresponding component in reactant pattern
  //     if there is no such mapping, ignore it because it might be a component of a
  //     new molecule instance
  for (auto prod_pat_it: prod_pattern_mapping) {
    vertex_descriptor_t prod_desc = prod_pat_it.first;
    const Node& prod_comp = products_index[prod_pat_it.first];
    if (!prod_comp.is_mol) {
      // product->pattern
      vertex_descriptor_t pat_desc = prod_pat_it.second;

      //   find component corresponding to reactant pattern in reactants_graph
      //     this mapping must exist because we matched the pattern graph to reactants graph
      auto pat_reac_it = pattern_reactant_mapping.find(pat_desc);
      assert(pat_reac_it != pattern_reactant_mapping.end() && "Mapping should exist?");

      vertex_descriptor_t reac_desc = pat_reac_it->second;

      //   manipulate state and or bond
      Node& reac_comp = reactants_index[reac_desc];
      assert(!reac_comp.is_mol);
      const ComponentInstance& prod_ci = *prod_comp.component;
      ComponentInstance& reac_ci = *reac_comp.component;

      // update state
      if (prod_ci.state_is_set() && prod_ci.state_id != reac_ci.state_id) {
        reac_ci.state_id = prod_ci.state_id;
      }

      // and bond
      // TODO: not handling !? yet
      if (prod_ci.bond_value != reac_ci.bond_value) {
        // orig: !+
        if (reac_ci.bond_value == BOND_VALUE_ANY) {
          // new: (no bond)
          if (prod_ci.bond_value == BOND_VALUE_NO_BOND) {
            bonds_to_remove.insert(UnorderedPair(
                reac_desc,
                get_bond_target(reactants_graph, reactants_index, reac_desc)
            ));
          }
          // new: !1
          else if (prod_ci.bond_has_numeric_value()) {
            assert(false && "Cannot change bond from !+ to !1");
          }
          else {
            assert(false); // no other option
          }
        }
        // orig: (no bond)
        else if (reac_ci.bond_value == BOND_VALUE_NO_BOND) {
          // new: !1
          if (prod_ci.bond_has_numeric_value()) {

            vertex_descriptor_t target_reac_desc = get_new_bond_target(
                reactants_graph,
                pattern_reactant_mapping,
                prod_pattern_mapping,
                products_graph,
                products_index,
                prod_desc
            );

            bonds_to_add.insert(UnorderedPair(reac_desc, target_reac_desc));
          }
          // new: !+
          else if (prod_ci.bond_value == BOND_VALUE_ANY){
            assert(false && "Cannot change bond from to !+");
          }
          else {
            assert(false);
          }
        }
        // orig: !1
        else if (reac_ci.bond_has_numeric_value()) {
          // new: (no bond)
          if (prod_ci.bond_value == BOND_VALUE_NO_BOND) {
            bonds_to_remove.insert(UnorderedPair(
                reac_desc,
                get_bond_target(reactants_graph, reactants_index, reac_desc)
            ));
          }
          // new: !2
          else if (prod_ci.bond_has_numeric_value()) {
            assert(prod_ci.bond_value != reac_ci.bond_value);
            // remove original one
            bonds_to_remove.insert(UnorderedPair(
                reac_desc,
                get_bond_target(reactants_graph, reactants_index, reac_desc)
            ));

            // and create a new one
            vertex_descriptor_t target_reac_desc = get_new_bond_target(
                reactants_graph,
                pattern_reactant_mapping,
                prod_pattern_mapping,
                products_graph,
                products_index,
                prod_desc
            );
            bonds_to_add.insert(UnorderedPair(reac_desc, target_reac_desc));
          }
          // new: !+
          else if (prod_ci.bond_value == BOND_VALUE_ANY){
            assert(false && "Cannot change bond from !1 to !+");
          }
          else {
            assert(false);
          }
        }
      } // if (prod_ci.bond_value != reac_ci.bond_value)
    } // if (!prod_comp.is_mol)
  } // for (auto prod_pat_it: prod_pattern_mapping)

  // remove bonds
  for (auto p: bonds_to_remove) {
    boost::remove_edge(p.first, p.second, reactants_graph);
  }

  // add bonds
  for (auto p: bonds_to_add) {
    boost::add_edge(p.first, p.second, reactants_graph);
  }

  // TODO: needed when molecules are added
  //
  // for each molecule in product graph:
  //   if there is no mapping to reactant pattern graph:
  //     check that all components that have state are connected and defined in the product graph
  //     create a new molecule instance in the reactants_graph
  //     for each component of the molecule in product graph:
  //       add it and connect to molecule
  //       if it has a bond:
  //         if does the second target already exist:
  //           create a bond
  //
}


static void convert_graph_component_to_cplx_inst(
    Graph& graph,
    const vector<int>& graph_components,
    const int graph_component_index,
    CplxInstance& cplx
) {
  VertexNameMap index = boost::get(boost::vertex_name, graph);

  map<UnorderedPair, int> bonds;

  // for each molecule instance
  typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> mol_it;
  for (mol_it = boost::vertices(graph); mol_it.first != mol_it.second; ++mol_it.first) {
    Graph::vertex_descriptor mol_desc = *mol_it.first;

    if (graph_components[mol_desc] != graph_component_index) {
      continue;
    }

    const Node& mol = index[mol_desc];
    if (!mol.is_mol) {
      continue;
    }

    cplx.mol_instances.push_back(*mol.mol);
    MolInstance& mi = cplx.mol_instances.back();

    // we will recreate components because they might have changed
    mi.component_instances.clear();

    // for each of its components
    boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
    for (boost::tie(ei,edge_end) = boost::out_edges(mol_desc, graph); ei != edge_end; ++ei) {
      Graph::edge_descriptor e_mol_comp = *ei;
      Graph::vertex_descriptor comp_desc = boost::target(e_mol_comp, graph);

      const Node& comp = index[comp_desc];
      assert(!comp.is_mol && "Only a component may be connected to a molecule.");

      mi.component_instances.push_back(*comp.component); // we use state as it its
      ComponentInstance& compi = mi.component_instances.back();

      // we need to set bonds
      Graph::vertex_descriptor bound_comp_desc = get_bond_target(graph, index, comp_desc, false);
      if (bound_comp_desc == TARGET_NOT_FOUND) {
        compi.bond_value = BOND_VALUE_NO_BOND;
      }
      else {
        UnorderedPair bond(comp_desc, bound_comp_desc);
        auto it = bonds.find(bond);
        if (it != bonds.end()) {
          compi.bond_value = it->second;
        }
        else {
          int next_bond_index = bonds.size() + 1; // we start from 1
          bonds[bond] = next_bond_index;
          compi.bond_value = next_bond_index;
        }
      }
    }

  } // for each vertex


  cplx.finalize();
}


static void create_products_from_reactants_graph(
    const BNGData* bng_data,
    Graph& reactants_graph,
    vector<CplxInstance>& created_products
) {
  created_products.clear();

  // the output of the algorithm is recorded in the component property map comp,
  // which will contain numbers giving the component number assigned to each vertex
  // the resulting vector is indexed by vertex_descriptor
  vector<int> comp(boost::num_vertices(reactants_graph));
  int num_components = boost::connected_components(
      reactants_graph,
      boost::make_iterator_property_map(
          comp.begin(), boost::get(boost::vertex_index, reactants_graph), comp[0]
      )
  );

  for (int i = 0; i < num_components; i++) {
    created_products.push_back(CplxInstance(bng_data));
    CplxInstance& new_cplx = created_products.back();
    convert_graph_component_to_cplx_inst(reactants_graph, comp, i, new_cplx);
  }
}

void RxnRule::create_products_for_complex_rxn(
    const vector<const CplxInstance*>& input_reactants,
    vector<CplxInstance>& created_products
) const {

  assert(mol_instances_are_fully_maintained && "Assuming this for now");

  assert(input_reactants.size() == reactants.size());
  assert(input_reactants.size() == 1 || input_reactants.size() == 2);

  // merge input reactant graphs
  Graph reactants_graph = input_reactants[0]->get_graph();
  if (input_reactants.size() == 2) {
    merge_graphs(reactants_graph, input_reactants[1]->get_graph());
  }

  // merge reactant patterns
  // this can be precomputed
  Graph pattern_graph = reactants[0].get_graph();

#ifdef DEBUG_CPLX_MATCHING
  cout << "Base pattern:\n";
  dump_graph(reactants[0].get_graph());
#endif

  if (reactants.size() == 2) {
    #ifdef DEBUG_CPLX_MATCHING
      cout << "Merging with:\n";
      dump_graph(reactants[1].get_graph());
    #endif
    merge_graphs(pattern_graph, reactants[1].get_graph());
  }

  // compute mapping reactant pattern -> reactant
#ifdef DEBUG_CPLX_MATCHING
  cout << "Pattern:\n";
  dump_graph(pattern_graph);

  cout << "Reactants:\n";
  dump_graph(reactants_graph);
#endif

  VertexMappingVector pattern_reactant_mappings;
  get_subgraph_isomorphism_mappings(
      pattern_graph, // pattern
      reactants_graph, // actual reactant
      false, // do not stop with first match
      pattern_reactant_mappings
  );
  assert(pattern_reactant_mappings.size() != 0 &&
      "Did not find a match of patterns onto reaction.");
  assert(pattern_reactant_mappings.size() == 1 &&
      "We do not support multiple matches yet and then there must be at least one match.");


  // create graph from products
  // this can be precomputed
  Graph products_graph;
  if (products.size() > 0) {
    products_graph = products[0].get_graph();
  }
  for (size_t i = 1; i < products.size(); i++) {
    merge_graphs(products_graph, products[i].get_graph());
  }
#ifdef DEBUG_CPLX_MATCHING
  cout << "Products:\n";
  dump_graph(products_graph);
#endif

  // compute mapping reactant patterns -> product patterns
  // this can be precomputed
  // nodes in the reaction products graph have an extra option that they can ignore the state
  // find all mappings and compute score for each of those mappings so that we get the best match
  // -> each matching state is a positive point

  VertexMapping product_pattern_mapping;
  // TODO: boost subgraph won't find a mapping of one of the reactants is not present in products,
  // we will need to write this manually, hopefully the graphs won't be too large
  // still must avoid exponential complexity...
  find_best_product_to_pattern_mapping(products_graph, reactants_graph, product_pattern_mapping);

#ifdef DEBUG_CPLX_MATCHING
  dump_graph_mapping(product_pattern_mapping);
#endif

  // manipulate nodes using information about products
  apply_rxn_on_reactants_graph(
      reactants_graph,
      pattern_reactant_mappings[0],
      pattern_graph,
      product_pattern_mapping,
      products_graph
  );

#ifdef DEBUG_CPLX_MATCHING
  cout << "\nReactants after applying rxn:\n";
  dump_graph(products_graph);
#endif

  // and finally create products, each disconnected graph in the result is a
  // separate complex instance
  create_products_from_reactants_graph(bng_data, reactants_graph, created_products);

#ifdef DEBUG_CPLX_MATCHING
  cout << "Resulting products:\n";
  for (auto& c: created_products) {
    c.dump(false);
    cout << "\n";
  }
#endif

}


bool RxnRule::is_cplx_reactant_on_both_sides_of_rxn(const uint index) const {
  assert(is_finalized());
  for (const CplxIndexPair& cplx_index_pair: cplx_mapping) {
    if (index == cplx_index_pair.reactant_index) {
      return true;
    }
  }
  return false;
}


bool RxnRule::is_cplx_product_on_both_sides_of_rxn(const uint index) const {
  assert(is_finalized());
  for (const CplxIndexPair& cplx_index_pair: cplx_mapping) {
    if (index == cplx_index_pair.product_index) {
      return true;
    }
  }
  return false;
}


bool RxnRule::find_assigned_mol_reactant_for_product(const CplxMolIndex& product_cmi, CplxMolIndex& reactant_cmi) const {
  // this is not a time critical search
  for (const CMIndexPair& cmi_pair: mol_mapping) {
    if (product_cmi == cmi_pair.product_cmi) {
      reactant_cmi = cmi_pair.reactant_cmi;
      return true;
    }
  }
  return false;
}


// Finds a matching already not assigned product,
// the product must be
// NOTE: BNGL2.pl provides more detailed reporting, see tests N220, N230-N232
bool RxnRule::find_most_fitting_unassigned_mol_product(const CplxMolIndex& reactant_cmi, CplxMolIndex& best_product_cmi) const {
  const MolInstance& reactant_mol_inst = get_mol_reactant(reactant_cmi);

  int best_score = -1;
  CplxMolIndex best_cmi;

  for (uint complex_index = 0; complex_index < products.size(); complex_index++) {
    for (uint molecule_index = 0; molecule_index < products[complex_index].mol_instances.size(); molecule_index++) {

      CplxMolIndex product_cmi(complex_index, molecule_index);

      // must have the same molecule type
      const MolInstance& product_mol_inst = get_mol_product(product_cmi);
      if (reactant_mol_inst.mol_type_id != product_mol_inst.mol_type_id) {
        continue;
      }

      // and must not be assigned
      CplxMolIndex found_cmi_ignored;
      bool found = find_assigned_mol_reactant_for_product(product_cmi, found_cmi_ignored);
      if (found) {
        continue;
      }

      // ok, this product was not mapped yet
      int num_explicitly_listed_components_in_reactant = 0;
      int num_same_explicitly_listed_components = 0;
      int num_same_component_states = 0; // when the components are expl. listed

      assert(reactant_mol_inst.component_instances.size() == product_mol_inst.component_instances.size());
      for (uint i = 0; i < reactant_mol_inst.component_instances.size(); i++) {

        const ComponentInstance& reactant_comp_inst = reactant_mol_inst.component_instances[i];
        const ComponentInstance& product_comp_inst = product_mol_inst.component_instances[i];

        // component must be explicitly listed on both sides to be considered
        if (reactant_comp_inst.explicitly_listed_in_pattern) {
          num_explicitly_listed_components_in_reactant++;

          // the same component must be explicitly listed and
          // if state is specified, it must be set on both sides
          if (product_comp_inst.explicitly_listed_in_pattern &&
              reactant_comp_inst.state_is_set() == product_comp_inst.state_is_set())  {

            num_same_explicitly_listed_components++;
            if (reactant_comp_inst.state_id == product_comp_inst.state_id) {
              num_same_component_states++;
            }
          }
        }
        else if (product_comp_inst.explicitly_listed_in_pattern) {
          // component is not listed in reactants, just in products
          num_same_explicitly_listed_components--;
        }
      }

      // all listed components match?
      if (num_explicitly_listed_components_in_reactant == num_same_explicitly_listed_components) {
        if (num_same_component_states > best_score) {
          best_score = num_same_component_states;
          best_cmi = product_cmi;
        }
      }
    }
  }

  if (best_score != -1) {
    best_product_cmi = best_cmi;
    return true;
  }
  else {
    return false;
  }
}


// check if it makes sense to compute molecule_mapping at all
bool RxnRule::has_same_mols_in_reactants_and_products() const {
  map<mol_type_id_t, int> reactant_types, product_types;

  for (const CplxInstance& ci: reactants) {
    for (const MolInstance& mi: ci.mol_instances) {
      if (reactant_types.count(mi.mol_type_id)) {
        reactant_types[mi.mol_type_id]++;
      }
      else {
        reactant_types[mi.mol_type_id] = 1;
      }
    }
  }

  for (const CplxInstance& ci: products) {
    for (const MolInstance& mi: ci.mol_instances) {
      if (product_types.count(mi.mol_type_id)) {
        product_types[mi.mol_type_id]++;
      }
      else {
        product_types[mi.mol_type_id] = 1;
      }
    }
  }

  return reactant_types == product_types;
}


bool RxnRule::find_assigned_cplx_reactant_for_product(const uint product_index, uint& reactant_index) const {
  // this is not a time critical search
  for (const CplxIndexPair& cplx_index_pair: cplx_mapping) {
    if (product_index == cplx_index_pair.product_index) {
      reactant_index = cplx_index_pair.reactant_index;
      return true;
    }
  }
  return false;
}


// a matching reactant and product must be identical
void RxnRule::compute_cplx_reactants_products_mapping() {

  cplx_mapping.clear();

  for (uint ri = 0; ri < reactants.size(); ri++) {
    for (uint pi = 0; pi < products.size(); pi++) {
      uint index_ignored;
      if (find_assigned_cplx_reactant_for_product(pi, index_ignored)) {
        // already used
        continue;
      }

      if (reactants[ri].matches_fully(products[pi], true)) {
        cplx_mapping.push_back(CplxIndexPair(ri, pi));
        // reactant was mapped, continue with the next reactant
        break;
      }
    }
  }
}


bool RxnRule::compute_mol_reactants_products_mapping(MolInstance& not_matching_mol_inst, CplxMolIndex& not_matching_cmi) {
  mol_mapping.clear();

  mol_instances_are_fully_maintained = has_same_mols_in_reactants_and_products();

  for (uint complex_index = 0; complex_index < reactants.size(); complex_index++) {
    for (uint molecule_index = 0; molecule_index < reactants[complex_index].mol_instances.size(); molecule_index++) {

      CplxMolIndex reactant_cmi = CplxMolIndex(complex_index, molecule_index);
      CplxMolIndex product_cmi;
      bool found = find_most_fitting_unassigned_mol_product(reactant_cmi, product_cmi);

      if (found) {
        mol_mapping.push_back(CMIndexPair(reactant_cmi, product_cmi));
      }
      else if (mol_instances_are_fully_maintained) {
        // reporting error only if there should be a full match
        not_matching_mol_inst = get_mol_reactant(reactant_cmi);
        not_matching_cmi = reactant_cmi;

        return false;
      }
    }
  }

  return true;
}


bool RxnRule::compute_reactants_products_mapping() {

  compute_cplx_reactants_products_mapping();

  MolInstance not_matching_mol_inst_ignored;
  CplxMolIndex not_matching_cmi_ignored;
  bool ok = compute_mol_reactants_products_mapping(not_matching_mol_inst_ignored, not_matching_cmi_ignored);
  assert(ok);
  return ok;
}


bool RxnRule::compute_reactants_products_mapping_w_error_output(const BNGData& bng_data, std::ostream& out) {

  compute_cplx_reactants_products_mapping();

  // NOTE: we might need to direct the molecule mapping using cplx mapping,
  // but let's see later

  MolInstance not_matching_mol_inst;
  CplxMolIndex not_matching_cmi;
  bool ok = compute_mol_reactants_products_mapping(not_matching_mol_inst, not_matching_cmi);
  if (!ok) {
    out << "Did not find a matching molecule in products for reactant molecule ";
    out << not_matching_mol_inst.to_str(bng_data, true);
    out << " listed as complex " << not_matching_cmi.cplx_index << " and molecule " << not_matching_cmi.mol_index << ".";
  }
  return ok;
}


void RxnRule::move_products_that_are_also_reactants_to_be_the_first_products() {

  // for each reactant (from the end since we want the products to be ordered in the same way)
  for (int pi = products.size() - 1; pi > 0; pi--) {
    uint ri;
    bool found = find_assigned_cplx_reactant_for_product(pi, ri);

    if (found) {
      // move product to the front
      CplxInstance prod = products[pi];
      products.erase(products.begin() + pi);
      products.insert(products.begin(), prod);

      // update mapping
      compute_reactants_products_mapping();
    }
  }
}


bool RxnRule::species_can_be_reactant(const species_id_t id, const SpeciesContainer& all_species) {

  // check caches first
  if (species_applicable_as_reactants.count(id) != 0) {
    return true;
  }
  if (species_not_applicable_as_reactants.count(id) != 0) {
    return false;
  }

  // need to find out
  const CplxInstance& inst = all_species.get_as_cplx_instance(id);

  // at least one should match
  bool matches = false;
  for (const CplxInstance& reactant: reactants) {
    if (reactant.matches_pattern(inst, true)) {
      matches = true;
      break;
    }
    else {
      matches = false;
    }
  }

  if (matches) {
    species_applicable_as_reactants.insert_unique(id);
  }
  else {
    species_not_applicable_as_reactants.insert_unique(id);
  }

  return matches;
}


bool RxnRule::species_is_both_bimol_reactants(const species_id_t id, const SpeciesContainer& all_species) {

  if (!is_bimol()) {
    return false;
  }

  // check if the species can be a reactant at all
  if (!species_can_be_reactant(id, all_species)) {
    return false;
  }

  // then the reactants must be identical (this can be precomputed)
  bool res = reactants[0].matches_pattern(reactants[1], true);
  assert(res == reactants[1].matches_pattern(reactants[0], true) && "Pattern identity must be bijective");
  return res;
}


bool RxnRule::update_variable_rxn_rate(const float_t current_time, const RxnClass* requester) {
  if (!may_update_rxn_rate()) {
    return false;
  }
  assert(!variable_rates.empty());
  assert(next_variable_rate_index < variable_rates.size());

  if (variable_rates[next_variable_rate_index].time > current_time
      && !cmp_eq(variable_rates[next_variable_rate_index].time, current_time) ) {
    return false;
  }

  // find which time to use - the highest but still smaller than the following one
  size_t current_index = next_variable_rate_index;
  while (current_index < variable_rates.size() &&
          (current_time > variable_rates[current_index + 1].time ||
           cmp_eq(current_time, variable_rates[current_index + 1].time)
        )
  ) {
    current_index++;
  }

  // should we use the last entry?
  if (current_index == variable_rates.size() &&
      next_variable_rate_index == current_index - 1) {
    current_index = next_variable_rate_index;
  }

  // current_time >= time for next change
  rate_constant = variable_rates[current_index].rate_constant;
  next_variable_rate_index = current_index + 1;

  // notify parents that update is needed
  for (RxnClass* user: rxn_classes_where_used) {
    // do not call update on the class that called us
    if (user != requester) {
      user->update_rxn_rates_if_needed(current_time);
    }
  }

  return true;
}


std::string RxnRule::to_str(const BNGData& bng_data) const {
  stringstream ss;
  ss << name << ": ";

  ss << complex_instance_vector_to_str(bng_data, reactants);
  ss << " -> ";
  ss << complex_instance_vector_to_str(bng_data, products);

  ss << " " << rate_constant;

  return ss.str();
}


std::string RxnRule::complex_instance_vector_to_str(const BNGData& bng_data, const CplxInstanceVector& complexes) const {
  stringstream ss;
  for (size_t i = 0; i < complexes.size(); i++) {
    ss << complexes[i].to_str(bng_data, is_surf_rxn());

    if (i != complexes.size() - 1) {
      ss << " + ";
    }
  }
  return ss.str();
}


std::string RxnRule::reactants_to_str(const BNGData& bng_data) const {
  return complex_instance_vector_to_str(bng_data, reactants);
}


std::string RxnRule::products_to_str(const BNGData& bng_data) const {
  return complex_instance_vector_to_str(bng_data, products);
}


void RxnRule::dump_complex_instance_vector(
    const BNGData& bng_data, const CplxInstanceVector& complexes,
    const std::string ind) const {

  for (size_t i = 0; i < complexes.size(); i++) {
    cout << ind << "CplxInstance " << i << ":\n";
    complexes[i].dump(true, ind + "  ");
    cout << "\n";
  }
}


void RxnRule::dump(const BNGData& bng_data, const bool for_diff, const std::string ind) const {
  if (!for_diff) {
    cout << ind << to_str(bng_data);
  }
  else {
    cout << ind << "name: " << name << "\n";
    cout << ind << "id: " << id << "\n";

    cout << ind << "type: ";
    switch (type) {
      case RxnType::Standard:
        cout << "Standard";
        break;
      case RxnType::AbsorbRegionBorder:
        cout << "AbsorbRegionBorder";
        break;
      case RxnType::Reflect:
        cout << "Reflect";
        break;
      case RxnType::Transparent:
        cout << "Transparent";
        break;
      default:
        assert(false);
    }
    cout << "\n";

    cout << ind << "rate_constant: " << rate_constant << "\n";
    cout << ind << "variable_rates.size: " << variable_rates.size() << "\n";

    if (!variable_rates.empty()) {
      for (size_t i = 0; i < variable_rates.size(); i++) {
        cout << ind + "  " << "t: " << variable_rates[i].time << ", r: " << variable_rates[i].rate_constant << "\n";
      }
    }

    cout << ind << "mol_instances_are_fully_maintained: " << mol_instances_are_fully_maintained << "\n";

    cout << ind << "reactants:\n";
    dump_complex_instance_vector(bng_data, reactants, ind);
    cout << ind << "products:\n";
    dump_complex_instance_vector(bng_data, products, ind);
  }
}


} /* namespace BNG */
