/*
 * rule.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#include <iostream>
#include <sstream>
#include <vector>


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
    if (first < b.first) {
      return true;
    }
    else if (first == b.first) {
      return second < b.second;
    }
    else {
      return false;
    }
  }

  vertex_descriptor_t first;
  vertex_descriptor_t second;
};


static void merge_graphs(Graph& srcdst, const Graph& src) {
  boost::copy_graph(src, srcdst);
}


// appends found pathways into the pathways array
void RxnRule::define_rxn_pathways_for_specific_reactants(
    SpeciesContainer& all_species,
    const BNGConfig& bng_config,
    const species_id_t reactant_a_species_id,
    const species_id_t reactant_b_species_id,
    const float_t pb_factor,
    RxnClassPathwayVector& pathways
) const {

  if (is_simple()) {
    std::vector<species_id_t> product_species;
    for (const CplxInstance& product: products) {
      // simple product is deterministic
      // simple species are defined mcell3 mode but in BNG mode they
      // may not have been created (they are based on molecule types)
      species_id_t species_id = all_species.find_full_match(product);
      if (species_id == SPECIES_ID_INVALID) {
        species_id = all_species.add(Species(product, *bng_data, bng_config));
        assert(species_id != SPECIES_ID_INVALID);
      }
      product_species.push_back(species_id);
    }

    float_t prob;
    if (cmp_eq(get_rate_constant(), FLT_GIGANTIC)) {
      // special surface reactions are not scaled
      prob = get_rate_constant();
    }
    else {
      prob = get_rate_constant() * pb_factor;
    }
    pathways.push_back(RxnClassPathway(id, prob, product_species));
  }
  else {
    // TODO LATER: we also need to maintain the IDs of the elementary molecules
    // but such a thing is not supported at all yet
    //
    // when the rule modifies an elementary molecule, but there is multiple
    // elementary molecules that match the reaction pattern such as in:
    // complex A(a,b~X).A(a,b~Y) and rule A(a) -> A(a!1).B(b!1).
    // there the rule can be applied to one of the distinct elementary molecules
    //
    // also having identical components inside of an elementary molecule
    // may lead to nondeterminism:
    // complex A(a,a~X) and rule A(a) -> A(a~Y),
    // there the rule can be applied to one of the distinct components

    // for now, we support reaction with a single outcome
    // and since we will depend on the input (once we will be maintaining IDs
    // of the elementary molecules), we do this computation on the fly
    vector<const CplxInstance*> reactants;
    // downcast, CplxInstance is sufficient for product computation
    reactants.push_back(dynamic_cast<const CplxInstance*>(&all_species.get(reactant_a_species_id)));
    if (is_bimol()) {
      reactants.push_back(dynamic_cast<const CplxInstance*>(&all_species.get(reactant_b_species_id)));
    }

    // we might get multiple matches on reactant(s), the numbers
    // of matches multiplied give us the total number of variants
    // a single random number then is used to choose a single variant
    vector<vector<CplxInstance>> product_sets;
    create_products_for_complex_rxn(
        reactants,
        product_sets
    );

    for (const auto& product_cplxs: product_sets) {
      // define the products as species
      std::vector<species_id_t> product_species;
      for (const CplxInstance& product: product_cplxs) {
        species_id_t species_id = all_species.find_or_add(Species(product, *bng_data, bng_config));
        assert(species_id != SPECIES_ID_INVALID);
        product_species.push_back(species_id);
      }

      assert(pb_factor != 0);
      assert(!cmp_eq(get_rate_constant(), FLT_GIGANTIC));

      // the probability is divided by the number of mapping of pattern onto pattern
      // because so many more product sets we will get
      float_t prob = get_rate_constant() * pb_factor / num_patterns_onto_patterns_mapping;

      pathways.push_back(RxnClassPathway(id, prob, product_species));
    }
  }
}

void RxnRule::create_patterns_graph() {
  patterns_graph.clear();
  // create graph for reactant patterns
  patterns_graph = reactants[0].get_graph();
  if (reactants.size() == 2) {
    merge_graphs(patterns_graph, reactants[1].get_graph());
  }
}


void RxnRule::create_products_graph() {
  // mark each of the original product complexes with the index of the product
  // we have our own copies and it is easier to mark the source graphs before merging
  for (size_t i = 0; i < products.size(); i++) {
    Graph& graph = products[i].get_graph();
    VertexNameMap index = boost::get(boost::vertex_name, graph);
    typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> it;
    for (it = boost::vertices(graph); it.first != it.second; ++it.first) {
      vertex_descriptor_t desc = *it.first;
      Node& node = index[desc];
      node.product_index = i;
    }
  }

  // create graph from products
  if (products.size() > 0) {
    products_graph = products[0].get_graph();
  }
  for (size_t i = 1; i < products.size(); i++) {
    merge_graphs(products_graph, products[i].get_graph());
  }
}


void RxnRule::compute_num_patterns_onto_patterns_mapping() {
  assert(boost::num_vertices(patterns_graph) != 0 &&
      "Patterns graph must have been initialized");

  if (is_unimol() && products.size() == 1) {
    // single reactant and single product ->
    // we are changing states
    VertexMappingVector pat_to_pat_mappings;
    get_subgraph_isomorphism_mappings(
        patterns_graph,
        patterns_graph,
        false,
        pat_to_pat_mappings
    );
    num_patterns_onto_patterns_mapping = pat_to_pat_mappings.size();
  }
  else {
    // multiple reactants or multiple products ->
    // we are creating bonds

    // prepare graph from reactants, only they will have all
    // components present
    CplxInstanceVector full_reactants;
    for (const CplxInstance& reac: reactants) {

      // make a copy
      CplxInstance full_reac = reac;

      // add missing components
      for (MolInstance& mi: full_reac.mol_instances) {
        mi.insert_missing_components_as_any_state_pattern(*bng_data);
      }
      // reinitialize
      full_reac.finalize();
      full_reactants.push_back(full_reac);
    }

    Graph full_patterns_graph;
    // create graph for reactant patterns
    full_patterns_graph = full_reactants[0].get_graph();
    if (full_reactants.size() == 2) {
      merge_graphs(full_patterns_graph, full_reactants[1].get_graph());
    }

    VertexMappingVector pat_to_pat_mappings;
    get_subgraph_isomorphism_mappings(
        patterns_graph,
        full_patterns_graph,
        false,
        pat_to_pat_mappings
    );

    num_patterns_onto_patterns_mapping = pat_to_pat_mappings.size();
  }
}


void RxnRule::finalize() {
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

  create_patterns_graph();
  create_products_graph();
  compute_num_patterns_onto_patterns_mapping();
  compute_reactants_products_mapping();

  // for MCell3 compatibility, updates mapping when needed
  move_products_that_are_also_reactants_to_be_the_first_products();

  set_finalized();
}


struct MolCompInfo {
  MolCompInfo(const vertex_descriptor_t desc_, const MolInstance* mi_)
    : desc(desc_), already_matched(false), mi(mi_), ci(nullptr) {
  }
  MolCompInfo(const vertex_descriptor_t desc_, const ComponentInstance* ci_)
    : desc(desc_), already_matched(false), mi(nullptr), ci(ci_) {
  }

  const MolInstance* get_mi() const {
    assert(mi != nullptr);
    return mi;
  }

  const ComponentInstance* get_ci() const {
    assert(ci != nullptr);
    return ci;
  }

  vertex_descriptor_t desc;

  vector<int> matching_score;
  bool already_matched;
private:
  const MolInstance* mi;
  const ComponentInstance* ci;
};


static void get_all_mol_instances_from_graph(
    Graph& graph,
    vector<MolCompInfo>& res
) {
  VertexNameMap index = boost::get(boost::vertex_name, graph);
  typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> it;
  for (it = boost::vertices(graph); it.first != it.second; ++it.first) {
    vertex_descriptor_t desc = *it.first;
    const Node& mi_node = index[desc];
    if (mi_node.is_mol) {
      res.push_back(MolCompInfo(desc, mi_node.mol));
    }
  }
}


static void get_all_component_instances_of_mol_from_graph(
    Graph& graph,
    const vertex_descriptor_t mol_desc, // specifies molecule whose components we are collecting
    vector<MolCompInfo>& res
) {
  VertexNameMap index = boost::get(boost::vertex_name, graph);

  boost::graph_traits<Graph>::out_edge_iterator pat_ei, pat_edge_end;
  for (boost::tie(pat_ei, pat_edge_end) = boost::out_edges(mol_desc, graph); pat_ei != pat_edge_end; ++pat_ei) {
    vertex_descriptor_t connected_node_desc = boost::target(*pat_ei, graph);
    const Node& mi_node = index[connected_node_desc];
    assert(!mi_node.is_mol);
    res.push_back(MolCompInfo(connected_node_desc, mi_node.component));
  }
}


static int get_component_instance_matching_score(
    const ComponentInstance& ci1,
    const ComponentInstance& ci2
) {
  if (ci1.component_type_id != ci2.component_type_id) {
    return 0; // not a match
  }

  int res = 1;
  if (ci1.state_id == ci2.state_id) {
    res++;
  }

  // NOTE: bond matching needs to be improved
  if (ci1.bond_value == ci2.bond_value) {
    res++;
  }

  return res;
}


static int get_mol_instance_matching_score(const MolInstance& pat, const MolInstance& prod) {
  if (pat.mol_type_id != prod.mol_type_id) {
    return -1;
  }
  int res = 0;

  // add a lot of points if we have the same number of components
  if (pat.component_instances.size() == prod.component_instances.size()) {
    res += pat.component_instances.size();
  }

  // and add points for each matching component, they are ordered according to
  // the MolType
  // TODO: probably use the same approach for matching components that is done later
  for (size_t i = 0; i < pat.component_instances.size(); i++) {
    const ComponentInstance& pat_compi = pat.component_instances[i];

    if (i < prod.component_instances.size()) {
      res += get_component_instance_matching_score(pat_compi, prod.component_instances[i]);
    }
  }

  return res;
}


// returns true if the next mapping was found and indices were set
static bool get_best_not_matched_mapping(
    const vector<MolCompInfo>& patterns,
    const vector<MolCompInfo>& products,
    size_t& highest_pat_index,
    size_t& highest_prod_index
) {
  int best_score = -1;

  for (size_t pat_i = 0; pat_i < patterns.size(); pat_i++) {
    if (!patterns[pat_i].already_matched) {
      for (size_t prod_i = 0; prod_i < patterns[pat_i].matching_score.size(); prod_i++) {
        if (!products[prod_i].already_matched && patterns[pat_i].matching_score[prod_i] > best_score) {
          best_score = patterns[pat_i].matching_score[prod_i];
          highest_pat_index = pat_i;
          highest_prod_index = prod_i;
        }
      }
    }
  }

  return best_score != -1;
}


static void select_best_mapping(
    vector<MolCompInfo>& patterns,
    vector<MolCompInfo>& products,
    VertexMapping& prod_reac_mapping
) {
  // sort the score and select the best 'global' match
  // find the globally highest score
  // TODO: check why we are not using the loop induction var.
  for (size_t i = 0; i < patterns.size(); i++) {
    size_t highest_pat_index;
    size_t highest_prod_index;
    bool found = get_best_not_matched_mapping(
        patterns, products, highest_pat_index, highest_prod_index);
    if (!found) {
      break;
    }

    // define this mapping, it is in the product->pattern direction
    vertex_descriptor_t prod_desc = products[highest_prod_index].desc;
    assert(prod_reac_mapping.count(prod_desc) == 0);
    prod_reac_mapping[prod_desc] = patterns[highest_pat_index].desc;

    // and mark that we used these molecule instances
    patterns[highest_pat_index].already_matched = true;
    products[highest_prod_index].already_matched = true;
  }
}


static void find_best_product_to_pattern_mapping(
    Graph& products_graph,
    Graph& patterns_graph,
    VertexMapping& prod_reac_mapping
) {

  prod_reac_mapping.clear();

  // prepare arrays of patterns and products
  vector<MolCompInfo> pattern_mols;
  get_all_mol_instances_from_graph(patterns_graph, pattern_mols);
  vector<MolCompInfo> product_mols;
  get_all_mol_instances_from_graph(products_graph, product_mols);


  // compute matching score for each pair of patterns and products
  for (MolCompInfo& pat: pattern_mols) {
    pat.matching_score.resize(product_mols.size());
    for (size_t i = 0; i < product_mols.size(); i++) {
      MolCompInfo& prod = product_mols[i];
      // TODO: deal also with bonds, when a molecule instance is connected to
      // certain molecules and the product as well, this should give it higher score
      pat.matching_score[i] = get_mol_instance_matching_score(*pat.get_mi(), *prod.get_mi());
    }
  }

  // get best mapping for molecule instances
  select_best_mapping(
      pattern_mols,
      product_mols,
      prod_reac_mapping
  );

  // ok, we matched molecules, now we need to match components as well
  // we are adding to the prod_reac_mapping, so we will make a copy to iterate over it
  // NOTE: in other places we use just order, we must use a similar approach
  // probably just define some sorting...
  VertexMapping mol_prod_reac_mapping = prod_reac_mapping;
  for (auto pair_it: mol_prod_reac_mapping) {
    vector<MolCompInfo> pattern_comps;
    get_all_component_instances_of_mol_from_graph(patterns_graph, pair_it.second, pattern_comps);
    vector<MolCompInfo> product_comps;
    get_all_component_instances_of_mol_from_graph(products_graph, pair_it.first, product_comps);

    // compute matching score for each pair of patterns and products
    for (MolCompInfo& pat: pattern_comps) {
      pat.matching_score.resize(product_comps.size());
      for (size_t i = 0; i < product_comps.size(); i++) {
        MolCompInfo& prod = product_comps[i];
        pat.matching_score[i] = get_component_instance_matching_score(*pat.get_ci(), *prod.get_ci());
      }
    }

    // get best mapping for component instances
    select_best_mapping(
        pattern_comps,
        product_comps,
        prod_reac_mapping
    );
  }
}


// vertex pointed to by reac_desc
// must be a component and must have a single bond
// to another component and also a bond to its molecule instance
// returns -1 when the component is not connected to another component,
// asserts in debug mode when must_exist is set to true and target was not found
//
// the release_assert calls are here because boost's library produced
// weird errors with release build, they are cheap anyway so they are kept there
const vertex_descriptor_t TARGET_NOT_FOUND = -1;
static vertex_descriptor_t get_bond_target(
    Graph& graph,
    const vertex_descriptor_t desc,
    const bool must_exist = true
) {
  #ifdef DEBUG_CPLX_RXNS
    cout << "get_bond_target:\n";
    dump_graph(graph);
    cout << "desc:" << (int)desc << "\n";
  #endif

  bool comp_found = false;
  bool mol_found = false;

  vertex_descriptor_t res = TARGET_NOT_FOUND;
  VertexNameMap index = boost::get(boost::vertex_name, graph);

  boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
  for (boost::tie(ei,edge_end) = boost::out_edges(desc, graph); ei != edge_end; ++ei) {
    Graph::edge_descriptor e_desc = *ei;
    Graph::vertex_descriptor n_desc = boost::target(e_desc, graph);
    const Node& n = index[n_desc];

    #ifdef DEBUG_CPLX_RXNS
      cout << "  checking " << n << "\n";
    #endif

    if (n.is_mol) {
      release_assert(!mol_found);
      mol_found = true; // just for debug
    }
    else {
      release_assert(!comp_found);
      comp_found = true;
      res = n_desc;
    }
  }

  release_assert(mol_found && "Component must be connected to its molecule");
  if (must_exist) {
    release_assert(comp_found);
  }

  return res;
}

// returns TARGET_NOT_FOUND when there is no mapping from products to pattern
vertex_descriptor_t get_new_bond_target(
    const Graph& reactants_graph,
    const VertexMapping& pattern_reactant_mapping,
    const VertexMapping& prod_pattern_mapping,
    Graph& products_graph,
    const vertex_descriptor_t prod_desc
) {

  // to which node we should connect?
  vertex_descriptor_t prog_target_desc = get_bond_target(products_graph, prod_desc);
  #ifdef DEBUG_CPLX_RXNS
    cout << "For " << (int)prod_desc << " found " << (int)prog_target_desc << "\n";
  #endif
  release_assert(prog_target_desc != TARGET_NOT_FOUND);

  // to which reactant we will point?
  auto target_prod_pat_it = prod_pattern_mapping.find(prog_target_desc);
  if (target_prod_pat_it == prod_pattern_mapping.end()) {
    return TARGET_NOT_FOUND;
  }
  auto target_pat_reac_it = pattern_reactant_mapping.find(target_prod_pat_it->second);
  assert(target_pat_reac_it != pattern_reactant_mapping.end() && "Mapping must exist");

  return target_pat_reac_it->second;
}


static void mark_consumed_reactants(
    Graph& reactants_graph,
    const VertexMapping& pattern_reactant_mapping,
    const VertexMapping& product_pattern_mapping
) {
  // we need to remove anything that forms a disconnected graph and is
  // not matched by the pattern

  size_t num_vertices = boost::num_vertices(reactants_graph);
  // compute connected components
  vector<int> component_per_vertex(num_vertices);
  int num_components = boost::connected_components(
      reactants_graph,
      boost::make_iterator_property_map(
          component_per_vertex.begin(), boost::get(boost::vertex_index, reactants_graph), component_per_vertex[0]
      )
  );

  // which of the components are unused
  vector<bool> used_components(num_components);

  // we start from products, and go through all the mappings from the product
  // to the pattern
  for (auto prod_pat_it: product_pattern_mapping) {
    vertex_descriptor_t pat_desc = prod_pat_it.second;
    // ok, we know that a product maps to this pattern vertex,
    // there must a be mapping onto the reactant otherwise this rxn cound
    // not have been triggered
    auto pat_reac_it = pattern_reactant_mapping.find(pat_desc);
    assert(pat_reac_it != pattern_reactant_mapping.end());
    vertex_descriptor_t reac_desc = pat_reac_it->second;

    // and mark that we must keep this graph component, i.e. at least
    // one of the graph components is kept
    used_components[component_per_vertex[reac_desc]] = true;
  }

  // we must not remove any vertices here because it would break the pattern_reactant_mapping,
  // simply mark the nodes that we do not want to include
  VertexNameMap index = boost::get(boost::vertex_name, reactants_graph);

  // and now mark each vertex based on whether it belongs to a component that we will keep
  // this information is used later in convert_graph_component_to_product_cplx_inst
  typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> mol_it;
  for (mol_it = boost::vertices(reactants_graph); mol_it.first != mol_it.second; ++mol_it.first) {
    Graph::vertex_descriptor reac_desc = *mol_it.first;

    Node& mol = index[reac_desc];
    mol.used_in_rxn_product = used_components[component_per_vertex[reac_desc]];
  }
}


static void add_new_products(
    Graph& reactants_graph,
    const VertexMapping& pattern_reactant_mapping,
    const VertexMapping& product_pattern_mapping,
    Graph& products_graph,
    const VertexNameMap& products_index
) {
  // this set will contain a pair of product graph_descriptors for bonds between components
  set<UnorderedPair> product_bonds_to_add_to_reactants;

  // we also create a new mapping directly from all products onto reactants (bypassing patterns)
  VertexMapping product_reactant_mapping;

  // for each molecule in product graph:
  std::pair<vertex_iter_t, vertex_iter_t> mol_it;
  for (mol_it = boost::vertices(products_graph); mol_it.first != mol_it.second; ++mol_it.first) {
    vertex_descriptor_t prod_desc = *mol_it.first;
    auto prod_pat_it = product_pattern_mapping.find(prod_desc);
    const Node& prod_mol_node = products_index[prod_desc];

    // if we are dealing with a molecule and there is no mapping to reactant pattern graph
    if (prod_mol_node.is_mol && prod_pat_it == product_pattern_mapping.end()) {
      // create a new molecule instance in the reactants_graph
      // (the node points to a MolInst owned by products of this rxn)
      vertex_descriptor_t new_reac_desc = boost::add_vertex(prod_mol_node, reactants_graph);
      product_reactant_mapping[prod_desc] = new_reac_desc;

      // for each component of the molecule in product graph
      boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
      for (boost::tie(ei,edge_end) = boost::out_edges(prod_desc, products_graph); ei != edge_end; ++ei) {
        Graph::edge_descriptor e_desc = *ei;
        Graph::vertex_descriptor prod_comp_desc = boost::target(e_desc, products_graph);
        const Node& comp_node = products_index[prod_comp_desc];
        assert(!comp_node.is_mol);

        // add this component and connect it to the molecule
        vertex_descriptor_t new_comp_desc = boost::add_vertex(comp_node, reactants_graph);
        boost::add_edge(new_reac_desc, new_comp_desc, reactants_graph);
        product_reactant_mapping[prod_comp_desc] = new_comp_desc;

        vertex_descriptor_t prod_bond_target = get_bond_target(products_graph, prod_comp_desc, false);
        if (prod_bond_target != TARGET_NOT_FOUND) {
          // the target component might have not been created, so we must remember to make the bond later
          product_bonds_to_add_to_reactants.insert(UnorderedPair(prod_comp_desc, prod_bond_target));
        }
      }
    }
    else {
      // mapped molecule or mapped component - need to define direct mapping as well
      if (prod_pat_it != product_pattern_mapping.end()) {
        auto pat_reac_it = pattern_reactant_mapping.find(prod_pat_it->second);
        assert(pat_reac_it != pattern_reactant_mapping.end());
        product_reactant_mapping[prod_pat_it->first] = pat_reac_it->second;
      }
    }
  }

  // create bonds between new components
  for (auto p: product_bonds_to_add_to_reactants) {
    assert(product_reactant_mapping.count(p.first) != 0);
    assert(product_reactant_mapping.count(p.second) != 0);
    vertex_descriptor_t reac_comp_desc1 = product_reactant_mapping[p.first];
    vertex_descriptor_t reac_comp_desc2 = product_reactant_mapping[p.second];
    boost::add_edge(reac_comp_desc1, reac_comp_desc2, reactants_graph);
  }
}

// manipulate the reactants_graph according to how
// products should look like
static void apply_rxn_on_reactants_graph(
    Graph& reactants_graph,
    const VertexMapping& pattern_reactant_mapping,
    const Graph& pattern_graph,
    const VertexMapping& product_pattern_mapping,
    Graph& products_graph
) {
  set<vertex_descriptor_t> mol_instances_to_keep;

  VertexNameMap reactants_index = boost::get(boost::vertex_name, reactants_graph);
  VertexNameMap products_index = boost::get(boost::vertex_name, products_graph);

  // mark all disconnected graphs whose molecule instance is not mapped to
  // by products
  mark_consumed_reactants(reactants_graph, pattern_reactant_mapping, product_pattern_mapping);

  set<UnorderedPair> bonds_to_remove;
  set<UnorderedPair> bonds_to_add;

  //
  // for each component in product graph:
  //   find corresponding component in reactant pattern
  //     if there is no such mapping, ignore it because it might be a component of a
  //     new molecule instance
  for (auto prod_pat_it: product_pattern_mapping) {
    vertex_descriptor_t prod_desc = prod_pat_it.first;
    const Node& prod_node = products_index[prod_pat_it.first];
    // product->pattern
    vertex_descriptor_t pat_desc = prod_pat_it.second;

    // find component corresponding to reactant pattern in reactants_graph
    //   this mapping must exist because we matched the pattern graph to reactants graph
    auto pat_reac_it = pattern_reactant_mapping.find(pat_desc);
    assert(pat_reac_it != pattern_reactant_mapping.end() && "Mapping should exist?");

    vertex_descriptor_t reac_desc = pat_reac_it->second;

    // manipulate state and or bond
    Node& reac_node = reactants_index[reac_desc];

    // set product index marker
    assert(prod_node.product_index != INDEX_INVALID);
    reac_node.product_index = prod_node.product_index;

    if (!prod_node.is_mol) {
      assert(!reac_node.is_mol);
      const ComponentInstance& prod_ci = *prod_node.component;
      ComponentInstance& reac_ci = *reac_node.component;

      // update state
      if (prod_ci.state_is_set() && prod_ci.state_id != reac_ci.state_id) {
        reac_ci.state_id = prod_ci.state_id;
      }

      // and bond,
      // assuming that there will be no change for !?
      if (prod_ci.bond_value != reac_ci.bond_value &&
          prod_ci.bond_value != BOND_VALUE_ANY &&
          reac_ci.bond_value != BOND_VALUE_ANY
      ) {
        // orig: !+
        if (reac_ci.bond_value == BOND_VALUE_BOUND) {
          // new: (no bond)
          if (prod_ci.bond_value == BOND_VALUE_UNBOUND) {
            bonds_to_remove.insert(UnorderedPair(
                reac_desc,
                get_bond_target(reactants_graph, reac_desc)
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
        else if (reac_ci.bond_value == BOND_VALUE_UNBOUND) {
          // new: !1
          if (prod_ci.bond_has_numeric_value()) {

            vertex_descriptor_t target_reac_desc = get_new_bond_target(
                reactants_graph,
                pattern_reactant_mapping,
                product_pattern_mapping,
                products_graph,
                prod_desc
            );

            // target does not have to exist when this is a new product,
            // it will be added to the reactants graph later
            if (target_reac_desc != TARGET_NOT_FOUND) {
              bonds_to_add.insert(UnorderedPair(reac_desc, target_reac_desc));
            }
          }
          // new: !+
          else if (prod_ci.bond_value == BOND_VALUE_BOUND){
            assert(false && "Cannot change bond from to !+");
          }
          else {
            assert(false);
          }
        }
        // orig: !1
        else if (reac_ci.bond_has_numeric_value()) {
          // new: (no bond)
          if (prod_ci.bond_value == BOND_VALUE_UNBOUND) {
            bonds_to_remove.insert(UnorderedPair(
                reac_desc,
                get_bond_target(reactants_graph, reac_desc)
            ));
          }
          // new: !2
          else if (prod_ci.bond_has_numeric_value()) {
            assert(prod_ci.bond_value != reac_ci.bond_value);
            // remove original one
            bonds_to_remove.insert(UnorderedPair(
                reac_desc,
                get_bond_target(reactants_graph, reac_desc)
            ));

            // and create a new one
            vertex_descriptor_t target_reac_desc = get_new_bond_target(
                reactants_graph,
                pattern_reactant_mapping,
                product_pattern_mapping,
                products_graph,
                prod_desc
            );
            // target does not have to exist when this is a new product,
            // it will be added to the reactants graph later
            if (target_reac_desc != TARGET_NOT_FOUND) {
              bonds_to_add.insert(UnorderedPair(reac_desc, target_reac_desc));
            }
          }
          // new: !+
          else if (prod_ci.bond_value == BOND_VALUE_BOUND){
            // ignored - product has '+' and for the reactant, there is already some bond set
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

  // add new products for molecules and components present in the
  // products graph but missing in the pattern graph
  add_new_products(
      reactants_graph, pattern_reactant_mapping,
      product_pattern_mapping, products_graph, products_index
  );
}


// returns false if any of the encountered nodes has its
// used_in_rxn_product set to false which means that this is
// a reactant that was consumed by a reaction
static bool convert_graph_component_to_product_cplx_inst(
    Graph& graph,
    const vector<int>& graph_components,
    const int graph_component_index,
    CplxInstance& cplx,
    uint& product_index
) {
  product_index = INDEX_INVALID;

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

    // check molecule whether it has a marker that says which product it was
    // we need to maintain the ordering of products
    if (mol.product_index != INDEX_INVALID) {
      assert(product_index == INDEX_INVALID || product_index == mol.product_index);
      product_index = mol.product_index;
    }

    if (!mol.used_in_rxn_product) {
      // this is a consumed reactant
      return false;
    }

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
      Graph::vertex_descriptor bound_comp_desc = get_bond_target(graph, comp_desc, false);
      if (bound_comp_desc == TARGET_NOT_FOUND) {
        compi.bond_value = BOND_VALUE_UNBOUND;
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

  return true;
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
  vector<int> graph_components(boost::num_vertices(reactants_graph));
  int num_graph_components = boost::connected_components(
      reactants_graph,
      boost::make_iterator_property_map(
          graph_components.begin(), boost::get(boost::vertex_index, reactants_graph), graph_components[0]
      )
  );

  for (int i = 0; i < num_graph_components; i++) {
    CplxInstance product_cplx(bng_data);

    // some reactants may have been removed
    uint product_index;
    bool is_rxn_product =
        convert_graph_component_to_product_cplx_inst(reactants_graph, graph_components, i, product_cplx, product_index);

    if (is_rxn_product) {
      assert(product_index != INDEX_INVALID);

      // products must be ordered by their position in rxn products
      if (product_index >= created_products.size()) {
        created_products.resize(product_index + 1, CplxInstance(bng_data));
      }
      created_products[product_index] = product_cplx;
    }
  }
}


void RxnRule::create_products_for_complex_rxn(
    const vector<const CplxInstance*>& input_reactants,
    vector<vector<CplxInstance>>& created_product_sets
) const {
  // the result of this function is cached in rxnclass
  assert(input_reactants.size() == reactants.size());
  assert(input_reactants.size() == 1 || input_reactants.size() == 2);

  // merge input reactant graphs
  Graph reactants_graph = input_reactants[0]->get_graph();
  if (input_reactants.size() == 2) {
    merge_graphs(reactants_graph, input_reactants[1]->get_graph());
  }

  // compute mapping reactant pattern -> reactant
  VertexMappingVector pattern_reactant_mappings;
  get_subgraph_isomorphism_mappings(
      patterns_graph, // pattern
      reactants_graph, // actual reactant
      false, // do not stop with first match
      pattern_reactant_mappings
  );
  assert(pattern_reactant_mappings.size() != 0 &&
      "Did not find a match of patterns onto reaction.");

  release_assert(pattern_reactant_mappings.size() < MAX_PRODUCT_SETS_PER_RXN
      && "Encountered a huge number of potential product sets for a single reaction");

  // now, for each of the mappings, compute what different products we might get
  for (const VertexMapping& mapping: pattern_reactant_mappings) {

  #ifdef DEBUG_CPLX_MATCHING
    cout << "Products:\n";
    dump_graph(products_graph);
  #endif

    // we need to make a copy of the reactants because we will be modifying them
    vector<CplxInstance> input_reactants_copy;
    for (const CplxInstance* ci: input_reactants) {
      input_reactants_copy.push_back(*ci);
    }
    Graph reactants_graph_copy = input_reactants_copy[0].get_graph();
    if (input_reactants_copy.size() == 2) {
      merge_graphs(reactants_graph_copy, input_reactants_copy[1].get_graph());
    }

    // manipulate nodes using information about products
    apply_rxn_on_reactants_graph(
        reactants_graph_copy,
        mapping,
        patterns_graph,
        products_to_patterns_mapping,
        products_graph
    );

  #ifdef DEBUG_CPLX_MATCHING
    cout << "\nReactants after applying rxn:\n";
    dump_graph(reactants_graph_copy);
  #endif

    // and finally create products, each disconnected graph in the result is a
    // separate complex instance
    created_product_sets.push_back(vector<CplxInstance>());
    create_products_from_reactants_graph(bng_data, reactants_graph_copy, created_product_sets.back());

  #ifdef DEBUG_CPLX_MATCHING
    cout << "Resulting products:\n";
    for (auto& c: created_product_sets.back()) {
      c.dump(false);
      cout << "\n";
    }
  #endif
  }
}


bool RxnRule::is_cplx_reactant_on_both_sides_of_rxn(const uint index) const {
#ifdef MCELL4_DO_NOT_REUSE_REACTANT
  return false;
#endif
  assert(is_finalized());
  for (const CplxIndexPair& cplx_index_pair: simple_cplx_mapping) {
    if (index == cplx_index_pair.reactant_index) {
      return true;
    }
  }
  return false;
}


bool RxnRule::is_cplx_product_on_both_sides_of_rxn(const uint index) const {
#ifdef MCELL4_DO_NOT_REUSE_REACTANT
  return false;
#endif
  assert(is_finalized());
  for (const CplxIndexPair& cplx_index_pair: simple_cplx_mapping) {
    if (index == cplx_index_pair.product_index) {
      return true;
    }
  }
  return false;
}


bool RxnRule::get_assigned_simple_cplx_reactant_for_product(const uint product_index, uint& reactant_index) const {
  // this is not a time critical search
  for (const CplxIndexPair& cplx_index_pair: simple_cplx_mapping) {
    if (product_index == cplx_index_pair.product_index) {
      reactant_index = cplx_index_pair.reactant_index;
      return true;
    }
  }
  return false;
}


static size_t find_mol_instance_with_address(const CplxInstanceVector& cplx_vex, const MolInstance* mi_addr) {
  for (size_t i = 0; i < cplx_vex.size(); i++) {
    const CplxInstance& ci = cplx_vex[i];
    for (size_t k = 0; k < ci.mol_instances.size(); k++) {
      if (mi_addr == &ci.mol_instances[k]) {
        return i;
      }
    }
  }
  release_assert(false && "Molecule instance must be found.");
  return 0;
}


void RxnRule::compute_reactants_products_mapping() {

  // compute mapping reactant patterns -> product patterns
  // nodes in the reaction products graph have an extra option that they can ignore the state
  // find all mappings and compute score for each of those mappings so that we get the best match

  // boost subgraph won't find a mapping of one of the reactants is not present in products,
  // so we need to do this manually, we always stop branching at components of a single molecule
  // so the complexity can be handled in reasonable time because there won't be usually may identical
  // molecules
  find_best_product_to_pattern_mapping(
      products_graph,
      patterns_graph,
      products_to_patterns_mapping
  );

#ifdef DEBUG_CPLX_MATCHING
  dump_graph_mapping(products_to_patterns_mapping);
#endif

  VertexNameMap patterns_index = boost::get(boost::vertex_name, patterns_graph);
  VertexNameMap products_index = boost::get(boost::vertex_name, products_graph);

  // also compute simple_cplx_mapping
  simple_cplx_mapping.clear();
  // and set mol_instances_are_fully_maintained
  uint num_prod_molecule_instances = 0;
  uint num_mapped_molecule_instances = 0;
  // for each molecule instance
  typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> prod_mol_it;
  for (prod_mol_it = boost::vertices(products_graph); prod_mol_it.first != prod_mol_it.second; ++prod_mol_it.first) {
    Graph::vertex_descriptor prod_mol_desc = *prod_mol_it.first;

    const Node& prod_mol = products_index[prod_mol_desc];
    if (!prod_mol.is_mol) {
      continue;
    }

    num_prod_molecule_instances++;

    // is it mapped?
    auto map_it = products_to_patterns_mapping.find(prod_mol_desc);
    if (map_it == products_to_patterns_mapping.end()) {
      continue;
    }

    // count that we have a mapping for this molecule instance
    num_mapped_molecule_instances++;

    // the complex must be simple for the simple_cplx_mapping, i.e. have no edge
    // note that we could compute the full mapping but we need this just for mcell3
    // and code for complexes that are not simple would be longer
    boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
    boost::tie(ei,edge_end) = boost::out_edges(prod_mol_desc, products_graph);
    if (ei != edge_end) {
      continue;
    }

    // ok, we can finally define our mapping
    const Node& pat_mol = patterns_index[map_it->second];
    assert(pat_mol.is_mol);


    // the graph has just a pointer to the molecule instance and since the graph
    // should be independent on the reactions, we are not storing any indices there
    uint reac_cplx_index = find_mol_instance_with_address(reactants, pat_mol.mol);
    uint prod_cplx_index = find_mol_instance_with_address(products, prod_mol.mol);


    simple_cplx_mapping.push_back(CplxIndexPair(reac_cplx_index, prod_cplx_index));
  }

  // count number of molecules in pattern
  uint num_pat_molecule_instances = 0;
  std::pair<vertex_iter, vertex_iter> pat_mol_it;
  for (pat_mol_it = boost::vertices(patterns_graph); pat_mol_it.first != pat_mol_it.second; ++pat_mol_it.first) {
    Graph::vertex_descriptor pat_mol_desc = *pat_mol_it.first;
    const Node& pat_mol = patterns_index[pat_mol_desc];
    if (!pat_mol.is_mol) {
      continue;
    }

    num_pat_molecule_instances++;
  }

  mol_instances_are_fully_maintained =
      num_pat_molecule_instances == num_prod_molecule_instances &&
      num_prod_molecule_instances == num_mapped_molecule_instances;
}


bool RxnRule::check_components_mapping(
    const MolInstance& first_mi,
    const MolInstance& second_mi,
    const char* msg,
    std::ostream& out
) {
  // TODO: use computed component ordering from products_to_patterns_mapping
  // maybe store vertex descriptor along with molecule and component instances
  bool ok = true;
  for (size_t i = 0; i < first_mi.component_instances.size(); i++) {
    const ComponentInstance& first_compi = first_mi.component_instances[i];

    if (i >= second_mi.component_instances.size() ||
        second_mi.component_instances[i].component_type_id != first_compi.component_type_id
    ) {
      out <<
          "Molecule " << msg << ": Component(s) " <<
          bng_data->get_component_type(first_mi.component_instances[i].component_type_id).name <<
          " missing from molecule " <<
          second_mi.to_str(*bng_data) << ".";
      ok = false;
    }
  }

  return ok;
}


bool RxnRule::check_components_states(
    const MolInstance& prod_mi,
    const MolInstance& pat_mi,
    std::ostream& out
) {
  // TODO: use computed component ordering from products_to_patterns_mapping
  bool ok = true;
  for (size_t i = 0; i < prod_mi.component_instances.size(); i++) {
    const ComponentInstance& prod_compi = prod_mi.component_instances[i];

    if (i >= pat_mi.component_instances.size() ||
        pat_mi.component_instances[i].component_type_id == prod_compi.component_type_id
    ) {
      const ComponentInstance& pat_compi = pat_mi.component_instances[i];

      // ok, components have the same type, we need to check state
      if (pat_compi.state_is_set() && !prod_compi.state_is_set()) {
        out <<
            "Component " << bng_data->get_component_type(pat_compi.component_type_id).name <<
            " with state attribute defined in reactant pattern cannot map to component with undefined state attribute in product pattern," <<
            " error for pattern molecule " <<  pat_mi.to_str(*bng_data) << ".";
        ok = false;
      }
    }
  }

  return ok;
}

// called from semantic analyzer
bool RxnRule::check_reactants_products_mapping(std::ostream& out) {
  assert(is_finalized());

  VertexNameMap products_index = boost::get(boost::vertex_name, products_graph);
  VertexNameMap patterns_index = boost::get(boost::vertex_name, patterns_graph);

  bool ok = true;

  for (auto map_it: products_to_patterns_mapping) {
    const Node& prod_node = products_index[map_it.first];
    if (!prod_node.is_mol) {
      continue;
    }
    const MolInstance& prod_mi = *prod_node.mol;

    const Node& pat_node = patterns_index[map_it.second];
    assert(prod_node.is_mol);
    const MolInstance& pat_mi = *pat_node.mol;

    // check that this molecule uses the same components in both directions
    // the components are ordered according to the definition in MolType
    ok = ok && check_components_mapping(pat_mi, prod_mi, "created in reaction rule", out);
    ok = ok && check_components_mapping(prod_mi, pat_mi, "used as pattern in reaction rule", out);

    ok = ok && check_components_states(prod_mi, pat_mi, out);
  }

  return ok;
}


void RxnRule::move_products_that_are_also_reactants_to_be_the_first_products() {
#ifdef MCELL4_SORT_RXN_PRODUCTS_BY_NAME
  // this is a behavior of NFsim that products seem to be sorted by name,
  // let's assume by the first molecule of each complex
  // for now let's sort max 2 products
  if (products.size() == 2) {

    const string& name0 = products[0].to_str(*bng_data, false); // we do not care whether this is a surf or vol rxn
    const string& name1 = products[1].to_str(*bng_data, false);

    if (name0 > name1) {
      CplxInstance tmp = products[0];
      products[0] = products[1];
      products[1] = tmp;
    }

    // then we need to recompute the products_graph
    create_products_graph();

    // update mapping
    compute_reactants_products_mapping();
  }
  return;
#endif

  // for each reactant (from the end since we want the products to be ordered in the same way)
  for (int pi = products.size() - 1; pi > 0; pi--) {
    uint ri;
    bool found = get_assigned_simple_cplx_reactant_for_product(pi, ri);

    if (found) {
      // move product to the front
      CplxInstance prod = products[pi];
      products.erase(products.begin() + pi);
      products.insert(products.begin(), prod);

      // this swap does not seem to to call the CplxInstancer copy ctor, so
      // we must call graph update manually...
      // TODO: this is a bit weird, needs to be checked better, e.g. test mdl/2200 failed without this update
      for (size_t pi_update = 0; pi_update < products.size(); pi_update++) {
        products[pi_update].create_graph();
      }
      // then we need to recompute the products_graph
      create_products_graph();

      // update mapping
      compute_reactants_products_mapping();
    }
  }
}

// returns -1 if the species cannot be a reactant
// returns 0 if this is an unimol rxn or it is the first reactant of a bimol reaction
// returns 1 if this is the second reactant of a bimol reaction
int RxnRule::get_reactant_index(const CplxInstance& inst, const SpeciesContainer& all_species) {
  for (size_t i = 0; i < reactants.size(); i++) {
    const CplxInstance& reactant = reactants[i];
    if (inst.matches_pattern(reactant, true)) { // reactant is the pattern to be matched
      return i;
    }
  }
  return -1;
}


bool RxnRule::species_can_be_reactant(const species_id_t id, const SpeciesContainer& all_species) {

  // check caches first
  if (species_applicable_as_reactants.count(id) != 0) {
    return true;
  }
  if (species_not_applicable_as_reactants.count(id) != 0) {
    return false;
  }

  // need to check whether the complex instance can be a reactant
  const CplxInstance& inst = all_species.get_as_cplx_instance(id);

  int index = get_reactant_index(inst, all_species);
  assert(index >= -1 && index <= 1);
  bool matches = index != -1;

  if (matches) {
    species_applicable_as_reactants.insert_unique(id);
  }
  else {
    species_not_applicable_as_reactants.insert_unique(id);
  }

  return matches;
}


bool RxnRule::species_can_be_bimol_reactants(
    const species_id_t id1, const species_id_t id2, const SpeciesContainer& all_species
) {
  assert(is_bimol());

  // check whether either of the species can be reactant
  if (!species_can_be_reactant(id1, all_species)) {
    return false;
  }
  if (!species_can_be_reactant(id2, all_species)) {
    return false;
  }

  // need to find out whether both can be used at the same time
  const CplxInstance& inst1 = all_species.get_as_cplx_instance(id1);
  const CplxInstance& inst2 = all_species.get_as_cplx_instance(id2);

  // which reactants match id1?
  bool id1_matches[2] = {false, false};
  bool id2_matches[2] = {false, false};
  for (size_t i = 0; i < reactants.size(); i++) {
    id1_matches[i] = inst1.matches_pattern(reactants[i], true);
    id2_matches[i] = inst2.matches_pattern(reactants[i], true);
  }

  // there must be direct or crossed covering of the reactants, i.e.
  // id1 -> patA && id2 -> patB or
  // id1 -> patB && id2 -> patA
  return (id1_matches[0] && id2_matches[1]) || (id1_matches[1] && id2_matches[0]);
}

// returns true when the species matches both reactants
bool RxnRule::species_is_both_bimol_reactants(const species_id_t id, const SpeciesContainer& all_species) {

  if (!is_bimol()) {
    return false;
  }

  // check if the species can be a reactant at all
  if (!species_can_be_reactant(id, all_species)) {
    return false;
  }

  const CplxInstance& inst = all_species.get_as_cplx_instance(id);
  return
      inst.matches_pattern(reactants[0], true) &&
      inst.matches_pattern(reactants[1], true);
}


bool RxnRule::update_variable_rxn_rate(const float_t current_time, const RxnClass* requester) {
  if (!may_update_rxn_rate()) {
    return false;
  }
  assert(!base_variable_rates.empty());
  assert(next_variable_rate_index < base_variable_rates.size());

  if (base_variable_rates[next_variable_rate_index].time > current_time
      && !cmp_eq(base_variable_rates[next_variable_rate_index].time, current_time) ) {
    return false;
  }

  // find which time to use - the highest but still smaller than the following one
  size_t current_index = next_variable_rate_index;
  while (current_index < base_variable_rates.size() &&
          (current_time > base_variable_rates[current_index + 1].time ||
           cmp_eq(current_time, base_variable_rates[current_index + 1].time)
        )
  ) {
    current_index++;
  }

  // should we use the last entry?
  if (current_index == base_variable_rates.size() &&
      next_variable_rate_index == current_index - 1) {
    current_index = next_variable_rate_index;
  }

  // current_time >= time for next change
  base_rate_constant = base_variable_rates[current_index].rate_constant;
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


std::string RxnRule::to_str(const bool with_rate_constant, const bool with_name) const {
  stringstream ss;
  if (with_name) {
    ss << name << ": ";
  }

  ss << complex_instance_vector_to_str(reactants);
  ss << " -> ";
  ss << complex_instance_vector_to_str(products);

  if (with_rate_constant) {
    ss << " " << base_rate_constant;
  }

  ss << " (";
  ss << "id: " << id << ", ";
  ss << "pat on pat mappings: " << num_patterns_onto_patterns_mapping;
  ss << ")";

  return ss.str();
}


std::string RxnRule::complex_instance_vector_to_str(const CplxInstanceVector& complexes) const {
  stringstream ss;
  for (size_t i = 0; i < complexes.size(); i++) {
    ss << complexes[i].to_str(*bng_data, is_surf_rxn());

    if (i != complexes.size() - 1) {
      ss << " + ";
    }
  }
  return ss.str();
}


std::string RxnRule::reactants_to_str() const {
  return complex_instance_vector_to_str(reactants);
}


std::string RxnRule::products_to_str() const {
  return complex_instance_vector_to_str(products);
}


void RxnRule::dump_complex_instance_vector(
    const CplxInstanceVector& complexes,
    const std::string ind) const {

  for (size_t i = 0; i < complexes.size(); i++) {
    cout << ind << "CplxInstance " << i << ":\n";
    complexes[i].dump(true, ind + "  ");

    if (!is_simple()) {
      cout << ind + "  " << "  graph:\n";
      dump_graph(complexes[i].get_graph(), bng_data, ind + "  ");
    }
    cout << "\n";
  }
}


void RxnRule::dump(const bool for_diff, const std::string ind) const {
  if (!for_diff) {
    cout << ind << to_str();
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

    cout << ind << "base_rate_constant: " << base_rate_constant << "\n";
    cout << ind << "variable_rates.size: " << base_variable_rates.size() << "\n";

    if (!base_variable_rates.empty()) {
      for (size_t i = 0; i < base_variable_rates.size(); i++) {
        cout << ind + "  " << "t: " << base_variable_rates[i].time << ", r: " << base_variable_rates[i].rate_constant << "\n";
      }
    }

    cout << ind << "mol_instances_are_fully_maintained: " << mol_instances_are_fully_maintained << "\n";

    cout << ind << "reactants:\n";
    dump_complex_instance_vector(reactants, ind);
    cout << ind << "products:\n";
    dump_complex_instance_vector(products, ind);
  }
}


} /* namespace BNG */
