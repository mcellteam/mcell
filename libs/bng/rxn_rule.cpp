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
#include "bng/bngl_names.h"
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


static void set_graph_reactant_pattern_indices(Graph& g, const uint index_to_set) {
  VertexNameMap graph_index = boost::get(boost::vertex_name, g);
  typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter; // TODO: get rid of these typedefs
  std::pair<vertex_iter, vertex_iter> it;
  for (it = boost::vertices(g); it.first != it.second; ++it.first) {
    vertex_descriptor_t desc = *it.first;
    Node& node = graph_index[desc];
    node.reactant_pattern_index = index_to_set;
  }
}


static bool less_product_cplxs_by_rxn_rule_index(
    const ProductSpeciesPtrWIndices& a1, const ProductSpeciesPtrWIndices& a2
) {
  assert(!a1.rule_product_indices.empty() && !a2.rule_product_indices.empty());

  // each product description has a set of rxn product indices,
  // sort by the smallest of them
  return *a1.rule_product_indices.begin() < *a2.rule_product_indices.begin();
}


// appends found pathways into the pathways array
void RxnRule::define_rxn_pathways_for_specific_reactants(
    SpeciesContainer& all_species,
    const BNGConfig& bng_config,
    const species_id_t reactant_a_species_id,
    const species_id_t reactant_b_species_id,
    const float_t pb_factor,
    RxnClassPathwayVector& pathways
) {

  if (is_simple()) {
    RxnProductsVector product_species;
    for (size_t i = 0; i < products.size(); i++) {
      const Cplx& product = products[i];

      // simple product is deterministic
      // simple species are defined mcell3 mode but in BNG mode they
      // may not have been created (they are based on molecule types)
      species_id_t species_id = all_species.find_full_match(product);
      if (species_id == SPECIES_ID_INVALID) {
        // no need to make simple species removable
        Species* new_species = new Species(product, *bng_data, bng_config);
        species_id = all_species.find_or_add_delete_if_exist(new_species);
        assert(species_id != SPECIES_ID_INVALID);
      }
      product_species.push_back(ProductSpeciesIdWIndices(species_id, i));
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

    // prepare set of reactants
    vector<species_id_t> reactant_species;
    reactant_species.push_back(reactant_a_species_id);
    if (is_bimol()) {
      reactant_species.push_back(reactant_b_species_id);
    }

    // we might get multiple matches on reactant(s), the numbers
    // of matches multiplied give us the total number of variants
    // a single random number then is used to choose a single variant
    ProductSetsVector product_sets;
    create_products_for_complex_rxn(
        all_species,
        bng_config,
        reactant_species,
        pb_factor,
        pathways
    );

  }
}


void RxnRule::create_patterns_graph() {
  // mark nodes of reactants - we must be able to distinguish them when
  // applying this rule
  for (size_t i = 0; i < reactants.size(); i++) {
    Graph& graph = reactants[i].get_graph();
    set_graph_reactant_pattern_indices(graph, i);
  }

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


void RxnRule::finalize() {
  bool simple = true;

  // finalize all reactants and products
  for (Cplx& ci: reactants) {
    ci.finalize_cplx();
    simple = simple && ci.is_simple();
  }

  for (Cplx& ci: products) {
    ci.finalize_cplx();
    simple = simple && ci.is_simple();
  }

  if (simple) {
    set_flag(RXN_FLAG_SIMPLE);
  }

  create_patterns_graph();
  create_products_graph();
  compute_reactants_products_mapping();

  // for MCell3 compatibility, updates mapping when needed
  move_products_that_are_also_reactants_to_be_the_first_products();

  // set flag that tells us whether we have to do equivalence checks
  // when constructing sets of possible products
  if (matching_may_produce_multiple_identical_results()) {
    set_flag(RXN_FLAG_MAY_PRODUCE_MUTLIPLE_IDENTICAL_PRODUCTS);
  }

  set_finalized();
}


bool RxnRule::matching_may_produce_multiple_identical_results() const {

  Graph patterns_graph_no_indices = patterns_graph;

  // we must not have pattern indices when checking for symmetry
  set_graph_reactant_pattern_indices(patterns_graph_no_indices, INDEX_INVALID);

  // are there multiple matches of the patterns graph onto itself?
  VertexMappingVector mappings;
  get_subgraph_isomorphism_mappings(
      patterns_graph_no_indices,
      patterns_graph_no_indices,
      false, // do not stop with first match
      mappings
  );
  assert(mappings.size() >= 1);
  if (mappings.size() > 1) {
    return true;
  }

  // now check whether we can match multiple components in the same molecule
  return may_modify_more_than_one_identical_component();
}


bool RxnRule::may_modify_more_than_one_identical_component() const {

  for (const Cplx& reactant: reactants) {
    for (const ElemMol& mi: reactant.elem_mols) {
      const ElemMolType& mt = bng_data->get_elem_mol_type(mi.elem_mol_type_id);

      // for each component of the pattern
      for (const Component& ci: mi.components) {
        // how many times is this component type used in the molecule type template
        if (mt.get_component_uses_count(ci.component_type_id) > 1) {
          return true;
        }
      }
    }
  }

  return false;
}


struct MolCompInfo {
  MolCompInfo(const vertex_descriptor_t desc_, const ElemMol* em_)
    : desc(desc_), already_matched(false), em(em_), comp(nullptr) {
  }
  MolCompInfo(const vertex_descriptor_t desc_, const Component* comp_)
    : desc(desc_), already_matched(false), em(nullptr), comp(comp_) {
  }

  const ElemMol* get_em() const {
    assert(em != nullptr);
    return em;
  }

  const Component* get_comp() const {
    assert(comp != nullptr);
    return comp;
  }

  // graph node corresponding to this molecule or component
  vertex_descriptor_t desc;

  // if this is a node from a pattern this array contains matching score
  // from this pattern node onto ith product node (i is index to this array)
  vector<int> matching_score;

  // set to true if this pattern node has its product node assigned
  bool already_matched;

  // if this is an elementary molecule, this array contains information
  // to which elementary molecule types are we bound
  set<elem_mol_type_id_t> connected_em_types;

private:
  const ElemMol* em;
  const Component* comp;
};


// - vertex pointed to by desc must be a component and must have a single bond
//   to another component and also a bond to its molecule instance
// - returns -1 when the component is not connected to another component
// - asserts in debug mode when must_exist is set to true and target was not found
//
// - the release_assert calls are here because boost's library produced
//   weird errors with release build, they are cheap anyway so they are kept there
const vertex_descriptor_t TARGET_NOT_FOUND = -1;
static vertex_descriptor_t get_bond_target(
    Graph& graph,
    const vertex_descriptor_t desc,
    const bool must_exist = true,
    const bool looking_for_component = true
) {
  #ifdef DEBUG_CPLX_RXNS
    cout << "get_bond_target:\n";
    dump_graph(graph);
    cout << "desc:" << (int)desc << "\n";
  #endif

  bool comp_found = false;

  vertex_descriptor_t res_em = TARGET_NOT_FOUND;
  vertex_descriptor_t res_comp = TARGET_NOT_FOUND;
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
      release_assert(res_em == TARGET_NOT_FOUND);
      res_em = n_desc;
    }
    else {
      // TODO, see test tests/bngl/0967_vol_in_a_plus_surf_b_to_surf_c_cplx_reusing_mol_in_prod
      release_assert(res_comp == TARGET_NOT_FOUND &&
          "Used rule is not supported yet, probably due to finding only a partial match from reactant to product such as sb(b!1).sy(y!1) -> sc(c!1).sy(y!1)");
      res_comp = n_desc;
    }
  }

  release_assert(res_em != TARGET_NOT_FOUND && "Component must be connected to its molecule");
  release_assert(!must_exist || res_comp != TARGET_NOT_FOUND);

  return (looking_for_component) ? res_comp : res_em;
}


// inserts all elementary moelcules connected to our em_desc through bonds
static void get_all_connected_elem_mol_types(
    Graph& graph, const vertex_descriptor_t em_desc, set<elem_mol_type_id_t>& connected_em_types) {

  VertexNameMap index = boost::get(boost::vertex_name, graph);
  assert(index[em_desc].is_mol);

  #ifdef DEBUG_CPLX_RXNS
    cout << "get_all_connected_elem_mol_types:\n";
    dump_graph(graph);
  #endif

  // for each connected component
  boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
  for (boost::tie(ei,edge_end) = boost::out_edges(em_desc, graph); ei != edge_end; ++ei) {
    edge_descriptor_t e_desc = *ei;

    // component of our elem mol
    vertex_descriptor_t comp_desc = boost::target(e_desc, graph);
    assert(!index[comp_desc].is_mol);

    // component to which 'comp' is bound
    vertex_descriptor_t connected_comp_desc = get_bond_target(graph, comp_desc, false, true);

    if (connected_comp_desc == TARGET_NOT_FOUND) {
      // component has no bond
      continue;
    }

    // finally get the elem mol
    assert(!index[connected_comp_desc].is_mol);
    vertex_descriptor_t connected_em_desc = get_bond_target(graph, connected_comp_desc, true, false);
    assert(connected_em_desc != TARGET_NOT_FOUND);
    const Node& em = index[connected_em_desc];
    assert(em.is_mol);

    connected_em_types.insert(em.mol->elem_mol_type_id);
  }
}


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
      MolCompInfo info(desc, mi_node.mol);

      // collect connected elementary molecules
      get_all_connected_elem_mol_types(graph, desc, info.connected_em_types);

      res.push_back(info);
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
    const Component& ci1,
    const Component& ci2
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


bool are_replaceable_elem_mols_do_not_check_connections(
    const BNGData& bng_data, const ElemMol& em1, const ElemMol& em2) {

  // components must have the same names and states
  const ElemMolType& emt1 = bng_data.get_elem_mol_type(em1.elem_mol_type_id);
  const ElemMolType& emt2 = bng_data.get_elem_mol_type(em2.elem_mol_type_id);

  if (emt1.component_type_ids.size() != emt2.component_type_ids.size()) {
    return false;
  }

  // same names, in the same order, and the components have the same states
  for (size_t i = 0; i < emt1.component_type_ids.size(); i++) {
    const ComponentType& ct1 = bng_data.get_component_type(emt1.component_type_ids[i]);
    const ComponentType& ct2 = bng_data.get_component_type(emt2.component_type_ids[i]);

    if (ct1.name != ct2.name) {
      return false;
    }

    if (ct1.allowed_state_ids != ct2.allowed_state_ids) {
      return false;
    }
  }

  return true;
}


bool are_replaceable_elem_mols(const BNGData& bng_data, const MolCompInfo& em1_info, const MolCompInfo& em2_info) {
  // special case where we can replace e.g. A_vol(x, y) with A_surf(x, y) in rxn
  // B.A_vol -> B.A_surf

  const ElemMol& em1 = *em1_info.get_em();
  const ElemMol& em2 = *em2_info.get_em();

  // both molecules must be connected and must be connected to the same molecules
  if (!em1.has_bond() || !em2.has_bond()) {
    return false;
  }

  // TODO: maybe allow just change from surf to vol and vice versa?

  // we much check connections so that for instance elem mols in this rule are not matched:
  // sb(b!1).sy(y!1) -> sc(c!1).sz(y!1), we allow only
  // sb(b!1).sy(y!1) -> sb(c!1).sz(y!1)
  // might need to be relaxed if we would need to change type of multiple
  // elem mols in one rule that are connected to each other but for now let's require that the
  // surroundings are exactly the same
  if (em1_info.connected_em_types != em2_info.connected_em_types) {
    return false;
  }

  return are_replaceable_elem_mols_do_not_check_connections(bng_data, em1, em2);
}


static int get_elem_mol_matching_score(const BNGData& bng_data, const MolCompInfo& pat_info, const MolCompInfo& prod_info) {

  const ElemMol& pat = *pat_info.get_em();
  const ElemMol& prod = *prod_info.get_em();

  // elementary molecule type must be the same or must be replaceable
  if (pat.elem_mol_type_id != prod.elem_mol_type_id &&
      !are_replaceable_elem_mols(bng_data, pat_info, prod_info)) {
    return -1;
  }

  int res = 0;

  // lot of points if the elem mol type is the same
  if (pat.elem_mol_type_id == prod.elem_mol_type_id) {
    res += 10;
  }

  // add a lot of points if we have the same number of components
  if (pat.components.size() == prod.components.size()) {
    res += pat.components.size();
  }

  // and add points for each matching component, they are ordered according to
  // the MolType
  // TODO: probably use the same approach for matching components that is done later
  for (size_t i = 0; i < pat.components.size(); i++) {
    const Component& pat_compi = pat.components[i];

    if (i < prod.components.size()) {
      res += get_component_instance_matching_score(pat_compi, prod.components[i]);
    }
  }
  
  // add extra point if compartment is the same 
  if (pat.compartment_id == prod.compartment_id) {
    res++;
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


void find_best_product_to_pattern_mapping(
    const BNGData& bng_data,
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
      pat.matching_score[i] = get_elem_mol_matching_score(bng_data, pat, prod);
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
        pat.matching_score[i] = get_component_instance_matching_score(*pat.get_comp(), *prod.get_comp());
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


static void change_elem_mol_type_and_component_types(
    const BNGData& bng_data, ElemMol& reac_em, elem_mol_type_id_t target_elem_mol_type_id) {

  elem_mol_type_id_t orig_elem_mol_type_id = reac_em.elem_mol_type_id;
  // no need to check graph components because graph contains pointers to components owned by reac_em
  // and we are not changing any bonds
  reac_em.elem_mol_type_id = target_elem_mol_type_id;

  const ElemMolType& orig_emt = bng_data.get_elem_mol_type(orig_elem_mol_type_id);
  const ElemMolType& target_emt = bng_data.get_elem_mol_type(target_elem_mol_type_id);
  release_assert(reac_em.components.size() == orig_emt.component_type_ids.size());
  release_assert(reac_em.components.size() == target_emt.component_type_ids.size());
  for (size_t i = 0; i < reac_em.components.size(); i++) {
    Component& comp = reac_em.components[i];
    component_type_id_t target_comp_id = target_emt.component_type_ids[i];
    // debug check that component name matches (they must be ordered in the same way)
    assert(bng_data.get_component_type(comp.component_type_id).name ==
        bng_data.get_component_type(target_comp_id).name);

    comp.component_type_id = target_comp_id;
  }
}


// goes through all reactant elem mols and components in pattern_reactant_mapping
// and sets their modified_ordering_index
static void set_modified_ordering_index_for_all_reactants_matching_pattern(
  Graph& reactants_graph,
  const VertexMapping& pattern_reactant_mapping
) {
  VertexNameMap reactants_index = boost::get(boost::vertex_name, reactants_graph);

  for (auto pat_reac_pair: pattern_reactant_mapping) {
    vertex_descriptor_t reac_desc = pat_reac_pair.second;
    Node& reac_node = reactants_index[reac_desc];
    reac_node.modified_ordering_index = reac_node.ordering_index;
  }
}


// manipulate the reactants_graph according to how
// products should look like
static void apply_rxn_on_reactants_graph(
    const BNGData& bng_data,
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
  for (const auto& prod_pat_pair: product_pattern_mapping) {
    vertex_descriptor_t prod_desc = prod_pat_pair.first;
    const Node& prod_node = products_index[prod_pat_pair.first];
    // product->pattern
    vertex_descriptor_t pat_desc = prod_pat_pair.second;

    // find component corresponding to reactant pattern in reactants_graph
    //   this mapping must exist because we matched the pattern graph to reactants graph
    auto pat_reac_it = pattern_reactant_mapping.find(pat_desc);
    assert(pat_reac_it != pattern_reactant_mapping.end() && "Mapping must exist");

    vertex_descriptor_t reac_desc = pat_reac_it->second;

    // manipulate state and or bond
    Node& reac_node = reactants_index[reac_desc];

    // set product index marker
    assert(prod_node.product_index != INDEX_INVALID);
    reac_node.product_index = prod_node.product_index;

    if (!prod_node.is_mol) {
      assert(!reac_node.is_mol);
      const Component& prod_ci = *prod_node.component;
      Component& reac_ci = *reac_node.component;

      // update state
      if (prod_ci.state_is_set() && prod_ci.state_id != reac_ci.state_id) {
        reac_ci.state_id = prod_ci.state_id;
        reac_node.modified_ordering_index = reac_node.ordering_index;
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
            reac_node.modified_ordering_index = reac_node.ordering_index;
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

            reac_node.modified_ordering_index = reac_node.ordering_index;
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

            reac_node.modified_ordering_index = reac_node.ordering_index;
          }
          // new: !2
          else if (prod_ci.bond_has_numeric_value()) {
            assert(prod_ci.bond_value != reac_ci.bond_value);
            // remove original one
            vertex_descriptor_t orig_target_desc = get_bond_target(reactants_graph, reac_desc);
            bonds_to_remove.insert(UnorderedPair(
                reac_desc,
                orig_target_desc
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

            // a change occurs only when the target is different
            if (orig_target_desc != target_reac_desc) {
              reac_node.modified_ordering_index = reac_node.ordering_index;
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
    else {
      // elementary molecule (mol)
      assert(prod_node.mol != nullptr);
      const ElemMol& prod_em = *prod_node.mol;

      assert(reac_node.mol != nullptr);
      ElemMol& reac_em = *reac_node.mol;

      // in rules such as A(x!1).B(a!1) -> A(x!1).C(a!1),
      // the elem mol id might differ, inserted into pattern->product mapping by using
      // are_replaceable_elem_mols
      if (prod_em.elem_mol_type_id != reac_em.elem_mol_type_id) {
        // we need to change the elem mol type and also types of all components
        change_elem_mol_type_and_component_types(bng_data, reac_em, prod_em.elem_mol_type_id);

        // - also remember that we made a change here so that different products are not matched,
        //   see test bngl/0055 - rule A(bc!1).B(a!1) -> A(bc!1).C(a!1) is applied onto
        //   A(x,bc!2).B(a!2,a!3).A(bc!3,x) -> there must be two distinct products because reactant
        //   bonds 2 and 3 can be matched
        // - setting the reac_node.modified_ordering_index does not help because we won't get different
        //   results, we must mark e.g. the component whose bond changed, it is the one that is mapped
        //   by the pattern
        // - a solution here is to mark all reactant elementary molecules and components matched by the
        //   rule
        set_modified_ordering_index_for_all_reactants_matching_pattern(
            reactants_graph, pattern_reactant_mapping);
      }

      // update compartment if needed
      if (is_specific_compartment_id(prod_em.compartment_id)) {
        reac_em.compartment_id = prod_em.compartment_id;
      }
    }
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
    Species* species,
    std::set<uint>& product_indices
) {
  product_indices.clear();

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

    // check molecule whether it has a marker that says which product it was,
    // we need to maintain the ordering of products
    if (mol.product_index != INDEX_INVALID) {
      product_indices.insert(mol.product_index);
    }

    if (!mol.used_in_rxn_product) {
      // this is a consumed reactant
      return false;
    }

    if (!mol.is_mol) {
      continue;
    }

    species->elem_mols.push_back(*mol.mol);
    ElemMol& mi = species->elem_mols.back();

    // we will recreate components because they might have changed
    mi.components.clear();

    // for each of its components
    boost::graph_traits<Graph>::out_edge_iterator ei, edge_end;
    for (boost::tie(ei,edge_end) = boost::out_edges(mol_desc, graph); ei != edge_end; ++ei) {
      Graph::edge_descriptor e_mol_comp = *ei;
      Graph::vertex_descriptor comp_desc = boost::target(e_mol_comp, graph);

      const Node& comp = index[comp_desc];
      assert(!comp.is_mol && "Only a component may be connected to a molecule.");

      mi.components.push_back(*comp.component); // we use state as it its
      Component& compi = mi.components.back();

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

  return true;
}


// product compartment may be set only partially in the rxn rule, this sets them accordingly
void RxnRule::set_product_compartments(Species* product_species, std::set<uint>& product_indices) const {
  // first we need to sef surf/vol flags
  product_species->finalize_cplx();

  bool is_surf = product_species->is_surf();

  if (!is_surf) {
    // we might need to fix volume compartment, e.g. from s(c!1)@PM.c(s!1,c!2)@EC.c(s,c!2)@EC -> s(c)@PM + c(s)@CP
    // we get product c(s!1,c!2)@CP.c(s,c!2)@EC, we use the specific compartment in CP to set the compartment of all
    // elem mols

    // get compartment of the product if set
    compartment_id_t vol_prod_comp_id = COMPARTMENT_ID_NONE;
    for (uint i: product_indices) {
      compartment_id_t curr_prod_comp_id = products[i].get_primary_compartment_id();
      release_assert(vol_prod_comp_id == COMPARTMENT_ID_NONE || vol_prod_comp_id == curr_prod_comp_id);
      vol_prod_comp_id = curr_prod_comp_id;
    }

    if (is_specific_compartment_id(vol_prod_comp_id)) {
      product_species->set_compartment_id(vol_prod_comp_id);
    }
  }

  compartment_id_t primary_compartment_id = product_species->get_primary_compartment_id();

  // we will be setting compartments only when the primary compartment is set, otherwise the
  // compartment is set when a molecule is added to a partition
  // IN/OUT compartments are allowed but ignored because they were replaced with orientation
  if (primary_compartment_id != COMPARTMENT_ID_NONE && !is_in_out_compartment_id(primary_compartment_id)) {

    for (ElemMol& em: product_species->elem_mols) {
      if (is_surf) {
        if (em.is_surf()) {
          // all surface elementary molecules must get the primary surface compartment
          em.compartment_id = primary_compartment_id;
        }
        else {
          // compartments for volume elem mols must be set
          // TODO: there should be a check for this e.g. during initialization
          release_assert(em.compartment_id != COMPARTMENT_ID_NONE &&
              "Rxn product sets surface compartment but volume compartment is unset, this is not supported because "
              "it is not clear which compartment should the volume elementary molecule use, please check your reaction rules");
        }
      }
      else {
        // volume - if set, all elem mols must have the same compartment
        assert(em.is_vol());
        if (em.compartment_id == COMPARTMENT_ID_NONE) {
          em.compartment_id = primary_compartment_id;
        }
        else if (is_specific_compartment_id(em.compartment_id) &&
            !bng_data->get_compartment(em.compartment_id).is_3d) {
          // override 2D compartment by a volume compartment
          em.compartment_id = primary_compartment_id;
        }
        else {
          release_assert(em.compartment_id == primary_compartment_id && "Volume compartment mismatch in rxn product");
        }
      }
    }
  }
}


void RxnRule::create_products_from_reactants_graph(
    const BNGData* bng_data,
    Graph& reactants_graph,
    ProductCplxWIndicesVector& created_products
) const {
  created_products.clear();
  created_products.reserve(2);

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
    Species* product_species = new Species(*bng_data);

    // some reactants may have been removed,
    // also, products still may be connected even though that are disconnected on the right-hand side of the
    // rule - e.g. when rule
    // A(b!1).B(a!1) -> A(b) + B(a)
    // is applied on
    // A(b!1,c!2).B(a!1,c!3).C(a!2,B!3) - a single complex is a result
    std::set<uint> product_indices;
    bool is_rxn_product =
        convert_graph_component_to_product_cplx_inst(
            reactants_graph, graph_components, i, product_species, product_indices);

    if (is_rxn_product) {
      release_assert(!product_indices.empty());

      // to which rxn product this product corresponds?
      set_product_compartments(product_species, product_indices);

      // remember product with its indices
      created_products.push_back(ProductSpeciesPtrWIndices(product_species, product_indices));
    }
    else {
      delete product_species;
    }
  }
}


// goes through all nodes of the graph and sets a unique index to each of them
static void set_ordering_indices(Graph& g) {
  uint ordering_index = 0;
  VertexNameMap index = boost::get(boost::vertex_name, g);
  typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
  pair<vertex_iter, vertex_iter> it;
  for (it = boost::vertices(g); it.first != it.second; ++it.first) {
    Graph::vertex_descriptor desc = *it.first;
    Node& n = index[desc];
    n.ordering_index = ordering_index;
    ordering_index++;
  }
}


static bool is_graph_unique_wrt_modified_ordering(
    Graph& new_graph, vector<Graph>& distinct_product_graphs) {

  for (Graph& g: distinct_product_graphs) {
    VertexMappingVector mappings;
    get_subgraph_isomorphism_mappings(
        g, // existing graph
        new_graph,
        true, // stop with first match
        mappings
    );

    if (!mappings.empty()) {
      // already present in distinct_product_graphs
      return false;
    }
  }

  // not found
  return true;
}


static bool less_pattern_reactant_mappings(
    const VertexMapping& m1, const VertexMapping& m2
) {
  for (auto it1: m1) {
    auto it2 = m2.find(it1.first);
    release_assert(it2 != m2.end() && "All patern nodes must be mapped");

    if (it1.second != it2->second) {
      return it1.second < it2->second;
    }
  }
  // identical mapping is ok, 
  // some std::sort implementations compare the same elements
  return false; 
}


static void sort_mappings(VertexMappingVector& mappings) {
  // used to get different results on Linux and MacOS,
  // need to sort the mappings in some way,
  // there should not be many elements in each of those mappings,
  // so the comparison is not very efficient,
  // NOTE: maybe we will need to optimize this

  std::sort(
      mappings.begin(),
      mappings.end(),
      less_pattern_reactant_mappings
  );
}


void RxnRule::create_products_for_complex_rxn(
    SpeciesContainer& all_species,
    const BNGConfig& bng_config,
    const std::vector<species_id_t>& reactant_species,
    const float_t pb_factor,
    RxnClassPathwayVector& pathways
) {

  // the result of this function is cached in rxnclass
  assert(reactant_species.size() == reactants.size());
  assert(reactant_species.size() == 1 || reactant_species.size() == 2);

  // input reactants are always created from canonicalized species therefore the next time we recreate
  // them they will be the same
  vector<const Cplx*> input_reactants;
  // downcast, Cplx is sufficient for product computation
  input_reactants.push_back(dynamic_cast<const Cplx*>(&all_species.get(reactant_species[0])));
  if (is_bimol()) {
    input_reactants.push_back(dynamic_cast<const Cplx*>(&all_species.get(reactant_species[1])));
  }

  if (bng_config.notifications.bng_verbosity_level >= 1) {
    cout << "Creating products for complex rxn " << name <<
        " for reactant(s) " << input_reactants[0]->to_str();
    if (input_reactants.size() == 2) {
      cout << " + " << input_reactants[1]->to_str() << "\n";
    }
  }

  Graph reactants_graph;
  bool use_symmetric_reactants_graph = false;
  Graph symmetric_reactants_graph; // needed for cases when both patterns can be matched onto both reactants

  if (is_bimol()) {
    // figure out which reactant pattern matches which input reactant
    vector<pair<uint, uint>> reac_indices;
    get_bimol_reactant_indices(
        reactant_species[0], reactant_species[1], all_species,
        reac_indices
    );
    release_assert(reac_indices.size() == 1 || reac_indices.size() == 2);
    release_assert(reac_indices[0].first != reac_indices[0].second);

    // and mark graph coming from each reactant so that during matching it is clear that
    // we did not match one input reactant with both patterns
    reactants_graph = input_reactants[0]->get_graph();
    set_graph_reactant_pattern_indices(reactants_graph, reac_indices[0].first);

    Graph reac1_graph = input_reactants[1]->get_graph();
    set_graph_reactant_pattern_indices(reac1_graph, reac_indices[0].second);
    merge_graphs(reactants_graph, reac1_graph);

    if (reac_indices.size() == 2) {
      // if both patterns match both reactants, we must also evaluate the symmetric variant
      use_symmetric_reactants_graph = true;

      release_assert(reac_indices[1].first != reac_indices[1].second);

      symmetric_reactants_graph = input_reactants[0]->get_graph();
      set_graph_reactant_pattern_indices(symmetric_reactants_graph, reac_indices[1].first);

      Graph symmetric_reac1_graph = input_reactants[1]->get_graph();
      set_graph_reactant_pattern_indices(symmetric_reac1_graph, reac_indices[1].second);
      merge_graphs(symmetric_reactants_graph, symmetric_reac1_graph);
    }
  }
  else {
    reactants_graph = input_reactants[0]->get_graph();
    set_graph_reactant_pattern_indices(reactants_graph, 0); // only one reactant
  }

  // compute mapping reactant pattern -> reactant
  VertexMappingVector pattern_reactant_mappings;
  get_subgraph_isomorphism_mappings(
      patterns_graph, // pattern
      reactants_graph, // actual reactant
      false, // do not stop with first match
      pattern_reactant_mappings
  );

  if (use_symmetric_reactants_graph) {
    // we need to evaluate case when both patterns match both reactants
    VertexMappingVector symmetric_pattern_reactant_mappings;
    get_subgraph_isomorphism_mappings(
        patterns_graph, // pattern
        symmetric_reactants_graph, // actual reactant
        false, // do not stop with first match
        symmetric_pattern_reactant_mappings
    );
    pattern_reactant_mappings.insert(
        pattern_reactant_mappings.end(),
        symmetric_pattern_reactant_mappings.begin(),
        symmetric_pattern_reactant_mappings.end()
    );
  }

  // and append mappings for cases when patterns match both reactants

  assert(pattern_reactant_mappings.size() != 0 &&
      "Did not find a match of patterns onto reaction.");

  release_assert(pattern_reactant_mappings.size() < MAX_PRODUCT_SETS_PER_RXN
      && "Encountered a huge number of potential product sets for a single reaction");

  if (bng_config.notifications.bng_verbosity_level >= 1) {
    cout << "  - found " << pattern_reactant_mappings.size() << " potential products";
    cout.flush();
  }

  // sort the mappings to make sure we get identical results everywhere
  sort_mappings(pattern_reactant_mappings);

  // decide whether only cache the possible products
  if (pattern_reactant_mappings.size() > MAX_IMMEDIATELLY_COMPUTED_PRODUCT_SETS_PER_RXN &&
      !has_flag(RXN_FLAG_MAY_PRODUCE_MUTLIPLE_IDENTICAL_PRODUCTS)) {

    for (const VertexMapping& mapping: pattern_reactant_mappings) {
      // define the products as species
      RxnProductsVector product_species;

      // we cannot sort the products because we do not know what they are
      // FIXME: we should unify this... - changing MAX_IMMEDIATELLY_COMPUTED_PRODUCT_SETS_PER_RXN
      // may change results

      // the probability is divided by the number of mapping of pattern onto pattern
      // because so many more product sets we will get
      float_t prob = get_rate_constant() * pb_factor;

      pathways.push_back(RxnClassPathway(id, prob, mapping));
    }
    return;
  }

  vector<vector<Cplx>> input_reactants_copies;
  vector<Graph> distinct_product_graphs;
#ifndef NDEBUG
  vector<Graph> debug_distinct_product_graphs;
#endif

  // now, for each of the mappings, compute what different products we might get
  for (const VertexMapping& mapping: pattern_reactant_mappings) {


  #ifdef DEBUG_CPLX_MATCHING
    cout << "Products:\n";
    dump_graph(products_graph);
  #endif

    // we need to make a copy of the reactants because we will be modifying them
    // a new graph will have its ordering indices cleared
    input_reactants_copies.push_back(vector<Cplx>());
    vector<Cplx>& input_reactants_copy = input_reactants_copies.back();
    for (const Cplx* ci: input_reactants) {
      input_reactants_copy.push_back(*ci);
    }
    Graph reactants_graph_copy = input_reactants_copy[0].get_graph();
    if (input_reactants_copy.size() == 2) {
      merge_graphs(reactants_graph_copy, input_reactants_copy[1].get_graph());
    }

    set_ordering_indices(reactants_graph_copy);

    // manipulate nodes using information about products
    apply_rxn_on_reactants_graph(
        *bng_data,
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

    if (has_flag(RXN_FLAG_MAY_PRODUCE_MUTLIPLE_IDENTICAL_PRODUCTS)) {
      // we must verify that we don't have this product yet,
      // this is quite time expensive
      // TODO: some optimization is needed
      if (is_graph_unique_wrt_modified_ordering(reactants_graph_copy, distinct_product_graphs)) {
        distinct_product_graphs.push_back(reactants_graph_copy);
      }
    }
    else {
      // each match is a unique product because the patterns are not symmetrical and
      // there are no multiple components of the same name
      distinct_product_graphs.push_back(reactants_graph_copy);

    #ifndef NDEBUG
      // the assumption above should be ok but to be sure let's check it
      if (is_graph_unique_wrt_modified_ordering(reactants_graph_copy, debug_distinct_product_graphs)) {
        debug_distinct_product_graphs.push_back(reactants_graph_copy);
      }
    #endif
    }
  }

#ifndef NDEBUG
  if (!has_flag(RXN_FLAG_MAY_PRODUCE_MUTLIPLE_IDENTICAL_PRODUCTS)) {
    assert(distinct_product_graphs.size() == debug_distinct_product_graphs.size());
  }
#endif

  if (bng_config.notifications.bng_verbosity_level >= 1) {
    cout << ", of it " << distinct_product_graphs.size() << " unique products\n";
  }

  ProductSetsVector created_product_sets;
  for (Graph& product_graph: distinct_product_graphs) {
    // and finally create products, each disconnected graph in the result is a
    // separate complex

    created_product_sets.push_back(ProductCplxWIndicesVector());
    create_products_from_reactants_graph(bng_data, product_graph, created_product_sets.back());

  #ifdef DEBUG_CPLX_MATCHING
    cout << "Resulting products:\n";
    for (auto& c: created_product_sets.back()) {
      c.product_species->dump();
      cout << "\n";
    }
  #endif

  }

  // and convert resulting complexes represented by Species* to species ids
  for (ProductCplxWIndicesVector& product_cplxs: created_product_sets) {
    // define the products as species
    RxnProductsVector product_species;

    // sort products by the rxn rule product indices (if applicable)
    sort(product_cplxs.begin(), product_cplxs.end(), less_product_cplxs_by_rxn_rule_index);

    // iterating over map sorted by product indices
    for (ProductSpeciesPtrWIndices& product_w_indices: product_cplxs) {
      // need to transform cplx into species id, the possibly new species will be removable
      product_w_indices.product_species->finalize_species(bng_config);
      species_id_t species_id = all_species.find_or_add_delete_if_exist(
          product_w_indices.product_species, true);

      assert(species_id != SPECIES_ID_INVALID);
      product_species.push_back(
          ProductSpeciesIdWIndices(species_id, product_w_indices.rule_product_indices));
    }

    assert(pb_factor != 0);
    assert(!cmp_eq(get_rate_constant(), FLT_GIGANTIC));

    // the probability is divided by the number of mapping of pattern onto pattern
    // because so many more product sets we will get
    float_t prob = get_rate_constant() * pb_factor;

    pathways.push_back(RxnClassPathway(id, prob, product_species));
  }
}


// TODO: very similar code as in the function above, must be unified
void RxnRule::define_rxn_pathway_using_mapping(
  SpeciesContainer& all_species,
  const BNGConfig& bng_config,
  const std::vector<species_id_t>& reactant_species,
  RxnClassPathway& pathway
) const {
  assert(!pathway.products_are_defined);
  assert(!pathway.rule_mapping_onto_reactants.empty());

  // we need to make a copy of the reactants because we will be modifying them
  // a new graph will have its ordering indices cleared
  vector<Cplx> input_reactants_copy;
  input_reactants_copy.reserve(2);
  for (species_id_t s_id: reactant_species) {
    input_reactants_copy.push_back(all_species.get(s_id));
  }

  // and create a graph from them
  // because we creating it from canonical species, we know for sure that this is the same
  // graph as was used when mapping was computed
  Graph reactants_graph = input_reactants_copy[0].get_graph();
  if (input_reactants_copy.size() == 2) {
    merge_graphs(reactants_graph, input_reactants_copy[1].get_graph());
  }

  set_ordering_indices(reactants_graph);

  // manipulate nodes using information about products and the precomputed mapping
  apply_rxn_on_reactants_graph(
      *bng_data,
      reactants_graph,
      pathway.rule_mapping_onto_reactants,
      patterns_graph,
      products_to_patterns_mapping,
      products_graph
  );

  // convert graph with products into cplx instances
  ProductCplxWIndicesVector product_cplxs;
  create_products_from_reactants_graph(bng_data, reactants_graph, product_cplxs);

  // sort products by the rxn rule product indices (if applicable)
  sort(product_cplxs.begin(), product_cplxs.end(), less_product_cplxs_by_rxn_rule_index);

  // and cplx instances into species
  // we are not setting the resulting compartment, neither orientation
  for (ProductSpeciesPtrWIndices& product_w_indices: product_cplxs) {
    // need to transform cplx into species id, the possibly new species will be removable
    product_w_indices.product_species->finalize_species(bng_config);
    species_id_t species_id = all_species.find_or_add_delete_if_exist(
        product_w_indices.product_species, true);

    assert(species_id != SPECIES_ID_INVALID);
    pathway.product_species_w_indices.push_back(
        ProductSpeciesIdWIndices(species_id, product_w_indices.rule_product_indices)
    );
  }

  pathway.rule_mapping_onto_reactants.clear();
  pathway.products_are_defined = true;
}


bool RxnRule::is_simple_cplx_reactant_on_both_sides_of_rxn_w_identical_compartments(const uint index) const {
#ifdef MCELL4_DO_NOT_REUSE_REACTANT
  return false;
#endif
  assert(is_finalized());
  for (const CplxIndexPair& cplx_index_pair: pat_prod_cplx_mapping) {
    if (index == cplx_index_pair.reactant_index && cplx_index_pair.is_simple_mapping) {
      return true;
    }
  }
  return false;
}


// ignores compartments
bool RxnRule::is_simple_cplx_product_on_both_sides_of_rxn_w_identical_compartments(const uint index) const {
#ifdef MCELL4_DO_NOT_REUSE_REACTANT
  return false;
#endif
  assert(is_finalized());
  for (const CplxIndexPair& cplx_index_pair: pat_prod_cplx_mapping) {
    if (index == cplx_index_pair.product_index && cplx_index_pair.is_simple_mapping) {
      return true;
    }
  }
  return false;
}


bool RxnRule::get_assigned_cplx_reactant_for_product(
    const uint product_index, const bool only_simple, uint& reactant_index) const {
  // this is not a time critical search
  for (const CplxIndexPair& cplx_index_pair: pat_prod_cplx_mapping) {
    if (only_simple && !cplx_index_pair.is_simple_mapping) {
      continue;
    }

    if (product_index == cplx_index_pair.product_index) {
      reactant_index = cplx_index_pair.reactant_index;
      return true;
    }
  }
  return false;
}


static size_t find_mol_instance_with_address(const CplxVector& cplx_vex, const ElemMol* mi_addr) {
  for (size_t i = 0; i < cplx_vex.size(); i++) {
    const Cplx& ci = cplx_vex[i];
    for (size_t k = 0; k < ci.elem_mols.size(); k++) {
      if (mi_addr == &ci.elem_mols[k]) {
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
      *bng_data,
      products_graph,
      patterns_graph,
      products_to_patterns_mapping
  );

#ifdef DEBUG_CPLX_MATCHING
  dump_graph_mapping(products_to_patterns_mapping);
#endif

  VertexNameMap patterns_index = boost::get(boost::vertex_name, patterns_graph);
  VertexNameMap products_index = boost::get(boost::vertex_name, products_graph);

  // also compute pat_prod_cplx_mapping
  pat_prod_cplx_mapping.clear();
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
    bool is_simple_mapping = ei == edge_end;

    const Node& pat_mol = patterns_index[map_it->second];
    assert(pat_mol.is_mol);

    // the molecules must have the same compartment,
    // the search above in find_best_product_to_pattern_mapping must ignore compartments
    if (pat_mol.mol->compartment_id != prod_mol.mol->compartment_id) {
      continue;
    }

    // the graph has just a pointer to the molecule instance and since the graph
    // should be independent on the reactions, we are not storing any indices there
    uint reac_cplx_index = find_mol_instance_with_address(reactants, pat_mol.mol);
    uint prod_cplx_index = find_mol_instance_with_address(products, prod_mol.mol);

    // we can finally define our mapping
    pat_prod_cplx_mapping.push_back(CplxIndexPair(reac_cplx_index, prod_cplx_index, is_simple_mapping));
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
    const ElemMol& first_mi,
    const ElemMol& second_mi,
    const char* msg,
    std::ostream& out
) {
  // TODO: use computed component ordering from products_to_patterns_mapping
  // maybe store vertex descriptor along with molecule and component instances
  bool ok = true;
  for (size_t i = 0; i < first_mi.components.size(); i++) {
    const Component& first_compi = first_mi.components[i];

    if (i >= second_mi.components.size() ||
        second_mi.components[i].component_type_id != first_compi.component_type_id
    ) {
      out <<
          "Molecule " << msg << ": Component(s) " <<
          bng_data->get_component_type(first_mi.components[i].component_type_id).name <<
          " missing from molecule " <<
          second_mi.to_str(*bng_data) << ".";
      ok = false;
    }
  }

  return ok;
}


bool RxnRule::check_components_states(
    const ElemMol& prod_mi,
    const ElemMol& pat_mi,
    std::ostream& out
) {
  // TODO: use computed component ordering from products_to_patterns_mapping
  bool ok = true;
  for (size_t i = 0; i < prod_mi.components.size(); i++) {
    const Component& prod_compi = prod_mi.components[i];

    if (i >= pat_mi.components.size() ||
        pat_mi.components[i].component_type_id == prod_compi.component_type_id
    ) {
      const Component& pat_compi = pat_mi.components[i];

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

// called from semantic analyzer when rule is analyzed
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
    const ElemMol& prod_mi = *prod_node.mol;

    const Node& pat_node = patterns_index[map_it.second];
    assert(prod_node.is_mol);
    const ElemMol& pat_mi = *pat_node.mol;

    // allow special case for A_vol -> A_surf
    // TODO: maybe will need to be more strict
    if (pat_mi.elem_mol_type_id != prod_mi.elem_mol_type_id &&
        are_replaceable_elem_mols_do_not_check_connections(*bng_data, pat_mi, prod_mi)) {
      // ok
      continue;
    }
    else {
      // check that this molecule uses the same components in both directions
      // the components are ordered according to the definition in MolType
      ok = ok && check_components_mapping(pat_mi, prod_mi, "created in reaction rule", out);
      ok = ok && check_components_mapping(prod_mi, pat_mi, "used as pattern in reaction rule", out);

      ok = ok && check_components_states(prod_mi, pat_mi, out);
    }
  }

  return ok;
}


void RxnRule::move_products_that_are_also_reactants_to_be_the_first_products() {
#if defined(MCELL4_SORT_RXN_PRODUCTS_BY_NAME) || defined(MCELL4_SORT_RXN_PRODUCTS_BY_NAME_REV) || defined(MCELL4_SORT_RXN_PRODUCTS_BY_LENGTH_DESC)
  // this is a behavior of NFsim that products seem to be sorted by name,
  // let's assume by the first molecule of each complex
  // for now let's sort max 2 products
  if (products.size() == 2) {

    const string& name0 = products[0].to_str(false); // we do not care whether this is a surf or vol rxn
    const string& name1 = products[1].to_str(false);

#if defined(MCELL4_SORT_RXN_PRODUCTS_BY_NAME)
    if (name0 > name1) {
#elif defined(MCELL4_SORT_RXN_PRODUCTS_BY_NAME_REV)
    if (name0 < name1) {
#elif defined(MCELL4_SORT_RXN_PRODUCTS_BY_LENGTH_DESC)
    if (name0.size() < name1.size()) {
#endif
      Cplx tmp = products[0];
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
    // doing this only for simple complexes to keep compatibility with MCell3
    bool found = get_assigned_cplx_reactant_for_product(pi, true, ri);

    if (found) {
      // move product to the front
      Cplx prod = products[pi];
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

void RxnRule::get_reactant_indices(
    const species_id_t species_id, const SpeciesContainer& all_species,
    std::vector<uint>& indices
) {
  // update cache
  if (!species_can_be_reactant(species_id, all_species)) {
    return;
  }

  for (uint i = 0; i < reactants.size(); i++) {
    if (species_applicable_as_reactant[i].count(species_id)) {
      indices.push_back(i);
    }
  }
}

// stores indices where reactant can be applied into array indices (0 or 1)
void RxnRule::get_reactant_indices_uncached(
    const Cplx& cplx, const SpeciesContainer& all_species,
    std::vector<uint>& indices) const {

  for (size_t i = 0; i < reactants.size(); i++) {
    const Cplx& reactant = reactants[i];
    if (cplx.matches_pattern(reactant, true)) { // reactant is the pattern to be matched
      indices.push_back(i);
    }
  }
}


// returns the possible mappings using which patterns can be matched onto a reactant
// TODO: replace with get_reactant_indices
void RxnRule::get_bimol_reactant_indices(
    const species_id_t id1, const species_id_t id2,
    const SpeciesContainer& all_species,
    std::vector<std::pair<uint, uint>>& reac_indices) {

  assert(is_bimol());

  // update cache to be sure
  species_can_be_reactant(id1, all_species);
  species_can_be_reactant(id2, all_species);

  bool id1_matches[2] = {false, false};
  bool id2_matches[2] = {false, false};
  for (size_t i = 0; i < reactants.size(); i++) {
    id1_matches[i] = species_applicable_as_reactant[i].count(id1) != 0;
    id2_matches[i] = species_applicable_as_reactant[i].count(id2) != 0;
  }

  if (id1_matches[0] && id2_matches[1]) {
    // order fits
    reac_indices.push_back(make_pair(0u, 1u));
  }

  if (id1_matches[1] && id2_matches[0]) {
    // order is switched (may be true even if order fits)
    reac_indices.push_back(make_pair(1u, 0u));
  }
}


bool RxnRule::species_can_be_reactant(const species_id_t id, const SpeciesContainer& all_species) {

  // check caches first - if we processed this species
  if (species_applicable_as_any_reactant.count(id) != 0) {
    return true;
  }
  if (species_not_applicable_as_any_reactant.count(id) != 0) {
    return false;
  }

  // need to check whether the complex instance can be a reactant
  const Cplx& inst = all_species.get_as_cplx(id);

  std::vector<uint> indices;
  get_reactant_indices_uncached(inst, all_species, indices);

  if (!indices.empty()) {
    species_applicable_as_any_reactant.insert_unique(id);
    for (uint i: indices) {
      species_applicable_as_reactant[i].insert_unique(id);
    }

    return true;
  }
  else {
    species_not_applicable_as_any_reactant.insert_unique(id);
    species_not_applicable_as_reactant[0].insert_unique(id);
    if (is_bimol()) {
      species_not_applicable_as_reactant[1].insert_unique(id);
    }

    return false;
  }
}


bool RxnRule::species_can_be_bimol_reactants(
    const species_id_t id1, const species_id_t id2, const SpeciesContainer& all_species,
    uint* assigned_index1, uint* assigned_index2, bool* both_match_both_patterns
) {
  assert(is_bimol());

  // check whether either of the species can be reactant, also initializes cache
  if (!species_can_be_reactant(id1, all_species)) {
    return false;
  }
  if (!species_can_be_reactant(id2, all_species)) {
    return false;
  }

  // which reactants match id1?
  bool id1_matches[2] = {false, false};
  bool id2_matches[2] = {false, false};
  for (size_t i = 0; i < reactants.size(); i++) {
    id1_matches[i] = species_applicable_as_reactant[i].count(id1) != 0;
    id2_matches[i] = species_applicable_as_reactant[i].count(id2) != 0;
  }

  // there must be direct or crossed covering of the reactants, i.e.
  // id1 -> patA && id2 -> patB or
  // id1 -> patB && id2 -> patA
  bool res = (id1_matches[0] && id2_matches[1]) || (id1_matches[1] && id2_matches[0]);

  if (res && (assigned_index1 != nullptr || assigned_index2 != nullptr)) {
    uint index1;
    uint index2;
    if ((id1_matches[0] && id2_matches[1])) {
      index1 = 0;
      index2 = 1;
    }
    else {
      index1 = 1;
      index2 = 0;
    }
    if (assigned_index1 != nullptr) {
      *assigned_index1 = index1;
    }
    if (assigned_index2 != nullptr) {
      *assigned_index2 = index2;
    }
    if (both_match_both_patterns != nullptr) {
      *both_match_both_patterns =
          id1_matches[0] && id1_matches[1] &&
          id2_matches[0] && id2_matches[1];
    }
  }

  return res;
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

  return
      species_applicable_as_reactant[0].count(id) != 0 &&
      species_applicable_as_reactant[1].count(id) != 0;
}


void RxnRule::remove_species_id_references(const species_id_t id) {
  // this rxn rule might not have seen this species
  species_applicable_as_any_reactant.erase(id);
  species_applicable_as_reactant[0].erase(id);
  species_applicable_as_reactant[1].erase(id);
  species_not_applicable_as_any_reactant.erase(id);
  species_not_applicable_as_reactant[0].erase(id);
  species_not_applicable_as_reactant[1].erase(id);
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


bool RxnRule::update_rxn_rate(const float_t new_rate) {

  // skip if rate is the same
  if (base_rate_constant == new_rate) {
    return false;
  }

  // update the rate
  base_rate_constant = new_rate;

  // notify parents that update is needed
  for (RxnClass* user: rxn_classes_where_used) {
    // NOTE: this can be done more efficiently if needed,
    // no need to recompute products again
    user->init_rxn_pathways_and_rates(true);
  }

  return true;
}


std::string RxnRule::to_str(const bool with_rate_constant, const bool with_name, const bool with_id) const {
  stringstream ss;
  if (with_name && name != "") {
    ss << name << ": ";
  }

  ss << cplx_vector_to_str(reactants);
  ss << " -> ";
  if (!products.empty()) {
    ss << cplx_vector_to_str(products);
  }
  else {
    ss << COMPLEX_Null;
  }

  if (with_rate_constant) {
    ss << " " << base_rate_constant;
  }

  if (with_id) {
    ss << " (";
    ss << "id: " << id;
    ss << ")";
  }

  return ss.str();
}


std::string RxnRule::cplx_vector_to_str(const CplxVector& complexes) const {
  stringstream ss;
  for (size_t i = 0; i < complexes.size(); i++) {
    ss << complexes[i].to_str(!is_vol_rxn());

    if (i != complexes.size() - 1) {
      ss << " + ";
    }
  }
  return ss.str();
}


std::string RxnRule::reactants_to_str() const {
  return cplx_vector_to_str(reactants);
}


std::string RxnRule::products_to_str() const {
  return cplx_vector_to_str(products);
}


void RxnRule::dump_cplx_vector(
    const CplxVector& complexes,
    const std::string ind) const {

  for (size_t i = 0; i < complexes.size(); i++) {
    cout << ind << "Cplx " << i << ":\n";
    complexes[i].dump(true, ind + "  ");

    if (!is_simple()) {
      cout << ind + "  " << "  graph:\n";
      dump_graph(complexes[i].get_graph(), bng_data, ind + "  ");
    }
    cout << "\n";
  }
}


void RxnRule::dump(
    const bool for_diff, const std::string ind, std::ostream& out_reaction_rules) const {
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
      size_t max = base_variable_rates.size();
      #ifdef NDEBUG
        if (max > 10) {
          cout << ind << "Printing only the first 10 out of " << max << "\n";
          max = 10;
        }
      #endif
      for (size_t i = 0; i < max; i++) {
        cout << ind + "  " << "t: " << base_variable_rates[i].time << ", r: " << base_variable_rates[i].rate_constant << "\n";
      }
    }

    cout << ind << "mol_instances_are_fully_maintained: " << mol_instances_are_fully_maintained << "\n";

    cout << ind << "reactants:\n";
    dump_cplx_vector(reactants, ind);
    cout << ind << "products:\n";
    dump_cplx_vector(products, ind);
  }
}


} /* namespace BNG */
