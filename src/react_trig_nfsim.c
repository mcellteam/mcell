#include "config.h"

//#include "hashmap.h"
#include "map_c.h"
#include "logging.h"
#include "mcell_structs.h"
#include "nfsim_func.h"
#include "react.h"
#include "react_nfsim.h"
#include "react_util.h"
#include "sym_table.h"
#include <stdlib.h>
#include <string.h>
//#include "lru.h"

map_t reaction_map = NULL;
map_t reaction_preliminary_map = NULL;


void clear_maps() {
  hashmap_clear(reaction_map);
  hashmap_clear(reaction_preliminary_map);
}

// struct rxn *rx;

unsigned long lhash(const char *keystring) {
  unsigned long key = crc32((unsigned char *)(keystring), strlen(keystring));

  /* Robert Jenkins' 32 bit Mix Function */
  key += (key << 12);
  key ^= (key >> 22);
  key += (key << 4);
  key ^= (key >> 9);
  key += (key << 10);
  key ^= (key >> 2);
  key += (key << 7);
  key ^= (key >> 12);

  /* Knuth's Multiplicative Method */
  key = (key >> 3) * 2654435761;
  return key;
}

queryOptions initializeNFSimQueryForBimolecularReactions(
    struct graph_data *graph1, struct graph_data *graph2, const char *onlyActive) {
  // constant settings

  static const char *optionKeys[2] = {"numReactants", "onlyActive"};
  static char *optionValues[2];
  optionValues[0] = (char *)"2";
  optionValues[1] = (char *)onlyActive;

  static const int optionSeeds[2] = {1, 1};
  static char *speciesArray[2];
  // initialize speciesArray with the string we are going to query
  speciesArray[0] = graph1->graph_pattern;

  if (graph2)
    speciesArray[1] = graph2->graph_pattern;

  // copy these settings to the options object
  queryOptions options;
  options.initKeys = speciesArray;
  options.initValues = optionSeeds;
  // we can potentially query just for the bimolecular reactions a single
  // reactant is involved in
  // the only active flag would need to be off though
  if (graph2)
    options.numOfInitElements = 2;
  else
    options.numOfInitElements = 1;

  options.optionKeys = optionKeys;
  options.optionValues = optionValues;
  options.numOfOptions = 2;
  return options;
}

/**********************************************************************
 *
 * This function creates a queryOptions object for designing an NFSim
 * experiment query
 *
 * In: The abstract molecule whose nfsim status we are going to verify
 *
 * Out: the queryOptions object we will use to interact with nfsim
 *
 **********************************************************************/
queryOptions
initializeNFSimQueryForUnimolecularReactions(struct abstract_molecule *am) {
  // constant settings
  static const char *optionKeys[1] = {"numReactants"};
  static char *optionValues[1] = {(char *)"1"};
  static const int optionSeeds[1] = {1};
  static char *speciesArray[1];
  // initialize speciesArray with the string we are going to query
  // const char** speciesArray = CHECKED_MALLOC_ARRAY(char*, 1, "string array of
  // patterns");
  speciesArray[0] = am->graph_data->graph_pattern;
  // copy these settings to the options object
  queryOptions options;
  options.initKeys = speciesArray;
  options.initValues = optionSeeds;
  options.numOfInitElements = 1;

  options.optionKeys = optionKeys;
  options.optionValues = optionValues;
  options.numOfOptions = 1;
  // free(speciesArray);
  return options;
}

/*************************************************************************
   In: hashA - hash value for first molecule
       hashB - hash value for second molecule
       reacA - species of first molecule
       reacB - species of second molecule
   Out: 1 if any reaction exists naming the two specified reactants, 0
       otherwise.
   Note: This is a quick test used to determine which per-species lists to
   traverse when checking for mol-mol collisions.
*************************************************************************/

long rxnfound = 0;
long rxnmissed = 0;


int trigger_bimolecular_preliminary_nfsim(struct abstract_molecule *reacA,
                                          struct abstract_molecule *reacB) {

  if (reaction_preliminary_map == NULL)
    reaction_preliminary_map = hashmap_new();

  bool *isValidReaction = NULL;
  unsigned long reactionHash = reacA->graph_data->graph_pattern_hash +
                               reacB->graph_data->graph_pattern_hash;
  int error = hashmap_get_nohash(reaction_preliminary_map, reactionHash,
                                 reactionHash, (void **)(&isValidReaction));
  // error = find_in_cache(reaction_key, rx);

  // XXX: it might be worth it to return the rx object since we already queried
  // it
  if (error == MAP_OK) {
    if (isValidReaction != NULL)
      return 1;

    rxnfound++;

    return 0;
  }

  rxnmissed++;

  queryOptions options = initializeNFSimQueryForBimolecularReactions(
      reacA->graph_data, reacB->graph_data, "1");

  void *results = mapvectormap_create();
  initAndQueryByNumReactant_c(options, results);
  // if we have reactions great! although we can't decide what to do with them
  // yet. bummer.

  if (mapvectormap_size(results) > 0) {
    isValidReaction = CHECKED_MALLOC_ARRAY(bool, 1, "this is a valid reaction");
    *isValidReaction = true;
    mapvectormap_delete(results);

    error = hashmap_put_nohash(reaction_preliminary_map, reactionHash,
                               reactionHash, isValidReaction);
    return 1;
  } else {
    // if we know there's no reactions there's no need to check again later
    mapvectormap_delete(results);
    error = hashmap_put_nohash(reaction_preliminary_map, reactionHash,
                               reactionHash, NULL);
    return 0;
  }
}

/***********

********/
int trigger_bimolecular_nfsim(struct volume *state,
                              struct abstract_molecule *reacA,
                              struct abstract_molecule *reacB, short orientA,
                              short orientB, struct rxn **matching_rxns) {

  int error;
  int num_matching_rxns = 0;

  if (reaction_map == NULL)
    reaction_map = hashmap_new();

  // memset(&reaction_key[0], 0, sizeof(reaction_key));
  struct rxn *rx = NULL;

  unsigned long reactionHash = reacA->graph_data->graph_pattern_hash +
                               reacB->graph_data->graph_pattern_hash;
  // sprintf(reaction_key,"%lu",reacA->graph_pattern_hash +
  // reacB->graph_pattern_hash);
  // mcell_log("reaction_key %s %s %s",reacA->graph_pattern,
  // reacB->graph_pattern, reaction_key);
  // A + B ->C is the same as B +A -> C
  // if (strcmp(reacA->graph_pattern, reacB->graph_pattern) < 0)
  //    sprintf(reaction_key,"%s-%s",reacA->graph_pattern,reacB->graph_pattern);
  // else
  //    sprintf(reaction_key,"%s-%s",reacB->graph_pattern,reacA->graph_pattern);

  error = hashmap_get_nohash(reaction_map, reactionHash, reactionHash,
                             (void **)(&rx));
  // error = find_in_cache(reaction_key, rx);

  if (error == MAP_OK) {
    // if(error != -1){

    if (rx != NULL) {
      int result = process_bimolecular(reacA, reacB, rx, orientA, orientB,
                                       matching_rxns, num_matching_rxns);

      if (result == 1)
        num_matching_rxns++;
    }
    return num_matching_rxns;
  }

  // mcell_log("+++++ %s %s %s",reacA->graph_pattern, reacB->graph_pattern,
  // reaction_key);
  queryOptions options = initializeNFSimQueryForBimolecularReactions(
      reacA->graph_data, reacB->graph_data, "1");
  // reset, init, query the nfsim system
  void *results = mapvectormap_create();
  initAndQueryByNumReactant_c(options, results);

  // XXX: it would probably would be more natural to make a method that queries
  // all the reactions
  // associated with a given species
  if (mapvectormap_size(results) > 0) {
    rx = new_reaction();
    initializeNFSimReaction(state, rx, 2, results, reacA, reacB);

    /*int result = process_bimolecular(reacA, reacB, rx, orientA, orientB,*/
    /*                matching_rxns, num_matching_rxns, 0);*/
    int result = process_bimolecular(reacA, reacB, rx, orientA, orientB,
                                     matching_rxns, num_matching_rxns);

    if (result == 1)
      num_matching_rxns++;
  }
  // store value in hashmap

  error = hashmap_put_nohash(reaction_map, reactionHash, reactionHash, rx);
  // add_to_cache(reaction_key, rx);

  // CLEANUP
  // delete_reactantQueryResults(query2);
  mapvectormap_delete(results);
  return num_matching_rxns;
}

/*


is_surface: true if there is a surface reactant
*/
int adjust_rates_nfsim(struct volume *state, struct rxn *rx, bool is_surface) {
  double pb_factor = 0.0;
  // int max_num_surf_products = set_product_geometries(path, rx, prod);
  pb_factor = compute_pb_factor(
      state->time_unit, state->length_unit, state->grid_density,
      state->rx_radius_3d /
          state->r_length_unit, // transform back to a unitless scale
      &state->rxn_flags,
      &state->create_shared_walls_info_flag, rx,
      is_surface); // max_num_surf_products);
  rx->pb_factor = pb_factor;
  // mcell_log("!!!pb_factor %.10e",pb_factor);

  // JJT: balance out rate (code from scale_rxn_probabilities)
  // extracted because we only want to change the rate for one path

  double rate;
  for (int i = 0; i < rx->n_pathways; i++) {

    // if (!rx->rates || !rx->rates[i]) {
    rate = pb_factor * rx->cum_probs[i];
    mcell_log("!!%s %.10e %.10e", rx->external_reaction_data[i].reaction_name,
              rx->cum_probs[i], rate);
    //} else
    //  rate = 0.0;
    rx->cum_probs[i] = rate;
  }

  return 0;

  /*if (scale_rxn_probabilities(&state->reaction_prob_limit_flag, state->notify,
                          path, rx, pb_factor))
    return 1;*/
}

int initializeNFSimReaction(struct volume *state, struct rxn *r,
                            int n_reactants, void *results,
                            struct abstract_molecule *reacA,
                            struct abstract_molecule *reacB) {

  char **resultKeys = mapvectormap_getKeys(results);
  void *headComplex = mapvectormap_get(results, resultKeys[0]);

  // we know that it only contains one result
  int keysize = mapvectormap_size(results);

  int headNumAssociatedReactions = mapvector_size(headComplex);
  state->n_NFSimPReactions += headNumAssociatedReactions;

  r->cum_probs = CHECKED_MALLOC_ARRAY(double, headNumAssociatedReactions,
                                      "cumulative probabilities");
  r->external_reaction_data = CHECKED_MALLOC_ARRAY(
      struct external_reaction_datastruct, headNumAssociatedReactions,
      "external reaction names");
  r->n_pathways = headNumAssociatedReactions;
  r->product_idx = CHECKED_MALLOC_ARRAY(u_int, headNumAssociatedReactions + 1,
                                        "the different possible products");
  r->product_idx_aux = CHECKED_MALLOC_ARRAY(int, headNumAssociatedReactions + 1,
                                            "the different possible products");
  // r->product_graph_pattern =
  // CHECKED_MALLOC_ARRAY(char**,query2.numOfAssociatedReactions[0],
  //                                  "graph patterns of the possible
  //                                  products");

  r->reactant_graph_data =
      CHECKED_MALLOC_ARRAY(struct graph_data *, n_reactants,
                           "graph patterns of the possible products");

  r->reactant_graph_data[0] = reacA->graph_data;
  if (reacB) {
    r->reactant_graph_data[1] = reacB->graph_data;
  }

  r->product_graph_data =
      CHECKED_MALLOC_ARRAY(struct graph_data **, headNumAssociatedReactions,
                           "graph patterns of the possible products");

  r->nfsim_players =
      CHECKED_MALLOC_ARRAY(struct species **, headNumAssociatedReactions,
                           "products associated to each path");

  r->nfsim_geometries =
      CHECKED_MALLOC_ARRAY(short *, headNumAssociatedReactions,
                           "geometries associated to each path");

  memset(r->nfsim_players, 0,
         sizeof(struct species **) * headNumAssociatedReactions);

  r->n_reactants = n_reactants;

  // XXX:do we really have to go over all of them or do we need to filter out
  // repeats?
  // for (int i=0; i<query2.numOfResults; i++){
  void *pathInformation;
  for (int path = 0; path < headNumAssociatedReactions; path++) {
    pathInformation = mapvector_get(headComplex, path);
    r->cum_probs[path] = atof(map_get(pathInformation, "rate"));
    r->external_reaction_data[path].reaction_name =
        strdup(map_get(pathInformation, "name"));
    r->external_reaction_data[path].resample = 0;
    if (strcmp(map_get(pathInformation, "resample"), "true") == 0)
      r->external_reaction_data[path].resample = 1;
    r->product_idx[path] = 0;
    r->product_idx_aux[path] = -1;
  }

  // XXX: is this necessary?
  // temporarily initialize reaction list with nfsim_players
  r->players = CHECKED_MALLOC_ARRAY(struct species *, n_reactants,
                                    "reaction players array");
  r->geometries =
      CHECKED_MALLOC_ARRAY(short, n_reactants, "geometry players array");

  r->players[0] = reacA->properties;
  r->geometries[0] = 0;
  if (reacB != NULL) {
    r->players[1] = reacB->properties;
    r->geometries[1] = 0;
  }

  // nfsim diffusion function, depending on whether the user wants us to do
  // this.
  initialize_rxn_diffusion_functions(r);

  bool orientation_flag1 = 0, orientation_flag2 = 0;
  int reactantOrientation1, reactantOrientation2;

  calculate_reactant_orientation(reacA, reacB, &orientation_flag1,
                                 &orientation_flag2, &reactantOrientation1,
                                 &reactantOrientation2);
  if (orientation_flag1) {
    r->geometries[0] = reactantOrientation1;
  }
  if (reacB && orientation_flag2) {
    r->geometries[1] = reactantOrientation2;
  }

  // adjust reaction probabilities
  if (reacB != NULL)
    mcell_log("++++ %s %s", reacA->graph_data->graph_pattern,
              reacB->graph_data->graph_pattern);
  else
    mcell_log("---- %s ", reacA->graph_data->graph_pattern);
  adjust_rates_nfsim(state, r, orientation_flag1 & orientation_flag2);

  // calculate cummulative probabilities
  for (int n_pathway = 1; n_pathway < r->n_pathways; ++n_pathway) {
    r->cum_probs[n_pathway] += r->cum_probs[n_pathway - 1];
  }

  // assert(r->cum_probs[r->n_pathways-1] <= 1.0);

  if (r->n_pathways > 0)
    r->min_noreaction_p = r->max_fixed_p = r->cum_probs[r->n_pathways - 1];
  else
    r->min_noreaction_p = r->max_fixed_p = 1.0;

  // cleanup
  for (int i = 0; i < keysize; i++)
    free(resultKeys[i]);

  free(resultKeys);

  return 0;
}

struct rxn *pick_unimolecular_reaction_nfsim(struct volume *state,
                                             struct abstract_molecule *am) {

  int error;
  struct rxn *rx = NULL;
  if (reaction_map == NULL)
    reaction_map = hashmap_new();

  // memset(&reaction_key[0], 0, sizeof(reaction_key));
  // sprintf(reaction_key,"%s",am->graph_pattern);
  // sprintf(reaction_key,"%lu",am->graph_pattern_hash);

  // check in the hashmap in case this is a reaction we have encountered before
  error =
      hashmap_get_nohash(reaction_map, am->graph_data->graph_pattern_hash,
                         am->graph_data->graph_pattern_hash, (void **)(&rx));
  // error = find_in_cache(reaction_key, rx);

  if (error == MAP_OK) {
    return rx;
  }

  // if(error != -1)
  //    return rx;

  // otherwise build the object
  queryOptions options = initializeNFSimQueryForUnimolecularReactions(am);
  // reset, init, query the nfsim system
  void *results = mapvectormap_create();
  initAndQueryByNumReactant_c(options, results);

  // XXX: it would probably would be more natural to  make a method that queries
  // all the reactions
  // associated with a given species
  if (mapvectormap_size(results) > 0) {
    rx = new_reaction();
    initializeNFSimReaction(state, rx, 1, results, am, NULL);
  }

  // store newly created reaction in the hashmap
  error = hashmap_put_nohash(reaction_map, am->graph_data->graph_pattern_hash,
                             am->graph_data->graph_pattern_hash, rx);
  // add_to_cache(reaction_key, rx);

  // CLEANUP
  mapvectormap_delete(results);
  return rx;
  // delete_reactantQueryResults(query2);
}
