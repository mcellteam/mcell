#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "count_util.h"
#include "grid_util.h"
#include "logging.h"
#include "mcell_reactions.h"
#include "nfsim_func.h"
#include "react.h"
#include "react_nfsim.h"
#include "react_util.h"
#include "rng.h"
#include "util.h"
#include "vol_util.h"
#include "wall_util.h"

static queryOptions
initializeNFSimQueryforUnimolecularFiring(struct abstract_molecule *am,
                                          const char *external_path);

static bool calculate_nfsim_reactivity(struct graph_data *);
static void calculate_nfsim_diffusion_derived_data(struct volume *state,
                                                   struct graph_data *data);

void set_nfsim_product_geometries(struct pathway *path, struct rxn *rx,
                                  int prod_orientation, int prod_index) {
  if ((prod_orientation + path->orientation1) *
              (prod_orientation - path->orientation1) ==
          0 &&
      prod_orientation * path->orientation1 != 0) {
    if (prod_orientation == path->orientation1)
      rx->geometries[prod_index] = 1;
    else
      rx->geometries[prod_index] = -1;
  } else if (rx->n_reactants > 1 &&
             (prod_orientation + path->orientation2) *
                     (prod_orientation - path->orientation2) ==
                 0 &&
             prod_orientation * path->orientation2 != 0) {
    if (prod_orientation == path->orientation2)
      rx->geometries[prod_index] = 2;
    else
      rx->geometries[prod_index] = -2;
  } else {
    rx->geometries[prod_index] = 1;
  }
}

queryOptions initializeNFSimQueryNoFiring(struct abstract_molecule *am) {
  // constant settings
  queryOptions options;

  static const char *optionKeys[1] = {"systemQuery"};

  options.optionValues = CHECKED_MALLOC_ARRAY(char *, 1, "option array");

  options.optionValues[0] = strdup("complex");

  // initialize speciesArray with the string we are going to query
  char **speciesArray =
      CHECKED_MALLOC_ARRAY(char *, 1, "string array of patterns");
  speciesArray[0] = am->graph_data->graph_pattern;

  static const int optionSeeds[1] = {1};

  // copy these settings to the options object
  options.initKeys = speciesArray;
  options.initValues = optionSeeds;
  options.numOfInitElements = 1;

  // experiment design: query for reactions with 1 reactant
  options.optionKeys = optionKeys;
  // options.optionValues = optionValues;
  options.numOfOptions = 1;
  return options;
}

/**********************************************************************
 *
 * This function creates a queryOptions object for designing an NFSim
 * experiment query after a unimolecular reaction fires.
 *
 * In: The abstract molecule whose nfsim status we are going to verify
 *
 * Out: the queryOptions object we will use to interact with nfsim
 *
 **********************************************************************/
queryOptions
initializeNFSimQueryforUnimolecularFiring(struct abstract_molecule *am,
                                          const char *external_path) {
  // constant settings
  queryOptions options;

  static const char *optionKeys[2] = {"systemQuery", "reaction"};

  options.optionValues = CHECKED_MALLOC_ARRAY(char *, 2, "option array");

  options.optionValues[0] = strdup("complex");
  options.optionValues[1] = strdup(external_path);

  // initialize speciesArray with the string we are going to query
  char **speciesArray =
      CHECKED_MALLOC_ARRAY(char *, 1, "string array of patterns");
  speciesArray[0] = am->graph_data->graph_pattern;

  static const int optionSeeds[1] = {1};

  // copy these settings to the options object
  options.initKeys = speciesArray;
  options.initValues = optionSeeds;
  options.numOfInitElements = 1;

  // experiment design: query for reactions with 1 reactant
  options.optionKeys = optionKeys;
  // options.optionValues = optionValues;
  options.numOfOptions = 2;
  return options;
}

queryOptions
initializeNFSimQueryforBimolecularFiring(struct abstract_molecule *am,
                                         struct abstract_molecule *am2,
                                         const char *external_path) {
  // constant settings
  queryOptions options;

  static const char *optionKeys[2] = {"systemQuery", "reaction"};

  options.optionValues = CHECKED_MALLOC_ARRAY(char *, 2, "option array");

  // options.optionValues[0] = CHECKED_MALLOC_ARRAY(char, strlen("complex"),
  // "reaction which we will step over");
  // strcpy(options.optionValues[0], "complex");

  // optionValues[0] = "complex";

  // options.optionValues[1] = CHECKED_MALLOC_ARRAY(char, strlen(external_path),
  // "reaction which we will step over");
  // strcpy(options.optionValues[1], external_path);

  options.optionValues[0] = strdup("complex");
  options.optionValues[1] = strdup(external_path);

  // initialize speciesArray with the string we are going to query
  char **speciesArray =
      CHECKED_MALLOC_ARRAY(char *, 2, "string array of patterns");
  speciesArray[0] = am->graph_data->graph_pattern;
  speciesArray[1] = am2->graph_data->graph_pattern;

  static const int optionSeeds[2] = {1, 1};

  // copy these settings to the options object
  options.initKeys = speciesArray;
  options.initValues = optionSeeds;
  options.numOfInitElements = 2;

  // experiment design: query for reactions with 1 reactant
  options.optionKeys = optionKeys;
  // options.optionValues = optionValues;
  options.numOfOptions = 2;
  return options;
}

static void find_objects(struct geom_object* current_parent,
                  const char* name1, struct geom_object** obj1,
                  const char* name2, struct geom_object** obj2) {

  struct geom_object* curr = current_parent;
  while (curr != NULL) {
    // did we find our object?
    if (*obj1 == NULL && curr->last_name != NULL && strcmp(curr->last_name, name1) == 0) {
      *obj1 = curr;
    }
    if (*obj2 == NULL && curr->last_name != NULL && strcmp(curr->last_name, name2) == 0) {
      *obj2 = curr;
    }
    if (*obj1 != NULL && *obj2 != NULL) {
      // terminate search
      return;
    }
    // check children
    find_objects(curr->first_child, name1, obj1, name2, obj2);

    // next object in list
    curr = curr->next;
  }
}

static bool has_region(struct geom_object* obj, const char* reg_name) {

  for (struct region_list *r = obj->regions; r != NULL; r = r->next) {
    if (r->reg != NULL && strcmp(r->reg->region_last_name, reg_name) == 0) {
      return true;
    }
  }
  return false;
}

int prepare_reaction_nfsim(struct volume *world, struct rxn *rx, void *results,
                           int path, struct abstract_molecule *reac,
                           struct abstract_molecule *reac2) {

  const char *product_pattern = NULL;
  void *individualResult;
  int numOfResults = mapvector_size(results);

  for (int productIdx = 0; productIdx < numOfResults; productIdx++) {
    individualResult = mapvector_get(results, productIdx);
    product_pattern = map_get(individualResult,
                              "label"); // results->results[productIdx].label;
    constructNauty_c(product_pattern, 1);
    // TODO: we are ignoring optimizing for overlaps for now
    /*if(strcmp(reac->graph_pattern, product_pattern) == 0 && overlapFlag == 0){
      overlapFlag = 1;
    }
    else{
        rx->product_idx_aux[path]++;

    }*/
  }
  rx->product_idx_aux[path] = numOfResults;
  if (rx->nfsim_players[path] == NULL) {
    rx->nfsim_players[path] = CHECKED_MALLOC_ARRAY(
        struct species *, numOfResults, "reaction players array");
    rx->nfsim_geometries[path] = CHECKED_MALLOC_ARRAY(
        short, numOfResults, "geometries associated to this path");
  }

  rx->product_graph_data[path] = CHECKED_MALLOC_ARRAY(
      struct graph_data *, rx->product_idx_aux[path],
      "graph patterns for products that have been added to the system");
  int counter = 0;
  const char *diffusion;
  rx->external_reaction_data[path].products = numOfResults;
  for (int productIdx = 0; productIdx < numOfResults; productIdx++) {
    individualResult = mapvector_get(results, productIdx);

    product_pattern = map_get(individualResult,
                              "label"); // results->results[productIdx].label;

    // query  graph_pattern hashmap instead of recreating stuff
    unsigned long graph_hash = lhash(product_pattern);
    int error =
        get_graph_data(graph_hash, &rx->product_graph_data[path][counter]);
    if (error != 0) {
      rx->product_graph_data[path][counter] = CHECKED_MALLOC_ARRAY(
          struct graph_data, 1, "graph pattern for a single path");
      rx->product_graph_data[path][counter]->graph_pattern =
          strdup(product_pattern);
      rx->product_graph_data[path][counter]->graph_pattern_hash =
          lhash(product_pattern);
      diffusion = map_get(individualResult, "diffusion_function");
      if (diffusion) {
        rx->product_graph_data[path][counter]->graph_diffusion =
            atof(diffusion);
        calculate_nfsim_diffusion_derived_data(
            world, rx->product_graph_data[path][counter]);
      } else {
        rx->product_graph_data[path][counter]->graph_diffusion = -1;
        rx->product_graph_data[path][counter]->time_step = -1;
        rx->product_graph_data[path][counter]->space_step = -1;
      }

      calculate_nfsim_reactivity(rx->product_graph_data[path][counter]);

      store_graph_data(graph_hash, rx->product_graph_data[path][counter]);
      world->n_NFSimSpecies += 1;
    }
    counter++;
    //}
  }

  // recalculate the position of all players
  int num_players = rx->n_reactants;
  int kk = rx->n_pathways;
  if (kk <= RX_SPECIAL)
    kk = 1;

  for (int n_pathway = 0; n_pathway < kk; n_pathway++) {
    int k = rx->product_idx_aux[n_pathway] + rx->n_reactants;
    rx->product_idx[n_pathway] = num_players;
    num_players += k;
  }
  rx->product_idx[kk] = num_players;

  // we will be recreating the players and geometries arrays. this might not be
  // the most efficient approach but this is because we don't know the total
  // number of products before each path is fired individually in nfsim
  // XXX: maybe fire them manually and just get the number of products per path
  // even if we don't store path information?
  if (rx->players != NULL)
    free(rx->players);
  free(rx->geometries);

  rx->players = CHECKED_MALLOC_ARRAY(struct species *, num_players,
                                     "reaction players array");
  rx->geometries =
      CHECKED_MALLOC_ARRAY(short, num_players, "reaction geometries array");

  memset(rx->players, 0, sizeof(struct species *) * num_players);
  memset(rx->geometries, 0, sizeof(short) * num_players);

  // create pathway that will form players
  struct pathway *pathp = (struct pathway *)CHECKED_MALLOC_STRUCT(
      struct pathway, "reaction pathway");
  if (pathp == NULL) {
    return -1;
  }
  memset(pathp, 0, sizeof(struct pathway));

  /* Scan reactants, copying into the new pathway */
  int num_vol_mols = 0;
  int num_surface_mols = 0;
  int all_3d = 1;
  /*int complex_type = 0;*/
  int reactant_idx = 0;
  int oriented_count = 0;
  /*int num_complex_reactants = 0;*/
  bool orientation_flag1 = 0, orientation_flag2 = 0;
  int reactantOrientation1, reactantOrientation2, productOrientation;

  // obtain orientation information
  calculate_reactant_orientation(reac, reac2, &orientation_flag1,
                                 &orientation_flag2, &reactantOrientation1,
                                 &reactantOrientation2);
  struct sym_entry *nfsim_molecule = reac->properties->sym;

  // create first reactant entry
  rx->geometries[0] = reactantOrientation1;
  struct mcell_species *reactants = mcell_add_to_species_list(
      nfsim_molecule, orientation_flag1, reactantOrientation1, NULL);

  // second reactant entry
  if (reac2 != NULL) {
    rx->geometries[1] = reactantOrientation2;
    nfsim_molecule = reac2->properties->sym;
    reactants = mcell_add_to_species_list(nfsim_molecule, orientation_flag2,
                                          reactantOrientation2, reactants);
  }
  // we will be storing prodcut compartment information in here
  compartmentStruct *compartmentInfoArray =
      CHECKED_MALLOC_ARRAY(compartmentStruct, numOfResults,
                           "Creating array of compartment information");

  for (int i = 0; i < numOfResults; i++) {
    // what compartment is the species now in
    individualResult = mapvector_get(results, i);
    compartmentInfoArray[i] = getCompartmentInformation_c(map_get(
        individualResult,
        "compartment")); // getCompartmentInformation_c(results->results[i].compartment);
  }

  struct species *nfsim_molecule_template;

  // get the mcell species proxy we will be using
  if (compartmentInfoArray[0].spatialDimensions == 2)
    nfsim_molecule_template = world->global_nfsim_surface;
  else
    nfsim_molecule_template = world->global_nfsim_volume;

  // calculate orientation information if its not a vol vol reaciton
  if (orientation_flag2) {
    orientation_flag1 = true;
    individualResult = mapvector_get(results, 0);
    const char* outside = compartmentInfoArray[0].outside;
    // originalCompartment is in fact the taget compartment where we should place out product
    const char* originalCompartment = map_get(individualResult, "originalCompartment");

    bool originalCompartmentEmpty = strcmp(originalCompartment, "") == 0;

    if (!originalCompartmentEmpty && strcmp(outside, originalCompartment) == 0) {
      // outside and original are the same, not completely sure what to do here,
      // keeping original implementation
      productOrientation = -1;
    }
    else if (!originalCompartmentEmpty && strcmp(outside, "") != 0) {
      // 1) find object with "outside or "originalCompartment" name in world->root_instance,
      // then check its regions
      // 2) a region of its object is either "ALL" or its surface
      // 3) if we go from surface into the object -> orientation == -1,
      //    all other cases are kept as they were implemented originally because we
      //    do not have other information about the hierarchy of objects
      struct geom_object* objOutside = NULL;
      struct geom_object* objOrigCompartment = NULL;
      find_objects(world->root_instance, outside, &objOutside, originalCompartment, &objOrigCompartment);

      // if we found only one of these objects, check its regions
      if (objOutside != NULL && objOrigCompartment != NULL) {
        // both are objects, default behavior should be ok
        mcell_warn(
            "Handling orientation of NFsim reaction where both compartents %s (orig) and %s (outside) are objects is not supported yet, "
            "using default orientation 1 (outside)", originalCompartment, outside);
        productOrientation = 1;
      } else if (objOutside == NULL && objOrigCompartment == NULL) {
        mcell_error(
            "NFsim reaction returned compartments %s (orig) and %s (outside) that were not identified as object, this is not supported yet.",
            originalCompartment, outside
        );
        exit(1);
      } else {
        if (objOrigCompartment != NULL) {
          if (has_region(objOrigCompartment, outside)) {
            // we go from outside (e.g. PM) inside (e.g. CP)
            productOrientation = -1;
          } else {
            // relationship is not known...
            mcell_warn(
                "Handling orientation of NFsim reaction where both compartent %s (orig) and %s (outside) are not in object-region relationship is not supported yet, "
                "using default orientation 1 (outside)", originalCompartment, outside);
            productOrientation = 1;
          }
        } else {
          if (has_region(objOutside, originalCompartment)) {
            // we go from inside (e.g. CP) outside (e.g. PM)
            productOrientation = 1;
          } else {
            // relationship is not known but we can assume that we go onto surface because the original compartment object is NULL,
            // therefore it is a region not an object,
            // selecting 1 because we need to stay on the 'outside' side of the surface
            productOrientation = 1;
          }
        }
      }
    } else {
      // outside is empty, we get this type of information from reactions like this:
      // MemA@PM -> Mem@PM+A@EC  -  we need to go outside from membrane PM to EC
      productOrientation = 1;
    }

  } else {
    orientation_flag1 = false;
    productOrientation = 0;
  }

  struct mcell_species *products =
      mcell_add_to_species_list(nfsim_molecule_template->sym, orientation_flag1,
                                productOrientation, NULL);
  rx->nfsim_geometries[path][0] = productOrientation;
  // if theres more than one product
  for (int i = 1; i < numOfResults; i++) {
    // compartmentInfo =
    // getCompartmentInformation_c(results->results[i].compartment);
    if (compartmentInfoArray[i].spatialDimensions == 2)
      nfsim_molecule_template = world->global_nfsim_surface;
    else
      nfsim_molecule_template = world->global_nfsim_volume;

    if (orientation_flag2) {

      individualResult = mapvector_get(results, i);
      if (strcmp(compartmentInfoArray[i].outside,
                 map_get(individualResult, "originalCompartment")) ==
          0) { // results->results[i].originalCompartment) == 0){
        productOrientation = -1;
      } else {
        productOrientation = 1;
      }

    } else {
      productOrientation = 0;
    }

    rx->nfsim_geometries[path][i] = productOrientation;
    products = mcell_add_to_species_list(nfsim_molecule_template->sym,
                                         orientation_flag1, productOrientation,
                                         products);
  }

  // free up compartment struct helpers
  for (int i = 0; i < numOfResults; i++) {
    delete_compartmentStructs(compartmentInfoArray[i]);
  }
  free(compartmentInfoArray);

  // create out pathway
  if (extract_reactants(pathp, reactants, &reactant_idx, &num_vol_mols,
                        &num_surface_mols, &all_3d, &oriented_count) != 0) {
    return -1;
  }

  int num_surf_products = 0;
  /*int num_complex_products = 0;*/
  int bidirectional = 0;

  if (extract_products(world->notify, pathp, products, &num_surf_products,
                       bidirectional, all_3d) == MCELL_FAIL) {
    return MCELL_FAIL;
  }

  mcell_delete_species_list(reactants);
  mcell_delete_species_list(products);

  int k = 0;
  struct product *prod = NULL;

  // store this nfsim path information
  counter = 0;

  // XXX: it might be necessary to clear out the stuff in pathp at the very end
  for (prod = pathp->product_head; counter < numOfResults;) {
    rx->nfsim_players[path][counter] = prod->prod;
    ++counter;
    if (counter < numOfResults)
      prod = prod->next;
  }
  prod->next = NULL;

  if (rx->players != NULL) {
    free(rx->players);
    free(rx->geometries);
  }

  rx->players = CHECKED_MALLOC_ARRAY(struct species *, num_players,
                                     "reaction players array");
  rx->geometries =
      CHECKED_MALLOC_ARRAY(short, num_players, "reaction geometries array");

  memset(rx->players, 0, sizeof(struct species *) * num_players);
  memset(rx->geometries, 0, sizeof(short) * num_players);

  // recreate players array
  rx->players[0] = pathp->reactant1;

  // if its a bimolecular reaction
  if (reac2 != NULL)
    rx->players[1] = pathp->reactant2;

  for (int n_pathway = 0; n_pathway < rx->n_pathways; n_pathway++) {
    k = rx->product_idx[n_pathway] + rx->n_reactants;
    counter = 0;
    for (counter = 0; counter < rx->product_idx_aux[n_pathway]; counter++) {
      // XXX: right now we are ignoring recycled species which is inefficient

      // if (recycled1 == 0 && prod->prod == pathp->reactant1) {
      //  recycled1 = 1;
      //  kk = rx->product_idx[path] + 0;
      //}
      // else {
      kk = k;
      k++;
      //}
      // kk = rx->product_idx[path] + 0;
      rx->players[kk] = rx->nfsim_players[n_pathway][counter];
      set_nfsim_product_geometries(
          pathp, rx, rx->nfsim_geometries[n_pathway][counter], kk);
    }
  } /* end for (n_pathway = 0, ...) */

  init_reaction_info(rx);
  rx->info[path].pathname = NULL;

  // adjust reaction probabilities
  // adjust_rates_nfsim(world, rx, pathp);

  // cleanup path information
  struct product *tmp = pathp->product_head;
  struct product *tmp2 = NULL;
  while (tmp != NULL) {
    tmp2 = tmp->next;
    free(tmp);
    tmp = tmp2;
  }
  free(pathp);

  return MCELL_SUCCESS;
}

void free_reaction_nfsim(struct rxn *rx, int path) {
  free(rx->nfsim_players[path]);
  free(rx->nfsim_geometries[path]);

  // i actually need to iterate over all products
  /*for(int i=0; i < rx->external_reaction_data[path].products;i++){
    free(rx->product_graph_data[path][i]->graph_pattern);
    free(rx->product_graph_data[path][i]);
  }*/
  free(rx->product_graph_data[path]);

  rx->nfsim_players[path] = NULL;
  rx->nfsim_geometries[path] = NULL;
  rx->product_graph_data[path] = NULL;
}

//int outcome_unimolecular_nfsim(struct volume *world, struct rxn *rx, int path,
//                               struct abstract_molecule *reac, double t) {
//  int result = RX_A_OK;
//
//  if (rx->product_idx_aux[path] == -1) {
//    mcell_log("uni restart %s\n", reac->graph_data->graph_pattern);
//    queryOptions options = initializeNFSimQueryforUnimolecularFiring(
//        reac, rx->external_reaction_data[path].reaction_name);
//
//    void *results = mapvector_create();
//
//    initAndQuerySystemStatus_c(options, results);
//
//    constructNauty_c(reac->graph_data->graph_pattern, -1);
//
//    // fill in the rxn react structure with the appropiate information
//    world->n_NFSimReactions += 1;
//    prepare_reaction_nfsim(world, rx, results, path, reac, NULL);
//  }
//  return result;
//}

/*
Calculate the space_Step and time_step associated to this molecule based on
right now it only considers global time_steps and space_steps
*/
void calculate_nfsim_diffusion_derived_data(struct volume *state,
                                            struct graph_data *data) {

  double global_time_unit = state->time_unit;

  const char *compartment1 =
      extractSpeciesCompartmentFromNauty_c(data->graph_pattern);
  compartmentStruct reactantCompartmentInfo1 =
      getCompartmentInformation_c(compartment1);
  // free string allocated in extractSpeciesCompartmentFromNauty_c
  free((char*)compartment1);

  if (!distinguishable(state->space_step, 0, EPS_C)) // Global timestep
  {
    data->space_step =
        sqrt(4.0 * 1.0e8 * data->graph_diffusion * global_time_unit) *
        state->r_length_unit;
    data->time_step = 1.0;
  } else /* Global spacestep */
  {
    double space_step = state->space_step * state->length_unit;
    if (reactantCompartmentInfo1.spatialDimensions == 2) {
      data->time_step =
          space_step * space_step /
          (MY_PI * 1.0e8 * data->graph_diffusion * global_time_unit);
    } else {
      data->time_step =
          space_step * space_step * MY_PI /
          (16.0 * 1.0e8 * data->graph_diffusion * global_time_unit);
    }
    data->space_step = sqrt(4.0 * 1.0e8 * data->graph_diffusion *
                            data->time_step * global_time_unit) *
                       state->r_length_unit;
  }

  // free strings allocated in getCompartmentInformation_c
  free(reactantCompartmentInfo1.name);
  free(reactantCompartmentInfo1.outside);
}

bool calculate_nfsim_reactivity(struct graph_data *graph) {
  graph->flags = 0;

  queryOptions options =
      initializeNFSimQueryForBimolecularReactions(graph, NULL, "0");
  void *results = mapvectormap_create();
  initAndQueryByNumReactant_c(options, results);

  bool dimensionalityFlag = true;
  if (mapvectormap_size(results) > 0) {
    char **resultKeys = mapvectormap_getKeys(results);
    // we know that it only contains one result
    void *headComplex = mapvectormap_get(results, resultKeys[0]);

    int resultSize = mapvectormap_size(results);
    for (int i = 0; i < resultSize; i++) {
      free(resultKeys[i]);
    }
    free(resultKeys);

    int headNumAssociatedReactions = mapvector_size(headComplex);
    void *pathInformation;

    for (int path = 0; path < headNumAssociatedReactions; path++) {
      pathInformation = mapvector_get(headComplex, path);
      const char *dimensionality =
          map_get(pathInformation, "reactionDimensionality");
      if (!dimensionality) {
        dimensionalityFlag = false;
        break;
      }
      if (strcmp(dimensionality, "VOLSURF") == 0) {
        graph->flags |= CAN_VOLSURF;
      }
      if (strcmp(dimensionality, "VOLVOL") == 0) {
        graph->flags |= CAN_VOLVOL;
      }
      if (strcmp(dimensionality, "SURFSURF") == 0) {
        graph->flags |= CAN_SURFSURF;
      }
    }
  }

  // cleanup
  mapvectormap_delete(results);

  if (dimensionalityFlag) {
    return true;
  } else {
    graph->flags = -1;
    return false;
  }
}

/**
Doesn't fire any reactions, it just queries NFSim for the properties of a given
graph pattern
used for initialization and copies them to the reac->graph_data object
**/
void properties_nfsim(struct volume *world, struct abstract_molecule *reac) {

  // XXX: this could be stored in a hashmap for initialization purposes if it
  // becomes oto much of an
  // issue to query everytime. we are not doi9ng it right now since this
  // function is only called
  // when initially releasing particles. moreover it might be worth it to have a
  // hashmap containing
  // graph_data properties instead of it being stored in every particle...

  // initialize system with only one molecule
  queryOptions options = initializeNFSimQueryNoFiring(reac);
  void *results = mapvector_create();
  initAndQuerySystemStatus_c(options, results);
  // get the first result since we are only querying information for one
  // molecule
  void *individualResult = mapvector_get(results, 0);
  const char *result = map_get(individualResult, "diffusion_function");

  if (result) {
    reac->graph_data->graph_diffusion = atof(result);
    calculate_nfsim_diffusion_derived_data(world, reac->graph_data);
  } else {
    reac->graph_data->graph_diffusion = -1;
    reac->graph_data->space_step = -1;
    reac->graph_data->time_step = -1;
    reac->get_diffusion = get_standard_diffusion;
    reac->get_space_step = get_standard_space_step;
    reac->get_time_step = get_standard_time_step;
  }
  mapvector_delete(results);

  // now lets get information about the reactionality of this reactant
  if (calculate_nfsim_reactivity(reac->graph_data)) {
    reac->get_flags = get_nfsim_flags;
  } else {
    reac->get_flags = get_standard_flags;
  }

  free(options.optionValues[0]);
  free(options.optionValues);
  free(options.initKeys);
}

int outcome_nfsim(struct volume *world, struct rxn *rx, int path,
                  struct abstract_molecule *reac,
                  struct abstract_molecule *reac2, double t) {
  int result = RX_A_OK;
  queryOptions options;
  // if this is a reaction rule that can be mapped in different ways we need to
  // resample
  if (rx->product_idx_aux[path] != -1 &&
      rx->external_reaction_data[path].resample == 1) {
    free_reaction_nfsim(rx, path);
    rx->product_idx_aux[path] = -1;
  }
  // if we don't have previous information about this path then build up the
  // rxn structure

  if (rx->product_idx_aux[path] == -1) {

    world->n_NFSimReactions += 1;

    if (reac2 == NULL)
      options = initializeNFSimQueryforUnimolecularFiring(
          reac, rx->external_reaction_data[path].reaction_name);
    else {
      options = initializeNFSimQueryforBimolecularFiring(
          reac, reac2, rx->external_reaction_data[path].reaction_name);
    }

    void *results = mapvector_create();

    initAndQuerySystemStatus_c(options, results);

    // queryResults results = initAndQuerySystemStatus_c(options);

    // frees up the option object
    free(options.optionValues[0]);
    free(options.optionValues[1]);
    free(options.optionValues);
    free(options.initKeys);

    // fill in the rxn react structure with the appropiate information
    prepare_reaction_nfsim(world, rx, results, path, reac, reac2);
    // frees up the query result
    mapvector_delete(results);

  }
  // otherwise just update populations
  else {
    for (int i = 0; i < rx->product_idx_aux[path]; i++) {
      constructNauty_c(rx->product_graph_data[path][i]->graph_pattern, 1);
    }
    // and clean the info information for good measure
    rx->info[path].pathname = NULL;
  }
  // also decrease reactant populations
  constructNauty_c(reac->graph_data->graph_pattern, -1);
  if (reac2 != NULL)
    constructNauty_c(reac2->graph_data->graph_pattern, -1);

  return result;
}
