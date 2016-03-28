#include "config.h"

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "logging.h"
#include "rng.h"
#include "util.h"
#include "grid_util.h"
#include "count_util.h"
#include "react.h"
#include "vol_util.h"
#include "macromolecule.h"
#include "wall_util.h"
#include "react_util.h"
#include "mcell_reactions.h"


 static queryOptions initializeNFSimQueryforUnimolecularFiring(struct abstract_molecule *am,
                                                           const char* external_path);


/**********************************************************************
 *
 * This function creates a queryOptions object for designing an NFSim experiment query after a 
 * unimolecular reaction fires
 * In: The abstract molecule whose nfsim status we are going to verify
 *
 * Out: the queryOptions object we will use to interact with nfsim
 *
 **********************************************************************/
 queryOptions initializeNFSimQueryforUnimolecularFiring(struct abstract_molecule *am,
                                                           const char* external_path){
    //constant settings
    queryOptions options;

    static const char* optionKeys[2]  = {"systemQuery","reaction"};

    options.optionValues = CHECKED_MALLOC_ARRAY(char*, 2, "option array");

    options.optionValues[0] = strdup("complex");
    options.optionValues[1] = strdup(external_path);


    //initialize speciesArray with the string we are going to query
    const char** speciesArray = CHECKED_MALLOC_ARRAY(char*, 1, "string array of patterns");
    speciesArray[0] = am->graph_pattern;

    static const int optionSeeds[1]= {1};


    //copy these settings to the options object
    options.initKeys = speciesArray;
    options.initValues = optionSeeds;
    options.numOfInitElements = 1;

    //experiment design: query for reactions with 1 reactant
    options.optionKeys = optionKeys;
    //options.optionValues = optionValues;
    options.numOfOptions = 2;
    return options;
}

 queryOptions initializeNFSimQueryforBimolecularFiring(struct abstract_molecule *am, struct abstract_molecule *am2,
                                                           const char* external_path){
    //constant settings
    queryOptions options;

    static const char* optionKeys[2]  = {"systemQuery","reaction"};

    options.optionValues = CHECKED_MALLOC_ARRAY(char*, 2, "option array");

    //options.optionValues[0] = CHECKED_MALLOC_ARRAY(char, strlen("complex"), "reaction which we will step over");
    //strcpy(options.optionValues[0], "complex");

    //optionValues[0] = "complex";

    //options.optionValues[1] = CHECKED_MALLOC_ARRAY(char, strlen(external_path), "reaction which we will step over");
    //strcpy(options.optionValues[1], external_path);

    options.optionValues[0] = strdup("complex");
    options.optionValues[1] = strdup(external_path);

    //initialize speciesArray with the string we are going to query
    const char** speciesArray = CHECKED_MALLOC_ARRAY(char*, 2, "string array of patterns");
    speciesArray[0] = am->graph_pattern;
    speciesArray[1] = am2->graph_pattern;

    static const int optionSeeds[2]= {1, 1};


    //copy these settings to the options object
    options.initKeys = speciesArray;
    options.initValues = optionSeeds;
    options.numOfInitElements = 2;

    //experiment design: query for reactions with 1 reactant
    options.optionKeys = optionKeys;
    //options.optionValues = optionValues;
    options.numOfOptions = 2;
    return options;
}


int prepare_reaction_nfsim(struct volume *world, struct rxn* rx, queryResults* results, 
  int path, struct abstract_molecule *reac, struct abstract_molecule *reac2)
{
  //mcell_log("+++++++++");
  //mcell_log("%s",reac->graph_pattern);
  //if(reac2 != NULL)
  //mcell_log("%s",reac2->graph_pattern);
  //mcell_log("--------");
  char* product_pattern = NULL;
  int overlapFlag = 0;


  for(int productIdx = 0; productIdx < results->numOfResults; productIdx++){
    product_pattern = results->results[productIdx].label;
    constructNauty_c(product_pattern, 1);
    //TODO: we are ignoring optimizing for overlaps for now
    /*if(strcmp(reac->graph_pattern, product_pattern) == 0 && overlapFlag == 0){
      overlapFlag = 1;
    }
    else{
        rx->product_idx_aux[path]++;
      
    }*/

  }
  rx->product_idx_aux[path] = results->numOfResults;
  if(rx->nfsim_players[path] == NULL){
    rx->nfsim_players[path] = CHECKED_MALLOC_ARRAY(struct species*, results->numOfResults,
                               "reaction players array");
  }


  rx->product_graph_pattern[path] = CHECKED_MALLOC_ARRAY(char*,rx->product_idx_aux[path],
                                      "graph patterns for products that have been added to the system");
  overlapFlag = 0;
  int counter = 0;

  for(int productIdx = 0; productIdx < results->numOfResults; productIdx++){
    product_pattern = results->results[productIdx].label;
    /*if(strcmp(reac->graph_pattern, product_pattern) == 0 && overlapFlag == 0){
      overlapFlag = 1;
    }*/
    //else{
      rx->product_graph_pattern[path][counter] = product_pattern;
      counter++;
    //}
  }


  //recalculate the position of all players
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



  //create pathway that will form players
  struct pathway *pathp = (struct pathway *)CHECKED_MALLOC_STRUCT(
      struct pathway, "reaction pathway");
  if (pathp == NULL) {
    return -1;
  }
  memset(pathp, 0, sizeof(pathp));

  /* Scan reactants, copying into the new pathway */
  int num_vol_mols = 0;
  int num_surface_mols = 0;
  int all_3d = 1;
  int complex_type = 0;
  int reactant_idx = 0;
  int oriented_count = 0;
  int num_complex_reactants = 0;
  bool orientation_flag;
  int orientation;
  //obtain the generic vol nfsim reactant
  struct sym_table* nfsim_molecule = reac->properties->sym;
  
  // if its a volume molecule
  if(reac->flags & ON_GRID == 0){
    orientation_flag = false;
    orientation = 0;

  }
  else{
    //else if its a surface molecule
    orientation_flag = true;
    orientation = 1;
  }

  struct mcell_species *reactants = mcell_add_to_species_list(nfsim_molecule, orientation_flag, orientation, 0, NULL);

  if(reac2 != NULL){
      nfsim_molecule = reac2->properties->sym;
      if(reac2->flags & ON_GRID == 0 && !orientation_flag){
        orientation_flag = false;
        orientation = 0;
      }
      else{
        //else if its a surface molecule
        orientation_flag = true;
        orientation = 1;
      }
      reactants = mcell_add_to_species_list(nfsim_molecule, orientation_flag, orientation, 0, reactants);    

  }

  compartmentStruct* compartmentInfoArray =CHECKED_MALLOC_ARRAY(compartmentStruct, results->numOfResults,
                               "Creating array of compartment information");


  for(int i =0; i < results->numOfResults; i++){
    // what compartment is the species now in
    compartmentInfoArray[i] = getCompartmentInformation_c(results->results[i].compartment);
  }


  struct species* nfsim_molecule_template;

  if(compartmentInfoArray[0].spatialDimensions == 2)
    nfsim_molecule_template = world->global_nfsim_surface;
  else
    nfsim_molecule_template = world->global_nfsim_volume;

  if(nfsim_molecule_template->flags & ON_GRID == 0  && !orientation_flag){
    orientation_flag = false;
    orientation = 0;
  }
  else{
    //else if its a surface molecule
    orientation_flag = true;
    //if we are moving towards the inside of a membrane
    if(strcmp(compartmentInfoArray[0].outside, 
              results->results[0].originalCompartment) == 0){
      orientation = -1;
    }
    else{
      orientation = 1;
    }
  }

  struct mcell_species *products =
        mcell_add_to_species_list(nfsim_molecule_template->sym, orientation_flag, orientation, 0, NULL);

  //FIXME: im leaking information here
  for(int i =1; i <results->numOfResults; i++){
    //compartmentInfo = getCompartmentInformation_c(results->results[i].compartment);
    if(compartmentInfoArray[i].spatialDimensions == 2)
      nfsim_molecule_template = world->global_nfsim_surface;
    else
      nfsim_molecule_template = world->global_nfsim_volume;

    if(nfsim_molecule_template->flags & ON_GRID == 0  && !orientation_flag){
      orientation_flag = false;
      orientation = 0;
    }
    else{
      //else if its a surface molecule
      orientation_flag = true;
      if(strcmp(compartmentInfoArray[i].outside, 
                results->results[i].originalCompartment) == 0){
        orientation = -1;
      }
      else{
        orientation = 1;
      }
    }

    products = mcell_add_to_species_list(nfsim_molecule_template->sym, orientation_flag, orientation, 0, products);
  }

  //free up compartment struct helpers
  for(int i =0; i < results->numOfResults; i++){
    delete_compartmentStructs(compartmentInfoArray[i]);
  }
  //create out pathway
  if (extract_reactants(pathp, reactants, &reactant_idx, &num_vol_mols,
                        &num_surface_mols, &num_complex_reactants, &all_3d,
                        &oriented_count, &complex_type) != 0) {
    return -1;
  }

  


  int num_surf_products = 0;
  int num_complex_products = 0;
  int bidirectional = 0;
  
  if (extract_products(world->notify, pathp, products, &num_surf_products,
                       &num_complex_products, bidirectional, complex_type,
                       all_3d) == MCELL_FAIL) {
    return MCELL_FAIL;
  }

  mcell_delete_species_list(reactants);
  mcell_delete_species_list(products);

  int recycled1 = 0;
  int k = 0;
  struct product *prod = NULL;

  //store this nfsim path information
  counter = 0;

  //XXX: it might be necessary to clear out the stuff in pathp at the very end
  for (prod = pathp->product_head; counter < results->numOfResults;) {
    rx->nfsim_players[path][counter] = prod->prod;
    ++counter;
    if(counter < results->numOfResults)
      prod=prod->next;
  }
  prod->next = NULL;

  if(rx->players != NULL)
    free(rx->players);
    free(rx->geometries);

  rx->players = CHECKED_MALLOC_ARRAY(struct species *, num_players,
                               "reaction players array");
  rx->geometries = CHECKED_MALLOC_ARRAY(short, num_players,
                                        "reaction geometries array");


  memset(rx->players, 0, sizeof(struct species*) * num_players);
  memset(rx->geometries, 0, sizeof(short) * num_players);

  rx->players[0] = pathp->reactant1;

  //if its a bimolecular reaction
  if(reac2 != NULL)
    rx->players[1] = pathp->reactant2;

  for (int n_pathway = 0; n_pathway < rx->n_pathways; n_pathway++) {
    k = rx->product_idx[n_pathway] + rx->n_reactants;
    counter = 0;
  for(counter=0;counter<rx->product_idx_aux[n_pathway];counter++){
    //XXX: right now we are ignoring recycled species which is inneficient
    //if (recycled1 == 0 && prod->prod == pathp->reactant1) {
    //  recycled1 = 1;
    //  kk = rx->product_idx[path] + 0;
    //} 
    //else {
      kk = k;
      k++;
    //}
    //kk = rx->product_idx[path] + 0;
    rx->players[kk] = rx->nfsim_players[n_pathway][counter];

  }
  //k = rx->product_idx[n_pathway];
  } /* end for (n_pathway = 0, ...) */


  init_reaction_info(rx);
  rx->info[path].pathname = NULL;

  //adjust reaction probabilities
  //adjust_rates_nfsim(world, rx, pathp);

  struct product *tmp = pathp->product_head;
  struct product *tmp2 = NULL;
  while(tmp != NULL){
    tmp2 = tmp->next;
    free(tmp);
    tmp = tmp2;
  }
  //free(pathp->product_head);
  free(pathp);


  return MCELL_SUCCESS;

} 

int outcome_unimolecular_nfsim(struct volume *world, struct rxn *rx, int path,
                         struct abstract_molecule *reac, double t){
    int result = RX_A_OK;

    if(rx->product_idx_aux[path] == -1){
      mcell_log("uni restart %s\n",reac->graph_pattern);
      queryOptions options = initializeNFSimQueryforUnimolecularFiring(reac, 
                                          rx->external_reaction_names[path]);

      queryResults results = initAndQuerySystemStatus_c(options);

      constructNauty_c(reac->graph_pattern, -1);
      //fill in the rxn react structure with the appropiate information
      prepare_reaction_nfsim(world, rx, &results, path, reac, NULL);
    } 
    return result;
}

int outcome_nfsim(struct volume *world, struct rxn *rx, int path,
                         struct abstract_molecule *reac, struct abstract_molecule *reac2, double t){
    int result = RX_A_OK;
    queryOptions options;
    //if we don't have previous information about this path then build up the rxn structure
    if(rx->product_idx_aux[path] == -1){
      if(reac2 == NULL)
        options = initializeNFSimQueryforUnimolecularFiring(reac, 
                                            rx->external_reaction_names[path]);
      else{
        options = initializeNFSimQueryforBimolecularFiring(reac, reac2,
                                            rx->external_reaction_names[path]);


      }


      queryResults results = initAndQuerySystemStatus_c(options);

      //frees up the option object
      free(options.optionValues[0]);
      free(options.optionValues[1]);
      free(options.optionValues);
      free(options.initKeys);

      //fill in the rxn react structure with the appropiate information
      prepare_reaction_nfsim(world, rx, &results, path, reac, reac2);

      //frees up the query result
      //for(int i=0; i <results.numOfResults; i++)
      //  free(results.results[i]);
      free(results.results);

    }
    //otherwise just update populations
    else{
      for(int i=0; i< rx->product_idx_aux[path]; i++){
        constructNauty_c(rx->product_graph_pattern[path][i],1);
      }
      //and clean the info information for good measure
      rx->info[path].pathname = NULL;

    }
    //also decrease reactant populations
    constructNauty_c(reac->graph_pattern, -1);
    if(reac2 != NULL)
      constructNauty_c(reac2->graph_pattern, -1);


    return result;
}
