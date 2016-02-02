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

    options.optionValues[0] = CHECKED_MALLOC_ARRAY(char, strlen("complex"), "reaction which we will step over");
    strcpy(options.optionValues[0], "complex");

    //optionValues[0] = "complex";

    options.optionValues[1] = CHECKED_MALLOC_ARRAY(char, strlen(external_path), "reaction which we will step over");
    strcpy(options.optionValues[1], external_path);

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


int outcome_unimolecular_nfsim(struct volume *world, struct rxn *rx, int path,
                         struct abstract_molecule *reac, double t){
    int result = RX_A_OK;

    if(rx->product_idx_aux[path] == 0){
      queryOptions options = initializeNFSimQueryforUnimolecularFiring(reac, 
                                          rx->external_reaction_names[path]);

      queryResults results = initAndQuerySystemStatus_c(options);

      constructNauty_c(reac->graph_pattern, -1);
      //determine if reactant is part of product
      int overlapFlag = 0;
      char* product_pattern;

      for(int productIdx = 0; productIdx < results.numOfResults; productIdx++){
        product_pattern = results.results[productIdx];
        constructNauty_c(product_pattern, 1);
        //TODO: we are ignoring optimizing for overlaps for now
        /*if(strcmp(reac->graph_pattern, product_pattern) == 0 && overlapFlag == 0){
          overlapFlag = 1;
        }
        else{*/

          ++rx->product_idx_aux[path];
        //}

      }
      rx->product_graph_pattern[path] = CHECKED_MALLOC_ARRAY(char*,rx->product_idx_aux[path],
                                          "graph patterns for products that have been added to the system");
      overlapFlag = 0;
      int counter = 0;

      for(int productIdx = 0; productIdx < results.numOfResults; productIdx++){
        product_pattern = results.results[productIdx];

        if(strcmp(reac->graph_pattern, product_pattern) == 0 && overlapFlag == 0){
          overlapFlag = 1;
        }
        else{
          rx->product_graph_pattern[path][counter] = product_pattern;
          counter++;
        }
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
      memset(pathp, 0, sizeof(struct pathway));

      /* Scan reactants, copying into the new pathway */
      int num_vol_mols = 0;
      int num_surface_mols = 0;
      int all_3d = 1;
      int complex_type = 0;
      int reactant_idx = 0;
      int oriented_count = 0;
      int num_complex_reactants = 0;

      //obtain the generic vol nfsim reactant
      struct sym_table* nfsim_vol_molecule = reac->properties->sym;
      struct mcell_species *reactants = mcell_add_to_species_list(nfsim_vol_molecule, true, 1, 0, NULL);

      //create list of products and add it to a list
      struct mcell_species *products =
          mcell_add_to_species_list(nfsim_vol_molecule, true, -1, 0, NULL);

      for(int i =1; i <results.numOfResults; i++){
          products = mcell_add_to_species_list(nfsim_vol_molecule, true, -1, 0, products);
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


      //recreate player array
      //struct species ** players = CHECKED_MALLOC_ARRAY(struct species *, num_players,
      //                             "reaction players array");



      int recycled1 = 0;
      int k = 0;
      struct product *prod = NULL;
      //XXX:I actually need to keep the previous path information if im going to 
      //want not having to query the damn reaction everytime but lets leave it
      //like this for now
      rx->players = CHECKED_MALLOC_ARRAY(struct species *, num_players,
                                   "reaction players array");

      rx->players[0] = pathp->reactant1;

      //normally we would iterate over all pathways but right now we only care about one of them
      //for (int n_pathway = 0; path != NULL; n_pathway++, path = path->next) {

      k = rx->product_idx[path] + rx->n_reactants;
  
      for (prod = pathp->product_head; prod != NULL; prod = prod->next) {
        //if (recycled1 == 0 && prod->prod == pathp->reactant1) {
        //  recycled1 = 1;
        //  kk = rx->product_idx[path] + 0;
        //} 
        //else {
          kk = k;
          k++;
        //}
        //kk = rx->product_idx[path] + 0;
        rx->players[kk] = prod->prod;

      }
      k = rx->product_idx[path];
      //if (recycled1 == 0)
        rx->players[k] = NULL;
    //} /* end for (n_pathway = 0, ...) */

      //set reactants and products associated with each 
      //players[0] = path->reactant1;

      init_reaction_info(rx);
      rx->info[path].pathname = NULL;

    } 
    return result;
}