#include "config.h"

#include <string.h>
#include <stdlib.h>
#include "logging.h"
#include "mcell_structs.h"
#include "react.h"
#include "hashmap.h"
#include "react_util.h"
//#include "lru.h"

map_t reaction_map = NULL;
char reaction_key[300];

struct rxn *rx;



unsigned long
lhash(const char *keystring)
{
    unsigned long key = crc32((unsigned char*)(keystring), strlen(keystring));

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


 queryOptions initializeNFSimQueryForBimolecularReactions(struct abstract_molecule *am, 
                                                          struct abstract_molecule* am2)
 {
    //constant settings

    static const char* optionKeys[1]  = {"numReactants"};
    static const char* optionValues[1] = {"2"};
    static const int optionSeeds[2]= {1,1};
    static const char** speciesArray[2];
    //initialize speciesArray with the string we are going to query
    speciesArray[0] = am->graph_data->graph_pattern;
    speciesArray[1] = am2->graph_data->graph_pattern;
    
    //copy these settings to the options object
    queryOptions options;
    options.initKeys = speciesArray;
    options.initValues = optionSeeds;
    options.numOfInitElements = 2;

    options.optionKeys = optionKeys;
    options.optionValues = optionValues;
    options.numOfOptions = 1;
    return options;
}

/**********************************************************************
 *
 * This function creates a queryOptions object for designing an NFSim experiment query
 *
 * In: The abstract molecule whose nfsim status we are going to verify
 *
 * Out: the queryOptions object we will use to interact with nfsim
 *
 **********************************************************************/
 queryOptions initializeNFSimQueryForUnimolecularReactions(struct abstract_molecule *am){
    //constant settings
    static const char* optionKeys[1]  = {"numReactants"};
    static const char* optionValues[1] = {"1"};
    static const int optionSeeds[1]= {1};
    static char** speciesArray[1];
    //initialize speciesArray with the string we are going to query
    //const char** speciesArray = CHECKED_MALLOC_ARRAY(char*, 1, "string array of patterns");
    speciesArray[0] = am->graph_data->graph_pattern;
    //copy these settings to the options object
    queryOptions options;
    options.initKeys = speciesArray;
    options.initValues = optionSeeds;
    options.numOfInitElements = 1;

    options.optionKeys = optionKeys;
    options.optionValues = optionValues;
    options.numOfOptions = 1;
    //free(speciesArray);
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
int trigger_bimolecular_preliminary_nfsim(struct abstract_molecule *reacA,
                                    struct abstract_molecule *reacB) {

    if (reaction_map == NULL)
        reaction_map = hashmap_new();

    unsigned long reactionHash = reacA->graph_data->graph_pattern_hash + reacB->graph_data->graph_pattern_hash;
    int error = hashmap_get_nohash(reaction_map, reactionHash, reactionHash, (void**)(&rx));
    //error = find_in_cache(reaction_key, rx);

    //XXX: it might be worth it to return the rx object since we already queried it

    if(error == MAP_OK){
        if (rx != NULL)
            return 1;
        return 0;
    }
    //if it doesn't exist in the hashmap yet we have to ask nfsim
    queryOptions options = initializeNFSimQueryForBimolecularReactions(reacA, reacB);

    

    reactantQueryResults query2 = initAndQueryByNumReactant_c(options);
    //if we have reactions great! although we can't decide what to do with them yet. bummer.
    if(query2.numOfResults > 0){
        return 1;
    }
    else{
        //if we know there's no reactions there's no need to check again later
        error = hashmap_put_nohash(reaction_map, reactionHash, reactionHash, NULL);
        return 0;
    }
    
}

/***********

********/
int trigger_bimolecular_nfsim(struct volume* state, struct abstract_molecule *reacA,
                        struct abstract_molecule *reacB,short orientA,
                        short orientB,  struct rxn **matching_rxns){


    int error;
    int num_matching_rxns = 0; 

    if (reaction_map == NULL)
        reaction_map = hashmap_new();

    //memset(&reaction_key[0], 0, sizeof(reaction_key));
    rx = NULL;

    unsigned long reactionHash = reacA->graph_data->graph_pattern_hash + reacB->graph_data->graph_pattern_hash;
    //sprintf(reaction_key,"%lu",reacA->graph_pattern_hash + reacB->graph_pattern_hash);
    //mcell_log("reaction_key %s %s %s",reacA->graph_pattern, reacB->graph_pattern, reaction_key);
    //A + B ->C is the same as B +A -> C
    //if (strcmp(reacA->graph_pattern, reacB->graph_pattern) < 0)
    //    sprintf(reaction_key,"%s-%s",reacA->graph_pattern,reacB->graph_pattern);
    //else
    //    sprintf(reaction_key,"%s-%s",reacB->graph_pattern,reacA->graph_pattern);
    

    error = hashmap_get_nohash(reaction_map, reactionHash, reactionHash, (void**)(&rx));
    //error = find_in_cache(reaction_key, rx);

    if(error == MAP_OK){
    //if(error != -1){

        if(rx != NULL){
            int result = process_bimolecular(reacA, reacB, rx, orientA, orientB,
                            matching_rxns, num_matching_rxns, 0);

            if(result == 1)
                num_matching_rxns++;
        }
        return num_matching_rxns;
    }

    //mcell_log("+++++ %s %s %s",reacA->graph_pattern, reacB->graph_pattern, reaction_key);
    queryOptions options = initializeNFSimQueryForBimolecularReactions(reacA, reacB);
    //reset, init, query the nfsim system
    reactantQueryResults query2 = initAndQueryByNumReactant_c(options);

    //XXX: it would probably would be more natural to make a method that queries all the reactions
    // associated with a given species
    if(query2.numOfResults > 0){
        rx = new_reaction();
        initializeNFSimReaction(state, rx, 2, query2, reacA, reacB);

        int result = process_bimolecular(reacA, reacB, rx, orientA, orientB,
                        matching_rxns, num_matching_rxns, 0);

        if(result == 1)
            num_matching_rxns++;
    }
    //store value in hashmap

    error = hashmap_put_nohash(reaction_map, reactionHash, reactionHash, rx);
    //add_to_cache(reaction_key, rx);

    //CLEANUP
    delete_reactantQueryResults(query2);
    return num_matching_rxns;



}

/*


is_surface: true if there is a surface reactant
*/
int adjust_rates_nfsim(struct volume* state, struct rxn *rx, bool is_surface){
    double pb_factor = 0.0;
    //int max_num_surf_products = set_product_geometries(path, rx, prod);
    pb_factor = compute_pb_factor(
    state->time_unit, state->length_unit, state->grid_density,
    state->rx_radius_3d/state->r_length_unit,   //transform back to a unitless scale
    &state->rxn_flags,
    &state->create_shared_walls_info_flag,
    rx, is_surface); //max_num_surf_products);
    rx->pb_factor = pb_factor;
    //mcell_log("!!!pb_factor %.10e",pb_factor);

    //JJT: balance out rate (code from scale_rxn_probabilities)
    //extracted because we only want to change the rate for one path

    double rate;
    for(int i =0; i < rx->n_pathways; i++){

        if (!rx->rates || !rx->rates[i]) {
          rate = pb_factor * rx->cum_probs[i];
          //mcell_log("!!%.10e %.10e",rx->cum_probs[i],rate);
        } else
          rate = 0.0;
        rx->cum_probs[i] = rate;
    }

    return 0;

    /*if (scale_rxn_probabilities(&state->reaction_prob_limit_flag, state->notify,
                            path, rx, pb_factor))
      return 1;*/

}



int initializeNFSimReaction(struct volume *state,
                            struct rxn* r, int n_reactants, reactantQueryResults query2, 
                            struct abstract_molecule* reacA, struct abstract_molecule* reacB){
    r->cum_probs = CHECKED_MALLOC_ARRAY(double, query2.numOfAssociatedReactions[0],
                                      "cumulative probabilities");
    r->external_reaction_names = CHECKED_MALLOC_ARRAY(char*, query2.numOfAssociatedReactions[0],
                                      "external reaction names");
    r->n_pathways = query2.numOfAssociatedReactions[0];
    r->product_idx = CHECKED_MALLOC_ARRAY(u_int, query2.numOfAssociatedReactions[0]+1,
                                      "the different possible products");
    r->product_idx_aux = CHECKED_MALLOC_ARRAY(int, query2.numOfAssociatedReactions[0]+1,
                                      "the different possible products");
    //r->product_graph_pattern = CHECKED_MALLOC_ARRAY(char**,query2.numOfAssociatedReactions[0],
    //                                  "graph patterns of the possible products");
    
    r->reactant_graph_data = CHECKED_MALLOC_ARRAY(struct graph_data*, n_reactants,
                                      "graph patterns of the possible products");

    r->reactant_graph_data[0] = reacA->graph_data;
    if(reacB){
        r->reactant_graph_data[1] = reacB->graph_data;
    }

    r->product_graph_data = CHECKED_MALLOC_ARRAY(struct graph_data**,query2.numOfAssociatedReactions[0],
                                      "graph patterns of the possible products");

    r->nfsim_players = CHECKED_MALLOC_ARRAY(struct species**,query2.numOfAssociatedReactions[0],
                                      "products associated to each path");

    r->nfsim_geometries = CHECKED_MALLOC_ARRAY(short*,query2.numOfAssociatedReactions[0],
                                      "geometries associated to each path");

    memset(r->nfsim_players, 0, sizeof(struct species**)*query2.numOfAssociatedReactions[0]);

    r->n_reactants = n_reactants;

    //XXX:do we really have to go over all of them or do we need to filter out repeats?
    //for (int i=0; i<query2.numOfResults; i++){
    int reactionNameLength;
    for(int path=0;path<query2.numOfAssociatedReactions[0]; path++){
        r->cum_probs[path] = query2.associatedReactions[0].rates[path];
        r->external_reaction_names[path] = strdup(query2.associatedReactions[0].reactionNames[path]);
        r->product_idx[path] = 0;
        r->product_idx_aux[path] = -1;
    }


    //XXX: is this necessary?
    //temporarily initialize reaction list with nfsim_players
    r->players = CHECKED_MALLOC_ARRAY(struct species *, n_reactants,
                           "reaction players array");
    r->geometries = CHECKED_MALLOC_ARRAY(short, n_reactants,
                           "geometry players array");

    r->players[0] = reacA->properties;
    r->geometries[0] = 0;
    if(reacB != NULL){
    r->players[1] = reacB->properties;
    r->geometries[1] = 0;
    }

    //nfsim diffusion function, depending on whether the user wants us to do this.
    initialize_rxn_diffusion_functions(r);

    bool orientation_flag1, orientation_flag2 = 0;
    int reactantOrientation1, reactantOrientation2;

    calculate_reactant_orientation(reacA, reacB, &orientation_flag1, &orientation_flag2,
                                   &reactantOrientation1, &reactantOrientation2);
    if(orientation_flag1){
        rx->geometries[0] = reactantOrientation1;
    }
    if(orientation_flag2){
        rx->geometries[1] = reactantOrientation2;
    }

    //adjust reaction probabilities
    //if (reacB != NULL)
    //    mcell_log("++++ %s %s",reacA->graph_data->graph_pattern, reacB->graph_data->graph_pattern);
    //else
    //    mcell_log("---- %s ",reacA->graph_data->graph_pattern);
    adjust_rates_nfsim(state, r, orientation_flag1 & orientation_flag2);

    //calculate cummulative probabilities
    for (int n_pathway = 1; n_pathway < r->n_pathways; ++n_pathway){
      r->cum_probs[n_pathway] += r->cum_probs[n_pathway - 1];
    }

    if (r->n_pathways > 0)
        r->min_noreaction_p = r->max_fixed_p = r->cum_probs[r->n_pathways - 1];
    else
        r->min_noreaction_p = r->max_fixed_p = 1.0;


    return 0;
}

struct rxn *pick_unimolecular_reaction_nfsim(struct volume *state,
                                       struct abstract_molecule *am){

    int error;

    if (reaction_map == NULL)
        reaction_map = hashmap_new();

    //memset(&reaction_key[0], 0, sizeof(reaction_key));
    rx = NULL;
    //sprintf(reaction_key,"%s",am->graph_pattern);
    //sprintf(reaction_key,"%lu",am->graph_pattern_hash);

    
    //check in the hashmap in case this is a reaction we have encountered before
    error = hashmap_get_nohash(reaction_map, am->graph_data->graph_pattern_hash, am->graph_data->graph_pattern_hash, (void**)(&rx));
    //error = find_in_cache(reaction_key, rx);

    if(error == MAP_OK){
        return rx;
    }
    
    //if(error != -1)
    //    return rx;

    //otherwise build the object
    queryOptions options = initializeNFSimQueryForUnimolecularReactions(am);
    //reset, init, query the nfsim system
    reactantQueryResults query2 = initAndQueryByNumReactant_c(options);

    //XXX: it would probably would be more natural to make a method that queries all the reactions
    // associated with a given species
    if(query2.numOfResults > 0){
      rx = new_reaction();
      initializeNFSimReaction(state, rx, 1, query2, am, NULL);

    }
 
    //store newly created reaction in the hashmap
    error = hashmap_put_nohash(reaction_map, am->graph_data->graph_pattern_hash, am->graph_data->graph_pattern_hash, rx);
    //add_to_cache(reaction_key, rx);
 
    //CLEANUP
    delete_reactantQueryResults(query2);

    return rx;

}