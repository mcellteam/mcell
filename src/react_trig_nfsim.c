#include "config.h"

#include <string.h>
#include <stdlib.h>
#include "logging.h"
#include "mcell_structs.h"
#include "react.h"
#include "hashmap.h"

map_t reaction_map = NULL;
char reaction_key[300];
struct rxn *rx;


 queryOptions initializeNFSimQueryForBimolecularReactions(struct abstract_molecule *am, 
                                                          struct abstract_molecule* am2)
 {
    //constant settings

    static const char* optionKeys[1]  = {"numReactants"};
    static const char* optionValues[1] = {"2"};
    static const int optionSeeds[2]= {1,1};
    static char** speciesArray[2];
    //initialize speciesArray with the string we are going to query
    //const char** speciesArray = CHECKED_MALLOC_ARRAY(char*, 1, "string array of patterns");
    speciesArray[0] = am->graph_pattern;
    speciesArray[1] = am2->graph_pattern;
    
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
    speciesArray[0] = am->graph_pattern;
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
int trigger_bimolecular_preliminary_nfsim(struct volume_molecule *reacA,
                                    struct volume_molecule *reacB) {

    queryOptions options = initializeNFSimQueryForBimolecularReactions(reacA, reacB);

    //reset, init, query the nfsim system

    reactantQueryResults query2 = initAndQueryByNumReactant_c(options);

    if(query2.numOfResults > 0){
        return 1;
    }
    return 0;
}
/***********

********/
int trigger_bimolecular_nfsim(struct abstract_molecule *reacA,
                        struct abstract_molecule *reacB,short orientA,
                        short orientB,  struct rxn **matching_rxns){


    int error;
    int num_matching_rxns = 0; 

    if (reaction_map == NULL)
        reaction_map = hashmap_new();

    memset(&reaction_key[0], 0, sizeof(reaction_key));
    rx = NULL;
    sprintf(reaction_key,"%s-%s",reacA->graph_pattern,reacB->graph_pattern);
    error = hashmap_get(reaction_map, reaction_key, (void**)(&rx));
    if(error == MAP_OK){
        //mcell_log("???? %d %p %s \n",hashmap_hash(reaction_map, reaction_key), rx, reaction_key);

        if(rx != NULL){
            int result = process_bimolecular(reacA, reacB, rx, orientA, orientB,
                            matching_rxns, num_matching_rxns, 0);

            if(result == 1)
                num_matching_rxns++;
        }
        return num_matching_rxns;
    }
    
    queryOptions options = initializeNFSimQueryForBimolecularReactions(reacA, reacB);
    //reset, init, query the nfsim system
    reactantQueryResults query2 = initAndQueryByNumReactant_c(options);

    //XXX: it would probably would be more natural to make a method that queries all the reactions
    // associated with a given species
    if(query2.numOfResults > 0){
        rx = new_reaction();
        initializeNFSimReaction(rx, 2, query2);
        //XXX:I need to conserve path information. THIS IS A HACK
        rx->players = CHECKED_MALLOC_ARRAY(struct species *, 2,
                               "reaction players array");
        rx->geometries = CHECKED_MALLOC_ARRAY(short, 2,
                               "geometry players array");

        rx->players[0] = reacA->properties;
        rx->players[1] = reacB->properties;

        rx->geometries[0] = 0;
        rx->geometries[1] = 0;

        int result = process_bimolecular(reacA, reacB, rx, orientA, orientB,
                        matching_rxns, num_matching_rxns, 0);

        if(result == 1)
            num_matching_rxns++;
    }
    //store value in hashmap
    //mcell_log("!!!!!!!!!!! %d %p %s\n",hashmap_hash(reaction_map, reaction_key), rx, reaction_key);

    error = hashmap_put(reaction_map, reaction_key, rx);

    //CLEANUP
    delete_reactantQueryResults(query2);
    return num_matching_rxns;



}


int initializeNFSimReaction(struct rxn* r, int n_reactants, reactantQueryResults query2){
    r->cum_probs = CHECKED_MALLOC_ARRAY(double, query2.numOfAssociatedReactions[0],
                                      "cumulative probabilities");
    r->external_reaction_names = CHECKED_MALLOC_ARRAY(char*, query2.numOfAssociatedReactions[0],
                                      "external reaction names");
    r->n_pathways = query2.numOfAssociatedReactions[0];
    r->product_idx = CHECKED_MALLOC_ARRAY(u_int, query2.numOfAssociatedReactions[0]+1,
                                      "the different possible products");
    r->product_idx_aux = CHECKED_MALLOC_ARRAY(int, query2.numOfAssociatedReactions[0]+1,
                                      "the different possible products");
    r->product_graph_pattern = CHECKED_MALLOC_ARRAY(char**,query2.numOfAssociatedReactions[0],
                                      "graph patterns of the possible products");


    r->nfsim_players = CHECKED_MALLOC_ARRAY(struct species**,query2.numOfAssociatedReactions[0],
                                      "products associated to each path");
    memset(r->nfsim_players, 0, sizeof(struct species**)*query2.numOfAssociatedReactions[0]);

    r->n_reactants = n_reactants;

    //XXX:do we really have to go over all of them or do we need to filter out repeats?
    //for (int i=0; i<query2.numOfResults; i++){
    int reactionNameLength;
    for(int path=0;path<query2.numOfAssociatedReactions[0]; path++){
        //reactionNameLength = strlen(query2.associatedReactions[0].reactionNames[path]);
        r->cum_probs[path] = query2.associatedReactions[0].rates[path];
        //r->external_reaction_names[path] = CHECKED_MALLOC_ARRAY(char, reactionNameLength+1,
        //                                  "external reaction names");
        //strcpy(r->external_reaction_names[path], query2.associatedReactions[0].reactionNames[path]);
        r->external_reaction_names[path] = strdup(query2.associatedReactions[0].reactionNames[path]);
        r->product_idx[path] = 0;
        r->product_idx_aux[path] = -1;
    }

    //calculate cummulative probabilities
    for (int n_pathway = 1; n_pathway < r->n_pathways; ++n_pathway)
      r->cum_probs[n_pathway] += r->cum_probs[n_pathway - 1];

    if (r->n_pathways > 0)
        r->min_noreaction_p = r->max_fixed_p = r->cum_probs[r->n_pathways - 1];
    else
        r->min_noreaction_p = r->max_fixed_p = 1.0;
}

struct rxn *pick_unimolecular_reaction_nfsim(struct volume *state,
                                       struct abstract_molecule *am){

    int error;

    if (reaction_map == NULL)
        reaction_map = hashmap_new();

    memset(&reaction_key[0], 0, sizeof(reaction_key));
    rx = NULL;
    sprintf(reaction_key,"%s",am->graph_pattern);

    //check in the hashmap in case this is a reaction we have encountered before
    error = hashmap_get(reaction_map, reaction_key, (void**)(&rx));

    if(error == MAP_OK){
        return rx;
    }
    //otherwise build the object
    queryOptions options = initializeNFSimQueryForUnimolecularReactions(am);
    //reset, init, query the nfsim system
    reactantQueryResults query2 = initAndQueryByNumReactant_c(options);

    //XXX: it would probably would be more natural to make a method that queries all the reactions
    // associated with a given species
    if(query2.numOfResults > 0){
      rx = new_reaction();
      initializeNFSimReaction(rx, 1, query2);

    }
    //store newly created reaction in the hashmap
    error = hashmap_put(reaction_map, reaction_key, rx);

    //CLEANUP
    delete_reactantQueryResults(query2);

    return rx;

}