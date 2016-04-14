#include "nfsim_func.h"
#include "hashmap.h"

static map_t graph_reaction_map = NULL;


void initialize_diffusion_function(struct abstract_molecule* this){
    if(this->properties->flags & EXTERNAL_SPECIES){
        this->get_diffusion = get_nfsim_diffusion;
        this->get_space_step = get_nfsim_space_step;
        this->get_time_step = get_nfsim_time_step;
    }
    else{
        this->get_diffusion = get_standard_diffusion;
        this->get_space_step = get_standard_space_step;
        this->get_time_step = get_standard_time_step;

    }
}

double get_standard_diffusion(struct abstract_molecule* this){
    return this->properties->D;
}

double get_nfsim_diffusion(struct abstract_molecule* this){
    //nfsim returns diffusion -1 when the user didnt define any diffusion functions
    if(this->graph_data->graph_diffusion > 0)
        return this->graph_data->graph_diffusion;
    return get_standard_diffusion(this);
}

double get_standard_space_step(struct abstract_molecule* this){
    return this->properties->space_step;
}

double get_nfsim_space_step(struct abstract_molecule* this){
    //nfsim returns diffusion -1 when the user didnt define any diffusion functions
    if(this->graph_data->space_step >= 0)
        return this->graph_data->space_step;
    return get_standard_space_step(this);
}


double get_standard_time_step(struct abstract_molecule* this){
    return this->properties->time_step;
}

double get_nfsim_time_step(struct abstract_molecule* this){
    //nfsim returns diffusion -1 when the user didnt define any diffusion functions
    if(this->graph_data->time_step >= 0)
        return this->graph_data->time_step;
    return get_standard_time_step(this);
}


double rxn_get_standard_diffusion(struct rxn* this, int index){
    return this->players[index]->D;
}

double rxn_get_nfsim_diffusion(struct rxn* this, int index){
    return this->reactant_graph_data[index]->graph_diffusion;
}


void initialize_graph_hashmap(){
    graph_reaction_map = hashmap_new();
}

int get_graph_data(unsigned long graph_pattern_hash, struct graph_data* graph_data){
    return hashmap_get_nohash(graph_reaction_map, graph_pattern_hash, graph_pattern_hash, (void**)(&graph_data));

}

int store_graph_data(unsigned long graph_pattern_hash, struct graph_data* graph_data){
    return hashmap_put_nohash(graph_reaction_map, graph_pattern_hash, graph_pattern_hash, graph_data);
}
