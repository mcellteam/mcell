#include "nfsim_func.h"
#include "hashmap.h"

static map_t graph_reaction_map = NULL;

void initialize_diffusion_function(struct abstract_molecule* this);
void initialize_rxn_diffusion_functions(struct rxn* this);

double get_standard_diffusion(struct abstract_molecule* self);
double get_nfsim_diffusion(struct abstract_molecule* self);
double get_standard_space_step(struct abstract_molecule* self);
double get_nfsim_space_step(struct abstract_molecule* self);
double get_standard_time_step(struct abstract_molecule* self);
double get_nfsim_time_step(struct abstract_molecule* self);

double rxn_get_nfsim_diffusion(struct rxn*, int);
double rxn_get_standard_diffusion(struct rxn*, int);
double rxn_get_nfsim_space_step(struct rxn*, int);
double rxn_get_standard_time_step(struct rxn*, int);
double rxn_get_nfsim_time_step(struct rxn*, int);
double rxn_get_standard_space_step(struct rxn*, int);


void initialize_diffusion_function(struct abstract_molecule* this){
    //if nfsim provides spatial parameters...
    if(this->properties->flags & EXTERNAL_SPECIES && this->graph_data->graph_diffusion != -1.0){
        this->get_diffusion = get_nfsim_diffusion;
        this->get_space_step = get_nfsim_space_step;
        this->get_time_step = get_nfsim_time_step;
    }
    else{
        this->get_diffusion = get_standard_diffusion;
        this->get_space_step = get_standard_space_step;
        this->get_time_step = get_standard_time_step;
    }

    if(this->properties->flags & EXTERNAL_SPECIES && this->graph_data->flags >= 0){
        this->get_flags = get_nfsim_flags;
    }
    else{
        this->get_flags = get_standard_flags;
        
    }



}

void initialize_rxn_diffusion_functions(struct rxn* this){
    if(this->players[0]->flags & EXTERNAL_SPECIES){
        this->get_reactant_diffusion = rxn_get_nfsim_diffusion;
        this->get_reactant_space_step = rxn_get_nfsim_space_step;
        this->get_reactant_time_step = rxn_get_nfsim_time_step;
    }
    else{
        this->get_reactant_diffusion = rxn_get_standard_diffusion;
        this->get_reactant_space_step = rxn_get_standard_space_step;
        this->get_reactant_time_step = rxn_get_standard_time_step;
    }
}

u_int get_standard_flags(struct abstract_molecule* this){
    return this->properties->flags;
}

u_int get_nfsim_flags(struct abstract_molecule* this){
    if(this->graph_data->flags <0)
        return get_standard_flags(this);

    return (unsigned int) this->graph_data->flags;
}

double get_standard_diffusion(struct abstract_molecule* this){
    return this->properties->D;
}

double get_nfsim_diffusion(struct abstract_molecule* this){
    //nfsim returns diffusion -1 when the user didnt define any diffusion functions
    if(this->graph_data->graph_diffusion > 0){
        return this->graph_data->graph_diffusion;
    }
    return get_standard_diffusion(this);
}

double get_standard_space_step(struct abstract_molecule* this){
    return this->properties->space_step;
}

double get_nfsim_space_step(struct abstract_molecule* this){
    //nfsim returns diffusion -1 when the user didnt define any diffusion functions
    if(this->graph_data->graph_diffusion >= 0)
        return this->graph_data->space_step;
    return get_standard_space_step(this);
}


double get_standard_time_step(struct abstract_molecule* this){
    return this->properties->time_step;
}

double get_nfsim_time_step(struct abstract_molecule* this){
    //nfsim returns diffusion -1 when the user didnt define any diffusion functions
    if(this->graph_data->graph_diffusion >= 0)
        return this->graph_data->time_step;
    return get_standard_time_step(this);
}


double rxn_get_standard_diffusion(struct rxn* this, int index){
    return this->players[index]->D;
}

double rxn_get_nfsim_diffusion(struct rxn* this, int index){
    if(this->reactant_graph_data && this->reactant_graph_data[index]->graph_diffusion >=0){
        return this->reactant_graph_data[index]->graph_diffusion;
    }
    return rxn_get_standard_diffusion(this, index);
}


double rxn_get_standard_time_step(struct rxn* this, int index){
    return this->players[index]->time_step;
}

double rxn_get_nfsim_time_step(struct rxn* this, int index){
    if(this->reactant_graph_data && this->reactant_graph_data[index]->graph_diffusion >=0){
        return this->reactant_graph_data[index]->time_step;
    }
    return rxn_get_standard_time_step(this, index);
}


double rxn_get_standard_space_step(struct rxn* this, int index){
    return this->players[index]->space_step;
}


double rxn_get_nfsim_space_step(struct rxn* this, int index){
    if(this->reactant_graph_data && this->reactant_graph_data[index]->graph_diffusion >=0){
        return this->reactant_graph_data[index]->space_step;
    }
    return rxn_get_standard_space_step(this, index);
}


void initialize_graph_hashmap(){
    graph_reaction_map = hashmap_new();
}

int get_graph_data(unsigned long graph_pattern_hash, struct graph_data** graph_data){
    return hashmap_get_nohash(graph_reaction_map, graph_pattern_hash, graph_pattern_hash, (void**)(graph_data));

}

int store_graph_data(unsigned long graph_pattern_hash, struct graph_data* graph_data){
    return hashmap_put_nohash(graph_reaction_map, graph_pattern_hash, graph_pattern_hash, graph_data);
}
