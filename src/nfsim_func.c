#include "nfsim_func.h"
#include "hashmap.h"

static map_t graph_reaction_map = NULL;

double get_standard_diffusion(struct abstract_molecule* self){
    return self->properties->D;
}


double get_nfsim_diffusion(struct abstract_molecule* self){
    return 0.0;
}


void initialize_graph_hashmap(){
    graph_reaction_map = hashmap_new();
}

int get_graph_data(unsigned long graph_pattern_hash, struct graph_data* graph_data){
    return hashmap_get_nohash(graph_reaction_map, graph_pattern_hash, graph_pattern_hash, (void**)(&graph_data));

}

int store_graph_data(unsigned long graph_pattern_hash, struct graph_data* graph_data){
    hashmap_put_nohash(graph_reaction_map, graph_pattern_hash, graph_pattern_hash, graph_data);
}
