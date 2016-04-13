#ifndef NFSIM_FUNC
#define NFSIM_FUNC

#include "mcell_structs.h"

double get_standard_diffusion(struct abstract_molecule* self);
double get_nfsim_diffusion(struct abstract_molecule* self);

int get_graph_data(unsigned long graph_pattern_hash, struct graph_data* graph_data);
int store_graph_data(unsigned long graph_pattern_hash, struct graph_data* graph_data);
void initialize_graph_hashmap();

#endif
