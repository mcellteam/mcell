#ifndef NFSIM_FUNC
#define NFSIM_FUNC

#include "mcell_structs.h"

// typedef double (*get_reactant_diffusion)(int a, int b);

void initialize_diffusion_function(struct abstract_molecule *this);
void initialize_rxn_diffusion_functions(struct rxn *this);

double get_standard_diffusion(struct abstract_molecule *self);
double get_nfsim_diffusion(struct abstract_molecule *self);

double get_standard_space_step(struct abstract_molecule *self);
double get_nfsim_space_step(struct abstract_molecule *self);

double get_standard_time_step(struct abstract_molecule *self);
double get_nfsim_time_step(struct abstract_molecule *self);

double rxn_get_nfsim_diffusion(struct rxn *, int);
double rxn_get_standard_diffusion(struct rxn *, int);

double rxn_get_nfsim_space_step(struct rxn *, int);
double rxn_get_standard_time_step(struct rxn *, int);

double rxn_get_nfsim_time_step(struct rxn *, int);
double rxn_get_standard_space_step(struct rxn *, int);

void initialize_graph_hashmap();
int get_graph_data(unsigned long graph_pattern_hash,
                   struct graph_data **graph_data);
int store_graph_data(unsigned long graph_pattern_hash,
                     struct graph_data *graph_data);

u_int get_nfsim_flags(struct abstract_molecule *this);
u_int get_standard_flags(struct abstract_molecule *this);

#endif
