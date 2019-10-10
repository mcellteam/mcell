#ifndef NFSIM_FUNC
#define NFSIM_FUNC

#include "mcell_structs.h"

// typedef double (*get_reactant_diffusion)(int a, int b);

void initialize_diffusion_function(struct abstract_molecule *this_ptr);
void initialize_rxn_diffusion_functions(struct rxn *this_ptr);

double get_standard_diffusion(void *self);
double get_nfsim_diffusion(void *self);

double get_standard_space_step(void *self);
double get_nfsim_space_step(void *self);

double get_standard_time_step(void *self);
double get_nfsim_time_step(void *self);

double rxn_get_nfsim_diffusion(struct rxn *, int);
double rxn_get_standard_diffusion(struct rxn *, int);

double rxn_get_nfsim_space_step(struct rxn *, int);
double rxn_get_standard_time_step(struct rxn *, int);

double rxn_get_nfsim_time_step(struct rxn *, int);
double rxn_get_standard_space_step(struct rxn *, int);

void initialize_graph_hashmap(void);
int get_graph_data(unsigned long graph_pattern_hash,
                   struct graph_data **graph_data);
int store_graph_data(unsigned long graph_pattern_hash,
                     struct graph_data *graph_data);

u_int get_nfsim_flags(void *this_ptr);
u_int get_standard_flags(void *this_ptr);

#endif
