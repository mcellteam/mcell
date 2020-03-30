#include "nfsim_func.h"
#include "map_c.h"

static map_t graph_reaction_map = NULL;

void initialize_rxn_diffusion_functions(struct rxn *mol_ptr);

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

void initialize_diffusion_function(struct abstract_molecule *mol_ptr) {
  // if nfsim provides spatial parameters...
  if (mol_ptr->properties->flags & EXTERNAL_SPECIES &&
      mol_ptr->graph_data->graph_diffusion != -1.0) {
    mol_ptr->get_diffusion = &get_nfsim_diffusion;
    mol_ptr->get_space_step = &get_nfsim_space_step;
    mol_ptr->get_time_step = &get_nfsim_time_step;
  } else {
    mol_ptr->get_diffusion = &get_standard_diffusion;
    mol_ptr->get_space_step = &get_standard_space_step;
    mol_ptr->get_time_step = &get_standard_time_step;
  }

  if (mol_ptr->properties->flags & EXTERNAL_SPECIES &&
      mol_ptr->graph_data->flags >= 0) {
    mol_ptr->get_flags = &get_nfsim_flags;
  } else {
    mol_ptr->get_flags = &get_standard_flags;
  }
}

void initialize_rxn_diffusion_functions(struct rxn *mol_ptr) {
  if (mol_ptr->players[0]->flags & EXTERNAL_SPECIES) {
    mol_ptr->get_reactant_diffusion = rxn_get_nfsim_diffusion;
    mol_ptr->get_reactant_space_step = rxn_get_nfsim_space_step;
    mol_ptr->get_reactant_time_step = rxn_get_nfsim_time_step;
  } else {
    mol_ptr->get_reactant_diffusion = rxn_get_standard_diffusion;
    mol_ptr->get_reactant_space_step = rxn_get_standard_space_step;
    mol_ptr->get_reactant_time_step = rxn_get_standard_time_step;
  }
}

u_int get_standard_flags(void *mol_ptr) {
  struct abstract_molecule *am = (struct abstract_molecule *)mol_ptr;
  return am->properties->flags;
}

u_int get_nfsim_flags(void *mol_ptr) {
  struct abstract_molecule *am = (struct abstract_molecule *)mol_ptr;
  if (am->graph_data->flags < 0)
    return get_standard_flags(am);

  return (unsigned int)am->graph_data->flags;
}

double get_standard_diffusion(void *mol_ptr) {
  struct abstract_molecule *am = (struct abstract_molecule *)mol_ptr;
  return am->properties->D;
}

double get_nfsim_diffusion(void *mol_ptr) {
  // nfsim returns diffusion -1 when the user didnt define any diffusion
  // functions
  struct abstract_molecule *am = (struct abstract_molecule *)mol_ptr;
  if (am->graph_data->graph_diffusion > 0) {
    return am->graph_data->graph_diffusion;
  }
  return get_standard_diffusion(am);
}

double get_standard_space_step(void *mol_ptr) {
  struct abstract_molecule *am = (struct abstract_molecule *)mol_ptr;
  return am->properties->space_step;
}

double get_nfsim_space_step(void *mol_ptr) {
  // nfsim returns diffusion -1 when the user didnt define any diffusion
  // functions
  struct abstract_molecule *am = (struct abstract_molecule *)mol_ptr;
  if (am->graph_data->graph_diffusion >= 0)
    return am->graph_data->space_step;
  return get_standard_space_step(am);
}

double get_standard_time_step(void *mol_ptr) {
  struct abstract_molecule *am = (struct abstract_molecule *)mol_ptr;
  return am->properties->time_step;
}

double get_nfsim_time_step(void *mol_ptr) {
  // nfsim returns diffusion -1 when the user didnt define any diffusion
  // functions
  struct abstract_molecule *am = (struct abstract_molecule *)mol_ptr;
  if (am->graph_data->graph_diffusion >= 0)
    return am->graph_data->time_step;
  return get_standard_time_step(am);
}

double rxn_get_standard_diffusion(struct rxn *mol_ptr, int index) {
  return mol_ptr->players[index]->D;
}

double rxn_get_nfsim_diffusion(struct rxn *mol_ptr, int index) {
  if (mol_ptr->reactant_graph_data &&
      mol_ptr->reactant_graph_data[index]->graph_diffusion >= 0) {
    return mol_ptr->reactant_graph_data[index]->graph_diffusion;
  }
  return rxn_get_standard_diffusion(mol_ptr, index);
}

double rxn_get_standard_time_step(struct rxn *mol_ptr, int index) {
  return mol_ptr->players[index]->time_step;
}

double rxn_get_nfsim_time_step(struct rxn *mol_ptr, int index) {
  if (mol_ptr->reactant_graph_data &&
      mol_ptr->reactant_graph_data[index]->graph_diffusion >= 0) {
    return mol_ptr->reactant_graph_data[index]->time_step;
  }
  return rxn_get_standard_time_step(mol_ptr, index);
}

double rxn_get_standard_space_step(struct rxn *mol_ptr, int index) {
  return mol_ptr->players[index]->space_step;
}

double rxn_get_nfsim_space_step(struct rxn *mol_ptr, int index) {
  if (mol_ptr->reactant_graph_data &&
      mol_ptr->reactant_graph_data[index]->graph_diffusion >= 0) {
    return mol_ptr->reactant_graph_data[index]->space_step;
  }
  return rxn_get_standard_space_step(mol_ptr, index);
}

void initialize_graph_hashmap() { graph_reaction_map = hashmap_new(); }

int get_graph_data(unsigned long graph_pattern_hash,
                   struct graph_data **graph_data) {
  return hashmap_get_nohash(graph_reaction_map, graph_pattern_hash,
                            graph_pattern_hash, (void **)(graph_data));
}

int store_graph_data(unsigned long graph_pattern_hash,
                     struct graph_data *graph_data) {
  return hashmap_put_nohash(graph_reaction_map, graph_pattern_hash,
                            graph_pattern_hash, graph_data);
}
