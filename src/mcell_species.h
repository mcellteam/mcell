/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#pragma once

#include <stdbool.h>

//#include "mcell_engine.h"
#include "mcell_init.h"
#include "mcell_structs.h"

typedef struct sym_entry mcell_symbol;

struct mcell_species_spec {
  const char *name;
  double D;
  int is_2d;               // 3D = 0; 2D = 1
  double custom_time_step; // default is 0.0
  int target_only;         // default is 0
  double max_step_length;  // default is 0.0
  double space_step;
  /*
  * external species denotes a species whose participation in the system will 
  * be queried to an external algorithm (nfsim for now)
  */
  int external_species; //default is 0
};

struct mcell_species {
  struct mcell_species *next;
  struct sym_entry *mol_type;
  short orient_set;
  short orient;
};

struct mcell_species_list {
  struct mcell_species *mol_type_head;
  struct mcell_species *mol_type_tail;
};

MCELL_STATUS mcell_create_species(MCELL_STATE *state,
                                  struct mcell_species_spec *species,
                                  mcell_symbol **species_ptr);

struct mcell_species *
mcell_add_to_species_list(mcell_symbol *species_ptr, bool is_oriented,
                          int orientation, struct mcell_species *species_list);

void mcell_delete_species_list(struct mcell_species *species);

int new_mol_species(MCELL_STATE *state, const char *name, struct sym_entry **sym_ptr);
