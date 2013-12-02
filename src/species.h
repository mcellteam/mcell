/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by                                                      *
 * The Salk Institute for Biological Studies and                                   *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University                    *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

#include "libmcell.h"
#include "mdlparse_aux.h"

/* These are the functions used to create a new species and were adapted from
 * their original use in the parser. Now, they are used by both the parser and
 * the API (via mcell_create_species). The parser versions are just thin
 * wrappers around these. */

struct sym_table *new_mol_species(MCELL_STATE* state, char *name);
struct species *assemble_mol_species(MCELL_STATE* state,
                                     struct sym_table *sym,
                                     double D_ref,
                                     double D,
                                     int is_2d,
                                     double custom_time_step,
                                     int target_only,
                                     double max_step_length);
int add_to_species_list(struct mem_helper *species_list_mem,
                        struct species_list *list,
                        struct species *spec);
void finish_molecule(MCELL_STATE* state, struct species *mol);
void finish_molecules(MCELL_STATE* state,
                      struct species_list_item *mols);
