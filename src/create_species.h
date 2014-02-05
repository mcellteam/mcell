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

#ifndef CREATE_SPECIES_H
#define CREATE_SPECIES_H
#include "libmcell.h"

struct species_list_item
{
  struct species_list_item *next;
  struct species *spec;
};

struct species_list
{
  struct species_list_item *species_head;
  struct species_list_item *species_tail;
  int                       species_count;
};

// It might make more sense to put the following structs somewhere else.

struct species_opt_orient
{
  struct species_opt_orient *next;
  struct sym_table *mol_type;
  short orient_set;
  short orient;
  short is_subunit;
};

struct species_opt_orient_list
{
  struct species_opt_orient *mol_type_head;
  struct species_opt_orient *mol_type_tail;
};

/* These are the functions used to create a new species and were adapted from
 * their original use in the parser. Now, the parser versions are just thin
 * wrappers around these. */

// assemble_mol_species is used by both the parser and the API (via
// mcell_create_species).*/
struct species *assemble_mol_species(MCELL_STATE* state,
                                     struct sym_table *sym,
                                     double D_ref,
                                     double D,
                                     int is_2d,
                                     double custom_time_step,
                                     int target_only,
                                     double max_step_length);
// The following functions are *only* used by the parser
int add_to_species_list(struct mem_helper *species_list_mem,
                        struct species_list *list,
                        struct species *spec);
void print_species_summary(MCELL_STATE* state, struct species *mol);
void print_species_summaries(MCELL_STATE* state,
                             struct species_list_item *mols);
#endif
