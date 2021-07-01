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

#include "mcell_structs.h"

struct species *new_species(void);
struct geom_object *new_object(void);
struct release_pattern *new_release_pattern(void);
struct rxn *new_reaction(void);
struct rxn_pathname *new_reaction_pathname(void);
struct region *new_region(void);
struct file_stream *new_filestream(void);

int dump_symtab(struct sym_table_head *hashtab);

unsigned long hash(char const *sym);
struct sym_entry *retrieve_sym(char const *sym, struct sym_table_head *hashtab);
struct sym_entry *store_sym(char const *sym, enum symbol_type_t sym_type,
                            struct sym_table_head *hashtab, void *data);
struct sym_table_head *init_symtab(int size);
void destroy_symtab(struct sym_table_head *tab);
