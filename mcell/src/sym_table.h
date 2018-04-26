/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#pragma once

#include "mcell_structs.h"

#ifndef ISAAC64_H
/* These guys come in for free if we're using Jenkins' random numbers */
typedef uint32_t ub4;      /* unsigned 4-byte quantities */
typedef unsigned char ub1; /* unsigned 1-byte quantities */
#endif

struct species *new_species(void);
struct object *new_object(void);
struct release_pattern *new_release_pattern(void);
struct rxn *new_reaction(void);
struct rxn_pathname *new_reaction_pathname(void);
struct region *new_region(void);
struct file_stream *new_filestream(void);

ub4 jenkins_hash(ub1 *sym, ub4 length);
unsigned long hash(char const *sym);
struct sym_entry *retrieve_sym(char const *sym, struct sym_table_head *hashtab);
struct sym_entry *store_sym(char const *sym, enum symbol_type_t sym_type,
                            struct sym_table_head *hashtab, void *data);
struct sym_table_head *init_symtab(int size);
void destroy_symtab(struct sym_table_head *tab);
