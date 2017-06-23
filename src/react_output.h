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

/* Header file for reaction output routines */

extern int emergency_output_hook_enabled;

void install_emergency_output_hooks(struct volume *world);

int truncate_output_file(char *name, double start_value);

void add_trigger_output(struct volume *world, struct counter *c,
                        struct output_request *ear, int n, short flags);

int flush_reaction_output(struct volume *world);

int check_reaction_output_file(struct output_set *os);

int update_reaction_output(struct volume *world, struct output_block *block);

int write_reaction_output(struct volume *world, struct output_set *set);

struct output_expression *new_output_expr(struct mem_helper *oexpr_mem);
void set_oexpr_column(struct output_expression *oe, struct output_column *oc);
void learn_oexpr_flags(struct output_expression *oe);
struct output_expression *dupl_oexpr_tree(struct output_expression *root,
                                          struct mem_helper *oexpr_mem);
struct output_expression *first_oexpr_tree(struct output_expression *root);
struct output_expression *next_oexpr_tree(struct output_expression *leaf);
void eval_oexpr_tree(struct output_expression *root, int skip_const);
void oexpr_flood_convert(struct output_expression *root, char old_oper,
                         char new_oper);
char *oexpr_title(struct output_expression *root);
