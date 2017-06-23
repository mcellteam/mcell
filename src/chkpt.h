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

#include <stdio.h>

#include "mcell_structs.h"

/* header file for chkpt.c, MCell checkpointing functions */

int create_chkpt(struct volume *world, char const *filename);
int write_chkpt(struct volume *world, FILE *fs);
int read_chkpt(struct volume *world, FILE *fs);
void chkpt_signal_handler(int signo);

int set_checkpoint_state(struct volume *world);

double compute_scaled_time(struct volume *world, double real_time);

unsigned long long
count_items_in_scheduler(struct storage_list *storage_head);
