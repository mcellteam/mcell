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

#include <stdio.h>

#include "mcell_structs.h"

/* header file for chkpt.c, MCell checkpointing functions */

int create_chkpt(struct volume *world, char const *filename);
int write_chkpt(struct volume *world, FILE *fs);
int read_chkpt(struct volume *world, FILE *fs, bool only_time_and_iter);
void chkpt_signal_handler(int signo);

int set_checkpoint_state(struct volume *world);

double compute_scaled_time(struct volume *world, double real_time);

unsigned long long
count_items_in_scheduler(struct storage_list *storage_head);
