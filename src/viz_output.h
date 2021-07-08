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

#define ASCII_VIZ_EXTERNAL_SPECIES_NAME

#include "mcell_structs.h"

/* Header file for visualization output routines */

int update_frame_data_list(struct volume *world,
                           struct viz_output_block *vizblk);

int init_frame_data_list(struct volume *world, struct viz_output_block *vizblk);

int finalize_viz_output(struct volume *world, struct viz_output_block *vizblk);
