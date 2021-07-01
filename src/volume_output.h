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

int update_volume_output(struct volume *wrld, struct volume_output_item *vo);
int output_volume_output_item(struct volume *wrld, char const *filename,
                              struct volume_output_item *vo);
