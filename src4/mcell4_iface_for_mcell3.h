/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_MCELL4_IFACE_FOR_MCELL3_H_
#define SRC4_MCELL4_IFACE_FOR_MCELL3_H_

#include "mcell_structs_shared.h"

bool mcell4_convert_mcell3_volume(struct volume* s);
bool mcell4_run_simulation(const bool dump_initial_state, const bool dump_with_geometry = false);
void mcell4_convert_to_data_model(const bool only_for_viz);
void mcell4_delete_world();

#endif // SRC4_MCELL4_IFACE_FOR_MCELL3_H_
