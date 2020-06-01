/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#ifndef SRC4_MCELL4_IFACE_FOR_MCELL3_H_
#define SRC4_MCELL4_IFACE_FOR_MCELL3_H_

#include "mcell_structs.h"

bool mcell4_convert_mcell3_volume(struct volume* s);
bool mcell4_run_simulation(const bool dump_initial_state);
void mcell4_convert_to_data_model();
void mcell4_delete_world();

#endif // SRC4_MCELL4_IFACE_FOR_MCELL3_H_
