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

#include "mcell3_world_converter.h"
#include "world.h"
#include "datamodel_defines.h"

// holds global instance of world after conversion
MCell::MCell3WorldConverter g_converter;


bool mcell4_convert_mcell3_volume(volume* s) {
  return g_converter.convert(s);
}


bool mcell4_run_simulation(const bool dump_initial_state) {
  g_converter.world->run_simulation(dump_initial_state);
  return true;
}

void mcell4_convert_to_datamodel() {
  g_converter.world->export_visualization_datamodel(DEFAULT_DATAMODEL_FILENAME);
  mcell_log("Datamodel was exported to '%s'.", DEFAULT_DATAMODEL_FILENAME);
}

void mcell4_delete_world() {
  g_converter.reset();
}
