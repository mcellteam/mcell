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

#include "mcell3_world_converter.h"
#include "world.h"
#include "datamodel_defines.h"

// holds global instance of world after conversion
MCell::MCell3WorldConverter g_converter;


bool mcell4_convert_mcell3_volume(volume* s) {
  return g_converter.convert(s);
}


bool mcell4_run_simulation(const bool dump_initial_state, const bool dump_with_geometry) {
  g_converter.world->init_and_run_simulation(dump_initial_state, dump_with_geometry);
  return true;
}

void mcell4_convert_to_data_model(const bool only_for_viz) {
  const char* fname;
  if (only_for_viz) {
    fname = DEFAULT_DATA_MODEL_VIZ_FILENAME;
  }
  else {
    fname = DEFAULT_DATA_MODEL_FILENAME;
  }
  g_converter.world->export_data_model(fname, only_for_viz);
  mcell_log("Datamodel was exported to '%s'.", fname);
}

void mcell4_delete_world() {
  g_converter.reset();
}
