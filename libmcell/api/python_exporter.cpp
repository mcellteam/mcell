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

#include "api/python_exporter.h"
#include "api/python_export_constants.h"
#include "api/python_export_utils.h"
#include "api/model.h"
#include "world.h"
#include "src/util.h"

namespace MCell {
namespace API {

PythonExporter::PythonExporter(Model* model_) :
  model(model_) {

  assert(model != nullptr);
  world = model->get_world();
  assert(world != nullptr);
}


void PythonExporter::open_and_check_file(
    const std::string file_name, std::ofstream& out,
    const bool for_append,
    const bool bngl) {

  open_and_check_file_w_prefix(output_dir, file_name, out, for_append, bngl);
}


void PythonExporter::save_checkpoint(const std::string& output_dir_) {
  output_dir = output_dir_;
  if (output_dir.back() != BNG::PATH_SEPARATOR) {
    output_dir += BNG::PATH_SEPARATOR;
  }

  ::make_parent_dir(output_dir.c_str());

  PythonExportContext ctx;

  // parameters
  // - includes rng state

  // subsystem
  save_subsystem(ctx);

  // geometry

  // instantiation

  // observables

  // molecules
  // - volume
  // - surface

  // model
  // - config
  // - warnings
  // - notifications
}


void PythonExporter::save_subsystem(PythonExportContext& ctx) {
  std::ofstream out_subsystem;
  open_and_check_file(SUBSYSTEM, out_subsystem);
  out_subsystem << MCELL_IMPORT;

  model->Subsystem::export_to_python(out_subsystem, ctx);
  out_subsystem.close();
}

} // namespace API
} // namespace MCell

