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

using namespace std;

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

  save_subsystem(ctx);

  string geometry_objects_name = save_geometry(ctx);

  save_instantiation(ctx, geometry_objects_name);

  // observables

  // molecules
  // - volume
  // - surface

  // model
  // - config
  // - warnings
  // - notifications
}


std::string PythonExporter::save_subsystem(PythonExportContext& ctx) {
  ofstream out;
  open_and_check_file(SUBSYSTEM, out);
  out << MCELL_IMPORT;

  string res_name = model->Subsystem::export_to_python(out, ctx);

  out.close();
  return res_name;
}


std::string PythonExporter::save_geometry(PythonExportContext& ctx) {
  ofstream out;
  open_and_check_file(GEOMETRY, out);
  out << MCELL_IMPORT;

  string res_name = model->export_vec_geometry_objects(out, ctx, "");

  out.close();
  return res_name;
}


std::string PythonExporter::save_instantiation(PythonExportContext& ctx, const std::string& geometry_objects_name) {
  // prints out everything, even past releases
  // for checkpointing, we always need to fully finish the current iteration and then start the new one
  std::ofstream out;
  open_and_check_file(INSTANTIATION, out);
  out << MCELL_IMPORT;
  out << get_import_star(GEOMETRY);
  out << "\n";

  string release_sites_name = model->export_vec_release_sites(out, ctx, "");

  gen_ctor_call(out, INSTANTIATION, NAME_CLASS_INSTANTIATION);
  gen_param_id(out, NAME_RELEASE_SITES, release_sites_name, true);
  gen_param_id(out, NAME_GEOMETRY_OBJECTS, geometry_objects_name, false);
  out << CTOR_END;

  out.close();
  return INSTANTIATION;
}

} // namespace API
} // namespace MCell

