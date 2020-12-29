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
#include "api/rng_state.h"
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

  // parameters - later, once we will maintain the association,

  string subsystem_name = save_subsystem(ctx);

  string geometry_objects_name = save_geometry(ctx);

  string instantiation_name = save_instantiation(ctx, geometry_objects_name);

  // need to append - some config flag?
  string observables_name = save_observables(ctx);

  // molecules
  // - volume
  // - surface
  // rng state,
  // starting, ending iterations, ...
  // all generated variables are prefixed with state_
  std::map<std::string, std::string> config_variable_names;
  save_simulation_state(ctx, config_variable_names);

  // model
  // - config
  // - warnings
  // - notifications
  save_model(ctx, subsystem_name, instantiation_name, observables_name, config_variable_names);
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


std::string PythonExporter::save_observables(PythonExportContext& ctx) {
  ofstream out;
  open_and_check_file(OBSERVABLES, out);
  out << MCELL_IMPORT;
  out << get_import_star(SUBSYSTEM);
  out << get_import_star(GEOMETRY);
  out << get_import_star(INSTANTIATION);
  out << "\n";

  string res_name = model->Observables::export_to_python(out, ctx);

  out.close();
  return res_name;
}


// state_variable_names are indexed by Config class attribute names
void PythonExporter::save_simulation_state(
    PythonExportContext& ctx,
    std::map<std::string, std::string>& config_variable_names
) {
  ofstream out;
  open_and_check_file(SIMULATION_STATE, out);
  out << MCELL_IMPORT;
  out << get_import_star(SUBSYSTEM);
  out << get_import_star(GEOMETRY);
  out << "\n";

  // current iteration and time
  gen_assign(out, NAME_INITIAL_ITERATION, world->stats.get_current_iteration());
  config_variable_names[NAME_INITIAL_ITERATION] = NAME_INITIAL_ITERATION;

  gen_assign(out, NAME_INITIAL_TIME, world->stats.get_current_iteration() * world->config.time_unit);
  config_variable_names[NAME_INITIAL_TIME] = NAME_INITIAL_TIME;
  out << "\n";

  // molecules
  save_molecules(ctx, out);

  // rng state
  RngState rng_state = RngState(world->rng);
  config_variable_names[NAME_RNG_STATE] = rng_state.export_to_python(out, ctx);
}


void PythonExporter::save_molecules(PythonExportContext& ctx, std::ostream& out) {


  // prepare species map

  // prepare geometry objects map

  // for each molecule

    // vol

    // surf
}


std::string PythonExporter::save_model(
    PythonExportContext& ctx,
    const std::string& subsystem_name,
    const std::string& instantiation_name,
    const std::string& observables_name,
    const std::map<std::string, std::string>& config_variable_names) {

  // prints out everything, even past releases
  // for checkpointing, we always need to fully finish the current iteration and then start the new one
  std::ofstream out;
  open_and_check_file(MODEL, out);

  // imports
  out << INTERPRETER;
  out << BASE_MODEL_IMPORTS;
  out << MCELL_PATH_SETUP;
  out << "\n";
  out << MCELL_IMPORT;

  // TODO: version check, warning

  out << get_import(SUBSYSTEM);
  out << get_import(INSTANTIATION);
  out << get_import(OBSERVABLES);
  out << get_import(SIMULATION_STATE);
  out << "\n";

  // create model object
  gen_ctor_call(out, MODEL, NAME_CLASS_MODEL, false);
  out << "\n";

  // config, notifications, warnings
  gen_assign(out, MODEL, NAME_CONFIG, model->config.export_to_python(out, ctx));
  out << "\n";
  gen_assign(out, MODEL, NAME_NOTIFICATIONS, model->notifications.export_to_python(out, ctx));
  out << "\n";
  gen_assign(out, MODEL, NAME_WARNINGS, model->warnings.export_to_python(out, ctx));
  out << "\n";

  // subsystem
  string subsystem_prefix = S(SUBSYSTEM) + "." + subsystem_name + ".";
  gen_assign(out, MODEL, NAME_SPECIES, subsystem_prefix + NAME_SPECIES);
  gen_assign(out, MODEL, NAME_REACTION_RULES, subsystem_prefix + NAME_REACTION_RULES);
  gen_assign(out, MODEL, NAME_SURFACE_CLASSES, subsystem_prefix + NAME_SURFACE_CLASSES);
  gen_assign(out, MODEL, NAME_ELEMENTARY_MOLECULE_TYPES, subsystem_prefix + NAME_ELEMENTARY_MOLECULE_TYPES);
  out << "\n";

  // instantiation
  // TODO: do not generate past releases?
  string instantiation_prefix = S(INSTANTIATION) + "." + instantiation_name + ".";
  gen_assign(out, MODEL, NAME_RELEASE_SITES, instantiation_prefix + NAME_RELEASE_SITES);
  gen_assign(out, MODEL, NAME_GEOMETRY_OBJECTS, instantiation_prefix + NAME_GEOMETRY_OBJECTS);
  out << "\n";

  // observables
  string observables_prefix = S(OBSERVABLES) + "." + observables_name + ".";
  gen_assign(out, MODEL, NAME_VIZ_OUTPUTS, observables_prefix + NAME_VIZ_OUTPUTS);
  gen_assign(out, MODEL, NAME_COUNTS, observables_prefix + NAME_COUNTS);
  out << "\n";

  // checkpoint-specific config

  // - time step (explicit)?
  vector<string> config_vars = { NAME_INITIAL_ITERATION, NAME_INITIAL_TIME, NAME_RNG_STATE };
  for (string& var: config_vars) {
    auto it = config_variable_names.find(var);
    release_assert(it != config_variable_names.end());
    gen_assign(out, MODEL, NAME_CONFIG, it->first, S(SIMULATION_STATE) + "." + it->second);
  }

  // - append to observables


  // resume simulation

  // initialize
  // run
  // end

  out.close();
  return MODEL;
}


} // namespace API
} // namespace MCell

