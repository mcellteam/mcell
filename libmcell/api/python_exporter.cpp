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
#include "api/species.h"
#include "api/geometry_object.h"
#include "api/chkpt_vol_mol.h"
#include "api/chkpt_surf_mol.h"
#include "api/rng_state.h"
#include "api/checkpoint_signals.h"

#include "src/util.h"
#include "src4/world.h"
#include "src4/partition.h"
#include "src4/molecule.h"
#include "src4/custom_function_call_event.h"
#include "src4/mol_or_rxn_count_event.h"

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

  string subsystem_name = export_subsystem(ctx);

  string geometry_objects_name = export_geometry(ctx);

  string instantiation_name = export_instantiation(ctx, geometry_objects_name);

  // need to append - some config flag?
  string observables_name = export_observables(ctx);

  // molecules
  // - volume
  // - surface
  // rng state,
  // starting, ending iterations, ...
  // all generated variables are prefixed with state_
  std::map<std::string, std::string> config_variable_names;
  export_simulation_state(ctx, config_variable_names);

  // model
  // - config
  // - warnings
  // - notifications
  export_model(ctx, subsystem_name, instantiation_name, observables_name, config_variable_names);
}


std::string PythonExporter::export_subsystem(PythonExportContext& ctx) {
  ofstream out;
  open_and_check_file(SUBSYSTEM, out);
  out << IMPORT_MCELL_AS_M;

  string res_name = model->Subsystem::export_to_python(out, ctx);

  out.close();
  return res_name;
}


std::string PythonExporter::export_geometry(PythonExportContext& ctx) {
  ofstream out;
  open_and_check_file(GEOMETRY, out);
  out << IMPORT_MCELL_AS_M;
  out << get_import_star(SUBSYSTEM);

  string res_name = model->export_vec_geometry_objects(out, ctx, "");

  out.close();
  return res_name;
}


std::string PythonExporter::export_instantiation(PythonExportContext& ctx, const std::string& geometry_objects_name) {
  // prints out everything, even past releases
  // for checkpointing, we always need to fully finish the current iteration and then start the new one
  std::ofstream out;
  open_and_check_file(INSTANTIATION, out);
  out << IMPORT_MCELL_AS_M;
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


std::string PythonExporter::export_observables(PythonExportContext& ctx) {
  ofstream out;
  open_and_check_file(OBSERVABLES, out);
  out << IMPORT_MCELL_AS_M;
  out << get_import_star(SUBSYSTEM);
  out << get_import_star(GEOMETRY);
  out << get_import_star(INSTANTIATION);
  out << "\n";

  // we need to set the initial_reactions_count of reaction counts here
  for (auto& count: model->Observables::counts) {
    if (is_set(count->count_expression)) {
      assert(count->count_event != nullptr);
      assert(count->count_event->mol_rxn_count_items.size() == 1);

      for (const auto& count_term: count->count_event->mol_rxn_count_items[0].terms) {
        if (count_term.is_rxn_count()) {
          throw RuntimeError(
            "Checkpointing of Count objects that use expressions containing reaction counts is not supported yet.");
        }
      }
    }
    else {
      if (is_set(count->reaction_rule)) {
        // set value that will be exported
        count->initial_reactions_count_export_override = count->get_current_value();
      }
    }
  }

  string res_name = model->Observables::export_to_python(out, ctx);

  out.close();
  return res_name;
}


// state_variable_names are indexed by Config class attribute names
void PythonExporter::export_simulation_state(
    PythonExportContext& ctx,
    std::map<std::string, std::string>& config_variable_names
) {
  ofstream out;
  open_and_check_file(SIMULATION_STATE, out);
  out << IMPORT_MCELL_AS_M;
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
  export_molecules(out, ctx);

  // rng state
  RngState rng_state = RngState(world->rng);
  config_variable_names[NAME_INITIAL_RNG_STATE] = rng_state.export_to_python(out, ctx);
}


void PythonExporter::export_molecules(std::ostream& out, PythonExportContext& ctx) {

  // first export used species
  stringstream species_out;
  species_out << "# species used by molecules but not defined in subsystem\n";
  species_out << NAME_SPECIES << " = m.Vector" << NAME_CLASS_SPECIES << "([ ";
  int num_exported_species = 0;

  // prepare species map
  IdSpeciesMap id_species_map;
  for (const BNG::Species* species: world->get_all_species().get_species_vector()) {
    assert(id_species_map.count(species->id) == 0);

    if (species->get_num_instantiations() > 0) {
      std::shared_ptr<Species> subsystem_species = model->get_species_with_id(species->id);
      if (is_set(subsystem_species)) {
        // use existing object
        id_species_map[species->id] = subsystem_species;
      }
      else if (!species->is_reactive_surface()){
        // - reactive surfaces/surface classes must be defined in subsystem
        // - create a new object and export it directly so that all the newly used species are
        //   in the beginning of the simulation_state module file
        shared_ptr<API::Species> new_species = make_shared<API::Species>(species->name);
        std::string name = new_species->export_to_python(out, ctx);
        id_species_map[species->id] = new_species;

        // add it to the list of species to be used
        if (num_exported_species != 0 && num_exported_species % 8 == 0) {
          species_out << "\n  ";
        }
        species_out << name << ", ";
      }
    }
  }
  species_out << "\n])\n\n";
  out << species_out.str();

  // prepare geometry objects map
  IdGeometryObjectMap id_geometry_object_map;
  for (const auto& obj: model->geometry_objects) {
    assert(obj->geometry_object_id != GEOMETRY_OBJECT_ID_INVALID);
    id_geometry_object_map[obj->geometry_object_id] = obj;
  }

  // for each partition
  stringstream dummy_out;
  out << NAME_CHECKPOINTED_MOLECULES << " = [\n";

  for (const MCell::Partition& p: world->get_partitions()) {
    // for each molecule
    for (const MCell::Molecule& m: p.get_molecules()) {
      if (m.is_defunct()) {
        continue;
      }
      assert(id_species_map.count(m.species_id) != 0);

      // vol
      PythonExportContext dummy_ctx; // we do not want to use the pointer comparison checks in export
      if (m.is_vol()) {
        ChkptVolMol vm = ChkptVolMol(
            m, id_species_map, world->config.time_unit, world->config.length_unit);
        out << IND4 << vm.export_to_python(dummy_out, ctx) << ",\n";
      }
      else {
        assert(m.is_surf());
        ChkptSurfMol sm = ChkptSurfMol(
            m, id_species_map, world->config.time_unit, world->config.length_unit,
            p, id_geometry_object_map);
        out << IND4 << sm.export_to_python(dummy_out, ctx) << ",\n";
      }
    }
  }
  assert(dummy_out.str() == "" && "No other code must be exported.");

  out << "]\n\n";
}


std::string PythonExporter::export_model(
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
  out << IMPORT_SYS_OS;
  out << "\n";
  out << MCELL_PATH_SETUP;
  out << "\n";
  out << IMPORT_MCELL_AS_M;

  // TODO: version check, warning
  out << make_section_comment("import model and saved simulation state");
  out << get_import(SUBSYSTEM);
  out << get_import(INSTANTIATION);
  out << get_import(OBSERVABLES);
  out << get_import(SIMULATION_STATE);
  out << "\n";

  // create model object
  out << make_section_comment("model setup");

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

  out << make_section_comment("saved simulation state and checkpoint config");
  // - time step (explicit)?
  vector<string> config_vars = { NAME_INITIAL_ITERATION, NAME_INITIAL_TIME, NAME_INITIAL_RNG_STATE };
  for (string& var: config_vars) {
    auto it = config_variable_names.find(var);
    release_assert(it != config_variable_names.end());
    gen_assign(out, MODEL, NAME_CONFIG, it->first, S(SIMULATION_STATE) + "." + it->second);
  }
  gen_assign(out, MODEL, NAME_CONFIG, NAME_APPEND_TO_COUNT_OUTPUT_DATA, true);

  out << "# internal type VectorSpecies does not provide operator += yet\n";
  out << "for s in " << SIMULATION_STATE << "." << NAME_SPECIES << ":";
  out << IND4 << MODEL << "." << NAME_SPECIES << ".append(s)\n";
  gen_assign(out, MODEL, NAME_CHECKPOINTED_MOLECULES, S(SIMULATION_STATE) + "." + NAME_CHECKPOINTED_MOLECULES);
  out << "\n";

  out << make_section_comment("resume simulation");
  gen_method_call(out, MODEL, NAME_INITIALIZE);
  out << "\n";

  export_checkpoint_iterations(out);

  gen_method_call(
      out, MODEL, NAME_RUN_ITERATIONS,
      S(MODEL) + "." + NAME_CONFIG + "." + NAME_TOTAL_ITERATIONS + " - " + S(SIMULATION_STATE) + "." + NAME_INITIAL_ITERATION);
  gen_method_call(out, MODEL, NAME_END_SIMULATION);
  out.close();
  return MODEL;
}


void PythonExporter::export_checkpoint_iterations(std::ostream& out) {

  vector<BaseEvent*> checkpoint_events;
  world->scheduler.get_all_events_with_type_index(
      EVENT_TYPE_INDEX_CALL_START_ITERATION_CHECKPOINT, checkpoint_events);

  for (auto be: checkpoint_events) {
    auto ce = dynamic_cast<CustomFunctionCallEvent<CheckpointSaveEventContext>*>(be);

    bool arg2_nondefault = !ce->return_from_run_n_iterations;
    bool arg3_nondefault = !ce->function_arg.append_it_to_dir; // with custom dir, this is false

    out << MODEL << "." << NAME_SCHEDULE_CHECKPOINT << "(\n";
    gen_param(out, NAME_ITERATION, ce->event_time, arg2_nondefault);

    if (arg2_nondefault) {
      gen_param(out, NAME_ITERATION, !ce->return_from_run_n_iterations, arg3_nondefault);
    }

    if (arg2_nondefault) {
      gen_param(out, NAME_ITERATION, ce->function_arg.dir_prefix, true);
    }
    out << ")\n";
  }
}

} // namespace API
} // namespace MCell

