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

#include "model.h"

#include <string>

#include "world.h"
#include "mcell4_converter.h"
#include "api_utils.h"

#include "species.h"
#include "reaction_rule.h"
#include "release_site.h"
#include "geometry_object.h"
#include "viz_output.h"
#include "count.h"

using namespace std;

namespace MCell {
namespace API {


Model::~Model() {
  delete world;
}


void Model::add_subsystem(std::shared_ptr<Subsystem> subsystem) {
  append_vec_to_vec(species, subsystem->species);
  append_vec_to_vec(reaction_rules, subsystem->reaction_rules);
  append_vec_to_vec(surface_classes, subsystem->surface_classes);
}


void Model::add_instantiation_data(std::shared_ptr<InstantiationData> instantiation_data) {
  append_vec_to_vec(release_sites, instantiation_data->release_sites);
  append_vec_to_vec(geometry_objects, instantiation_data->geometry_objects);
}


void Model::add_observables(std::shared_ptr<Observables> observables) {
  append_vec_to_vec(viz_outputs, observables->viz_outputs, true);
  append_vec_to_vec(counts, observables->counts, true);
}


void Model::initialize() {
  if (world != nullptr) {
    throw RuntimeError("Model.initialize() can be called only once");
  }

  world = new World();

  // semantic checks are done during conversion
  MCell4Converter converter;

  converter.convert(this, world);

  // set that all used objects were initialized
  vec_set_initialized(species);
  vec_set_initialized(reaction_rules);
  vec_set_initialized(release_sites);
  vec_set_initialized(geometry_objects);
  vec_set_initialized(viz_outputs);
  vec_set_initialized(counts);

  world->init_simulation();
}


void Model::run_iterations(const long iterations) {
  if (world == nullptr) {
    throw RuntimeError("Model was not initialized, call Model.initialize() first");
  }
  uint output_frequency = World::determine_output_frequency(iterations);
  world->run_n_iterations(iterations, output_frequency, false);
}


void Model::end_simulation(const bool print_final_report) {
  // the first argument specifies that the last mol/rxn count and viz_output events will be run
  world->end_simulation(true, print_final_report);
}


std::string Model::to_str(const std::string ind) const {
  std::stringstream ss;
  ss << "Model" << ": " <<
      "\n" << ind + "  " <<
      //"config=" << "(" << ((config != nullptr) ? config->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      //"warnings=" << "(" << ((warnings != nullptr) ? warnings->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      //"notifications=" << "(" << ((notifications != nullptr) ? notifications->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "reaction_rules=" << vec_ptr_to_str(reaction_rules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "species=" << vec_ptr_to_str(species, ind + "  ") << ", " << "\n" << ind + "  " <<
      "release_sites=" << vec_ptr_to_str(release_sites, ind + "  ") << ", " << "\n" << ind + "  " <<
      "geometry_objects=" << vec_ptr_to_str(geometry_objects, ind + "  ");
  return ss.str();
}


void Model::dump_internal_state() {
  world->dump();
}

void Model::export_data_model(const std::string& file) {
  if (is_set(file)) {
    world->export_data_model(file);
  }
  else {
    // use the first viz_output
    if (viz_outputs.empty()) {
      throw ValueError(
          S("Method ") + NAME_EXPORT_DATA_MODEL + " of " + NAME_CLASS_MODEL + " requires a file argument when there are no instances of " +
          NAME_CLASS_VIZ_OUTPUT + " present."
      );
    }

    if (!is_set(viz_outputs[0]->filename_prefix)) {
      throw ValueError(
          S("Method ") + NAME_EXPORT_DATA_MODEL + " of " + NAME_CLASS_MODEL + ": the first VizOutput instance does not have its " +
          NAME_FILENAME_PREFIX + " set."
      );
    }

    world->export_data_model_to_dir(viz_outputs[0]->filename_prefix);
  }

}


std::vector<int> Model::get_molecule_ids(std::shared_ptr<Species> species) {
  std::vector<int> res;
  assert(false && "TODO");
  return res;
}


std::shared_ptr<Molecule> Model::get_molecule(const int id) {
  std::shared_ptr<Molecule> res;
  assert(false && "TODO");
  return res;
}


void Model::dump() const {
  cout << to_str();
}

}
}

