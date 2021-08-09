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

#include "model.h"

#include <string>

#include "api/mcell4_converter.h"
#include "api/python_exporter.h"
#include "api/api_utils.h"
#include "api/molecule.h"

#include "api/species.h"
#include "api/reaction_rule.h"
#include "api/release_site.h"
#include "api/geometry_object.h"
#include "api/viz_output.h"
#include "api/count.h"
#include "api/wall.h"
#include "api/wall_wall_hit_info.h"
#include "api/checkpoint_signals.h"

#include "world.h"
#include "scheduler.h"
#include "diffuse_react_event.h"
#include "release_event.h"
#include "rxn_utils.inl"
#include "molecule.h"
#include "viz_output_event.h"
#include "custom_function_call_event.h"
#include "bngl_exporter.h"

#include "bng/rxn_class.h"


using namespace std;

namespace MCell {
namespace API {

class RngState;

Model::~Model() {
  unset_checkpoint_signals(this); // may be safely called multiple times
  delete world;
}


void Model::add_subsystem(std::shared_ptr<Subsystem> subsystem) {
  error_if_initialized(NAME_CLASS_SUBSYSTEM);

  append_vec_to_vec_canonical_name(elementary_molecule_types, subsystem->elementary_molecule_types);
  append_vec_to_vec_canonical_name(species, subsystem->species);
  append_vec_to_vec(surface_classes, subsystem->surface_classes);
  append_vec_to_vec_canonical_name(reaction_rules, subsystem->reaction_rules);
}


void Model::add_instantiation(std::shared_ptr<Instantiation> instantiation) {
  error_if_initialized(NAME_CLASS_INSTANTIATION);

  append_vec_to_vec(release_sites, instantiation->release_sites);
  append_vec_to_vec(geometry_objects, instantiation->geometry_objects);
}


void Model::add_observables(std::shared_ptr<Observables> observables) {
  error_if_initialized(NAME_CLASS_OBSERVABLES);

  append_vec_to_vec(viz_outputs, observables->viz_outputs, true);
  append_vec_to_vec(counts, observables->counts, true);
}


void Model::initialize(const bool print_copyright) {
  if (world != nullptr) {
    throw RuntimeError("Model.initialize() can be called only once");
  }

  // semantic checks of the following classes must be called manually because
  // their contents may have changed and check_semantics is automatically called only from their ctor
  config.check_semantics();
  warnings.check_semantics();
  notifications.check_semantics();

  // first add species superclasses
  std::vector<std::shared_ptr<Species>> superspecies = { AllMolecules, AllVolumeMolecules, AllSurfaceMolecules };
  species.insert(species.begin(), superspecies.begin(), superspecies.end());

  // Species objects might have created their own ElementaryMoleculeType
  // objects, we must unify it first (required for python export)
  unify_and_register_elementary_molecule_types();

  world = new World(callbacks);

  // semantic checks are done during conversion
  MCell4Converter converter(this, world);

  converter.convert_before_init();

  // set that all used objects were initialized
  vec_set_initialized(species);
  vec_set_initialized(reaction_rules);
  vec_set_initialized(release_sites);
  vec_set_initialized(geometry_objects);
  vec_set_initialized(viz_outputs);
  vec_set_initialized(counts);

  world->init_simulation(world->config.get_simulation_start_time());

  // convert also rng state checkpointed molecules,
  // must be done after world initialization
  converter.convert_after_init();

  initialize_introspection(this);

  // set action for SIGUSR1 and SIGUSR2
  set_checkpoint_signals(this);

  initialized = true;

  if (print_copyright) {
    cout <<
      "Copyright (C) 2006-2021 by\n"
      "  The National Center for Multiscale Modeling of Biological Systems,\n"
      "  The Salk Institute for Biological Studies, and\n"
      "  Pittsburgh Supercomputing Center, Carnegie Mellon University,\n"
      "\n"
      "**********************************************************************\n"
      "MCell development is supported by the NIGMS-funded (P41GM103712)\n"
      "National Center for Multiscale Modeling of Biological Systems (MMBioS).\n"
      "Please acknowledge MCell in your publications.\n"
      "**********************************************************************\n\n";
  }
}


uint64_t Model::run_iterations(const double iterations) {
  // release GIL before calling into long-running C++ code, 
  // necessary for other Python threads to be allowed to run (e.g. a timer)
  py::gil_scoped_release release;

  if (world == nullptr) {
    throw RuntimeError("Model was not initialized, call Model.initialize() first");
  }
  return world->run_n_iterations(iterations, false);
}


void Model::end_simulation(const bool print_final_report) {
  // the first argument specifies that the last mol/rxn count and viz_output events will be run
  // cannot call callbacks anymore (no diffusion event is executed)
  world->end_simulation(print_final_report);

  // unset action for SIGUSR1 and SIGUSR2
  unset_checkpoint_signals(this);
}


void Model::dump_internal_state(const bool with_geometry) {
  world->dump(with_geometry);
}


void Model::export_geometry(const std::string& output_files_prefix) {
  if (world == nullptr) {
    throw RuntimeError(S("Model must be initialized before a call to ") + NAME_EXPORT_DATA_MODEL + ".");
  }

  if (is_set(output_files_prefix)) {
    world->export_geometry_to_obj(output_files_prefix);
  }
  else {
    stringstream prefix;
    prefix <<
        get_first_viz_output_files_prefix(NAME_EXPORT_GEOMETRY) << ".geom." <<
        VizOutputEvent::iterations_to_string(world->stats.get_current_iteration(), config.total_iterations);

    world->export_geometry_to_obj(prefix.str());
  }
}


void Model::export_data_model_viz_or_full(
    const std::string& file,
    const bool only_for_visualization,
    const char* method_name) {

  if (world == nullptr) {
    throw RuntimeError(S("Model must be initialized before a call to ") + NAME_EXPORT_DATA_MODEL + ".");
  }

  if (is_set(file)) {
    world->export_data_model(file, only_for_visualization);
  }
  else {
    world->export_data_model_to_dir(get_first_viz_output_files_prefix(method_name), only_for_visualization);
  }
}


void Model::release_molecules(std::shared_ptr<ReleaseSite> release_site) {
  if (!initialized) {
    throw RuntimeError(S("Model must be initialized before calling ") + NAME_RELEASE_MOLECULES + ".");
  }

  // check that time is now or in the future
  double iteration_start_time = world->stats.get_current_iteration() * world->config.time_unit;
  if (release_site->release_time < iteration_start_time) {
    throw ValueError("Cannot release molecules for time " + to_string(release_site->release_time) +
        " before the start time of the current iteration " + to_string(iteration_start_time) + ".");
  }

  if (is_set(release_site->release_pattern)) {
    throw ValueError(S("Cannot release molecules with a release pattern, method ") + NAME_RELEASE_MOLECULES +
        " may be used only for immediate releases.");
  }

  // convert to a ReleaseEvent
  MCell4Converter converter(this, world);
  MCell::ReleaseEvent* rel_event = converter.convert_single_release_event(release_site);
  assert(rel_event != nullptr);

  // TODO: we must improve handling of cases when the release is in the future
  rel_event->event_time = release_site->release_time / world->config.time_unit;

  // figure out whether the DiffuseAndReactEvent is running
  MCell::BaseEvent* current_event = world->scheduler.get_event_being_executed();

  MCell::DiffuseReactEvent* diffuse_event = nullptr;
  if (current_event != nullptr && current_event->type_index == EVENT_TYPE_INDEX_DIFFUSE_REACT) {
    diffuse_event = dynamic_cast<MCell::DiffuseReactEvent*>(current_event);
    assert(diffuse_event != nullptr);
  }
  
  // and execute the release
  rel_event->release_immediatelly(diffuse_event);
  delete rel_event;
}


std::vector<int> Model::run_reaction(
    std::shared_ptr<ReactionRule> reaction_rule,
    const std::vector<int> reactant_ids,
    const double time) {

  // need to release lock because we might be calling callbacks that are locking back
  py::gil_scoped_release release;

  if (!initialized) {
    throw RuntimeError(S("Model must be initialized before calling ") + NAME_RUN_REACTION + ".");
  }

  if (reaction_rule->reactants.size() != reactant_ids.size()) {
    throw RuntimeError("Reaction expects " + to_string(reaction_rule->reactants.size()) +
        " reactants but " + to_string(reactant_ids.size()) + " reactants were provided.");
  }

  if (reaction_rule->fwd_rxn_rule_id == BNG::RXN_RULE_ID_INVALID) {
    throw RuntimeError("Reaction rule is not present in model and was not initialized.");
  }

  if (reaction_rule->rev_rxn_rule_id != BNG::RXN_RULE_ID_INVALID) {
    throw RuntimeError(S("Method ") + NAME_RUN_REACTION + " can be used only with irreversible reactions.");
  }

  const BNG::RxnRule* rxn = world->get_all_rxns().get(reaction_rule->fwd_rxn_rule_id);

  if (!rxn->is_unimol()) {
    throw RuntimeError(S("Method ") + NAME_RUN_REACTION + " currently supports only unimolecular reactions.");
  }

  MCell::Partition& p = world->get_partition(PARTITION_ID_INITIAL);
  molecule_id_t id1 = reactant_ids[0];
  if (!p.does_molecule_exist(id1)) {
    throw RuntimeError("Molecule with id " + to_string(id1) + " does not exist.");
  }
  MCell::Molecule& m1 = p.get_m(id1);
  if (m1.is_defunct()) {
    throw RuntimeError("Molecule with id " + to_string(id1) + " was removed.");
  }

  std::vector<int> res;
  if (rxn->is_unimol()) {
    // check if the requested rxn is applicable for our reactant
    // also determine pathway index
    BNG::RxnClass* rxn_class = world->bng_engine.get_all_rxns().get_unimol_rxn_class(m1.species_id);
    rxn_class->init_rxn_pathways_and_rates(); // initialize if needed
    BNG::rxn_class_pathway_index_t index = 0;
    while (index < (BNG::rxn_class_pathway_index_t)rxn_class->get_num_pathways() &&
        rxn_class->get_rxn_for_pathway(index)->id != rxn->id) {
      index++;
    }

    if (index >= (BNG::rxn_class_pathway_index_t)rxn_class->get_num_pathways()) {
      const BNG::Species& species = world->get_all_species().get(m1.species_id);
      throw RuntimeError("Reaction rule " + reaction_rule->to_bngl_str() +
          " cannot be applied on molecule with species " + species.name);
    }

    // if we are in a callback, we are probably in a diffuse react event
    BaseEvent* current_event = world->scheduler.get_event_being_executed();
    DiffuseReactEvent* diffuse_react_event;
    bool using_temporary_event;
    if (current_event != nullptr && current_event->type_index == EVENT_TYPE_INDEX_DIFFUSE_REACT) {
      diffuse_react_event = dynamic_cast<DiffuseReactEvent*>(current_event);
      using_temporary_event = false;
    }
    else {
      // otherwise we will instantiate a new event
      diffuse_react_event = new DiffuseReactEvent(world);
      using_temporary_event = true;
    }

    MoleculeIdsVector product_ids;
    diffuse_react_event->outcome_unimolecular(p, m1, time / world->config.time_unit, rxn_class, index, &product_ids);

    if (using_temporary_event) {
      delete diffuse_react_event;
    }

    res.insert(res.begin(), product_ids.begin(), product_ids.end());
  }
  else {
    // TODO
    release_assert(false);
  }

  return res;
}


void Model::add_vertex_move(
    std::shared_ptr<GeometryObject> object, const int vertex_index, const std::vector<double> displacement
) {
  // - currently, it is not expected that the user will have access to the scheduled vertex moves
  // - later we can use the object id to determine the partition

  object->check_is_initialized();

  release_assert(
      object->first_vertex_index + vertex_index <
      world->get_partition(PARTITION_ID_INITIAL).get_geometry_vertex_count()
  );

  if (displacement.size() != 3) {
    throw ValueError(S("Argument ") + NAME_DISPLACEMENT + " must be a list of exactly three floats.");
  }
  Vec3 disp_vec = Vec3(displacement);

  vertex_moves.push_back(
      VertexMoveInfo(
          PARTITION_ID_INITIAL,
          object->geometry_object_id,
          object->get_partition_vertex_index(vertex_index),
          disp_vec * Vec3(world->config.rcp_length_unit) // convert to internal units
      )
  );
}


std::vector<std::shared_ptr<WallWallHitInfo>> Model::apply_vertex_moves(
    const bool collect_wall_wall_hits,
    const bool randomize_order) {

  // run the actual vertex update
  std::set<GeometryObjectWallUnorderedPair> colliding_walls;
  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  // may update vertex_move displacements due to collisions
  p.apply_vertex_moves(randomize_order, vertex_moves, colliding_walls);

  std::vector<std::shared_ptr<WallWallHitInfo>> res;

  // convert information on processed hits
  if (collect_wall_wall_hits) {
    for (const auto& wall_pair: colliding_walls) {
      auto obj1 = get_geometry_object_with_id(wall_pair.geometry_object_id1);
      int obj_wall_index1 = obj1->get_object_wall_index(wall_pair.wall_index1);
      auto obj2 = get_geometry_object_with_id(wall_pair.geometry_object_id2);
      int obj_wall_index2 = obj2->get_object_wall_index(wall_pair.wall_index2);

      auto pair = make_shared<WallWallHitInfo>(DefaultCtorArgType());
      pair->wall1 = get_wall(obj1, obj_wall_index1);
      pair->wall2 = get_wall(obj2, obj_wall_index2);

      res.push_back(pair);
    }
  }

  // also update the API copy of the geometry
  for (const VertexMoveInfo& move_info: vertex_moves) {
    // get the new vertex coordinates
    assert(move_info.partition_id == PARTITION_ID_INITIAL);
    const Vec3& pos = p.get_geometry_vertex(move_info.vertex_index);

    // API object whose vertex we need to update
    std::shared_ptr<API::GeometryObject> obj =
        get_geometry_object_with_id(move_info.geometry_object_id);

    // API object's vertex index
    int obj_vertex_index = obj->get_object_vertex_index(move_info.vertex_index);

    Vec3 pos_um = pos * Vec3(world->config.length_unit);
    auto& pos_to_update = obj->vertex_list[obj_vertex_index];
    pos_to_update[0] = pos_um.x;
    pos_to_update[1] = pos_um.y;
    pos_to_update[2] = pos_um.z;
  }

  vertex_moves.clear();

  return res;
}


void Model::register_mol_wall_hit_callback(
    const std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)> function,
    py::object context,
    std::shared_ptr<GeometryObject> object,
    std::shared_ptr<Species> species
) {
  if (!initialized) {
    throw RuntimeError("Model must be initialized before registering callbacks.");
  }

  geometry_object_id_t geometry_object_id = GEOMETRY_OBJECT_ID_INVALID;
  if (is_set(object)) {
    if (object->geometry_object_id == GEOMETRY_OBJECT_ID_INVALID) {
      throw RuntimeError("Geometry object " + object->name + " is not present in model.");
    }
    geometry_object_id = object->geometry_object_id;
  }

  species_id_t species_id = SPECIES_ID_INVALID;
  if (is_set(species)) {
    if (species->species_id == SPECIES_ID_INVALID) {
      throw RuntimeError("Species object " + object->name + " is not present in model.");
    }
    species_id = species->species_id;
  }

  callbacks.register_mol_wall_hit_callback(function, context, geometry_object_id, species_id);
}


void Model::register_reaction_callback(
    const std::function<void(std::shared_ptr<ReactionInfo>, py::object)> function,
    py::object context,
    std::shared_ptr<ReactionRule> reaction_rule
) {
  if (!initialized) {
    throw RuntimeError("Model must be initialized before registering callbacks.");
  }

  if (reaction_rule->is_reversible()) {
    throw RuntimeError(S("Reaction callback cannot be registered for reversible reactions. ") +
        "Split the reaction rule's forward and reverese direction into separate " +
        NAME_REACTION_RULE + " objects.");
  }

  BNG::rxn_rule_id_t rxn_id = reaction_rule->fwd_rxn_rule_id;
  if (rxn_id == BNG::RXN_RULE_ID_INVALID) {
    throw RuntimeError(S(NAME_REACTION_RULE) + reaction_rule->name + " with its BNGL representation " +
        reaction_rule->to_bngl_str() + " is not present in model.");
  }

  callbacks.register_rxn_callback(function, context, rxn_id);
}


void Model::load_bngl(
    const std::string& file_name,
    const std::string& observables_path_or_file,
    std::shared_ptr<Region> default_release_region,
    const std::map<std::string, double>& parameter_overrides,
    const CountOutputFormat observables_output_format) {

  BNG::BNGData bng_data;

  int num_errors = BNG::parse_bngl_file(file_name, bng_data, parameter_overrides);
  if (num_errors != 0) {
    throw RuntimeError("Could not parse BNGL file " + file_name + ".");
  }

  // now convert everything we parsed into the API classes so that the user can
  // inspect or manipulate it if needed
  convert_bng_data_to_subsystem_data(bng_data);

  // needs subsystem data created in the last step
  convert_bng_data_to_instantiation(bng_data, default_release_region);

  convert_bng_data_to_observables_data(
      bng_data,
      get_count_output_format(observables_output_format, observables_path_or_file),
      observables_path_or_file);
}


void Model::export_to_bngl(const std::string& file_name,
    const BNGSimulationMethod simulation_method) {
  if (!initialized) {
    throw RuntimeError("Model must be initialized for BNGL export.");
  }

  if (simulation_method != BNGSimulationMethod::NONE &&
      simulation_method != BNGSimulationMethod::ODE &&
      simulation_method != BNGSimulationMethod::SSA &&
      simulation_method != BNGSimulationMethod::PLA &&
      simulation_method != BNGSimulationMethod::NF) {
    throw ValueError(S("Invalid ") + NAME_SIMULATION_METHOD + " argument value.");
  }

  MCell::BNGLExporter bngl_exporter;
  string err_msg = bngl_exporter.export_to_bngl(world, file_name, simulation_method);

  if (err_msg != "") {
    throw RuntimeError("BNGL export failed: " + err_msg);
  }
}


void Model::save_checkpoint(const std::string& custom_dir) {
  if (!initialized) {
    throw RuntimeError(S("Model must be initialized for ") + NAME_SAVE_CHECKPOINT + ".");
  }

  // prepare output directory name
  string dir;
  if (is_set(custom_dir)) {
    dir = custom_dir;
  }
  else {
    // TODO: move the VizOutputEvent::iterations_to_string to api_utils
    dir =
        world->config.get_default_checkpoint_dir_prefix() +
        VizOutputEvent::iterations_to_string(world->stats.get_current_iteration(), config.total_iterations) +
        BNG::PATH_SEPARATOR;
  }

  PythonExporter exporter(this);
  exporter.save_checkpoint(dir);
}


// can be called asynchronously at any point of simulation, only single-threaded
// cannot be called from signal handlers
void Model::schedule_checkpoint(
    const uint64_t iteration,
    const bool continue_simulation,
    const std::string& custom_dir) {

  // NOTE: accessing current iteration value asynchronously,
  // it is a 64-bit integer so write is done atomically with a single instruction
  uint64_t current_it = world->stats.get_current_iteration();
  if (iteration != 0 && iteration < current_it) {
    throw RuntimeError(S("Method ") + NAME_SCHEDULE_CHECKPOINT + " - cannot schedule checkpoint to the past.");
  }

  // same as above - 'initialized' uses a simple type, therefore writes are atomic
  if (!initialized) {
    throw RuntimeError(S("Method ") + NAME_SCHEDULE_CHECKPOINT + " cannot be called before model initialization.");
  }

  // prepare context for callback
  CheckpointSaveEventContext ctx;
  ctx.model = this;

  if (is_set(custom_dir)) {
    ctx.dir_prefix = custom_dir;
    ctx.append_it_to_dir = false;
  }
  else {
    ctx.dir_prefix = world->config.get_default_checkpoint_dir_prefix();
    ctx.append_it_to_dir = true;
  }

  world->schedule_checkpoint_event(iteration, continue_simulation, ctx);
}


// overrides from derived classes Subsystem, Instantiation, and Observables,
// in .cpp because implementation in .h file would need too many headers to be included
void Model::add_species(std::shared_ptr<Species> s) {
  error_if_initialized(NAME_CLASS_SPECIES);
  Subsystem::add_species(s);
}

void Model::add_reaction_rule(std::shared_ptr<ReactionRule> r) {
  error_if_initialized(NAME_CLASS_REACTION_RULE);
  Subsystem::add_reaction_rule(r);
}

void Model::add_surface_class(std::shared_ptr<SurfaceClass> sc) {
  error_if_initialized(NAME_CLASS_SURFACE_CLASS);
  Subsystem::add_surface_class(sc);
}

void Model::add_release_site(std::shared_ptr<ReleaseSite> s) {
  error_if_initialized(NAME_CLASS_RELEASE_SITE);
  Instantiation::add_release_site(s);
}

void Model::add_geometry_object(std::shared_ptr<GeometryObject> o) {
  error_if_initialized(NAME_CLASS_GEOMETRY_OBJECT);
  Instantiation::add_geometry_object(o);
}

void Model::add_viz_output(std::shared_ptr<VizOutput> viz_output) {
  error_if_initialized(NAME_CLASS_VIZ_OUTPUT);
  Observables::add_viz_output(viz_output);
};

void Model::add_count(std::shared_ptr<Count> count) {
  error_if_initialized(NAME_CLASS_OBSERVABLES);
  Observables::add_count(count);
};


std::shared_ptr<GeometryObject> Model::get_geometry_object_with_id(const geometry_object_id_t id) {

  // map for caching
  auto it = geometry_object_id_to_obj_cache.find(id);
  if (it != geometry_object_id_to_obj_cache.end()) {
    return it->second;
  }

  for (auto o: geometry_objects) {
    if (o->geometry_object_id == id) {
      geometry_object_id_to_obj_cache[o->geometry_object_id] = o;
      return o;
    }
  }
  return std::shared_ptr<GeometryObject>(nullptr);
}


std::shared_ptr<ReactionRule> Model::get_reaction_rule_with_fwd_id(const BNG::rxn_rule_id_t id) {
  // not very efficient, we may need some caching/map later
  for (auto r: reaction_rules) {
    if (r->fwd_rxn_rule_id == id) {
      return r;
    }
  }
  return std::shared_ptr<ReactionRule>(nullptr);
}

void Model::dump() const {
  cout << to_str();
}

}
}

