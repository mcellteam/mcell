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

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_model.h"
#include "api/model.h"
#include "api/base_chkpt_mol.h"
#include "api/complex.h"
#include "api/api_config.h"
#include "api/count.h"
#include "api/elementary_molecule_type.h"
#include "api/geometry_object.h"
#include "api/instantiation.h"
#include "api/mol_wall_hit_info.h"
#include "api/molecule.h"
#include "api/notifications.h"
#include "api/observables.h"
#include "api/reaction_info.h"
#include "api/reaction_rule.h"
#include "api/region.h"
#include "api/release_site.h"
#include "api/species.h"
#include "api/subsystem.h"
#include "api/surface_class.h"
#include "api/viz_output.h"
#include "api/wall.h"
#include "api/wall_wall_hit_info.h"
#include "api/warnings.h"

namespace MCell {
namespace API {

bool GenModel::__eq__(const Model& other) const {
  return
    config == other.config &&
    warnings == other.warnings &&
    notifications == other.notifications &&
    vec_ptr_eq(species, other.species) &&
    vec_ptr_eq(reaction_rules, other.reaction_rules) &&
    vec_ptr_eq(surface_classes, other.surface_classes) &&
    vec_ptr_eq(elementary_molecule_types, other.elementary_molecule_types) &&
    vec_ptr_eq(release_sites, other.release_sites) &&
    vec_ptr_eq(geometry_objects, other.geometry_objects) &&
    vec_ptr_eq(checkpointed_molecules, other.checkpointed_molecules) &&
    vec_ptr_eq(viz_outputs, other.viz_outputs) &&
    vec_ptr_eq(counts, other.counts);
}

bool GenModel::eq_nonarray_attributes(const Model& other, const bool ignore_name) const {
  return
    config == other.config &&
    warnings == other.warnings &&
    notifications == other.notifications &&
    true /*species*/ &&
    true /*reaction_rules*/ &&
    true /*surface_classes*/ &&
    true /*elementary_molecule_types*/ &&
    true /*release_sites*/ &&
    true /*geometry_objects*/ &&
    true /*checkpointed_molecules*/ &&
    true /*viz_outputs*/ &&
    true /*counts*/;
}

std::string GenModel::to_str(const std::string ind) const {
  #if 0 // not generated correctly yet
  std::stringstream ss;
  ss << "Model" << ": " <<
      "\n" << ind + "  " << "config=" << "(" << ((config != nullptr) ? config->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "warnings=" << "(" << ((warnings != nullptr) ? warnings->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "notifications=" << "(" << ((notifications != nullptr) ? notifications->to_str(ind + "  ") : "null" ) << ")" << ", " << "\n" << ind + "  " <<
      "species=" << vec_ptr_to_str(species, ind + "  ") << ", " << "\n" << ind + "  " <<
      "reaction_rules=" << vec_ptr_to_str(reaction_rules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_classes=" << vec_ptr_to_str(surface_classes, ind + "  ") << ", " << "\n" << ind + "  " <<
      "elementary_molecule_types=" << vec_ptr_to_str(elementary_molecule_types, ind + "  ") << ", " << "\n" << ind + "  " <<
      "release_sites=" << vec_ptr_to_str(release_sites, ind + "  ") << ", " << "\n" << ind + "  " <<
      "geometry_objects=" << vec_ptr_to_str(geometry_objects, ind + "  ") << ", " << "\n" << ind + "  " <<
      "checkpointed_molecules=" << vec_ptr_to_str(checkpointed_molecules, ind + "  ") << ", " << "\n" << ind + "  " <<
      "viz_outputs=" << vec_ptr_to_str(viz_outputs, ind + "  ") << ", " << "\n" << ind + "  " <<
      "counts=" << vec_ptr_to_str(counts, ind + "  ");
  return ss.str();
  #else
  return "";
  #endif
}

py::class_<Model> define_pybinding_Model(py::module& m) {
  return py::class_<Model, std::shared_ptr<Model>>(m, "Model")
      .def(
          py::init<
          >()
      )
      .def("__str__", &Model::to_str, py::arg("ind") = std::string(""))
      .def("__eq__", &Model::__eq__, py::arg("other"))
      .def("initialize", &Model::initialize, py::arg("print_copyright") = true, "Initializes model, initialization blocks most of changes to \ncontained components. \n\n- print_copyright\n")
      .def("run_iterations", &Model::run_iterations, py::arg("iterations"), "Runs specified number of iterations. Returns the number of iterations\nexecuted (it might be less than the requested number of iterations when \na checkpoint was scheduled). \n\n- iterations: Number of iterations to run. Value is truncated to an integer.\n\n")
      .def("end_simulation", &Model::end_simulation, py::arg("print_final_report") = true, "Generates the last visualization and reaction output (if they were defined), then\nflushes all buffers and optionally prints simulation report. \nBuffers are also flushed when the Model object is destroyed.   \n\n- print_final_report\n")
      .def("add_subsystem", &Model::add_subsystem, py::arg("subsystem"), "- subsystem\n")
      .def("add_instantiation", &Model::add_instantiation, py::arg("instantiation"), "- instantiation\n")
      .def("add_observables", &Model::add_observables, py::arg("observables"), "- observables\n")
      .def("dump_internal_state", &Model::dump_internal_state, "Prints out the simulation engine's internal state, mainly for debugging.")
      .def("export_data_model", &Model::export_data_model, py::arg("file") = STR_UNSET, "If file is not set, then uses the first VizOutput to determine the target directory \nand creates name using the current iteration. Fails if argument file is not set and there is no VizOutput.\nMust be called after initialization.\nAlways exports the current state, i.e. with the current . \nEvents (ReleaseSites and VizOutputs) with scheduled time other than zero cannot be imported correectly yet.  \n\n- file\n")
      .def("export_viz_data_model", &Model::export_viz_data_model, py::arg("file") = STR_UNSET, "Same as export_data_model, only the created data model will contain only information required for visualization in CellBlender. This makes the loading ofthemodel by CellBlender faster and also allows to avoid potential compatibility issues.\n- file\n")
      .def("release_molecules", &Model::release_molecules, py::arg("release_site"), "Performs immediate release based on the definition of the release site argument.\nThe ReleaseSite.release_time must not be in the past and should be withing the current iteration.\nThe ReleaseEvent must not use a release_pattern because this is an immediate release and it is not \nscheduled into the global scheduler. \n\n- release_site\n")
      .def("run_reaction", &Model::run_reaction, py::arg("reaction_rule"), py::arg("reactant_ids"), py::arg("time"), "Run a single reaction on reactants. Callbacks will be called if they are registered for the given reaction.\nReturns a list of product IDs.\nNote: only unimolecular reactions are currently supported.\n\n- reaction_rule: Reaction rule to run.\n\n- reactant_ids: The number of reactants for a unimolecular reaction must be 1 and for a bimolecular reaction must be 2.\nReactants for a bimolecular reaction do not have to be listed in the same order as in the reaction rule definition. \n\n\n- time: Precise time in seconds when this reaction occurs. Important to know for how long the products\nwill be diffused when they are created in a middle of a time step. \n\n\n")
      .def("add_vertex_move", &Model::add_vertex_move, py::arg("object"), py::arg("vertex_index"), py::arg("displacement"), "Adds a displacement for given object's vertex, only stored until apply_vertex_moves is called\n- object: Object whose vertex will be changed\n\n- vertex_index: Index of vertex in object's vertex list that will be changed\n\n- displacement: Change of vertex coordinates (in um), will be added to the current coordinates of the vertex,\nmust contain exactly three floating point values.\n\n\n")
      .def("apply_vertex_moves", &Model::apply_vertex_moves, py::arg("collect_wall_wall_hits") = false, "Applies all the vertex moves specified with add_vertex_move call.\nWalls of different objects are checked against collisions and move the maximal way so that they do not \noverlap. (the current pllementation is a bit basic and may not work 100% correctly) \nWhen collect_wall_wall_hits is True, a list of wall pairs that collided is returned,\nwhen collect_wall_wall_hits is False, and empty list is returned.\n\n- collect_wall_wall_hits: When set to True, a list of wall pairs that collided is returned,\notherwise an empty list is returned.\n\n\n")
      .def("register_mol_wall_hit_callback", &Model::register_mol_wall_hit_callback, py::arg("function"), py::arg("context"), py::arg("object") = nullptr, py::arg("species") = nullptr, "There can be currently only a single wall hit callback registered.\n- function: Callback function to be called. \nIt must have two arguments MolWallHitInfo and context.\n\n\n- context: Context passed to the callback function, the callback function can store\ninformation to this object. Some context must be always passed, even when \nit is a useless python object. \n\n\n- object: Only hits of this object will be reported, any object hit is reported when not set.\n\n- species: Only hits of molecules of this species will be reported, any species hit is reported when not set.\n\n")
      .def("register_reaction_callback", &Model::register_reaction_callback, py::arg("function"), py::arg("context"), py::arg("reaction_rule"), "Defines a function to be called when a reaction was processed.\nIt is allowed to do state modifications except for removing reacting molecules, \nthey will be removed automatically after return from this callback.  \n\n- function: Callback function to be called. \nIt must have two arguments ReactionInfo and context.\nCalled when it is decided that the reaction will happen.\nAfter return the reaction proceeds as it would without a callback. \n\n\n- context: Context passed to the callback function, the callback function can store\ninformation to this object. Some context must be always passed, even when \nit is a useless python object. \n\n\n- reaction_rule: The callback function will be called whenever is this reaction rule applied.\n\n")
      .def("load_bngl", &Model::load_bngl, py::arg("file_name"), py::arg("observables_files_prefix") = "", py::arg("default_release_region") = nullptr, py::arg("parameter_overrides") = std::map<std::string, float_t>(), "Loads sections: molecule types, reaction rules, seed species, and observables from a BNGL file\nand creates objects in the current model according to it.\nAll elementary molecule types used in the seed species section must be defined in subsystem.\nIf an item in the seed species section does not have its compartment set,\nthe argument default_region must be set and the molecules are released into or onto the \ndefault_region. \n\n- file_name\n- observables_files_prefix: Prefix to be used when creating files with observable values.\n\n- default_release_region\n- parameter_overrides\n")
      .def("export_to_bngl", &Model::export_to_bngl, py::arg("file_name"), "Exports all defined species, reaction rules and applicable observables\nas a BNGL file. \nLimited currrently to exactly one volume compartment and volume reactions.\n\n- file_name: Output file name.\n\n")
      .def("save_checkpoint", &Model::save_checkpoint, py::arg("custom_dir") = STR_UNSET, "Saves current model state as checkpoint. \nThe default directory structure is checkpoints/seed_<SEED>/it_<ITERATION>,\nit can be changed by setting 'custom_dir'.\nIf used during an iteration, schedules an event for the end of the current iteration\nthat saves the checkpoint (effectively calls 'checkpoint_after_iteration(0, False, custom_dir)'.  \n\n- custom_dir: Sets custom directory where the checkpoint will be stored. \nThe default is 'checkpoints/seed_<SEED>/it_<ITERATION>'. \n\n\n")
      .def("schedule_checkpoint", &Model::schedule_checkpoint, py::arg("iteration") = 0, py::arg("continue_simulation") = false, py::arg("custom_dir") = STR_UNSET, "Schedules checkpoint save that will occur when an iteration is started  \nright before any other events scheduled for the given iteration are executed.\nCan be called asynchronously at any time after initialization.\n\n- iteration: Specifies iteration number when the checkpoint save will occur. \nPlease note that iterations are counted from 0.\nTo schedule a checkpoint for the closest time as possible, keep the default value 0,\nthis will schedule checkpoint for the beginning of the iteration with number current iteration + 1.  \nIf calling schedule_checkpoint from a different thread (e.g. by using threading.Timer), \nit is highly recommended to keep the default value 0 or choose some time that will be \nfor sure in the future.\n\n\n- continue_simulation: When false, saving the checkpoint means that we want to terminate the simulation \nright after the save, the currently running function Model.run_iterations\ndoes not simulate any following iterations and execution returns from this function\nto execute the next statement which is usually 'model.end_simulation()'.\nWhen true, the checkpoint is just saved and simulation continues uninterrupted.\n      \n\n\n- custom_dir: Sets custom directory where the checkpoint will be stored. \nThe default is 'checkpoints/seed_<SEED>/it_<ITERATION>'. \n      \n\n")
      .def("add_species", &Model::add_species, py::arg("s"), "- s\n")
      .def("find_species", &Model::find_species, py::arg("name"), "- name\n")
      .def("add_reaction_rule", &Model::add_reaction_rule, py::arg("r"), "- r\n")
      .def("find_reaction_rule", &Model::find_reaction_rule, py::arg("name"), "- name\n")
      .def("add_surface_class", &Model::add_surface_class, py::arg("sc"), "- sc\n")
      .def("find_surface_class", &Model::find_surface_class, py::arg("name"), "- name\n")
      .def("add_elementary_molecule_type", &Model::add_elementary_molecule_type, py::arg("mt"), "- mt\n")
      .def("find_elementary_molecule_type", &Model::find_elementary_molecule_type, py::arg("name"), "- name\n")
      .def("load_bngl_molecule_types_and_reaction_rules", &Model::load_bngl_molecule_types_and_reaction_rules, py::arg("file_name"), py::arg("parameter_overrides") = std::map<std::string, float_t>(), "Parses a BNGL file and only reads molecule types and\nreaction rules sections, e.g. ignores observables. \nParameter values are evaluated and the result value \nis directly used.  \nCompartments names are stored in rxn rules as strings because\ncompartments belong to geometry objects and the subsystem is independent\non specific geometry.\nHowever they must be defined on initialization.\n \n\n- file_name\n- parameter_overrides\n")
      .def("add_release_site", &Model::add_release_site, py::arg("s"), "Adds a reference to the release site s to the list of release sites.\n- s\n")
      .def("find_release_site", &Model::find_release_site, py::arg("name"), "Finds a release site by its name, returns None if no such release site is present.\n- name\n")
      .def("add_geometry_object", &Model::add_geometry_object, py::arg("o"), "Adds a reference to the geometry object o to the list of geometry objects.\n- o\n")
      .def("find_geometry_object", &Model::find_geometry_object, py::arg("name"), "Finds a geometry object by its name, returns None if no such geometry object is present.\n- name\n")
      .def("find_volume_compartment_object", &Model::find_volume_compartment_object, py::arg("name"), "Finds a geometry object by its name, the geometry object must be a BNGL compartment.\nReturns None if no such geometry object is present.\n\n- name\n")
      .def("find_surface_compartment_object", &Model::find_surface_compartment_object, py::arg("name"), "Finds a geometry object that is a BNGL compartment and its surface name is name.\nReturns None if no such geometry object is present.\n\n- name\n")
      .def("load_bngl_seed_species", &Model::load_bngl_seed_species, py::arg("file_name"), py::arg("default_release_region") = nullptr, py::arg("parameter_overrides") = std::map<std::string, float_t>(), "Loads section seed species from a BNGL file and creates release sites according to it.\nAll elementary molecule types used in the seed species section must be already defined in subsystem.\nIf an item in the BNGL seed species section does not have its compartment set,\nthe argument default_region must be set and the molecules are then released into or onto the \ndefault_region. \nDoes not create geometry objects. \nAll compartments used in the loaded BNGL seed species section must exist in the model before \nmodel intialization.\n  \n\n- file_name: Path to the BNGL file.\n\n- default_release_region: Used for seed species that have no compartments specified.\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n")
      .def("add_viz_output", &Model::add_viz_output, py::arg("viz_output"), "- viz_output\n")
      .def("add_count", &Model::add_count, py::arg("count"), "- count\n")
      .def("find_count", &Model::find_count, py::arg("name"), "- name\n")
      .def("load_bngl_observables", &Model::load_bngl_observables, py::arg("file_name"), py::arg("output_files_prefix") = "", py::arg("parameter_overrides") = std::map<std::string, float_t>(), "Loads section observables from a BNGL file and creates Count objects according to it.\nAll elementary molecule types used in the seed species section must be defined in subsystem.\n\n- file_name: BNGL file name.\n\n- output_files_prefix: Prefix to be used when creating files with observable values.\n\n- parameter_overrides\n")
      .def("get_molecule_ids", &Model::get_molecule_ids, py::arg("pattern") = nullptr, "Returns a list of ids of molecules.\nIf the arguments pattern is not set, the list of all molecule ids is returned.  \nIf the argument pattern is set, the list of all molecule ids whose species match \nthe pattern is returned. Matching of patterns with compartments works exactly in the \nsame was as in observables.\n\n- pattern\n")
      .def("get_molecule", &Model::get_molecule, py::arg("id"), "Returns a molecule from the simulated environment, None if the molecule does not exist\n- id\n")
      .def("get_species_name", &Model::get_species_name, py::arg("species_id"), "Returns a string representing canonical species name in the BNGL format. \n\n- species_id\n")
      .def("get_vertex", &Model::get_vertex, py::arg("object"), py::arg("vertex_index"), "Returns coordinates of a vertex.\n- object\n- vertex_index: This is the index of the vertex in object's walls (wall_list).\n\n")
      .def("get_wall", &Model::get_wall, py::arg("object"), py::arg("wall_index"), "Returns information about a wall belonging to a given object.\n- object\n- wall_index: This is the index of the wall in object's walls (wall_list).\n\n")
      .def("get_vertex_unit_normal", &Model::get_vertex_unit_normal, py::arg("object"), py::arg("vertex_index"), "Returns sum of all wall normals that use this vertex converted to a unit vector of length 1um.\nThis represents the unit vector pointing outwards from the vertex.\n\n- object\n- vertex_index: This is the index of the vertex in object's vertex_list.\n\n")
      .def("get_wall_unit_normal", &Model::get_wall_unit_normal, py::arg("object"), py::arg("wall_index"), "Returns wall normal converted to a unit vector of length 1um.\n- object\n- wall_index: This is the index of the vertex in object's walls (wall_list).\n\n")
      .def("dump", &Model::dump)
      .def_property("config", &Model::get_config, &Model::set_config)
      .def_property("warnings", &Model::get_warnings, &Model::set_warnings)
      .def_property("notifications", &Model::get_notifications, &Model::set_notifications)
      .def_property("species", &Model::get_species, &Model::set_species)
      .def_property("reaction_rules", &Model::get_reaction_rules, &Model::set_reaction_rules)
      .def_property("surface_classes", &Model::get_surface_classes, &Model::set_surface_classes)
      .def_property("elementary_molecule_types", &Model::get_elementary_molecule_types, &Model::set_elementary_molecule_types, "Used mainly when a BNGL file is loaded, if BNGL species is defined through \nPython API, this array is populated automatically \n")
      .def_property("release_sites", &Model::get_release_sites, &Model::set_release_sites, "List of release sites to be included in the model.  \n")
      .def_property("geometry_objects", &Model::get_geometry_objects, &Model::set_geometry_objects, "List of geometry objects to be included in the model.  \n")
      .def_property("checkpointed_molecules", &Model::get_checkpointed_molecules, &Model::set_checkpointed_molecules, "Used when resuming simulation from a checkpoint.\n")
      .def_property("viz_outputs", &Model::get_viz_outputs, &Model::set_viz_outputs)
      .def_property("counts", &Model::get_counts, &Model::set_counts)
    ;
}

std::string GenModel::export_to_python(std::ostream& out, PythonExportContext& ctx) {
# if 0 // not to be used
  std::string exported_name = "model";

  bool str_export = export_as_string_without_newlines();
  std::string nl = "";
  std::string ind = " ";
  std::stringstream ss;
  if (!str_export) {
    nl = "\n";
    ind = "    ";
    ss << exported_name << " = ";
  }
  ss << "m.Model(" << nl;
  if (species != std::vector<std::shared_ptr<Species>>() && !skip_vectors_export()) {
    ss << ind << "species = " << export_vec_species(out, ctx, exported_name) << "," << nl;
  }
  if (reaction_rules != std::vector<std::shared_ptr<ReactionRule>>() && !skip_vectors_export()) {
    ss << ind << "reaction_rules = " << export_vec_reaction_rules(out, ctx, exported_name) << "," << nl;
  }
  if (surface_classes != std::vector<std::shared_ptr<SurfaceClass>>() && !skip_vectors_export()) {
    ss << ind << "surface_classes = " << export_vec_surface_classes(out, ctx, exported_name) << "," << nl;
  }
  if (elementary_molecule_types != std::vector<std::shared_ptr<ElementaryMoleculeType>>() && !skip_vectors_export()) {
    ss << ind << "elementary_molecule_types = " << export_vec_elementary_molecule_types(out, ctx, exported_name) << "," << nl;
  }
  if (release_sites != std::vector<std::shared_ptr<ReleaseSite>>() && !skip_vectors_export()) {
    ss << ind << "release_sites = " << export_vec_release_sites(out, ctx, exported_name) << "," << nl;
  }
  if (geometry_objects != std::vector<std::shared_ptr<GeometryObject>>() && !skip_vectors_export()) {
    ss << ind << "geometry_objects = " << export_vec_geometry_objects(out, ctx, exported_name) << "," << nl;
  }
  if (checkpointed_molecules != std::vector<std::shared_ptr<BaseChkptMol>>() && !skip_vectors_export()) {
    ss << ind << "checkpointed_molecules = " << export_vec_checkpointed_molecules(out, ctx, exported_name) << "," << nl;
  }
  if (viz_outputs != std::vector<std::shared_ptr<VizOutput>>() && !skip_vectors_export()) {
    ss << ind << "viz_outputs = " << export_vec_viz_outputs(out, ctx, exported_name) << "," << nl;
  }
  if (counts != std::vector<std::shared_ptr<Count>>() && !skip_vectors_export()) {
    ss << ind << "counts = " << export_vec_counts(out, ctx, exported_name) << "," << nl;
  }
  if (config != Config()) {
    ss << ind << "config = " << config.export_to_python(out, ctx) << "," << nl;
  }
  if (warnings != Warnings()) {
    ss << ind << "warnings = " << warnings.export_to_python(out, ctx) << "," << nl;
  }
  if (notifications != Notifications()) {
    ss << ind << "notifications = " << notifications.export_to_python(out, ctx) << "," << nl;
  }
  ss << ")" << nl << nl;
  if (!str_export) {
    out << ss.str();
    return exported_name;
  }
  else {
    return ss.str();
  }
#else // # if 0
  assert(false);
  return "";
#endif
}

std::string GenModel::export_vec_species(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_species";
  }
  else {
    exported_name = "species";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < species.size(); i++) {
    const auto& item = species[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenModel::export_vec_reaction_rules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_reaction_rules";
  }
  else {
    exported_name = "reaction_rules";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < reaction_rules.size(); i++) {
    const auto& item = reaction_rules[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenModel::export_vec_surface_classes(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_surface_classes";
  }
  else {
    exported_name = "surface_classes";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < surface_classes.size(); i++) {
    const auto& item = surface_classes[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenModel::export_vec_elementary_molecule_types(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_elementary_molecule_types";
  }
  else {
    exported_name = "elementary_molecule_types";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < elementary_molecule_types.size(); i++) {
    const auto& item = elementary_molecule_types[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenModel::export_vec_release_sites(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_release_sites";
  }
  else {
    exported_name = "release_sites";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < release_sites.size(); i++) {
    const auto& item = release_sites[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenModel::export_vec_geometry_objects(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_geometry_objects";
  }
  else {
    exported_name = "geometry_objects";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < geometry_objects.size(); i++) {
    const auto& item = geometry_objects[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenModel::export_vec_checkpointed_molecules(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_checkpointed_molecules";
  }
  else {
    exported_name = "checkpointed_molecules";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < checkpointed_molecules.size(); i++) {
    const auto& item = checkpointed_molecules[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenModel::export_vec_viz_outputs(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_viz_outputs";
  }
  else {
    exported_name = "viz_outputs";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < viz_outputs.size(); i++) {
    const auto& item = viz_outputs[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

std::string GenModel::export_vec_counts(std::ostream& out, PythonExportContext& ctx, const std::string& parent_name) {
  // prints vector into 'out' and returns name of the generated list
  std::stringstream ss;
  std::string exported_name;
  if (parent_name != ""){
    exported_name = parent_name+ "_counts";
  }
  else {
    exported_name = "counts";
  }

  ss << exported_name << " = [\n";
  for (size_t i = 0; i < counts.size(); i++) {
    const auto& item = counts[i];
    if (i == 0) {
      ss << "    ";
    }
    else if (i % 16 == 0) {
      ss << "\n  ";
    }
    if (!item->skip_python_export()) {
      std::string name = item->export_to_python(out, ctx);
      ss << name << ", ";
    }
  }
  ss << "\n]\n\n";
  out << ss.str();
  return exported_name;
}

} // namespace API
} // namespace MCell

