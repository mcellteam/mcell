/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <sstream>
#include "api/pybind11_stl_include.h"
#include "api/python_export_utils.h"
#include "gen_model.h"
#include "api/model.h"
#include "api/base_chkpt_mol.h"
#include "api/color.h"
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

std::shared_ptr<Model> GenModel::copy_model() const {
  std::shared_ptr<Model> res = std::make_shared<Model>(DefaultCtorArgType());
  res->config = config;
  res->warnings = warnings;
  res->notifications = notifications;
  res->species = species;
  res->reaction_rules = reaction_rules;
  res->surface_classes = surface_classes;
  res->elementary_molecule_types = elementary_molecule_types;
  res->release_sites = release_sites;
  res->geometry_objects = geometry_objects;
  res->checkpointed_molecules = checkpointed_molecules;
  res->viz_outputs = viz_outputs;
  res->counts = counts;

  return res;
}

std::shared_ptr<Model> GenModel::deepcopy_model(py::dict) const {
  std::shared_ptr<Model> res = std::make_shared<Model>(DefaultCtorArgType());
  res->config = config;
  res->warnings = warnings;
  res->notifications = notifications;
  for (const auto& item: species) {
    res->species.push_back((is_set(item)) ? item->deepcopy_species() : nullptr);
  }
  for (const auto& item: reaction_rules) {
    res->reaction_rules.push_back((is_set(item)) ? item->deepcopy_reaction_rule() : nullptr);
  }
  for (const auto& item: surface_classes) {
    res->surface_classes.push_back((is_set(item)) ? item->deepcopy_surface_class() : nullptr);
  }
  for (const auto& item: elementary_molecule_types) {
    res->elementary_molecule_types.push_back((is_set(item)) ? item->deepcopy_elementary_molecule_type() : nullptr);
  }
  for (const auto& item: release_sites) {
    res->release_sites.push_back((is_set(item)) ? item->deepcopy_release_site() : nullptr);
  }
  for (const auto& item: geometry_objects) {
    res->geometry_objects.push_back((is_set(item)) ? item->deepcopy_geometry_object() : nullptr);
  }
  for (const auto& item: checkpointed_molecules) {
    res->checkpointed_molecules.push_back((is_set(item)) ? item->deepcopy_base_chkpt_mol() : nullptr);
  }
  for (const auto& item: viz_outputs) {
    res->viz_outputs.push_back((is_set(item)) ? item->deepcopy_viz_output() : nullptr);
  }
  for (const auto& item: counts) {
    res->counts.push_back((is_set(item)) ? item->deepcopy_count() : nullptr);
  }

  return res;
}

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

std::string GenModel::to_str(const bool all_details, const std::string ind) const {
  std::stringstream ss;
  ss << "Model" << ": " <<
      "\n" << ind + "  " << "config=" << "(" << config.to_str(all_details, ind + "  ") << ")" << ", " << "\n" << ind + "  " <<
      "warnings=" << "(" << warnings.to_str(all_details, ind + "  ") << ")" << ", " << "\n" << ind + "  " <<
      "notifications=" << "(" << notifications.to_str(all_details, ind + "  ") << ")" << ", " << "\n" << ind + "  " <<
      "species=" << vec_ptr_to_str(species, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "reaction_rules=" << vec_ptr_to_str(reaction_rules, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "surface_classes=" << vec_ptr_to_str(surface_classes, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "elementary_molecule_types=" << vec_ptr_to_str(elementary_molecule_types, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "release_sites=" << vec_ptr_to_str(release_sites, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "geometry_objects=" << vec_ptr_to_str(geometry_objects, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "checkpointed_molecules=" << vec_ptr_to_str(checkpointed_molecules, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "viz_outputs=" << vec_ptr_to_str(viz_outputs, all_details, ind + "  ") << ", " << "\n" << ind + "  " <<
      "counts=" << vec_ptr_to_str(counts, all_details, ind + "  ");
  return ss.str();
}

py::class_<Model> define_pybinding_Model(py::module& m) {
  return py::class_<Model, std::shared_ptr<Model>>(m, "Model", "This is the main class that is used to assemble all simulation input \nand configuration. It also provides methods to do initialization,\nrun simulation, and introspect the running simulation.\n")
      .def(
          py::init<
          >()
      )
      .def("__copy__", &Model::copy_model)
      .def("__deepcopy__", &Model::deepcopy_model, py::arg("memo"))
      .def("__str__", &Model::to_str, py::arg("all_details") = false, py::arg("ind") = std::string(""))
      .def("__eq__", &Model::__eq__, py::arg("other"))
      .def("initialize", &Model::initialize, py::arg("print_copyright") = true, "Initializes model, initialization blocks most of changes to \ncontained components. \n\n- print_copyright: Prints information about MCell.\n\n")
      .def("run_iterations", &Model::run_iterations, py::arg("iterations"), "Runs specified number of iterations. Returns the number of iterations\nexecuted (it might be less than the requested number of iterations when \na checkpoint was scheduled). \n\n- iterations: Number of iterations to run. Value is truncated to an integer.\n\n")
      .def("end_simulation", &Model::end_simulation, py::arg("print_final_report") = true, "Generates the last visualization and reaction output (if they are included \nin the model), then flushes all buffers and optionally prints simulation report. \nBuffers are also flushed when the Model object is destroyed such as when Ctrl-C\nis pressed during simulation.   \n\n- print_final_report: Print information on simulation time and counts of selected events.\n\n")
      .def("add_subsystem", &Model::add_subsystem, py::arg("subsystem"), "Adds all components of a Subsystem object to the model.\n- subsystem\n")
      .def("add_instantiation", &Model::add_instantiation, py::arg("instantiation"), "Adds all components of an Instantiation object to the model.\n- instantiation\n")
      .def("add_observables", &Model::add_observables, py::arg("observables"), "Adds all counts and viz outputs of an Observables object to the model.\n- observables\n")
      .def("dump_internal_state", &Model::dump_internal_state, py::arg("with_geometry") = false, "Prints out the simulation engine's internal state, mainly for debugging.\n- with_geometry: Include geometry in the dump.\n\n")
      .def("export_data_model", &Model::export_data_model, py::arg("file") = STR_UNSET, "Exports the current state of the model into a data model JSON format.\nDoes not export state of molecules.\nMust be called after model initialization.\nAlways exports the current state, i.e. with the current geometry and reaction rates. \nEvents (ReleaseSites and VizOutputs) with scheduled time other than zero are not exported correctly yet.  \n\n- file: If file is not set, then uses the first VizOutput to determine the target directory \nand creates name using the current iteration. Fails if argument file is not set and \nthere is no VizOutput in the model.\n\n\n")
      .def("export_viz_data_model", &Model::export_viz_data_model, py::arg("file") = STR_UNSET, "Same as export_data_model, only the created data model will contain only information required for visualization\nin CellBlender. This makes the loading of the model by CellBlender faster and also allows to avoid potential\ncompatibility issues.\nMust be called after model initialization.\n\n- file: Optional path to the output data model file.\n\n")
      .def("export_geometry", &Model::export_geometry, py::arg("output_files_prefix") = STR_UNSET, "Exports model geometry as Wavefront OBJ format. \nMust be called after model initialization.\nDoes not export material colors (yet).\n\n- output_files_prefix: Optional prefix for .obj and .mtl files that will be created on export. \nIf output_files_prefix is not set, then uses the first VizOutput to determine the target directory \nand creates names using the current iteration. Fails if argument output_files_prefix is not set and \nthere is no VizOutput in the model.\n\n\n")
      .def("release_molecules", &Model::release_molecules, py::arg("release_site"), "Performs immediate release of molecules based on the definition of the release site argument.\nThe ReleaseSite.release_time must not be in the past and must be within the current iteration \nmeaning that the time must be greater or equal iteration * time_step and less than (iteration + 1) * time_step.\nThe ReleaseEvent must not use a release_pattern because this is an immediate release and it is not \nscheduled into the global scheduler.\n\n- release_site\n")
      .def("run_reaction", &Model::run_reaction, py::arg("reaction_rule"), py::arg("reactant_ids"), py::arg("time"), "Run a single reaction on reactants. Callbacks will be called if they are registered for the given reaction.\nReturns a list of product IDs.\nNote: only unimolecular reactions are currently supported.\n\n- reaction_rule: Reaction rule to run.\n\n- reactant_ids: The number of reactants for a unimolecular reaction must be 1 and for a bimolecular reaction must be 2.\nReactants for a bimolecular reaction do not have to be listed in the same order as in the reaction rule definition. \n\n\n- time: Precise time in seconds when this reaction occurs. Important to know for how long the products\nwill be diffused when they are created in a middle of a time step. \n\n\n")
      .def("add_vertex_move", &Model::add_vertex_move, py::arg("object"), py::arg("vertex_index"), py::arg("displacement"), "Appends information about a displacement for given object's vertex into an internal list of vertex moves. \nTo do the actual geometry change, call Model.apply_vertex_moves.\nThe reason why we first need to collect all changes and then apply them all at the same time is for performance\nreasons. \n\n- object: Object whose vertex will be changed.\n\n- vertex_index: Index of vertex in object's vertex list that will be changed.\n\n- displacement: Change of vertex coordinates [x, y, z] (in um) that will be added to the current \ncoordinates of the vertex.\n\n\n")
      .def("apply_vertex_moves", &Model::apply_vertex_moves, py::arg("collect_wall_wall_hits") = false, "Applies all the vertex moves specified with Model.add_vertex_move call.\nWalls of different objects are checked against collisions and move the maximal way so that they do not \noverlap.\nThe API representation (GeometryObject) is not updated, only the internal MCell data are changed.\nNote: It is not supported yet to move two objects that woudl collide at the same time.  \nWhen collect_wall_wall_hits is True, a list of wall pairs that collided is returned,\nwhen collect_wall_wall_hits is False, and empty list is returned.\n\n- collect_wall_wall_hits: When set to True, a list of wall pairs that collided is returned,\notherwise an empty list is returned.\n\n\n")
      .def("register_mol_wall_hit_callback", &Model::register_mol_wall_hit_callback, py::arg("function"), py::arg("context"), py::arg("object") = nullptr, py::arg("species") = nullptr, "Register a callback for event when a molecule hits a wall. \nMay be called only after model initialization because it internally uses geometry object\nand species ids that are set during the initialization. \n\n- function: Callback function to be called. \nThe function must have two arguments MolWallHitInfo and context.\nDo not modify the received MolWallHitInfo object since it may be reused for other \nwall hit callbacks (e.g. when the first callback is for a specific geometry object and \nthe second callback is for any geometry object). \nThe context object (py::object type argument) is on the other hand provided \nto be modified and one can for instance use it to count the number of hits.. \n\n\n- context: Context passed to the callback function, the callback function can store\ninformation to this object. Some context must be always passed, even when \nit is a useless python object. \n\n\n- object: Only hits of this object will be reported, any object hit is reported when not set.\n\n- species: Only hits of molecules of this species will be reported, any hit of volume molecules of \nany species is reported when this argument is not set.\nSets an internal flag for this species to make sure that the species id does not change \nduring simulation.           \n\n\n")
      .def("register_reaction_callback", &Model::register_reaction_callback, py::arg("function"), py::arg("context"), py::arg("reaction_rule"), "Defines a function to be called when a reaction was processed.\nIt is allowed to do state modifications except for removing reacting molecules, \nthey will be removed automatically after return from this callback. \nUnlimited number of reaction callbacks is allowed. \nMay be called only after model initialization because it internally uses \nreaction rule ids that are set during the initialization. \n\n- function: Callback function to be called. \nThe function must have two arguments ReactionInfo and context.\nCalled right after a reaction occured but before the reactants were removed.\nAfter return the reaction proceeds and reactants are removed (unless they were kept\nby the reaction such as with reaction A + B -> A + C).\n\n\n- context: Context passed to the callback function, the callback function can store\ninformation to this object. Some context must be always passed, even when \nit is a useless python object. \n\n\n- reaction_rule: The callback function will be called whenever this reaction rule is applied.\n\n")
      .def("load_bngl", &Model::load_bngl, py::arg("file_name"), py::arg("observables_path_or_file") = "", py::arg("default_release_region") = nullptr, py::arg("parameter_overrides") = std::map<std::string, double>(), py::arg("observables_output_format") = CountOutputFormat::AUTOMATIC_FROM_EXTENSION, "Loads sections: molecule types, reaction rules, seed species, and observables from a BNGL file\nand creates objects in the current model according to it.\nAll elementary molecule types used in the seed species section must be defined in subsystem.\nIf an item in the seed species section does not have its compartment set,\nthe argument default_region must be set and the molecules are released into or onto the \ndefault_region. \n\n- file_name: Path to the BNGL file to be loaded.\n\n- observables_path_or_file: Directory prefix or file name where observable values will be stored.\nIf a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty \nstring is used,\neach observable gets its own file and the output file format for created Count \nobjects is CountOutputFormat.DAT.\nIf a file has a .gdat extension such as \n'./react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat', all observable are stored in this \nfile and the output file format for created Count objects is CountOutputFormat.GDAT.\nMust not be empty when observables_output_format is explicitly set to CountOutputFormat.GDAT.\n\n\n- default_release_region: Used as region for releases for seed species that have no compartments specified.\n\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n- observables_output_format: Selection of output format. Default setting uses automatic detection\nbased on contents of the 'observables_path_or_file' attribute.\n\n\n")
      .def("export_to_bngl", &Model::export_to_bngl, py::arg("file_name"), py::arg("simulation_method") = BNGSimulationMethod::ODE, "Exports all defined species, reaction rules and applicable observables\nas a BNGL file that can be then loaded by MCell4 or BioNetGen. \nThe resulting file should be validated that it produces expected results. \nMany MCell features cannot be exported into BNGL and when such a feature is \nencountered the export fails with a RuntimeError exception.\nHowever, the export code tries to export as much as possible and one can catch\nthe RuntimeError exception and use the possibly incomplete BNGL file anyway.   \n\n- file_name: Output file name.\n\n- simulation_method: Selection of the BioNetGen simulation method. \nSelects BioNetGen action to run with the selected simulation method.\nFor BNGSimulationMethod.NF the export is limited to a single volume and\na single surface and the enerated rates use volume and surface area so that \nsimulation with NFSim produces corect results. \n\n\n")
      .def("save_checkpoint", &Model::save_checkpoint, py::arg("custom_dir") = STR_UNSET, "Saves current model state as checkpoint. \nThe default directory structure is checkpoints/seed_<SEED>/it_<ITERATION>,\nit can be changed by setting 'custom_dir'.\nIf used during an iteration such as in a callback, an event is scheduled for the  \nbeginning of the next iteration. This scheduled event saves the checkpoint.  \n\n- custom_dir: Sets custom directory where the checkpoint will be stored. \nThe default is 'checkpoints/seed_<SEED>/it_<ITERATION>'. \n\n\n")
      .def("schedule_checkpoint", &Model::schedule_checkpoint, py::arg("iteration") = 0, py::arg("continue_simulation") = false, py::arg("custom_dir") = STR_UNSET, "Schedules checkpoint save event that will occur when an iteration is started.  \nThis means that it will be executed right before any other events scheduled for \nthe given iteration are executed.\nCan be called asynchronously at any time after initialization.\n\n- iteration: Specifies iteration number when the checkpoint save will occur. \nPlease note that iterations are counted from 0.\nTo schedule a checkpoint for the closest time as possible, keep the default value 0,\nthis will schedule checkpoint for the beginning of the iteration with number current iteration + 1.  \nIf calling schedule_checkpoint from a different thread (e.g. by using threading.Timer), \nit is highly recommended to keep the default value 0 or choose some time that will be \nfor sure in the future.\n\n\n- continue_simulation: When false, saving the checkpoint means that we want to terminate the simulation \nright after the save. The currently running function Model.run_iterations\nwill not simulate any following iterations and execution will return from this function\nto execute the next statement which is usually 'model.end_simulation()'.\nWhen true, the checkpoint is saved and simulation continues uninterrupted.\n      \n\n\n- custom_dir: Sets custom directory where the checkpoint will be stored. \nThe default is 'checkpoints/seed_<SEED>/it_<ITERATION>'. \n\n\n")
      .def("add_species", &Model::add_species, py::arg("s"), "Add a reference to a Species object to the species list.\n- s\n")
      .def("find_species", &Model::find_species, py::arg("name"), "Find a Species object using name in the species list. \nReturns None if no such species is found.\n\n- name\n")
      .def("add_reaction_rule", &Model::add_reaction_rule, py::arg("r"), "Add a reference to a ReactionRule object to the reaction_rules list.\n- r\n")
      .def("find_reaction_rule", &Model::find_reaction_rule, py::arg("name"), "Find a ReactionRule object using name in the reaction_rules list. \nReturns None if no such reaction rule is found.\n\n- name\n")
      .def("add_surface_class", &Model::add_surface_class, py::arg("sc"), "Add a reference to a SurfaceClass object to the surface_classes list.\n- sc\n")
      .def("find_surface_class", &Model::find_surface_class, py::arg("name"), "Find a SurfaceClass object using name in the surface_classes list. \nReturns None if no such surface class is found.\n\n- name\n")
      .def("add_elementary_molecule_type", &Model::add_elementary_molecule_type, py::arg("mt"), "Add a reference to an ElementaryMoleculeType object to the elementary_molecule_types list.\n- mt\n")
      .def("find_elementary_molecule_type", &Model::find_elementary_molecule_type, py::arg("name"), "Find an ElementaryMoleculeType object using name in the elementary_molecule_types list. \nReturns None if no such elementary molecule type is found.\n\n- name\n")
      .def("load_bngl_molecule_types_and_reaction_rules", &Model::load_bngl_molecule_types_and_reaction_rules, py::arg("file_name"), py::arg("parameter_overrides") = std::map<std::string, double>(), "Parses a BNGL file, only reads molecule types and reaction rules sections, \ni.e. ignores observables and seed species. \nParameter values are evaluated and the result value is directly used.  \nCompartments names are stored in rxn rules as strings because compartments belong \nto geometry objects and the subsystem is independent on specific geometry.\nHowever, the compartments and their objects must be defined before initialization.\n\n- file_name: Path to the BNGL file to be loaded.\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n")
      .def("add_release_site", &Model::add_release_site, py::arg("s"), "Adds a reference to the release site s to the list of release sites.\n- s\n")
      .def("find_release_site", &Model::find_release_site, py::arg("name"), "Finds a release site by its name, returns None if no such release site is present.\n- name\n")
      .def("add_geometry_object", &Model::add_geometry_object, py::arg("o"), "Adds a reference to the geometry object o to the list of geometry objects.\n- o\n")
      .def("find_geometry_object", &Model::find_geometry_object, py::arg("name"), "Finds a geometry object by its name, returns None if no such geometry object is present.\n- name\n")
      .def("find_volume_compartment_object", &Model::find_volume_compartment_object, py::arg("name"), "Finds a geometry object by its name, the geometry object must be a BNGL compartment.\nReturns None if no such geometry object is present.\n\n- name\n")
      .def("find_surface_compartment_object", &Model::find_surface_compartment_object, py::arg("name"), "Finds a geometry object that is a BNGL compartment and its surface name is name.\nReturns None if no such geometry object is present.\n\n- name\n")
      .def("load_bngl_compartments_and_seed_species", &Model::load_bngl_compartments_and_seed_species, py::arg("file_name"), py::arg("default_release_region") = nullptr, py::arg("parameter_overrides") = std::map<std::string, double>(), "First loads section compartments and for each 3D compartment that does not \nalready exist as a geometry object in this Instantiation object, creates a \nbox with compartment's volume and also sets its 2D (membrane) compartment name.\nWhen multiple identical geometry objects are added to the final Model object, \nonly one copy is left so one can merge multiple Instantiation objects that created \ncompartments assuming that their volume is the same.        \nThen loads section seed species from a BNGL file and creates release sites according to it.\nAll elementary molecule types used in the seed species section must be already defined in subsystem.\nIf an item in the BNGL seed species section does not have its compartment set,\nthe argument default_region must be set and the molecules are then released into or onto the \ndefault_region. \n\n- file_name: Path to the BNGL file.\n\n- default_release_region: Used as region for releases for seed species that have no compartments specified.\n\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n")
      .def("add_viz_output", &Model::add_viz_output, py::arg("viz_output"), "Adds a reference to the viz_output object to the list of visualization output specifications.\n- viz_output\n")
      .def("add_count", &Model::add_count, py::arg("count"), "Adds a reference to the count object to the list of count specifications.\n- count\n")
      .def("find_count", &Model::find_count, py::arg("name"), "Finds a count object by its name, returns None if no such count is present.\n- name\n")
      .def("load_bngl_observables", &Model::load_bngl_observables, py::arg("file_name"), py::arg("observables_path_or_file") = "", py::arg("parameter_overrides") = std::map<std::string, double>(), py::arg("observables_output_format") = CountOutputFormat::AUTOMATIC_FROM_EXTENSION, "Loads section observables from a BNGL file and creates Count objects according to it.\nAll elementary molecule types used in the seed species section must be defined in subsystem.\n\n- file_name: Path to the BNGL file.\n\n- observables_path_or_file: Directory prefix or file name where observable values will be stored.\nIf a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty \nstring is used, each observable gets its own file and the output file format for created Count \nobjects is CountOutputFormat.DAT.\nIf a file has a .gdat extension such as \n'./react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat', all observable are stored in this \nfile and the output file format for created Count objects is CountOutputFormat.GDAT.\nMust not be empty when observables_output_format is explicitly set to CountOutputFormat.GDAT.\n\n\n- parameter_overrides: For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,\nits value is ignored and instead value parameter_overrides[k] is used.\n\n\n- observables_output_format: Selection of output format. Default setting uses automatic detection\nbased on contents of the 'observables_path_or_file' attribute.\n             \n\n\n")
      .def("get_molecule_ids", &Model::get_molecule_ids, py::arg("pattern") = nullptr, "Returns a list of ids of molecules.\nIf the arguments pattern is not set, the list of all molecule ids is returned.  \nIf the argument pattern is set, the list of all molecule ids whose species match \nthe pattern is returned. \n\n- pattern: BNGL pattern to select molecules based on their species, might use compartments.\n\n")
      .def("get_molecule", &Model::get_molecule, py::arg("id"), "Returns a information on a molecule from the simulated environment, \nNone if the molecule does not exist.\n\n- id: Unique id of the molecule to be retrieved.\n\n")
      .def("get_species_name", &Model::get_species_name, py::arg("species_id"), "Returns a string representing canonical species name in the BNGL format.\n\n- species_id: Id of the species.\n\n")
      .def("get_vertex", &Model::get_vertex, py::arg("object"), py::arg("vertex_index"), "Returns coordinates of a vertex.\n- object\n- vertex_index: This is the index of the vertex in the geometry object's walls (wall_list).\n\n")
      .def("get_wall", &Model::get_wall, py::arg("object"), py::arg("wall_index"), "Returns information about a wall belonging to a given object.\n- object: Geometry object whose wall to retrieve.\n\n- wall_index: This is the index of the wall in the geometry object's walls (wall_list).\n\n")
      .def("get_vertex_unit_normal", &Model::get_vertex_unit_normal, py::arg("object"), py::arg("vertex_index"), "Returns sum of all wall normals that use this vertex converted to a unit vector of \nlength 1 um (micrometer).\nThis represents the unit vector pointing outwards from the vertex.\n\n- object: Geometry object whose vertex to retrieve.\n\n- vertex_index: This is the index of the vertex in the geometry object's vertex_list.\n\n")
      .def("get_wall_unit_normal", &Model::get_wall_unit_normal, py::arg("object"), py::arg("wall_index"), "Returns wall normal converted to a unit vector of length 1um.\n- object: Geometry object whose wall's normal to retrieve.\n\n- wall_index: This is the index of the vertex in the geometry object's walls (wall_list).\n\n")
      .def("get_wall_color", &Model::get_wall_color, py::arg("object"), py::arg("wall_index"), "Returns color of a wall.\n- object: Geometry object whose wall's color to retrieve.\n\n- wall_index: This is the index of the vertex in the geometry object's walls (wall_list).\n\n")
      .def("set_wall_color", &Model::set_wall_color, py::arg("object"), py::arg("wall_index"), py::arg("color"), "Sets color of a wall.\n- object: Geometry object whose wall's color to retrieve.\n\n- wall_index: This is the index of the vertex in the geometry object's walls (wall_list).\n\n- color: Color to be set.\n\n")
      .def("dump", &Model::dump)
      .def_property("config", &Model::get_config, &Model::set_config, "Simulation configuration.")
      .def_property("warnings", &Model::get_warnings, &Model::set_warnings, "Configuration on how to report warnings.")
      .def_property("notifications", &Model::get_notifications, &Model::set_notifications, "Configuration on how to report certain notifications.")
      .def_property("species", &Model::get_species, &Model::set_species, py::return_value_policy::reference, "List of species to be included in the model for initialization.\nUsed usually only for simple species (species that are defined using a\nsingle molecule type without components such as 'A').\nOther species may be created inside simulation  \n")
      .def_property("reaction_rules", &Model::get_reaction_rules, &Model::set_reaction_rules, py::return_value_policy::reference)
      .def_property("surface_classes", &Model::get_surface_classes, &Model::set_surface_classes, py::return_value_policy::reference)
      .def_property("elementary_molecule_types", &Model::get_elementary_molecule_types, &Model::set_elementary_molecule_types, py::return_value_policy::reference, "Contains list of elementary molecule types with their diffusion constants and other information. \nPopulated when a BNGL file is loaded and also on initialization from Species objects present in \nthe species list.\n")
      .def_property("release_sites", &Model::get_release_sites, &Model::set_release_sites, py::return_value_policy::reference, "List of release sites to be included in the model.  \n")
      .def_property("geometry_objects", &Model::get_geometry_objects, &Model::set_geometry_objects, py::return_value_policy::reference, "List of geometry objects to be included in the model.  \n")
      .def_property("checkpointed_molecules", &Model::get_checkpointed_molecules, &Model::set_checkpointed_molecules, py::return_value_policy::reference, "Used when resuming simulation from a checkpoint.\n")
      .def_property("viz_outputs", &Model::get_viz_outputs, &Model::set_viz_outputs, py::return_value_policy::reference, "List of visualization outputs to be included in the model.\nThere is usually just one VizOutput object.   \n")
      .def_property("counts", &Model::get_counts, &Model::set_counts, py::return_value_policy::reference, "List of counts to be included in the model.\n")
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

