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

#ifndef API_GEN_NAMES
#define API_GEN_NAMES


namespace MCell {
namespace API {

const char* const NAME_CLASS_COMPLEX_INSTANCE = "ComplexInstance";
const char* const NAME_CLASS_COMPONENT_INSTANCE = "ComponentInstance";
const char* const NAME_CLASS_COMPONENT_TYPE = "ComponentType";
const char* const NAME_CLASS_CONFIG = "Config";
const char* const NAME_CLASS_COUNT = "Count";
const char* const NAME_CLASS_COUNT_TERM = "CountTerm";
const char* const NAME_CLASS_GEOMETRY_OBJECT = "GeometryObject";
const char* const NAME_CLASS_INSTANTIATION_DATA = "InstantiationData";
const char* const NAME_CLASS_MODEL = "Model";
const char* const NAME_CLASS_MOLECULE_INSTANCE = "MoleculeInstance";
const char* const NAME_CLASS_MOLECULE_TYPE = "MoleculeType";
const char* const NAME_CLASS_NOTIFICATIONS = "Notifications";
const char* const NAME_CLASS_OBSERVABLES = "Observables";
const char* const NAME_CLASS_REACTION_RULE = "ReactionRule";
const char* const NAME_CLASS_RELEASE_SITE = "ReleaseSite";
const char* const NAME_CLASS_SPECIES = "Species";
const char* const NAME_CLASS_SUBSYSTEM = "Subsystem";
const char* const NAME_CLASS_SURFACE_REGION = "SurfaceRegion";
const char* const NAME_CLASS_VIZ_OUTPUT = "VizOutput";
const char* const NAME_CLASS_WARNINGS = "Warnings";

const char* const NAME___ADD__ = "__add__";
const char* const NAME___SUB__ = "__sub__";
const char* const NAME_ADD_COUNT = "add_count";
const char* const NAME_ADD_GEOMETRY_OBJECT = "add_geometry_object";
const char* const NAME_ADD_INSTANTIATION_DATA = "add_instantiation_data";
const char* const NAME_ADD_OBSERVABLES = "add_observables";
const char* const NAME_ADD_REACTION_RULE = "add_reaction_rule";
const char* const NAME_ADD_RELEASE_SITE = "add_release_site";
const char* const NAME_ADD_SPECIES = "add_species";
const char* const NAME_ADD_SUBSYSTEM = "add_subsystem";
const char* const NAME_ADD_VIZ_OUTPUT = "add_viz_output";
const char* const NAME_BOND = "bond";
const char* const NAME_CENTER_MOLECULES_ON_GRID = "center_molecules_on_grid";
const char* const NAME_COMPONENT_TYPE = "component_type";
const char* const NAME_COMPONENTS = "components";
const char* const NAME_CONFIG = "config";
const char* const NAME_COUNT = "count";
const char* const NAME_COUNT_EXPRESSION = "count_expression";
const char* const NAME_COUNTS = "counts";
const char* const NAME_DEGENERATE_POLYGONS = "degenerate_polygons";
const char* const NAME_DIFFUSION_CONSTANT_2D = "diffusion_constant_2d";
const char* const NAME_DIFFUSION_CONSTANT_3D = "diffusion_constant_3d";
const char* const NAME_DIFFUSION_CONSTANT_REPORT = "diffusion_constant_report";
const char* const NAME_DUMP_INTERNAL_STATE = "dump_internal_state";
const char* const NAME_ELEMENT_CONNECTIONS = "element_connections";
const char* const NAME_END_SIMULATION = "end_simulation";
const char* const NAME_EVERY_N_TIMESTEPS = "every_n_timesteps";
const char* const NAME_FILENAME = "filename";
const char* const NAME_FILENAME_PREFIX = "filename_prefix";
const char* const NAME_FINAL_SUMMARY = "final_summary";
const char* const NAME_FIND_GEOMETRY_OBJECT = "find_geometry_object";
const char* const NAME_FIND_REACTION_RULE = "find_reaction_rule";
const char* const NAME_FIND_RELEASE_SITE = "find_release_site";
const char* const NAME_FIND_SPECIES = "find_species";
const char* const NAME_FWD_RATE = "fwd_rate";
const char* const NAME_GEOMETRY_OBJECTS = "geometry_objects";
const char* const NAME_HIGH_REACTION_PROBABILITY = "high_reaction_probability";
const char* const NAME_INITIAL_ORIENTATION = "initial_orientation";
const char* const NAME_INITIALIZE = "initialize";
const char* const NAME_INST = "inst";
const char* const NAME_INSTANTIATION_DATA = "instantiation_data";
const char* const NAME_INTERACTION_RADIUS = "interaction_radius";
const char* const NAME_ITERATION_REPORT = "iteration_report";
const char* const NAME_ITERATIONS = "iterations";
const char* const NAME_LEFT_NODE = "left_node";
const char* const NAME_LIFETIME_THRESHOLD = "lifetime_threshold";
const char* const NAME_LIFETIME_TOO_SHORT = "lifetime_too_short";
const char* const NAME_LOCATION = "location";
const char* const NAME_MICROSCOPIC_REVERSIBILITY = "microscopic_reversibility";
const char* const NAME_MISSED_REACTIONS = "missed_reactions";
const char* const NAME_MISSED_REACTIONS_THRESHOLD = "missed_reactions_threshold";
const char* const NAME_MISSING_SURFACE_ORIENTATION = "missing_surface_orientation";
const char* const NAME_MODE = "mode";
const char* const NAME_MOLECULE_COLLISION_REPORT = "molecule_collision_report";
const char* const NAME_MOLECULE_INSTANCES = "molecule_instances";
const char* const NAME_MOLECULE_TYPE = "molecule_type";
const char* const NAME_NAME = "name";
const char* const NAME_NEGATIVE_DIFFUSION_CONSTANT = "negative_diffusion_constant";
const char* const NAME_NEGATIVE_REACTION_RATE = "negative_reaction_rate";
const char* const NAME_NODE_TYPE = "node_type";
const char* const NAME_NOTIFICATIONS = "notifications";
const char* const NAME_NUMBER_TO_RELEASE = "number_to_release";
const char* const NAME_O = "o";
const char* const NAME_OBSERVABLES = "observables";
const char* const NAME_OP2 = "op2";
const char* const NAME_ORIENTATION = "orientation";
const char* const NAME_PARTITION_DIMENSION = "partition_dimension";
const char* const NAME_PRINT_FINAL_REPORT = "print_final_report";
const char* const NAME_PROBABILITY_REPORT = "probability_report";
const char* const NAME_PRODUCTS = "products";
const char* const NAME_PROGRESS_REPORT = "progress_report";
const char* const NAME_R = "r";
const char* const NAME_REACTANTS = "reactants";
const char* const NAME_REACTION_RULE = "reaction_rule";
const char* const NAME_REACTION_RULES = "reaction_rules";
const char* const NAME_REGION = "region";
const char* const NAME_RELEASE_EVENT_REPORT = "release_event_report";
const char* const NAME_RELEASE_PROBABILITY = "release_probability";
const char* const NAME_RELEASE_SITES = "release_sites";
const char* const NAME_REV_RATE = "rev_rate";
const char* const NAME_RIGHT_NODE = "right_node";
const char* const NAME_RUN_ITERATIONS = "run_iterations";
const char* const NAME_S = "s";
const char* const NAME_SEED = "seed";
const char* const NAME_SHAPE = "shape";
const char* const NAME_SITE_DIAMETER = "site_diameter";
const char* const NAME_SITE_RADIUS = "site_radius";
const char* const NAME_SPECIES = "species";
const char* const NAME_SPECIES_LIST = "species_list";
const char* const NAME_STATE = "state";
const char* const NAME_STATES = "states";
const char* const NAME_SUBPARTITION_DIMENSION = "subpartition_dimension";
const char* const NAME_SUBSYSTEM = "subsystem";
const char* const NAME_SURFACE_GRID_DENSITY = "surface_grid_density";
const char* const NAME_SURFACE_REGIONS = "surface_regions";
const char* const NAME_TIME_STEP = "time_step";
const char* const NAME_TOTAL_ITERATIONS_HINT = "total_iterations_hint";
const char* const NAME_USELESS_VOLUME_ORIENTATION = "useless_volume_orientation";
const char* const NAME_VACANCY_SEARCH_DISTANCE = "vacancy_search_distance";
const char* const NAME_VARYING_PROBABILITY_REPORT = "varying_probability_report";
const char* const NAME_VERTEX_LIST = "vertex_list";
const char* const NAME_VIZ_OUTPUT = "viz_output";
const char* const NAME_VIZ_OUTPUTS = "viz_outputs";
const char* const NAME_WARNINGS = "warnings";

const char* const NAME_ENUM_WARNING_LEVEL = "WarningLevel";
const char* const NAME_ENUM_ORIENTATION = "Orientation";
const char* const NAME_ENUM_NOTIFICATION = "Notification";
const char* const NAME_ENUM_EXPR_NODE_TYPE = "ExprNodeType";
const char* const NAME_ENUM_VIZ_MODE = "VizMode";
const char* const NAME_ENUM_SHAPE = "Shape";

const char* const NAME_ADD = "Add";
const char* const NAME_ASCII = "Ascii";
const char* const NAME_BRIEF = "Brief";
const char* const NAME_CELLBLENDER = "Cellblender";
const char* const NAME_DOWN = "Down";
const char* const NAME_ERROR = "Error";
const char* const NAME_FULL = "Full";
const char* const NAME_IGNORE = "Ignore";
const char* const NAME_LEAF = "Leaf";
const char* const NAME_NONE = "None";
const char* const NAME_NOT_SET = "NotSet";
const char* const NAME_SPHERICAL = "Spherical";
const char* const NAME_SUB = "Sub";
const char* const NAME_UNSET = "Unset";
const char* const NAME_UP = "Up";
const char* const NAME_WARNING = "Warning";

const char* const NAME_BOND_BOUND = "BOND_BOUND";
const char* const NAME_BOND_UNBOUND = "BOND_UNBOUND";
const char* const NAME_PARTITION_EDGE_EXTRA_MARGIN_UM = "PARTITION_EDGE_EXTRA_MARGIN_UM";
const char* const NAME_STATE_UNSET = "STATE_UNSET";
const char* const NAME_STATE_UNSET_INT = "STATE_UNSET_INT";

} // namespace API
} // namespace MCell

#endif // API_GEN_NAMES
