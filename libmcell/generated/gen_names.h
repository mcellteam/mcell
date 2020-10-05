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

const char* const NAME_CLASS_BNGL_UTILS = "bngl_utils";
const char* const NAME_CLASS_COMPLEX_INSTANCE = "ComplexInstance";
const char* const NAME_CLASS_COMPONENT_INSTANCE = "ComponentInstance";
const char* const NAME_CLASS_COMPONENT_TYPE = "ComponentType";
const char* const NAME_CLASS_CONFIG = "Config";
const char* const NAME_CLASS_COUNT = "Count";
const char* const NAME_CLASS_COUNT_TERM = "CountTerm";
const char* const NAME_CLASS_ELEMENTARY_MOLECULE_INSTANCE = "ElementaryMoleculeInstance";
const char* const NAME_CLASS_ELEMENTARY_MOLECULE_TYPE = "ElementaryMoleculeType";
const char* const NAME_CLASS_GEOMETRY_UTILS = "geometry_utils";
const char* const NAME_CLASS_GEOMETRY_OBJECT = "GeometryObject";
const char* const NAME_CLASS_INITIAL_SURFACE_RELEASE = "InitialSurfaceRelease";
const char* const NAME_CLASS_INSTANTIATION_DATA = "InstantiationData";
const char* const NAME_CLASS_MODEL = "Model";
const char* const NAME_CLASS_MOLECULE = "Molecule";
const char* const NAME_CLASS_MOLECULE_RELEASE_INFO = "MoleculeReleaseInfo";
const char* const NAME_CLASS_NOTIFICATIONS = "Notifications";
const char* const NAME_CLASS_OBSERVABLES = "Observables";
const char* const NAME_CLASS_REACTION_RULE = "ReactionRule";
const char* const NAME_CLASS_REGION = "Region";
const char* const NAME_CLASS_RELEASE_PATTERN = "ReleasePattern";
const char* const NAME_CLASS_RELEASE_SITE = "ReleaseSite";
const char* const NAME_CLASS_SPECIES = "Species";
const char* const NAME_CLASS_SUBSYSTEM = "Subsystem";
const char* const NAME_CLASS_SURFACE_CLASS = "SurfaceClass";
const char* const NAME_CLASS_SURFACE_PROPERTY = "SurfaceProperty";
const char* const NAME_CLASS_SURFACE_REGION = "SurfaceRegion";
const char* const NAME_CLASS_VIZ_OUTPUT = "VizOutput";
const char* const NAME_CLASS_WALL_HIT_INFO = "WallHitInfo";
const char* const NAME_CLASS_WARNINGS = "Warnings";

const char* const NAME___ADD__ = "__add__";
const char* const NAME___MUL__ = "__mul__";
const char* const NAME___SUB__ = "__sub__";
const char* const NAME_ADD_COUNT = "add_count";
const char* const NAME_ADD_ELEMENTARY_MOLECULE_TYPE = "add_elementary_molecule_type";
const char* const NAME_ADD_GEOMETRY_OBJECT = "add_geometry_object";
const char* const NAME_ADD_INSTANTIATION_DATA = "add_instantiation_data";
const char* const NAME_ADD_OBSERVABLES = "add_observables";
const char* const NAME_ADD_REACTION_RULE = "add_reaction_rule";
const char* const NAME_ADD_RELEASE_SITE = "add_release_site";
const char* const NAME_ADD_SPECIES = "add_species";
const char* const NAME_ADD_SUBSYSTEM = "add_subsystem";
const char* const NAME_ADD_SURFACE_CLASS = "add_surface_class";
const char* const NAME_ADD_VERTEX_MOVE = "add_vertex_move";
const char* const NAME_ADD_VIZ_OUTPUT = "add_viz_output";
const char* const NAME_AFFECTED_SPECIES = "affected_species";
const char* const NAME_ALL_SPECIES = "all_species";
const char* const NAME_APPLY_VERTEX_MOVES = "apply_vertex_moves";
const char* const NAME_AS_SPECIES = "as_species";
const char* const NAME_BNG_VERBOSITY_LEVEL = "bng_verbosity_level";
const char* const NAME_BNGL_SPECIES = "bngl_species";
const char* const NAME_BOND = "bond";
const char* const NAME_CENTER_MOLECULES_ON_GRID = "center_molecules_on_grid";
const char* const NAME_CHECK_OVERLAPPED_WALLS = "check_overlapped_walls";
const char* const NAME_COMPLEX_INSTANCE = "complex_instance";
const char* const NAME_COMPONENT_TYPE = "component_type";
const char* const NAME_COMPONENTS = "components";
const char* const NAME_CONCENTRATION = "concentration";
const char* const NAME_CONFIG = "config";
const char* const NAME_CONTEXT = "context";
const char* const NAME_COUNT = "count";
const char* const NAME_COUNT_EXPRESSION = "count_expression";
const char* const NAME_COUNTS = "counts";
const char* const NAME_CREATE_BOX = "create_box";
const char* const NAME_CUSTOM_SPACE_STEP = "custom_space_step";
const char* const NAME_CUSTOM_TIME_STEP = "custom_time_step";
const char* const NAME_DEFAULT_RELEASE_REGION = "default_release_region";
const char* const NAME_DEGENERATE_POLYGONS = "degenerate_polygons";
const char* const NAME_DENSITY = "density";
const char* const NAME_DIFFUSION_CONSTANT_2D = "diffusion_constant_2d";
const char* const NAME_DIFFUSION_CONSTANT_3D = "diffusion_constant_3d";
const char* const NAME_DIFFUSION_CONSTANT_REPORT = "diffusion_constant_report";
const char* const NAME_DISPLACEMENT = "displacement";
const char* const NAME_DUMP_INTERNAL_STATE = "dump_internal_state";
const char* const NAME_EDGE_LENGTH = "edge_length";
const char* const NAME_ELEMENT_CONNECTIONS = "element_connections";
const char* const NAME_ELEMENTARY_MOLECULE_INSTANCES = "elementary_molecule_instances";
const char* const NAME_ELEMENTARY_MOLECULE_TYPE = "elementary_molecule_type";
const char* const NAME_ELEMENTARY_MOLECULE_TYPES = "elementary_molecule_types";
const char* const NAME_END_SIMULATION = "end_simulation";
const char* const NAME_EVERY_N_TIMESTEPS = "every_n_timesteps";
const char* const NAME_EXPORT_DATA_MODEL = "export_data_model";
const char* const NAME_EXPORT_VIZ_DATA_MODEL = "export_viz_data_model";
const char* const NAME_FILE = "file";
const char* const NAME_FILE_NAME = "file_name";
const char* const NAME_FINAL_SUMMARY = "final_summary";
const char* const NAME_FIND_ELEMENTARY_MOLECULE_TYPE = "find_elementary_molecule_type";
const char* const NAME_FIND_GEOMETRY_OBJECT = "find_geometry_object";
const char* const NAME_FIND_REACTION_RULE = "find_reaction_rule";
const char* const NAME_FIND_RELEASE_SITE = "find_release_site";
const char* const NAME_FIND_SPECIES = "find_species";
const char* const NAME_FIND_SURFACE_CLASS = "find_surface_class";
const char* const NAME_FUNCTION = "function";
const char* const NAME_FWD_RATE = "fwd_rate";
const char* const NAME_GEOMETRY_OBJECT_ID = "geometry_object_id";
const char* const NAME_GEOMETRY_OBJECTS = "geometry_objects";
const char* const NAME_GET_MOLECULE = "get_molecule";
const char* const NAME_GET_MOLECULE_IDS = "get_molecule_ids";
const char* const NAME_HIGH_REACTION_PROBABILITY = "high_reaction_probability";
const char* const NAME_ID = "id";
const char* const NAME_INDEX = "index";
const char* const NAME_INITIAL_PARTITION_ORIGIN = "initial_partition_origin";
const char* const NAME_INITIAL_SURFACE_RELEASES = "initial_surface_releases";
const char* const NAME_INITIALIZE = "initialize";
const char* const NAME_INST = "inst";
const char* const NAME_INSTANTIATION_DATA = "instantiation_data";
const char* const NAME_INTERACTION_RADIUS = "interaction_radius";
const char* const NAME_ITERATION_REPORT = "iteration_report";
const char* const NAME_ITERATIONS = "iterations";
const char* const NAME_LEFT_NODE = "left_node";
const char* const NAME_LIFETIME_THRESHOLD = "lifetime_threshold";
const char* const NAME_LIFETIME_TOO_SHORT = "lifetime_too_short";
const char* const NAME_LOAD_BNGL = "load_bngl";
const char* const NAME_LOAD_BNGL_MOLECULE_TYPES_AND_REACTION_RULES = "load_bngl_molecule_types_and_reaction_rules";
const char* const NAME_LOAD_BNGL_OBSERVABLES = "load_bngl_observables";
const char* const NAME_LOAD_BNGL_PARAMETERS = "load_bngl_parameters";
const char* const NAME_LOAD_BNGL_SEED_SPECIES = "load_bngl_seed_species";
const char* const NAME_LOCATION = "location";
const char* const NAME_MISSED_REACTIONS = "missed_reactions";
const char* const NAME_MISSED_REACTIONS_THRESHOLD = "missed_reactions_threshold";
const char* const NAME_MISSING_SURFACE_ORIENTATION = "missing_surface_orientation";
const char* const NAME_MODE = "mode";
const char* const NAME_MOLECULE_COLLISION_REPORT = "molecule_collision_report";
const char* const NAME_MOLECULE_ID = "molecule_id";
const char* const NAME_MOLECULE_LIST = "molecule_list";
const char* const NAME_MOLECULES_PATTERN = "molecules_pattern";
const char* const NAME_MT = "mt";
const char* const NAME_MULTIPLIER = "multiplier";
const char* const NAME_NAME = "name";
const char* const NAME_NEGATIVE_DIFFUSION_CONSTANT = "negative_diffusion_constant";
const char* const NAME_NEGATIVE_REACTION_RATE = "negative_reaction_rate";
const char* const NAME_NODE_TYPE = "node_type";
const char* const NAME_NOTIFICATIONS = "notifications";
const char* const NAME_NUMBER_OF_TRAINS = "number_of_trains";
const char* const NAME_NUMBER_TO_RELEASE = "number_to_release";
const char* const NAME_O = "o";
const char* const NAME_OBJECT = "object";
const char* const NAME_OBSERVABLES = "observables";
const char* const NAME_OBSERVABLES_FILES_PREFIX = "observables_files_prefix";
const char* const NAME_OP2 = "op2";
const char* const NAME_ORIENTATION = "orientation";
const char* const NAME_OTHER = "other";
const char* const NAME_OUTPUT_FILES_PREFIX = "output_files_prefix";
const char* const NAME_PARAMETER_OVERRIDES = "parameter_overrides";
const char* const NAME_PARTITION_DIMENSION = "partition_dimension";
const char* const NAME_POS = "pos";
const char* const NAME_POS3D = "pos3d";
const char* const NAME_POS_BEFORE_HIT = "pos_before_hit";
const char* const NAME_PRINT_FINAL_REPORT = "print_final_report";
const char* const NAME_PROBABILITY_REPORT = "probability_report";
const char* const NAME_PRODUCTS = "products";
const char* const NAME_PROGRESS_REPORT = "progress_report";
const char* const NAME_PROPERTIES = "properties";
const char* const NAME_R = "r";
const char* const NAME_REACTANTS = "reactants";
const char* const NAME_REACTION_RULE = "reaction_rule";
const char* const NAME_REACTION_RULES = "reaction_rules";
const char* const NAME_REGION = "region";
const char* const NAME_REGISTER_WALL_HIT_CALLBACK = "register_wall_hit_callback";
const char* const NAME_RELEASE_EVENT_REPORT = "release_event_report";
const char* const NAME_RELEASE_INTERVAL = "release_interval";
const char* const NAME_RELEASE_PATTERN = "release_pattern";
const char* const NAME_RELEASE_PROBABILITY = "release_probability";
const char* const NAME_RELEASE_SITES = "release_sites";
const char* const NAME_RELEASE_TIME = "release_time";
const char* const NAME_REMOVE = "remove";
const char* const NAME_REV_NAME = "rev_name";
const char* const NAME_REV_RATE = "rev_rate";
const char* const NAME_RIGHT_NODE = "right_node";
const char* const NAME_RUN_ITERATIONS = "run_iterations";
const char* const NAME_RXN_AND_SPECIES_REPORT = "rxn_and_species_report";
const char* const NAME_S = "s";
const char* const NAME_SC = "sc";
const char* const NAME_SEED = "seed";
const char* const NAME_SHAPE = "shape";
const char* const NAME_SITE_DIAMETER = "site_diameter";
const char* const NAME_SITE_RADIUS = "site_radius";
const char* const NAME_SPECIES = "species";
const char* const NAME_SPECIES_LIST = "species_list";
const char* const NAME_SPECIES_PATTERN = "species_pattern";
const char* const NAME_STATE = "state";
const char* const NAME_STATES = "states";
const char* const NAME_SUBPARTITION_DIMENSION = "subpartition_dimension";
const char* const NAME_SUBSYSTEM = "subsystem";
const char* const NAME_SURFACE_CLASS = "surface_class";
const char* const NAME_SURFACE_CLASSES = "surface_classes";
const char* const NAME_SURFACE_GRID_DENSITY = "surface_grid_density";
const char* const NAME_SURFACE_REGIONS = "surface_regions";
const char* const NAME_TARGET_ONLY = "target_only";
const char* const NAME_TIME = "time";
const char* const NAME_TIME_BEFORE_HIT = "time_before_hit";
const char* const NAME_TIME_STEP = "time_step";
const char* const NAME_TO_BNGL_STR = "to_bngl_str";
const char* const NAME_TOTAL_ITERATIONS_HINT = "total_iterations_hint";
const char* const NAME_TRAIN_DURATION = "train_duration";
const char* const NAME_TRAIN_INTERVAL = "train_interval";
const char* const NAME_TYPE = "type";
const char* const NAME_USELESS_VOLUME_ORIENTATION = "useless_volume_orientation";
const char* const NAME_VACANCY_SEARCH_DISTANCE = "vacancy_search_distance";
const char* const NAME_VARIABLE_RATE = "variable_rate";
const char* const NAME_VARYING_PROBABILITY_REPORT = "varying_probability_report";
const char* const NAME_VERTEX_LIST = "vertex_list";
const char* const NAME_VIZ_OUTPUT = "viz_output";
const char* const NAME_VIZ_OUTPUTS = "viz_outputs";
const char* const NAME_WALL_ID = "wall_id";
const char* const NAME_WALL_INDICES = "wall_indices";
const char* const NAME_WARNINGS = "warnings";

const char* const NAME_ENUM_EXPR_NODE_TYPE = "ExprNodeType";
const char* const NAME_ENUM_NOTIFICATION = "Notification";
const char* const NAME_ENUM_ORIENTATION = "Orientation";
const char* const NAME_ENUM_REGION_NODE_TYPE = "RegionNodeType";
const char* const NAME_ENUM_SHAPE = "Shape";
const char* const NAME_ENUM_SURFACE_PROPERTY_TYPE = "SurfacePropertyType";
const char* const NAME_ENUM_VIZ_MODE = "VizMode";
const char* const NAME_ENUM_WARNING_LEVEL = "WarningLevel";

const char* const NAME_EV_ABSORPTIVE = "ABSORPTIVE";
const char* const NAME_EV_ADD = "ADD";
const char* const NAME_EV_ANY = "ANY";
const char* const NAME_EV_ASCII = "ASCII";
const char* const NAME_EV_BRIEF = "BRIEF";
const char* const NAME_EV_CELLBLENDER = "CELLBLENDER";
const char* const NAME_EV_DIFFERENCE = "DIFFERENCE";
const char* const NAME_EV_DOWN = "DOWN";
const char* const NAME_EV_ERROR = "ERROR";
const char* const NAME_EV_FULL = "FULL";
const char* const NAME_EV_IGNORE = "IGNORE";
const char* const NAME_EV_INTERSECT = "INTERSECT";
const char* const NAME_EV_LEAF = "LEAF";
const char* const NAME_EV_LEAF_GEOMETRY_OBJECT = "LEAF_GEOMETRY_OBJECT";
const char* const NAME_EV_LEAF_SURFACE_REGION = "LEAF_SURFACE_REGION";
const char* const NAME_EV_LIST = "LIST";
const char* const NAME_EV_NONE = "NONE";
const char* const NAME_EV_NOT_SET = "NOT_SET";
const char* const NAME_EV_REFLECTIVE = "REFLECTIVE";
const char* const NAME_EV_REGION_EXPR = "REGION_EXPR";
const char* const NAME_EV_SPHERICAL = "SPHERICAL";
const char* const NAME_EV_SUB = "SUB";
const char* const NAME_EV_TRANSPARENT = "TRANSPARENT";
const char* const NAME_EV_UNION = "UNION";
const char* const NAME_EV_UNSET = "UNSET";
const char* const NAME_EV_UP = "UP";
const char* const NAME_EV_WARNING = "WARNING";

const char* const NAME_CV_ALL_MOLECULES = "ALL_MOLECULES";
const char* const NAME_CV_ALL_SURFACE_MOLECULES = "ALL_SURFACE_MOLECULES";
const char* const NAME_CV_ALL_VOLUME_MOLECULES = "ALL_VOLUME_MOLECULES";
const char* const NAME_CV_AllMolecules = "AllMolecules";
const char* const NAME_CV_AllSurfaceMolecules = "AllSurfaceMolecules";
const char* const NAME_CV_AllVolumeMolecules = "AllVolumeMolecules";
const char* const NAME_CV_BOND_ANY = "BOND_ANY";
const char* const NAME_CV_BOND_BOUND = "BOND_BOUND";
const char* const NAME_CV_BOND_UNBOUND = "BOND_UNBOUND";
const char* const NAME_CV_DEFAULT_COUNT_BUFFER_SIZE = "DEFAULT_COUNT_BUFFER_SIZE";
const char* const NAME_CV_MOLECULE_ID_INVALID = "MOLECULE_ID_INVALID";
const char* const NAME_CV_NUMBER_OF_TRAINS_UNLIMITED = "NUMBER_OF_TRAINS_UNLIMITED";
const char* const NAME_CV_PARTITION_EDGE_EXTRA_MARGIN_UM = "PARTITION_EDGE_EXTRA_MARGIN_UM";
const char* const NAME_CV_STATE_UNSET = "STATE_UNSET";
const char* const NAME_CV_STATE_UNSET_INT = "STATE_UNSET_INT";
const char* const NAME_CV_TIME_INFINITY = "TIME_INFINITY";

} // namespace API
} // namespace MCell

#endif // API_GEN_NAMES

