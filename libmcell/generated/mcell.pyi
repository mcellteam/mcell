from typing import List, Dict, Callable, Any
from enum import Enum

INT32_MAX = 2147483647 # do not use this constant in your code

FLT_MAX = 3.40282346638528859812e+38 # do not use this constant in your code

# "forward" declarations to make the type hints valid
class BaseChkptMol():
    pass
class ChkptSurfMol():
    pass
class ChkptVolMol():
    pass
class Color():
    pass
class Complex():
    pass
class Component():
    pass
class ComponentType():
    pass
class Config():
    pass
class Count():
    pass
class CountTerm():
    pass
class ElementaryMolecule():
    pass
class ElementaryMoleculeType():
    pass
class GeometryObject():
    pass
class InitialSurfaceRelease():
    pass
class Instantiation():
    pass
class Introspection():
    pass
class Model():
    pass
class MolWallHitInfo():
    pass
class Molecule():
    pass
class MoleculeReleaseInfo():
    pass
class Notifications():
    pass
class Observables():
    pass
class ReactionInfo():
    pass
class ReactionRule():
    pass
class Region():
    pass
class ReleasePattern():
    pass
class ReleaseSite():
    pass
class RngState():
    pass
class Species():
    pass
class Subsystem():
    pass
class SurfaceClass():
    pass
class SurfaceProperty():
    pass
class SurfaceRegion():
    pass
class VizOutput():
    pass
class Wall():
    pass
class WallWallHitInfo():
    pass
class Warnings():
    pass
class bngl_utils():
    pass
class data_utils():
    pass
class geometry_utils():
    pass
class run_utils():
    pass

class Vec3():
    def __init__(self, x : float = 0, y : float = 0, z : float = 0):
        self.x = x
        self.y = y
        self.z = z

class Vec2():
    def __init__(self, u : float = 0, v : float = 0):
        self.u = u
        self.v = v

class Orientation(Enum):
    DOWN = -1
    NONE = 0
    UP = 1
    NOT_SET = 2
    ANY = 3
    DEFAULT = 4

class Notification(Enum):
    NONE = 0
    BRIEF = 1
    FULL = 2

class WarningLevel(Enum):
    IGNORE = 0
    WARNING = 1
    ERROR = 2

class VizMode(Enum):
    ASCII = 0
    CELLBLENDER_V1 = 1
    CELLBLENDER = 2

class Shape(Enum):
    UNSET = 0
    SPHERICAL = 1
    REGION_EXPR = 2
    LIST = 3
    COMPARTMENT = 4

class SurfacePropertyType(Enum):
    UNSET = 0
    REACTIVE = 1
    REFLECTIVE = 2
    TRANSPARENT = 3
    ABSORPTIVE = 4
    CONCENTRATION_CLAMP = 5
    FLUX_CLAMP = 6

class ExprNodeType(Enum):
    UNSET = 0
    LEAF = 1
    ADD = 2
    SUB = 3

class RegionNodeType(Enum):
    UNSET = 0
    LEAF_GEOMETRY_OBJECT = 1
    LEAF_SURFACE_REGION = 2
    UNION = 3
    DIFFERENCE = 4
    INTERSECT = 5

class ReactionType(Enum):
    UNSET = 0
    UNIMOL_VOLUME = 1
    UNIMOL_SURFACE = 2
    VOLUME_VOLUME = 3
    VOLUME_SURFACE = 4
    SURFACE_SURFACE = 5

class MoleculeType(Enum):
    UNSET = 0
    VOLUME = 1
    SURFACE = 2

class BNGSimulationMethod(Enum):
    NONE = 0
    ODE = 1
    SSA = 2
    PLA = 3
    NF = 4

class CountOutputFormat(Enum):
    UNSET = 0
    AUTOMATIC_FROM_EXTENSION = 1
    DAT = 2
    GDAT = 3



STATE_UNSET = 'STATE_UNSET'
STATE_UNSET_INT = -1
BOND_UNBOUND = -1
BOND_BOUND = -2
BOND_ANY = -3
PARTITION_EDGE_EXTRA_MARGIN_UM = 0.01
DEFAULT_COUNT_BUFFER_SIZE = 100
ALL_MOLECULES = 'ALL_MOLECULES'
ALL_VOLUME_MOLECULES = 'ALL_VOLUME_MOLECULES'
ALL_SURFACE_MOLECULES = 'ALL_SURFACE_MOLECULES'
DEFAULT_CHECKPOINTS_DIR = 'checkpoints'
DEFAULT_SEED_DIR_PREFIX = 'seed_'
DEFAULT_SEED_DIR_DIGITS = 5
DEFAULT_ITERATION_DIR_PREFIX = 'it_'
ID_INVALID = -1
NUMBER_OF_TRAINS_UNLIMITED = -1
TIME_INFINITY = 1e140
INT_UNSET = INT32_MAX
FLT_UNSET = FLT_MAX
RNG_SIZE = 256


class BaseChkptMol():
    def __init__(
            self,
            id : int,
            species : Species,
            diffusion_time : float,
            birthday : float,
            flags : int,
            unimol_rx_time : float = None
        ):
        self.id = id
        self.species = species
        self.diffusion_time = diffusion_time
        self.birthday = birthday
        self.flags = flags
        self.unimol_rx_time = unimol_rx_time


class ChkptSurfMol():
    def __init__(
            self,
            pos : Vec2,
            orientation : Orientation,
            geometry_object : GeometryObject,
            wall_index : int,
            grid_tile_index : int,
            id : int,
            species : Species,
            diffusion_time : float,
            birthday : float,
            flags : int,
            unimol_rx_time : float = None
        ):
        self.pos = pos
        self.orientation = orientation
        self.geometry_object = geometry_object
        self.wall_index = wall_index
        self.grid_tile_index = grid_tile_index
        self.id = id
        self.species = species
        self.diffusion_time = diffusion_time
        self.birthday = birthday
        self.flags = flags
        self.unimol_rx_time = unimol_rx_time


class ChkptVolMol():
    def __init__(
            self,
            pos : Vec3,
            id : int,
            species : Species,
            diffusion_time : float,
            birthday : float,
            flags : int,
            unimol_rx_time : float = None
        ):
        self.pos = pos
        self.id = id
        self.species = species
        self.diffusion_time = diffusion_time
        self.birthday = birthday
        self.flags = flags
        self.unimol_rx_time = unimol_rx_time


class Color():
    def __init__(
            self,
            red : float = None,
            green : float = None,
            blue : float = None,
            alpha : float = 1,
            rgba : int = 0
        ):
        self.red = red
        self.green = green
        self.blue = blue
        self.alpha = alpha
        self.rgba = rgba


class Complex():
    def __init__(
            self,
            name : str = None,
            elementary_molecules : List[ElementaryMolecule] = None,
            orientation : Orientation = Orientation.DEFAULT,
            compartment_name : str = None
        ):
        self.name = name
        self.elementary_molecules = elementary_molecules
        self.orientation = orientation
        self.compartment_name = compartment_name


    def to_bngl_str(
            self,
        ) -> 'str':
        pass

    def as_species(
            self,
        ) -> 'Species':
        pass

class Component():
    def __init__(
            self,
            component_type : ComponentType,
            state : str = STATE_UNSET,
            bond : int = BOND_UNBOUND
        ):
        self.component_type = component_type
        self.state = state
        self.bond = bond


    def to_bngl_str(
            self,
        ) -> 'str':
        pass

class ComponentType():
    def __init__(
            self,
            name : str,
            states : List[str] = None
        ):
        self.name = name
        self.states = states


    def inst(
            self,
            state : str = STATE_UNSET,
            bond : int = BOND_UNBOUND
        ) -> 'Component':
        pass

    def to_bngl_str(
            self,
        ) -> 'str':
        pass

class Config():
    def __init__(
            self,
            seed : int = 1,
            time_step : float = 1e-6,
            use_bng_units : bool = False,
            surface_grid_density : float = 10000,
            interaction_radius : float = None,
            intermembrane_interaction_radius : float = None,
            vacancy_search_distance : float = 10,
            center_molecules_on_grid : bool = False,
            partition_dimension : float = 10,
            initial_partition_origin : List[float] = None,
            subpartition_dimension : float = 0.5,
            total_iterations : float = 1000000,
            check_overlapped_walls : bool = True,
            reaction_class_cleanup_periodicity : int = 500,
            species_cleanup_periodicity : int = 10000,
            molecules_order_random_shuffle_periodicity : int = 10000,
            sort_molecules : bool = False,
            memory_limit_gb : int = -1,
            initial_iteration : int = 0,
            initial_time : float = 0,
            initial_rng_state : RngState = None,
            append_to_count_output_data : bool = False,
            continue_after_sigalrm : bool = False
        ):
        self.seed = seed
        self.time_step = time_step
        self.use_bng_units = use_bng_units
        self.surface_grid_density = surface_grid_density
        self.interaction_radius = interaction_radius
        self.intermembrane_interaction_radius = intermembrane_interaction_radius
        self.vacancy_search_distance = vacancy_search_distance
        self.center_molecules_on_grid = center_molecules_on_grid
        self.partition_dimension = partition_dimension
        self.initial_partition_origin = initial_partition_origin
        self.subpartition_dimension = subpartition_dimension
        self.total_iterations = total_iterations
        self.check_overlapped_walls = check_overlapped_walls
        self.reaction_class_cleanup_periodicity = reaction_class_cleanup_periodicity
        self.species_cleanup_periodicity = species_cleanup_periodicity
        self.molecules_order_random_shuffle_periodicity = molecules_order_random_shuffle_periodicity
        self.sort_molecules = sort_molecules
        self.memory_limit_gb = memory_limit_gb
        self.initial_iteration = initial_iteration
        self.initial_time = initial_time
        self.initial_rng_state = initial_rng_state
        self.append_to_count_output_data = append_to_count_output_data
        self.continue_after_sigalrm = continue_after_sigalrm


class Count():
    def __init__(
            self,
            name : str = None,
            file_name : str = None,
            expression : CountTerm = None,
            multiplier : float = 1,
            every_n_timesteps : float = 1,
            output_format : CountOutputFormat = CountOutputFormat.AUTOMATIC_FROM_EXTENSION
        ):
        self.name = name
        self.file_name = file_name
        self.expression = expression
        self.multiplier = multiplier
        self.every_n_timesteps = every_n_timesteps
        self.output_format = output_format


    def get_current_value(
            self,
        ) -> 'float':
        pass

class CountTerm():
    def __init__(
            self,
            species_pattern : Complex = None,
            molecules_pattern : Complex = None,
            reaction_rule : ReactionRule = None,
            region : Region = None,
            node_type : ExprNodeType = ExprNodeType.LEAF,
            left_node : CountTerm = None,
            right_node : CountTerm = None,
            initial_reactions_count : int = 0
        ):
        self.species_pattern = species_pattern
        self.molecules_pattern = molecules_pattern
        self.reaction_rule = reaction_rule
        self.region = region
        self.node_type = node_type
        self.left_node = left_node
        self.right_node = right_node
        self.initial_reactions_count = initial_reactions_count


    def __add__(
            self,
            op2 : CountTerm
        ) -> 'CountTerm':
        pass

    def __sub__(
            self,
            op2 : CountTerm
        ) -> 'CountTerm':
        pass

class ElementaryMolecule():
    def __init__(
            self,
            elementary_molecule_type : ElementaryMoleculeType,
            components : List[Component] = None,
            compartment_name : str = None
        ):
        self.elementary_molecule_type = elementary_molecule_type
        self.components = components
        self.compartment_name = compartment_name


    def to_bngl_str(
            self,
            with_compartment : bool = True
        ) -> 'str':
        pass

class ElementaryMoleculeType():
    def __init__(
            self,
            name : str,
            components : List[ComponentType] = None,
            diffusion_constant_2d : float = None,
            diffusion_constant_3d : float = None,
            custom_time_step : float = None,
            custom_space_step : float = None,
            target_only : bool = False
        ):
        self.name = name
        self.components = components
        self.diffusion_constant_2d = diffusion_constant_2d
        self.diffusion_constant_3d = diffusion_constant_3d
        self.custom_time_step = custom_time_step
        self.custom_space_step = custom_space_step
        self.target_only = target_only


    def inst(
            self,
            components : List[Component] = None,
            compartment_name : str = None
        ) -> 'ElementaryMolecule':
        pass

    def to_bngl_str(
            self,
        ) -> 'str':
        pass

class GeometryObject():
    def __init__(
            self,
            name : str,
            vertex_list : List[List[float]],
            wall_list : List[List[int]],
            is_bngl_compartment : bool = False,
            surface_compartment_name : str = None,
            surface_regions : List[SurfaceRegion] = None,
            surface_class : SurfaceClass = None,
            initial_surface_releases : List[InitialSurfaceRelease] = None,
            initial_color : Color = None,
            node_type : RegionNodeType = RegionNodeType.UNSET,
            left_node : Region = None,
            right_node : Region = None
        ):
        self.name = name
        self.vertex_list = vertex_list
        self.wall_list = wall_list
        self.is_bngl_compartment = is_bngl_compartment
        self.surface_compartment_name = surface_compartment_name
        self.surface_regions = surface_regions
        self.surface_class = surface_class
        self.initial_surface_releases = initial_surface_releases
        self.initial_color = initial_color
        self.node_type = node_type
        self.left_node = left_node
        self.right_node = right_node


    def translate(
            self,
            move : List[float]
        ) -> None:
        pass

    def __add__(
            self,
            other : Region
        ) -> 'Region':
        pass

    def __sub__(
            self,
            other : Region
        ) -> 'Region':
        pass

    def __mul__(
            self,
            other : Region
        ) -> 'Region':
        pass

class InitialSurfaceRelease():
    def __init__(
            self,
            complex : Complex,
            number_to_release : int = None,
            density : float = None
        ):
        self.complex = complex
        self.number_to_release = number_to_release
        self.density = density


class Instantiation():
    def __init__(
            self,
            release_sites : List[ReleaseSite] = None,
            geometry_objects : List[GeometryObject] = None,
            checkpointed_molecules : List[BaseChkptMol] = None
        ):
        self.release_sites = release_sites
        self.geometry_objects = geometry_objects
        self.checkpointed_molecules = checkpointed_molecules


    def add_release_site(
            self,
            s : ReleaseSite
        ) -> None:
        pass

    def find_release_site(
            self,
            name : str
        ) -> 'ReleaseSite':
        pass

    def add_geometry_object(
            self,
            o : GeometryObject
        ) -> None:
        pass

    def find_geometry_object(
            self,
            name : str
        ) -> 'GeometryObject':
        pass

    def find_volume_compartment_object(
            self,
            name : str
        ) -> 'GeometryObject':
        pass

    def find_surface_compartment_object(
            self,
            name : str
        ) -> 'GeometryObject':
        pass

    def load_bngl_compartments_and_seed_species(
            self,
            file_name : str,
            default_release_region : Region = None,
            parameter_overrides : Dict[str, float] = None
        ) -> None:
        pass

class Introspection():
    def __init__(
            self,
        ):
        pass

    def get_molecule_ids(
            self,
            pattern : Complex = None
        ) -> 'List[int]':
        pass

    def get_molecule(
            self,
            id : int
        ) -> 'Molecule':
        pass

    def get_species_name(
            self,
            species_id : int
        ) -> 'str':
        pass

    def get_vertex(
            self,
            object : GeometryObject,
            vertex_index : int
        ) -> 'Vec3':
        pass

    def get_wall(
            self,
            object : GeometryObject,
            wall_index : int
        ) -> 'Wall':
        pass

    def get_vertex_unit_normal(
            self,
            object : GeometryObject,
            vertex_index : int
        ) -> 'Vec3':
        pass

    def get_wall_unit_normal(
            self,
            object : GeometryObject,
            wall_index : int
        ) -> 'Vec3':
        pass

    def get_wall_color(
            self,
            object : GeometryObject,
            wall_index : int
        ) -> 'Color':
        pass

    def set_wall_color(
            self,
            object : GeometryObject,
            wall_index : int,
            color : Color
        ) -> None:
        pass

class Model():
    def __init__(
            self,
            config : Config = Config(),
            warnings : Warnings = Warnings(),
            notifications : Notifications = Notifications(),
            species : List[Species] = None,
            reaction_rules : List[ReactionRule] = None,
            surface_classes : List[SurfaceClass] = None,
            elementary_molecule_types : List[ElementaryMoleculeType] = None,
            release_sites : List[ReleaseSite] = None,
            geometry_objects : List[GeometryObject] = None,
            checkpointed_molecules : List[BaseChkptMol] = None,
            viz_outputs : List[VizOutput] = None,
            counts : List[Count] = None
        ):
        self.config = config
        self.warnings = warnings
        self.notifications = notifications
        self.species = species
        self.reaction_rules = reaction_rules
        self.surface_classes = surface_classes
        self.elementary_molecule_types = elementary_molecule_types
        self.release_sites = release_sites
        self.geometry_objects = geometry_objects
        self.checkpointed_molecules = checkpointed_molecules
        self.viz_outputs = viz_outputs
        self.counts = counts


    def initialize(
            self,
            print_copyright : bool = True
        ) -> None:
        pass

    def run_iterations(
            self,
            iterations : float
        ) -> 'int':
        pass

    def end_simulation(
            self,
            print_final_report : bool = True
        ) -> None:
        pass

    def add_subsystem(
            self,
            subsystem : Subsystem
        ) -> None:
        pass

    def add_instantiation(
            self,
            instantiation : Instantiation
        ) -> None:
        pass

    def add_observables(
            self,
            observables : Observables
        ) -> None:
        pass

    def dump_internal_state(
            self,
            with_geometry : bool = False
        ) -> None:
        pass

    def export_data_model(
            self,
            file : str = None
        ) -> None:
        pass

    def export_viz_data_model(
            self,
            file : str = None
        ) -> None:
        pass

    def export_geometry(
            self,
            output_files_prefix : str = None
        ) -> None:
        pass

    def release_molecules(
            self,
            release_site : ReleaseSite
        ) -> None:
        pass

    def run_reaction(
            self,
            reaction_rule : ReactionRule,
            reactant_ids : List[int],
            time : float
        ) -> 'List[int]':
        pass

    def add_vertex_move(
            self,
            object : GeometryObject,
            vertex_index : int,
            displacement : List[float]
        ) -> None:
        pass

    def apply_vertex_moves(
            self,
            collect_wall_wall_hits : bool = False
        ) -> 'List[WallWallHitInfo]':
        pass

    def register_mol_wall_hit_callback(
            self,
            function : Callable, # std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)>,
            context : Any, # py::object,
            object : GeometryObject = None,
            species : Species = None
        ) -> None:
        pass

    def register_reaction_callback(
            self,
            function : Callable, # std::function<void(std::shared_ptr<ReactionInfo>, py::object)>,
            context : Any, # py::object,
            reaction_rule : ReactionRule
        ) -> None:
        pass

    def load_bngl(
            self,
            file_name : str,
            observables_path_or_file : str = '',
            default_release_region : Region = None,
            parameter_overrides : Dict[str, float] = None,
            observables_output_format : CountOutputFormat = CountOutputFormat.AUTOMATIC_FROM_EXTENSION
        ) -> None:
        pass

    def export_to_bngl(
            self,
            file_name : str,
            simulation_method : BNGSimulationMethod = BNGSimulationMethod.ODE
        ) -> None:
        pass

    def save_checkpoint(
            self,
            custom_dir : str = None
        ) -> None:
        pass

    def schedule_checkpoint(
            self,
            iteration : int = 0,
            continue_simulation : bool = False,
            custom_dir : str = None
        ) -> None:
        pass

    def add_species(
            self,
            s : Species
        ) -> None:
        pass

    def find_species(
            self,
            name : str
        ) -> 'Species':
        pass

    def add_reaction_rule(
            self,
            r : ReactionRule
        ) -> None:
        pass

    def find_reaction_rule(
            self,
            name : str
        ) -> 'ReactionRule':
        pass

    def add_surface_class(
            self,
            sc : SurfaceClass
        ) -> None:
        pass

    def find_surface_class(
            self,
            name : str
        ) -> 'SurfaceClass':
        pass

    def add_elementary_molecule_type(
            self,
            mt : ElementaryMoleculeType
        ) -> None:
        pass

    def find_elementary_molecule_type(
            self,
            name : str
        ) -> 'ElementaryMoleculeType':
        pass

    def load_bngl_molecule_types_and_reaction_rules(
            self,
            file_name : str,
            parameter_overrides : Dict[str, float] = None
        ) -> None:
        pass

    def add_release_site(
            self,
            s : ReleaseSite
        ) -> None:
        pass

    def find_release_site(
            self,
            name : str
        ) -> 'ReleaseSite':
        pass

    def add_geometry_object(
            self,
            o : GeometryObject
        ) -> None:
        pass

    def find_geometry_object(
            self,
            name : str
        ) -> 'GeometryObject':
        pass

    def find_volume_compartment_object(
            self,
            name : str
        ) -> 'GeometryObject':
        pass

    def find_surface_compartment_object(
            self,
            name : str
        ) -> 'GeometryObject':
        pass

    def load_bngl_compartments_and_seed_species(
            self,
            file_name : str,
            default_release_region : Region = None,
            parameter_overrides : Dict[str, float] = None
        ) -> None:
        pass

    def add_viz_output(
            self,
            viz_output : VizOutput
        ) -> None:
        pass

    def add_count(
            self,
            count : Count
        ) -> None:
        pass

    def find_count(
            self,
            name : str
        ) -> 'Count':
        pass

    def load_bngl_observables(
            self,
            file_name : str,
            observables_path_or_file : str = '',
            parameter_overrides : Dict[str, float] = None,
            observables_output_format : CountOutputFormat = CountOutputFormat.AUTOMATIC_FROM_EXTENSION
        ) -> None:
        pass

    def get_molecule_ids(
            self,
            pattern : Complex = None
        ) -> 'List[int]':
        pass

    def get_molecule(
            self,
            id : int
        ) -> 'Molecule':
        pass

    def get_species_name(
            self,
            species_id : int
        ) -> 'str':
        pass

    def get_vertex(
            self,
            object : GeometryObject,
            vertex_index : int
        ) -> 'Vec3':
        pass

    def get_wall(
            self,
            object : GeometryObject,
            wall_index : int
        ) -> 'Wall':
        pass

    def get_vertex_unit_normal(
            self,
            object : GeometryObject,
            vertex_index : int
        ) -> 'Vec3':
        pass

    def get_wall_unit_normal(
            self,
            object : GeometryObject,
            wall_index : int
        ) -> 'Vec3':
        pass

    def get_wall_color(
            self,
            object : GeometryObject,
            wall_index : int
        ) -> 'Color':
        pass

    def set_wall_color(
            self,
            object : GeometryObject,
            wall_index : int,
            color : Color
        ) -> None:
        pass

class MolWallHitInfo():
    def __init__(
            self,
            molecule_id : int,
            geometry_object : GeometryObject,
            wall_index : int,
            time : float,
            pos3d : Vec3,
            time_before_hit : float,
            pos3d_before_hit : Vec3
        ):
        self.molecule_id = molecule_id
        self.geometry_object = geometry_object
        self.wall_index = wall_index
        self.time = time
        self.pos3d = pos3d
        self.time_before_hit = time_before_hit
        self.pos3d_before_hit = pos3d_before_hit


class Molecule():
    def __init__(
            self,
            id : int = ID_INVALID,
            type : MoleculeType = MoleculeType.UNSET,
            species_id : int = ID_INVALID,
            pos3d : Vec3 = None,
            orientation : Orientation = Orientation.NOT_SET,
            pos2d : Vec2 = None,
            geometry_object : GeometryObject = None,
            wall_index : int = -1
        ):
        self.id = id
        self.type = type
        self.species_id = species_id
        self.pos3d = pos3d
        self.orientation = orientation
        self.pos2d = pos2d
        self.geometry_object = geometry_object
        self.wall_index = wall_index


    def remove(
            self,
        ) -> None:
        pass

class MoleculeReleaseInfo():
    def __init__(
            self,
            complex : Complex,
            location : List[float]
        ):
        self.complex = complex
        self.location = location


class Notifications():
    def __init__(
            self,
            bng_verbosity_level : int = 0,
            rxn_and_species_report : bool = False,
            simulation_stats_every_n_iterations : int = 0,
            rxn_probability_changed : bool = True,
            iteration_report : bool = True
        ):
        self.bng_verbosity_level = bng_verbosity_level
        self.rxn_and_species_report = rxn_and_species_report
        self.simulation_stats_every_n_iterations = simulation_stats_every_n_iterations
        self.rxn_probability_changed = rxn_probability_changed
        self.iteration_report = iteration_report


class Observables():
    def __init__(
            self,
            viz_outputs : List[VizOutput] = None,
            counts : List[Count] = None
        ):
        self.viz_outputs = viz_outputs
        self.counts = counts


    def add_viz_output(
            self,
            viz_output : VizOutput
        ) -> None:
        pass

    def add_count(
            self,
            count : Count
        ) -> None:
        pass

    def find_count(
            self,
            name : str
        ) -> 'Count':
        pass

    def load_bngl_observables(
            self,
            file_name : str,
            observables_path_or_file : str = '',
            parameter_overrides : Dict[str, float] = None,
            observables_output_format : CountOutputFormat = CountOutputFormat.AUTOMATIC_FROM_EXTENSION
        ) -> None:
        pass

class ReactionInfo():
    def __init__(
            self,
            type : ReactionType,
            reactant_ids : List[int],
            product_ids : List[int],
            reaction_rule : ReactionRule,
            time : float,
            pos3d : Vec3,
            geometry_object : GeometryObject = None,
            wall_index : int = -1,
            pos2d : Vec2 = None
        ):
        self.type = type
        self.reactant_ids = reactant_ids
        self.product_ids = product_ids
        self.reaction_rule = reaction_rule
        self.time = time
        self.pos3d = pos3d
        self.geometry_object = geometry_object
        self.wall_index = wall_index
        self.pos2d = pos2d


class ReactionRule():
    def __init__(
            self,
            name : str = None,
            reactants : List[Complex] = None,
            products : List[Complex] = None,
            fwd_rate : float = None,
            rev_name : str = None,
            rev_rate : float = None,
            variable_rate : List[List[float]] = None,
            is_intermembrane_surface_reaction : bool = False
        ):
        self.name = name
        self.reactants = reactants
        self.products = products
        self.fwd_rate = fwd_rate
        self.rev_name = rev_name
        self.rev_rate = rev_rate
        self.variable_rate = variable_rate
        self.is_intermembrane_surface_reaction = is_intermembrane_surface_reaction


    def to_bngl_str(
            self,
        ) -> 'str':
        pass

class Region():
    def __init__(
            self,
            node_type : RegionNodeType = RegionNodeType.UNSET,
            left_node : Region = None,
            right_node : Region = None
        ):
        self.node_type = node_type
        self.left_node = left_node
        self.right_node = right_node


    def __add__(
            self,
            other : Region
        ) -> 'Region':
        pass

    def __sub__(
            self,
            other : Region
        ) -> 'Region':
        pass

    def __mul__(
            self,
            other : Region
        ) -> 'Region':
        pass

class ReleasePattern():
    def __init__(
            self,
            name : str = None,
            release_interval : float = TIME_INFINITY,
            train_duration : float = TIME_INFINITY,
            train_interval : float = TIME_INFINITY,
            number_of_trains : int = 1
        ):
        self.name = name
        self.release_interval = release_interval
        self.train_duration = train_duration
        self.train_interval = train_interval
        self.number_of_trains = number_of_trains


class ReleaseSite():
    def __init__(
            self,
            name : str,
            complex : Complex = None,
            molecule_list : List[MoleculeReleaseInfo] = None,
            release_time : float = 0,
            release_pattern : ReleasePattern = None,
            shape : Shape = Shape.UNSET,
            region : Region = None,
            location : List[float] = None,
            site_diameter : float = 0,
            site_radius : float = None,
            number_to_release : float = None,
            density : float = None,
            concentration : float = None,
            release_probability : float = 1
        ):
        self.name = name
        self.complex = complex
        self.molecule_list = molecule_list
        self.release_time = release_time
        self.release_pattern = release_pattern
        self.shape = shape
        self.region = region
        self.location = location
        self.site_diameter = site_diameter
        self.site_radius = site_radius
        self.number_to_release = number_to_release
        self.density = density
        self.concentration = concentration
        self.release_probability = release_probability


class RngState():
    def __init__(
            self,
            randcnt : int,
            aa : int,
            bb : int,
            cc : int,
            randslr : List[int],
            mm : List[int],
            rngblocks : int
        ):
        self.randcnt = randcnt
        self.aa = aa
        self.bb = bb
        self.cc = cc
        self.randslr = randslr
        self.mm = mm
        self.rngblocks = rngblocks


class Species():
    def __init__(
            self,
            name : str = None,
            diffusion_constant_2d : float = None,
            diffusion_constant_3d : float = None,
            custom_time_step : float = None,
            custom_space_step : float = None,
            target_only : bool = False,
            elementary_molecules : List[ElementaryMolecule] = None,
            orientation : Orientation = Orientation.DEFAULT,
            compartment_name : str = None
        ):
        self.name = name
        self.diffusion_constant_2d = diffusion_constant_2d
        self.diffusion_constant_3d = diffusion_constant_3d
        self.custom_time_step = custom_time_step
        self.custom_space_step = custom_space_step
        self.target_only = target_only
        self.name = name
        self.elementary_molecules = elementary_molecules
        self.orientation = orientation
        self.compartment_name = compartment_name


    def inst(
            self,
            orientation : Orientation = Orientation.DEFAULT,
            compartment_name : str = None
        ) -> 'Complex':
        pass

    def to_bngl_str(
            self,
        ) -> 'str':
        pass

    def as_species(
            self,
        ) -> 'Species':
        pass

class Subsystem():
    def __init__(
            self,
            species : List[Species] = None,
            reaction_rules : List[ReactionRule] = None,
            surface_classes : List[SurfaceClass] = None,
            elementary_molecule_types : List[ElementaryMoleculeType] = None
        ):
        self.species = species
        self.reaction_rules = reaction_rules
        self.surface_classes = surface_classes
        self.elementary_molecule_types = elementary_molecule_types


    def add_species(
            self,
            s : Species
        ) -> None:
        pass

    def find_species(
            self,
            name : str
        ) -> 'Species':
        pass

    def add_reaction_rule(
            self,
            r : ReactionRule
        ) -> None:
        pass

    def find_reaction_rule(
            self,
            name : str
        ) -> 'ReactionRule':
        pass

    def add_surface_class(
            self,
            sc : SurfaceClass
        ) -> None:
        pass

    def find_surface_class(
            self,
            name : str
        ) -> 'SurfaceClass':
        pass

    def add_elementary_molecule_type(
            self,
            mt : ElementaryMoleculeType
        ) -> None:
        pass

    def find_elementary_molecule_type(
            self,
            name : str
        ) -> 'ElementaryMoleculeType':
        pass

    def load_bngl_molecule_types_and_reaction_rules(
            self,
            file_name : str,
            parameter_overrides : Dict[str, float] = None
        ) -> None:
        pass

class SurfaceClass():
    def __init__(
            self,
            name : str,
            properties : List[SurfaceProperty] = None,
            type : SurfacePropertyType = SurfacePropertyType.UNSET,
            affected_complex_pattern : Complex = None,
            concentration : float = None
        ):
        self.name = name
        self.properties = properties
        self.type = type
        self.affected_complex_pattern = affected_complex_pattern
        self.concentration = concentration


class SurfaceProperty():
    def __init__(
            self,
            type : SurfacePropertyType = SurfacePropertyType.UNSET,
            affected_complex_pattern : Complex = None,
            concentration : float = None
        ):
        self.type = type
        self.affected_complex_pattern = affected_complex_pattern
        self.concentration = concentration


class SurfaceRegion():
    def __init__(
            self,
            name : str,
            wall_indices : List[int],
            surface_class : SurfaceClass = None,
            initial_surface_releases : List[InitialSurfaceRelease] = None,
            initial_color : Color = None,
            node_type : RegionNodeType = RegionNodeType.UNSET,
            left_node : Region = None,
            right_node : Region = None
        ):
        self.name = name
        self.wall_indices = wall_indices
        self.surface_class = surface_class
        self.initial_surface_releases = initial_surface_releases
        self.initial_color = initial_color
        self.node_type = node_type
        self.left_node = left_node
        self.right_node = right_node


    def __add__(
            self,
            other : Region
        ) -> 'Region':
        pass

    def __sub__(
            self,
            other : Region
        ) -> 'Region':
        pass

    def __mul__(
            self,
            other : Region
        ) -> 'Region':
        pass

class VizOutput():
    def __init__(
            self,
            output_files_prefix : str,
            species_list : List[Species] = None,
            mode : VizMode = VizMode.ASCII,
            every_n_timesteps : float = 1
        ):
        self.output_files_prefix = output_files_prefix
        self.species_list = species_list
        self.mode = mode
        self.every_n_timesteps = every_n_timesteps


class Wall():
    def __init__(
            self,
            geometry_object : GeometryObject,
            wall_index : int,
            vertices : List[Vec3],
            area : float,
            unit_normal : Vec3,
            is_movable : bool = True
        ):
        self.geometry_object = geometry_object
        self.wall_index = wall_index
        self.vertices = vertices
        self.area = area
        self.unit_normal = unit_normal
        self.is_movable = is_movable


class WallWallHitInfo():
    def __init__(
            self,
            wall1 : Wall,
            wall2 : Wall
        ):
        self.wall1 = wall1
        self.wall2 = wall2


class Warnings():
    def __init__(
            self,
            high_reaction_probability : WarningLevel = WarningLevel.IGNORE,
            molecule_placement_failure : WarningLevel = WarningLevel.ERROR
        ):
        self.high_reaction_probability = high_reaction_probability
        self.molecule_placement_failure = molecule_placement_failure


class bngl_utils():
    def __init__(
            self,
        ):
        pass

    def load_bngl_parameters(
            self,
            file_name : str,
            parameter_overrides : Dict[str, float] = None
        ) -> 'Dict[str, float]':
        pass

class data_utils():
    def __init__(
            self,
        ):
        pass

    def load_dat_file(
            self,
            file_name : str
        ) -> 'List[List[float]]':
        pass

class geometry_utils():
    def __init__(
            self,
        ):
        pass

    def create_box(
            self,
            name : str,
            edge_dimension : float = None,
            xyz_dimensions : List[float] = None
        ) -> 'GeometryObject':
        pass

    def create_icosphere(
            self,
            name : str,
            radius : float,
            subdivisions : int
        ) -> 'GeometryObject':
        pass

class run_utils():
    def __init__(
            self,
        ):
        pass

    def get_last_checkpoint_dir(
            self,
            seed : int
        ) -> 'str':
        pass

    def remove_cwd(
            self,
            paths : List[str]
        ) -> 'List[str]':
        pass

AllMolecules = Species('ALL_MOLECULES')
AllVolumeMolecules = Species('ALL_VOLUME_MOLECULES')
AllSurfaceMolecules = Species('ALL_SURFACE_MOLECULES')
