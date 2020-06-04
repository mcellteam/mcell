from typing import List
from enum import Enum

class Vec3():
    def __init__(self, x : float = 0, y : float = 0, z : float = 0):
        self.x = x
        self.y = y
        self.z = z

class Orientation(Enum):
    DOWN = -1
    NONE = 0
    UP = 1
    NOT_SET = 2
    ANY = 3

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
    CELLBLENDER = 1

class Shape(Enum):
    UNSET = 0
    SPHERICAL = 1
    REGION_EXPR = 2

class SurfacePropertyType(Enum):
    UNSET = 0
    REFLECTIVE = 1
    TRANSPARENT = 2
    ABSORPTIVE = 3

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



STATE_UNSET = 'state_unset'
STATE_UNSET_INT = -1
BOND_UNBOUND = 0
BOND_BOUND = -1
PARTITION_EDGE_EXTRA_MARGIN_UM = 0.01
DEFAULT_COUNT_BUFFER_SIZE = 10000
ALL_MOLECULES = 'ALL_MOLECULES'
ALL_VOLUME_MOLECULES = 'ALL_VOLUME_MOLECULES'
ALL_SURFACE_MOLECULES = 'ALL_SURFACE_MOLECULES'


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
        ) -> 'ComponentInstance':
        pass

    def inst(
            self,
            state : int = STATE_UNSET_INT,
            bond : int = BOND_UNBOUND
        ) -> 'ComponentInstance':
        pass

class ComponentInstance():
    def __init__(
            self,
            component_type : ComponentType,
            state : str = STATE_UNSET,
            bond : int = BOND_UNBOUND
        ):
        self.component_type = component_type
        self.state = state
        self.bond = bond


class ElementaryMoleculeType():
    def __init__(
            self,
            name : str,
            components : List[ComponentType] = None,
            diffusion_constant_2d : float = None,
            diffusion_constant_3d : float = None
        ):
        self.name = name
        self.components = components
        self.diffusion_constant_2d = diffusion_constant_2d
        self.diffusion_constant_3d = diffusion_constant_3d


    def inst(
            self,
            components : List[ComponentInstance]
        ) -> 'ElementaryMoleculeInstance':
        pass

class ElementaryMoleculeInstance():
    def __init__(
            self,
            elementary_molecule_type : ElementaryMoleculeType,
            components : List[ComponentInstance] = None
        ):
        self.elementary_molecule_type = elementary_molecule_type
        self.components = components


class ComplexInstance():
    def __init__(
            self,
            molecule_instances : List[ElementaryMoleculeInstance] = None,
            orientation : Orientation = Orientation.NONE
        ):
        self.molecule_instances = molecule_instances
        self.orientation = orientation


class Species():
    def __init__(
            self,
            name : str,
            diffusion_constant_2d : float = None,
            diffusion_constant_3d : float = None,
            molecule_instances : List[ElementaryMoleculeInstance] = None,
            orientation : Orientation = Orientation.NONE
        ):
        self.name = name
        self.diffusion_constant_2d = diffusion_constant_2d
        self.diffusion_constant_3d = diffusion_constant_3d
        self.molecule_instances = molecule_instances
        self.orientation = orientation


    def inst(
            self,
            orientation : Orientation = Orientation.NOT_SET
        ) -> 'ComplexInstance':
        pass

class SurfaceProperty():
    def __init__(
            self,
            type : SurfacePropertyType = SurfacePropertyType.UNSET,
            affected_species : Species = None,
            orientation : Orientation = Orientation.NOT_SET
        ):
        self.type = type
        self.affected_species = affected_species
        self.orientation = orientation


class SurfaceClass():
    def __init__(
            self,
            name : str,
            properties : List[SurfaceProperty] = None,
            type : SurfacePropertyType = SurfacePropertyType.UNSET,
            affected_species : Species = None,
            orientation : Orientation = Orientation.NOT_SET
        ):
        self.name = name
        self.properties = properties
        self.type = type
        self.affected_species = affected_species
        self.orientation = orientation


class ReactionRule():
    def __init__(
            self,
            name : str = None,
            reactants : List[ComplexInstance] = None,
            products : List[ComplexInstance] = None,
            fwd_rate : float = None,
            rev_name : str = None,
            rev_rate : float = None,
            variable_rate : List[List[float]] = None
        ):
        self.name = name
        self.reactants = reactants
        self.products = products
        self.fwd_rate = fwd_rate
        self.rev_name = rev_name
        self.rev_rate = rev_rate
        self.variable_rate = variable_rate


class Subsystem():
    def __init__(
            self,
            species : List[Species] = None,
            reaction_rules : List[ReactionRule] = None,
            surface_classes : List[SurfaceClass] = None
        ):
        self.species = species
        self.reaction_rules = reaction_rules
        self.surface_classes = surface_classes


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

class Region():
    def __init__(
            self,
            node_type : RegionNodeType = RegionNodeType.UNSET,
            left_node : 'Region' = None,
            right_node : 'Region' = None
        ):
        self.node_type = node_type
        self.left_node = left_node
        self.right_node = right_node


    def __add__(
            self,
            other : 'Region'
        ) -> 'Region':
        pass

    def __sub__(
            self,
            other : 'Region'
        ) -> 'Region':
        pass

    def __mul__(
            self,
            other : 'Region'
        ) -> 'Region':
        pass

class SurfaceRegion():
    def __init__(
            self,
            name : str,
            wall_indices : List[int],
            surface_class : SurfaceClass = None,
            node_type : RegionNodeType = RegionNodeType.UNSET,
            left_node : Region = None,
            right_node : Region = None
        ):
        self.name = name
        self.wall_indices = wall_indices
        self.surface_class = surface_class
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

class GeometryObject():
    def __init__(
            self,
            name : str,
            vertex_list : List[List[float]],
            element_connections : List[List[int]],
            surface_regions : List[SurfaceRegion] = None,
            surface_class : SurfaceClass = None,
            node_type : RegionNodeType = RegionNodeType.UNSET,
            left_node : Region = None,
            right_node : Region = None
        ):
        self.name = name
        self.vertex_list = vertex_list
        self.element_connections = element_connections
        self.surface_regions = surface_regions
        self.surface_class = surface_class
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

class ReleaseSite():
    def __init__(
            self,
            name : str,
            species : Species,
            orientation : Orientation = Orientation.NONE,
            shape : Shape = Shape.UNSET,
            region : Region = None,
            location : Vec3 = None,
            site_diameter : float = 0,
            site_radius : float = None,
            number_to_release : int = None,
            density : float = None,
            release_probability : float = None
        ):
        self.name = name
        self.species = species
        self.orientation = orientation
        self.shape = shape
        self.region = region
        self.location = location
        self.site_diameter = site_diameter
        self.site_radius = site_radius
        self.number_to_release = number_to_release
        self.density = density
        self.release_probability = release_probability


class InstantiationData():
    def __init__(
            self,
            release_sites : List[ReleaseSite] = None,
            geometry_objects : List[GeometryObject] = None
        ):
        self.release_sites = release_sites
        self.geometry_objects = geometry_objects


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

class VizOutput():
    def __init__(
            self,
            filename_prefix : str,
            species_list : List[Species] = None,
            mode : VizMode = VizMode.ASCII,
            every_n_timesteps : int = 1
        ):
        self.filename_prefix = filename_prefix
        self.species_list = species_list
        self.mode = mode
        self.every_n_timesteps = every_n_timesteps


class CountTerm():
    def __init__(
            self,
            species : Species = None,
            reaction_rule : ReactionRule = None,
            region : Region = None,
            orientation : Orientation = Orientation.NOT_SET,
            node_type : ExprNodeType = ExprNodeType.LEAF,
            left_node : 'CountTerm' = None,
            right_node : 'CountTerm' = None
        ):
        self.species = species
        self.reaction_rule = reaction_rule
        self.region = region
        self.orientation = orientation
        self.node_type = node_type
        self.left_node = left_node
        self.right_node = right_node


    def __add__(
            self,
            op2 : 'CountTerm'
        ) -> 'CountTerm':
        pass

    def __sub__(
            self,
            op2 : 'CountTerm'
        ) -> 'CountTerm':
        pass

class Count():
    def __init__(
            self,
            filename : str,
            count_expression : CountTerm = None,
            every_n_timesteps : int = 1,
            species : Species = None,
            reaction_rule : ReactionRule = None,
            region : Region = None,
            orientation : Orientation = Orientation.NOT_SET,
            node_type : ExprNodeType = ExprNodeType.LEAF,
            left_node : CountTerm = None,
            right_node : CountTerm = None
        ):
        self.filename = filename
        self.count_expression = count_expression
        self.every_n_timesteps = every_n_timesteps
        self.species = species
        self.reaction_rule = reaction_rule
        self.region = region
        self.orientation = orientation
        self.node_type = node_type
        self.left_node = left_node
        self.right_node = right_node


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

class Config():
    def __init__(
            self,
            seed : int = 1,
            time_step : float = 1e-6,
            surface_grid_density : float = 10000,
            interaction_radius : float = None,
            vacancy_search_distance : float = 10,
            center_molecules_on_grid : bool = False,
            partition_dimension : float = 10,
            subpartition_dimension : float = 0.5,
            total_iterations_hint : int = 1000000
        ):
        self.seed = seed
        self.time_step = time_step
        self.surface_grid_density = surface_grid_density
        self.interaction_radius = interaction_radius
        self.vacancy_search_distance = vacancy_search_distance
        self.center_molecules_on_grid = center_molecules_on_grid
        self.partition_dimension = partition_dimension
        self.subpartition_dimension = subpartition_dimension
        self.total_iterations_hint = total_iterations_hint


class Notifications():
    def __init__(
            self,
            probability_report : bool = True,
            diffusion_constant_report : Notification = Notification.BRIEF,
            final_summary : bool = True,
            iteration_report : bool = True,
            varying_probability_report : bool = True,
            progress_report : bool = True,
            release_event_report : bool = True,
            molecule_collision_report : bool = True
        ):
        self.probability_report = probability_report
        self.diffusion_constant_report = diffusion_constant_report
        self.final_summary = final_summary
        self.iteration_report = iteration_report
        self.varying_probability_report = varying_probability_report
        self.progress_report = progress_report
        self.release_event_report = release_event_report
        self.molecule_collision_report = molecule_collision_report


class Warnings():
    def __init__(
            self,
            molecule_collision_report : WarningLevel = WarningLevel.WARNING,
            degenerate_polygons : WarningLevel = WarningLevel.WARNING,
            negative_diffusion_constant : WarningLevel = WarningLevel.WARNING,
            missing_surface_orientation : WarningLevel = WarningLevel.ERROR,
            negative_reaction_rate : WarningLevel = WarningLevel.WARNING,
            useless_volume_orientation : WarningLevel = WarningLevel.WARNING,
            high_reaction_probability : WarningLevel = WarningLevel.IGNORE,
            lifetime_too_short : WarningLevel = WarningLevel.WARNING,
            lifetime_threshold : float = 50,
            missed_reactions : WarningLevel = WarningLevel.WARNING,
            missed_reactions_threshold : float = 0.00100000004749745
        ):
        self.molecule_collision_report = molecule_collision_report
        self.degenerate_polygons = degenerate_polygons
        self.negative_diffusion_constant = negative_diffusion_constant
        self.missing_surface_orientation = missing_surface_orientation
        self.negative_reaction_rate = negative_reaction_rate
        self.useless_volume_orientation = useless_volume_orientation
        self.high_reaction_probability = high_reaction_probability
        self.lifetime_too_short = lifetime_too_short
        self.lifetime_threshold = lifetime_threshold
        self.missed_reactions = missed_reactions
        self.missed_reactions_threshold = missed_reactions_threshold


class Model():
    def __init__(
            self,
            config : Config = Config(),
            warnings : Warnings = Warnings(),
            notifications : Notifications = Notifications(),
            species : List[Species] = None,
            reaction_rules : List[ReactionRule] = None,
            surface_classes : List[SurfaceClass] = None,
            release_sites : List[ReleaseSite] = None,
            geometry_objects : List[GeometryObject] = None,
            viz_outputs : List[VizOutput] = None,
            counts : List[Count] = None
        ):
        self.config = config
        self.warnings = warnings
        self.notifications = notifications
        self.species = species
        self.reaction_rules = reaction_rules
        self.surface_classes = surface_classes
        self.release_sites = release_sites
        self.geometry_objects = geometry_objects
        self.viz_outputs = viz_outputs
        self.counts = counts


    def initialize(
            self,
        ) -> None:
        pass

    def run_iterations(
            self,
            iterations : int
        ) -> None:
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

    def add_instantiation_data(
            self,
            instantiation_data : InstantiationData
        ) -> None:
        pass

    def add_observables(
            self,
            observables : Observables
        ) -> None:
        pass

    def dump_internal_state(
            self,
        ) -> None:
        pass

    def export_data_model(
            self,
            file : str = None
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

AllMolecules = Species('ALL_MOLECULES')
AllVolumeMolecules = Species('ALL_VOLUME_MOLECULES')
AllSurfaceMolecules = Species('ALL_SURFACE_MOLECULES')
