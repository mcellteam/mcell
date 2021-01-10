**************************************
MCell 4 Python API Enums and Constants
**************************************

Orientation
===========


  | Orientation of a Complex.

* | **DOWN** = -1
* | **NONE** = 0
* | **UP** = 1
* | **NOT_SET** = 2
* | **ANY** = 3
* | **DEFAULT** = 4
  | DEFAULT means NONE for volume complexes and UP for surface complexes.


Notification
============

* | **NONE** = 0
* | **BRIEF** = 1
* | **FULL** = 2

WarningLevel
============

* | **IGNORE** = 0
  | Do something sensible and continue silently

* | **WARNING** = 1
  | Do something sensible but emit a warning message

* | **ERROR** = 2
  | Treat the warning as an error and stop


VizMode
=======

* | **ASCII** = 0
* | **CELLBLENDER** = 1

Shape
=====

* | **UNSET** = 0
* | **SPHERICAL** = 1
* | **REGION_EXPR** = 2
* | **LIST** = 3
* | **COMPARTMENT** = 4

SurfacePropertyType
===================

* | **UNSET** = 0
* | **REFLECTIVE** = 1
* | **TRANSPARENT** = 2
* | **ABSORPTIVE** = 3
* | **CONCENTRATION_CLAMP** = 4
  | Clamps concentration at a surface by periodically releasing molecules that correspond
  | to the wall being a transparent boundary to area with given concentration, 
  | and by absorbing all molecules that hit this surface.

* | **FLUX_CLAMP** = 5
  | Clamps flux at a surface by periodically releasing molecules that correspond
  | to the wall being a transparent boundary to area with given concentration. 
  | The clamped surface reflects these molecules.


ExprNodeType
============


  | Used internally to represent expression trees.

* | **UNSET** = 0
* | **LEAF** = 1
* | **ADD** = 2
* | **SUB** = 3

RegionNodeType
==============


  | Used internally to represent region trees.

* | **UNSET** = 0
* | **LEAF_GEOMETRY_OBJECT** = 1
* | **LEAF_SURFACE_REGION** = 2
* | **UNION** = 3
* | **DIFFERENCE** = 4
* | **INTERSECT** = 5

ReactionType
============


  | Used in reaction callbacks.

* | **UNSET** = 0
* | **UNIMOL_VOLUME** = 1
* | **UNIMOL_SURFACE** = 2
* | **VOLUME_VOLUME** = 3
* | **VOLUME_SURFACE** = 4
* | **SURFACE_SURFACE** = 5

MoleculeType
============


  | Used in molecule introspection and internally in checkpointing.

* | **UNSET** = 0
* | **VOLUME** = 1
* | **SURFACE** = 2



Constants
=========

* | **STATE_UNSET**: str = STATE_UNSET
* | **STATE_UNSET_INT**: int = -1
* | **BOND_UNBOUND**: int = -1
  | Represents cases when a component must not be bound in a pattern

* | **BOND_BOUND**: int = -2
  | Represents bond type !+ in a pattern

* | **BOND_ANY**: int = -3
  | Represents bond type !? in a pattern

* | **PARTITION_EDGE_EXTRA_MARGIN_UM**: float = 0.01
  | Internal constant used to match partition setup when comparing models against their MDL variant.

* | **DEFAULT_COUNT_BUFFER_SIZE**: int = 10000
  | Internal constant used to initialize buffer size for mol and rxn counts.

* | **ALL_MOLECULES**: str = ALL_MOLECULES
* | **ALL_VOLUME_MOLECULES**: str = ALL_VOLUME_MOLECULES
* | **ALL_SURFACE_MOLECULES**: str = ALL_SURFACE_MOLECULES
* | **DEFAULT_CHECKPOINTS_DIR**: str = checkpoints
* | **DEFAULT_SEED_DIR_PREFIX**: str = seed_
* | **DEFAULT_SEED_DIR_DIGITS**: int = 5
* | **DEFAULT_ITERATION_DIR_PREFIX**: str = it_
* | **AllMolecules**: Species = AllMolecules
* | **AllVolumeMolecules**: Species = AllVolumeMolecules
* | **AllSurfaceMolecules**: Species = AllSurfaceMolecules
* | **ID_INVALID**: int = -1
* | **NUMBER_OF_TRAINS_UNLIMITED**: int = -1
* | **TIME_INFINITY**: float = 1e140
* | **INT_UNSET**: int = INT32_MAX
  | This is a special integer value that means that an argument was not set, 
  | its value is 2147483647.

* | **FLT_UNSET**: float = FLT_MAX
  | This is a special floating point value that means that an argument was not set, 
  | its value is 3.40282346638528859812e+38F.

* | **RNG_SIZE**: int = 256
  | Size of arrays of



**************************
MCell 4 Python API Classes
**************************

BaseChkptMol
============

All times are in us (microseconds).

Attributes:
***********
* | **id**: int

* | **species**: Species

* | **diffusion_time**: float

* | **birthday**: float

* | **flags**: int

* | **unimol_rx_time**: float = None

ChkptSurfMol
============

Attributes:
***********
* | **pos**: Vec2

* | **orientation**: Orientation

* | **geometry_object**: GeometryObject

* | **wall_index**: int

* | **grid_tile_index**: int

* | **id**: int

* | **species**: Species

* | **diffusion_time**: float

* | **birthday**: float

* | **flags**: int

* | **unimol_rx_time**: float = None

ChkptVolMol
===========

Attributes:
***********
* | **pos**: Vec3

* | **id**: int

* | **species**: Species

* | **diffusion_time**: float

* | **birthday**: float

* | **flags**: int

* | **unimol_rx_time**: float = None

Complex
=======

This class represents a complex molecule composed of molecule instances.
It is either defined using a BNGL string or using a list of elementary molecule instances.
On top of that, orientation may be defined.
This class is used as argument in cases where either a fully qualified instance or a pattern 
can be provided such as in observable Count.  
Comparison operator __eq__ first converts complexes to their canonical representation and 
then does comparison so for instance m.Complex('A(b!1).B(a!1)') == m.Complex('B(a!2).A(b!2)').

Attributes:
***********
* | **name**: str = None
  | When set, this complex instance is initialized from a BNGL string passed as this argument, 
  | the string is parsed during model initialization so the molecule types it uses
  | don't have to be defined before initialization.

* | **elementary_molecules**: List[ElementaryMolecule] = None
  | Individual molecule instances contained in the complex.

* | **orientation**: Orientation = Orientation.DEFAULT
  | Specifies orientation of a molecule. 
  | When Orientation.DEFAULT if kept then during model initialization is
  | 'orientation' set to Orientation.NONE for volume complexes and to 
  | Orientation.UP for surface complexes.
  | Ignored by derived class Species.

* | **compartment_name**: str = None
  | Specifies compartment name of this Complex. 
  | Only one of 'orientation' and 'compartment_name' can be set. 
  | If a 2D/surface compartment is specified, the complex must be a surface complex and 
  | orientation is set to Orientation.UP.
  | If a 3D/volume compartment is specified, the complex must be a volume complex and
  | orientation is set to Orientation.NONE. 
  | Ignored by derived class Species.


Methods:
*********
* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation


* | **as_species**

   * | return type: Species


  | Returns a Species object based on this Complex. All species-specific 
  | attributes are set to their default values and 'name' is set to value returned by 
  | 'to_bngl_str()'.



Component
=========

Instance of a component belonging to a molecule instance.
A component instance may have its state set.
It is also used to connect molecule instance in a complex instance.

Attributes:
***********
* | **component_type**: ComponentType

* | **state**: str = STATE_UNSET

* | **bond**: int = BOND_UNBOUND


Methods:
*********
* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation.



ComponentType
=============

Attributes:
***********
* | **name**: str

* | **states**: List[str] = None


Methods:
*********
* | **inst**

   * | state: str = STATE_UNSET
   * | bond: int = BOND_UNBOUND
   * | return type: Component


* | **inst**

   * | state: int = STATE_UNSET_INT
   * | bond: int = BOND_UNBOUND
   * | return type: Component


* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation.



Config
======

Attributes:
***********
* | **seed**: int = 1

* | **time_step**: float = 1e-6
  | Default value is 1us, in seconds

* | **surface_grid_density**: float = 10000

* | **interaction_radius**: float = None
  | Diffusing volume molecules will interact with each other when
  | they get within N microns of each other. The default is
  | 1/sqrt(PI \* Sigma_s) where Sigma_s is the surface grid density 
  | (default or user-specified).

* | **intermembrane_interaction_radius**: float = None
  | Diffusing surface molecules will interact with surface molecules on other
  | walls when they get within N microns of each other. The default is
  | 1/sqrt(PI \* Sigma_s) where Sigma_s is the surface grid density 
  | (default or user-specified).

* | **vacancy_search_distance**: float = 10
  | Normally, a reaction will not proceed on a surface unless there
  | is room to place all products on the single grid element where
  | the reaction is initiated. By increasing r from its default value
  | of 0, one can specify how far from the reaction’s location, in microns, the
  | reaction can place its products. To be useful, r must
  | be larger than the longest axis of the grid element on the triangle
  | in question. The reaction will then proceed if there is room to
  | place its products within a radius r, and will place those products as 
  | close as possible to the place where the reaction occurs
  | (deterministically, so small-scale directional bias is possible).

* | **center_molecules_on_grid**: bool = False

* | **initial_partition_origin**: List[float] = None
  | Optional placement of the partition 0 placement, specifies the left, lower and front 
  | point. If not set, value -partition_dimension/2 is used for each of the dimensions 
  | placing the center of the partition to (0, 0, 0).

* | **partition_dimension**: float = 10

* | **subpartition_dimension**: float = 0.5

* | **total_iterations**: float = 1000000
  | Required for checkpointing so that the checkpointed model has information on
  | the intended total number of iterations. 
  | Also used when generating visualization data files and also for other reporting uses. 
  | Value is truncated to an integer.

* | **check_overlapped_walls**: bool = True
  | Enables check for overlapped walls. Overlapping walls can cause issues during 
  | simulation such as a molecule escaping closed geometry when it hits two walls 
  | that overlap.

* | **reaction_class_cleanup_periodicity**: int = 500
  | Reaction class cleanup removes computed reaction classes for inactive species from memory.
  | This provides faster reaction lookup faster but when the same reaction class is 
  | needed again, it must be recomputed.

* | **species_cleanup_periodicity**: int = 10000
  | Species cleanup removes inactive species from memory. It removes also all reaction classes 
  | that reference it.
  | This provides faster addition of new species lookup faster but when the species is 
  | needed again, it must be recomputed.

* | **sort_molecules**: bool = False
  | Enables sorting of molecules for diffusion, this may improve cache locality.
  | Produces different results when enabled.

* | **memory_limit_gb**: int = -1
  | Sets memory limit in GB for simulation run. 
  | When this limit is hit, all buffers are flushed and simulation is terminated with an error.

* | **initial_iteration**: int = 0
  | Initial iteration, used when resuming a checkpoint.

* | **initial_time**: float = 0
  | Initial time in us, used when resuming a checkpoint.
  | Will be truncated to be a multiple of time step.

* | **initial_rng_state**: RngState = None
  | Used for checkpointing, may contain state of the random number generator to be set 
  | after initialization right before the first event is started. 
  | When not set, the set 'seed' value is used to initialize the random number generator.

* | **append_to_count_output_data**: bool = False
  | Used for checkpointing, instead of creating new files for Count observables data, 
  | new values are appended to the existing files. If such files do not exist, new files are
  | created.

* | **continue_after_sigalrm**: bool = False
  | MCell registers a SIGALRM signal handler. When SIGALRM signal is received and 
  | continue_after_sigalrm is False, checkpoint is stored and simulation is terminated. 
  | When continue_after_sigalrm is True, checkpoint is stored and simulation continues.

Count
=====

Attributes:
***********
* | **name**: str = None
  | Name of a count may be specified when one needs to search for them later. 
  | Also when the count is created while loading a BNGL file, its name is set.

* | **file_name**: str = None
  | File name with an optional path must be set. It is not deduced automatically.

* | **count_expression**: CountTerm = None
  | The count expression must be composed only from CountTerm objects that are added or 
  | subtracted.

* | **multiplier**: float = 1
  | In some cases it might be useful to multiply the whole count by a constant to get 
  | for instance concentration. The count_expression is not an arbitrary expression
  | and such multiplication can be done through this attribute.

* | **every_n_timesteps**: float = 1
  | Value is truncated (floored) to an integer.
  | If value is set to 0, this Count is used only on-demand through calls to its
  | get_current_value method.

* | **species_pattern**: Complex = None
  | Count the number of molecules that match the given complex instance pattern.
  | Counts each molecule exactly once. 
  | If the pattern has a compartment set, this specifies the counted region.

* | **molecules_pattern**: Complex = None
  | Count the number of matches of the given pattern on molecules.
  | The observable will count a molecule every time it matches the pattern.
  | When the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a 
  | molecule A(a!1).A(a!1,b!2).B(a!2) will be counted twice because the 
  | pattern may match in two different ways. 
  | If the pattern has a compartment set, this specifies the counted region.

* | **reaction_rule**: ReactionRule = None
  | Count the number of reactions that occurred since the start of the simulation.

* | **region**: Region = None
  | Only a GeometryObject or SurfaceRegion can be passed as the region argument, 
  | compound regions (created with +, -, \*) are not supproted yet.   
  | Cannot be set when 'species_pattern' or 'molecules_pattern' has a  
  | compartment specified.
  | If pattern compartment is not specified and 'region' is left 'unset', 
  | counting is done in the whole world.

* | **node_type**: ExprNodeType = ExprNodeType.LEAF
  | Internal, used to represent an expression

* | **left_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the right operand

* | **initial_reactions_count**: int = 0
  | Used for checkpointing, allows to set initial count of reactions that occurred.
  | Ignored when molecules are counted.


Methods:
*********
* | **get_current_value**

   * | return type: float


  | Returns the current value for this count. Cannot be used to count reactions.
  | The model must be initialized with this Count present as one of the observables.


* | **__add__**

   * | op2: CountTerm
   * | return type: CountTerm


* | **__sub__**

   * | op2: CountTerm
   * | return type: CountTerm



CountTerm
=========

Attributes:
***********
* | **species_pattern**: Complex = None
  | Count the number of molecules that match the given complex instance pattern.
  | Counts each molecule exactly once. 
  | If the pattern has a compartment set, this specifies the counted region.

* | **molecules_pattern**: Complex = None
  | Count the number of matches of the given pattern on molecules.
  | The observable will count a molecule every time it matches the pattern.
  | When the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a 
  | molecule A(a!1).A(a!1,b!2).B(a!2) will be counted twice because the 
  | pattern may match in two different ways. 
  | If the pattern has a compartment set, this specifies the counted region.

* | **reaction_rule**: ReactionRule = None
  | Count the number of reactions that occurred since the start of the simulation.

* | **region**: Region = None
  | Only a GeometryObject or SurfaceRegion can be passed as the region argument, 
  | compound regions (created with +, -, \*) are not supproted yet.   
  | Cannot be set when 'species_pattern' or 'molecules_pattern' has a  
  | compartment specified.
  | If pattern compartment is not specified and 'region' is left 'unset', 
  | counting is done in the whole world.

* | **node_type**: ExprNodeType = ExprNodeType.LEAF
  | Internal, used to represent an expression

* | **left_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the right operand

* | **initial_reactions_count**: int = 0
  | Used for checkpointing, allows to set initial count of reactions that occurred.
  | Ignored when molecules are counted.


Methods:
*********
* | **__add__**

   * | op2: CountTerm
   * | return type: CountTerm


* | **__sub__**

   * | op2: CountTerm
   * | return type: CountTerm



ElementaryMolecule
==================

Attributes:
***********
* | **elementary_molecule_type**: ElementaryMoleculeType

* | **components**: List[Component] = None


Methods:
*********
* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation



ElementaryMoleculeType
======================

Attributes:
***********
* | **name**: str

* | **components**: List[ComponentType] = None

* | **diffusion_constant_2d**: float = None
  | This molecule is constrained to a surface and diffuses with diffusion constant D.

* | **diffusion_constant_3d**: float = None
  | This molecule diffuses in space with diffusion constant D. D can be zero, in which case the molecule doesn’t move. The units of D are cm 2 /s.

* | **custom_time_step**: float = None
  | This molecule should take timesteps of length t (in seconds). Use either this or custom_time_step.

* | **custom_space_step**: float = None
  | This molecule should take steps of average length L (in microns). Use either this or custom_time_step.

* | **target_only**: bool = False
  | This molecule will not initiate reactions when it runs into other molecules. This
  | setting can speed up simulations when applied to a molecule at high concentrations 
  | that reacts with a molecule at low concentrations (it is more efficient for
  | the low-concentration molecule to trigger the reactions). This directive does
  | not affect unimolecular reactions.


Methods:
*********
* | **inst**

   * | components: List[Component] = None
   * | return type: ElementaryMolecule


* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation.



GeometryObject
==============

Attributes:
***********
* | **name**: str
  | Name of the object. Also represents BNGL compartment name if 'is_bngl_compartment' is True.

* | **vertex_list**: List[List[float]]
  | List of [x,y,z] triplets specifying positions of individual vertices.
  | Equivalent to List[Vec3] however, defining a constructor Vec3(List[float]) then 
  | tries to convert all lists of floats to Vec3

* | **wall_list**: List[List[int]]
  | List of [a,b,c] triplets specifying each wall, individual values are indices into the vertex list.
  | Equivalent to List[IVec3].

* | **is_bngl_compartment**: bool = False

* | **surface_compartment_name**: str = None

* | **surface_regions**: List[SurfaceRegion] = None

* | **surface_class**: SurfaceClass = None
  | Surface class for the whole object's surface. It is applied to the whole surface of this object 
  | except for those surface regions that have their specific surface class set explicitly.

* | **initial_surface_releases**: List[InitialSurfaceRelease] = None
  | Equivalent to MDL's MODIFY_SURFACE_REGIONS/MOLECULE_DENSITY or MOLECULE_NUMBER,
  | each item defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.

* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods:
*********
* | **translate**

   * | move: Vec3

  | Move object by a specified vector, must be done before model initialization.


* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of thwo regions


* | **__sub__**

   * | other: Region
   * | return type: Region


* | **__mul__**

   * | other: Region
   * | return type: Region



InitialSurfaceRelease
=====================

Defines molecules to be released onto a SurfaceRegion right when simulation starts

Attributes:
***********
* | **complex**: Complex

* | **number_to_release**: int = None
  | Number of molecules to be released onto a region,
  | only one of number_to_release and density can be set.

* | **density**: float = None
  | Density of molecules to be released onto a region,
  | only one of number_to_release and density can be set.

Instantiation
=============

Attributes:
***********
* | **release_sites**: List[ReleaseSite] = None

* | **geometry_objects**: List[GeometryObject] = None

* | **checkpointed_molecules**: List[BaseChkptMol] = None
  | Used when resuming simulation from a checkpoint.


Methods:
*********
* | **add_release_site**

   * | s: ReleaseSite

  | Makes a copy of the release site


* | **find_release_site**

   * | name: str
   * | return type: ReleaseSite


* | **add_geometry_object**

   * | o: GeometryObject

  | Makes a copy of the geometry object, in the future we will probably add some transformations


* | **find_geometry_object**

   * | name: str
   * | return type: GeometryObject


* | **find_volume_compartment**

   * | name: str
   * | return type: GeometryObject


* | **find_surface_compartment**

   * | name: str
   * | return type: GeometryObject


* | **load_bngl_seed_species**

   * | file_name: str
   * | subsystem: Subsystem
   * | default_release_region: Region = None
     | Used for seed species that have no compartments specified

   * | parameter_overrides: Dict[str, float] = None

  | Loads section seed species from a BNGL file and creates release sites according to it.
  | All elementary molecule types used in the seed species section must be already defined in subsystem.
  | If an item in the BNGL seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region.



Introspection
=============

This class is used only as a base class to Model, it is not provided through API. Provides methods to introspect simulation state.


Methods:
*********
* | **get_molecule_ids**

   * | species: Species = None
   * | return type: List[int]


  | Returns a list of ids of molecules of given Species existing in the simulated environment,
  | if the argument species is not set, returns list of all molecules.


* | **get_molecule**

   * | id: int
   * | return type: Molecule


  | Returns a molecule from the simulated environment, None if the molecule does not exist


* | **get_species_name**

   * | species_id: int
   * | return type: str


  | Returns a string representing canonical species name in the BNGL format.


* | **get_vertex**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in object's walls (wall_list).

   * | return type: Vec3


  | Returns coordinates of a vertex.


* | **get_wall**

   * | object: GeometryObject
   * | wall_index: int
     | This is the index of the wall in object's walls (wall_list).

   * | return type: Wall


  | Returns information about a wall belonging to a given object.


* | **get_vertex_unit_normal**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in object's vertex_list.

   * | return type: Vec3


  | Returns sum of all wall normals that use this vertex converted to a unit vector of length 1um.
  | This represents the unit vector pointing outwards from the vertex.


* | **get_wall_unit_normal**

   * | object: GeometryObject
   * | wall_index: int
     | This is the index of the vertex in object's walls (wall_list).

   * | return type: Vec3


  | Returns wall normal converted to a unit vector of length 1um.



Model
=====

Attributes:
***********
* | **config**: Config = Config()

* | **warnings**: Warnings = Warnings()

* | **notifications**: Notifications = Notifications()

* | **species**: List[Species] = None

* | **reaction_rules**: List[ReactionRule] = None

* | **surface_classes**: List[SurfaceClass] = None

* | **elementary_molecule_types**: List[ElementaryMoleculeType] = None
  | Used mainly when a BNGL file is loaded, if BNGL species is defined through 
  | Python API, this array is populated automatically

* | **release_sites**: List[ReleaseSite] = None

* | **geometry_objects**: List[GeometryObject] = None

* | **checkpointed_molecules**: List[BaseChkptMol] = None
  | Used when resuming simulation from a checkpoint.

* | **viz_outputs**: List[VizOutput] = None

* | **counts**: List[Count] = None


Methods:
*********
* | **initialize**


  | Initializes model, initialization blocks most of changes to 
  | contained components (the attributes


* | **run_iterations**

   * | iterations: float
     | Number of iterations to run. Value is truncated to an integer.

   * | return type: int


  | Runs specified number of iterations. Returns the number of iterations
  | executed (it might be less than the requested number of iterations when 
  | a checkpoint was scheduled).


* | **end_simulation**

   * | print_final_report: bool = True

  | Generates the last visualization and reaction output (if they were defined), then
  | flushes all buffers and optionally prints simulation report. 
  | Buffers are also flushed when the Model object is destroyed.


* | **add_subsystem**

   * | subsystem: Subsystem

* | **add_instantiation**

   * | instantiation: Instantiation

* | **add_observables**

   * | observables: Observables

* | **dump_internal_state**


  | Prints out the simulation engine's internal state, mainly for debugging.


* | **export_data_model**

   * | file: str = None

  | If file is not set, then uses the first VizOutput to determine the target directory 
  | and creates name using the current iteration. Fails if argument file is not set and there is no VizOutput.
  | Must be called after initialization.
  | Always exports the current state, i.e. with the current . 
  | Events (ReleaseSites and VizOutputs) with scheduled time other than zero cannot be imported correectly yet.


* | **export_viz_data_model**

   * | file: str = None

  | Same as export_data_model, only the created data model will contain only information required for visualization in CellBlender. This makes the loading ofthemodel by CellBlender faster and also allows to avoid potential compatibility issues.


* | **release_molecules**

   * | release_site: ReleaseSite

  | Performs immediate release based on the definition of the release site argument.
  | The ReleaseSite.release_time must not be in the past and should be withing the current iteration.
  | The ReleaseEvent must not use a release_pattern because this is an immediate release and it is not 
  | scheduled into the global scheduler.


* | **run_reaction**

   * | reaction_rule: ReactionRule
     | Reaction rule to run.

   * | reactant_ids: List[int]
     | The number of reactants for a unimolecular reaction must be 1 and for a bimolecular reaction must be 2.
     | Reactants for a bimolecular reaction do not have to be listed in the same order as in the reaction rule definition.

   * | time: float
     | Precise time in seconds when this reaction occurs. Important to know for how long the products
     | will be diffused when they are created in a middle of a time step.

   * | return type: List[int]


  | Run a single reaction on reactants. Callbacks will be called if they are registered for the given reaction.
  | Returns a list of product IDs.
  | Note\: only unimolecular reactions are currently supported.


* | **add_vertex_move**

   * | object: GeometryObject
     | Object whose vertex will be changed

   * | vertex_index: int
     | Index of vertex in object's vertex list that will be changed

   * | displacement: Vec3
     | Change of vertex coordinates (in um), will be added to the current coordinates of the vertex


  | Adds a displacement for given object's vertex, only stored until apply_vertex_moves is called


* | **apply_vertex_moves**

   * | collect_wall_wall_hits: bool = False
     | When set to True, a list of wall pairs that collided is returned,
     | otherwise an empty list is returned.

   * | return type: List[WallWallHitInfo]


  | Applies all the vertex moves specified with add_vertex_move call.
  | Walls of different objects are checked against collisions and move the maximal way so that they do not 
  | overlap. (the current pllementation is a bit basic and may not work 100% correctly) 
  | When collect_wall_wall_hits is True, a list of wall pairs that collided is returned,
  | when collect_wall_wall_hits is False, and empty list is returned.


* | **register_mol_wall_hit_callback**

   * | function: Callable, # std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)>
     | Callback function to be called. 
     | It must have two arguments MolWallHitInfo and context.

   * | context: Any, # py::object
     | Context passed to the callback function, the callback function can store
     | information to this object. Some context must be always passed, even when 
     | it is a useless python object.

   * | object: GeometryObject = None
     | Only hits of this object will be reported, any object hit is reported when not set.

   * | species: Species = None
     | Only hits of molecules of this species will be reported, any species hit is reported when not set.


  | There can be currently only a single wall hit callback registered.


* | **register_reaction_callback**

   * | function: Callable, # std::function<void(std::shared_ptr<ReactionInfo>, py::object)>
     | Callback function to be called. 
     | It must have two arguments ReactionInfo and context.
     | Called when it is decided that the reaction will happen.
     | After return the reaction proceeds as it would without a callback.

   * | context: Any, # py::object
     | Context passed to the callback function, the callback function can store
     | information to this object. Some context must be always passed, even when 
     | it is a useless python object.

   * | reaction_rule: ReactionRule
     | The callback function will be called whenever is this reaction rule applied.


  | Defines a function to be called when a reaction was processed.
  | It is allowed to do state modifications except for removing reacting molecules, 
  | they will be removed automatically after return from this callback.


* | **load_bngl**

   * | file_name: str
   * | observables_files_prefix: str = ''
     | Prefix to be used when creating files with observable values.

   * | default_release_region: Region = None
   * | parameter_overrides: Dict[str, float] = None

  | Loads sections\: molecule types, reaction rules, seed species, and observables from a BNGL file
  | and creates objects in the current model according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.
  | If an item in the seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region.


* | **export_to_bngl**

   * | file_name: str
     | Output file name.


  | Exports all defined species, reaction rules and applicable observables
  | as a BNGL file. 
  | Limited currrently to exactly one volume compartment and volume reactions.


* | **save_checkpoint**

   * | custom_dir: str = None
     | Sets custom directory where the checkpoint will be stored. 
     | The default is 'checkpoints/seed_<SEED>/it_<ITERATION>'.


  | Saves current model state as checkpoint. 
  | The default directory structure is checkpoints/seed_<SEED>/it_<ITERATION>,
  | it can be changed by setting 'custom_dir'.
  | If used during an iteration, schedules an event for the end of the current iteration
  | that saves the checkpoint (effectively calls 'checkpoint_after_iteration(0, False, custom_dir)'.


* | **schedule_checkpoint**

   * | iteration: int = 0
     | Specifies iteration number when the checkpoint save will occur. 
     | Please note that iterations are counted from 0, i.e. a call run_iterations(3)
     | runs iterations 0, 1, 2.
     | To schedule a checkpoint for the closest time as possible, keep the default value 0. 
     | If calling schedule_checkpoint from a different thread (e.g. by using threading.Timer), 
     | it is highly recommended to keep the default value 0 or choose some time that will be 
     | for sure in the future.

   * | continue_simulation: bool = False
     | When false, saving the checkpoint means that we want to terminate the simulation 
     | right after the save, the currently running function Model.run_iterations
     | does not simulate any following iterations and execution returns from this function
     | to execute the next statement which is usually 'model.end_simulation()'.
     | When true, the checkpoint is just saved and simulation continues uninterrupted.

   * | custom_dir: str = None
     | Sets custom directory where the checkpoint will be stored. 
     | The default is 'checkpoints/seed_<SEED>/it_<ITERATION>'.


  | Schedules checkpoint save that will occur when an iteration is started  
  | right before any other events scheduled for the given iteration are executed.
  | Can be called asynchronously at any time after initialization.


* | **add_species**

   * | s: Species

* | **find_species**

   * | name: str
   * | return type: Species


* | **add_reaction_rule**

   * | r: ReactionRule

* | **find_reaction_rule**

   * | name: str
   * | return type: ReactionRule


* | **add_surface_class**

   * | sc: SurfaceClass

* | **find_surface_class**

   * | name: str
   * | return type: SurfaceClass


* | **add_elementary_molecule_type**

   * | mt: ElementaryMoleculeType

* | **find_elementary_molecule_type**

   * | name: str
   * | return type: ElementaryMoleculeType


* | **load_bngl_molecule_types_and_reaction_rules**

   * | file_name: str
   * | parameter_overrides: Dict[str, float] = None

  | Parses a BNGL file and only reads molecule types and
  | reaction rules sections, e.g. ignores observables. 
  | Parameter values are evaluated and the result value 
  | is directly used.  
  | Compartments names are stored in rxn rules as strings because
  | compartments belong to geometry objects and the subsystem is independent
  | on specific geometry.
  | However they must be defined on initialization.


* | **add_release_site**

   * | s: ReleaseSite

  | Makes a copy of the release site


* | **find_release_site**

   * | name: str
   * | return type: ReleaseSite


* | **add_geometry_object**

   * | o: GeometryObject

  | Makes a copy of the geometry object, in the future we will probably add some transformations


* | **find_geometry_object**

   * | name: str
   * | return type: GeometryObject


* | **find_volume_compartment**

   * | name: str
   * | return type: GeometryObject


* | **find_surface_compartment**

   * | name: str
   * | return type: GeometryObject


* | **load_bngl_seed_species**

   * | file_name: str
   * | subsystem: Subsystem
   * | default_release_region: Region = None
     | Used for seed species that have no compartments specified

   * | parameter_overrides: Dict[str, float] = None

  | Loads section seed species from a BNGL file and creates release sites according to it.
  | All elementary molecule types used in the seed species section must be already defined in subsystem.
  | If an item in the BNGL seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region.


* | **add_viz_output**

   * | viz_output: VizOutput

* | **add_count**

   * | count: Count

* | **find_count**

   * | name: str
   * | return type: Count


* | **load_bngl_observables**

   * | file_name: str
     | BNGL file name.

   * | subsystem: Subsystem
   * | output_files_prefix: str = ''
     | Prefix to be used when creating files with observable values.

   * | parameter_overrides: Dict[str, float] = None

  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.


* | **get_molecule_ids**

   * | species: Species = None
   * | return type: List[int]


  | Returns a list of ids of molecules of given Species existing in the simulated environment,
  | if the argument species is not set, returns list of all molecules.


* | **get_molecule**

   * | id: int
   * | return type: Molecule


  | Returns a molecule from the simulated environment, None if the molecule does not exist


* | **get_species_name**

   * | species_id: int
   * | return type: str


  | Returns a string representing canonical species name in the BNGL format.


* | **get_vertex**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in object's walls (wall_list).

   * | return type: Vec3


  | Returns coordinates of a vertex.


* | **get_wall**

   * | object: GeometryObject
   * | wall_index: int
     | This is the index of the wall in object's walls (wall_list).

   * | return type: Wall


  | Returns information about a wall belonging to a given object.


* | **get_vertex_unit_normal**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in object's vertex_list.

   * | return type: Vec3


  | Returns sum of all wall normals that use this vertex converted to a unit vector of length 1um.
  | This represents the unit vector pointing outwards from the vertex.


* | **get_wall_unit_normal**

   * | object: GeometryObject
   * | wall_index: int
     | This is the index of the vertex in object's walls (wall_list).

   * | return type: Vec3


  | Returns wall normal converted to a unit vector of length 1um.



MolWallHitInfo
==============

Attributes:
***********
* | **molecule_id**: int

* | **geometry_object**: GeometryObject
  | Object that was hit.

* | **wall_index**: int
  | Index of wall belonging to the geometry_object.

* | **time**: float
  | Time of the hit

* | **pos3d**: Vec3
  | Position of the hit

* | **time_before_hit**: float
  | Time when the molecule started to diffuse towards the hit wall. 
  | It is either the start of the molecule's diffusion or 
  | if a wall was hit later then the time of last wall hit.

* | **pos3d_before_hit**: Vec3
  | Position of the molecule at time_before_hit

Molecule
========

This is a Python representation of a molecule obtained from Model 
during simulation.

Attributes:
***********
* | **id**: int = ID_INVALID
  | Unique id of this molecule

* | **type**: MoleculeType = MoleculeType.UNSET

* | **species_id**: int = ID_INVALID
  | Species id of this molecule.
  | The id value is only temporary and can be invalidated by simulating an iteration.

* | **pos3d**: Vec3 = None
  | Contains position of a molecule in 3D space.

* | **orientation**: Orientation = Orientation.NOT_SET
  | Contains orientation for surface molecule. Volume molecules 
  | have always orientation set to Orientation.NONE.

* | **pos2d**: Vec2 = None
  | Set only for surface molecules.

* | **geometry_object**: GeometryObject = None
  | Set only for surface molecules.
  | Object on whose surface is the molecule located.

* | **wall_index**: int = -1
  | Set only for surface molecules.
  | Index of wall belonging to the geometry_object where is the 
  | molecule located.


Methods:
*********
* | **remove**


  | Removes this molecule from simulation. Any subsequent modifications
  | of this object won't have any effect.



MoleculeReleaseInfo
===================

Attributes:
***********
* | **complex**: Complex
  | Complex instance defining the molecule that will be released.
  | Orientation of the complex instance is used to define orientation of the released molecule,
  | when Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and
  | surface molecules are released with Orientation.UP.
  | Compartment must not be set because this specific release definition states the location.

* | **location**: List[float]
  | 3D position where the molecule will be released. 
  | If a molecule has a 2D diffusion constant, it will be
  | placed on the surface closest to the coordinate given. 
  | Argument must have exactly three floating point values.

Notifications
=============

Attributes:
***********
* | **bng_verbosity_level**: int = 0
  | Sets verbosity level that enables printouts of extra information on BioNetGen 
  | species and rules created and used during simulation.

* | **rxn_and_species_report**: bool = True
  | Simulation generates files rxn_report_SEED.txt species_report_SEED.txt that contain
  | details on reaction classes and species that were created based on reaction rules.

* | **simulation_stats_every_n_iterations**: int = 0
  | When set to a value other than 0, internal simulation stats will be printed.

Observables
===========

Neither VizOutput, nor Count have name, therefore there are no find_* methods.

Attributes:
***********
* | **viz_outputs**: List[VizOutput] = None

* | **counts**: List[Count] = None


Methods:
*********
* | **add_viz_output**

   * | viz_output: VizOutput

* | **add_count**

   * | count: Count

* | **find_count**

   * | name: str
   * | return type: Count


* | **load_bngl_observables**

   * | file_name: str
     | BNGL file name.

   * | subsystem: Subsystem
   * | output_files_prefix: str = ''
     | Prefix to be used when creating files with observable values.

   * | parameter_overrides: Dict[str, float] = None

  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.



ReactionInfo
============

Data structure passed to a reaction callback.

Attributes:
***********
* | **type**: ReactionType
  | Specifies whether the reaction is unimolecular or bimolecular and
  | also provides information in reactant types.

* | **reactant_ids**: List[int]
  | IDs of the reacting molecules, contains 1 ID for a unimolecular reaction, 2 IDs for a bimolecular reaction.
  | For a bimolecular reaction, the first ID is always the molecule that was diffused and the second one 
  | is the molecule that was hit.
  | IDs can be used to obtain location of the molecules. The position of the first molecule obtained through 
  | model.get_molecule() is the position of the diffusing molecule before the collision.
  | All the reactants are removed after return from this callback, unless they are kept by the reaction such as in A + B -> A + C.

* | **product_ids**: List[int]
  | IDs of reaction product molecules. They already exist in the simulated system together with reactants, however reactants 
  | will be removed after return from this callback.

* | **reaction_rule**: ReactionRule
  | Reaction rule of the reaction.

* | **time**: float
  | Time of the reaction

* | **pos3d**: Vec3
  | Specifies where reaction occured in the 3d space, specific meaning depends on the reaction type\:
  | - unimolecular reaction - position of the reacting molecule,
  | - volume-volume or surface-surface reaction - position of the first reactant,
  | - volume-surface reaction - position where the volume molecule hit the wall with the surface molecule.

* | **geometry_object**: GeometryObject = None
  | Set only for surface reactions.
  | Object on whose surface where the reaction occured.

* | **wall_index**: int = -1
  | Set only for surface reactions.
  | Index of wall belonging to the geometry_object where the reaction occured, 
  | i.e. where the volume molecule hit the wall with a surface molecule or
  | wall where the diffusing surface reactant reacted.

* | **pos2d**: Vec2 = None
  | Set only for surface reactions.
  | Specifies where reaction occured in the 2d UV coordinates defined by the wall where the reaction occured, 
  | specific meaning depends on the reaction type\:
  | - unimolecular reaction - position of the reacting molecule,
  | - volume-surface and surface-surface reaction - position of the second reactant.

ReactionRule
============

Attributes:
***********
* | **name**: str = None
  | Name of the reaction. If this is a reversible reaction, then it is the name of the 
  | reaction in forward direction.

* | **reactants**: List[Complex] = None

* | **products**: List[Complex] = None

* | **fwd_rate**: float = None
  | Rates have following units\: unimolecular [s^-1], volume bimolecular [M^-1\*s^-1], 
  | The units of the reaction rate for uni- and bimolecular reactions are
  |   \* [s^-1] for unimolecular reactions,
  |   \* [M^-1\*s^-1] for bimolecular reactions between either two volume molecules, a volume molecule 
  |                 and a surface (molecule), 
  |   \* [um^2\*N^-1\*s^-1] bimolecular reactions between two surface molecules on the same surface, and
  |   \* [N^-1\*s^-1] bimolecular reactions between two surface molecules on different objects 
  |     (this is a highly experimental feature and the unit will likely change in the future, 
  |      not sure if probability is computed correctly, it works the way that the surface molecule 
  |      is first diffused and then a potential collisions within the distance of Config.intermembrane_interaction_radius
  |      are evaluated). 
  | Here, M is the molarity of the solution and N the number of reactants.
  | May be changed after model initialization. 
  | Setting of value is ignored if the rate does not change. 
  | If the new value differs from previous, updates all information related 
  | to the new rate including recomputation of reaction times for molecules if this is a
  | unimolecular reaction.

* | **rev_name**: str = None
  | Name of the reaction in reverse direction.

* | **rev_rate**: float = None
  | Reverse reactions rate, reaction is unidirectional when not specified.
  | May be changed after model initialization, in the case behaves the same was as for 
  | changing the 'fwd_rate'.

* | **variable_rate**: List[List[float]] = None
  | Variable rate is applicable only for irreversible reactions. Members fwd_rate and rev_rate 
  | must not be set. The array passed as this argument must have as its items a pair of floats (time, rate).

* | **is_intermembrane_surface_reaction**: bool = False
  | Experimental, see addintinal explanation in 'fwd' rate.
  | Then set to true, this is a special type of surface-surface reaction that 
  | allows for two surface molecules to react when they are on different geometrical objects. 
  | This support is limited for now, the reaction rule must be in the form of A + B -> C + D 
  | where all reactants and products must be surface molecules and 
  | their orientation must be 'any' (default).


Methods:
*********
* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to the reaction rule's BNGL representation, does not contain rates.



Region
======

Represents region construted from 1 or more multiple, usually unnamed?

Attributes:
***********
* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods:
*********
* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of thwo regions


* | **__sub__**

   * | other: Region
   * | return type: Region


* | **__mul__**

   * | other: Region
   * | return type: Region



ReleasePattern
==============

Attributes:
***********
* | **name**: str = None
  | Name of the release pattern

* | **release_interval**: float = TIME_INFINITY
  | During a train of releases, release molecules after every t seconds. 
  | Default is to release only once.

* | **train_duration**: float = TIME_INFINITY
  | The train of releases lasts for t seconds before turning off. 
  | Default is to never turn off.

* | **train_interval**: float = TIME_INFINITY
  | A new train of releases happens every t seconds. 
  | Default is to never have a new train. 
  | The train interval must not be shorter than the train duration.

* | **number_of_trains**: int = 1
  | Repeat the release process for n trains of releases. Default is one train.
  | For unlimited number of trains use constant NUMBER_OF_TRAINS_UNLIMITED.

ReleaseSite
===========

Attributes:
***********
* | **name**: str
  | Name of the release site

* | **complex**: Complex = None
  | Defines the species of the molecule that will be released. Not used for the LIST shape. 
  | Must be set when molecule_list is empty and unset when molecule_list is not empty.
  | Orientation of the complex instance is used to define orientation of the released molecule,
  | when Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and
  | surface molecules are released with Orientation.UP.
  | When compartment is specified, this sets shape to Shape.COMPARTMENT and the molecules are released 
  | into the compartment.

* | **molecule_list**: List[MoleculeReleaseInfo] = None
  | Used for LIST shape release mode. 
  | Only one of number_to_release, density, concentration or molecule_list can be set.

* | **release_time**: float = 0
  | Specifies time in seconds when the release event is executed.
  | In case when a release pattern is used, this is the time of the first release.      
  | Equivalent to MDL's RELEASE_PATTERN command DELAY.

* | **release_pattern**: ReleasePattern = None
  | Use the release pattern to define schedule of releases. 
  | The default is to release the specified number of molecules at the set release_time.

* | **shape**: Shape = Shape.UNSET
  | Set automatically when

* | **region**: Region = None
  | Sets shape to Shape.REGION_EXPR.

* | **location**: Vec3 = None

* | **site_diameter**: float = 0
  | For a geometrical release site, this releases molecules uniformly within
  | a radius r. Not used for releases on regions.
  | Usually required for Shape.List type of releases.

* | **site_radius**: float = None
  | For a geometrical release site, this releases molecules uniformly within
  | a radius r. Not used for releases on regions.

* | **number_to_release**: float = None
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Value is truncated (floored) to an integer.

* | **density**: float = None
  | Only one of number_to_release, density, concentration or molecule_list can be set.

* | **concentration**: float = None
  | Only one of number_to_release, density, concentration or molecule_list can be set.

* | **release_probability**: float = None

RngState
========

Internal checkpointing structure holding state of the random number generator.

Attributes:
***********
* | **randcnt**: int

* | **aa**: int

* | **bb**: int

* | **cc**: int

* | **randslr**: List[int]
  | Must contain RNG_SIZE items.

* | **mm**: List[int]
  | Must contain RNG_SIZE items.

* | **rngblocks**: int
  | Must contain RNG_SIZE items.

Species
=======

There are three ways how to use this class\:
1) definition of simple species - in this case 'name' is 
a single identifier and at least 'diffusion_constant_2d' or 
'diffusion_constant_3d' must be provided.
Example\: m.Species('A', diffusion_constant_3d=1e-6). 
Such a definition must be added to subsystem or model so that  
during model initialization this species is transformed to MCell 
representation and an ElementaryMoleculeType 'A' with a given 
diffusion constant is created as well.
2) full definition of complex species - in this case the 
inherited attribute 'elementary_molecules' from Complex
is used as a definition of the complex and this gives information 
on diffusion constants of the used elementary molecules.
Example\: m.Species(elementary_molecules=[ei1, ei2]). 
Such a definition must be added to subsystem or model.   
3) declaration of species - in this case only 'name' in the form of 
an BNGL string is provided. The complex instance specified by the name 
must be fully qualified (i.e. all components are present and those 
components that have a state have their state set).
No information on diffusion constants and other properties of 
used elementary molecules is provided, it must be provided elsewhere.
Example\: m.Species('A(b!1).B(a!1)').
This is a common form of usage when reaction rules are provided in a BNGL file.
Such declaration does no need to be added to subsystem or model.
This form is used as argument in cases where a fully qualified instance  
must be provided such as in molecule releases.

Attributes:
***********
* | **name**: str = None
  | Name of the species in the BNGL format. 
  | One must either specify 'name' or 'elementary_molecules' 
  | (inherited from Complex). This argument 'name' is parsed during model 
  | initialization.

* | **diffusion_constant_2d**: float = None
  | This molecule is constrained to a surface and diffuses with diffusion constant D.

* | **diffusion_constant_3d**: float = None
  | This molecule diffuses in space with diffusion constant D. D can be zero, in which case the molecule doesn’t move. The units of D are cm 2 /s.

* | **custom_time_step**: float = None
  | This molecule should take timesteps of length t (in seconds). Use either this or custom_time_step.

* | **custom_space_step**: float = None
  | This molecule should take steps of average length L (in microns). Use either this or custom_time_step.

* | **target_only**: bool = False
  | This molecule will not initiate reactions when it runs into other molecules. This
  | setting can speed up simulations when applied to a molecule at high concentrations 
  | that reacts with a molecule at low concentrations (it is more efficient for
  | the low-concentration molecule to trigger the reactions). This directive does
  | not affect unimolecular reactions.

* | **name**: str = None
  | When set, this complex instance is initialized from a BNGL string passed as this argument, 
  | the string is parsed during model initialization so the molecule types it uses
  | don't have to be defined before initialization.

* | **elementary_molecules**: List[ElementaryMolecule] = None
  | Individual molecule instances contained in the complex.

* | **orientation**: Orientation = Orientation.DEFAULT
  | Specifies orientation of a molecule. 
  | When Orientation.DEFAULT if kept then during model initialization is
  | 'orientation' set to Orientation.NONE for volume complexes and to 
  | Orientation.UP for surface complexes.
  | Ignored by derived class Species.

* | **compartment_name**: str = None
  | Specifies compartment name of this Complex. 
  | Only one of 'orientation' and 'compartment_name' can be set. 
  | If a 2D/surface compartment is specified, the complex must be a surface complex and 
  | orientation is set to Orientation.UP.
  | If a 3D/volume compartment is specified, the complex must be a volume complex and
  | orientation is set to Orientation.NONE. 
  | Ignored by derived class Species.


Methods:
*********
* | **inst**

   * | orientation: Orientation = Orientation.DEFAULT
     | Maximum one of orientation or compartment_name can be set, not both.

   * | compartment_name: str = None
     | Maximum one of orientation or compartment_name can be set, not both.

   * | return type: Complex


  | Creates a Complex of this Species. Can be currently used only for simple species, i.e. those that
  | have a single molecule instance and no components.


* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation


* | **as_species**

   * | return type: Species


  | Returns a Species object based on this Complex. All species-specific 
  | attributes are set to their default values and 'name' is set to value returned by 
  | 'to_bngl_str()'.



Subsystem
=========

Attributes:
***********
* | **species**: List[Species] = None

* | **reaction_rules**: List[ReactionRule] = None

* | **surface_classes**: List[SurfaceClass] = None

* | **elementary_molecule_types**: List[ElementaryMoleculeType] = None
  | Used mainly when a BNGL file is loaded, if BNGL species is defined through 
  | Python API, this array is populated automatically


Methods:
*********
* | **add_species**

   * | s: Species

* | **find_species**

   * | name: str
   * | return type: Species


* | **add_reaction_rule**

   * | r: ReactionRule

* | **find_reaction_rule**

   * | name: str
   * | return type: ReactionRule


* | **add_surface_class**

   * | sc: SurfaceClass

* | **find_surface_class**

   * | name: str
   * | return type: SurfaceClass


* | **add_elementary_molecule_type**

   * | mt: ElementaryMoleculeType

* | **find_elementary_molecule_type**

   * | name: str
   * | return type: ElementaryMoleculeType


* | **load_bngl_molecule_types_and_reaction_rules**

   * | file_name: str
   * | parameter_overrides: Dict[str, float] = None

  | Parses a BNGL file and only reads molecule types and
  | reaction rules sections, e.g. ignores observables. 
  | Parameter values are evaluated and the result value 
  | is directly used.  
  | Compartments names are stored in rxn rules as strings because
  | compartments belong to geometry objects and the subsystem is independent
  | on specific geometry.
  | However they must be defined on initialization.



SurfaceClass
============

Defining a surface class allows surfaces to behave like species (in a limited way).

Attributes:
***********
* | **name**: str
  | Name of the surface class

* | **properties**: List[SurfaceProperty] = None
  | A surface class can either have a list of properties or just one property.
  | In the usual case of having one property, one can set the attributes 
  | type, affected_species, etc. inherited from SurfaceProperty directly.

* | **type**: SurfacePropertyType = SurfacePropertyType.UNSET
  | Must be set.

* | **affected_complex_pattern**: Complex = None
  | A complex pattern with optional orientation must be set.
  | Default orientation means that the pattern matches any orientation.
  | For concentration or flux clamp the orientation specifies on which side  
  | will be the concentration held 
  | (UP is front or outside, DOWN is back or inside, and DEFAULT, ANY or NONE is on both sides).
  | The complex pattern must not have any compartment.

* | **concentration**: float = None
  | Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or 
  | SurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side 
  | of the wall that has this concentration or flux clamped.

SurfaceProperty
===============

Attributes:
***********
* | **type**: SurfacePropertyType = SurfacePropertyType.UNSET
  | Must be set.

* | **affected_complex_pattern**: Complex = None
  | A complex pattern with optional orientation must be set.
  | Default orientation means that the pattern matches any orientation.
  | For concentration or flux clamp the orientation specifies on which side  
  | will be the concentration held 
  | (UP is front or outside, DOWN is back or inside, and DEFAULT, ANY or NONE is on both sides).
  | The complex pattern must not have any compartment.

* | **concentration**: float = None
  | Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or 
  | SurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side 
  | of the wall that has this concentration or flux clamped.

SurfaceRegion
=============

Surface region  in MDL, however a new class Region was instroduced in MCell4 so it was renamed 
to avoid confusion.

Attributes:
***********
* | **name**: str

* | **wall_indices**: List[int]
  | Surface region must be a part of a GeometryObject, items in this list are indices to 
  | its wall_list array

* | **surface_class**: SurfaceClass = None
  | Has higher priority than the parent geometry object's surface class.

* | **initial_surface_releases**: List[InitialSurfaceRelease] = None
  | Equivalent to MDL's MODIFY_SURFACE_REGIONS/MOLECULE_DENSITY or MOLECULE_NUMBER,
  | each item defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.

* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods:
*********
* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of thwo regions


* | **__sub__**

   * | other: Region
   * | return type: Region


* | **__mul__**

   * | other: Region
   * | return type: Region



VizOutput
=========

Attributes:
***********
* | **output_files_prefix**: str

* | **species_list**: List[Species] = None
  | Specifies a list of species to be visualized, when empty, all_species will be generated.

* | **mode**: VizMode = VizMode.ASCII

* | **every_n_timesteps**: float = 1
  | Value is truncated (floored) to an integer.
  | Value 0 means that the viz output is ran only once at iteration 0.

Wall
====

This is a Python representation of a molecule obtained from Model 
during simulation.

Attributes:
***********
* | **geometry_object**: GeometryObject
  | Object to which this wall belongs.

* | **wall_index**: int
  | Index of this wall in the object to which this wall belongs.

* | **vertices**: List[Vec3]
  | Vertices of the triangle that represents this wall.

* | **area**: float

* | **unit_normal**: Vec3
  | Normal of this wall with unit length of 1 um.
  | To get just the unit vector, not the whole wall, there is also method Model.get_wall_unit_normal.

* | **is_movable**: bool = True
  | If True, whis wall can be moved through Model.apply_vertex_moves,
  | if False, wall moves are ignored.

WallWallHitInfo
===============

Attributes:
***********
* | **wall1**: Wall

* | **wall2**: Wall

Warnings
========

This is a placeholder for future warnings settings. Empty for now.

bngl_utils
==========


Methods:
*********
* | **load_bngl_parameters**

   * | file_name: str
   * | parameter_overrides: Dict[str, float] = None
   * | return type: Dict[str, float]



geometry_utils
==============


Methods:
*********
* | **create_box**

   * | name: str
     | Name of the created geometry object

   * | edge_length: float
     | Specifies length of each edge of the box.

   * | return type: GeometryObject


  | Creates a GeometryObject whose center is at (0, 0, 0).



run_utils
=========


Methods:
*********
* | **get_last_checkpoint_dir**

   * | seed: int
   * | return type: str


  | Searches the directory checkpoints for the last checkpoint for the given 
  | parameters and returns the directory name if such a directory exists. 
  | Returns empty string if no checkpoint directory was found.
  | Currently supports only the seed argument.


* | **remove_cwd**

   * | paths: List[str]
   * | return type: List[str]


  | Removes all directory names items pointing to the current working directory from a list and 
  | returns a new list.



