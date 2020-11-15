**************************
MCell 4 Python API Classes
**************************
Complex
=======

This class represents a complex molecule composed of molecule instances.
It is either defined using a BNGL string or using a list of elementary molecule instances.
On top of that, orientation may be defined.
This class is used as argument in cases where either a fully qualified instance or a pattern 
can be provided such as in observable Count.  


Attributes
**********
* | **name**: str = None
  | When set, this complex instance is initialized from a BNGL string passed as this argument, 
  | the string is parsed during model initialization so the molecule types it uses
  | don't have to be defined before initialization.
  | 

* | **elementary_molecule_instances**: List[ElementaryMoleculeInstance] = None
  | Individual molecule instances contained in the complex.

* | **orientation**: Orientation = Orientation.DEFAULT
  | Specifies orientation of a molecule. 
  | When Orientation.DEFAULT if kept then during model initialization is
  | 'orientation' set to Orientation.NONE for volume complexes and to 
  | Orientation.UP for surface complexes.
  | Ignored by derived class Species.
  | 

* | **compartment_name**: str = None
  | Specifies compartment name of this Complex. 
  | Only one of 'orientation' and 'compartment_name' can be set. 
  | If a 2D/surface compartment is specified, the complex must be a surface complex and 
  | orientation is set to Orientation.UP.
  | If a 3D/volume compartment is specified, the complex must be a volume complex and
  | orientation is set to Orientation.NONE. 
  | Ignored by derived class Species.
  | 


Methods
********
* | **to_bngl_str**
  | Creates a string that corresponds to its BNGL representation

   * return type: str

* | **as_species**
  | Returns a Species object based on this Complex. All species-specific 
  | attributes are set to their default values and 'name' is set to value returned by 
  | 'to_bngl_str()'.
  | 

   * return type: Species


ComponentInstance
=================

Instance of a component belonging to a molecule instance.
A component instance may have its state set.
It is also used to connect molecule instance in a complex instance.


Attributes
**********
* | **component_type**: ComponentType

* | **state**: str = STATE_UNSET

* | **bond**: int = BOND_UNBOUND


Methods
********
* | **to_bngl_str**
  | Creates a string that corresponds to its BNGL representation

   * return type: str


ComponentType
=============

Attributes
**********
* | **name**: str

* | **states**: List[str] = None


Methods
********
* | **inst**
   * return type: ComponentInstance
   * | state: str = STATE_UNSET   * | bond: int = BOND_UNBOUND
* | **inst**
   * return type: ComponentInstance
   * | state: int = STATE_UNSET_INT   * | bond: int = BOND_UNBOUND

Config
======

Attributes
**********
* | **seed**: int = 1

* | **time_step**: float = 1e-6
  | Default value is 1us, in seconds

* | **surface_grid_density**: float = 10000

* | **interaction_radius**: float = None
  | Diffusing volume molecules will interact with each other when
  | they get within N microns of each other. The default is
  | 1/sqrt(PI * Sigma_s) where Sigma_s is the surface grid density 
  | (default or userspecified)
  | 

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
  | 

* | **center_molecules_on_grid**: bool = False

* | **initial_partition_origin**: List[float] = None
  | Optional placement of the partition 0 placement, specifies the left, lower and front 
  | point. If not set, value -partition_dimension/2 is used for each of the dimensions 
  | placing the center of the partition to (0, 0, 0).   
  | 

* | **partition_dimension**: float = 10

* | **subpartition_dimension**: float = 0.5

* | **total_iterations_hint**: float = 1000000
  | Estimated value of total iterations, used when generating visualization data 
  | files and also for other reporting uses. Value is truncated to an integer.
  | 

* | **check_overlapped_walls**: bool = True
  | Enables check for overlapped walls. Overlapping walls can cause issues during 
  | simulation such as a molecule escaping closed geometry when it hits two walls 
  | that overlap. 
  | 


Methods
********

Count
=====

Attributes
**********
* | **file_name**: str

* | **count_expression**: CountTerm = None
  | The count expression must be composed only from CountTerm objects that are added or 
  | subtracted.  
  | 

* | **multiplier**: float = 1
  | In some cases it might be useful to multiply the whole count by a constant to get 
  | for instance concentration. The count_expression is not an arbitrary expression
  | and such multiplication can be done through this attribute.  
  | 

* | **every_n_timesteps**: float = 1
  | Value is truncated (floored) to an integer.
  | 

* | **species_pattern**: Complex = None
  | Count the number of molecules that match the given complex instance pattern.
  | Counts each molecule exactly once. 
  | If the pattern has a compartment set, this specifies the counted region.  
  | 

* | **molecules_pattern**: Complex = None
  | Count the number of matches of the given pattern on molecules.
  | The observable will count a molecule every time it matches the pattern.
  | When the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a 
  | molecule A(a!1).A(a!1,b!2).B(a!2) will be counted twice because the 
  | pattern may match in two different ways. 
  | If the pattern has a compartment set, this specifies the counted region.  
  | 

* | **reaction_rule**: ReactionRule = None

* | **region**: Region = None
  | Only a GeometryObject or SurfaceRegion can be passed as the region argument, 
  | compound regions (created with +, -, *) are not supproted yet.   
  | Cannot be set when 'species_pattern' or 'molecules_pattern' has a  
  | compartment specified.
  | If pattern compartment is not specified and 'region' is left 'unset', 
  | counting is done in the whole world.
  | 

* | **node_type**: ExprNodeType = ExprNodeType.LEAF
  | Internal, used to represent an expression

* | **left_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods
********
* | **__add__**
   * return type: CountTerm
   * | op2: CountTerm
* | **__sub__**
   * return type: CountTerm
   * | op2: CountTerm

CountTerm
=========

Attributes
**********
* | **species_pattern**: Complex = None
  | Count the number of molecules that match the given complex instance pattern.
  | Counts each molecule exactly once. 
  | If the pattern has a compartment set, this specifies the counted region.  
  | 

* | **molecules_pattern**: Complex = None
  | Count the number of matches of the given pattern on molecules.
  | The observable will count a molecule every time it matches the pattern.
  | When the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a 
  | molecule A(a!1).A(a!1,b!2).B(a!2) will be counted twice because the 
  | pattern may match in two different ways. 
  | If the pattern has a compartment set, this specifies the counted region.  
  | 

* | **reaction_rule**: ReactionRule = None

* | **region**: Region = None
  | Only a GeometryObject or SurfaceRegion can be passed as the region argument, 
  | compound regions (created with +, -, *) are not supproted yet.   
  | Cannot be set when 'species_pattern' or 'molecules_pattern' has a  
  | compartment specified.
  | If pattern compartment is not specified and 'region' is left 'unset', 
  | counting is done in the whole world.
  | 

* | **node_type**: ExprNodeType = ExprNodeType.LEAF
  | Internal, used to represent an expression

* | **left_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods
********
* | **__add__**
   * return type: CountTerm
   * | op2: CountTerm
* | **__sub__**
   * return type: CountTerm
   * | op2: CountTerm

ElementaryMoleculeInstance
==========================

Attributes
**********
* | **elementary_molecule_type**: ElementaryMoleculeType

* | **components**: List[ComponentInstance] = None


Methods
********
* | **to_bngl_str**
  | Creates a string that corresponds to its BNGL representation

   * return type: str


ElementaryMoleculeType
======================

Attributes
**********
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
  | 


Methods
********
* | **inst**
   * return type: ElementaryMoleculeInstance
   * | components: List[ComponentInstance] = None

GeometryObject
==============

Attributes
**********
* | **name**: str
  | Name of the object. Also represents BNGL compartment name if 'is_bngl_compartment' is True.
  | 

* | **vertex_list**: List[List[float]]
  | List of [x,y,z] triplets specifying positions of individual vertices.
  | Equivalent to List[Vec3] however, defining a constructor Vec3(List[float]) then 
  | tries to convert all lists of floats to Vec3 
  |   
  | 

* | **wall_list**: List[List[int]]
  | List of [a,b,c] triplets specifying each wall, individual values are indices into the vertex list.
  | Equivalent to List[IVec3]. 
  | 

* | **is_bngl_compartment**: bool = False

* | **surface_compartment_name**: str = None

* | **surface_regions**: List[SurfaceRegion] = None

* | **surface_class**: SurfaceClass = None
  | Surface class for the whole object's surface. It is applied to the whole surface of this object 
  | except for those surface regions that have their specific surface class set explicitly.
  | 

* | **initial_surface_releases**: List[InitialSurfaceRelease] = None
  | Equivalent to MDL's MODIFY_SURFACE_REGIONS/MOLECULE_DENSITY or MOLECULE_NUMBER,
  | each item defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.
  | 

* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.
  | 

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods
********
* | **translate**
  | Move object by a specified vector, must be done before model initialization.

   * | move: Vec3
* | **__add__**
  | Computes union of thwo regions

   * return type: Region
   * | other: Region
* | **__sub__**
   * return type: Region
   * | other: Region
* | **__mul__**
   * return type: Region
   * | other: Region

InitialSurfaceRelease
=====================

Defines molecules to be released onto a SurfaceRegion right when simulation starts

Attributes
**********
* | **complex**: Complex

* | **number_to_release**: int = None
  | Number of molecules to be released onto a region,
  | only one of number_to_release and density can be set.
  | 

* | **density**: float = None
  | Density of molecules to be released onto a region,
  | only one of number_to_release and density can be set.
  | 


Methods
********

InstantiationData
=================

Attributes
**********
* | **release_sites**: List[ReleaseSite] = None

* | **geometry_objects**: List[GeometryObject] = None


Methods
********
* | **add_release_site**
  | Makes a copy of the release site

   * | s: ReleaseSite
* | **find_release_site**
   * return type: ReleaseSite
   * | name: str
* | **add_geometry_object**
  | Makes a copy of the geometry object, in the future we will probably add some transformations

   * | o: GeometryObject
* | **find_geometry_object**
   * return type: GeometryObject
   * | name: str
* | **find_volume_compartment**
   * return type: GeometryObject
   * | name: str
* | **find_surface_compartment**
   * return type: GeometryObject
   * | name: str
* | **load_bngl_seed_species**
  | Loads section seed species from a BNGL file and creates release sites according to it.
  | All elementary molecule types used in the seed species section must be already defined in subsystem.
  | If an item in the BNGL seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region. 
  | 

   * | file_name: str   * | subsystem: Subsystem   * | default_release_region: Region = None     | Loads section seed species from a BNGL file and creates release sites according to it.
     | All elementary molecule types used in the seed species section must be already defined in subsystem.
     | If an item in the BNGL seed species section does not have its compartment set,
     | the argument default_region must be set and the molecules are released into or onto the 
     | default_region. 
     | 

   * | parameter_overrides: Dict[str, float] = None     | Loads section seed species from a BNGL file and creates release sites according to it.
     | All elementary molecule types used in the seed species section must be already defined in subsystem.
     | If an item in the BNGL seed species section does not have its compartment set,
     | the argument default_region must be set and the molecules are released into or onto the 
     | default_region. 
     | 



Model
=====

Attributes
**********
* | **config**: Config = Config()

* | **warnings**: Warnings = Warnings()

* | **notifications**: Notifications = Notifications()

* | **species**: List[Species] = None

* | **reaction_rules**: List[ReactionRule] = None

* | **surface_classes**: List[SurfaceClass] = None

* | **elementary_molecule_types**: List[ElementaryMoleculeType] = None
  | Used mainly when a BNGL file is loaded, if BNGL species is defined through 
  | Python API, this array is populated automatically 
  | 

* | **release_sites**: List[ReleaseSite] = None

* | **geometry_objects**: List[GeometryObject] = None

* | **viz_outputs**: List[VizOutput] = None

* | **counts**: List[Count] = None


Methods
********
* | **initialize**

* | **run_iterations**
   * | iterations: float
* | **end_simulation**
  | Generates the last visualization and reaction output (if they were defined), then
  | flushes all buffers and optionally prints simulation report. 
  | Buffers are also flushed when the Model object is destroyed.   
  | 

   * | print_final_report: bool = True     | Generates the last visualization and reaction output (if they were defined), then
     | flushes all buffers and optionally prints simulation report. 
     | Buffers are also flushed when the Model object is destroyed.   
     | 


* | **add_subsystem**
   * | subsystem: Subsystem
* | **add_instantiation_data**
   * | instantiation_data: InstantiationData
* | **add_observables**
   * | observables: Observables
* | **dump_internal_state**
  | Prints out the simulation engine's internal state, mainly for debugging.


* | **export_data_model**
  | If file is not set, then uses the first VizOutput to determine the target directory 
  | and creates name using the current iteration. Fails if argument file is not set and there is no VizOutput.
  | Must be called after initialization.
  | Always exports the current state, i.e. with the current . 
  | Events (ReleaseSites and VizOutputs) with scheduled time other than zero cannot be imported correectly yet.  
  | 

   * | file: str = None     | If file is not set, then uses the first VizOutput to determine the target directory 
     | and creates name using the current iteration. Fails if argument file is not set and there is no VizOutput.
     | Must be called after initialization.
     | Always exports the current state, i.e. with the current . 
     | Events (ReleaseSites and VizOutputs) with scheduled time other than zero cannot be imported correectly yet.  
     | 


* | **export_viz_data_model**
  | Same as export_data_model, only the created data model will contain only information required for visualization in CellBlender. This makes the loading ofthemodel by CellBlender faster and also allows to avoid potential compatibility issues.

   * | file: str = None     | Same as export_data_model, only the created data model will contain only information required for visualization in CellBlender. This makes the loading ofthemodel by CellBlender faster and also allows to avoid potential compatibility issues.


* | **get_molecule_ids**
  | Returns a list of ids of molecules of given Species existing in the simulated environment,
  | if the argument species is not set, returns list of all molecules.      
  | 

   * return type: List[int]
   * | species: Species = None     | Returns a list of ids of molecules of given Species existing in the simulated environment,
     | if the argument species is not set, returns list of all molecules.      
     | 


* | **get_molecule**
  | Returns a molecule from the simulated environment, None if the molecule does not exist

   * return type: Molecule
   * | id: int
* | **get_vertex**
  | Returns coordinates of a vertex.

   * return type: Vec3
   * | object: GeometryObject   * | vertex_index: int
* | **get_wall**
  | Returns information about a wall belonging to a given object.

   * return type: Wall
   * | object: GeometryObject   * | wall_index: int
* | **get_vertex_unit_normal**
  | Returns sum of all wall normals that use this vertex converted to a unit vector of length 1um.
  | This represents the unit vector pointing outwards from the vertex.
  | 

   * return type: Vec3
   * | object: GeometryObject   * | vertex_index: int
* | **get_wall_unit_normal**
  | Returns wall normal converted to a unit vector of length 1um.

   * return type: Vec3
   * | object: GeometryObject   * | wall_index: int
* | **add_vertex_move**
  | Adds a displacement for given object's vertex, only stored until apply_vertex_moves is called

   * | object: GeometryObject   * | vertex_index: int   * | displacement: Vec3
* | **apply_vertex_moves**
  | Applies all the vertex moves specified with add_vertex_move call.
  | Walls of different objects are checked against collisions and move the maximal way so that they do not 
  | overlap. (the current pllementation is a bit basic and may not work 100% correctly) 
  | When collect_wall_wall_hits is True, a list of wall pairs that collided is returned,
  | when collect_wall_wall_hits is False, and empty list is returned.
  | 

   * return type: List[WallWallHitInfo]
   * | collect_wall_wall_hits: bool = False     | Applies all the vertex moves specified with add_vertex_move call.
     | Walls of different objects are checked against collisions and move the maximal way so that they do not 
     | overlap. (the current pllementation is a bit basic and may not work 100% correctly) 
     | When collect_wall_wall_hits is True, a list of wall pairs that collided is returned,
     | when collect_wall_wall_hits is False, and empty list is returned.
     | 


* | **register_mol_wall_hit_callback**
  | There can be currently only a single wall hit callback registered.

   * | function: std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)>   * | context: py::object   * | object: GeometryObject = None     | There can be currently only a single wall hit callback registered.

   * | species: Species = None     | There can be currently only a single wall hit callback registered.


* | **load_bngl**
  | Loads sections\: molecule types, reaction rules, seed species, and observables from a BNGL file
  | and creates objects in the current model according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.
  | If an item in the seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region. 
  | 

   * | file_name: str   * | observables_files_prefix: str =      | Loads sections\: molecule types, reaction rules, seed species, and observables from a BNGL file
     | and creates objects in the current model according to it.
     | All elementary molecule types used in the seed species section must be defined in subsystem.
     | If an item in the seed species section does not have its compartment set,
     | the argument default_region must be set and the molecules are released into or onto the 
     | default_region. 
     | 

   * | default_release_region: Region = None     | Loads sections\: molecule types, reaction rules, seed species, and observables from a BNGL file
     | and creates objects in the current model according to it.
     | All elementary molecule types used in the seed species section must be defined in subsystem.
     | If an item in the seed species section does not have its compartment set,
     | the argument default_region must be set and the molecules are released into or onto the 
     | default_region. 
     | 

   * | parameter_overrides: Dict[str, float] = None     | Loads sections\: molecule types, reaction rules, seed species, and observables from a BNGL file
     | and creates objects in the current model according to it.
     | All elementary molecule types used in the seed species section must be defined in subsystem.
     | If an item in the seed species section does not have its compartment set,
     | the argument default_region must be set and the molecules are released into or onto the 
     | default_region. 
     | 


* | **export_to_bngl**
  | Exports all defined species, reaction rules and applicable observables
  | as a BNGL file. 
  | Limited currrently to exactly one volume compartment and volume reactions.
  | 

   * | file_name: str
* | **add_species**
   * | s: Species
* | **find_species**
   * return type: Species
   * | name: str
* | **add_reaction_rule**
   * | r: ReactionRule
* | **find_reaction_rule**
   * return type: ReactionRule
   * | name: str
* | **add_surface_class**
   * | sc: SurfaceClass
* | **find_surface_class**
   * return type: SurfaceClass
   * | name: str
* | **add_elementary_molecule_type**
   * | mt: ElementaryMoleculeType
* | **find_elementary_molecule_type**
   * return type: ElementaryMoleculeType
   * | name: str
* | **load_bngl_molecule_types_and_reaction_rules**
  | Parses a BNGL file and only reads molecule types and
  | reaction rules sections, e.g. ignores observables. 
  | Parameter values are evaluated and the result value 
  | is directly used.  
  | Compartments names are stored in rxn rules as strings because
  | compartments belong to geometry objects and the subsystem is independent
  | on specific geometry.
  | However they must be defined on initialization.
  |  
  | 

   * | file_name: str   * | parameter_overrides: Dict[str, float] = None     | Parses a BNGL file and only reads molecule types and
     | reaction rules sections, e.g. ignores observables. 
     | Parameter values are evaluated and the result value 
     | is directly used.  
     | Compartments names are stored in rxn rules as strings because
     | compartments belong to geometry objects and the subsystem is independent
     | on specific geometry.
     | However they must be defined on initialization.
     |  
     | 


* | **add_release_site**
  | Makes a copy of the release site

   * | s: ReleaseSite
* | **find_release_site**
   * return type: ReleaseSite
   * | name: str
* | **add_geometry_object**
  | Makes a copy of the geometry object, in the future we will probably add some transformations

   * | o: GeometryObject
* | **find_geometry_object**
   * return type: GeometryObject
   * | name: str
* | **find_volume_compartment**
   * return type: GeometryObject
   * | name: str
* | **find_surface_compartment**
   * return type: GeometryObject
   * | name: str
* | **load_bngl_seed_species**
  | Loads section seed species from a BNGL file and creates release sites according to it.
  | All elementary molecule types used in the seed species section must be already defined in subsystem.
  | If an item in the BNGL seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region. 
  | 

   * | file_name: str   * | subsystem: Subsystem   * | default_release_region: Region = None     | Loads section seed species from a BNGL file and creates release sites according to it.
     | All elementary molecule types used in the seed species section must be already defined in subsystem.
     | If an item in the BNGL seed species section does not have its compartment set,
     | the argument default_region must be set and the molecules are released into or onto the 
     | default_region. 
     | 

   * | parameter_overrides: Dict[str, float] = None     | Loads section seed species from a BNGL file and creates release sites according to it.
     | All elementary molecule types used in the seed species section must be already defined in subsystem.
     | If an item in the BNGL seed species section does not have its compartment set,
     | the argument default_region must be set and the molecules are released into or onto the 
     | default_region. 
     | 


* | **add_viz_output**
   * | viz_output: VizOutput
* | **add_count**
   * | count: Count
* | **load_bngl_observables**
  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.
  | 

   * | file_name: str   * | subsystem: Subsystem   * | output_files_prefix: str =      | Loads section observables from a BNGL file and creates Count objects according to it.
     | All elementary molecule types used in the seed species section must be defined in subsystem.
     | 

   * | parameter_overrides: Dict[str, float] = None     | Loads section observables from a BNGL file and creates Count objects according to it.
     | All elementary molecule types used in the seed species section must be defined in subsystem.
     | 



MolWallHitInfo
==============

Attributes
**********
* | **molecule_id**: int

* | **geometry_object**: GeometryObject
  | Object that was hit.

* | **wall_index**: int
  | Index of wall belonging to the geometry_object.

* | **time**: float
  | Time of the hit

* | **pos**: Vec3
  | Position of the hit

* | **time_before_hit**: float
  | Time when the molecule started to diffuse towards the hit wall. 
  | It is either the start of the molecule's diffusion or 
  | if a wall was hit later then the time of last wall hit.
  |   
  | 

* | **pos_before_hit**: Vec3
  | Position of the molecule at time_before_hit


Methods
********

Molecule
========

This is a Python representation of a molecule obtained from Model 
during simulation.


Attributes
**********
* | **id**: int = MOLECULE_ID_INVALID
  | Unique id of this molecule

* | **species**: Species = None

* | **pos3d**: Vec3 = None
  | TODO - Right now, contains only position of this is a volume molecule 
  | Contains position in space both for surface and volume molecules,
  | it won't be possible to change it for surface molecules.
  |  
  | 

* | **orientation**: Orientation = Orientation.NOT_SET
  | Contains orientation for surface molecule. Volume molecules 
  | have always orientation set to Orientation.NONE.
  | 


Methods
********
* | **remove**
  | Removes this molecule from simulation. Any subsequent modifications
  | of this object won't have any effect.  
  | 



MoleculeReleaseInfo
===================

Attributes
**********
* | **complex**: Complex
  | Complex instance defining the molecule that will be released.
  | Orientation of the complex instance is used to define orientation of the released molecule,
  | when Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and
  | surface molecules are released with Orientation.UP.
  | Compartment must not be set because this specific release definition states the location.  
  | 

* | **location**: List[float]
  | 3D position where the molecule will be released. 
  | If a molecule has a 2D diffusion constant, it will be
  | placed on the surface closest to the coordinate given. 
  | Argument must have exactly three floating point values.
  |   
  | 


Methods
********

Notifications
=============

Attributes
**********
* | **bng_verbosity_level**: int = 0
  | Sets verbosity level that enables printouts of extra information on BioNetGen 
  | species and rules created and used during simulation.
  | 

* | **rxn_and_species_report**: bool = True
  | Simulation generates files rxn_report_SEED.txt species_report_SEED.txt that contain
  | details on reaction classes and species that were created based on reaction rules.   
  | 

* | **probability_report**: bool = True

* | **diffusion_constant_report**: Notification = Notification.BRIEF

* | **final_summary**: bool = True

* | **iteration_report**: bool = True

* | **varying_probability_report**: bool = True
  | Related to changing rxn probabilities at runtime

* | **progress_report**: bool = True

* | **release_event_report**: bool = True

* | **molecule_collision_report**: bool = True


Methods
********

Observables
===========

Neither VizOutput, nor Count have name, therefore there are no find_* methods.


Attributes
**********
* | **viz_outputs**: List[VizOutput] = None

* | **counts**: List[Count] = None


Methods
********
* | **add_viz_output**
   * | viz_output: VizOutput
* | **add_count**
   * | count: Count
* | **load_bngl_observables**
  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.
  | 

   * | file_name: str   * | subsystem: Subsystem   * | output_files_prefix: str =      | Loads section observables from a BNGL file and creates Count objects according to it.
     | All elementary molecule types used in the seed species section must be defined in subsystem.
     | 

   * | parameter_overrides: Dict[str, float] = None     | Loads section observables from a BNGL file and creates Count objects according to it.
     | All elementary molecule types used in the seed species section must be defined in subsystem.
     | 



ReactionRule
============

Attributes
**********
* | **name**: str = None
  | Name of the reaction. If this is a reversible reaction, then it is the name of the 
  | reaction in forward direction.
  | 

* | **reactants**: List[Complex] = None

* | **products**: List[Complex] = None

* | **fwd_rate**: float = None

* | **rev_name**: str = None
  | Name of the reaction in reverse direction.

* | **rev_rate**: float = None
  | Reverse reactions rate, reaction is unidirectional when not specified

* | **variable_rate**: List[List[float]] = None
  | Variable rate is applicable only for irreversible reactions. Members fwd_rate and rev_rate 
  | must not be set. The array passed as this argument must have as its items a pair of floats (time, rate).     
  | 


Methods
********
* | **to_bngl_str**
  | Creates a string that corresponds to the reaction rule's BNGL representation, does not contain rates.

   * return type: str


Region
======

Represents region construted from 1 or more multiple, usually unnamed?

Attributes
**********
* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.
  | 

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods
********
* | **__add__**
  | Computes union of thwo regions

   * return type: Region
   * | other: Region
* | **__sub__**
   * return type: Region
   * | other: Region
* | **__mul__**
   * return type: Region
   * | other: Region

ReleasePattern
==============

Attributes
**********
* | **name**: str = None
  | Name of the release pattern

* | **release_interval**: float = TIME_INFINITY
  | During a train of releases, release molecules after every t seconds. 
  | Default is to release only once.
  | 

* | **train_duration**: float = TIME_INFINITY
  | The train of releases lasts for t seconds before turning off. 
  | Default is to never turn off.
  | 

* | **train_interval**: float = TIME_INFINITY
  | A new train of releases happens every t seconds. 
  | Default is to never have a new train. 
  | The train interval must not be shorter than the train duration.
  | 

* | **number_of_trains**: int = 1
  | Repeat the release process for n trains of releases. Default is one train.
  | For unlimited number of trains use constant NUMBER_OF_TRAINS_UNLIMITED.
  | 


Methods
********

ReleaseSite
===========

Attributes
**********
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
  | 

* | **molecule_list**: List[MoleculeReleaseInfo] = None
  | Used for LIST shape release mode. 
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | 

* | **release_time**: float = 0
  | Specifies time in seconds when the release event is executed.
  | In case when a release pattern is used, this is the time of the first release.      
  | Equivalent to MDL's RELEASE_PATTERN command DELAY.
  | 

* | **release_pattern**: ReleasePattern = None
  | Use the release pattern to define schedule of releases. 
  | The default is to release the specified number of molecules at the set release_time. 
  | 

* | **shape**: Shape = Shape.UNSET
  | Set automatically when 
  | 

* | **region**: Region = None
  | Sets shape to Shape.REGION_EXPR.

* | **location**: Vec3 = None

* | **site_diameter**: float = 0
  | For a geometrical release site, this releases molecules uniformly within
  | a radius r. Not used for releases on regions.
  | Usually required for Shape.List type of releases.
  | 

* | **site_radius**: float = None
  | For a geometrical release site, this releases molecules uniformly within
  | a radius r. Not used for releases on regions.
  | 

* | **number_to_release**: float = None
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Value is truncated (floored) to an integer.
  | 

* | **density**: float = None
  | Only one of number_to_release, density, concentration or molecule_list can be set.

* | **concentration**: float = None
  | Only one of number_to_release, density, concentration or molecule_list can be set.

* | **release_probability**: float = None


Methods
********

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
inherited attribute 'elementary_molecule_instances' from Complex
is used as a definition of the complex and this gives information 
on diffusion constants of the used elementary molecules.
Example\: m.Species(elementary_molecule_instances=[ei1, ei2]). 
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


Attributes
**********
* | **name**: str = None
  | Name of the species in the BNGL format. 
  | One must either specify 'name' or 'elementary_molecule_instances' 
  | (inherited from Complex). This argument 'name' is parsed during model 
  | initialization.    
  | 

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
  | 

* | **name**: str = None
  | When set, this complex instance is initialized from a BNGL string passed as this argument, 
  | the string is parsed during model initialization so the molecule types it uses
  | don't have to be defined before initialization.
  | 

* | **elementary_molecule_instances**: List[ElementaryMoleculeInstance] = None
  | Individual molecule instances contained in the complex.

* | **orientation**: Orientation = Orientation.DEFAULT
  | Specifies orientation of a molecule. 
  | When Orientation.DEFAULT if kept then during model initialization is
  | 'orientation' set to Orientation.NONE for volume complexes and to 
  | Orientation.UP for surface complexes.
  | Ignored by derived class Species.
  | 

* | **compartment_name**: str = None
  | Specifies compartment name of this Complex. 
  | Only one of 'orientation' and 'compartment_name' can be set. 
  | If a 2D/surface compartment is specified, the complex must be a surface complex and 
  | orientation is set to Orientation.UP.
  | If a 3D/volume compartment is specified, the complex must be a volume complex and
  | orientation is set to Orientation.NONE. 
  | Ignored by derived class Species.
  | 


Methods
********
* | **inst**
  | Creates a Complex of this Species. Can be currently used only for simple species, i.e. those that
  | have a single molecule instance and no components.
  | 

   * return type: Complex
   * | orientation: Orientation = Orientation.DEFAULT     | Creates a Complex of this Species. Can be currently used only for simple species, i.e. those that
     | have a single molecule instance and no components.
     | 

   * | compartment_name: str = None     | Creates a Complex of this Species. Can be currently used only for simple species, i.e. those that
     | have a single molecule instance and no components.
     | 


* | **to_bngl_str**
  | Creates a string that corresponds to its BNGL representation

   * return type: str

* | **as_species**
  | Returns a Species object based on this Complex. All species-specific 
  | attributes are set to their default values and 'name' is set to value returned by 
  | 'to_bngl_str()'.
  | 

   * return type: Species


Subsystem
=========

Attributes
**********
* | **species**: List[Species] = None

* | **reaction_rules**: List[ReactionRule] = None

* | **surface_classes**: List[SurfaceClass] = None

* | **elementary_molecule_types**: List[ElementaryMoleculeType] = None
  | Used mainly when a BNGL file is loaded, if BNGL species is defined through 
  | Python API, this array is populated automatically 
  | 


Methods
********
* | **add_species**
   * | s: Species
* | **find_species**
   * return type: Species
   * | name: str
* | **add_reaction_rule**
   * | r: ReactionRule
* | **find_reaction_rule**
   * return type: ReactionRule
   * | name: str
* | **add_surface_class**
   * | sc: SurfaceClass
* | **find_surface_class**
   * return type: SurfaceClass
   * | name: str
* | **add_elementary_molecule_type**
   * | mt: ElementaryMoleculeType
* | **find_elementary_molecule_type**
   * return type: ElementaryMoleculeType
   * | name: str
* | **load_bngl_molecule_types_and_reaction_rules**
  | Parses a BNGL file and only reads molecule types and
  | reaction rules sections, e.g. ignores observables. 
  | Parameter values are evaluated and the result value 
  | is directly used.  
  | Compartments names are stored in rxn rules as strings because
  | compartments belong to geometry objects and the subsystem is independent
  | on specific geometry.
  | However they must be defined on initialization.
  |  
  | 

   * | file_name: str   * | parameter_overrides: Dict[str, float] = None     | Parses a BNGL file and only reads molecule types and
     | reaction rules sections, e.g. ignores observables. 
     | Parameter values are evaluated and the result value 
     | is directly used.  
     | Compartments names are stored in rxn rules as strings because
     | compartments belong to geometry objects and the subsystem is independent
     | on specific geometry.
     | However they must be defined on initialization.
     |  
     | 



SurfaceClass
============

Defining a surface class allows surfaces to behave like species (in a limited way).

Attributes
**********
* | **name**: str
  | Name of the surface class

* | **properties**: List[SurfaceProperty] = None
  | A surface class can either have a list of properties or just one property.
  | In the usual case of having one property, one can set the attributes 
  | type, affected_species, etc. inherited from SurfaceProperty directly.
  | 

* | **type**: SurfacePropertyType = SurfacePropertyType.UNSET
  | Must be set.

* | **affected_complex_pattern**: Complex = None
  | A complex pattern with optional orientation must be set.
  | Default orientation means that the pattern matches any orientation.
  | For concentration or flux clamp the orientation specifies on which side  
  | will be the concentration held 
  | (UP is front or outside, DOWN is back or inside, and DEFAULT, ANY or NONE is on both sides).
  | The complex pattern must not have any compartment.
  | 

* | **concentration**: float = None
  | Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or 
  | SurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side 
  | of the wall that has this concentration or flux clamped.
  | 


Methods
********

SurfaceProperty
===============

Attributes
**********
* | **type**: SurfacePropertyType = SurfacePropertyType.UNSET
  | Must be set.

* | **affected_complex_pattern**: Complex = None
  | A complex pattern with optional orientation must be set.
  | Default orientation means that the pattern matches any orientation.
  | For concentration or flux clamp the orientation specifies on which side  
  | will be the concentration held 
  | (UP is front or outside, DOWN is back or inside, and DEFAULT, ANY or NONE is on both sides).
  | The complex pattern must not have any compartment.
  | 

* | **concentration**: float = None
  | Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or 
  | SurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side 
  | of the wall that has this concentration or flux clamped.
  | 


Methods
********

SurfaceRegion
=============

Surface region  in MDL, however a new class Region was instroduced in MCell4 so it was renamed 
to avoid confusion.


Attributes
**********
* | **name**: str

* | **wall_indices**: List[int]
  | Surface region must be a part of a GeometryObject, items in this list are indices to 
  | its wall_list array
  | 

* | **surface_class**: SurfaceClass = None
  | Has higher priority than the parent geometry object's surface class.

* | **initial_surface_releases**: List[InitialSurfaceRelease] = None
  | Equivalent to MDL's MODIFY_SURFACE_REGIONS/MOLECULE_DENSITY or MOLECULE_NUMBER,
  | each item defines either density or number of molecules to be released on this surface 
  | regions when simulation starts. 
  | 

* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.
  | 

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods
********
* | **__add__**
  | Computes union of thwo regions

   * return type: Region
   * | other: Region
* | **__sub__**
   * return type: Region
   * | other: Region
* | **__mul__**
   * return type: Region
   * | other: Region

VizOutput
=========

Attributes
**********
* | **output_files_prefix**: str

* | **species_list**: List[Species] = None
  | When empty and all_species is false, empty files will be generated

* | **all_species**: bool = False
  | Visualize all species

* | **mode**: VizMode = VizMode.ASCII

* | **every_n_timesteps**: float = 1
  | Value is truncated (floored) to an integer.
  | 


Methods
********

Wall
====

This is a Python representation of a molecule obtained from Model 
during simulation.


Attributes
**********
* | **geometry_object**: GeometryObject
  | Object to which this wall belongs.

* | **wall_index**: int
  | Index of this wall in the object to which this wall belongs.

* | **vertices**: List[Vec3]
  | Vertices of the triangle that represents this wall.

* | **area**: float

* | **normal**: Vec3
  | Normal of this wall, no guarantees on the length of this vector are given.
  | To get a unit vector, use Model.get_wall_unit_normal instead. 
  | 

* | **is_movable**: bool = True
  | If True, whis wall can be moved through Model.apply_vertex_moves,
  | if False, wall moves are ignored. 
  | 


Methods
********

WallWallHitInfo
===============

Attributes
**********
* | **wall1**: Wall

* | **wall2**: Wall


Methods
********

Warnings
========

Attributes
**********
* | **molecule_collision_report**: WarningLevel = WarningLevel.WARNING

* | **degenerate_polygons**: WarningLevel = WarningLevel.WARNING

* | **negative_diffusion_constant**: WarningLevel = WarningLevel.WARNING

* | **missing_surface_orientation**: WarningLevel = WarningLevel.ERROR

* | **negative_reaction_rate**: WarningLevel = WarningLevel.WARNING

* | **useless_volume_orientation**: WarningLevel = WarningLevel.WARNING

* | **high_reaction_probability**: WarningLevel = WarningLevel.IGNORE

* | **lifetime_too_short**: WarningLevel = WarningLevel.WARNING

* | **lifetime_threshold**: float = 50

* | **missed_reactions**: WarningLevel = WarningLevel.WARNING

* | **missed_reactions_threshold**: float = 0.00100000004749745


Methods
********

bngl_utils
==========

Attributes
**********

Methods
********
* | **load_bngl_parameters**
   * return type: Dict[str, float]
   * | file_name: str   * | parameter_overrides: Dict[str, float] = None

geometry_utils
==============

Attributes
**********

Methods
********
* | **create_box**
  | Creates a GeometryObject whose center is at (0, 0, 0).

   * return type: GeometryObject
   * | name: str   * | edge_length: float

