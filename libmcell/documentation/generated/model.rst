*****
Model
*****
Model
=====

This is the main class that is used to assemble all simulation input 
and configuration. It also provides methods to do initialization,
run simulation, and introspect the running simulation.

Attributes:
***********
* | **config**: Config = Config()
  | Simulation configuration.

* | **warnings**: Warnings = Warnings()
  | Configuration on how to report warnings.

* | **notifications**: Notifications = Notifications()
  | Configuration on how to report certain notifications.

* | **species**: List[Species] = None

* | **reaction_rules**: List[ReactionRule] = None

* | **surface_classes**: List[SurfaceClass] = None

* | **elementary_molecule_types**: List[ElementaryMoleculeType] = None
  | Used mainly when a BNGL file is loaded, if BNGL species is defined through 
  | Python API, this array is populated automatically

* | **release_sites**: List[ReleaseSite] = None
  | List of release sites to be included in the model.

* | **geometry_objects**: List[GeometryObject] = None
  | List of geometry objects to be included in the model.

* | **checkpointed_molecules**: List[BaseChkptMol] = None
  | Used when resuming simulation from a checkpoint.

* | **viz_outputs**: List[VizOutput] = None
  | List of visualization outputs to be included in the model.
  | There is usually just one VizOutput object.

* | **counts**: List[Count] = None
  | List of counts to be included in the model.


Methods:
*********
* | **initialize**

   * | print_copyright: bool = True
     | Prints information about MCell.


  | Initializes model, initialization blocks most of changes to 
  | contained components.


* | **run_iterations**

   * | iterations: float
     | Number of iterations to run. Value is truncated to an integer.

   * | return type: int


  | Runs specified number of iterations. Returns the number of iterations
  | executed (it might be less than the requested number of iterations when 
  | a checkpoint was scheduled).


* | **end_simulation**

   * | print_final_report: bool = True
     | Print information on simulation time and counts of selected events.


  | Generates the last visualization and reaction output (if they are included 
  | in the model), then flushes all buffers and optionally prints simulation report. 
  | Buffers are also flushed when the Model object is destroyed such as when Ctrl-C
  | is pressed during simulation.


* | **add_subsystem**

   * | subsystem: Subsystem

  | Adds all components of a Subsystem object to the model.


* | **add_instantiation**

   * | instantiation: Instantiation

  | Adds all components of an Instantiation object to the model.


* | **add_observables**

   * | observables: Observables

  | Adds all counts and viz outputs of an Observables object to the model.


* | **dump_internal_state**


  | Prints out the simulation engine's internal state, mainly for debugging.


* | **export_data_model**

   * | file: str = None

  | Exports the current state of the model into a data model JSON format.
  | Does not export state of molecules.
  | Must be called after model initialization.
  | If file is not set, then uses the first VizOutput to determine the target directory 
  | and creates name using the current iteration. Fails if argument file is not set and there is no VizOutput.
  | Always exports the current state, i.e. with the current geometry and reaction rates. 
  | Events (ReleaseSites and VizOutputs) with scheduled time other than zero are not exported correctly yet.


* | **export_viz_data_model**

   * | file: str = None
     | Optional path to the output data model file.


  | Same as export_data_model, only the created data model will contain only information required for visualization in CellBlender. This makes the loading ofthemodel by CellBlender faster and also allows to avoid potential compatibility issues.


* | **release_molecules**

   * | release_site: ReleaseSite

  | Performs immediate release of molecules based on the definition of the release site argument.
  | The ReleaseSite.release_time must not be in the past and must be within the current iteration 
  | meaning that the time must be greater or equal iteration \* time_step and less than (iteration + 1) \* time_step.
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
     | Object whose vertex will be changed.

   * | vertex_index: int
     | Index of vertex in object's vertex list that will be changed.

   * | displacement: List[float]
     | Change of vertex coordinates [x, y, z] (in um) that will be added to the current 
     | coordinates of the vertex.


  | Appends information about a displacement for given object's vertex into an internal list of vertex moves. 
  | To do the actual geometry change, call Model.apply_vertex_moves.
  | The reason why we first need to collect all changes and then apply them all at the same time is for performance
  | reasons.


* | **apply_vertex_moves**

   * | collect_wall_wall_hits: bool = False
     | When set to True, a list of wall pairs that collided is returned,
     | otherwise an empty list is returned.

   * | return type: List[WallWallHitInfo]


  | Applies all the vertex moves specified with Model.add_vertex_move call.
  | Walls of different objects are checked against collisions and move the maximal way so that they do not 
  | overlap.
  | Note\: It is not supported yet to move two objects that woudl collide at the same time.  
  | When collect_wall_wall_hits is True, a list of wall pairs that collided is returned,
  | when collect_wall_wall_hits is False, and empty list is returned.


* | **register_mol_wall_hit_callback**

   * | function: Callable, # std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)>
     | Callback function to be called. 
     | The function must have two arguments MolWallHitInfo and context.

   * | context: Any, # py::object
     | Context passed to the callback function, the callback function can store
     | information to this object. Some context must be always passed, even when 
     | it is a useless python object.

   * | object: GeometryObject = None
     | Only hits of this object will be reported, any object hit is reported when not set.

   * | species: Species = None
     | Only hits of molecules of this species will be reported, any species hit is reported when not set.


  | Register a callback for event when a molecule hits a wall. 
  | Note\: There can be currently only a single wall hit callback registered.


* | **register_reaction_callback**

   * | function: Callable, # std::function<void(std::shared_ptr<ReactionInfo>, py::object)>
     | Callback function to be called. 
     | The function must have two arguments ReactionInfo and context.
     | Called right after a reaction occured but before the reactants were removed.
     | After return the reaction proceeds and reactants are removed (unless they were kept
     | by the reaction such as with reaction A + B -> A + C).

   * | context: Any, # py::object
     | Context passed to the callback function, the callback function can store
     | information to this object. Some context must be always passed, even when 
     | it is a useless python object.

   * | reaction_rule: ReactionRule
     | The callback function will be called whenever this reaction rule is applied.


  | Defines a function to be called when a reaction was processed.
  | It is allowed to do state modifications except for removing reacting molecules, 
  | they will be removed automatically after return from this callback. 
  | Unlimited number of reaction callbacks is allowed.


* | **load_bngl**

   * | file_name: str
   * | observables_files_prefix: str = ''
     | Prefix to be used when creating files where observable values are stored during simulation.

   * | default_release_region: Region = None
     | Used as region for releases for seed species that have no compartments specified.

   * | parameter_overrides: Dict[str, float] = None
     | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
     | its value is ignored and instead value parameter_overrides[k] is used.


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
  | as a BNGL file that can be then loaded by MCell4 or BioNetGen. 
  | Note\: Limited currrently to exactly one volume compartment and volume reactions.


* | **save_checkpoint**

   * | custom_dir: str = None
     | Sets custom directory where the checkpoint will be stored. 
     | The default is 'checkpoints/seed_<SEED>/it_<ITERATION>'.


  | Saves current model state as checkpoint. 
  | The default directory structure is checkpoints/seed_<SEED>/it_<ITERATION>,
  | it can be changed by setting 'custom_dir'.
  | If used during an iteration such as in a callback, an event is scheduled for the  
  | beginning of the next iteration. This scheduled event saves the checkpoint.


* | **schedule_checkpoint**

   * | iteration: int = 0
     | Specifies iteration number when the checkpoint save will occur. 
     | Please note that iterations are counted from 0.
     | To schedule a checkpoint for the closest time as possible, keep the default value 0,
     | this will schedule checkpoint for the beginning of the iteration with number current iteration + 1.  
     | If calling schedule_checkpoint from a different thread (e.g. by using threading.Timer), 
     | it is highly recommended to keep the default value 0 or choose some time that will be 
     | for sure in the future.

   * | continue_simulation: bool = False
     | When false, saving the checkpoint means that we want to terminate the simulation 
     | right after the save. The currently running function Model.run_iterations
     | will not simulate any following iterations and execution will return from this function
     | to execute the next statement which is usually 'model.end_simulation()'.
     | When true, the checkpoint is saved and simulation continues uninterrupted.

   * | custom_dir: str = None
     | Sets custom directory where the checkpoint will be stored. 
     | The default is 'checkpoints/seed_<SEED>/it_<ITERATION>'.


  | Schedules checkpoint save event that will occur when an iteration is started.  
  | This means that it will be executed right before any other events scheduled for 
  | the given iteration are executed.
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

  | Adds a reference to the release site s to the list of release sites.


* | **find_release_site**

   * | name: str
   * | return type: ReleaseSite


  | Finds a release site by its name, returns None if no such release site is present.


* | **add_geometry_object**

   * | o: GeometryObject

  | Adds a reference to the geometry object o to the list of geometry objects.


* | **find_geometry_object**

   * | name: str
   * | return type: GeometryObject


  | Finds a geometry object by its name, returns None if no such geometry object is present.


* | **find_volume_compartment_object**

   * | name: str
   * | return type: GeometryObject


  | Finds a geometry object by its name, the geometry object must be a BNGL compartment.
  | Returns None if no such geometry object is present.


* | **find_surface_compartment_object**

   * | name: str
   * | return type: GeometryObject


  | Finds a geometry object that is a BNGL compartment and its surface name is name.
  | Returns None if no such geometry object is present.


* | **load_bngl_seed_species**

   * | file_name: str
     | Path to the BNGL file.

   * | default_release_region: Region = None
     | Used as region for releases for seed species that have no compartments specified.

   * | parameter_overrides: Dict[str, float] = None
     | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
     | its value is ignored and instead value parameter_overrides[k] is used.


  | Loads section seed species from a BNGL file and creates release sites according to it.
  | All elementary molecule types used in the seed species section must be already defined in subsystem.
  | If an item in the BNGL seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are then released into or onto the 
  | default_region. 
  | Does not create geometry objects. 
  | All compartments used in the loaded BNGL seed species section must exist in the model before 
  | model intialization.


* | **add_viz_output**

   * | viz_output: VizOutput

  | Adds a reference to the viz_output object to the list of visualization output specifications.


* | **add_count**

   * | count: Count

  | Adds a reference to the count object to the list of count specifications.


* | **find_count**

   * | name: str
   * | return type: Count


  | Finds a count object by its name, returns None if no such count is present.


* | **load_bngl_observables**

   * | file_name: str
     | Path to the BNGL file.

   * | output_files_prefix: str = ''
     | Prefix to be used when creating files with observable values.
     | The usual value is './react_data/seed_' + str(SEED).zfill(5) + '/'.

   * | parameter_overrides: Dict[str, float] = None
     | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
     | its value is ignored and instead value parameter_overrides[k] is used.


  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.


* | **get_molecule_ids**

   * | pattern: Complex = None
     | BNGL pattern to select molecules based on their species, might use compartments.

   * | return type: List[int]


  | Returns a list of ids of molecules.
  | If the arguments pattern is not set, the list of all molecule ids is returned.  
  | If the argument pattern is set, the list of all molecule ids whose species match 
  | the pattern is returned.


* | **get_molecule**

   * | id: int
     | Unique id of the molecule to be retrieved.

   * | return type: Molecule


  | Returns a information on a molecule from the simulated environment, 
  | None if the molecule does not exist.


* | **get_species_name**

   * | species_id: int
     | Id of the species.

   * | return type: str


  | Returns a string representing canonical species name in the BNGL format.


* | **get_vertex**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | return type: Vec3


  | Returns coordinates of a vertex.


* | **get_wall**

   * | object: GeometryObject
     | Geometry object whose wall to retrieve.

   * | wall_index: int
     | This is the index of the wall in the geometry object's walls (wall_list).

   * | return type: Wall


  | Returns information about a wall belonging to a given object.


* | **get_vertex_unit_normal**

   * | object: GeometryObject
     | Geometry object whose vertex to retrieve.

   * | vertex_index: int
     | This is the index of the vertex in the geometry object's vertex_list.

   * | return type: Vec3


  | Returns sum of all wall normals that use this vertex converted to a unit vector of length 1um.
  | This represents the unit vector pointing outwards from the vertex.


* | **get_wall_unit_normal**

   * | object: GeometryObject
     | Geometry object whose wall's normal to retrieve.

   * | wall_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | return type: Vec3


  | Returns wall normal converted to a unit vector of length 1um.


