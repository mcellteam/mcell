.. _api-model:

*****
Model
*****
Model
=====

This is the main class that is used to assemble all simulation input 
and configuration. It also provides methods to do initialization,
run simulation, and introspect the running simulation.

Example: `1400_rel_site_for_each_it/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1400_rel_site_for_each_it/model.py>`_ 

Attributes:
***********
* | **config**: Config = Config()
  | Simulation configuration.

* | **warnings**: Warnings = Warnings()
  | Configuration on how to report warnings.

* | **notifications**: Notifications = Notifications()
  | Configuration on how to report certain notifications.

* | **species**: List[Species] = None
  | List of species to be included in the model for initialization.
  | Used usually only for simple species (species that are defined using a
  | single molecule type without components such as 'A').
  | Other species may be created inside simulation

* | **reaction_rules**: List[ReactionRule] = None

* | **surface_classes**: List[SurfaceClass] = None

* | **elementary_molecule_types**: List[ElementaryMoleculeType] = None
  | Contains list of elementary molecule types with their diffusion constants and other information. 
  | Populated when a BNGL file is loaded and also on initialization from Species objects present in 
  | the species list.

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

   * | with_geometry: bool = False
     | Include geometry in the dump.


  | Prints out the simulation engine's internal state, mainly for debugging.


* | **export_data_model**

   * | file: str = None
     | If file is not set, then uses the first VizOutput to determine the target directory 
     | and creates name using the current iteration. Fails if argument file is not set and 
     | there is no VizOutput in the model.


  | Exports the current state of the model into a data model JSON format.
  | Does not export state of molecules.
  | Must be called after model initialization.
  | Always exports the current state, i.e. with the current geometry and reaction rates. 
  | Events (ReleaseSites and VizOutputs) with scheduled time other than zero are not exported correctly yet.


* | **export_viz_data_model**

   * | file: str = None
     | Optional path to the output data model file.


  | Same as export_data_model, only the created data model will contain only information required for visualization
  | in CellBlender. This makes the loading of the model by CellBlender faster and also allows to avoid potential
  | compatibility issues.
  | Must be called after model initialization.

  | Example: `1520_sphere_collision/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1520_sphere_collision/model.py>`_ 


* | **export_geometry**

   * | output_files_prefix: str = None
     | Optional prefix for .obj and .mtl files that will be created on export. 
     | If output_files_prefix is not set, then uses the first VizOutput to determine the target directory 
     | and creates names using the current iteration. Fails if argument output_files_prefix is not set and 
     | there is no VizOutput in the model.


  | Exports model geometry as Wavefront OBJ format. 
  | Must be called after model initialization.
  | Does not export material colors (yet).


* | **release_molecules**

   * | release_site: ReleaseSite

  | Performs immediate release of molecules based on the definition of the release site argument.
  | The ReleaseSite.release_time must not be in the past and must be within the current iteration 
  | meaning that the time must be greater or equal iteration \* time_step and less than (iteration + 1) \* time_step.
  | The ReleaseEvent must not use a release_pattern because this is an immediate release and it is not 
  | scheduled into the global scheduler.

  | Example: `2300_immediate_release/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2300_immediate_release/model.py>`_ 


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

  | Example: `1850_run_unimol_rxn_in_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1850_run_unimol_rxn_in_callback/model.py>`_ 


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

  | Example: `1510_tetrahedron_box_collision_moving_3_verts/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1510_tetrahedron_box_collision_moving_3_verts/model.py>`_ 


* | **apply_vertex_moves**

   * | collect_wall_wall_hits: bool = False
     | When set to True, a list of wall pairs that collided is returned,
     | otherwise an empty list is returned.

   * | return type: List[WallWallHitInfo]


  | Applies all the vertex moves specified with Model.add_vertex_move call.
  | Walls of different objects are checked against collisions and move the maximal way so that they do not 
  | overlap.
  | The API representation (GeometryObject) is not updated, only the internal MCell data are changed.
  | Note\: It is not supported yet to move two objects that woudl collide at the same time.  
  | When collect_wall_wall_hits is True, a list of wall pairs that collided is returned,
  | when collect_wall_wall_hits is False, and empty list is returned.

  | Example: `1510_tetrahedron_box_collision_moving_3_verts/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1510_tetrahedron_box_collision_moving_3_verts/model.py>`_ 


* | **register_mol_wall_hit_callback**

   * | function: Callable, # std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)>
     | Callback function to be called. 
     | The function must have two arguments MolWallHitInfo and context.
     | Do not modify the received MolWallHitInfo object since it may be reused for other 
     | wall hit callbacks (e.g. when the first callback is for a specific geometry object and 
     | the second callback is for any geometry object). 
     | The context object (py::object type argument) is on the other hand provided 
     | to be modified and one can for instance use it to count the number of hits..

   * | context: Any, # py::object
     | Context passed to the callback function, the callback function can store
     | information to this object. Some context must be always passed, even when 
     | it is a useless python object.

   * | object: GeometryObject = None
     | Only hits of this object will be reported, any object hit is reported when not set.

   * | species: Species = None
     | Only hits of molecules of this species will be reported, any hit of volume molecules of 
     | any species is reported when this argument is not set.
     | Sets an internal flag for this species to make sure that the species id does not change 
     | during simulation.


  | Register a callback for event when a molecule hits a wall. 
  | May be called only after model initialization because it internally uses geometry object
  | and species ids that are set during the initialization.

  | Example: `1300_wall_hit_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1300_wall_hit_callback/model.py>`_ 


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
  | May be called only after model initialization because it internally uses 
  | reaction rule ids that are set during the initialization.

  | Example: `1800_vol_rxn_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1800_vol_rxn_callback/model.py>`_ 


* | **load_bngl**

   * | file_name: str
     | Path to the BNGL file to be loaded.

   * | observables_path_or_file: str = ''
     | Directory prefix or file name where observable values will be stored.
     | If a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty 
     | string is used,
     | each observable gets its own file and the output file format for created Count 
     | objects is CountOutputFormat.DAT.
     | If a file has a .gdat extension such as 
     | './react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat', all observable are stored in this 
     | file and the output file format for created Count objects is CountOutputFormat.GDAT.
     | Must not be empty when observables_output_format is explicitly set to CountOutputFormat.GDAT.

   * | default_release_region: Region = None
     | Used as region for releases for seed species that have no compartments specified.

   * | parameter_overrides: Dict[str, float] = None
     | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
     | its value is ignored and instead value parameter_overrides[k] is used.

   * | observables_output_format: CountOutputFormat = CountOutputFormat.AUTOMATIC_FROM_EXTENSION
     | Selection of output format. Default setting uses automatic detection
     | based on contents of the 'observables_path_or_file' attribute.


  | Loads sections\: molecule types, reaction rules, seed species, and observables from a BNGL file
  | and creates objects in the current model according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.
  | If an item in the seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region.

  | Example: `1400_rel_site_for_each_it/model.py <https://github.com/mcellteam/mcell_tests/blob/master/pymcell4/1400_rel_site_for_each_it/model.py>`_ 


* | **export_to_bngl**

   * | file_name: str
     | Output file name.

   * | simulation_method: BNGSimulationMethod = BNGSimulationMethod.ODE
     | Selection of the BioNetGen simulation method. 
     | Selects BioNetGen action to run with the selected simulation method.
     | For BNGSimulationMethod.NF the export is limited to a single volume and
     | a single surface and the enerated rates use volume and surface area so that 
     | simulation with NFSim produces corect results.


  | Exports all defined species, reaction rules and applicable observables
  | as a BNGL file that can be then loaded by MCell4 or BioNetGen. 
  | The resulting file should be validated that it produces expected results. 
  | Many MCell features cannot be exported into BNGL and when such a feature is 
  | encountered the export fails with a RuntimeError exception.
  | However, the export code tries to export as much as possible and one can catch
  | the RuntimeError exception and use the possibly incomplete BNGL file anyway.


* | **save_checkpoint**

   * | custom_dir: str = None
     | Sets custom directory where the checkpoint will be stored. 
     | The default is 'checkpoints/seed_<SEED>/it_<ITERATION>'.


  | Saves current model state as checkpoint. 
  | The default directory structure is checkpoints/seed_<SEED>/it_<ITERATION>,
  | it can be changed by setting 'custom_dir'.
  | If used during an iteration such as in a callback, an event is scheduled for the  
  | beginning of the next iteration. This scheduled event saves the checkpoint.

  | Example: `2700_save_checkpoint_rxn_in_box/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/2700_save_checkpoint_rxn_in_box/model.py>`_ 


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

  | Example: `2760_schedule_checkpoint_async_w_timer/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/nutmeg4_pymcell4/2760_schedule_checkpoint_async_w_timer/model.py>`_ 


* | **add_species**

   * | s: Species

  | Add a reference to a Species object to the species list.


* | **find_species**

   * | name: str
   * | return type: Species


  | Find a Species object using name in the species list. 
  | Returns None if no such species is found.


* | **add_reaction_rule**

   * | r: ReactionRule

  | Add a reference to a ReactionRule object to the reaction_rules list.


* | **find_reaction_rule**

   * | name: str
   * | return type: ReactionRule


  | Find a ReactionRule object using name in the reaction_rules list. 
  | Returns None if no such reaction rule is found.


* | **add_surface_class**

   * | sc: SurfaceClass

  | Add a reference to a SurfaceClass object to the surface_classes list.


* | **find_surface_class**

   * | name: str
   * | return type: SurfaceClass


  | Find a SurfaceClass object using name in the surface_classes list. 
  | Returns None if no such surface class is found.


* | **add_elementary_molecule_type**

   * | mt: ElementaryMoleculeType

  | Add a reference to an ElementaryMoleculeType object to the elementary_molecule_types list.


* | **find_elementary_molecule_type**

   * | name: str
   * | return type: ElementaryMoleculeType


  | Find an ElementaryMoleculeType object using name in the elementary_molecule_types list. 
  | Returns None if no such elementary molecule type is found.


* | **load_bngl_molecule_types_and_reaction_rules**

   * | file_name: str
     | Path to the BNGL file to be loaded.

   * | parameter_overrides: Dict[str, float] = None
     | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
     | its value is ignored and instead value parameter_overrides[k] is used.


  | Parses a BNGL file, only reads molecule types and reaction rules sections, 
  | i.e. ignores observables and seed species. 
  | Parameter values are evaluated and the result value is directly used.  
  | Compartments names are stored in rxn rules as strings because compartments belong 
  | to geometry objects and the subsystem is independent on specific geometry.
  | However, the compartments and their objects must be defined before initialization.

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 


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


* | **load_bngl_compartments_and_seed_species**

   * | file_name: str
     | Path to the BNGL file.

   * | default_release_region: Region = None
     | Used as region for releases for seed species that have no compartments specified.

   * | parameter_overrides: Dict[str, float] = None
     | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
     | its value is ignored and instead value parameter_overrides[k] is used.


  | First loads section compartments and for each 3D compartment that does not 
  | already exist as a geometry object in this Instantiation object, creates a 
  | box with compartment's volume and also sets its 2D (membrane) compartment name.
  | When multiple identical geometry objects are added to the final Model object, 
  | only one copy is left so one can merge multiple Instantiation objects that created 
  | compartments assuming that their volume is the same.        
  | Then loads section seed species from a BNGL file and creates release sites according to it.
  | All elementary molecule types used in the seed species section must be already defined in subsystem.
  | If an item in the BNGL seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are then released into or onto the 
  | default_region.

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 


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

   * | observables_path_or_file: str = ''
     | Directory prefix or file name where observable values will be stored.
     | If a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty 
     | string is used, each observable gets its own file and the output file format for created Count 
     | objects is CountOutputFormat.DAT.
     | If a file has a .gdat extension such as 
     | './react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat', all observable are stored in this 
     | file and the output file format for created Count objects is CountOutputFormat.GDAT.
     | Must not be empty when observables_output_format is explicitly set to CountOutputFormat.GDAT.

   * | parameter_overrides: Dict[str, float] = None
     | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
     | its value is ignored and instead value parameter_overrides[k] is used.

   * | observables_output_format: CountOutputFormat = CountOutputFormat.AUTOMATIC_FROM_EXTENSION
     | Selection of output format. Default setting uses automatic detection
     | based on contents of the 'observables_path_or_file' attribute.


  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 


* | **get_molecule_ids**

   * | pattern: Complex = None
     | BNGL pattern to select molecules based on their species, might use compartments.

   * | return type: List[int]


  | Returns a list of ids of molecules.
  | If the arguments pattern is not set, the list of all molecule ids is returned.  
  | If the argument pattern is set, the list of all molecule ids whose species match 
  | the pattern is returned.

  | Example: `1910_get_molecule_ids_w_pattern/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1910_get_molecule_ids_w_pattern/model.py>`_ 


* | **get_molecule**

   * | id: int
     | Unique id of the molecule to be retrieved.

   * | return type: Molecule


  | Returns a information on a molecule from the simulated environment, 
  | None if the molecule does not exist.

  | Example: `1900_molecule_introspection/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1900_molecule_introspection/model.py>`_ 


* | **get_species_name**

   * | species_id: int
     | Id of the species.

   * | return type: str


  | Returns a string representing canonical species name in the BNGL format.

  | Example: `1850_run_unimol_rxn_in_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1850_run_unimol_rxn_in_callback/model.py>`_ 


* | **get_vertex**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | return type: Vec3


  | Returns coordinates of a vertex.

  | Example: `1340_get_vertex/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1340_get_vertex/model.py>`_ 


* | **get_wall**

   * | object: GeometryObject
     | Geometry object whose wall to retrieve.

   * | wall_index: int
     | This is the index of the wall in the geometry object's walls (wall_list).

   * | return type: Wall


  | Returns information about a wall belonging to a given object.

  | Example: `1330_get_wall/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1330_get_wall/model.py>`_ 


* | **get_vertex_unit_normal**

   * | object: GeometryObject
     | Geometry object whose vertex to retrieve.

   * | vertex_index: int
     | This is the index of the vertex in the geometry object's vertex_list.

   * | return type: Vec3


  | Returns sum of all wall normals that use this vertex converted to a unit vector of 
  | length 1 um (micrometer).
  | This represents the unit vector pointing outwards from the vertex.

  | Example: `1320_get_vertex_unit_normal/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1320_get_vertex_unit_normal/model.py>`_ 


* | **get_wall_unit_normal**

   * | object: GeometryObject
     | Geometry object whose wall's normal to retrieve.

   * | wall_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | return type: Vec3


  | Returns wall normal converted to a unit vector of length 1um.

  | Example: `1310_get_wall_unit_normal/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1310_get_wall_unit_normal/model.py>`_ 


* | **get_wall_color**

   * | object: GeometryObject
     | Geometry object whose wall's color to retrieve.

   * | wall_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | return type: Color


  | Returns color of a wall.


* | **set_wall_color**

   * | object: GeometryObject
     | Geometry object whose wall's color to retrieve.

   * | wall_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | color: Color
     | Color to be set.


  | Sets color of a wall.



