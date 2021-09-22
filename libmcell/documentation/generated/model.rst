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
.. _Model__config:

config: Config
--------------

  | Simulation configuration.
  | - default argument value in constructor: Config()

.. _Model__warnings:

warnings: Warnings
------------------

  | Configuration on how to report warnings.
  | - default argument value in constructor: Warnings()

.. _Model__notifications:

notifications: Notifications
----------------------------

  | Configuration on how to report certain notifications.
  | - default argument value in constructor: Notifications()

.. _Model__species:

species: List[Species]
----------------------

  | List of species to be included in the model for initialization.
  | Used usually only for simple species (species that are defined using a
  | single molecule type without components such as 'A').
  | Other species may be created inside simulation
  | - default argument value in constructor: None

.. _Model__reaction_rules:

reaction_rules: List[ReactionRule]
----------------------------------

  | - default argument value in constructor: None

.. _Model__surface_classes:

surface_classes: List[SurfaceClass]
-----------------------------------

  | - default argument value in constructor: None

.. _Model__elementary_molecule_types:

elementary_molecule_types: List[ElementaryMoleculeType]
-------------------------------------------------------

  | Contains list of elementary molecule types with their diffusion constants and other information. 
  | Populated when a BNGL file is loaded and also on initialization from Species objects present in 
  | the species list.
  | - default argument value in constructor: None

.. _Model__release_sites:

release_sites: List[ReleaseSite]
--------------------------------

  | List of release sites to be included in the model.
  | - default argument value in constructor: None

.. _Model__geometry_objects:

geometry_objects: List[GeometryObject]
--------------------------------------

  | List of geometry objects to be included in the model.
  | - default argument value in constructor: None

.. _Model__checkpointed_molecules:

checkpointed_molecules: List[BaseChkptMol]
------------------------------------------

  | Used when resuming simulation from a checkpoint.
  | - default argument value in constructor: None

.. _Model__viz_outputs:

viz_outputs: List[VizOutput]
----------------------------

  | List of visualization outputs to be included in the model.
  | There is usually just one VizOutput object.
  | - default argument value in constructor: None

.. _Model__counts:

counts: List[Count]
-------------------

  | List of counts to be included in the model.
  | - default argument value in constructor: None


Methods:
*********
.. _Model__initialize:

initialize (print_copyright: bool=True)
---------------------------------------


  | Initializes model, initialization blocks most of changes to 
  | contained components.

* | print_copyright: bool = True
  | Prints information about MCell.


.. _Model__run_iterations:

run_iterations (iterations: float) -> int
-----------------------------------------


  | Runs specified number of iterations. Returns the number of iterations
  | executed (it might be less than the requested number of iterations when 
  | a checkpoint was scheduled).

* | iterations: float
  | Number of iterations to run. Value is truncated to an integer.


.. _Model__end_simulation:

end_simulation (print_final_report: bool=True)
----------------------------------------------


  | Generates the last visualization and reaction output (if they are included 
  | in the model), then flushes all buffers and optionally prints simulation report. 
  | Buffers are also flushed when the Model object is destroyed such as when Ctrl-C
  | is pressed during simulation.

* | print_final_report: bool = True
  | Print information on simulation time and counts of selected events.


.. _Model__add_subsystem:

add_subsystem (subsystem: Subsystem)
------------------------------------


  | Adds all components of a Subsystem object to the model.

* | subsystem: Subsystem

.. _Model__add_instantiation:

add_instantiation (instantiation: Instantiation)
------------------------------------------------


  | Adds all components of an Instantiation object to the model.

* | instantiation: Instantiation

.. _Model__add_observables:

add_observables (observables: Observables)
------------------------------------------


  | Adds all counts and viz outputs of an Observables object to the model.

* | observables: Observables

.. _Model__dump_internal_state:

dump_internal_state (with_geometry: bool=False)
-----------------------------------------------


  | Prints out the simulation engine's internal state, mainly for debugging.

* | with_geometry: bool = False
  | Include geometry in the dump.


.. _Model__export_data_model:

export_data_model (file: str=None)
----------------------------------


  | Exports the current state of the model into a data model JSON format.
  | Does not export state of molecules.
  | Must be called after model initialization.
  | Always exports the current state, i.e. with the current geometry and reaction rates. 
  | Events (ReleaseSites and VizOutputs) with scheduled time other than zero are not exported correctly yet.

* | file: str = None
  | If file is not set, then uses the first VizOutput to determine the target directory 
  | and creates name using the current iteration. Fails if argument file is not set and 
  | there is no VizOutput in the model.


.. _Model__export_viz_data_model:

export_viz_data_model (file: str=None)
--------------------------------------


  | Same as export_data_model, only the created data model will contain only information required for visualization
  | in CellBlender. This makes the loading of the model by CellBlender faster and also allows to avoid potential
  | compatibility issues.
  | Must be called after model initialization.

* | file: str = None
  | Optional path to the output data model file.

  | Example: `1520_sphere_collision/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1520_sphere_collision/model.py>`_ 


.. _Model__export_geometry:

export_geometry (output_files_prefix: str=None)
-----------------------------------------------


  | Exports model geometry as Wavefront OBJ format. 
  | Must be called after model initialization.
  | Does not export material colors (yet).

* | output_files_prefix: str = None
  | Optional prefix for .obj and .mtl files that will be created on export. 
  | If output_files_prefix is not set, then uses the first VizOutput to determine the target directory 
  | and creates names using the current iteration. Fails if argument output_files_prefix is not set and 
  | there is no VizOutput in the model.


.. _Model__release_molecules:

release_molecules (release_site: ReleaseSite)
---------------------------------------------


  | Performs immediate release of molecules based on the definition of the release site argument.
  | The ReleaseSite.release_time must not be in the past and must be within the current iteration 
  | meaning that the time must be greater or equal iteration \* time_step and less than (iteration + 1) \* time_step.
  | The ReleaseEvent must not use a release_pattern because this is an immediate release and it is not 
  | scheduled into the global scheduler.

* | release_site: ReleaseSite
  | Example: `2300_immediate_release/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2300_immediate_release/model.py>`_ 


.. _Model__run_reaction:

run_reaction (reaction_rule: ReactionRule, reactant_ids: List[int], time: float) -> List[int]
---------------------------------------------------------------------------------------------


  | Run a single reaction on reactants. Callbacks will be called if they are registered for the given reaction.
  | Returns a list of product IDs.
  | Note\: only unimolecular reactions are currently supported.

* | reaction_rule: ReactionRule
  | Reaction rule to run.

* | reactant_ids: List[int]
  | The number of reactants for a unimolecular reaction must be 1 and for a bimolecular reaction must be 2.
  | Reactants for a bimolecular reaction do not have to be listed in the same order as in the reaction rule definition.

* | time: float
  | Precise time in seconds when this reaction occurs. Important to know for how long the products
  | will be diffused when they are created in a middle of a time step.

  | Example: `1850_run_unimol_rxn_in_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1850_run_unimol_rxn_in_callback/model.py>`_ 


.. _Model__add_vertex_move:

add_vertex_move (object: GeometryObject, vertex_index: int, displacement: List[float])
--------------------------------------------------------------------------------------


  | Appends information about a displacement for given object's vertex into an internal list of vertex moves. 
  | To do the actual geometry change, call Model.apply_vertex_moves.
  | The reason why we first need to collect all changes and then apply them all at the same time is for performance
  | reasons.

* | object: GeometryObject
  | Object whose vertex will be changed.

* | vertex_index: int
  | Index of vertex in object's vertex list that will be changed.

* | displacement: List[float]
  | Change of vertex coordinates [x, y, z] (in um) that will be added to the current 
  | coordinates of the vertex.

  | Example: `1510_tetrahedron_box_collision_moving_3_verts/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1510_tetrahedron_box_collision_moving_3_verts/model.py>`_ 


.. _Model__apply_vertex_moves:

apply_vertex_moves (collect_wall_wall_hits: bool=False, randomize_order: bool=True) -> List[WallWallHitInfo]
------------------------------------------------------------------------------------------------------------


  | Applies all the vertex moves specified with Model.add_vertex_move call.
  | 
  | All affected vertices are first divided based on to which geometery object they belong. 
  | Then each object is manipulated one by one. 
  | 
  | During vertex moves, collisions are checked\:
  | a) When a moved vertex hits a wall of another object, it is stopped at the wall.
  | b) When a second object's vertex would end up inside the moved object, the vertex move 
  | that would cause it is canceled (its displacement set to 0) because finding the maximum 
  | distance we can move is too computationally expensive. To minimize the impact of this 
  | cancellation, the vertices should be moved only by a small distance.
  | 
  | Applying vertex moves also takes paired molecules into account\: 
  | When moves are applied to an object, all moved molecules that are paired are collected.
  | For each of the paired molecules, we collect displacements for each 
  | of the vertices of the 'primary' wall where this molecule is located (that were provided by the user 
  | through add_vertex_move, and were possibly truncated due to collisions).
  | Then we find the second wall where the second molecule of the pair is located.
  | For each of the vertices of all 'secondary' walls, we collect a list of displacements
  | that move the vertices of 'primary' walls. 
  | Then, an average displacement is computed for each vertex, and these average displacements
  | are used to move the 'secondary' walls.
  | When a 'primary' wall collides, its displacement is clamped or canceled. This is true even if 
  | it collides with a 'secondary' wall that would be otherwise moved. So, the displacement of the 
  | 'primary' wall will mostly just pull the 'secondary' wall, not push. Therefore it is needed 
  | that both objects are active and pull each other. 
  | 
  | This process is well commented in MCell code\: 
  | `partition.cpp <https://github.com/mcellteam/mcell/blob/master/src4/partition.cpp>`_ in functions
  | apply_vertex_moves, apply_vertex_moves_per_object, and move_walls_with_paired_molecules. 
  |      
  | When argument collect_wall_wall_hits is True, a list of wall pairs that collided is returned,
  | when collect_wall_wall_hits is False, an empty list is returned.

* | collect_wall_wall_hits: bool = False
  | When set to True, a list of wall pairs that collided is returned,
  | otherwise an empty list is returned.

* | randomize_order: bool = True
  | When set to True (default), the ordering of the vertex move list created by add_vertex_move
  | calls is randomized. This allows to avoid any bias in the resulting positions of surface
  | molecules.  
  | However, the individual vertex moves are then sorted by the object to which the vertex belongs
  | and the moves are applied object by object for correctness. Setting this to True also radomizes the 
  | order of objects to which the vertex moves are applied.

  | Examples: `1510_tetrahedron_box_collision_moving_3_verts/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1510_tetrahedron_box_collision_moving_3_verts/model.py>`_ `3200_sphere_collision_against_each_other/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/3200_sphere_collision_against_each_other/model.py>`_ `3150_dyn_vert_intramembrane_rxns_and_paired_mols/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/3150_dyn_vert_intramembrane_rxns_and_paired_mols/model.py>`_ 


.. _Model__pair_molecules:

pair_molecules (id1: int, id2: int)
-----------------------------------


  | Sets that two surface molecules are paired. Paired molecules bind walls together
  | and when one wall is moved, the wall that is bound through a paired molecule is moved as well.
  | Throws exception if the molecule ids are not surface molecules.
  | Throws exception if the molecules are on the same object.  
  | Throws exception if any of the molecules is already paired.
  | May be called only after model initialization.

* | id1: int
* | id2: int
  | Examples: `2900_pair_unpair_molecules/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/nutmeg4_pymcell4/2900_pair_unpair_molecules/model.py>`_ `3160_dyn_vert_paired_mols_box_box/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/3160_dyn_vert_paired_mols_box_box/model.py>`_ 


.. _Model__unpair_molecules:

unpair_molecules (id1: int, id2: int)
-------------------------------------


  | Sets that two surface molecules are not paired. 
  | Throws exception if the molecule ids are not surface molecules. 
  | Throws exception if the molecules are not paired together.
  | May be called only after model initialization.

* | id1: int
* | id2: int
  | Example: `2900_pair_unpair_molecules/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/nutmeg4_pymcell4/2900_pair_unpair_molecules/model.py>`_ 


.. _Model__get_paired_molecule:

get_paired_molecule (id: int) -> int
------------------------------------


  | Return id of the molecule to which the molecule with 'id' is paired.
  | Returns ID_INVALID (-1) when the molecule is not paired.
  | May be called only after model initialization.

* | id: int
  | Example: `2900_pair_unpair_molecules/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/nutmeg4_pymcell4/2900_pair_unpair_molecules/model.py>`_ 


.. _Model__get_paired_molecules:

get_paired_molecules () -> Dict[uint32, uint32]
-----------------------------------------------


  | Returns a dictionary that contains all molecules that are paired.
  | Molecule ids are keys and the value associated with the key is the second paired molecule.
  | The returned dictionary is a copy and any changes made to it are ignored by MCell.
  | Note\: The reason why uint32 is used as the base type for the dictionary but type int is used
  | everywhere else for molecule ids is only for performance reasons.

  | Example: `3170_get_paired_molecules/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/3170_get_paired_molecules/model.py>`_ 


.. _Model__register_mol_wall_hit_callback:

register_mol_wall_hit_callback (function: Callable, # std::function<void(std::shared_ptr<MolWallHitInfo>, py::object)>, context: Any, # py::object, object: GeometryObject=None, species: Species=None)
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  | Register a callback for event when a molecule hits a wall. 
  | May be called only after model initialization because it internally uses geometry object
  | and species ids that are set during the initialization.

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

  | Example: `1300_wall_hit_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1300_wall_hit_callback/model.py>`_ 


.. _Model__register_reaction_callback:

register_reaction_callback (function: Callable, # std::function<void(std::shared_ptr<ReactionInfo>, py::object)>, context: Any, # py::object, reaction_rule: ReactionRule)
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  | Defines a function to be called when a reaction was processed.
  | It is allowed to do state modifications except for removing reacting molecules, 
  | they will be removed automatically after return from this callback. 
  | Unlimited number of reaction callbacks is allowed. 
  | May be called only after model initialization because it internally uses 
  | reaction rule ids that are set during the initialization.

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

  | Example: `1800_vol_rxn_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1800_vol_rxn_callback/model.py>`_ 


.. _Model__load_bngl:

load_bngl (file_name: str, observables_path_or_file: str=None, default_release_region: Region=None, parameter_overrides: Dict[str, float]=None, observables_output_format: CountOutputFormat=CountOutputFormat.AUTOMATIC_FROM_EXTENSION)
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  | Loads sections\: molecule types, reaction rules, seed species, and observables from a BNGL file
  | and creates objects in the current model according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.
  | If an item in the seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region.

* | file_name: str
  | Path to the BNGL file to be loaded.

* | observables_path_or_file: str = None
  | Directory prefix or file name where observable values will be stored.
  | If a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty 
  | string/unset is used, each observable gets its own file and the output file format for created Count 
  | objects is CountOutputFormat.DAT.
  | When not set, this path is used: './react_data/seed_' + str(model.config.seed).zfill(5) + '/'.
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

  | Example: `1400_rel_site_for_each_it/model.py <https://github.com/mcellteam/mcell_tests/blob/master/pymcell4/1400_rel_site_for_each_it/model.py>`_ 


.. _Model__export_to_bngl:

export_to_bngl (file_name: str, simulation_method: BNGSimulationMethod=BNGSimulationMethod.ODE)
-----------------------------------------------------------------------------------------------


  | Exports all defined species, reaction rules and applicable observables
  | as a BNGL file that can be then loaded by MCell4 or BioNetGen. 
  | The resulting file should be validated that it produces expected results. 
  | Many MCell features cannot be exported into BNGL and when such a feature is 
  | encountered the export fails with a RuntimeError exception.
  | However, the export code tries to export as much as possible and one can catch
  | the RuntimeError exception and use the possibly incomplete BNGL file anyway.

* | file_name: str
  | Output file name.

* | simulation_method: BNGSimulationMethod = BNGSimulationMethod.ODE
  | Selection of the BioNetGen simulation method. 
  | Selects BioNetGen action to run with the selected simulation method.
  | For BNGSimulationMethod.NF the export is limited to a single volume and
  | a single surface and the enerated rates use volume and surface area so that 
  | simulation with NFSim produces corect results.


.. _Model__save_checkpoint:

save_checkpoint (custom_dir: str=None)
--------------------------------------


  | Saves current model state as checkpoint. 
  | The default directory structure is checkpoints/seed_<SEED>/it_<ITERATION>,
  | it can be changed by setting 'custom_dir'.
  | If used during an iteration such as in a callback, an event is scheduled for the  
  | beginning of the next iteration. This scheduled event saves the checkpoint.

* | custom_dir: str = None
  | Sets custom directory where the checkpoint will be stored. 
  | The default is 'checkpoints/seed_<SEED>/it_<ITERATION>'.

  | Example: `2700_save_checkpoint_rxn_in_box/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/2700_save_checkpoint_rxn_in_box/model.py>`_ 


.. _Model__schedule_checkpoint:

schedule_checkpoint (iteration: int=0, continue_simulation: bool=False, custom_dir: str=None)
---------------------------------------------------------------------------------------------


  | Schedules checkpoint save event that will occur when an iteration is started.  
  | This means that it will be executed right before any other events scheduled for 
  | the given iteration are executed.
  | Can be called asynchronously at any time after initialization.

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

  | Example: `2760_schedule_checkpoint_async_w_timer/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/nutmeg4_pymcell4/2760_schedule_checkpoint_async_w_timer/model.py>`_ 


.. _Model__add_species:

add_species (s: Species)
------------------------


  | Add a reference to a Species object to the species list.

* | s: Species

.. _Model__find_species:

find_species (name: str) -> Species
-----------------------------------


  | Find a Species object using name in the species list. 
  | Returns None if no such species is found.

* | name: str

.. _Model__add_reaction_rule:

add_reaction_rule (r: ReactionRule)
-----------------------------------


  | Add a reference to a ReactionRule object to the reaction_rules list.

* | r: ReactionRule

.. _Model__find_reaction_rule:

find_reaction_rule (name: str) -> ReactionRule
----------------------------------------------


  | Find a ReactionRule object using name in the reaction_rules list. 
  | Returns None if no such reaction rule is found.

* | name: str

.. _Model__add_surface_class:

add_surface_class (sc: SurfaceClass)
------------------------------------


  | Add a reference to a SurfaceClass object to the surface_classes list.

* | sc: SurfaceClass

.. _Model__find_surface_class:

find_surface_class (name: str) -> SurfaceClass
----------------------------------------------


  | Find a SurfaceClass object using name in the surface_classes list. 
  | Returns None if no such surface class is found.

* | name: str

.. _Model__add_elementary_molecule_type:

add_elementary_molecule_type (mt: ElementaryMoleculeType)
---------------------------------------------------------


  | Add a reference to an ElementaryMoleculeType object to the elementary_molecule_types list.

* | mt: ElementaryMoleculeType

.. _Model__find_elementary_molecule_type:

find_elementary_molecule_type (name: str) -> ElementaryMoleculeType
-------------------------------------------------------------------


  | Find an ElementaryMoleculeType object using name in the elementary_molecule_types list. 
  | Returns None if no such elementary molecule type is found.

* | name: str

.. _Model__load_bngl_molecule_types_and_reaction_rules:

load_bngl_molecule_types_and_reaction_rules (file_name: str, parameter_overrides: Dict[str, float]=None)
--------------------------------------------------------------------------------------------------------


  | Parses a BNGL file, only reads molecule types and reaction rules sections, 
  | i.e. ignores observables and seed species. 
  | Parameter values are evaluated and the result value is directly used.  
  | Compartments names are stored in rxn rules as strings because compartments belong 
  | to geometry objects and the subsystem is independent on specific geometry.
  | However, the compartments and their objects must be defined before initialization.

* | file_name: str
  | Path to the BNGL file to be loaded.

* | parameter_overrides: Dict[str, float] = None
  | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
  | its value is ignored and instead value parameter_overrides[k] is used.

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 


.. _Model__add_release_site:

add_release_site (s: ReleaseSite)
---------------------------------


  | Adds a reference to the release site s to the list of release sites.

* | s: ReleaseSite

.. _Model__find_release_site:

find_release_site (name: str) -> ReleaseSite
--------------------------------------------


  | Finds a release site by its name, returns None if no such release site is present.

* | name: str

.. _Model__add_geometry_object:

add_geometry_object (o: GeometryObject)
---------------------------------------


  | Adds a reference to the geometry object o to the list of geometry objects.

* | o: GeometryObject

.. _Model__find_geometry_object:

find_geometry_object (name: str) -> GeometryObject
--------------------------------------------------


  | Finds a geometry object by its name, returns None if no such geometry object is present.

* | name: str

.. _Model__find_volume_compartment_object:

find_volume_compartment_object (name: str) -> GeometryObject
------------------------------------------------------------


  | Finds a geometry object by its name, the geometry object must be a BNGL compartment.
  | Returns None if no such geometry object is present.

* | name: str

.. _Model__find_surface_compartment_object:

find_surface_compartment_object (name: str) -> GeometryObject
-------------------------------------------------------------


  | Finds a geometry object that is a BNGL compartment and its surface name is name.
  | Returns None if no such geometry object is present.

* | name: str

.. _Model__load_bngl_compartments_and_seed_species:

load_bngl_compartments_and_seed_species (file_name: str, default_release_region: Region=None, parameter_overrides: Dict[str, float]=None)
-----------------------------------------------------------------------------------------------------------------------------------------


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

* | file_name: str
  | Path to the BNGL file.

* | default_release_region: Region = None
  | Used as region for releases for seed species that have no compartments specified.

* | parameter_overrides: Dict[str, float] = None
  | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
  | its value is ignored and instead value parameter_overrides[k] is used.

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 


.. _Model__add_viz_output:

add_viz_output (viz_output: VizOutput)
--------------------------------------


  | Adds a reference to the viz_output object to the list of visualization output specifications.

* | viz_output: VizOutput

.. _Model__add_count:

add_count (count: Count)
------------------------


  | Adds a reference to the count object to the list of count specifications.

* | count: Count

.. _Model__find_count:

find_count (name: str) -> Count
-------------------------------


  | Finds a count object by its name, returns None if no such count is present.

* | name: str

.. _Model__load_bngl_observables:

load_bngl_observables (file_name: str, observables_path_or_file: str=None, parameter_overrides: Dict[str, float]=None, observables_output_format: CountOutputFormat=CountOutputFormat.AUTOMATIC_FROM_EXTENSION)
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.

* | file_name: str
  | Path to the BNGL file.

* | observables_path_or_file: str = None
  | Directory prefix or file name where observable values will be stored.
  | If a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty 
  | string/unset is used, each observable gets its own file and the output file format for created Count 
  | objects is CountOutputFormat.DAT.
  | When not set, this path is used: './react_data/seed_' + str(model.config.seed).zfill(5) + '/'.
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

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 


.. _Model__get_molecule_ids:

get_molecule_ids (pattern: Complex=None) -> List[int]
-----------------------------------------------------


  | Returns a list of ids of molecules.
  | If the arguments pattern is not set, the list of all molecule ids is returned.  
  | If the argument pattern is set, the list of all molecule ids whose species match 
  | the pattern is returned.

* | pattern: Complex = None
  | BNGL pattern to select molecules based on their species, might use compartments.

  | Example: `1910_get_molecule_ids_w_pattern/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1910_get_molecule_ids_w_pattern/model.py>`_ 


.. _Model__get_molecule:

get_molecule (id: int) -> Molecule
----------------------------------


  | Returns a information on a molecule from the simulated environment, 
  | None if the molecule does not exist.

* | id: int
  | Unique id of the molecule to be retrieved.

  | Example: `1900_molecule_introspection/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1900_molecule_introspection/model.py>`_ 


.. _Model__get_species_name:

get_species_name (species_id: int) -> str
-----------------------------------------


  | Returns a string representing canonical species name in the BNGL format.

* | species_id: int
  | Id of the species.

  | Example: `1850_run_unimol_rxn_in_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1850_run_unimol_rxn_in_callback/model.py>`_ 


.. _Model__get_vertex:

get_vertex (object: GeometryObject, vertex_index: int) -> List[float]
---------------------------------------------------------------------


  | Returns coordinates of a vertex.

* | object: GeometryObject
* | vertex_index: int
  | This is the index of the vertex in the geometry object's walls (wall_list).

  | Example: `1340_get_vertex/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1340_get_vertex/model.py>`_ 


.. _Model__get_wall:

get_wall (object: GeometryObject, wall_index: int) -> Wall
----------------------------------------------------------


  | Returns information about a wall belonging to a given object.

* | object: GeometryObject
  | Geometry object whose wall to retrieve.

* | wall_index: int
  | This is the index of the wall in the geometry object's walls (wall_list).

  | Example: `1330_get_wall/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1330_get_wall/model.py>`_ 


.. _Model__get_vertex_unit_normal:

get_vertex_unit_normal (object: GeometryObject, vertex_index: int) -> List[float]
---------------------------------------------------------------------------------


  | Returns sum of all wall normals that use this vertex converted to a unit vector of 
  | length 1 um (micrometer).
  | This represents the unit vector pointing outwards from the vertex.

* | object: GeometryObject
  | Geometry object whose vertex to retrieve.

* | vertex_index: int
  | This is the index of the vertex in the geometry object's vertex_list.

  | Example: `1320_get_vertex_unit_normal/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1320_get_vertex_unit_normal/model.py>`_ 


.. _Model__get_wall_unit_normal:

get_wall_unit_normal (object: GeometryObject, wall_index: int) -> List[float]
-----------------------------------------------------------------------------


  | Returns wall normal converted to a unit vector of length 1um.

* | object: GeometryObject
  | Geometry object whose wall's normal to retrieve.

* | wall_index: int
  | This is the index of the vertex in the geometry object's walls (wall_list).

  | Example: `1310_get_wall_unit_normal/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1310_get_wall_unit_normal/model.py>`_ 


.. _Model__get_wall_color:

get_wall_color (object: GeometryObject, wall_index: int) -> Color
-----------------------------------------------------------------


  | Returns color of a wall.

* | object: GeometryObject
  | Geometry object whose wall's color to retrieve.

* | wall_index: int
  | This is the index of the vertex in the geometry object's walls (wall_list).


.. _Model__set_wall_color:

set_wall_color (object: GeometryObject, wall_index: int, color: Color)
----------------------------------------------------------------------


  | Sets color of a wall.

* | object: GeometryObject
  | Geometry object whose wall's color to retrieve.

* | wall_index: int
  | This is the index of the vertex in the geometry object's walls (wall_list).

* | color: Color
  | Color to be set.



