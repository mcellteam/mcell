.. _api-simulation_setup:

****************
Simulation setup
****************
Config
======

Class holds simulation configuration.

Attributes:
***********
.. _Config__seed:

seed: int
---------

  | Random generator seed value.
  | - default argument value in constructor: 1

.. _Config__time_step:

time_step: float
----------------

  | Set the simulation time step to time_step seconds. 1e-6 (1us) is a common value. 
  | One can set the time steps taken by individual molecules, but this 
  | time step is still used as a default.
  | - default argument value in constructor: 1e-6

.. _Config__use_bng_units:

use_bng_units: bool
-------------------

  | When False (default), MCell uses traditional MCell units for bimolecular reaction rates are:
  |  \* [M^-1\*s^-1] for bimolecular reactions between either two volume molecules, a volume molecule 
  |                and a surface (molecule), 
  |  \* [um^2\*N^-1\*s^-1] bimolecular reactions between two surface molecules on the same surface.
  | When True, BioNetGen units for bimolecular reaction rates are:
  |  \* [um^3\*N^-1\*s^-1] for any bimolecular reactions. Surface-surface reaction rate conversion assumes 10nm membrane thickness
  | BioNetGen units are compatible with BioNetGen's ODE, SSA, and PLA solvers given that seed species 
  | is copy number (N), these units are not compatible with NFSim. 
  | No other units are affected by this setting.
  | - default argument value in constructor: False

.. _Config__surface_grid_density:

surface_grid_density: float
---------------------------

  | Tile all surfaces so that they can hold molecules at N different positions per square micron.
  | - default argument value in constructor: 10000

.. _Config__interaction_radius:

interaction_radius: float
-------------------------

  | Diffusing volume molecules will interact with each other when
  | they get within N microns of each other. The default is
  | 1/sqrt(PI \* Sigma_s) where Sigma_s is the surface grid density 
  | (default or user-specified).
  | - default argument value in constructor: None

.. _Config__intermembrane_interaction_radius:

intermembrane_interaction_radius: float
---------------------------------------

  | Diffusing surface molecules will interact with surface molecules on other
  | walls when they get within N microns of each other. The default is
  | 1/sqrt(PI \* Sigma_s) where Sigma_s is the surface grid density 
  | (default or user-specified). 
  | When unset, the default value is computed as: 
  | 1.0 / sqrt_f(MY_PI \* surface_grid_density).
  | - default argument value in constructor: None

  | Example: `3000_intermembrane_rxns/customization.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/3000_intermembrane_rxns/customization.py>`_ 


.. _Config__vacancy_search_distance:

vacancy_search_distance: float
------------------------------

  | Rather internal, there is usually no need to change this value.
  | Used in dynamic geometry (see Model.apply_vertex_moves()). 
  | When a wall moves or its dimensions change, this is the maximum search distance 
  | use when looking onto which tiles place the molecules on this wall. 
  | If no empty tile is found within this distance, simulation fails.
  | - default argument value in constructor: 10

  | Example: `1200_dyn_vert_tetrahedron_vol_mol_multiple_changes/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1200_dyn_vert_tetrahedron_vol_mol_multiple_changes/model.py>`_ 


.. _Config__center_molecules_on_grid:

center_molecules_on_grid: bool
------------------------------

  | If set to True, then all molecules on a surface will be
  | located exactly at the center of their grid element. If False, the
  | molecules will be randomly located when placed, and reactions
  | will take place at the location of the target (or the site of impact
  | in the case of 3D molecule/surface reactions).
  | - default argument value in constructor: False

  | Example: `1210_dyn_vert_tetrahedron_surf_mol_multiple_changes/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1210_dyn_vert_tetrahedron_surf_mol_multiple_changes/model.py>`_ 


.. _Config__partition_dimension:

partition_dimension: float
--------------------------

  | All the simulated 3d space is placed in a partition. The partition is a cube and 
  | this partition_dimension specifies the length of its edge in um.
  | - default argument value in constructor: 10

  | Example: `1100_point_release/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1100_point_release/model.py>`_ 


.. _Config__initial_partition_origin:

initial_partition_origin: List[float]
-------------------------------------

  | Optional placement of the initial partition in um, specifies the left, lower front 
  | point. If not set, value -partition_dimension/2 is used for each of the dimensions 
  | placing the center of the partition to (0, 0, 0).
  | - default argument value in constructor: None

.. _Config__subpartition_dimension:

subpartition_dimension: float
-----------------------------

  | Subpartition are spatial division of 3D space used to accelerate collision checking.
  | In general, partitions should be chosen to avoid having too many surfaces and molecules
  | in one subpartition. 
  | If there are few surfaces and/or molecules in a subvolume, it is advantageous to have the 
  | subvolume as large as possible. Crossing partition boundaries takes a small amount of time, 
  | so it is rarely useful to have partitions more finely spaced than the average diffusion distance 
  | of the faster-moving molecules in the simulation.
  | - default argument value in constructor: 0.5

  | Example: `2000_bngl_a_plus_b_to_c_partitioning/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2000_bngl_a_plus_b_to_c_partitioning/model.py>`_ 


.. _Config__total_iterations:

total_iterations: float
-----------------------

  | Required for checkpointing so that the checkpointed model has information on
  | the intended total number of iterations. 
  | Also used when generating visualization data files and also for other reporting uses. 
  | Value is truncated to an integer.
  | - default argument value in constructor: 1000000

.. _Config__check_overlapped_walls:

check_overlapped_walls: bool
----------------------------

  | Enables check for overlapped walls. Overlapping walls can cause issues during 
  | simulation such as a molecule escaping closed geometry when it hits two walls 
  | that overlap.
  | - default argument value in constructor: True

.. _Config__reaction_class_cleanup_periodicity:

reaction_class_cleanup_periodicity: int
---------------------------------------

  | Reaction class cleanup removes computed reaction classes for inactive species from memory.
  | This provides faster reaction lookup faster but when the same reaction class is 
  | needed again, it must be recomputed.
  | - default argument value in constructor: 500

  | Example: `2701_concentration_based_rxn_rate_cleanup_check/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2701_concentration_based_rxn_rate_cleanup_check/model.py>`_ 


.. _Config__species_cleanup_periodicity:

species_cleanup_periodicity: int
--------------------------------

  | Species cleanup removes inactive species from memory. It removes also all reaction classes 
  | that reference it.
  | This provides faster addition of new species lookup faster but when the species is 
  | needed again, it must be recomputed.
  | - default argument value in constructor: 10000

  | Example: `2701_concentration_based_rxn_rate_cleanup_check/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2701_concentration_based_rxn_rate_cleanup_check/model.py>`_ 


.. _Config__molecules_order_random_shuffle_periodicity:

molecules_order_random_shuffle_periodicity: int
-----------------------------------------------

  | Randomly shuffle the order in which molecules are simulated.
  | This helps to overcome potential biases that may occur when 
  | molecules are ordered e.g. by their species when simulation starts. 
  | The first shuffling occurs at this iteration, i.e. no shuffle is done at iteration 0.
  | Setting this parameter to 0 disables the shuffling.
  | - default argument value in constructor: 10000

.. _Config__sort_molecules:

sort_molecules: bool
--------------------

  | Enables sorting of molecules for diffusion, this may improve cache locality and provide 
  | slightly better performance. 
  | Produces different results for the same seed when enabled because molecules are simulated 
  | in a different order.
  | - default argument value in constructor: False

.. _Config__memory_limit_gb:

memory_limit_gb: int
--------------------

  | Sets memory limit in GB for simulation run. 
  | When this limit is hit, all buffers are flushed and simulation is terminated with an error.
  | - default argument value in constructor: -1

  | Example: `0200_memory_limit/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/nutmeg4_pymcell4/0200_memory_limit/model.py>`_ 


.. _Config__initial_iteration:

initial_iteration: int
----------------------

  | Initial iteration, used when resuming a checkpoint.
  | - default argument value in constructor: 0

.. _Config__initial_time:

initial_time: float
-------------------

  | Initial time in us, used when resuming a checkpoint.
  | Will be truncated to be a multiple of time step.
  | - default argument value in constructor: 0

.. _Config__initial_rng_state:

initial_rng_state: RngState
---------------------------

  | Used for checkpointing, may contain state of the random number generator to be set 
  | after initialization right before the first event is started. 
  | When not set, the set 'seed' value is used to initialize the random number generator.
  | - default argument value in constructor: None

.. _Config__append_to_count_output_data:

append_to_count_output_data: bool
---------------------------------

  | Used for checkpointing, instead of creating new files for Count observables data, 
  | new values are appended to the existing files. If such files do not exist, new files are
  | created.
  | - default argument value in constructor: False

.. _Config__continue_after_sigalrm:

continue_after_sigalrm: bool
----------------------------

  | MCell registers a SIGALRM signal handler. When SIGALRM signal is received and 
  | continue_after_sigalrm is False, checkpoint is stored and simulation is terminated. 
  | When continue_after_sigalrm is True, checkpoint is stored and simulation continues.
  | SIGALRM is not supported on Windows.
  | - default argument value in constructor: False

  | Example: `2785_schedule_checkpoint_async_w_sigalrm_continue/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/nutmeg4_pymcell4/2785_schedule_checkpoint_async_w_sigalrm_continue/model.py>`_ 


Notifications
=============

Attributes:
***********
.. _Notifications__bng_verbosity_level:

bng_verbosity_level: int
------------------------

  | Sets verbosity level that enables printouts of extra information on BioNetGen 
  | species and rules created and used during simulation.
  | - default argument value in constructor: 0

.. _Notifications__rxn_and_species_report:

rxn_and_species_report: bool
----------------------------

  | When set to True, simulation generates files rxn_report_SEED.txt, and 
  | species_report_SEED.txt that contain details on reaction classes and species 
  | that were created based on reaction rules.
  | - default argument value in constructor: False

.. _Notifications__simulation_stats_every_n_iterations:

simulation_stats_every_n_iterations: int
----------------------------------------

  | When set to a value other than 0, internal simulation stats will be printed.
  | - default argument value in constructor: 0

.. _Notifications__rxn_probability_changed:

rxn_probability_changed: bool
-----------------------------

  | When True, information that a reaction's probability has changed is printed during simulation.
  | - default argument value in constructor: True

.. _Notifications__iteration_report:

iteration_report: bool
----------------------

  | When True, a running report of how many iterations have completed, chosen based 
  | on the total number of iterations, will be printed during simulation.
  | - default argument value in constructor: True

.. _Notifications__wall_overlap_report:

wall_overlap_report: bool
-------------------------

  | When True, information on wall overlaps will be printed.
  | - default argument value in constructor: False

Warnings
========

This class contains warnings settings. For now it contains only one configurable 
warning.

Attributes:
***********
.. _Warnings__high_reaction_probability:

high_reaction_probability: WarningLevel
---------------------------------------

  | Print a warning when a bimolecular reaction probability is over 0.5 but less or equal than 1.
  | Warning when probability is greater than 1 is always printed.
  | Cannot be set to WarningLevel.ERROR.
  | - default argument value in constructor: WarningLevel.IGNORE

  | Example: `0615_bimol_rxn_prob_over_05_less_1_warning_disabled/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/nutmeg4_pymcell4/0615_bimol_rxn_prob_over_05_less_1_warning_disabled/model.py>`_ 


.. _Warnings__molecule_placement_failure:

molecule_placement_failure: WarningLevel
----------------------------------------

  | Print a warning or end with an error when a release of a molecule fails.
  | - default argument value in constructor: WarningLevel.ERROR

