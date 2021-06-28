.. _api-simulation_setup:

****************
Simulation setup
****************
Config
======

Class holds simulation configuration.

Attributes:
***********
* | **seed**: int = 1
  | Random generator seed value.

* | **time_step**: float = 1e-6
  | Set the simulation time step to time_step seconds. 1e-6 (1us) is a common value. 
  | One can set the time steps taken by individual molecules, but this 
  | time step is still used as a default.

* | **surface_grid_density**: float = 10000
  | Tile all surfaces so that they can hold molecules at N different positions per square micron.

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
  | When unset, the default value is computed as: 
  | 1.0 / sqrt_f(MY_PI \* surface_grid_density).

  | Example: `3000_intermembrane_rxns/customization.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/3000_intermembrane_rxns/customization.py>`_ 


* | **vacancy_search_distance**: float = 10
  | Rather internal, there is usually no need to change this value.
  | Used in dynamic geometry (see Model.apply_vertex_moves()). 
  | When a wall moves or its dimensions change, this is the maximum search distance 
  | use when looking onto which tiles place the molecules on this wall. 
  | If no empty tile is found within this distance, simulation fails.

  | Example: `1200_dyn_vert_tetrahedron_vol_mol_multiple_changes/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/1200_dyn_vert_tetrahedron_vol_mol_multiple_changes/model.py>`_ 


* | **center_molecules_on_grid**: bool = False
  | If set to True, then all molecules on a surface will be
  | located exactly at the center of their grid element. If False, the
  | molecules will be randomly located when placed, and reactions
  | will take place at the location of the target (or the site of impact
  | in the case of 3D molecule/surface reactions).

  | Example: `1210_dyn_vert_tetrahedron_surf_mol_multiple_changes/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/1210_dyn_vert_tetrahedron_surf_mol_multiple_changes/model.py>`_ 


* | **partition_dimension**: float = 10
  | All the simulated 3d space is placed in a partition. The partition is a cube and 
  | this partition_dimension specifies the length of its edge in um.

  | Example: `1100_point_release/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/1100_point_release/model.py>`_ 


* | **initial_partition_origin**: List[float] = None
  | Optional placement of the initial partition in um, specifies the left, lower front 
  | point. If not set, value -partition_dimension/2 is used for each of the dimensions 
  | placing the center of the partition to (0, 0, 0).

* | **subpartition_dimension**: float = 0.5
  | Subpartition are spatial division of 3D space used to accelerate collision checking.
  | In general, partitions should be chosen to avoid having too many surfaces and molecules
  | in one subpartition. 
  | If there are few surfaces and/or molecules in a subvolume, it is advantageous to have the 
  | subvolume as large as possible. Crossing partition boundaries takes a small amount of time, 
  | so it is rarely useful to have partitions more finely spaced than the average diffusion distance 
  | of the faster-moving molecules in the simulation.

  | Example: `2000_bngl_a_plus_b_to_c_partitioning/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/2000_bngl_a_plus_b_to_c_partitioning/model.py>`_ 


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

  | Example: `2701_concentration_based_rxn_rate_cleanup_check/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/2701_concentration_based_rxn_rate_cleanup_check/model.py>`_ 


* | **species_cleanup_periodicity**: int = 10000
  | Species cleanup removes inactive species from memory. It removes also all reaction classes 
  | that reference it.
  | This provides faster addition of new species lookup faster but when the species is 
  | needed again, it must be recomputed.

  | Example: `2701_concentration_based_rxn_rate_cleanup_check/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/2701_concentration_based_rxn_rate_cleanup_check/model.py>`_ 


* | **molecules_order_random_shuffle_periodicity**: int = 10000
  | Randomly shuffle the order in which molecules are simulated.
  | This helps to overcome potential biases that may occur when 
  | molecules are ordered e.g. by their species when simulation starts. 
  | The first shuffling occurs at this iteration, i.e. no shuffle is done at iteration 0.
  | Setting this parameter to 0 disables the shuffling.

* | **sort_molecules**: bool = False
  | Enables sorting of molecules for diffusion, this may improve cache locality and provide 
  | slightly better performance. 
  | Produces different results for the same seed when enabled because molecules are simulated 
  | in a different order.

* | **memory_limit_gb**: int = -1
  | Sets memory limit in GB for simulation run. 
  | When this limit is hit, all buffers are flushed and simulation is terminated with an error.

  | Example: `0200_memory_limit/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/nutmeg4_pymcell4/0200_memory_limit/model.py>`_ 


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
  | SIGALRM is not supported on Windows.

  | Example: `2785_schedule_checkpoint_async_w_sigalrm_continue/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/nutmeg4_pymcell4/2785_schedule_checkpoint_async_w_sigalrm_continue/model.py>`_ 


Notifications
=============

Attributes:
***********
* | **bng_verbosity_level**: int = 0
  | Sets verbosity level that enables printouts of extra information on BioNetGen 
  | species and rules created and used during simulation.

* | **rxn_and_species_report**: bool = False
  | When set to True, simulation generates files rxn_report_SEED.txt, and 
  | species_report_SEED.txt that contain details on reaction classes and species 
  | that were created based on reaction rules.

* | **simulation_stats_every_n_iterations**: int = 0
  | When set to a value other than 0, internal simulation stats will be printed.

* | **rxn_probability_changed**: bool = True
  | When True, information that a reaction's probability has changed is printed during simulation.

* | **iteration_report**: bool = True
  | When True, a running report of how many iterations have completed, chosen based 
  | on the total number of iterations, will be printed during simulation.

Warnings
========

This class contains warnings settings. For now it contains only one configurable 
warning.

Attributes:
***********
* | **high_reaction_probability**: WarningLevel = WarningLevel.IGNORE
  | Print a warning when a bimolecular reaction probability is over 0.5 but less or equal than 1.
  | Warning when probability is greater than 1 is always printed.
  | Cannot be set to WarningLevel.ERROR.

  | Example: `0615_bimol_rxn_prob_over_05_less_1_warning_disabled/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/nutmeg4_pymcell4/0615_bimol_rxn_prob_over_05_less_1_warning_disabled/model.py>`_ 


