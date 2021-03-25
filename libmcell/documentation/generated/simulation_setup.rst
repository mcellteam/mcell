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

* | **rxn_probability_changed**: bool = True
  | When True, information that a reaction's probability has changed is printed during simulation.

Warnings
========

This class contains warnings settings. For now it contains only one configurable 
warning.

Attributes:
***********
* | **high_reaction_probability**: WarningLevel = WarningLevel.WARNING
  | Print a warning when a bimolecular reaction probability is over 0.5 but less or equal than 1.
  | Warning when probability is greater than 1 is always printed.
  | Cannot be set to WarningLevel.ERROR.

