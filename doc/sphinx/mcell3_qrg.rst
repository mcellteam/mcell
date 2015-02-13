.. title:: MCell Quick Reference Guide

In this document, the main text is in a sans-serif font. Command-line entries,
MDL file commands, and code is in a ``fixed-width`` font.  Values that must be
supplied by the user are in an *italicized sans-serif* font.

.. _running_mcell3:

Running MCell3
================

MCell3 runs on the command line. The format is

    ``mcell3`` *options* *filename*

By default, MCell3 sends informational messages, such as simulation progress,
to stdout (which will normally appear on the screen); error messages are sent
to stderr (which will also normally appear on the screen). Results of
simulations are written to files and do not appear as MCell3 is running.

A brief summary of MCell3 optional command-line arguments is given below.

.. tabularcolumns:: |l|L|

+-------------------------------------+---------------------------------------+
| **Argument**                        | **Explanation**                       |
+=====================================+=======================================+
| ``-seed`` *N*                       | Start with random number seed *N*     |
|                                     | instead of 1 (the default).           |
+-------------------------------------+---------------------------------------+
| ``-iterations`` *N*                 | Run the simulation with *N* timesteps |
|                                     | (overrides the value in the mdl file) |
+-------------------------------------+---------------------------------------+
| ``-help``                           | Print out a basic help screen.        |
+-------------------------------------+---------------------------------------+
| ``-logfile`` *filename*             | Send messages to *filename* instead   |
|                                     | of stdout/stderr                      |
+-------------------------------------+---------------------------------------+
| ``-errfile`` *filename*             | Send error messages to *filename*     |
|                                     | instead of stderr                     |
+-------------------------------------+---------------------------------------+
| ``-logfreq`` *N*                    | Print out a message when every *N*    |
|                                     | iterations have finished.             |
+-------------------------------------+---------------------------------------+
| ``-checkpoint_infile`` *filename*   | Use *filename* as a checkpoint file   |
|                                     | for the current simulation (overrides |
|                                     | any value in the MDL file).           |
+-------------------------------------+---------------------------------------+
| ``-with_checks``                    | Accepts 'yes' or 'no' as              |
|                                     | argument. If enabled, MCell will      |
|                                     | perform extended model checks. This   |
|                                     | option is enabled by default. Please  |
|                                     | note that model checking may have a   |
|                                     | noticeable run-time overhead. Thus,   |
|                                     | it may be advantageous to turn this   |
|                                     | option off after a given model has    |
|                                     | been checked once.                    |
+-------------------------------------+---------------------------------------+

.. _mdl_overview:

Model Description Language overview
=====================================

MCell3 runs simulations that are specified in *model description language*
(MDL) format. These files typically have the extension .mdl, but are not
required to. A MDL file is a text file with commands separated by whitespace.
The nature and type of whitespace (space, tab, newline) is unimportant to
MCell3. You are thus free to use whitespace to clarify the contents of the MDL
file.

.. _mdl_structure:

The structure of an MDL file
--------------------------------

Commands fall into five general groups, which usually should be given in the
order presented below. Although this is not always required, there are some
commands (e.g. defining a molecule) that must be used before others (e.g.
defining a reaction that uses that molecule). The order below should always be
safe:

#. Initialization. These commands set global parameters such as the
   time-step, spatial partitioning, and duration of the simulation.
#. Molecule definitions. These commands specify the names and diffusion
   constants of molecules in the simulation.
#. Reaction definitions. These commands specify the reactions that can
   occur between molecules and the rate at which those reactions occur.
#. Geometry specification. These commands describe the membranes and
   other boundaries within which the simulation occurs, plus where in
   the world to place molecules initially.
#. Output specification. These commands specify what data should be
   output as the simulation is running; this can include graphical
   snapshots of the simulation in progress, as well as lists of numbers
   of molecules or reactions as a function of time.

In addition, there are utility commands--defining variables and including other
MDL files--that can appear nearly anywhere.

.. _how_to_use:

How to use this document
----------------------------

This document gives a brief description of every valid MCell3 command.
Commands can be specified one after another; it is often convenient to put
commands on separate lines but this is not necessary.

Some commands have a scope, delimited by ``{`` and ``}`` (braces).  Within
these braces, a different set of commands become available. In this document,
each set of commands is given a different title. For example, commands given
within a ``DEFINE_MOLECULE`` block receive the title Define Molecule Commands.

.. _mdl_commands:

MDL commands
==============

.. _init_commands:

Initialization commands
---------------------------

The following initialization commands are required in every MDL file.

.. tabularcolumns:: |l|L|

+--------------------------------+----------------------------------------------------+
| **Command**                    | **Explanation**                                    |
+================================+====================================================+
| ``TIME_STEP =`` *t*            | Set the simulation time step to *t* seconds.       |
| :index:`\ <single:TIME_STEP>`  | ``1e-6`` is a common value. Later commands can     |
|                                | change the time steps taken by individual          |
|                                | molecules, but this time step is still used by all |
|                                | output statements.                                 |
+--------------------------------+----------------------------------------------------+
| ``ITERATIONS =`` *N*           | Run the simulation for *N* iterations.             |
| :index:`\ <single:ITERATIONS>` |                                                    |
+--------------------------------+----------------------------------------------------+

The following initialization commands are optional.

.. math::

.. tabularcolumns:: |p{7cm}|p{9cm}|

.. cssclass:: longtable

+------------------------------------------------+----------------------------------------------------+
| **Command**                                    | **Explanation**                                    |
+================================================+====================================================+
| ``TIME_STEP_MAX =`` *t*                        | MCell3 will move longer than the specified         |
| :index:`\ <single:TIME_STEP_MAX>`              | simulation time step if it seems safe. This        |
|                                                | command makes sure that the longest possible time  |
|                                                | step is no longer than *t* seconds, even if MCell3 |
|                                                | thinks a longer step would be safe. The default is |
|                                                | no limit.                                          |
+------------------------------------------------+----------------------------------------------------+
| ``SPACE_STEP =`` *N*                           | Have all diffusing molecules take time steps of    |
| :index:`\ <single:SPACE_STEP>`                 | different duration, chosen so that the mean        |
|                                                | diffusion distance is *N* microns for each         |
|                                                | molecule. By default, all molecules move the same  |
|                                                | time step.                                         |
+------------------------------------------------+----------------------------------------------------+
| ``CHECKPOINT_INFILE = "`` *filename* ``"``     | Start the simulation using the conditions          |
| :index:`\ <single:CHECKPOINT_INFILE>`          | specified in the checkpoint file *filename*. This  |
|                                                | will start at the time that the saved simulation   |
|                                                | left off, and will use molecules stored in the     |
|                                                | specified file instead of surface molecule         |
|                                                | densities/numbers specified in the MDL file.       |
|                                                | Release sites can add new molecules if the release |
|                                                | time is after the time the simulation starts.      |
+------------------------------------------------+----------------------------------------------------+
| ``CHECKPOINT_OUTFILE = "`` *filename* ``"``    | Save the state of the simulation when              |
| :index:`\ <single:CHECKPOINT_OUTFILE>`         | ``CHECKPOINT_ITERATIONS`` (described below) is     |
|                                                | reached, and stop.                                 |
+------------------------------------------------+----------------------------------------------------+
| ``CHECKPOINT_REALTIME =``                      | Create a checkpoint file after a specified period  |
| *time*  *exit_policy*                          | of *time* has elapsed. The *time* should be set to |
| :index:`\ <single:CHECKPOINT_REALTIME>`        | integer values separated by colons like this: 1:30 |
|                                                | (one minute, thirty seconds) or 1:5:0:0 (one day,  |
|                                                | five hours). Possible units and formatting are as  |
|                                                | follows: *days:hours:minutes:seconds*,             |
|                                                | *hours:minutes:seconds*, *minutes:seconds*, and    |
|                                                | *seconds*. The *exit_policy* is optional and can   |
|                                                | be set to ``EXIT`` or ``NOEXIT``. If set to        |
|                                                | ``EXIT``, then the simulation will stop after the  |
|                                                | checkpoint file is created. If set to ``NOEXIT``,  |
|                                                | then it will continue running. ``EXIT`` is the     |
|                                                | default.                                           |
+------------------------------------------------+----------------------------------------------------+
| ``CHECKPOINT_ITERATIONS =`` *N*                | Used with ``CHECKPOINT_OUTFILE``. This specifies   |
| :index:`\ <single:CHECKPOINT_ITERATIONS>`      | how many iterations to run before stopping and     |
|                                                | writing the checkpoint file. If *N* is larger than |
|                                                | ``ITERATIONS``, the simulation will terminate      |
|                                                | normally after the maximum amount of iterations as |
|                                                | specified by ``ITERATIONS`` has been reached.      |
+------------------------------------------------+----------------------------------------------------+
| ``SURFACE_GRID_DENSITY =`` *N*                 | Tile all surfaces so that they can hold molecules  |
| :index:`\ <single:SURFACE_GRID_DENSITY>`       | at *N* different positions per square micron. The  |
|                                                | default is 10000. For backwards compatibility,     |
|                                                | ``EFFECTOR_GRID_DENSITY`` works also.              |
+------------------------------------------------+----------------------------------------------------+
| ``INTERACTION_RADIUS =`` *N*                   | Diffusing volume molecules will interact with each |
| :index:`\ <single:INTERACTION_RADIUS>`         | other when they get within *N* microns of each     |
|                                                | other. The default is                              |
|                                                | :math:`1/\sqrt{\pi\cdot\sigma_s}` where            |
|                                                | :math:`\sigma_s` is the surface grid density       |
|                                                | (default or user-specified).                       |
+------------------------------------------------+----------------------------------------------------+
| ``PARTITION_`` *D* ``= [`` *list* ``]``        | Subdivide the *D* th axis of space, where *D* is   |
| :index:`\ <single:PARTITION_X>`                | ``X``, ``Y``, or ``Z``, at the boundaries given in |
| :index:`\ <single:PARTITION_Y>`                | *list* (in microns). In future versions, MCell3    |
| :index:`\ <single:PARTITION_Z>`                | will further subdivide space if it is              |
|                                                | computationally advantageous. By default, each     |
|                                                | axis will be split into between five and fifteen   |
|                                                | equal partitions. If you do not explicitly         |
|                                                | partition all three axes, MCell3 is likely to      |
|                                                | ignore your request and perform automatic          |
|                                                | partitioning. The spacing between adjacent         |
|                                                | partitions must be larger than the                 |
|                                                | ``INTERACTION_RADIUS``.                            |
+------------------------------------------------+----------------------------------------------------+
| ``RADIAL_DIRECTIONS =`` *N*                    | Specifies how many different directions to put in  |
| :index:`\ <single:RADIAL_DIRECTIONS>`          | the look-up table. The default is sensible. Don't  |
|                                                | use this unless you know what you're doing.        |
|                                                | Instead of a number, you can specify               |
|                                                | ``FULLY_RANDOM`` to generate the directions        |
|                                                | directly from double precision numbers (but this   |
|                                                | is slower).                                        |
+------------------------------------------------+----------------------------------------------------+
| ``RADIAL_SUBDIVISIONS =`` *N*                  | Specifies how many distances to put in the         |
| :index:`\ <single:RADIAL_SUBDIVISIONS>`        | diffusion look-up table. Again, the default is     |
|                                                | sensible. ``FULLY_RANDOM`` is not implemented      |
|                                                | here.                                              |
+------------------------------------------------+----------------------------------------------------+
| ``ACCURATE_3D_REACTIONS =`` *boolean*          | Specifies which method to use for computing 3D     |
| :index:`\ <single:ACCURATE_3D_REACTIONS>`      | molecule-molecule interactions. If *boolean* is    |
|                                                | ``TRUE``, then molecules will look through         |
|                                                | partition boundaries for potential interacting     |
|                                                | partners--this is slower but more accurate. If     |
|                                                | *boolean* is ``FALSE``, then molecule interaction  |
|                                                | disks will be clipped at partition boundaries and  |
|                                                | probabilities adjusted to get the correct rate--   |
|                                                | this is faster but can be less accurate. The       |
|                                                | default is ``TRUE``.                               |
+------------------------------------------------+----------------------------------------------------+
| ``CENTER_MOLECULES_ON_GRID =`` *boolean*       | If *boolean* is set to ``TRUE``, then all          |
| :index:`\ <single:CENTER_MOLECULES_ON_GRID>`   | molecules on a surface will be located exactly at  |
|                                                | the center of their grid element. If ``FALSE``,    |
|                                                | the molecules will be randomly located when        |
|                                                | placed, and reactions will take place at the       |
|                                                | location of the target (or the site of impact in   |
|                                                | the case of 3D molecule/surface reactions). The    |
|                                                | default is ``FALSE.``                              |
+------------------------------------------------+----------------------------------------------------+
| ``VACANCY_SEARCH_DISTANCE =`` *r*              | Normally, a reaction will not proceed on a surface |
| :index:`\ <single:VACANCY_SEARCH_DISTANCE>`    | unless there is room to place all products on the  |
|                                                | single grid element where the reaction is          |
|                                                | initiated. By increasing *r* from its default      |
|                                                | value of 0, one can specify how far from the       |
|                                                | reaction's location, in microns, the reaction can  |
|                                                | place its products. To be useful, *r* must be      |
|                                                | larger than the longest axis of the grid element   |
|                                                | on the triangle in question. The reaction will     |
|                                                | then proceed if there is room to place its         |
|                                                | products within a radius *r*, and will place those |
|                                                | products as close as possible to the place where   |
|                                                | the reaction occurs (deterministically, so small-  |
|                                                | scale directional bias is possible).               |
+------------------------------------------------+----------------------------------------------------+
| ``MICROSCOPIC_REVERSIBILITY =`` *value*        | If *value* is set to ``OFF``, then binding-        |
| :index:`\ <single:MICROSCOPIC_REVERSIBILITY>`  | unbinding reactions between molecules will be      |
|                                                | somewhat more efficient but may not be accurate if |
|                                                | the probability of binding is high (close to 1).   |
|                                                | If ``ON``, a more computationally demanding        |
|                                                | routine will be used to make sure binding-         |
|                                                | unbinding is more similar in both directions. If   |
|                                                | *value* is set to ``SURFACE_ONLY`` or              |
|                                                | ``VOLUME_ONLY``, the more accurate routines will   |
|                                                | be used only for reactions at surfaces or only for |
|                                                | those in the volume. ``OFF`` is the default.       |
+------------------------------------------------+----------------------------------------------------+
| | ``NOTIFICATIONS``                            | This block of commands lets you set the            |
| | ``{``                                        | informational messages that MCell3 generates. The  |
| |   *notification commands*                    | block can appear multiple times and applies to all |
| | ``}``                                        | MDL below it in the file. It can appear anywhere   |
|                                                | at the top level (but not inside other blocks).    |
| :index:`\ <single:NOTIFICATIONS>`              |                                                    |
+------------------------------------------------+----------------------------------------------------+
| | ``WARNINGS``                                 | This block of commands lets you control how MCell3 |
| | ``{``                                        | handles warnings---whether it generates a warning  |
| |   *warning policy commands*                  | and continues, silently handles the condition, or  |
| | ``}``                                        | generates an error and quits. The block can appear |
|                                                | multiple times and applies to all MDL below it in  |
| :index:`\ <single:WARNINGS>`                   | the file. It can appear anywhere at the top level  |
|                                                | (but not inside other blocks).                     |
+------------------------------------------------+----------------------------------------------------+

The following commands can be given in a notifications block; in each case,
setting the notification policy to ``OFF`` will prevent any informational
output regarding that aspect of the simulation. This will not affect warnings.

.. math::

.. tabularcolumns:: |l|L|

+----------------------------------------------------+---------------------------------------+
| **Notification Command**                           | **Explanation**                       |
+====================================================+=======================================+
|  ``BOX_TRIANGULATION_REPORT =`` *policy*           | If *policy* is ``ON``, MCell3 will    |
|  :index:`\ <single:BOX_TRIANGULATION_REPORT>`      | report how many triangles are         |
|                                                    | generated from each box object.       |
|                                                    | Default is ``OFF``.                   |
+----------------------------------------------------+---------------------------------------+
| ``DIFFUSION_CONSTANT_REPORT =`` *policy*           | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:DIFFUSION_CONSTANT_REPORT>`      | report four measures of the diffusion |
|                                                    | constant for each molecule. If        |
|                                                    | *policy* is ``BRIEF``, MCell3 will    |
|                                                    | report just one measure (average      |
|                                                    | diffusion distance per step) for each |
|                                                    | molecule. Default is ``BRIEF``.       |
+----------------------------------------------------+---------------------------------------+
| ``FILE_OUTPUT_REPORT =`` *policy*                  | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:FILE_OUTPUT_REPORT>`             | report every time reaction data is    |
|                                                    | written to disk. Default is ``OFF``.  |
+----------------------------------------------------+---------------------------------------+
| ``FINAL_SUMMARY =`` *policy*                       | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:FINAL_SUMMARY>`                  | give some information about the CPU   |
|                                                    | time used and some of the internal    |
|                                                    | events. Default is ``ON``.            |
+----------------------------------------------------+---------------------------------------+
| ``ITERATION_REPORT =`` *policy*                    | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:ITERATION_REPORT>`               | provide a running report of how many  |
|                                                    | iterations have completed, chosen     |
|                                                    | based on the total number of          |
|                                                    | iterations. If *policy* is an integer |
|                                                    | value, MCell3 will report each time   |
|                                                    | that number of iterations have        |
|                                                    | elapsed. Default is ``ON``.           |
+----------------------------------------------------+---------------------------------------+
| ``PARTITION_LOCATION_REPORT =`` *policy*           | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:PARTITION_LOCATION_REPORT>`      | print out the locations of the        |
|                                                    | partitions used for the simulation.   |
|                                                    | Default is ``OFF``.                   |
+----------------------------------------------------+---------------------------------------+
| ``PROBABILITY_REPORT =`` *policy*                  | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:PROBABILITY_REPORT>`             | print out the reaction probabilities  |
|                                                    | for each reaction (except special     |
|                                                    | internal surface reactions such as    |
|                                                    | absorptive surfaces). Default is      |
|                                                    | ``ON``. This will reset the reporting |
|                                                    | threshold to a probability of zero.   |
+----------------------------------------------------+---------------------------------------+
| ``PROBABILITY_REPORT_THRESHOLD =`` *p*             | MCell3 will print out the             |
| :index:`\ <single:PROBABILITY_REPORT_THRESHOLD>`   | probabilities for every reaction with |
|                                                    | probability greater than or equal to  |
|                                                    | *p*. This will override the policy    |
|                                                    | for probability reports.              |
+----------------------------------------------------+---------------------------------------+
| ``VARYING_PROBABILITY_REPORT =`` *policy*          | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:VARYING_PROBABILITY_REPORT>`     | print out the reaction probabilities  |
|                                                    | when a time- varying reaction updates |
|                                                    | its reaction rate (regardless of the  |
|                                                    | old or new probability). Default is   |
|                                                    | ``ON``.                               |
+----------------------------------------------------+---------------------------------------+
| ``PROGRESS_REPORT =`` *policy*                     | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:PROGRESS_REPORT>`                | print out messages indicating which   |
|                                                    | part of the simulation process is     |
|                                                    | underway (initializing, running,      |
|                                                    | etc.). Default is ``ON``.             |
+----------------------------------------------------+---------------------------------------+
| ``RELEASE_EVENT_REPORT =`` *policy*                | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:RELEASE_EVENT_REPORT>`           | print out a message every time        |
|                                                    | molecules are released through a      |
|                                                    | release site (indicating how many     |
|                                                    | molecules of which type were released |
|                                                    | and the iteration on which they were  |
|                                                    | released). Default is ``ON``.         |
+----------------------------------------------------+---------------------------------------+
| ``MOLECULE_COLLISION_REPORT =`` *policy*           | If *policy* is ``ON``, MCell3 will    |
| :index:`\ <single:MOLECULE_COLLISION_REPORT>`      | print, for each reaction type, the    |
|                                                    | number of bimolecular or trimolecular |
|                                                    | collisions that occured between       |
|                                                    | reactants during reactions. Default   |
|                                                    | is ``OFF``.                           |
+----------------------------------------------------+---------------------------------------+
| ``ALL_NOTIFICATIONS =`` *policy*                   | Set all notification policies to the  |
| :index:`\ <single:ALL_NOTIFICATIONS>`              | same value (``ON`` or ``OFF``). This  |
|                                                    | overrides the existing probability    |
|                                                    | report threshold, if there is one.    |
+----------------------------------------------------+---------------------------------------+

The following commands can be given in a warnings block. Setting the warning
policy to ``IGNORED`` will prevent any output and the condition will be handled
as best it can. ``WARNING`` will give a warning message, but the problem will
be handled and the simulation will continue.  Setting to ``ERROR`` will
generate an error and the simulation will stop. This will not affect
notification policies.

.. math::

+--------------------------------------------------+-----------------------------------------+
| **Warning Policy Command**                       | **Explanation**                         |
+==================================================+=========================================+
| | ``DEGENERATE_POLYGONS =`` *policy*             | Degenerate polygons are polygons with   |
| | :index:`\ <single:DEGENERATE_POLYGONS>`        | zero area and must be removed for the   |
|                                                  | simulation to run. The default policy   |
|                                                  | is ``WARNING``.                         |
+--------------------------------------------------+-----------------------------------------+
| ``HIGH_REACTION_PROBABILITY =`` *policy*         | Generate warnings or errors if reaction |
| :index:`\ <single:HIGH_REACTION_PROBABILITY>`    | probabilities exceed a certain          |
|                                                  | threshold. The default policy is        |
|                                                  | ``IGNORED``. The warnings or errors     |
|                                                  | will be generated both at parse time    |
|                                                  | and during run-time if there are time   |
|                                                  | varying reaction rates that exceed the  |
|                                                  | threshold.                              |
+--------------------------------------------------+-----------------------------------------+
| ``HIGH_PROBABILITY_THRESHOLD =`` *p*             | If the policy is to generate warnings   |
| :index:`\ <single:HIGH_PROBABILITY_THRESHOLD>`   | or errors on high probability           |
|                                                  | reactions, have them generated when the |
|                                                  | probability equals or exceeds *p*. The  |
|                                                  | default value is 1.0.                   |
+--------------------------------------------------+-----------------------------------------+
| ``LIFETIME_TOO_SHORT =`` *policy*                | Generate warnings if molecules have     |
| :index:`\ <single:LIFETIME_TOO_SHORT>`           | short lifetimes (which could affect the |
|                                                  | accuracy of the simulation). This       |
|                                                  | warning occurs after the simulation has |
|                                                  | ended, so ``ERROR``. is not a valid     |
|                                                  | option. The default policy is           |
|                                                  | ``WARNING``.                            |
+--------------------------------------------------+-----------------------------------------+
| ``LIFETIME_THRESHOLD =`` *n*                     | If the policy is to generate a warning  |
| :index:`\ <single:LIFETIME_THRESHOLD>`           | if molecules have short lifetimes, then |
|                                                  | generate warnings on molecules that     |
|                                                  | have an average lifetime of less than   |
|                                                  | *n* iterations. The default value is    |
|                                                  | 50.                                     |
+--------------------------------------------------+-----------------------------------------+
| ``MISSED_REACTIONS =`` *policy*                  | Generate errors or warnings if there    |
| :index:`\ <single:MISSED_REACTIONS>`             | are missed reactions (which usually is  |
|                                                  | a consequence of an overly high         |
|                                                  | reaction probability). This warning     |
|                                                  | occurs after the simulation has ended,  |
|                                                  | so ``ERROR``. is not a valid option.    |
|                                                  | The default policy is ``WARNING``.      |
+--------------------------------------------------+-----------------------------------------+
| ``MISSED_REACTION_THRESHOLD =`` *f*              | If the policy is to generate a warning  |
| :index:`\ <single:MISSED_REACTION_THRESHOLD>`    | if there are missed reactions, then     |
|                                                  | generate a warning for each reaction    |
|                                                  | where a fraction of at least *f* of     |
|                                                  | reactions were missed. The default      |
|                                                  | value is :math:`10^{-3}`.               |
+--------------------------------------------------+-----------------------------------------+
| ``NEGATIVE_DIFFUSION_CONSTANT =`` *policy*       | Diffusion constants cannot be negative, |
| :index:`\ <single:NEGATIVE_DIFFUSION_CONSTANT>`  | and will be set to zero if they are.    |
|                                                  | The default policy is ``WARNING``.      |
+--------------------------------------------------+-----------------------------------------+
| ``MISSING_SURFACE_ORIENTATION =`` *policy*       | Generate errors or warnings if a        |
| :index:`\ <single:MISSING_SURFACE_ORIENTATION>`  | molecule is placed on a surface or      |
|                                                  | reactions occur at a surface without a  |
|                                                  | specified orientation---the code will   |
|                                                  | assume you mean that there is no        |
|                                                  | orientation in the warning or silent    |
|                                                  | cases. To avoid triggering this         |
|                                                  | condition, if you want to have no       |
|                                                  | orientation, you must specify it        |
|                                                  | explicitly with ``',`` or ``,'`` or     |
|                                                  | ``;``. The default policy is ``ERROR``. |
+--------------------------------------------------+-----------------------------------------+
| ``NEGATIVE_REACTION_RATE =`` *policy*            | Reaction rate constants cannot be       |
| :index:`\ <single:NEGATIVE_REACTION_RATE>`       | negative, and will be set to zero if    |
|                                                  | they are. The default policy is         |
|                                                  | ``WARNING``.                            |
+--------------------------------------------------+-----------------------------------------+
| ``USELESS_VOLUME_ORIENTATION =`` *policy*        | Generate errors or warnings if a        |
| :index:`\ <single:USELESS_VOLUME_ORIENTATION>`   | molecule is placed in a volume or       |
|                                                  | reactions occur in free space but an    |
|                                                  | orientation is specified anyway---      |
|                                                  | there is no way to impose orientation   |
|                                                  | so the marks will be ignored. The       |
|                                                  | default policy is ``WARNING``.          |
+--------------------------------------------------+-----------------------------------------+
| ``ALL_WARNINGS =`` *policy*                      | Set all warning policies to the same    |
| :index:`\ <single:ALL_WARNINGS>`                 | value (``IGNORED``, ``WARNING`` or      |
|                                                  | ``ERROR``). If ``ERROR`` is not a valid |
|                                                  | choice, the policy will be set to       |
|                                                  | ``WARNING`` instead.                    |
+--------------------------------------------------+-----------------------------------------+

.. _molecule_def_commands:

Molecule definition commands
--------------------------------

All molecules must be defined by name in a ``DEFINE_MOLECULES`` block. The
names must be unique in the entire simulation (that is, unique within their own
MDL file and any included MDL files that make up the whole simulation).

A define molecule block can be one of the following:

.. tabularcolumns:: |p{7cm}|p{9cm}|

+---------------------------------+-------------------------------------------+
| **Command**                     | **Explanation**                           |
+=================================+===========================================+
| | ``DEFINE_MOLECULE`` *name*    | Define a single molecule called *name*.   |
| | ``{``                         | The molecule's properties are specified   |
| |    *define molecule commands* | by commands inside braces.                |
| | ``}``                         |                                           |
+---------------------------------+-------------------------------------------+

.. tabularcolumns:: |p{7cm}|p{9cm}|

+-----------------------------------------------------+--------------------------+
| **Command**                                         | **Explanation**          |
+=====================================================+==========================+
| | ``DEFINE_MOLECULES``                              | Define a series of       |
| | ``{``                                             | molecules by name. Each  |
| |    *nameA* ``{`` *define molecule commands* ``}`` | molecule's properties are|
| |    *nameB* ``{`` *define molecule commands* ``}`` | specified by commands    |
| |    *  ...*                                        | inside braces.           |
| | ``}``                                             |                          |
+-----------------------------------------------------+--------------------------+

Each molecule must have a diffusion constant set using one of the following
commands:

.. math::

.. tabularcolumns:: |p{7cm}|p{9cm}|

+-------------------------------------------+-------------------------------------------+
| **Define Molecule Command**               | **Explanation**                           |
+===========================================+===========================================+
| | ``DIFFUSION_CONSTANT =`` *D*            | This molecule diffuses in space with      |
| | :index:`\ <single:DIFFUSION_CONSTANT>`  | diffusion constant *D*. *D* can be zero,  |
|                                           | in which case the molecule doesn't        |
|                                           | move. Synonyms for this command are       |
|                                           | ``DIFFUSION_CONSTANT_3D`` and ``D_3D``.   |
|                                           | The units of *D* are :math:`cm^2/s`.      |
+-------------------------------------------+-------------------------------------------+
| ``DIFFUSION_CONSTANT_2D =`` *D*           | This molecule is constrained to a surface |
| :index:`\ <single:DIFFUSION_CONSTANT_2D>` | and diffuses with diffusion constant *D*. |
|                                           | ``D_2D`` is a synonym for this command.   |
+-------------------------------------------+-------------------------------------------+

The following optional commands can be applied to each molecule (and must
appear in this order, and after the diffusion constant is set):

.. math::

.. tabularcolumns:: |p{5cm}|p{11cm}|

+-------------------------------------------+---------------------------------------------+
| **Define Molecule Command**               | **Explanation**                             |
+===========================================+=============================================+
| | ``CUSTOM_TIME_STEP =`` *t*              | This molecule should take timesteps of      |
| | :index:`\ <single:CUSTOM_TIME_STEP>`    | length *t* (in seconds). Use either this or |
|                                           | ``CUSTOM_SPACE_STEP``, not both.            |
+-------------------------------------------+---------------------------------------------+
| ``CUSTOM_SPACE_STEP =`` *L*               | This molecule should take steps of average  |
| :index:`\ <single:CUSTOM_SPACE_STEP>`     | length *L* (in microns). If you use this    |
|                                           | directive, do not set ``CUSTOM_TIME_STEP``. |
|                                           | Providing a ``CUSTOM_SPACE_STEP`` for a     |
|                                           | molecule overrides a potentially present    |
|                                           | global ``SPACE_STEP`` for this particular   |
|                                           | molecule.                                   |
+-------------------------------------------+---------------------------------------------+
| ``TARGET_ONLY``                           | This molecule will not initiate reactions   |
| :index:`\ <single:TARGET_ONLY>`           | when it runs into other molecules. This     |
|                                           | setting can speed up simulations when       |
|                                           | applied to a molecule at high               |
|                                           | concentrations that reacts with a molecule  |
|                                           | at low concentrations (it is more efficient |
|                                           | for the low-concentration molecule to       |
|                                           | trigger the reactions). This directive does |
|                                           | not affect unimolecular reactions.          |
+-------------------------------------------+---------------------------------------------+
| ``MAXIMUM_STEP_LENGTH =`` *L*             | This molecule should never step farther     |
| :index:`\ <single:MAXIMUM_STEP_LENGTH>`   | than length *L* (in microns) during a       |
|                                           | single timestep. This can be used to speed  |
|                                           | up simulations by enforcing a certain       |
|                                           | maximum step length for molecules such as   |
|                                           | molecular motors on a surface without       |
|                                           | having to reduce the global timestep        |
|                                           | unnecessarily. Please use this keyword with |
|                                           | care since it may give rise to a            |
|                                           | non-equilibrium distribution of the given   |
|                                           | molecule and also cause deviations from     |
|                                           | mass action kinetics.                       |
+-------------------------------------------+---------------------------------------------+

.. _rxn_def_commands:

Reaction definition commands
--------------------------------

All reactions must be defined inside a reaction definition block:

.. tabularcolumns:: |p{6cm}|p{10cm}|

+-----------------------------+---------------------------------------------+
| **Command**                 | **Explanation**                             |
+=============================+=============================================+
|  | ``DEFINE_REACTIONS``     | Define a series of reactions inside braces. |
|  | ``{``                    |                                             |
|  | *  reaction commands*    |                                             |
|  | ``}``                    |                                             |
|                             |                                             |
+-----------------------------+---------------------------------------------+

Reactions are specified using arrow notation:

.. tabularcolumns:: |p{6cm}|p{10cm}|

+------------------------------+----------------------------------------------+
| **Reaction Command**         | **Explanation**                              |
+==============================+==============================================+
| *reactants* ``->``           | Define a reaction that occurs between one,   |
| *products* ``[``\ *rate*\    | two or three *reactants* (names of           |
| ``]``                        | molecules, separated by ``+``) and produces  |
|                              | an arbitrary number of *products* (also      |
|                              | separated by ``+``), with a specified        |
|                              | *rate*. If a molecule is in the *reactants*  |
|                              | list and not in the *products* list, it is   |
|                              | destroyed in the reaction. *rate* can either |
|                              | be a literal number or a filename, in        |
|                              | quotes, that contains two columns: the       |
|                              | second is the rate, while the first is the   |
|                              | time at which that rate should start being   |
|                              | used. This allows variable reaction rates.   |
|                              | If you do not want products, use the         |
|                              | ``NULL`` keyword as a placeholder.           |
+------------------------------+----------------------------------------------+
| *reactants* ``->``           | As above, and call the reaction *name* so it |
| *products* ``[``\ *rate*\    | can be referred to by count statements.      |
| ``]:``\ *name*               |                                              |
+------------------------------+----------------------------------------------+

The units of the reaction *rate* for uni- and bimolecular reactions are

-  [:math:`s^{-1}` ] for unimolecular reactions,
-  [:math:`M^{-1}s^{-1}` ] for bimolecular reactions between
   either two volume molecules or a volume molecule and a surface
   (molecule), and
-  [:math:`{\mu}m^2N^{-1}s^{-1}`] for bimolecular reactions
   between two surface molecules.

Here, M is the molarity of the solution and N the number of reactants.

This notation is perhaps best explained through examples. In the most basic
form, reactants and products are just the names of molecules, separated by
``+``:

.. math::

+------------------------------+----------------------------------------------+
| **Example**                  | **Explanation**                              |
+==============================+==============================================+
| ``A -> B [100]``             | Molecule ``A`` changes into molecule ``B``   |
|                              | at a rate of :math:`100 s^{-1}`.             |
+------------------------------+----------------------------------------------+
| ``A -> A + B [100]``         | Molecule ``A`` emits molecules of ``B`` at a |
|                              | rate of :math:`100 s^{-1}`.                  |
+------------------------------+----------------------------------------------+
| ``A -> NULL [100]``          | Molecule ``A`` is destroyed at a rate of     |
|                              | :math:`100 s^{-1}`.                          |
+------------------------------+----------------------------------------------+
| ``A + B -> A [1e6]``         | Molecule ``A`` destroys molecule ``B`` at a  |
|                              | rate of :math:`10^6M^{-1}s^{-1}`.            |
+------------------------------+----------------------------------------------+
| ``A + B -> A + C [1e6]``     | Molecule ``A`` catalytically converts ``B``  |
|                              | to ``C`` at a rate of                        |
|                              | :math:`10^6M^{-1}s^{-1}`                     |
+------------------------------+----------------------------------------------+
| ``A + B -> A + B + C [1e6]`` | Collision of ``A`` and ``B`` catalytically   |
|                              | generates ``C`` at a rate of                 |
|                              | :math:`10^6M^{-1}s^{-1}`.                    |
+------------------------------+----------------------------------------------+

Reactions can take place on surfaces or involve molecules contained therein
(surface molecules). Surfaces possess a front and a back side defined by the
direction of the surface normal which points from the back toward the front.
Surface molecules have an orientation in the form of a top and a bottom domain
and are positioned on surfaces with their top domain either on the surfaces'
front or back side, or top-front and top-back for short.

Reactions that explicitly involve surfaces are said to occur with an absolute
orientation regarding the surface. When reactions involving surface molecules
take place in the absence of explicit surfaces they are said to occur without
an absolute orientation. Below, we will illustrate both cases.

.. _rxn_wo_absolute_orient:

Reactions without absolute orientation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For reactions without an absolute orientation, the reaction specification lists
the required relative orientation of the reactants and products. This allows
one to write general reactions that do not depend on the way in which molecules
are inserted into surfaces, i.e., either top-front or top-back.

The two possible orientations are specified by ``'`` and ``,`` (apostrophe and
comma) after the molecule's name. Hence, a surface-bound molecule ``B`` can
have the orientations ``B'`` and ``B,``. The table below provides a few example
reactions

.. math::

.. tabularcolumns:: |l|L|

+-------------------------------+---------------------------------------------+
| **Example**                   | **Explanation**                             |
+===============================+=============================================+
| ``B' -> B, [10]``             | Molecule ``B`` flips (changes its           |
|                               | orientation) at a rate of :math:`10 s^{-1}` |
+-------------------------------+---------------------------------------------+
| ``B' -> B' + A' + C,[10]``    | Molecule ``B`` emits molecules of ``A`` on  |
|                               | the side it's pointing to and emits ``C``   |
|                               | on the other side, at a rate of             |
|                               | :math:`10 s^{-1}`                           |
+-------------------------------+---------------------------------------------+
| ``B, -> B, + A, + C' [10]``   | This specifies exactly the same reaction as |
|                               | above. ``B`` and ``A`` end up with the same |
|                               | orientation, while ``C`` has opposite       |
|                               | orientation.                                |
+-------------------------------+---------------------------------------------+

The best way to keep the relationships straight is to draw a "before" picture
with each reactant facing the direction of the tick mark, and an "after"
picture with each product facing in the direction of the tick mark. Clearly,
inverting this picture by flipping all tick marks results in the same reaction.
One can thus use tick marks that are consistent with ones mental picture.

Below are additional reaction examples involving a molecule ``A`` diffusing in
3D and surface molecules ``B`` and ``C``:

.. math::

.. tabularcolumns:: |l|L|

+-------------------------+---------------------------------------------------+
| **Example**             | **Explanation**                                   |
+=========================+===================================================+
| ``A' + B' -> C' [1e5]`` | Molecule ``A`` binds to ``B`` if it is on the     |
|                         | side that ``B`` is pointing to, producing a ``C`` |
|                         | facing the same way as ``B``, at a rate of        |
|                         | :math:`10^5M^{-1}s^{-1}`.                         |
+-------------------------+---------------------------------------------------+
| ``A, + B, -> C, [1e5]`` | The same reaction again---everything occurs on    |
|                         | the same side, but we wrote it on the bottom this |
|                         | time.                                             |
+-------------------------+---------------------------------------------------+
| ``A' + B, -> C' [1e5]`` | Molecule ``A`` binds when it hits the opposite    |
|                         | side of ``B``, producing a ``C`` facing the       |
|                         | opposite way as ``B`` (i.e. towards the side      |
|                         | ``A`` came from), at a rate of                    |
|                         | :math:`10^5M^{-1}s^{-1}`.                         |
+-------------------------+---------------------------------------------------+
| ``A, + B' -> C, [1e5]`` | Same as above.                                    |
+-------------------------+---------------------------------------------------+

So far, all examples have used the first orientation class, specified with
``'`` and ``,``. The second orientation class is specified by ``''`` and
``,,``. The third is ``'''`` and ``,,,`` and so on.  Molecules in different
orientation classes do not pay attention to each other's orientation. In a
reaction with orientation, every molecule must be explicitly given an
orientation class otherwise an error is generated. This behavior can be
adjusted to generate warnings or no messages instead; in this case, molecules
without an orientation class act without regard to orientation. Several
examples follow:

.. math::

.. tabularcolumns:: |l|L|

+------------------------------------+----------------------------------------+
| **Example**                        | **Explanation**                        |
+====================================+========================================+
| ``A'' + B, -> C' [1e5]``           | Molecule ``A`` binds to either side of |
|                                    | ``B`` (since they are in different     |
|                                    | orientation classes); this produces a  |
|                                    | ``C`` facing the opposite way as       |
|                                    | ``B``, at a rate of                    |
|                                    | :math:`10^5M^{-1}s^{-1}`.              |
+------------------------------------+----------------------------------------+
| ``A,, + B, -> C' [1e5]``           | This is the same reaction - since      |
|                                    | ``A`` is the only molecule in the      |
|                                    | second orientation class, it doesn't   |
|                                    | matter which way we specify things.    |
+------------------------------------+----------------------------------------+
| ``A,, + B' -> C, [1e5]``           | Same again--``B`` and ``C`` still have |
|                                    | opposite orientations.                 |
+------------------------------------+----------------------------------------+
| ``A, + B' -> C,, ``[1e5]``         | Molecule ``A`` hits the opposite side  |
|                                    | of ``B`` and produces ``C`` that is    |
|                                    | equally likely to point either way, at |
|                                    | a rate of :math:`10^5M^{-1}s^{-1}`.    |
+------------------------------------+----------------------------------------+
| ``A, + B' -> C'' [1e5]``           | Same as above, since ``C`` is still    |
|                                    | not in the same orientation class as   |
|                                    | the others.                            |
+------------------------------------+----------------------------------------+
| ``A' + B'' -> A, + B''' [1e5]``    | Molecule ``A`` hits molecule ``B`` on  |
|                                    | either side; ``A`` keeps traveling     |
|                                    | (goes to the other side) and ``B``     |
|                                    | tumbles to a random orientation, at a  |
|                                    | rate of :math:`10^5M^{-1}s^{-1}`       |
+------------------------------------+----------------------------------------+
| ``A' + B'' -> C''' + D'''' [1e5]`` | ``A`` and ``B`` react in any           |
|                                    | orientation and produce ``C`` and      |
|                                    | ``D`` in random orientations. All      |
|                                    | orientation classes are different, so  |
|                                    | there are no geometrical constraints   |
|                                    | here.                                  |
+------------------------------------+----------------------------------------+

There are more examples of how one would use this syntax to model well-known
biological reactions at the end of this document in section
:ref:`example_models`.

.. _rxns_w_absolute_orient:

Reactions with absolute orientation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reactions can specify an absolute orientation with respect to the surface on
which they take place via including a surface class specification in the
reaction definition. The general form for defining reactions with absolute
orientations is accomplished via the "@" character as shown below

.. tabularcolumns:: |p{8cm}|p{8cm}|

+-------------------------------------+---------------------------------------------------------------+
| **Reaction Command**                | **Explanation**                                               |
+=====================================+===============================================================+
| *reactants* ``@`` *surf_class_name* | Define a reaction that occurs between one or two oriented     |
| ``->`` *products* ``[`` *rate*      | *reactants* (names of molecules, separated by ``+``) on a set |
| ``]``                               | of surface regions identified by *surf_class_name.* The       |
|                                     | reaction produces an arbitrary number of oriented *products*  |
|                                     | (also separated by ``+``), with a specified *rate*. If a      |
|                                     | molecule is in the *reactants* list and not in the *products* |
|                                     | list, it is destroyed in the reaction. The rate can also be a |
|                                     | filename, in quotes, that contains two columns: the second is |
|                                     | the rate, while the first is the time at which that rate      |
|                                     | should start being used. This allows variable reaction rates. |
|                                     | If you do not want products, use the ``NULL`` keyword as a    |
|                                     | placeholder.                                                  |
+-------------------------------------+---------------------------------------------------------------+
| *reactants* ``@`` *surf_class_name* | As above, and call the reaction *name* so it can be referred  |
| ``->`` *products* ``[`` *rate*      | to by count statements.                                       |
| ``]:`` *name*                       |                                                               |
+-------------------------------------+---------------------------------------------------------------+

A reaction defined in this way takes place on all surface regions which specify
``SURFACE_CLASS`` *= surf_class_name.* The relative orientation of reactants
and products is specified as explained in :ref:`rxn_wo_absolute_orient` but now
the reaction takes place with respect to the orientation given for
*surf_class_name* indicating the front or back of the selected surface regions.
Please note that all reactants have to be listed to the left of
*surf_class_name* and no surface class specifications can occur on the product
side of the reaction definition. Furthermore, for bi-molecular reactions at
least one of the two reactants has to be a surface molecule.

The table below lists several examples of oriented reactions involving a
surface class *surf*, a 3D molecule ``A``, and surface molecules ``B`` and
``C``.

.. math::

.. tabularcolumns:: |l|L|

+----------------------------------+--------------------------------------------+
| **Example**                      | **Explanation**                            |
+==================================+============================================+
| ``A' + B' @ surf' -> C, [1e5]``  | The reaction affects surface molecules     |
|                                  | ``B`` located on surface regions           |
|                                  | identified by surface class ``surf`` which |
|                                  | have their top domain at the front of      |
|                                  | the surface. ``B`` reacts with ``A``       |
|                                  | approaching from the front at a rate of    |
|                                  | :math:`10^5M^{-1}s^{-1}` to yield          |
|                                  | surface molecule ``C`` whose orientation   |
|                                  | is flipped with respect to ``B``, i.e.,    |
|                                  | ``C`` has its top domain aligned to the    |
|                                  | back of the surface regions.               |
+----------------------------------+--------------------------------------------+
| ``A' + B, @ surf' -> C, [1e5]``  | Same as above, but ``B`` now has its top   |
|                                  | domain at the back of the surface and      |
|                                  | reaction product ``C`` assumes the same    |
|                                  | orientation.                               |
+----------------------------------+--------------------------------------------+
| ``A,, + B, @ surf' -> C' [1e5]`` | Since ``A`` is in an orientation class     |
|                                  | different from both ``B`` and ``surf``,    |
|                                  | ``A`` can react from both sides. ``B``     |
|                                  | has its top domain at the back of the      |
|                                  | surface and the reaction product ``C``     |
|                                  | has its orientation flipped, i.e., its     |
|                                  | top domain is at the front of the          |
|                                  | surface.                                   |
+----------------------------------+--------------------------------------------+
| ``A' + B' @ surf' -> C,, [1e5]`` | Same as the the first reaction, but        |
|                                  | since product ``C`` is in a orientation    |
|                                  | class different from either ``A``,         |
|                                  | ``B``, and ``surf``, its orientation is    |
|                                  | random with respect to the surface         |
|                                  | regions, i.e., its top domain can be       |
|                                  | either on the front or back.               |
+----------------------------------+--------------------------------------------+

Tick marks add, so that ``',`` and ``,'`` mean no orientation. Reactions will
occur from either orientation when given reactants with no orientation, and
products will orient randomly. A semicolon, ``;``, can be used instead of two
opposite tick marks. Orientations can also be specified numerically inside
``{}`` after the molecule name. For example, ``A{1}`` and ``A{-1}`` are
synonyms for ``A'`` and ``A,`` and ``A{0}`` is a synonym for ``A;.``

There are several variants of the normal reaction arrow ->. One can use an
arbitrary number of dashes in the arrow, i.e., ``->,`` ``-->,`` and ``------>``
all mean the same thing. In addition, the following arrows have different
meanings:

.. math::

.. tabularcolumns:: |l|L|

+-----------------------+-----------------------------------------------------+
| **Reaction Arrow**    | **Explanation**                                     |
+=======================+=====================================================+
| ``->``                | A unidirectional reaction going from reactants (on  |
|                       | the left) to products (on the right).               |
+-----------------------+-----------------------------------------------------+
| ``<->``               | A bidirectional reaction going in either direction; |
|                       | at most two molecule names can appear on each side. |
|                       | A rate must be given for each direction using the   |
|                       | notation :math:`[>k_{+}, <k_{-}]`, where            |
|                       | :math:`k_{+}` is the forward rate constant and      |
|                       | :math:`k_{-}` is the backward rate constant.        |
+-----------------------+-----------------------------------------------------+
| *reactant* ``--``     | This specifies a catalytic reaction where           |
| *catalyst* ``->``     | *reactant* is converted to *products* in the        |
| *products*            | presence of *catalyst*. This is the same as the     |
|                       | reaction *catalyst* + *reactant* -> *catalyst* +    |
|                       | *products*. Presently, there can only be one        |
|                       | reactant.                                           |
+-----------------------+-----------------------------------------------------+
| *reactant* ``<-``     | A bidirectional catalytic reaction. There can only  |
| *catalyst* ``->``     | be one reactant and one product.                    |
| *product*             |                                                     |
+-----------------------+-----------------------------------------------------+

Finally, a few special cases deserve particular mention

*  For catalytic reactions, if a catalyst is a surface class, the latter is not copied to the list of products, i.e.:

  * ``A' — SURF' -> C, [rate]`` is equivalent to
  * ``A'  @ SURF' -> C, [rate]``

*  Reversible reactions of the form  ``A' @ SURF' <--> C, [>rate1,<rate2]  `` or ``A' <-- SURF'--> C, [>rate1,<rate2]``   are equivalent to the following two reactions:

  * ``A' @ SURF' -> C, [rate1]``
  * ``C, @ SURF' -> A' [rate2]``

.. _trimolecular_rxns:

Trimolecular reactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the conventional unimolecular and bimolecular reaction syntax,
users can also specify trimolecular reactions between arbitrary combinations of
volume and surface molecules, i.e., reactions of the form ``A + B + C ->
products`` with ``A``, ``B``, and ``C`` either volume or surface molecules. As
for regular unimolecular and bimolecular reactions, the presence of surface
molecules in a trimolecular reaction requires the addition of tick marks to
specify their proper orientation. Please note that the trimolecular reaction
syntax does not allow for the presence of an additional surface class specifier
via the ``@`` syntax. The ability to formulate trimolecular reactions within
MCell3 is targeted toward users who wish to use MCell3 to simulate ODE based
models which may contain such trimolecular terms. Please note that since
intermediate species are not explicitly treated, trimolecular reactions are
only approximations to the true underlying microscopic reaction mechanism and
faithfully represent the latter only over a limited parameter range. In
general, it is preferable to describe models using elementary reaction
mechanisms via unimolecular and bimolecular reactions.

Below are a few examples of trimolecular reactions involving volume molecules
``A``, ``B``, ``C``, ``D``, ``E``,  and ``F``.

.. math::

.. tabularcolumns:: |l|L|

+-----------------------------------+-----------------------------------------+
| **Example**                       | **Explanation**                         |
+===================================+=========================================+
| ``A + B + C -> D [1e12]``         | Volume molecules ``A``, ``B`` and ``C`` |
|                                   | react to yield product ``D``, at a rate |
|                                   | of :math:`10^{12}M^{-2}s^{-1}`.         |
+-----------------------------------+-----------------------------------------+
| ``A + B + C -> D + E + F [1e11]`` | Volume molecule ``A``, ``B`` and ``C``  |
|                                   | react to yield the three volume         |
|                                   | products ``D``, ``E`` and ``F`` at a    |
|                                   | rate of :math:`10^{11}M^{-2}s^{-1}`.    |
+-----------------------------------+-----------------------------------------+

The following table shows several examples involving a mixture of volume
molecules ``A``, ``B``, ``C`` , ``D`` and surface molecules ``S``, ``R``,
``T``, and ``U``

.. math::

.. tabularcolumns:: |l|L|

+-----------------------------------------+-----------------------------------+
| **Example**                             | **Explanation**                   |
+=========================================+===================================+
| ``A' + B' + S, -> D' [1e12]``           | Volume molecules ``A`` and ``B``  |
|                                         | both react with the bottom of     |
|                                         | surface molecule ``S`` to yield   |
|                                         | volume product ``D`` which is     |
|                                         | released toward the same side     |
|                                         | from which ``A`` and ``B`` came   |
|                                         | from at a rate of                 |
|                                         | :math:`10^{12}M^{-2}s^{-1}`.      |
+-----------------------------------------+-----------------------------------+
| ``A, + B, + S' -> D, [1e12]``           | This reaction is identical to the |
|                                         | previous one.                     |
+-----------------------------------------+-----------------------------------+
| ``A, + B, + S' -> A' + B' + S' [1e9]``  | This reaction describes the       |
|                                         | action of a surface bound         |
|                                         | symporter molecule ``S``.         |
|                                         | Molecules ``A`` and ``B`` bind to |
|                                         | the bottom of ``S`` which then    |
|                                         | re-releases ``A`` and ``B`` at    |
|                                         | its top domain. This reaction     |
|                                         | happens with a rate of            |
|                                         | :math:`10^9M^{-2}s^{-1}`.         |
+-----------------------------------------+-----------------------------------+
| ``A, + B' + S' -> A' + B, + S' [1e9]``  | This is similar to the previous   |
|                                         | reaction but ``S`` now acts as an |
|                                         | antiporter for ``A`` and ``B``.   |
+-----------------------------------------+-----------------------------------+
| ``A, + S' + R'' -> T'' [1e11]``         | In this reaction, volume molecule |
|                                         | ``A`` facilitates the             |
|                                         | dimerization of surface molecules |
|                                         | ``S`` and ``R``. ``A`` reacts     |
|                                         | with the bottom of ``S`` and      |
|                                         | ``R`` in arbitrary orientation to |
|                                         | produce a dimer ``T`` that is     |
|                                         | oriented like ``R``. The reaction |
|                                         | happens with a rate of            |
|                                         | :math:`10^{11}                    |
|                                         | {\mu}m^2N^{-1}M^{-1}s^{-1}`.      |
+-----------------------------------------+-----------------------------------+
| ``R, + S, + T'' -> T'' + U,,, [1e11]``  | Identically oriented surface      |
|                                         | molecules ``R`` and ``S``         |
|                                         | dimerize in the presence of       |
|                                         | surface molecule ``T`` which is   |
|                                         | oriented opposite to both ``R``   |
|                                         | and ``S``. The reaction           |
|                                         | regenerates ``T`` in its original |
|                                         | orientation and creates the dimer |
|                                         | ``U`` which can have an arbitrary |
|                                         | orientation. This reaction occurs |
|                                         | at a rate of :math:`10^{11}       |
|                                         | {\mu}m^4N^{-2}s^{-1}`.            |
+-----------------------------------------+-----------------------------------+

The units for the rates of trimolecular reactions depend on the reaction type
and are as below, where M is the molarity of the solution and N the number of
reactants.

-  [:math:`M^{-2}s^{-1}`] for trimolecular reactions between
   either three volume molecules or two volume molecule and a surface
   molecule,
-  [:math:`{\mu}m^2N^{-1}M^{-1}s^{-1}`] for trimolecular reactions between one
   volume molecule and two surface molecules, and
-  [:math:`{\mu}m^4N^{-2}s^{-1}` ] for trimolecular reactions
   involving three surface molecules.

.. _geom_def_commands:

Geometry definition commands
--------------------------------

.. _surf_props:

Surface properties
~~~~~~~~~~~~~~~~~~~~~~~~

MCell3 allows the user to specify properties of the surfaces of objects. For
example, one may wish to specify that a surface does not block the diffusion of
molecules. Each type of surface is defined by name, and each surface name must
be unique in the simulation and should not match any molecule names. Surface
properties are specified inside a surface definition block:

.. tabularcolumns:: |p{5cm}|p{11cm}|

+-------------------------------------+---------------------------------------+
| **Command**                         | **Explanation**                       |
+=====================================+=======================================+
|  | ``DEFINE_SURFACE_CLASS`` *name*  | Define a single surface type called   |
|  | ``{``                            | *name*. The properties are specified  |
|  | *  surface property commands*    | by zero or more commands inside       |
|  | ``}``                            | braces.                               |
+-------------------------------------+---------------------------------------+

.. tabularcolumns:: |p{7cm}|p{9cm}|

+------------------------------------------------------+----------------------+
| **Command**                                          | **Explanation**      |
+======================================================+======================+
|  | ``DEFINE_SURFACE_CLASSES``                        | Define a series of   |
|  | ``{``                                             | surface types by     |
|  | *  nameA* ``{`` *surface property commands* ``}`` | name.                |
|  | *  nameB* ``{`` *surface property commands* ``}`` |                      |
|  | *  ...*                                           |                      |
|  | ``}``                                             |                      |
+------------------------------------------------------+----------------------+

To define surface properties, use the following commands:

.. math:: 

.. tabularcolumns:: |p{6cm}|p{10cm}|

.. cssclass:: longtable

+--------------------------------------+---------------------------------------------+
| **Surface Property Command**         | **Explanation**                             |
+======================================+=============================================+
| ``REFLECTIVE =`` *name*              | If *name* refers to a volume molecule it is |
| :index:`\ <single:REFLECTIVE>`       | reflected by any surface with this surface  |
|                                      | class. This is the default behavior for     |
|                                      | volume molecules. If *name* refers to a     |
|                                      | surface molecule it is reflected by the     |
|                                      | border of the surface with this surface     |
|                                      | class. Tick marks on the *name* allow       |
|                                      | selective reflection of volume molecules    |
|                                      | from only the front or back of a surface or |
|                                      | selective reflection of surface molecules   |
|                                      | with only a certain orientation from the    |
|                                      | surface's border. Using the keyword         |
|                                      | ``ALL_MOLECULES`` for *name* has the effect |
|                                      | that all volume molecules are reflected by  |
|                                      | surfaces with this surface class and all    |
|                                      | surface molecules are reflected by the      |
|                                      | border of the surfaces with this surface    |
|                                      | class. Using the keyword                    |
|                                      | ``ALL_VOLUME_MOLECULES`` for the *name* has |
|                                      | the effect that all volume molecules are    |
|                                      | reflected by surfaces with this surface     |
|                                      | class. Using the keyword                    |
|                                      | ``ALL_SURFACE_MOLECULES`` has the effect    |
|                                      | that all surface molecules are reflected by |
|                                      | the border of the surface with this surface |
|                                      | class.                                      |
+--------------------------------------+---------------------------------------------+
| ``TRANSPARENT =`` *name*             | If *name* refers to a volume molecule it    |
| :index:`\ <single:TRANSPARENT>`      | passes through all surfaces with this       |
|                                      | surface class. If *name* refers to a        |
|                                      | surface molecule it passes through the      |
|                                      | border of the surface with this surface     |
|                                      | class. This is the default behavior for     |
|                                      | surface molecules. Tick marks on\ *name*    |
|                                      | allow the creation of one-way transparent   |
|                                      | surfaces for volume molecules or one-way    |
|                                      | transparent surface borders for surface     |
|                                      | molecules. To make a surface with this      |
|                                      | surface class transparent to all volume     |
|                                      | molecules, use ``ALL_VOLUME_MOLECULES`` for |
|                                      | *name*. To make a border of the surface     |
|                                      | with this surface class transparent to all  |
|                                      | surface molecules, use                      |
|                                      | ``ALL_SURFACE_MOLECULES`` for *name*. Using |
|                                      | the keyword ``ALL_MOLECULES`` for *name*    |
|                                      | has the effect that surfaces with this      |
|                                      | surface class are transparent to all volume |
|                                      | molecules and borders of the surfaces with  |
|                                      | this surface class are transparent to all   |
|                                      | surface molecules.                          |
+--------------------------------------+---------------------------------------------+
| ``ABSORPTIVE =`` *name*              | If *name* refers to a volume molecule it is |
| :index:`\ <single:ABSORPTIVE>`       | destroyed if it touches surfaces with this  |
|                                      | surface class. If *name* refers to a        |
|                                      | surface molecule it is destroyed if it      |
|                                      | touches the border of the surface with this |
|                                      | surface class. Tick marks on *name* allow   |
|                                      | destruction from only one side of the       |
|                                      | surface for volume molecules or selective   |
|                                      | destruction for surface molecules on the    |
|                                      | surfaces's border based on their            |
|                                      | orientation. To make a surface with this    |
|                                      | surface class absorptive to all volume      |
|                                      | molecules, ``ALL_VOLUME_MOLECULES`` can be  |
|                                      | used for *name*. To make a border of the    |
|                                      | surface with this surface class absorptive  |
|                                      | to all surface molecules,                   |
|                                      | ``ALL_SURFACE_MOLECULES`` can be used for   |
|                                      | *name*. Using the keyword ``ALL_MOLECULES`` |
|                                      | has the effect that surfaces with this      |
|                                      | surface class are absorptive for all volume |
|                                      | molecules and borders of the surfaces with  |
|                                      | this surface class are absorptive for all   |
|                                      | surface molecules.                          |
+--------------------------------------+---------------------------------------------+
| ``CLAMP_CONCENTRATION``              | The molecule called *name* is destroyed if  |
| *name* ``=`` *value*                 | it touches the surface (as if it had passed |
| :index:`\ <single:CLAMP_CONC>`       | through), and new molecules are created at  |
|                                      | the surface, as if molecules had passed     |
|                                      | through from the other side at a            |
|                                      | concentration *value* (units = M).          |
|                                      | Orientation marks may be used; in this      |
|                                      | case, the other side of the surface is      |
|                                      | reflective. Note that this command is only  |
|                                      | used to set the effective concentration of  |
|                                      | a volume molecule at a surface; it is not   |
|                                      | valid to specify a surface molecule. This   |
|                                      | command can be abbreviated as               |
|                                      | ``CLAMP_CONC``.                             |
+--------------------------------------+---------------------------------------------+
| | ``MOLECULE_DENSITY``               | Add the named molecules at the specified    |
| | ``{``                              | densities *D1*, *D2*, *...*, (units =       |
| | *  name1* ``=`` *D1*               | :math:`{\mu}m^{-2}`) to every surface with  |
| | *  name2* ``=`` *D2*               | this surface class. Use orientation marks   |
| | ``}``                              | after the name to specify the direction     |
|                                      | relative to the surface normal. For example,|
|                                      | ``A'`` specifies a molecule in the same     |
|                                      | orientation as the surface, while ``A,``    |
|                                      | specifies the opposite orientation. Using   |
|                                      | both marks indicates that the molecule      |
|                                      | should be assigned an orientation randomly. |
|                                      |                                             |
|                                      |                                             |
+--------------------------------------+---------------------------------------------+
| | ``MOLECULE_NUMBER``                | Add the exact numbers *N1*, *N2*, *...*, of |
| | ``{``                              | molecules onto any region that is made out  |
| | *  name1* ``=`` *N1*               | of this surface class. Note: this usage is  |
| | *  name2* ``=`` *N2*               | not recommended; it is better to add exact  |
| | ``}``                              | numbers of molecules to the region.         |
|                                      | Orientation marks after the name must be    |
|                                      | used to specify the direction the molecules |
|                                      | are facing.                                 |
|                                      |                                             |
|                                      |                                             |
+--------------------------------------+---------------------------------------------+

Note that surface normals are defined by the right-hand rule applied to the
vertices in order as listed (see section :ref:`geom_objs`). Box objects are
converted internally into triangles and the surface normals point outwards.

.. _geom_objs:

Geometrical objects
~~~~~~~~~~~~~~~~~~~~~~~~~~

Two types of geometrical objects are supported in MCell3. Objects can not have
coincident surfaces. Geometrical objects can be defined using:

.. tabularcolumns:: |p{5cm}|p{11cm}|

+---------------------------------+-------------------------------------------+
| **Command**                     | **Explanation**                           |
+=================================+===========================================+
|   | *name* ``BOX``              | This defines a box object called *name*.  |
|   | ``{``                       | The shape and position of the box is      |
|   | *  box commands*            | defined by . Optionally, additional       |
|   | *  region commands*         | commands can create regions and perform   |
|   | *  transformation commands* | geometrical transformations on the box.   |
|   | ``}``                       | Internally, a box is represented as a set |
|                                 | of triangles.                             |
+---------------------------------+-------------------------------------------+
 
.. tabularcolumns:: |p{5cm}|p{11cm}|

+---------------------------------+-------------------------------------------+
| **Command**                     | **Explanation**                           |
+=================================+===========================================+
|   | *name* ``POLYGON_LIST``     | This defines a polygon list object called |
|   | ``{``                       | *name*. Polygon list objects explicitly   |
|   | *  polygon commands*        | give their triangular surface elements.   |
|   | *  region commands*         |                                           |
|   | *  transformation commands* |                                           |
|   | ``}``                       |                                           |
+---------------------------------+-------------------------------------------+

A variety of optional commands can be used inside a geometrical object
definition block, after corners or vertex list / element connections are
specified, to modify the basic composition of the object and its surface
properties. These are described below. Geometrical transformations are
described later, in section :ref:`geom_trans`.

.. tabularcolumns:: |p{7cm}|p{9cm}|

+-----------------------------------+------------------------------------------+
| **Box Command**                   | **Explanation**                          |
+===================================+==========================================+
| ``CORNERS = [`` *x1* ``,`` *y1*   | The box object has corners as specified. |
| ``,`` *z1* ``],[`` *x2* ``,``     | The first coordinates should be less than|
| *y2* ``,`` *z2* ``]``             | the second set of coordinates, although  |
|                                   | MCell3 may fix it if you do it           |
|                                   | incorrectly.                             |
+-----------------------------------+------------------------------------------+
| ``ASPECT_RATIO =`` *a*            | Make sure that the ratio of the long to  |
|                                   | short side of each triangle making up the|
|                                   | box is no more than *a*. The smallest    |
|                                   | allowed value is 2. The default is to not|
|                                   | care about triangle shape.               |
+-----------------------------------+------------------------------------------+

.. tabularcolumns:: |p{5cm}|p{11cm}|

+--------------------------------------------+---------------------------------------------------------------------------+
| **Polygon Command**                        | **Explanation**                                                           |
+============================================+===========================================================================+
| | ``VERTEX_LIST``                          | Specify the vertices of the triangles inside a polygon list object        |
| | ``{``                                    | inside braces. Each vertex is given by its triple                         |
| | ``  [`` *x0* ``,`` *y0* ``,`` *z0* ``]`` | ``[`` *x* ``,`` *y* ``,`` *z* ``]``. This command must be given           |
| | ``  [`` *x1* ``,`` *y1* ``,`` *z1* ``]`` | before the ``ELEMENT_CONNECTIONS`` command.                               |
| | ``  `` *...*                             |                                                                           |
| | ``}``                                    |                                                                           |
+--------------------------------------------+---------------------------------------------------------------------------+
| | ``ELEMENT_CONNECTIONS``                  | Specify the triangles by vertex indices. The vertices are numbered from   |
| | ``{``                                    | ``0`` upwards in the order they were given in the vertex list. The        |
| | ``  [`` *a0* ``,`` *b0* ``,`` *c0* ``]`` | direction of the surface normal is determined by the right-hand rule      |
| | ``  [`` *a1* ``,`` *b1* ``,`` *c1* ``]`` | while following the vertices. Each triangle is given by a triple          |
| | ``  `` *...*                             | ``[`` *a* ``,`` *b* ``,`` *c* ``]`` of vertex numbers. This               |
| | ``}``                                    | command must be given after the ``VERTEX_LIST`` command.                  |
+--------------------------------------------+---------------------------------------------------------------------------+

.. tabularcolumns:: |p{5cm}|p{11cm}|

+-------------------------------------+---------------------------------------------------------------------------+
| **Region Command**                  | **Explanation**                                                           |
+=====================================+===========================================================================+
| | ``DEFINE_SURFACE_REGIONS``        | Define regions on the object. The extent of a region is given by the      |
| | ``{``                             | element specifier commands (at least one is required). Molecules can be   |
| |    *nameA*                        | added and surface properties can be set with the optional regional        |
| |    ``{``                          | surface commands. You can have an arbitrary number of regions on an       |
| |      *element specifier commands* | object, and they may overlap if you wish. Molecules added to overlapping  |
| |      *regional surface commands*  | regions accumulate. Triangles belonging to multiple regions inherit all   |
| |    ``}``                          | parent regions' surface properties. Users have to make sure that in case  |
| |    *name2* ``{`` *...* ``}``      | of overlapped regions their surface properties are compatible. Every      |
| |    *...*                          | ``BOX`` and ``POLYGON_LIST`` object has a pre-defined ``ALL`` region      |
| | ``}``                             | which consists of the entire object and has no special properties.        |
|                                     |                                                                           |
+-------------------------------------+---------------------------------------------------------------------------+
| | ``REMOVE_ELEMENTS``               | Remove the portion of the object specified by the element specifiers.     |
| | ``{``                             | You can think of this as a special type of region that defines the        |
| |    *element specifier commands*   | removed portions of the object. No real region exists on any part of the  |
| | ``}``                             | object that has been removed. You can use a list of element               |
|                                     | numbers/names instead of element specifiers if you wish, but you cannot   |
|                                     | mix a list of element numbers/names with the element specifier syntax.    |
|                                     | It is an error to remove all elements in an object or region.             |
|                                     |                                                                           |
+-------------------------------------+---------------------------------------------------------------------------+

.. tabularcolumns:: |p{8cm}|p{8cm}|

+-------------------------------------------------+-------------------------------------------------------------------------------------+
| **Element Specifier Command**                   | **Explanation**                                                                     |
+=================================================+=====================================================================================+
| ``INCLUDE_ELEMENTS = [`` *list* ``]``           | Include the elements specified by number or name. For polygon objects, these refer  |
|                                                 | to the triangles defined by the element connections, counting from zero upwards in  |
|                                                 | the order given. For boxes, the side names ``LEFT``, ``RIGHT``, ``FRONT``,          |
|                                                 | ``BACK``, ``BOTTOM``, and ``TOP`` can be used to refer to the sides, where          |
|                                                 | left/right corresponds to the x axis (left is lower x values), front/back to y, and |
|                                                 | bottom/top to z. ``ALL_ELEMENTS`` refers to the entire object. Numbers can be       |
|                                                 | specified individually (separated by commas) or in ranges with the format *N*       |
|                                                 | ``TO`` *M*. The two styles can be mixed (separated by commas).                      |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| ``EXCLUDE_ELEMENTS = [`` *list* ``]``           | Exclude the elements listed. If this is the first element specifier, assume that    |
|                                                 | all elements not listed are included. If not, subtract from the existing list.      |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| ``INCLUDE_REGION =`` *name*                     | Include the existing region on this object called *name* into this region, too.     |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| ``EXCLUDE_REGION =`` *name*                     | Exclude the existing region on this object called *name* from this new region.      |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| ``INCLUDE_PATCH=[`` *x1* ``,`` *y1* ``,`` *z1*  | This specifier is only valid on box objects, and the corners must define a          |
| ``],[`` *x2* ``,`` *y2* ``,`` *z2* ``]``        | rectangular patch that is on exactly one side of the box. The box will be divided   |
|                                                 | into triangles in such a way that this patch consists of separate triangles and     |
|                                                 | will form a region.                                                                 |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| ``EXCLUDE_PATCH=[`` *x1* ``,`` *y1* ``,`` *z1*  | Exclude the patch from this region.                                                 |
| ``],[`` *x2* ``,`` *y2* ``,`` *z2* ``]``        |                                                                                     |
+-------------------------------------------------+-------------------------------------------------------------------------------------+

Multiple element specifier commands can be used within the same region
definition statement. When combining multiple commands the resulting elements
list may depend on the order of these keywords. After element specifiers,
regions can specify a surface type and add extra molecules using:

.. tabularcolumns:: |p{5cm}|p{11cm}|

+------------------------------------+-------------------------------------------------------+
| **Regional Surface Command**       | **Explanation**                                       |
+====================================+=======================================================+
| | ``SURFACE_CLASS =`` *name*       | Set the surface type of this region to the previously |
|                                    | defined surface class called *name*.                  |
+------------------------------------+-------------------------------------------------------+
| ``MOLECULE_DENSITY {`` *...* ``}`` | This is the same as the Surface Property Command of   |
|                                    | the same name.                                        |
+------------------------------------+-------------------------------------------------------+
| ``MOLECULE_NUMBER {`` *...* ``}``  | This is the same as the Surface Property Command of   |
|                                    | the same name. Its usage is recommended here, as a    |
|                                    | regional surface command, rather than as a surface    |
|                                    | property command, so that the number of molecules is  |
|                                    | specified in the same place as the geometry, thus     |
|                                    | making the density easier to figure out.              |
+------------------------------------+-------------------------------------------------------+

.. _rel_objs:

Release objects
~~~~~~~~~~~~~~~~~~~~~

Release objects place molecules into the world. Release objects provide the
only means of placing molecules in a three dimensional space, but some release
shapes can place molecules on surfaces as well. Release objects are defined
using the following commands:

.. tabularcolumns:: |p{5cm}|p{11cm}|

+-------------------------------------------------+---------------------------------------------------------------------------+
| **Command**                                     | **Explanation**                                                           |
+=================================================+===========================================================================+
| | *name* ``RELEASE_SITE``                       | Create a release site called *name*. The shape and method of release is   |
| | ``{``                                         | specified by the release site commands. Optionally, geometrical           |
| |    *release site commands*                    | transformations can be applied also.                                      |                
| |    *transformation commands*                  |                                                                           |
| | ``}``                                         |                                                                           |
+-------------------------------------------------+---------------------------------------------------------------------------+
| *name* ``CUBIC_RELEASE_SITE {`` *...* ``}``     | Create a cubic release site called *name*. Molecules are released in a    |
|                                                 | box as specified by the radius. (This is the same as using the            |
|                                                 | ``SHAPE=CUBIC`` command inside ``RELEASE_SITE``.)                         |
+-------------------------------------------------+---------------------------------------------------------------------------+
| *name* ``SPHERICAL_RELEASE_SITE {`` *...* ``}`` | Create a spherical release site called *name*. Molecules are released     |
|                                                 | uniformly within the sphere depending on the defined radius of the        |
|                                                 | object. (This is the same as using the ``SHAPE=SPHERICAL`` command        |
|                                                 | inside ``RELEASE_SITE``.)                                                 |
|                                                 |                                                                           |
+-------------------------------------------------+---------------------------------------------------------------------------+
| *name* ``SPHERICAL_SHELL_SITE {``\ *...*\ ``}`` | Create a spherical shell release site called *name*. Molecules are        |
|                                                 | distributed on a spherical shell at the defined radius of the object.     |
|                                                 | For now, you must specify the number to distribute, not a concentration.  |
|                                                 | (This is the same as using the ``SHAPE=SPHERICAL_SHELL`` command inside   |
|                                                 | ``RELEASE_SITE``.)                                                        |
+-------------------------------------------------+---------------------------------------------------------------------------+
| | ``DEFINE_RELEASE_PATTERN`` *name*             | Define a new release pattern according to the commands given. A release   |
| | ``{``                                         | pattern must be defined for anything other than release at the beginning  |
| |   *release pattern commands*                  | of the simulation. Release patterns must be defined before they are       |
| | ``}``                                         | used. Multiple release sites can use the same pattern.                    |
|                                                 |                                                                           |
+-------------------------------------------------+---------------------------------------------------------------------------+

The following commands define where, what, and when a release object releases
molecules:

.. math::

.. tabularcolumns:: |p{6cm}|p{10cm}|

.. cssclass:: longtable

+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| **Release Site Command**                                  | **Explanation**                                                                 |
+===========================================================+=================================================================================+
| ``SHAPE =`` *geometry*                                    | Release molecules in the specified shape. Valid shapes are ``CUBIC``,           |
| :index:`\ <single:SHAPE>`                                 | ``SPHERICAL``, ``SPHERICAL_SHELL``, and ``LIST``; or the name of region(s) on   |
|                                                           | which to release. Each region must already be instantiated or be inside the     |
|                                                           | same ``OBJECT`` as the release site (see ``OBJECT`` command). Region names can  |
|                                                           | be combined with ``+`` to indicate release on both regions, ``-`` to indicate   |
|                                                           | the release occurs on the first and not the second, and ``*`` to indicate the   |
|                                                           | release occurs only where the two regions overlap. Parentheses may be used for  |
|                                                           | grouping. Volume molecules will be released in the volume bounded by the        |
|                                                           | regions (each region must be closed); surface molecules will be released on the |
|                                                           | surface (and regions need not be closed). If the region name is omitted and     |
|                                                           | only the name of a ``BOX`` or ``POLYGON_LIST`` object is specified, the         |
|                                                           | object's ``ALL`` region will be used.                                           |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| ``LOCATION = [`` *x* ``,`` *y* ``,`` *z* ``]``            | The release occurs centered at this location. Only used for geometrical shapes. |
| :index:`\ <single:LOCATION>`                              |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| ``MOLECULE =`` *name*                                     | The named molecule is the one that will be released. Not used for the ``LIST``  |
|                                                           | shape. You must specify an orientation if the molecule is a surface molecule.   |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| | ``MOLECULE_POSITIONS``                                  | The named molecules are added in the locations given. The molecule names        |
| | ``{``                                                   | must be followed by orientation marks if they have a 2D diffusion               |
| |    *name1* ``[``\ *x1*\ ``,``\ *y1*\ ``,``\ *z1*\ ``]`` | constant. If a molecule has a 2D diffusion constant, it will be placed          |
| |    *name2* ``[``\ *x2*\ ``,``\ *y2*\ ``,``\ *z2*\ ``]`` | on the surface closest to the coordinate given. This command is used for        |
| |    *...*                                                | the ``LIST`` shape only.                                                        |
| | ``}``                                                   |                                                                                 |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| | ``SITE_DIAMETER =`` *d*                                 | For a geometrical release site, this releases molecules uniformly within        |
| | ``SITE_RADIUS =`` *r*                                   | a diameter *d* or a radius *r*. Not used for releases on regions. With          |
|                                                           | the ``LIST`` shape, this is the distance that surface molecules search          |
|                                                           | for a surface before giving up; free molecules pay no attention to this         |
|                                                           | value for the ``LIST`` shape.                                                   |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| | ``SITE_DIAMETER = [`` *x* ``,`` *y* ``,`` *z* ``]``     | Release is asymmetric with a different diameters in different                   |
| | ``SITE_RADIUS = [`` *x* ``,`` *y* ``,`` *z* ``]``       | directions, as indicated by the vector. Not used for releases on regions        |
|                                                           | or with the ``LIST`` shape.                                                     |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| ``RELEASE_PROBABILITY =`` *p*                             | This release does not occur every time, but rather with probability *p*.        |
| :index:`\ <single:RELEASE_PROBABILITY>`                   | (If omitted, the default is to release without fail.) Either the whole          |
|                                                           | release occurs or none of it does; the probability does not apply               |
|                                                           | molecule-by-molecule. *p* must be in the interval [0, 1].                       |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| ``NUMBER_TO_RELEASE =`` *n*                               | Release *n* molecules. For releases on regions, *n* can be negative, and        |
| :index:`\ <single:NUMBER_TO_RELEASE>`                     | the release will then remove molecules of that type from the region. To         |
|                                                           | remove all molecules of a type, just make *n* large and negative. It is         |
|                                                           | unwise to both add and remove molecules on the same timestep---the order        |
|                                                           | of addition and removal is not defined in that case. This directive is          |
|                                                           | not used for the ``LIST`` shape, as every molecule is specified.                |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| | ``CONCENTRATION =`` *c*                                 | Release molecules at concentration *c* molar for volumes and *d*                |
| | ``DENSITY =`` *d*                                       | molecules per square micron for surfaces. Neither can be used for the           |
|                                                           | ``LIST`` shape; ``DENSITY`` is only valid for regions.                          |
|                                                           |                                                                                 |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| | ``GAUSSIAN_RELEASE_NUMBER``                             | Release molecules according to a Gaussian distribution with the                 |
| | ``{``                                                   | specified mean and standard deviation.                                          |
| |   ``MEAN_NUMBER =`` *n*                                 |                                                                                 |
| |   ``STANDARD_DEVIATION =`` *s*                          |                                                                                 |
| | ``}``                                                   |                                                                                 |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+
| ``RELEASE_PATTERN =`` *name*                              | Use the named release pattern instead of the default. The default is to         |
| :index:`\ <single:RELEASE_PATTERN>`                       | release the specified number of molecules at the beginning of the               |
|                                                           | simulation. If *name* is the name of a reaction pathway, the release            |
|                                                           | event will happen every time that reaction happens. The location will           |
|                                                           | then be relative to the site of the reaction, and the z-axis will be            |
|                                                           | rotated to align with the surface normal if the reaction was at a               |
|                                                           | surface. This is much slower than creating products within a reaction,          |
|                                                           | so only use it for special cases (e.g. synaptic vesicle release with a          |
|                                                           | random or very large number of neurotransmitter molecules).                     |
|                                                           |                                                                                 |
+-----------------------------------------------------------+---------------------------------------------------------------------------------+

Release patterns are defined as follows.

.. math::

.. tabularcolumns:: |p{6cm}|p{10cm}|

.. cssclass:: longtable

+------------------------------------------+------------------------------------------+
| **Release Pattern Command**              | **Explanation**                          |
+==========================================+==========================================+
| ``DELAY =`` *t*                          | The release pattern will start at time   |
| :index:`\ <single:DELAY>`                | *t*. (Default is to start at time zero.) |
+------------------------------------------+------------------------------------------+
| ``RELEASE_INTERVAL =`` *t*               | During a train of releases, release      |
| :index:`\ <single:RELEASE_INTERVAL>`     | molecules after every *t* seconds.       |
|                                          | Default is to release only once (*t* =   |
|                                          | :math:`{\infty}`).                       |
+------------------------------------------+------------------------------------------+
| ``TRAIN_DURATION =`` *t*                 | The train of releases lasts for *t*      |
| :index:`\ <single:TRAIN_DURATION>`       | seconds before turning off. Default is   |
|                                          | to never turn off (*t* =                 |
|                                          | :math:`{\infty}`).                       |
+------------------------------------------+------------------------------------------+
| ``TRAIN_INTERVAL =`` *t*                 | A new train of releases happens every    |
| :index:`\ <single:TRAIN_INTERVAL>`       | *t* seconds. Default is to never have a  |
|                                          | new train (*t* = :math:`{\infty}`). The  |
|                                          | train interval must not be shorter than  |
|                                          | the train duration.                      |
+------------------------------------------+------------------------------------------+
| ``NUMBER_OF_TRAINS =`` *n*               | Repeat the release process for *n*       |
| :index:`\ <single:NUMBER_OF_TRAINS>`     | trains of releases. Default is one       |
|                                          | train.                                   |
+------------------------------------------+------------------------------------------+
| ``NUMBER_OF_TRAINS = UNLIMITED``         | Repeat trains forever.                   |
| :index:`\ <single:NUMBER_OF_TRAINS>`     |                                          |
+------------------------------------------+------------------------------------------+

.. _inst_group_mod_objs:

Instantiation, grouping, and modification of objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An object is a box, polygon, release site, or a meta object which
contains other objects. Meta objects are defined and modified using

.. tabularcolumns:: |p{6cm}|p{10cm}|

+--------------------------------------------------+--------------------------------------------------------------------------+
| **Command**                                      | **Explanation**                                                          |
+==================================================+==========================================================================+
| | *name* ``OBJECT``                              | Define a new object called *name*. Inside the braces, list other objects |
| | ``{``                                          | one at a time to be added (see below).                                   |
| |    *object specifier commands*                 |                                                                          |
| |    *transformation commands*                   |                                                                          |
| | ``}``                                          |                                                                          |
|                                                  |                                                                          |
+--------------------------------------------------+--------------------------------------------------------------------------+
| ``INSTANTIATE`` *name* ``OBJECT {`` *...* ``}``  | Same as above, except we also insert the object into the world. A        |
|                                                  | simulation must have at least one ``INSTANTIATE``\ d object.             |
|                                                  |                                                                          |
|                                                  |                                                                          |
|                                                  |                                                                          |
|                                                  |                                                                          |
|                                                  |                                                                          |
+--------------------------------------------------+--------------------------------------------------------------------------+
| | ``MODIFY_SURFACE_REGIONS``                     | This modifies surface regions on existing objects via their name and     |
| | ``{``                                          | region name. Element lists may not be changed, but otherwise all         |
| |    *nameA* ``[`` *regA1* ``] {``               | regional surface commands are available. The full name must be given in  |
| |      *regional surface commands*               | the case of separate objects (using *name1.name2* to refer to            |
| |    ``}``                                       | objects inside meta objects). If an object is included in a meta object, |
| |    *nameB* ``[`` *regB1* ``] {`` *...* ``}``   | then has a surface region modified, and is included in another meta      |
| |    *  ...*                                     | object, the surface regions will differ in those the two meta objects.   |
| | ``}``                                          |                                                                          |
|                                                  |                                                                          |
+--------------------------------------------------+--------------------------------------------------------------------------+

You can define release sites, boxes, and polygon objects inside another
object, as well as placing previously defined objects into existing
ones:

.. tabularcolumns:: |p{6cm}|p{10cm}|

+-----------------------------------------+--------------------------------------------------------------------------------+
| **Object Specifier Command**            | **Explanation**                                                                |
+=========================================+================================================================================+
| | *newname* ``OBJECT`` *oldname*        | Add the existing object called *oldname* into the existing object and label it |
| | ``{``                                 | *newname*. You can add extra commands (e.g. transformation) inside the braces. |
| |    *transformation commands*          | The old and new names can be the same thing.  Thereafter, this object can be   |
| | ``}``                                 | referred to in the world as *name.newname*.                                    |
|                                         |                                                                                |
+-----------------------------------------+--------------------------------------------------------------------------------+
| *name* ``BOX {`` *...* ``}``            | Create a box inside the existing object (using the same syntax as              |
|                                         | previously defined).                                                           |
|                                         |                                                                                |
+-----------------------------------------+--------------------------------------------------------------------------------+
| *name* ``POLYGON_LIST {``\ *...*\ ``}`` | Create a polygon list object inside the existing object (using the same        |
|                                         | syntax as previously defined).                                                 |
|                                         |                                                                                |
+-----------------------------------------+--------------------------------------------------------------------------------+
| *name* ``RELEASE_SITE {``\ *...*\ ``}`` | Create a release site inside the existing object.                              |
|                                         |                                                                                |
+-----------------------------------------+--------------------------------------------------------------------------------+
| *newname* ``OBJECT {`` *...* ``}``      | Create an object inside the existing object.                                   |
|                                         |                                                                                |
+-----------------------------------------+--------------------------------------------------------------------------------+

.. _geom_trans:

Geometrical transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the end of the definition of a release object or geometrical object, or in
the block where an object is instantiated, it can be moved using the following
transformation commands (placed at the end of the block before the closing
brace).

.. math::

.. tabularcolumns:: |p{6cm}|p{10cm}|

+----------------------------------------------------+-------------------------------------------+
| **Transformation Command**                         | **Explanation**                           |
+====================================================+===========================================+
| ``TRANSLATE = [`` *x* ``,`` *y* ``,`` *z* ``]``    | Move the object by the specified vector.  |
| :index:`\ <single:TRANSLATE>`                      |                                           |
+----------------------------------------------------+-------------------------------------------+
| ``SCALE = [`` *x* ``,`` *y* ``,`` *z* ``]``        | Scale the object by multiplying each      |
| :index:`\ <single:SCALE>`                          | coordinate by the corresponding value in  |
|                                                    | the vector.                               |
+----------------------------------------------------+-------------------------------------------+
| ``ROTATE = [`` *x* ``,`` *y* ``,`` *z* ``] ,`` *A* | Rotate *A* degrees about the axis defined |
| :index:`\ <single:ROTATE>`                         | by the supplied vector.                   |
+----------------------------------------------------+-------------------------------------------+

.. _output_spec_commands:

Output specification commands
---------------------------------

There are two forms of output in MCell3, visualization output and count output.
Visualization output contains the molecules of the model in a form suitable for
visualization or analysis that requires knowledge of the precise location of
particles. Count output reports running totals of summary statistics such as
the total number of molecules of a certain type in the world, the number of
times a reaction has occurred inside some object in the world, and so on. Count
output can also be written when triggered by a specific event such as a
reaction taking place.

.. _viz_output:

Visualization Output
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tabularcolumns:: |p{6cm}|p{10cm}|

+----------------------------+------------------------------------------------+
| **Command**                | **Explanation**                                |
+============================+================================================+
| | ``VIZ_OUTPUT``           | Define a new visualization output block. MDL   |
| | ``{``                    | files can have multiple ``VIZ_OUTPUT`` blocks. |
| |    *viz output commands* |                                                |
| | ``}``                    |                                                |
+----------------------------+------------------------------------------------+

Each viz output block consists of the following commands:

.. math::

.. tabularcolumns:: |p{6cm}|p{10cm}|

+-----------------------------+-----------------------------------------------+
| **Viz Output Command**      | **Explanation**                               |
+=============================+===============================================+
| ``MODE`` = *viz_mode*       | Specifies the mode of the visualization       |
| :index:`\ <single:MODE>`    | output. The valid values are                  |
|                             | ``CELLBLENDER`` , ``ASCII`` , and ``NONE``.   |
|                             | Most users will want to use                   |
|                             | ``CELLBLENDER`` mode. ``ASCII`` mode will     |
|                             | export the data in a human-readable format.   |
|                             | ``NONE`` mode can be used as a way to disable |
|                             | viz data (without using comments).            |
+-----------------------------+-----------------------------------------------+
| ``FILENAME = "``            | Directory and filename prefix for all of the  |
| *filename_specifier* ``"``  | binary or ASCII data files.                   |
+-----------------------------+-----------------------------------------------+
| | ``MOLECULES``             | Defines molecules visualization data output   |
| | ``{``                     | block.                                        |
| |    *data output block*    |                                               |
| | ``}``                     |                                               |
|                             |                                               |
+-----------------------------+-----------------------------------------------+

Each data output block consists of the following commands:

.. tabularcolumns:: |p{6cm}|p{10cm}|

+--------------------------------------------------+--------------------------------------------------------------------------+
| **Data Output Block Command**                    | **Explanation**                                                          |
+==================================================+==========================================================================+
| | ``NAME_LIST``                                  | Defines a valid name list. The valid values are either names separated   |
| | ``{``                                          | by any type of whitespace, strings with wildcards (in quotes) that match |
| |    *name list commands*                        | names, or keywords defined below. All children of the named objects are  |
| | ``}``                                          | included by default. If this statement occurs in a ``MESHES`` block, the |
|                                                  | names should be names of objects; in a ``MOLECULES`` block they should   |
|                                                  | be names of molecules.                                                   |
|                                                  |                                                                          |
+--------------------------------------------------+--------------------------------------------------------------------------+
| | ``TIME_POINTS``                                | Defines what data should be output at what times. The data types are     |
| | ``{``                                          | given below and valid notations for *time_points_list* are               |
| |    *data type* ``@`` *time_points_list*        | ``[`` *time1* ``]``, or ``[`` *time1* ``,`` *time2* ``,`` … ``,``        |
| | ``}``                                          | *time_end* ``]``, or ``[`` *time1* ``,`` *time2* ``, [`` *time3*         |
|                                                  | ``TO`` *time_end* ``STEP`` *delta_time*\ ``]]``, or ``ALL_TIMES.``       |
|                                                  | Mutually exclusive with ``ITERATION_NUMBERS``.                           |
|                                                  |                                                                          |
+--------------------------------------------------+--------------------------------------------------------------------------+
| | ``ITERATION_NUMBERS``                          | Defines what data should be output at what iterations. The data types    |
| | ``{``                                          | are given below and valid notations for *iteration_numbers_list* are     |
| |    *data type* ``@`` *iterations_numbers_list* | ``[`` *iteration1* ``]``, or ``[`` *iteration1* ``,``                    |
| | ``}``                                          | *iteration2* ``,`` … ``,`` *iteration_end* ``]``, or                     |
|                                                  | ``[`` *iteration1* ``,`` *iteration2* ``, [`` *iteration3* ``TO``        |
|                                                  | *iteration_end* ``STEP`` *delta_iteration* ``]]``, or                    |
|                                                  | ``ALL_ITERATIONS.`` Mutually exclusive with ``TIME_POINTS``.             |
|                                                  |                                                                          |
+--------------------------------------------------+--------------------------------------------------------------------------+

The following name list commands for ``MOLECULES`` are available:

+----------------------------------------+------------------------------------+
| **Name list Commands** (``MOLECULES``) | **Explanation**                    |
+========================================+====================================+
| ``ALL_MOLECULES``                      | All molecule names should be       |
|                                        | included in the ``NAME_LIST``      |
|                                        | sub-block inside ``MOLECULES``     |
|                                        | block.                             |
+----------------------------------------+------------------------------------+

The following data type commands for ``MOLECULES`` are available:

+---------------------------------+-------------------------------------------+
| **Data types** (``MOLECULES``)  | **Explanation**                           |
+=================================+===========================================+
| ``POSITIONS`` or ``ALL_DATA``   | Molecule position information should be   |
|                                 | written at the specified time/iteration.  |
+---------------------------------+-------------------------------------------+

All of the keywords in the ``VIZ_OUTPUT`` block are optional except
``FILENAME.``

Examples of ``VIZ_OUTPUT`` statements are given below.

**Option #1 (time style):**

.. code-block:: none

     VIZ_OUTPUT {
      FILENAME = "viz_data/output_example"
      MOLECULES {
        NAME_LIST { ALL_MOLECULES /* or list of molecule names */ }
        TIME_POINTS { ALL_DATA @ ALL_TIMES }
      }
    }

**Option #2 (iterations style):**

.. code-block:: none

     VIZ_OUTPUT {
      FILENAME = "viz_data/output_example"
      MOLECULES {
        NAME_LIST { ALL_MOLECULES /* or list of molecule names */ }
        ITERATION_NUMBERS { ALL_DATA @ ALL_ITERATIONS }
      }
    }

Usual UNIX-style wildcards like "\*" and "?" are allowed in the
*name_list* but must be enclosed in quotes. For example in the case of
``MOLECULES`` the following ``NAME_LIST`` statements are all valid:

.. code-block:: none

    NAME_LIST{A B C1 C2 C3} 
    NAME_LIST{A B "C*"} 
    NAME_LIST{A B "C?"}

.. _rxn_data_output:

Reaction Data Output
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tabularcolumns:: |p{6cm}|p{10cm}|

+---------------------------------+-------------------------------------------+
| **Command**                     | **Explanation**                           |
+=================================+===========================================+
| | ``REACTION_DATA_OUTPUT``      | Define a new count data output block which|
| | ``{``                         | contains the commands below.  Each MDL    |
| |    *reaction output commands* | file can have multiple reaction data      |
| | ``}``                         | output blocks.                            |
+---------------------------------+-------------------------------------------+

Each reaction data output block consists of the following commands:

.. tabularcolumns:: |p{6cm}|p{10cm}|

.. cssclass:: longtable

+------------------------------------------------------+-------------------------------------------------------------------------+
| **Reaction Output Command**                          | **Explanation**                                                         |
+======================================================+=========================================================================+
| | ``OUTPUT_BUFFER_SIZE =`` *N*                       | Write output to disk after every *N* lines. The default is *N* =10000.  |
|                                                      | This command is optional, but must be first if it is used. The output   |
|                                                      | will also always be written when the simulation terminates, regardless  |
|                                                      | of *N*.                                                                 |
+------------------------------------------------------+-------------------------------------------------------------------------+
| ``STEP =`` *t*                                       | Output this block every *t* seconds. Exactly one of ``STEP`` or the     |
|                                                      | following two commands should be used. Triggered output ignores the     |
|                                                      | values specified, but some value must still be given.                   |
+------------------------------------------------------+-------------------------------------------------------------------------+
| ``TIME_LIST = [`` *list* ``]``                       | Output this block at the times specified in the list.                   |
+------------------------------------------------------+-------------------------------------------------------------------------+
| ``ITERATION_LIST = [`` *list* ``]``                  | Output this block at the iteration numbers specified in the list (i.e.  |
|                                                      | after that number of timesteps).                                        |
+------------------------------------------------------+-------------------------------------------------------------------------+
| ``HEADER =`` *setting*                               | Output blocks by default have no header but can optionally have a       |
|                                                      | header line that states the output (name of molecule, reaction, etc.)   |
|                                                      | in each column. This command can set the behavior of that header line;  |
|                                                      | it applies to all output files until the next ``HEADER`` line. A        |
|                                                      | *setting* of ``ON`` turns on the header line; ``OFF`` prevents any      |
|                                                      | header. A string, in quotes, will turn the header on and prepend the    |
|                                                      | string to the line; this is useful to add comment character(s). For     |
|                                                      | example, ``''//''`` would add a C++-style comment prefix to the line.   |
|                                                      | For ``TRIGGER`` statements (see below), the column label (plus comment  |
|                                                      | character if specified) is appended to each line of output when headers |
|                                                      | are on.                                                                 |
+------------------------------------------------------+-------------------------------------------------------------------------+
| ``SHOW_EXACT_TIME =`` *setting*                      | ``TRIGGER`` statements (see below) can report timing information more   |
|                                                      | precisely than by iteration. However, if only iteration timing is of    |
|                                                      | interest, this can be set ``OFF``. The default is ``ON``. It applies to |
|                                                      | all output files until the next ``SHOW_EXACT_TIME`` line.               |
+------------------------------------------------------+-------------------------------------------------------------------------+
| ``{`` *value* ``} => "`` *file* ``"``                | Output the value in braces to the filename in quotes. The first column  |
|                                                      | will be the time (in seconds) of the iteration unless the               |
|                                                      | ``ITERATION_LIST`` specifier is used, in which case the first column    |
|                                                      | will be the iteration number. For ``COUNT`` values, the second column   |
|                                                      | will be the value of the count; other possibilities appear later in     |
|                                                      | this document. This command, and the variants listed below, can be      |
|                                                      | repeated to send different output to many files. The output symbol      |
|                                                      | ``=>`` has several variants which are described below.                  |
+------------------------------------------------------+-------------------------------------------------------------------------+
| ``{`` *value* ``: "`` *name* ``" } => "`` *file*     | Output the value in braces with the column header string *name* to the  |
| ``"``                                                | filename *file*. Not valid if *value* is found using wildcards. Trigger |
|                                                      | outputs put this header in the rightmost column on each line; count     |
|                                                      | outputs put the name at the top of the appropriate column.              |
+------------------------------------------------------+-------------------------------------------------------------------------+
| ``{`` *value* ``,`` *value* ``,`` *...* ``} => "``   | For counts, output the list of values in braces, one to a column, in    |
| *file* ``"``                                         | the order listed. The first column will be the time/iteration number;   |
|                                                      | successive columns will be the values in the order listed. If headers   |
|                                                      | are on, each column header can be customized by specifying:             |
|                                                      | ``''``\ *name*\ ``''`` after the value. For triggers, all the specified |
|                                                      | events will be combined into one file.                                  |
+------------------------------------------------------+-------------------------------------------------------------------------+

The *value* specified in braces is either a ``TRIGGER`` statement, a ``COUNT``
statement, or a mathematical operation involving ``COUNT`` statements and
constants. Currently, MCell supports addition (+), subtraction (-),
multiplication (\*), and division (/) with the corresponding operators given in
parenthesis. Furthermore, expressions can be grouped using parenthesis. Hence,
the following is a valid *value* expression

.. code-block:: none

    { (COUNT[A,WORLD] + COUNT[B,WORLD]) * 3.0 }

Wildcards can be used to select multiple molecules or reactions by name, but in
this case mathematical operations cannot be used. The wildcards ``?`` and ``*``
can be used to match any single character and any sequence of characters,
respectively; internally, this will generate one count/trigger statement per
matching name. Having headers on is convenient in this case, so one can tell
which column (for ``COUNT`` statements) or row (for ``TRIGGER`` statements)
corresponds to which name.

If a simulation starts from a checkpoint file, it will add to any existing
output files. Otherwise, the output files will be overwritten if they already
exist.

``COUNT`` statements are either *instantaneous*, and give information about the
state of the model at the instant the count is output---the number of molecules
in a region, for example---or are *cumulative*, and count the number of events
that have occurred since the beginning of the simulation. Alternatively, they
can output the time and location of each reaction or molecular collision of the
type specified. In all cases, if a region or object is referred to, it should
be the fully qualified name starting with the name of the instantiated object.

The ``COUNT`` statements themselves have the following syntax:

.. math::

.. tabularcolumns:: |p{8cm}|p{8cm}|

.. cssclass:: longtable

+------------------------------------------------------------+---------------------------------------------------------------------------------------+
| **Count Statement**                                        | **Explanation**                                                                       |
+============================================================+=======================================================================================+
| ``COUNT[`` *name* ``, WORLD]``                             | Count molecules or reactions in the world. If *name* refers to a molecule, this is an |
|                                                            | instantaneous count of the number of copies that molecule in the world. If *name*     |
|                                                            | refers to a reaction, count how many times that reaction has occurred since the       |
|                                                            | beginning of the simulation. If ``''`` *name* ``''`` is in quotes, in this command or |
|                                                            | any of the following commands, the string in quotes can contain wildcards which will  |
|                                                            | be matched to molecule and reaction names and will be listed in alphabetical order.   |
|                                                            | It is usually a good idea when using wildcards to turn on headers so one can see      |
|                                                            | which column is which.                                                                |
+------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``COUNT[`` *name* ``,`` *object* ``]``                     | Count molecules or reactions inside the object called *object*. This must be an       |
|                                                            | instantiated object. For example, if you have instantiated an object called           |
|                                                            | ``my_world`` with a box called ``my_box`` inside it, *object* would be                |
|                                                            | ``my_world.my_box``. If you are counting surface molecules or reactions at a surface, |
|                                                            | only the ones that actually occur on *object* will be counted (not those inside which |
|                                                            | are on a different object). Molecules with a 3D diffusion constant will be counted    |
|                                                            | inside the object, but the object must be closed. All counts are instantaneous.       |
+------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``COUNT[`` *name* ``,`` *region* ``]``                     | Count molecules or reactions inside the named region. For a grid molecule, *name* can |
|                                                            | also specify its surface orientation and in such a case has to be enclosed in quotes, |
|                                                            | e.g., " ``A,`` ". The surface orientation may be given by an arbitrary number of      |
|                                                            | either ``'``, ``,`` or ``;`` *.* Mixing is not possible. Equivalently, the numerical  |
|                                                            | orientation specifiers ``{-1},{0},`` or ``{1}`` can be used. Clearly, the             |
|                                                            | specification ``A;`` or ``A{0}`` is equivalent to ``A`` since ``;`` and ``{0}`` both  |
|                                                            | specify no orientation. The named region must be referenced fully. E.g. if ``my_box`` |
|                                                            | (from above) has a region called ``my_region``, the name would be                     |
|                                                            | ``my_world.my_box[my_region]``. The count is instantaneous. As with the object        |
|                                                            | syntax, molecules and reactions on surfaces must be on the named region, while volume |
|                                                            | molecules and reactions must be inside.                                               |
+------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``COUNT[`` *name* ``,`` *region* ``, ALL_ENCLOSED]``       | Count all molecules or reactions that occur in the area enclosed by region (not       |
|                                                            | counting those that occur on the surface of the region). Imagine, for an example, two |
|                                                            | cubes "outer" and "inner" such that "inner" is completely inside "outer". This        |
|                                                            | statement written for "outer" cube will effectively count grid molecules *name* on    |
|                                                            | the surface of "inner" cube only. For a grid molecule, *name* can also specify its    |
|                                                            | surface orientation and in such a case has to be enclosed in quotes, e.g., " ``A,``   |
|                                                            | ". The surface orientation may be given by an arbitrary number of either ``'``, ``,`` |
|                                                            | or ``;`` *.* Mixing is not possible. Equivalently, the numerical orientation          |
|                                                            | specifiers ``{-1},{0},`` or ``{1}`` can be used. Clearly, the specification ``A;`` or |
|                                                            | ``A{0}`` is equivalent to ``A`` since ``;`` and ``{0}`` both specify no orientation.  |
|                                                            | This ``COUNT`` statement lets you count surface molecules contained on surfaces that  |
|                                                            | lie within a box, for example. This will work with object names as well as region     |
|                                                            | names, but the object or region must be closed. It is only useful for surface         |
|                                                            | molecules and reactions at surfaces; adding ``ALL_ENCLOSED`` is valid for volume      |
|                                                            | molecules and reactions, but ``ALL_ENCLOSED`` is the default behavior. The count is   |
|                                                            | instantaneous.                                                                        |
+------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``COUNT[`` *molecule* ``,`` *region* ``, ESTIMATE_CONC]``  | Currently this feature applies only to volume molecules. Estimate the concentration   |
|                                                            | of the volume molecule at that region, averaged since the beginning of the simulation |
|                                                            | (output has units of :math:`{\mu}m`). A single object can be used instead of a        |
|                                                            | region. The region/object does not need to be closed. To find the average             |
|                                                            | concentration during one count interval, let :math:`t_i` be the time of the *i* th    |
|                                                            | output, let :math:`t_j` be some earlier output, and let *c* (*t*) be the              |
|                                                            | concentration averaged up to time *t*. Then the average concentration between times   |
|                                                            | :math:`t_j` and :math:`t_i` is :math:`\bar{c}(t_j \rightarrow                         |
|                                                            | t_i)=\frac{t_i\bar{c}(t_i)-t_j\bar{c}(t_j)}{t_i-t_j}` . Note that this is the         |
|                                                            | concentration all around the surface, so if the molecule can only reach one side, the |
|                                                            | concentration on that side will be twice what is reported here. The command can be    |
|                                                            | given verbosely as ``ESTIMATE_CONCENTRATION``. The estimate is based on a cumulative  |
|                                                            | count.                                                                                |
+------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``COUNT[`` *molecule* ``,`` *region* ``,`` *hits* ``]``    | For a volume molecule output the number of times the named molecule has hit the named |
|                                                            | region (or object). For a surface molecule output the number of times the named       |
|                                                            | molecule hit the boundary of the named region. The *hits* specifier should be one of  |
|                                                            | ``FRONT_HITS``, ``BACK_HITS``, ``ALL_HITS``, ``FRONT_CROSSINGS``, ``BACK_CROSSINGS``, |
|                                                            | and ``ALL_CROSSINGS``. For a volume molecule the meaning of these specifiers is       |
|                                                            | obvious. For a surface molecule FRONT means inside out direction, and BACK means      |
|                                                            | outside in direction. The count is cumulative.                                        |
+------------------------------------------------------------+---------------------------------------------------------------------------------------+
| ``EXPRESSION[`` *expression* ``]``                         | Evaluate and output a mathematical expression. This can be mixed with ``COUNT``       |
|                                                            | statements but not with ``TRIGGER`` statements.                                       |
+------------------------------------------------------------+---------------------------------------------------------------------------------------+

Cumulative counts are reset when a simulation is started from a checkpoint.
This breaks ``ESTIMATE_CONC``, but the other cumulative counts can be recovered
by adding the last report before the checkpoint to the first one after the
checkpoint.

``TRIGGER`` statements output the time and location each time the number of
molecules changes or a reaction happens. Most ``COUNT`` statements have a
corresponding ``TRIGGER``, but ``TRIGGER`` statements are not compatible with
the ``WORLD`` or the ``ESTIMATE_CONC`` directives.  Within output statements
pointing to the same output file, there can only be ``TRIGGER`` commands, i.e.,
they cannot be mixed with ``COUNT`` or ``EXPRESSION`` statements.

``TRIGGER`` statements obey the following syntax:

.. tabularcolumns:: |p{8cm}|p{8cm}|

.. cssclass:: longtable

+--------------------------------+--------------------------------------------+
| **Trigger Statement**          | **Explanation**                            |
+================================+============================================+
| ``TRIGGER[`` *molecule* ``,``  | Generates output each time the number of   |
| *region* ``]``                 | molecules inside the specified region      |
|                                | changes. The output has the format         |
|                                | *iteration_time exact_time X Y Z           |
|                                | orientation number [name]* as described    |
|                                | below. The sixth column, *orientation*,    |
|                                | gives the molecule orientation, i.e., it   |
|                                | is 0 for volume molecules and +/-1 for     |
|                                | surface molecules according to their       |
|                                | orientation with respect to the surface    |
|                                | containing them. The seventh column,       |
|                                | *number*, can take on values of +/-1       |
|                                | depending on if the molecule was added or  |
|                                | removed, respectively, from the region. If |
|                                | ``HEADER`` is on, the eighth column lists  |
|                                | the molecule name.                         |
+--------------------------------+--------------------------------------------+
| ``TRIGGER[`` *reaction* ``,``  | Generates output each time the named       |
| *region* ``]``                 | reaction takes place inside the specified  |
|                                | closed region. The output has the format   |
|                                | *iteration_time exact_time X Y Z [name].*  |
|                                | The fields are described below. Note that  |
|                                | since reactions do not have an orientation |
|                                | and always occur one at a time the         |
|                                | *orientation* and *number* fields are      |
|                                | omitted. If ``HEADER`` is on, the sixth    |
|                                | column lists the reaction name.            |
+--------------------------------+--------------------------------------------+
| ``TRIGGER[`` *name* ``,``      | This is equivalent to specifying a list of |
| *object* ``]``                 | ``TRIGGER`` statements which consist of    |
|                                | all regions in that object.                |
+--------------------------------+--------------------------------------------+
| ``TRIGGER[`` *name* ``,``      | Generates output each time the named       |
| *region* ``,ALL_ENCLOSED]``    | reaction takes place, or number of named   |
|                                | molecules changes, inside the specified    |
|                                | closed region. As with ``COUNT``           |
|                                | statements, this is only useful for        |
|                                | surface reactions and molecules, and does  |
|                                | not include the surface of the named       |
|                                | region, only events wholly inside it. The  |
|                                | output has the format appropriate for      |
|                                | molecules or reactions.                    |
+--------------------------------+--------------------------------------------+
| ``TRIGGER[`` *molecule* ``,``  | For a volume molecule generates output     |
| *region* ``,`` *hits* ``]``    | each time the molecule hits or crosses the |
|                                | named region. For a surface molecule       |
|                                | generates output each time the molecule    |
|                                | hits or crosses the boundary of the named  |
|                                | region. The *hits* specifier should be one |
|                                | of ``FRONT_HITS``, ``BACK_HITS``,          |
|                                | ``ALL_HITS``, ``FRONT_CROSSINGS``,         |
|                                | ``BACK_CROSSINGS``, and ``ALL_CROSSINGS``. |
|                                | For a volume molecule the meaning of these |
|                                | specifiers is obvious. For a surface       |
|                                | molecule FRONT means inside out direction, |
|                                | and BACK means outside in direction. The   |
|                                | output has the format *iteration_time      |
|                                | exact_time X Y Z orientation [name].* The  |
|                                | *orientation* column can take on values of |
|                                | +/-1 depending on if the region (region    |
|                                | boundary) was hit or crossed from the      |
|                                | front or the back (inside out or outside   |
|                                | in), respectively; other columns are       |
|                                | described below. Note that the *number*    |
|                                | column is omitted. If ``HEADER`` is on,    |
|                                | the seventh column lists the molecule      |
|                                | name.                                      |
+--------------------------------+--------------------------------------------+

The output contains one row of data for each even that happened. The format of
the columns is:

+---------------------+-----------------+-------+-------+-------+-------------------+--------------+------------+
| *iteration_time*    | *exact_time*    | *X*   | *Y*   | *Z*   | *[orientation]*   | *[number]*   | *[name]*   |
+=====================+=================+=======+=======+=======+===================+==============+============+
|                     |                 |       |       |       |                   |              |            |
+---------------------+-----------------+-------+-------+-------+-------------------+--------------+------------+

*Iteration_time* is the time of the iteration during which the event happened,
or the iteration number if ``ITERATION_LIST`` was specified for the block.

*Exact_time* is the time at which the event was scheduled, between
*iteration_time* and the time of the next iteration. Since events within one
iteration are not ordered precisely, *exact_time* values will not always
increase. This column can be turned off by using the ``SHOW_EXACT_TIME=OFF``
directive inside the ``REACTION_DATA_OUTPUT`` block. These values are always
times, even if ``ITERATION_LIST`` is specified for the block.

*X*, *Y*, and *Z* are the coordinates at which the event took place.  Reactions
and hits always report their coordinates precisely. Volume molecules that
disappear at a surface will report their final position as slightly inside the
surface along their last trajectory (so that it is possible to tell which side
of the surface they were on); if they react with another volume molecule they
will report the position they reached when their interaction disk intersected
the target molecule, not the position of the target. Surface molecules diffuse
by hopping rather than ray-tracing, so when a surface molecule leaves a region
of interest, the position reported is the last position where the molecule was
located inside the region, not the boundary of the region where it crossed out
(and conversely, when entering, it's the first position where the molecule
stopped at the end of its time-step).

*Orientation* and *number* are only provided for certain types of triggers and
are described above.

*Name* is the name of the molecule or a user-defined string, and present in the
last column (6, 7, or 8 depending on which type of trigger is used) if headers
are on.

The following output symbols can be used in place of ``=>`` and give the
behaviors described below. All output symbols will create files if none exist.
No output symbols will create directories---if the files that are referred to
cannot be created as specified, MCell3 will quit with an error message. Output
may create empty files if the simulation ends without producing output (either
because of an error condition or because the simulation did not run long enough
to reach the time/iteration of any reaction data output).

+--------------------+--------------------------------------------------------+
| **Output Symbol**  | **Explanation**                                        |
+====================+========================================================+
| ``=>``             | If a checkpoint file is not used, overwrite the        |
|                    | existing file (with headers if requested). If a        |
|                    | checkpoint file is used, discard any of the output     |
|                    | file that appears to be a later time than the start of |
|                    | the current run, and append to the file from that      |
|                    | point. Headers are not written unless the file has to  |
|                    | be created or is empty to begin with. This command     |
|                    | generally does "what you expect"-after the simulation  |
|                    | has run, it will contain data from earlier in the      |
|                    | simulation that the current run, plus the data created |
|                    | in the current simulation. If you switch between       |
|                    | ``ITERATION_LIST`` and other output time specifiers,   |
|                    | this command won't know whether output is by time or   |
|                    | by iteration number, so don't use this command if you  |
|                    | switch from one to the other after checkpointing.      |
+--------------------+--------------------------------------------------------+
| ``>``              | Always overwrite the file, whether or not a checkpoint |
|                    | is used. If headers are requested, they will appear at |
|                    | the beginning of the file.                             |
+--------------------+--------------------------------------------------------+
| ``+>``             | Always create a new file, whether or not a checkpoint  |
|                    | is used. If a file of the given name already exists    |
|                    | and is not empty, MCell3 will print an error message   |
|                    | and exit. If headers are requested, they will appear   |
|                    | at the beginning of the file.                          |
+--------------------+--------------------------------------------------------+
| ``>>``             | Always append to an existing file without removing any |
|                    | previous data. Headers are only written if the file    |
|                    | starts out empty or has to be created.                 |
+--------------------+--------------------------------------------------------+
| ``>>>``            | Always append to an existing file without removing any |
|                    | previous data and if headers are requested, write them |
|                    | even into the middle of the file.                      |
+--------------------+--------------------------------------------------------+

.. _deprecated_viz_commands:

Deprecated Visualization Output Commands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend that new users use ``CELLBLENDER`` mode. However, this section is
for users who still want to use DReAMM visualization output.

Each viz output block consists of the following commands:

.. tabularcolumns:: |p{7cm}|p{9cm}|

+-----------------------------------------------+---------------------------------------------------------------------------------------+
| **Viz Output Command**                        | **Explanation**                                                                       |
+===============================================+=======================================================================================+
| ``MODE`` = *viz_mode*                         | Specifies the mode of the visualization output. The mode defines the directory        |
|                                               | structure and number of files comprising the visualization output. The valid values   |
|                                               | are ``DREAMM_V3`` and ``DREAMM_V3_GROUPED``.                                          |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| ``FILENAME = "`` *filename_specifier* ``"``   | Name of the master header file containing all information for DReAMM and references   |
|                                               | to the multiple binary data files.                                                    |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| ``VIZ_MOLECULE_FORMAT =`` *output_mode*       | Specifies if molecule positions are being output in binary or plain-text format.      |
|                                               | Valid choices are ``BINARY`` or ``ASCII.`` This command can only be used in           |
|                                               | ``DREAMM_V3`` mode. If this command is not explicitly specified it defaults to        |
|                                               | ``BINARY.``                                                                           |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| ``VIZ_MESH_FORMAT =`` *output_mode*           | Specifies if mesh positions are being output in binary or plain-text format. Valid    |
|                                               | choices are ``BINARY`` or ``ASCII.`` This command can only be used in ``DREAMM_V3``   |
|                                               | mode. If this command is not explicitly specified it defaults to ``BINARY``.          |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| | ``MESHES``                                  | Defines meshes visualization data output block.                                       |
| | ``{``                                       |                                                                                       |
| |    *data output block*                      |                                                                                       |
| | ``}``                                       |                                                                                       |
|                                               |                                                                                       |
+-----------------------------------------------+---------------------------------------------------------------------------------------+
| | ``MOLECULES``                               | Defines molecules visualization data output block.                                    |
| | ``{``                                       |                                                                                       |
| |    *data output block*                      |                                                                                       |
| | ``}``                                       |                                                                                       |
|                                               |                                                                                       |
+-----------------------------------------------+---------------------------------------------------------------------------------------+

The following name list commands for ``MESHES`` are available:

+-------------------------------------+---------------------------------------+
| **Name list Commands** (``MESHES``) | **Explanation**                       |
+=====================================+=======================================+
| ``ALL_MESHES``                      | All mesh object names should be       |
|                                     | included in the ``NAME_LIST`` sub-    |
|                                     | block inside a ``MESHES`` block.      |
|                                     | ``ALL_MESHES`` is equivalent to       |
|                                     | naming the top-level mesh object      |
|                                     | (assuming that only a single          |
|                                     | ``INSTANTIATE`` block is present).    |
+-------------------------------------+---------------------------------------+

The name list commands for ``MOLECULES`` are identical to those listed
previously.

The following data type commands for ``MESHES`` are available:

+-----------------------------+-----------------------------------------------+
| **Data types** (``MESHES``) | **Explanation**                               |
+=============================+===============================================+
| ``GEOMETRY``                | Mesh vertex and connectivity information      |
|                             | should be written at the specified            |
|                             | time/iteration.                               |
+-----------------------------+-----------------------------------------------+
| ``REGION_DATA``             | Mesh region information should be written at  |
|                             | the specified time/iteration.                 |
+-----------------------------+-----------------------------------------------+
| ``ALL_DATA``                | Equivalent to using both ``GEOMETRY`` and     |
|                             | ``REGION_DATA``.                              |
+-----------------------------+-----------------------------------------------+

+--------------------------------+--------------------------------------------+
| **Data types** (``MOLECULES``) | **Explanation**                            |
+================================+============================================+
| ``POSITIONS``                  | Molecule position information should be    |
|                                | written at the specified time/iteration.   |
+--------------------------------+--------------------------------------------+
| ``ORIENTATIONS``               | Molecule orientation information should be |
|                                | written at the specified time/iteration.   |
+--------------------------------+--------------------------------------------+
| ``ALL_DATA``                   | Equivalent to using both ``POSITIONS`` and |
|                                | ``ORIENTATION``.                           |
+--------------------------------+--------------------------------------------+

There are two possible visualization output file formats. ``DREAMM_V3``
mode is the default, and creates files in native ``DX`` format. This
mode is optimized for speed of visualization, but creates many
individual files. It has a directory structure with the top-level
directory given by adding ``_viz_data`` to the filename\ *, i.e.*
*filename*\ ``_viz_data``. For example if
``FILENAME = "./viz_data/diffusion_box"`` then the directory
``diffusion_box_viz_data`` will be created inside the ``./viz_data``
directory. Inside the *filename*\ ``_viz_data`` directory there is the
data directory called ``frame_data`` and three files:

-  *filename*\ ``.dx`` (the header file)
-  *filename*\ ``.iteration_numbers.bin``
-  *filename*\ ``.time_values.bin``

The directory ``frame_data`` contains a number of sub-directories named
by combining the word ``iteration_`` with the iteration number of the
simulation, such as ``iteration_0``, ``iteration_20``, etc. Each of
these iteration sub-directories by itself contains up to nine files:

-  ``meshes.dx`` (header file for meshes)
-  ``mesh_positions.bin``
-  ``mesh_states.bin`` (*optional*)
-  ``region_indices.bin``
-  ``surface_molecules.dx`` (header file for surface molecules)
-  ``surface_molecules_orientations.bin``
-  ``surface_molecules_positions.bin``
-  ``surface_molecules_states.bin`` (*optional*)
-  ``volume_molecules.dx`` (header file for volume molecules)
-  ``volume_molecules_orientations.bin``
-  ``volume_molecules_positions.bin``
-  ``volume_molecules_states.bin`` (*optional*)

Visualization data output for the ``DREAMM_V3_GROUPED`` mode is in
native ``DX`` format and includes one master header file and seven
binary data files, plus up to two optional data files if state values
are specified in the ``NAME_LIST`` blocks:

-  *filename*\ ``.dx`` (the master header file)
-  *filename*\ ``.mesh_positions.bin``
-  *filename*\ ``.mesh_states.bin`` (*optional*)
-  *filename*\ ``.region_indices.bin``
-  *filename*\ ``.molecule_positions.bin``
-  *filename*\ ``.molecule_orientations.bin``
-  *filename*\ ``.molecule_states.bin`` (*optional*)
-  *filename*\ ``.iteration_numbers.bin``
-  *filename*\ ``.time_values.bin``

Because the ``DREAMM_V3_GROUPED`` mode produces a small number of files, they
each may become very large. Hence, reading the files may be slow, but this mode
may be best for use on production (supercomputer) machines to avoid
transferring large number of files.

All of the keywords in the ``VIZ_OUTPUT`` block are optional except
``FILENAME.`` If the user does not specify the ``FILENAME`` keyword an error
message is printed and the simulation aborted. Some of the binary files for
both formats may be empty. For example, if no regions are defined the file
``region_indices.bin`` will be empty. Similarly, if no meshes or molecules are
defined the corresponding ``mesh_positions.bin`` or all molecules related
binary files will be empty. This avoids unintentional mixing of pre-existing
and new files that could result during several runs if incomplete file sets
were to be generated with the same names. In DReAMM, the user will only need to
point to the *filename*\ ``.dx`` file, and the data from the binary files will
be imported as needed for different frames. While using checkpointing in case
of the ``DREAMM_V3_GROUPED`` format the resulting visualization output files
add the checkpoint sequence number to their names, like *filename*\
``.mesh_positions.1.bin``. When checkpointing using the ``DREAMM_V3`` format,
new ``iteration_``\ *#* subdirectories holding the additional simulation output
will be created and the files *filename*\ ``.iteration_numbers.bin``\ *,
filename*\ ``.time_values.bin``\ *,* and *filename*\ ``.dx`` will be updated in
place to reflect these changes.

Examples of ``VIZ_OUTPUT`` statements are given below.

**Short-hand #1 (time style):**

.. code-block:: none

     VIZ_OUTPUT {
      FILENAME = "viz_data/output_example"
      MESHES {
        NAME_LIST { ALL_MESHES /* or list of object names */ }
        TIME_POINTS { ALL_DATA @ [0] }
      }
      MOLECULES {
        NAME_LIST { ALL_MOLECULES /* or list of molecule names */ }
        TIME_POINTS { ALL_DATA @ ALL_TIMES }
      }
    }

**Short-hand #2 (iterations style):**

.. code-block:: none

     VIZ_OUTPUT {
      FILENAME = "viz_data/output_example"
      MESHES {
        NAME_LIST { ALL_MESHES /* or list of object names */ }
        ITERATION_NUMBERS { ALL_DATA @ [0] }
      }
      MOLECULES {
        NAME_LIST { ALL_MOLECULES /* or list of molecule names */ }
        ITERATION_NUMBERS { ALL_DATA @ ALL_ITERATIONS }
      }
    }

**Expanded case:**

.. code-block:: none

     VIZ_OUTPUT {
      FILENAME = "viz_data/output_example"
      MESHES {
        NAME_LIST { ALL_MESHES /* or list of object names */ }
        TIME_POINTS {
          GEOMETRY @ [0]
          REGION_DATA @ [0]
        }
      }
      MOLECULES {
        NAME_LIST { ALL_MOLECULES /* or list of molecule names */ }
        TIME_POINTS {
          POSITIONS @ ALL_TIMES
          ORIENTATIONS @ ALL_TIMES
        }
      }
    }

Each ``MESHES / NAME_LIST`` statement may contain a single mesh object name or
multiple mesh object names with optional state values. It is left to the user
to avoid possible confusion arising from overlapping object trees within a
single master header file and its associated data files.

.. _util_commands:

Utility commands
--------------------

MCell3 understands the standard numeric operations ``+ - * /`` as well as the
following standard numerical functions:

.. math::

.. tabularcolumns:: |p{9cm}|p{7cm}|

+---------------------------------------------+-----------------------------------------------------------------------------+
| **Numerical Command**                       | **Explanation**                                                             |
+=============================================+=============================================================================+
| ``SQRT(`` *x* ``)``                         | Return the square root of *x*                                               |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``EXP(`` *x* ``)``                          | Return the value of *e* raised to the :math:`x^{th}` power                  |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``LOG(`` *x* ``)``                          | Return the natural logarithm of *x*                                         |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``LOG10(`` *x* ``)``                        | Return the base 10 logarithm of *x*                                         |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``SIN(`` *x* ``)``                          | Return the sine of *x*                                                      |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``COS(`` *x* ``)``                          | Return the cosine of *x*                                                    |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``TAN(`` *x* ``)``                          | Return the tangent of *x*                                                   |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``ASIN(`` *x* ``)``                         | Return the inverse sine of *x*                                              |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``ACOS(`` *x* ``)``                         | Return the inverse cosine of *x*                                            |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``ATAN(`` *x* ``)``                         | Return the inverse tangent of *x*                                           |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``ABS(`` *x* ``)``                          | Return the absolute value of *x*                                            |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``CEIL(`` *x* ``)``                         | Return the smallest integer at least as big as *x*                          |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``FLOOR(`` *x* ``)``                        | Return the largest integer at no bigger than *x*                            |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``MAX(`` *x* ``,`` *y* ``)``                | Return the larger of *x* and *y*                                            |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``MIN(`` *x* ``,`` *y* ``)``                | Return the smaller of *x* and *y*                                           |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``RAND_UNIFORM``                            | Return a random number uniformly distributed between 0 and 1                |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``RAND_GAUSSIAN``                           | Return a random number from a Gaussian distribution with mean 0 and         |
|                                             | standard deviation 1.                                                       |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``PI``                                      | The numeric value :math:`{\pi}=3.14159265358979323846`.                     |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``SEED``                                    | The value of the random number generator seed                               |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``printf(format,var1,...)``                 | ``printf`` works similar to C's printf statement given a ``format`` string  |
|                                             | and corresponding variables ``var1``, ... . Since MCell treats all defined  |
|                                             | variables as doubles only floating point formats should be used in the      |
|                                             | format string otherwise the results are undefined.                          |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``sprintf(out_string, format, var1, ...)``  | Same as printf but the result is written to a string variable               |
|                                             | ``out_string``.                                                             |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``fprintf(file_stream, format, var1, ...)`` | Same as printf but the result is written to a filestream object             |
|                                             | ``file_stream``.                                                            |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``file_stream = fopen(filename, mode)``     | Open a file ``filename`` using ``mode``. The returned filestream object     |
|                                             | ``file_stream`` can be used in ``fprintf``. Supported modes are "r" (read), |
|                                             | "w" (write), and "a" (append). ``filename`` should be a quoted string of    |
|                                             | the file location.                                                          |
+---------------------------------------------+-----------------------------------------------------------------------------+
| ``fclose(file_stream)``                     | Close the filestream object ``file_stream.``                                |
+---------------------------------------------+-----------------------------------------------------------------------------+

At any outer block in MCell3, one can define variables simply by assigning a
value to the name of the variable. E.g.  ``my_lucky_number=13`` would be a
valid (if unusual) way to define a variable. Variables can take numeric, array,
or string values. String values consist of text between double quotes. Strings
can be combined with the ``&`` operator. Array values are lists of numbers
inside brackets separated by commas, or starting and ending values plus a step
size, as exemplified below (note the double brackets):

.. code-block:: none

    my_lucky_number = 13
    my_favorite_array = [1,3,5,7,11,17]
    my_second_favorite_array = [[1.3 TO 2.75 STEP 0.331]]
    my_boring_string = "la la la, la la la" & ", la la"

To turn the random number generator seed into a string that you can use as part
of a filename, use

``sprintf(my_string_name,"%g",SEED)``.

If you want it to be a fixed width, e.g. 3 characters padded with zeros, use
the appropriate format string, e.g. ``"%03g"``.

MCell3 comments are delimited by ``/*`` and ``*/`` and can be nested. MDL files
can include other MDL files using the following syntax:

.. tabularcolumns:: |p{6cm}|p{10cm}|

+----------------------------------------+------------------------------------+
| **Command**                            | **Explanation**                    |
+========================================+====================================+
| ``INCLUDE_FILE = "`` *filename* ``"``  | Parse the text in *filename* as if |
|                                        | it were inserted into this MDL     |
|                                        | file at this point.                |
+----------------------------------------+------------------------------------+

Paths are relative to the location that MCell was run from, not relative to the
MDL file being parsed.

.. _tech_details:

Technical details affecting simulation speed and accuracy
===========================================================

.. _partitioning:

Partitioning
----------------

In future releases, MCell3 will automatically partition space to improve
execution speed. Currently, however, this must be performed manually. In
general, partitions should be chosen to avoid having too many surfaces and
molecules in one subvolume defined by the partitions. Molecules that are
specified as ``TARGET_ONLY`` or which do not interact with other molecules
diffusing in 3D need only have relatively few surfaces in one subvolume.

If there are few surfaces and/or molecules in a subvolume, it is advantageous
to have the subvolume as large as possible. Crossing partition boundaries takes
a small amount of time, so it is rarely useful to have partitions more finely
spaced than the average diffusion distance of the faster-moving molecules in
the simulation.

In cases where the diffusing molecules do not interact with each other, they
can safely take extended time-steps by measuring how far they are from things
they could interact with. In this case, the partitions with no surfaces should
be as large as possible. For example, a box works well with partitions just
inside its outer walls.

Finally, note that partition placement is not exact. The model is divided into
16384 possible partition boundaries, so partitions may shift by up to about one
part in twenty thousand of the size of the model. For instance, if the model
has a structure that is :math:`6{\mu}m` long, partitions may vary by about
:math:`0.0003{\mu}m`. Thus, do not place partitions too close to objects in
your model or they may not appear on the side you expect them to appear.

.. _avoid_coincident_meshes:

Avoid Coincident Meshes
---------------------------

Coincident meshes are mesh regions that overlap in space exactly. Coincident
meshes are problematic since they may lead to ambiguities which MCell can not
resolve properly. For example, how should MCell treat two coincident meshes
that are simultaneously transparent and reflective to a certain molecule. This
could happen when the face of a transparent counting box coincides exactly with
the face of a reflective bounding surface of the model. Problems can also arise
during checkpointing when MCell attempts to place surface molecules onto
available surfaces meshes. In general, coincident meshes (including counting
boxes or regions) should be avoided if at all possible, e.g., by separating
them by a small distance.

.. _mean_diff_distance:

Mean diffusion distance
---------------------------

Diffusion in MCell3 (and in earlier versions of MCell) is modeled as a series
of motions in a straight line. This is a good approximation around geometry
that is of a larger scale than the mean diffusion length for the time-step of
the molecule in question. For accurate results around intricate geometry, it
may be necessary to reduce the time step (or space step).

.. _rxn_probs_molecule_lifetimes:

Reaction probabilities and Molecule Lifetimes
-------------------------------------------------

MCell3 assigns a probability to each reaction in the simulation. These
probabilities are computed to match the bulk mass action reaction rates
specified in the MDL file. However, tracking of mass action behavior will
become impossible if some of the computed probabilities go beyond a value of
1.0 and simulations will then fail to reproduce the expected mass action
results. Importantly, even if the reaction probabilities that are computed and
displayed at start-up are all smaller than 1.0, internal correction factors can
further increase the actual probability values beyond a value of 1.0.
Therefore, MCell will output a warning if any reaction probability goes above
the value specified in ``HIGH_PROBABILITY_THRESHOLD`` which is 1.0 by default.
If warnings are given (and possibly even if they are not), one should reduce
the time step to lower the probabilities and see if the same results are
generated. If not, simulations should be run with shorter time steps in order
to avoid overly high probabilities.

Short lived reaction species may lead to inaccurate reaction equilibria both in
the volume and on surfaces. Experience has shown that a minimum lifetime of at
least 50 iterations is typically required to obtain reliable reaction
estimates. By default, MCell will warn if a species lives less than 50
iterations (see ``LIFETIME_THRESHOLD`` keyword) and users are strongly advised
to ensure that the lifetimes of molecules in their simulations are longer than
that.

Unimolecular reactions with half-lives of less than one time step are also not
perfectly accurate. Although unimolecular transitions will always occur at the
right rate, other molecules may not experience the right effective
concentration of each state, since a short-lifetime species may not be
converted to another species until the end of the time step after which many
other molecules may have had a chance to interact with it. Thus, the
shortest-lifetime species in a series of unimolecular transitions should not
have a half-life of less than approximately one time step if other molecules
can interact with that state.

.. _interaction_radii:

Interaction radii
---------------------

Bimolecular reactions occur within a distance specified by the
``INTERACTION_RADIUS`` command. In many cases, one may want to increase or
decrease this value. In particular, in order to get the right probability of
reaction, MCell3 increases the probability of reaction when near surfaces.

If ``ACCURATE_3D_REACTIONS`` is set to ``FALSE``, MCell3 also treats partition
boundaries as opaque and increases the probability of reaction rather than
looking for molecules on the other side of the partition.  This speeds
execution time but can lead to error, the reaction rate has approximately 1-2%
error if the average spacing between surfaces is at least 10 times the
interaction radius, and the reaction probabilities are 0.3 or less. For
example, if one has partitions spaced :math:`0.02{\mu}m` apart, simulation
accuracy will be poor with the default interaction radius of
:math:`0.01{\mu}m`. Thus, one might wish to specify
``INTERACTION_RADIUS=0.001``.

.. _placing_molecules_in_world:

Placing molecules in the world
----------------------------------

There are two ways to place molecules on surfaces: with a release site on a
region, and as part of the property of a surface or region. Release sites are
more flexible but slower; if you do not need the flexibility of release site
notation, you're better off defining a region and using the
``MOLECULE_DENSITY`` or ``MOLECULE_NUMBER`` commands to add molecules at
initialization.

All placement of molecules in volumes is done with release sites.  However, the
geometrical release sites (``CUBIC`` and ``SPHERICAL``) require less
computation to place each molecule. Thus, these should be used preferentially
for simple geometry. To release many particles at a one point, use a cubic
release site and set the diameter to 0. To release many particles at different
points, use the ``LIST`` release type.

.. _checkpointing_sims:

Checkpointing Simulations
===========================

MCell has the ability to checkpoint simulations, i.e., simulations can be
interrupted (checkpointed) and then restarted from where they left off.
Checkpointing can be used to divide long running simulations into shorter
segments or to change certain model parameters during a single simulation run
(see below). The basic MDL structure of a checkpointed simulation is as follows

.. code-block:: none

    ...
    CHECKPOINT_INFILE = "chkpt_in" 
    CHECKPOINT_OUTFILE = "chkpt_out" 
    CHECKPOINT_ITERATIONS = 10
    ...

Provided this input file, MCell will read the simulation state of a previous
simulation run stored in the file *chkpt_in*, continue to simulate for another
``CHECKPOINT_ITERATIONS`` iterations, and then save the simulation state to the
file *chkpt_out*. If ``CHECKPOINT_INFILE`` is omitted or the file *chkpt_in* is
absent (e.g., during the initial run of a checkpointed simulation) MCell will
start the simulation solely based on the information present in the provided
MDL file(s).

A checkpoint file contains all information needed by MCell to continue from a
previously checkpointed simulation. More specifically, it contains the current
iteration number, the current time, the state of the random number generator,
as well as the identities, orientations, and locations of all molecules in the
simulation. Hence, when MCell restarts from a previously saved checkpoint file
it will continue the simulation at the iteration and time given in the
checkpoint file until it reaches the last iteration given by ``ITERATION`` or
until the next checkpoint is due (controlled by ``CHECKPOINT_ITERATIONS``)
whatever happens first. All molecules specified in the checkpoint file will be
placed at their appropriate locations and orientations for surface molecules.
Since MCell checkpoint files contain the complete state of the random number
generator a checkpointed simulation will be identical to its uninterrupted
counterpart if no parameters are changed in between.

When restarting from a previously checkpointed simulation users may change
simulation parameters such as timestep, reaction rates, as well as mesh
geometries. Furthermore, new molecules may be added to the simulation. The
ability to change mesh geometries between checkpoints provides a limited
ability to simulate dynamic model geometries such as fusion pore opening.
However, when changing meshes between checkpoints several points need to be
kept in mind. Presently, checkpoint files only contain absolute molecule
positions and orientations for surface molecules without any reference to their
relative location inside the model geometry. For volume molecules this means
that it is up to the user to ensure that, e.g., molecule A remains inside a
certain region in the face of a changing mesh geometry. This is typically the
case for expanding meshes but may be problematic for shrinking ones. When
placing surface molecules from a previous run, MCell tries to 'snap' them to
the surface closest to their location given in the checkpoint file. Like for
volume molecules, it is up to the user to ensure that this leads to the
expected result. Finally, introducing new molecules at the start of a new
checkpoint has to be accomplished via a release pattern and the proper time
delay. Regular releases always happen at time 0 and will hence have no effect
during any but the initial checkpoint run.

Time Based Checkpointing
----------------------------------------

Instead of checkpointing at a specific iteration, one can alternatively create
a checkpoint at a set time with the ``CHECKPOINT_REALTIME`` command. The value
assigned to this is a series of integers separated by colons. The units and
formatting are illustrated below:

 - days:hours:minutes:seconds
 - hours:minutes:seconds
 - minutes:seconds
 - seconds

The basic MDL structure of a checkpointed simulation is as follows:

.. code-block:: none

    ...
    CHECKPOINT_INFILE = "chkpt_in" 
    CHECKPOINT_OUTFILE = "chkpt_out" 
    CHECKPOINT_REALTIME = 1:2:3:30 NOEXIT
    ...

In this example, 1:2:3:30 stands for 1 day, 2 hours, 3 minutes, and 30 seconds.
The ``NOEXIT`` command is optional, but, without it, MCell will exit after
creating the checkpoint file.

Checkpointing with SIGUSR1 and SIGUSR2
----------------------------------------

One can also create a checkpoint file as needed during a running simulation by
using the kill command with SIGUSR signals and MCell's PID.  With SIGUSR1,
MCell will create a checkpoint and continue running. With SIGUSR2, MCell will
create a checkpoint and end the simulation. As an example, if an MCell
simulation is running with PID 7984, then a checkpoint file could be created
like this::

  kill -SIGUSR1 7984

This feature is currently only available on Linux and Macs. 

The ps or top commands can be used to find MCell's PID.

.. _example_models:

Example models
================

.. _ligand_gated_ion_channel:

Ligand-gated ion channel
----------------------------

Below are a set of molecule definitions and reactions that specify an
ion channel that is gated by the binding of a single ligand.

.. code-block:: none

    DEFINE_MOLECULES {
      channel_unbound { D_2D=0 }
      channel_bound   { D_2D=0 }
      channel_open    { D_2D=0 }
      ligand          { D_3D=2e-8 }
      ion             { D_3D=3e-8 }
    }
    DEFINE_REACTIONS {
      channel_unbound' + ligand' -> channel_bound'             [1e7]
      channel_bound'             -> channel_unbound' + ligand' [2e2]
      channel_bound'             -> channel_open'              [5e2]
      channel_open'              -> channel_open'    + ion,    [8e4]
    }

We have defined a reaction where a ligand binds to one end of a channel
(presumably the extracellular face), which causes the channel to be in its
bound state. In that state it can either release the ligand or become open.
While open, it will emit ions on the other end (presumably the intracellular
face). This would be suitable if the ion concentration is much higher outside
than inside, or the membrane potential makes it highly favorable for the ion to
move inside, so that we don't have to worry about the reverse reaction. If
there is no electrical driving force, we might have to model ions both inside
and outside:

.. code-block:: none

    DEFINE_REACTIONS {
      channel_unbound' + ligand' -> channel_bound'             [1e7]
      channel_bound'             -> channel_unbound' + ligand' [2e2]
      channel_bound'             -> channel_open'              [5e2]
      channel_open'' + ion'      -> channel_open''     + ion,  [4e7]
    }

Here, the ion travels in either direction just as easily since it pays no
attention to the orientation of the channel. However, if there was a modest
driving force, traveling in might be easier than traveling out, which would be
reflected in the rates.

.. code-block:: none

    DEFINE_REACTIONS {
      channel_unbound' + ligand' -> channel_bound'             [1e7]
      channel_bound'             -> channel_unbound' + ligand' [2e2]
      channel_bound'             -> channel_open'              [5e2]
      channel_open' + ion'       -> channel_open'    + ion,    [4e8]
      channel_open' + ion,       -> channel_open'    + ion'    [1e8]
    }

In this case, the ion is four times as likely to travel from outside to inside
as inside to outside.

.. _example_bimolec_rxn:

Example bimolecular reaction
--------------------------------

Here's a complete MDL file that implements a simple bimolecular reaction that
should achieve equilibrium at 482 molecules of each species.

.. code-block:: none

    time_step = 1.0e-6

    TIME_STEP = time_step
    TIME_STEP_MAX = time_step
    ITERATIONS = 1e-2/time_step
    EFFECTOR_GRID_DENSITY = 10000
    INTERACTION_RADIUS = 0.001

    PARTITION_X = [ [-0.1 TO 0.1 STEP 0.01] ]
    PARTITION_Y = [ [-0.1 TO 0.1 STEP 0.01] ]
    PARTITION_Z = [ [-0.1 TO 0.1 STEP 0.01] ]

    DEFINE_MOLECULES
    {
      A { D_3D = 100e-8 }
      B { D_3D = 100e-8 }
      C { D_3D = 100e-8 }
    }

    /* Your basic reversible binding reaction */
    DEFINE_REACTIONS
    {
      A + B -> C [1e7]
      C -> A + B [1e3]
    }

    small_box BOX
    {
      CORNERS = [-0.1,-0.1,-0.1] , [0.1,0.1,0.1]
      /* REMOVE_ELEMENTS { TOP,LEFT } */  /* Could remove sides ... */
      /* REMOVE_ELEMENTS { INCLUDE_PATCH = [0.1,0,0] , [0.1,0.05,0.05] } /*... or patch*/
    }
    INSTANTIATE my_world OBJECT
    {
      A_release CUBIC_RELEASE_SITE {
        LOCATION=[0,0,0]
        MOLECULE=A
        NUMBER_TO_RELEASE=482
        SITE_DIAMETER=0.196
      }
      B_release CUBIC_RELEASE_SITE {
        LOCATION=[0,0,0]
        MOLECULE=B
        NUMBER_TO_RELEASE=482
        SITE_DIAMETER=0.196
      }
      C_release CUBIC_RELEASE_SITE {
        LOCATION=[0,0,0]
        MOLECULE=C
        NUMBER_TO_RELEASE=482
        SITE_DIAMETER=0.196
      }

     my_box OBJECT small_box {}
    }

    REACTION_DATA_OUTPUT
    {
      STEP = 1e-5
      { COUNT [A,WORLD] } => "eq_A.dat"
      { COUNT [B,WORLD] } => "eq_B.dat"
      { COUNT [C,WORLD] } => "eq_C.dat"
    }

.. _authors:

Authors
=========

The following authors have contributed to this document:

-  Tom Bartol
-  Jacob Czech
-  Markus Dittrich
-  Boris Kaminsky
-  Rex Kerr
-  Joel Stiles
