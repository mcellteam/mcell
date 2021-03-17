*******************
Enums and Constants
*******************

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
* | **REACTIVE** = 1
  | This surface class does not do anything by itself but it can be used as a reactant in 
  | reaction rules.

* | **REFLECTIVE** = 2
* | **TRANSPARENT** = 3
* | **ABSORPTIVE** = 4
* | **CONCENTRATION_CLAMP** = 5
  | Clamps concentration at a surface by periodically releasing molecules that correspond
  | to the wall being a transparent boundary to area with given concentration, 
  | and by absorbing all molecules that hit this surface.

* | **FLUX_CLAMP** = 6
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

* | **DEFAULT_COUNT_BUFFER_SIZE**: int = 100
  | Internal constant used to initialize buffer size for molecule and reaction counts.

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



