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
  | Value DEFAULT means NONE for volume complexes and UP for surface complexes.


Notification
============

* | **NONE** = 0
* | **BRIEF** = 1
* | **FULL** = 2

WarningLevel
============

* | **IGNORE** = 0
  | Do something sensible and continue silently.

* | **WARNING** = 1
  | Do something sensible but emit a warning message.

* | **ERROR** = 2
  | Treat the warning as an error and stop.


VizMode
=======

* | **ASCII** = 0
  | Readable molecule visualization output.

* | **CELLBLENDER_V1** = 1
  | Binary molecule visualization output used by MCell3, format v1.
  | Allows only limited length of species name (256 chars) and 
  | does not contain molecule IDs.

* | **CELLBLENDER** = 2
  | Binary molecule visualization output, format v2.


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
  | This surface class does not do anything by itself, but it can be used as a reactant in 
  | reaction rules.

* | **REFLECTIVE** = 2
  | If used as a surface property for a volume molecule it is reflected by any surface with
  | this surface class. This is the default behavior for volume molecules.
  | If used for a surface molecule it is reflected by the border of the
  | surface with this surface class. 
  | Setting orientation in affected_complex_pattern allows selective reflection of volume 
  | molecules from only the front or back of a surface or selective reflection of surface 
  | molecules with only a certain orientation from the surface’s border. 
  | Using m.ALL_MOLECULES as affected_complex_pattern has the effect that all 
  | volume molecules are reflected by surfaces with this surface class and all surface molecules 
  | are reflected by the border of the surfaces with this surface class. 
  | Using m.ALL_VOLUME_MOLECULES as affected_complex_pattern has the effect that all
  | volume molecules are reflected by surfaces with this surface class. 
  | Using m.ALL_SURFACE_MOLECULES as affected_complex_pattern has the effect that all
  | surface molecules are reflected by the border of the surface with this surface class.

* | **TRANSPARENT** = 3
  | If used as a surface property for a volume molecule it passes through all surfaces with
  | this surface class.  
  | If used for a surface molecule it passes through the border of the surface with this surface 
  | class. This is the default behavior for surface molecules.
  | Setting orientation in affected_complex_pattern allows selective transparency of volume 
  | molecules from only the front or back of a surface or selective transparency for surface 
  | molecules with only a certain orientation from the surface’s border. 
  | To make a surface with this surface class transparent to all volume molecules,
  | use m.ALL_VOLUME_MOLECULES for affected_complex_pattern. 
  | To make a border of the surface with this surface class transparent to all surface molecules,
  | use m.ALL_SURFACE_MOLECULES for the affected_complex_pattern. 
  | Using m.ALL_MOLECULES for affected_complex_pattern has the effect that surfaces with this surface class 
  | are transparent to all volume molecules and borders of the surfaces with this surface class are 
  | transparent to all surface molecules.

* | **ABSORPTIVE** = 4
  | If affected_complex_pattern refers to a volume molecule it is destroyed if it touches surfaces with this surface class. 
  | If affected_complex_pattern refers to a surface molecule it is destroyed if it touches the border of the surface with 
  | this surface class, i.e., it is allowed to release surface molecules on absorptive surfaces, they get destroyed only
  | when they touch the border of this surface. 
  | Tick marks on name allow destruction from only one side of the surface for volume molecules or selective destruction 
  | for surface molecules on the surfaces’s border based on their orientation. 
  | To make a surface with this surface class absorptive to all volume molecules, m.ALL_VOLUME_MOLECULES 
  | can be used for affected_complex_pattern. 
  | To make a border of the surface with this surface class absorptive to all surface molecules,
  | m.ALL_SURFACE_MOLECULES can be used for name. 
  | Using m.ALL_MOLECULES as affected_complex_pattern has the effect that surfaces with this surface
  | class are absorptive for all volume molecules and borders of the surfaces with this surface class 
  | are absorptive for all surface molecules.

* | **CONCENTRATION_CLAMP** = 5
  | Clamps concentration at a surface by periodically releasing molecules that correspond
  | to the wall being a transparent boundary to the area with given concentration, 
  | and by absorbing all molecules that hit this surface.  
  | The molecules matching affected_complex_pattern are destroyed if they touch the surface (as if they
  | had passed through), and new molecules are created at the surface, as if molecules had passed through 
  | from the other side at a concentration value (units = M). 
  | Orientation marks may be used; in this case, the other side of the surface is reflective. 
  | Note that this command is only used to set the effective concentration of a volume molecule at a surface; 
  | it is not valid to specify a surface molecule.

* | **FLUX_CLAMP** = 6
  | Clamps flux at a surface by periodically releasing molecules that correspond
  | to the wall being a transparent boundary to the area with given concentration. 
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
  | Represents cases when a component must not be bound in a pattern.

* | **BOND_BOUND**: int = -2
  | Represents bond type !+ in a pattern.

* | **BOND_ANY**: int = -3
  | Represents bond type !? in a pattern.

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
  | This is a special floating-point value that means that an argument was not set, 
  | its value is 3.40282346638528859812e+38F.

* | **RNG_SIZE**: int = 256
  | Size of arrays of



