.. _api-instantiation:

*************
Instantiation
*************
InitialSurfaceRelease
=====================

Defines molecules to be released onto a SurfaceRegion right when simulation starts

Attributes:
***********
.. _InitialSurfaceRelease__complex:

complex: Complex
----------------

  | Defines the species of the molecule that will be released.


.. _InitialSurfaceRelease__number_to_release:

number_to_release: int
----------------------

  | Number of molecules to be released onto a region,
  | only one of number_to_release and density can be set.
  | - default argument value in constructor: None

.. _InitialSurfaceRelease__density:

density: float
--------------

  | Density of molecules to be released onto a region,
  | only one of number_to_release and density can be set.
  | - default argument value in constructor: None

Instantiation
=============

Container used to hold instantiation-related model data. 
Instantiation is usually specific for each model, defines 
the geometry and initial setup of molecule releases.

Example: `pymcell4/1250_organelle_move <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1250_organelle_move>`_ 

Attributes:
***********
.. _Instantiation__release_sites:

release_sites: List[ReleaseSite]
--------------------------------

  | List of release sites to be included in the model.
  | - default argument value in constructor: None

.. _Instantiation__geometry_objects:

geometry_objects: List[GeometryObject]
--------------------------------------

  | List of geometry objects to be included in the model.
  | - default argument value in constructor: None

.. _Instantiation__checkpointed_molecules:

checkpointed_molecules: List[BaseChkptMol]
------------------------------------------

  | Used when resuming simulation from a checkpoint.
  | - default argument value in constructor: None


Methods:
*********
.. _Instantiation__add_release_site:

add_release_site (s: ReleaseSite)
---------------------------------


  | Adds a reference to the release site s to the list of release sites.

* | s: ReleaseSite

.. _Instantiation__find_release_site:

find_release_site (name: str) -> ReleaseSite
--------------------------------------------


  | Finds a release site by its name, returns None if no such release site is present.

* | name: str

.. _Instantiation__add_geometry_object:

add_geometry_object (o: GeometryObject)
---------------------------------------


  | Adds a reference to the geometry object o to the list of geometry objects.

* | o: GeometryObject

.. _Instantiation__find_geometry_object:

find_geometry_object (name: str) -> GeometryObject
--------------------------------------------------


  | Finds a geometry object by its name, returns None if no such geometry object is present.

* | name: str

.. _Instantiation__find_volume_compartment_object:

find_volume_compartment_object (name: str) -> GeometryObject
------------------------------------------------------------


  | Finds a geometry object by its name, the geometry object must be a BNGL compartment.
  | Returns None if no such geometry object is present.

* | name: str

.. _Instantiation__find_surface_compartment_object:

find_surface_compartment_object (name: str) -> GeometryObject
-------------------------------------------------------------


  | Finds a geometry object that is a BNGL compartment and its surface name is name.
  | Returns None if no such geometry object is present.

* | name: str

.. _Instantiation__load_bngl_compartments_and_seed_species:

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



MoleculeReleaseInfo
===================

Defines a pair (molecule, location). Used in ReleaseSite when its shape is Shape.LIST.

Attributes:
***********
.. _MoleculeReleaseInfo__complex:

complex: Complex
----------------

  | Complex instance defining the molecule that will be released.
  | Orientation of the complex instance is used to define orientation of the released molecule,
  | when Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and
  | surface molecules are released with Orientation.UP.
  | Compartment must not be set because this specific release definition states the location.


.. _MoleculeReleaseInfo__location:

location: List[float]
---------------------

  | 3D position where the molecule will be released. 
  | If a molecule has a 2D diffusion constant, it will be
  | placed on the surface closest to the coordinate given. 
  | Argument must have exactly three floating point values [x, y, z].


ReleasePattern
==============

Defines a release pattern that specifies repeating molecule releases. 
Can be used by a ReleaseSite.

Attributes:
***********
.. _ReleasePattern__name:

name: str
---------

  | Name of the release pattern.
  | - default argument value in constructor: None

.. _ReleasePattern__release_interval:

release_interval: float
-----------------------

  | During a train of releases, release molecules after every t seconds. 
  | Default is to release only once.
  | - default argument value in constructor: TIME_INFINITY

.. _ReleasePattern__train_duration:

train_duration: float
---------------------

  | The train of releases lasts for t seconds before turning off. 
  | Default is to never turn off.
  | - default argument value in constructor: TIME_INFINITY

.. _ReleasePattern__train_interval:

train_interval: float
---------------------

  | A new train of releases happens every t seconds. 
  | Default is to never have a new train. 
  | The train interval must not be shorter than the train duration.
  | - default argument value in constructor: TIME_INFINITY

.. _ReleasePattern__number_of_trains:

number_of_trains: int
---------------------

  | Repeat the release process for n trains of releases. Default is one train.
  | For unlimited number of trains use a constant NUMBER_OF_TRAINS_UNLIMITED.
  | - default argument value in constructor: 1

ReleaseSite
===========

Defines a release site that specifies where, when and how should molecules be released.

Example: `1100_point_release/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1100_point_release/model.py>`_ 

Attributes:
***********
.. _ReleaseSite__name:

name: str
---------

  | Name of the release site


.. _ReleaseSite__complex:

complex: Complex
----------------

  | Defines the species of the molecule that will be released. Not used for the LIST shape. 
  | Must be set when molecule_list is empty and unset when molecule_list is not empty.
  | Orientation of the complex instance is used to define orientation of the released molecule,
  | when Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and
  | surface molecules are released with Orientation.UP.
  | When compartment is specified and region is not set, this sets shape to Shape.COMPARTMENT and 
  | the molecules are released into the compartment.
  | When this is a release of volume molecules, and both compartment and region are set, 
  | this sets shape to Shape.REGION_EXPR and the target region is the intersection 
  | of the region and the compartment.
  | - default argument value in constructor: None

.. _ReleaseSite__molecule_list:

molecule_list: List[MoleculeReleaseInfo]
----------------------------------------

  | Used for LIST shape release mode. 
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | - default argument value in constructor: None

.. _ReleaseSite__release_time:

release_time: float
-------------------

  | Specifies time in seconds when the release event is executed.
  | In case when a release pattern is used, this is the time of the first release.      
  | Equivalent to MDL's RELEASE_PATTERN command DELAY.
  | - default argument value in constructor: 0

.. _ReleaseSite__release_pattern:

release_pattern: ReleasePattern
-------------------------------

  | Use the release pattern to define schedule of releases. 
  | The default is to release the specified number of molecules at the set release_time.
  | - default argument value in constructor: None

.. _ReleaseSite__shape:

shape: Shape
------------

  | Defines how the molecules shoudl be released. 
  | Set automatically for these cases to the following values\: 
  | region is set - Shape.REGION_EXPR,
  | region is not set and complex uses a compartment - Shape.COMPARTMENT,
  | molecule_list is set - Shape.LIST,
  | location is set - Shape.SPHERICAL.
  | - default argument value in constructor: Shape.UNSET

.. _ReleaseSite__region:

region: Region
--------------

  | Defines a volume or surface region where to release molecules. 
  | Setting it sets shape to Shape.REGION_EXPR. 
  | When this is a release of volume molecules, and both compartment and region are set, 
  | this sets shape to Shape.REGION_EXPR and the target region is the intersection 
  | of the region and the compartment.
  | - default argument value in constructor: None

.. _ReleaseSite__location:

location: List[float]
---------------------

  | Defines center of a sphere where to release molecules. 
  | Setting it sets shape to Shape.SPHERICAL.
  | - default argument value in constructor: None

.. _ReleaseSite__site_diameter:

site_diameter: float
--------------------

  | For a geometrical release site, this releases molecules uniformly within
  | a radius r computed as site_diameter/2. 
  | Used only when shape is Shape.SPHERICAL.
  | Maximum one of site_diameter or site_radius may be set.
  | - default argument value in constructor: 0

.. _ReleaseSite__site_radius:

site_radius: float
------------------

  | For a geometrical release site, this releases molecules uniformly within
  | a radius site_radius.
  | Used only when shape is Shape.SPHERICAL.
  | Maximum one of site_diameter or site_radius may be set.
  | - default argument value in constructor: None

.. _ReleaseSite__number_to_release:

number_to_release: float
------------------------

  | Sets number of molecules to release. Cannot be set when shape is Shape.LIST. 
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Value is truncated (floored) to an integer.
  | - default argument value in constructor: None

.. _ReleaseSite__density:

density: float
--------------

  | Unit is molecules per square micron (for surfaces). 
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Cannot be set when shape is Shape.LIST.
  | - default argument value in constructor: None

.. _ReleaseSite__concentration:

concentration: float
--------------------

  | Unit is molar (moles per liter) for volumes.
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Cannot be set when shape is Shape.LIST.
  | - default argument value in constructor: None

.. _ReleaseSite__release_probability:

release_probability: float
--------------------------

  | This release does not occur every time, but rather with probability p. 
  | Either the whole release occurs or none of it does; the probability does not 
  | apply molecule-by-molecule. release_probability must be in the interval [0, 1].
  | - default argument value in constructor: 1

