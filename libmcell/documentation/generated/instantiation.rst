*************
Instantiation
*************
InitialSurfaceRelease
=====================

Defines molecules to be released onto a SurfaceRegion right when simulation starts

Attributes:
***********
* | **complex**: Complex
  | Defines the species of the molecule that will be released.

* | **number_to_release**: int = None
  | Number of molecules to be released onto a region,
  | only one of number_to_release and density can be set.

* | **density**: float = None
  | Density of molecules to be released onto a region,
  | only one of number_to_release and density can be set.

Instantiation
=============

Container used to hold instantiation-related model data. 
Instantiation is usually specific for each model, defines 
the geometry and initial setup of molecule releases.

Attributes:
***********
* | **release_sites**: List[ReleaseSite] = None
  | List of release sites to be included in the model.

* | **geometry_objects**: List[GeometryObject] = None
  | List of geometry objects to be included in the model.

* | **checkpointed_molecules**: List[BaseChkptMol] = None
  | Used when resuming simulation from a checkpoint.


Methods:
*********
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
     | Used for seed species that have no compartments specified.

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



MoleculeReleaseInfo
===================

Defines a pair (molecule, location). Used in ReleaseSite when its shape is Shape.LIST.

Attributes:
***********
* | **complex**: Complex
  | Complex instance defining the molecule that will be released.
  | Orientation of the complex instance is used to define orientation of the released molecule,
  | when Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and
  | surface molecules are released with Orientation.UP.
  | Compartment must not be set because this specific release definition states the location.

* | **location**: List[float]
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
* | **name**: str = None
  | Name of the release pattern.

* | **release_interval**: float = TIME_INFINITY
  | During a train of releases, release molecules after every t seconds. 
  | Default is to release only once.

* | **train_duration**: float = TIME_INFINITY
  | The train of releases lasts for t seconds before turning off. 
  | Default is to never turn off.

* | **train_interval**: float = TIME_INFINITY
  | A new train of releases happens every t seconds. 
  | Default is to never have a new train. 
  | The train interval must not be shorter than the train duration.

* | **number_of_trains**: int = 1
  | Repeat the release process for n trains of releases. Default is one train.
  | For unlimited number of trains use a constant NUMBER_OF_TRAINS_UNLIMITED.

ReleaseSite
===========

Defines a release site that specifies where, when and how should molecules be released.

Attributes:
***********
* | **name**: str
  | Name of the release site

* | **complex**: Complex = None
  | Defines the species of the molecule that will be released. Not used for the LIST shape. 
  | Must be set when molecule_list is empty and unset when molecule_list is not empty.
  | Orientation of the complex instance is used to define orientation of the released molecule,
  | when Orientation.DEFAULT is set, volume molecules are released with Orientation.NONE and
  | surface molecules are released with Orientation.UP.
  | When compartment is specified, this sets shape to Shape.COMPARTMENT and the molecules are released 
  | into the compartment.

* | **molecule_list**: List[MoleculeReleaseInfo] = None
  | Used for LIST shape release mode. 
  | Only one of number_to_release, density, concentration or molecule_list can be set.

* | **release_time**: float = 0
  | Specifies time in seconds when the release event is executed.
  | In case when a release pattern is used, this is the time of the first release.      
  | Equivalent to MDL's RELEASE_PATTERN command DELAY.

* | **release_pattern**: ReleasePattern = None
  | Use the release pattern to define schedule of releases. 
  | The default is to release the specified number of molecules at the set release_time.

* | **shape**: Shape = Shape.UNSET
  | Defines how the molecules shoudl be released. 
  | Set automatically for these cases to the following values\: 
  | region is set - Shape.REGION_EXPR,
  | region is not set and complex uses a compartment - Shape.COMPARTMENT,
  | molecule_list is set - Shape.LIST,
  | location is set - Shape.SPHERICAL.

* | **region**: Region = None
  | Defines a volume or surface region where to release molecules. 
  | Setting it sets shape to Shape.REGION_EXPR.

* | **location**: List[float] = None
  | Defines center of a sphere where to release molecules. 
  | Setting it sets shape to Shape.SPHERICAL.

* | **site_diameter**: float = 0
  | For a geometrical release site, this releases molecules uniformly within
  | a radius r computed as site_diameter/2. 
  | Used only when shape is Shape.SPHERICAL.
  | Maximum one of site_diameter or site_radius may be set.

* | **site_radius**: float = None
  | For a geometrical release site, this releases molecules uniformly within
  | a radius site_radius.
  | Used only when shape is Shape.SPHERICAL.
  | Maximum one of site_diameter or site_radius may be set.

* | **number_to_release**: float = None
  | Sets number of molecules to release. Cannot be set when shape is Shape.LIST. 
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Value is truncated (floored) to an integer.

* | **density**: float = None
  | Unit is molecules per square micron (for surfaces). 
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Cannot be set when shape is Shape.LIST.

* | **concentration**: float = None
  | Unit is molar (moles per liter) for volumes.
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Cannot be set when shape is Shape.LIST.

