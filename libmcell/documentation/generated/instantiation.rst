*************
Instantiation
*************
InitialSurfaceRelease
=====================

Defines molecules to be released onto a SurfaceRegion right when simulation starts

Attributes:
***********
* | **complex**: Complex

* | **number_to_release**: int = None
  | Number of molecules to be released onto a region,
  | only one of number_to_release and density can be set.

* | **density**: float = None
  | Density of molecules to be released onto a region,
  | only one of number_to_release and density can be set.

Instantiation
=============

Attributes:
***********
* | **release_sites**: List[ReleaseSite] = None

* | **geometry_objects**: List[GeometryObject] = None

* | **checkpointed_molecules**: List[BaseChkptMol] = None
  | Used when resuming simulation from a checkpoint.


Methods:
*********
* | **add_release_site**

   * | s: ReleaseSite

  | Makes a copy of the release site


* | **find_release_site**

   * | name: str
   * | return type: ReleaseSite


* | **add_geometry_object**

   * | o: GeometryObject

  | Makes a copy of the geometry object, in the future we will probably add some transformations


* | **find_geometry_object**

   * | name: str
   * | return type: GeometryObject


* | **find_volume_compartment_object**

   * | name: str
   * | return type: GeometryObject


* | **find_surface_compartment_object**

   * | name: str
   * | return type: GeometryObject


* | **load_bngl_seed_species**

   * | file_name: str
   * | default_release_region: Region = None
     | Used for seed species that have no compartments specified

   * | parameter_overrides: Dict[str, float] = None

  | Loads section seed species from a BNGL file and creates release sites according to it.
  | All elementary molecule types used in the seed species section must be already defined in subsystem.
  | If an item in the BNGL seed species section does not have its compartment set,
  | the argument default_region must be set and the molecules are released into or onto the 
  | default_region.



MoleculeReleaseInfo
===================

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
  | Argument must have exactly three floating point values.

ReleasePattern
==============

Attributes:
***********
* | **name**: str = None
  | Name of the release pattern

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
  | For unlimited number of trains use constant NUMBER_OF_TRAINS_UNLIMITED.

ReleaseSite
===========

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
  | Set automatically when

* | **region**: Region = None
  | Sets shape to Shape.REGION_EXPR.

* | **location**: List[float] = None

* | **site_diameter**: float = 0
  | For a geometrical release site, this releases molecules uniformly within
  | a radius r. Not used for releases on regions.
  | Usually required for Shape.List type of releases.

* | **site_radius**: float = None
  | For a geometrical release site, this releases molecules uniformly within
  | a radius r. Not used for releases on regions.

* | **number_to_release**: float = None
  | Only one of number_to_release, density, concentration or molecule_list can be set.
  | Value is truncated (floored) to an integer.

* | **density**: float = None
  | Unit is molecules per square micron (for surfaces).
  | Only one of number_to_release, density, concentration or molecule_list can be set.

* | **concentration**: float = None
  | Unit is molar (moles per liter) for volumes.
  | Only one of number_to_release, density, concentration or molecule_list can be set.

* | **release_probability**: float = None

