********
Geometry
********
Region
======

Represents region construted from 1 or more multiple, usually unnamed?

Attributes:
***********
* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods:
*********
* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of thwo regions


* | **__sub__**

   * | other: Region
   * | return type: Region


* | **__mul__**

   * | other: Region
   * | return type: Region



GeometryObject
==============

Attributes:
***********
* | **name**: str
  | Name of the object. Also represents BNGL compartment name if 'is_bngl_compartment' is True.

* | **vertex_list**: List[List[float]]
  | List of [x,y,z] triplets specifying positions of individual vertices.
  | Equivalent to List[Vec3] however, defining a constructor Vec3(List[float]) then 
  | tries to convert all lists of floats to Vec3

* | **wall_list**: List[List[int]]
  | List of [a,b,c] triplets specifying each wall, individual values are indices into the vertex list.
  | Equivalent to List[IVec3].

* | **is_bngl_compartment**: bool = False

* | **surface_compartment_name**: str = None

* | **surface_regions**: List[SurfaceRegion] = None

* | **surface_class**: SurfaceClass = None
  | Surface class for the whole object's surface. It is applied to the whole surface of this object 
  | except for those surface regions that have their specific surface class set explicitly.

* | **initial_surface_releases**: List[InitialSurfaceRelease] = None
  | Equivalent to MDL's MODIFY_SURFACE_REGIONS/MOLECULE_DENSITY or MOLECULE_NUMBER,
  | each item defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.

* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods:
*********
* | **translate**

   * | move: List[float]

  | Move object by a specified vector, must be done before model initialization.


* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of thwo regions


* | **__sub__**

   * | other: Region
   * | return type: Region


* | **__mul__**

   * | other: Region
   * | return type: Region



SurfaceRegion
=============

Surface region  in MDL, however a new class Region was instroduced in MCell4 so it was renamed 
to avoid confusion.

Attributes:
***********
* | **name**: str

* | **wall_indices**: List[int]
  | Surface region must be a part of a GeometryObject, items in this list are indices to 
  | its wall_list array

* | **surface_class**: SurfaceClass = None
  | Has higher priority than the parent geometry object's surface class.

* | **initial_surface_releases**: List[InitialSurfaceRelease] = None
  | Equivalent to MDL's MODIFY_SURFACE_REGIONS/MOLECULE_DENSITY or MOLECULE_NUMBER,
  | each item defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.

* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, when node_type is not Leaf, this is the right operand


Methods:
*********
* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of thwo regions


* | **__sub__**

   * | other: Region
   * | return type: Region


* | **__mul__**

   * | other: Region
   * | return type: Region



