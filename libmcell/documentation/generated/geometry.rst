********
Geometry
********
GeometryObject
==============

Class represents geometry objects defined by triangular surface elements.

Example: `1330_get_wall/geometry.py <https://github.com/mcellteam/mcell_tests/tree/mcell4_dev/tests/pymcell4_positive/1330_get_wall/geometry.py>`_ 

Attributes:
***********
* | **name**: str
  | Name of the object. Also represents BNGL compartment name if 'is_bngl_compartment' is True.

* | **vertex_list**: List[List[float]]
  | List of [x,y,z] triplets specifying positions of individual vertices of each triangle.

* | **wall_list**: List[List[int]]
  | List of [a,b,c] triplets specifying each wall, individual values are indices into the 
  | vertex_list attribute.

* | **is_bngl_compartment**: bool = False
  | Set to True if this object represents a 3D BNGL compartment. 
  | Its name will be then the BNGL compartment name.

* | **surface_compartment_name**: str = None
  | When is_bngl_compartment is True, this attribute can be set to specify its 
  | membrane (2D) compartment name.

* | **surface_regions**: List[SurfaceRegion] = None
  | All surface regions associated with this geometry object.

* | **surface_class**: SurfaceClass = None
  | Surface class for the whole object's surface. It is applied to the whole surface of this object 
  | except for those surface regions that have their specific surface class set explicitly.

* | **initial_surface_releases**: List[InitialSurfaceRelease] = None
  | Each item in this list defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.

* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, do not use. When node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, do not use. When node_type is not Leaf, this is the right operand


Methods:
*********
* | **translate**

   * | move: List[float]
     | 3D vector [x, y, z] that will be added to each vertex of this object.


  | Move object by a specified vector. 
  | Cannot be called after model was initialized.

  | Example: `1400_object_translate/model.py <https://github.com/mcellteam/mcell_tests/tree/mcell4_dev/tests/pymcell4_positive/1400_object_translate/model.py>`_ 


* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of two regions, use with Python operator '+'.


* | **__sub__**

   * | other: Region
   * | return type: Region


  | Computes difference of two regions, use with Python operator '-'.


* | **__mul__**

   * | other: Region
   * | return type: Region


  | Computes intersection of two regions, use with Python operator '\*'.



Region
======

Represents region construted from 1 or more multiple, usually unnamed?

Attributes:
***********
* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, do not use. When node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, do not use. When node_type is not Leaf, this is the right operand


Methods:
*********
* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of two regions, use with Python operator '+'.


* | **__sub__**

   * | other: Region
   * | return type: Region


  | Computes difference of two regions, use with Python operator '-'.


* | **__mul__**

   * | other: Region
   * | return type: Region


  | Computes intersection of two regions, use with Python operator '\*'.



SurfaceRegion
=============

Defines a region on the object. The extent of a region is given by the wall_indices list. 
Molecules can be added and surface properties can be set with the optional regional surface commands. 
You can have an arbitrary number of regions on an object, and they may overlap if
you wish. Molecules added to overlapping regions accumulate. Triangles belonging to 
multiple regions inherit all parent regions’ surface properties. Users
have to make sure that in case of overlapped regions their surface properties
are compatible.

Example: `1700_linear_conc_gradient_w_conc_clamp/geometry.py <https://github.com/mcellteam/mcell_tests/tree/mcell4_dev/tests/pymcell4_positive/1700_linear_conc_gradient_w_conc_clamp/geometry.py>`_ 

Attributes:
***********
* | **name**: str
  | Name of this region.

* | **wall_indices**: List[int]
  | Surface region must be a part of a GeometryObject, items in this list are indices to 
  | its wall_list array.

* | **surface_class**: SurfaceClass = None
  | Optional surface class assigned to this surface region.
  | If not set, it is inherited from the parent heometry object's surface_class.

* | **initial_surface_releases**: List[InitialSurfaceRelease] = None
  | Each item of this list defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.

* | **node_type**: RegionNodeType = RegionNodeType.UNSET
  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.

* | **left_node**: Region = None
  | Internal, do not use. When node_type is not Leaf, this is the left operand

* | **right_node**: Region = None
  | Internal, do not use. When node_type is not Leaf, this is the right operand


Methods:
*********
* | **__add__**

   * | other: Region
   * | return type: Region


  | Computes union of two regions, use with Python operator '+'.


* | **__sub__**

   * | other: Region
   * | return type: Region


  | Computes difference of two regions, use with Python operator '-'.


* | **__mul__**

   * | other: Region
   * | return type: Region


  | Computes intersection of two regions, use with Python operator '\*'.



