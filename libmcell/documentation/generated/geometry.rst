.. _api-geometry:

********
Geometry
********
Color
=====

Represents color with alpha component.
Provides two means to set value, either red, green, blue and alpha, 
or rgba. If both color individual components and rgba are set in initialization,
the individual components are used.

Attributes:
***********
.. _Color__red:

red: float
----------

  | Red component in range 0-1.
  | - default argument value in constructor: None

.. _Color__green:

green: float
------------

  | Green component in range 0-1.
  | - default argument value in constructor: None

.. _Color__blue:

blue: float
-----------

  | Blue component in range 0-1.
  | - default argument value in constructor: None

.. _Color__alpha:

alpha: float
------------

  | Alpha component in range 0-1. 1 means nontransparent.
  | - default argument value in constructor: 1

.. _Color__rgba:

rgba: int
---------

  | This attribute provides an alternative way of defining colors by supplying a 
  | 32-bit unsigned integer representation of the color with an aplha channel. 
  | In hexadecimal notation the first 2 digits are value for red, second 2 digits are 
  | green, third 2 digits are blue and the last two digits are alpha. 
  | The range for each component is thus 0x0-0xFF (0-255). 
  | Example\: 0x0000ffcc represents the same color as rgba(0, 0, 100%, 80%).
  | All values are valid.
  | - default argument value in constructor: 0

GeometryObject
==============

Class represents geometry objects defined by triangular surface elements.

Example: `1330_get_wall/geometry.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1330_get_wall/geometry.py>`_ 

Attributes:
***********
.. _GeometryObject__name:

name: str
---------

  | Name of the object. Also represents BNGL compartment name if 'is_bngl_compartment' is True.


.. _GeometryObject__vertex_list:

vertex_list: List[List[float]]
------------------------------

  | List of [x,y,z] triplets specifying positions of individual vertices of each triangle.


.. _GeometryObject__wall_list:

wall_list: List[List[int]]
--------------------------

  | List of [a,b,c] triplets specifying each wall, individual values are indices into the 
  | vertex_list attribute.


.. _GeometryObject__is_bngl_compartment:

is_bngl_compartment: bool
-------------------------

  | Set to True if this object represents a 3D BNGL compartment. 
  | Its name will be then the BNGL compartment name.
  | - default argument value in constructor: False

.. _GeometryObject__surface_compartment_name:

surface_compartment_name: str
-----------------------------

  | When is_bngl_compartment is True, this attribute can be set to specify its 
  | membrane (2D) compartment name.
  | - default argument value in constructor: None

.. _GeometryObject__surface_regions:

surface_regions: List[SurfaceRegion]
------------------------------------

  | All surface regions associated with this geometry object.
  | - default argument value in constructor: None

.. _GeometryObject__surface_class:

surface_class: SurfaceClass
---------------------------

  | Surface class for the whole object's surface. It is applied to the whole surface of this object 
  | except for those surface regions that have their specific surface class set explicitly.
  | - default argument value in constructor: None

.. _GeometryObject__initial_surface_releases:

initial_surface_releases: List[InitialSurfaceRelease]
-----------------------------------------------------

  | Each item in this list defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.
  | - default argument value in constructor: None

.. _GeometryObject__initial_color:

initial_color: Color
--------------------

  | Initial color for this geometry object. If a surface region has its color set, its value 
  | is used for the walls of that surface region.
  | - default argument value in constructor: None

.. _GeometryObject__node_type:

node_type: RegionNodeType
-------------------------

  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.
  | - default argument value in constructor: RegionNodeType.UNSET

.. _GeometryObject__left_node:

left_node: Region
-----------------

  | Internal, do not use. When node_type is not Leaf, this is the left operand
  | - default argument value in constructor: None

.. _GeometryObject__right_node:

right_node: Region
------------------

  | Internal, do not use. When node_type is not Leaf, this is the right operand
  | - default argument value in constructor: None


Methods:
*********
.. _GeometryObject__translate:

translate (move: List[float])
-----------------------------


  | Move object by a specified vector. 
  | Cannot be called after model was initialized.

* | move: List[float]
  | 3D vector [x, y, z] that will be added to each vertex of this object.

  | Example: `1400_object_translate/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1400_object_translate/model.py>`_ 


.. _GeometryObject____add__:

__add__ (other: Region) -> Region
---------------------------------


  | Computes union of two regions, use with Python operator '+'.

* | other: Region

.. _GeometryObject____sub__:

__sub__ (other: Region) -> Region
---------------------------------


  | Computes difference of two regions, use with Python operator '-'.

* | other: Region

.. _GeometryObject____mul__:

__mul__ (other: Region) -> Region
---------------------------------


  | Computes intersection of two regions, use with Python operator '\*'.

* | other: Region


Region
======

Represents region construted from 1 or more multiple, usually unnamed?

Attributes:
***********
.. _Region__node_type:

node_type: RegionNodeType
-------------------------

  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.
  | - default argument value in constructor: RegionNodeType.UNSET

.. _Region__left_node:

left_node: Region
-----------------

  | Internal, do not use. When node_type is not Leaf, this is the left operand
  | - default argument value in constructor: None

.. _Region__right_node:

right_node: Region
------------------

  | Internal, do not use. When node_type is not Leaf, this is the right operand
  | - default argument value in constructor: None


Methods:
*********
.. _Region____add__:

__add__ (other: Region) -> Region
---------------------------------


  | Computes union of two regions, use with Python operator '+'.

* | other: Region

.. _Region____sub__:

__sub__ (other: Region) -> Region
---------------------------------


  | Computes difference of two regions, use with Python operator '-'.

* | other: Region

.. _Region____mul__:

__mul__ (other: Region) -> Region
---------------------------------


  | Computes intersection of two regions, use with Python operator '\*'.

* | other: Region


SurfaceRegion
=============

Defines a region on the object. The extent of a region is given by the wall_indices list. 
Molecules can be added and surface properties can be set with the optional regional surface commands. 
You can have an arbitrary number of regions on an object, and they may overlap if
you wish. Molecules added to overlapping regions accumulate. Triangles belonging to 
multiple regions inherit all parent regions’ surface properties. Users
have to make sure that in case of overlapped regions their surface properties
are compatible.

Example: `1700_linear_conc_gradient_w_conc_clamp/geometry.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1700_linear_conc_gradient_w_conc_clamp/geometry.py>`_ 

Attributes:
***********
.. _SurfaceRegion__name:

name: str
---------

  | Name of this region.


.. _SurfaceRegion__wall_indices:

wall_indices: List[int]
-----------------------

  | Surface region must be a part of a GeometryObject, items in this list are indices to 
  | its wall_list array.


.. _SurfaceRegion__surface_class:

surface_class: SurfaceClass
---------------------------

  | Optional surface class assigned to this surface region.
  | If not set, it is inherited from the parent geometry object's surface_class.
  | - default argument value in constructor: None

.. _SurfaceRegion__initial_surface_releases:

initial_surface_releases: List[InitialSurfaceRelease]
-----------------------------------------------------

  | Each item of this list defines either density or number of molecules to be released on this surface 
  | regions when simulation starts.
  | - default argument value in constructor: None

.. _SurfaceRegion__initial_color:

initial_color: Color
--------------------

  | Initial color for this specific surface region. If not set, color of the parent's GeometryObject is used.
  | - default argument value in constructor: None

.. _SurfaceRegion__node_type:

node_type: RegionNodeType
-------------------------

  | When this values is LeafGeometryObject, then this object is of class GeometryObject,
  | when LeafSurfaceRegion, then it is of class SurfaceRegion.
  | - default argument value in constructor: RegionNodeType.UNSET

.. _SurfaceRegion__left_node:

left_node: Region
-----------------

  | Internal, do not use. When node_type is not Leaf, this is the left operand
  | - default argument value in constructor: None

.. _SurfaceRegion__right_node:

right_node: Region
------------------

  | Internal, do not use. When node_type is not Leaf, this is the right operand
  | - default argument value in constructor: None


Methods:
*********
.. _SurfaceRegion____add__:

__add__ (other: Region) -> Region
---------------------------------


  | Computes union of two regions, use with Python operator '+'.

* | other: Region

.. _SurfaceRegion____sub__:

__sub__ (other: Region) -> Region
---------------------------------


  | Computes difference of two regions, use with Python operator '-'.

* | other: Region

.. _SurfaceRegion____mul__:

__mul__ (other: Region) -> Region
---------------------------------


  | Computes intersection of two regions, use with Python operator '\*'.

* | other: Region


