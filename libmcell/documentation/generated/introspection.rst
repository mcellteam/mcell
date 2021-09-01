.. _api-introspection:

*************
Introspection
*************
Introspection
=============

Only internal. This class is used only as a base class to Model, it is not provided through API. Defines interface to introspect simulation state.


Methods:
*********
.. _Introspection__get_molecule_ids:

get_molecule_ids (pattern: Complex=None) -> List[int]
-----------------------------------------------------


  | Returns a list of ids of molecules.
  | If the arguments pattern is not set, the list of all molecule ids is returned.  
  | If the argument pattern is set, the list of all molecule ids whose species match 
  | the pattern is returned.

* | pattern: Complex = None
  | BNGL pattern to select molecules based on their species, might use compartments.

  | Example: `1910_get_molecule_ids_w_pattern/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1910_get_molecule_ids_w_pattern/model.py>`_ 


.. _Introspection__get_molecule:

get_molecule (id: int) -> Molecule
----------------------------------


  | Returns a information on a molecule from the simulated environment, 
  | None if the molecule does not exist.

* | id: int
  | Unique id of the molecule to be retrieved.

  | Example: `1900_molecule_introspection/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1900_molecule_introspection/model.py>`_ 


.. _Introspection__get_species_name:

get_species_name (species_id: int) -> str
-----------------------------------------


  | Returns a string representing canonical species name in the BNGL format.

* | species_id: int
  | Id of the species.

  | Example: `1850_run_unimol_rxn_in_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1850_run_unimol_rxn_in_callback/model.py>`_ 


.. _Introspection__get_vertex:

get_vertex (object: GeometryObject, vertex_index: int) -> Vec3
--------------------------------------------------------------


  | Returns coordinates of a vertex.

* | object: GeometryObject
* | vertex_index: int
  | This is the index of the vertex in the geometry object's walls (wall_list).

  | Example: `1340_get_vertex/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1340_get_vertex/model.py>`_ 


.. _Introspection__get_wall:

get_wall (object: GeometryObject, wall_index: int) -> Wall
----------------------------------------------------------


  | Returns information about a wall belonging to a given object.

* | object: GeometryObject
  | Geometry object whose wall to retrieve.

* | wall_index: int
  | This is the index of the wall in the geometry object's walls (wall_list).

  | Example: `1330_get_wall/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1330_get_wall/model.py>`_ 


.. _Introspection__get_vertex_unit_normal:

get_vertex_unit_normal (object: GeometryObject, vertex_index: int) -> Vec3
--------------------------------------------------------------------------


  | Returns sum of all wall normals that use this vertex converted to a unit vector of 
  | length 1 um (micrometer).
  | This represents the unit vector pointing outwards from the vertex.

* | object: GeometryObject
  | Geometry object whose vertex to retrieve.

* | vertex_index: int
  | This is the index of the vertex in the geometry object's vertex_list.

  | Example: `1320_get_vertex_unit_normal/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1320_get_vertex_unit_normal/model.py>`_ 


.. _Introspection__get_wall_unit_normal:

get_wall_unit_normal (object: GeometryObject, wall_index: int) -> Vec3
----------------------------------------------------------------------


  | Returns wall normal converted to a unit vector of length 1um.

* | object: GeometryObject
  | Geometry object whose wall's normal to retrieve.

* | wall_index: int
  | This is the index of the vertex in the geometry object's walls (wall_list).

  | Example: `1310_get_wall_unit_normal/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1310_get_wall_unit_normal/model.py>`_ 


.. _Introspection__get_wall_color:

get_wall_color (object: GeometryObject, wall_index: int) -> Color
-----------------------------------------------------------------


  | Returns color of a wall.

* | object: GeometryObject
  | Geometry object whose wall's color to retrieve.

* | wall_index: int
  | This is the index of the vertex in the geometry object's walls (wall_list).


.. _Introspection__set_wall_color:

set_wall_color (object: GeometryObject, wall_index: int, color: Color)
----------------------------------------------------------------------


  | Sets color of a wall.

* | object: GeometryObject
  | Geometry object whose wall's color to retrieve.

* | wall_index: int
  | This is the index of the vertex in the geometry object's walls (wall_list).

* | color: Color
  | Color to be set.



Molecule
========

Representation of a molecule obtained from Model 
during simulation obtained through Model.get_molecule.
Changes through changing attributes of this object are not allowed except 
for complete removal of this molecule.

Example: `1900_molecule_introspection/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1900_molecule_introspection/model.py>`_ 

Attributes:
***********
.. _Molecule__id:

id: int
-------

  | Unique id of this molecule. MCell assigns this unique id to each created 
  | molecule. All reactions change ID of molecules even in reactions such as 
  | A\@CP -> A\@EC.
  | - default argument value in constructor: ID_INVALID

.. _Molecule__type:

type: MoleculeType
------------------

  | Type of this molecule, either volume or surface.
  | - default argument value in constructor: MoleculeType.UNSET

.. _Molecule__species_id:

species_id: int
---------------

  | Species id of this molecule.
  | The species_id value is only temporary. Species ids are created and removed as needed
  | automatically and if this species is removed, this particular species_id value 
  | won't be valid. This can happen when a following iteration is simulated.
  | - default argument value in constructor: ID_INVALID

.. _Molecule__pos3d:

pos3d: Vec3
-----------

  | Contains position of a molecule in 3D space.
  | - default argument value in constructor: None

.. _Molecule__orientation:

orientation: Orientation
------------------------

  | Contains orientation for surface molecule. Volume molecules 
  | have always orientation set to Orientation.NONE.
  | - default argument value in constructor: Orientation.NOT_SET

.. _Molecule__pos2d:

pos2d: Vec2
-----------

  | Set only for surface molecules. Position on a wall in UV coordinates 
  | relative to the triangle of the wall.
  | - default argument value in constructor: None

.. _Molecule__geometry_object:

geometry_object: GeometryObject
-------------------------------

  | Set only for surface molecules.
  | Is set to a reference to the geometry object on whose surface is the molecule located.
  | - default argument value in constructor: None

.. _Molecule__wall_index:

wall_index: int
---------------

  | Set only for surface molecules.
  | Index of wall belonging to the geometry_object where is the 
  | molecule located.
  | - default argument value in constructor: -1


Methods:
*********
.. _Molecule__remove:

remove ()
---------


  | Removes this molecule from simulation. Any subsequent modifications
  | of this molecule won't have any effect.

  | Example: `1920_molecule_remove/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1920_molecule_remove/model.py>`_ 



Wall
====

Constant representation of wall of a geometry object.
Changes through changing attributes of this object are not allowed
except for the attribute is_movable.

Example: `1330_get_wall/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1330_get_wall/model.py>`_ 

Attributes:
***********
.. _Wall__geometry_object:

geometry_object: GeometryObject
-------------------------------

  | Object to which this wall belongs.


.. _Wall__wall_index:

wall_index: int
---------------

  | Index of this wall in the object to which this wall belongs.


.. _Wall__vertices:

vertices: List[Vec3]
--------------------

  | Vertices of the triangle that represents this wall.


.. _Wall__area:

area: float
-----------

  | Area of the wall in um^2.


.. _Wall__unit_normal:

unit_normal: Vec3
-----------------

  | Normal of this wall with unit length of 1 um.
  | There is also a method Model.get_wall_unit_normal that allows to 
  | retrieve just the normal value without the need to prepare this 
  | whole Wall object.


.. _Wall__is_movable:

is_movable: bool
----------------

  | If True, whis wall can be moved through Model.apply_vertex_moves,
  | if False, wall moves are ignored. 
  | Can be set during simulation.
  | - default argument value in constructor: True

WallWallHitInfo
===============

This class is used in the return type of Model.apply_vertex_moves.
Contains pair of walls that collided.

Example: `1515_tetrahedron_box_collision_moving_3_w_wall_wall_hit/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1515_tetrahedron_box_collision_moving_3_w_wall_wall_hit/model.py>`_ 

Attributes:
***********
.. _WallWallHitInfo__wall1:

wall1: Wall
-----------

  | First colliding wall.


.. _WallWallHitInfo__wall2:

wall2: Wall
-----------

  | Second colliding wall.


