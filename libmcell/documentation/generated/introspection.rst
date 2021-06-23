.. _api-introspection:

*************
Introspection
*************
Introspection
=============

Only internal. This class is used only as a base class to Model, it is not provided through API. Defines interface to introspect simulation state.


Methods:
*********
* | **get_molecule_ids**

   * | pattern: Complex = None
     | BNGL pattern to select molecules based on their species, might use compartments.

   * | return type: List[int]


  | Returns a list of ids of molecules.
  | If the arguments pattern is not set, the list of all molecule ids is returned.  
  | If the argument pattern is set, the list of all molecule ids whose species match 
  | the pattern is returned.

  | Example: `1910_get_molecule_ids_w_pattern/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1910_get_molecule_ids_w_pattern/model.py>`_ 


* | **get_molecule**

   * | id: int
     | Unique id of the molecule to be retrieved.

   * | return type: Molecule


  | Returns a information on a molecule from the simulated environment, 
  | None if the molecule does not exist.

  | Example: `1900_molecule_introspection/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1900_molecule_introspection/model.py>`_ 


* | **get_species_name**

   * | species_id: int
     | Id of the species.

   * | return type: str


  | Returns a string representing canonical species name in the BNGL format.

  | Example: `1850_run_unimol_rxn_in_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1850_run_unimol_rxn_in_callback/model.py>`_ 


* | **get_vertex**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | return type: Vec3


  | Returns coordinates of a vertex.

  | Example: `1340_get_vertex/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1340_get_vertex/model.py>`_ 


* | **get_wall**

   * | object: GeometryObject
     | Geometry object whose wall to retrieve.

   * | wall_index: int
     | This is the index of the wall in the geometry object's walls (wall_list).

   * | return type: Wall


  | Returns information about a wall belonging to a given object.

  | Example: `1330_get_wall/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1330_get_wall/model.py>`_ 


* | **get_vertex_unit_normal**

   * | object: GeometryObject
     | Geometry object whose vertex to retrieve.

   * | vertex_index: int
     | This is the index of the vertex in the geometry object's vertex_list.

   * | return type: Vec3


  | Returns sum of all wall normals that use this vertex converted to a unit vector of 
  | length 1 um (micrometer).
  | This represents the unit vector pointing outwards from the vertex.

  | Example: `1320_get_vertex_unit_normal/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1320_get_vertex_unit_normal/model.py>`_ 


* | **get_wall_unit_normal**

   * | object: GeometryObject
     | Geometry object whose wall's normal to retrieve.

   * | wall_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | return type: Vec3


  | Returns wall normal converted to a unit vector of length 1um.

  | Example: `1310_get_wall_unit_normal/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1310_get_wall_unit_normal/model.py>`_ 


* | **get_wall_color**

   * | object: GeometryObject
     | Geometry object whose wall's color to retrieve.

   * | wall_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | return type: Color


  | Returns color of a wall.


* | **set_wall_color**

   * | object: GeometryObject
     | Geometry object whose wall's color to retrieve.

   * | wall_index: int
     | This is the index of the vertex in the geometry object's walls (wall_list).

   * | color: Color
     | Color to be set.


  | Sets color of a wall.



Molecule
========

Representation of a molecule obtained from Model 
during simulation obtained through Model.get_molecule.
Changes through changing attributes of this object are not allowed except 
for complete removal of this molecule.

Example: `1900_molecule_introspection/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1900_molecule_introspection/model.py>`_ 

Attributes:
***********
* | **id**: int = ID_INVALID
  | Unique id of this molecule. MCell assigns this unique id to each created 
  | molecule. All reactions change ID of molecules even in reactions such as 
  | A\@CP -> A\@EC.

* | **type**: MoleculeType = MoleculeType.UNSET
  | Type of this molecule, either volume or surface.

* | **species_id**: int = ID_INVALID
  | Species id of this molecule.
  | The species_id value is only temporary. Species ids are created and removed as needed
  | automatically and if this species is removed, this particular species_id value 
  | won't be valid. This can happen when a following iteration is simulated.

* | **pos3d**: Vec3 = None
  | Contains position of a molecule in 3D space.

* | **orientation**: Orientation = Orientation.NOT_SET
  | Contains orientation for surface molecule. Volume molecules 
  | have always orientation set to Orientation.NONE.

* | **pos2d**: Vec2 = None
  | Set only for surface molecules. Position on a wall in UV coordinates 
  | relative to the triangle of the wall.

* | **geometry_object**: GeometryObject = None
  | Set only for surface molecules.
  | Is set to a reference to the geometry object on whose surface is the molecule located.

* | **wall_index**: int = -1
  | Set only for surface molecules.
  | Index of wall belonging to the geometry_object where is the 
  | molecule located.


Methods:
*********
* | **remove**


  | Removes this molecule from simulation. Any subsequent modifications
  | of this molecule won't have any effect.

  | Example: `1920_molecule_remove/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1920_molecule_remove/model.py>`_ 



Wall
====

Constant representation of wall of a geometry object.
Changes through changing attributes of this object are not allowed
except for the attribute is_movable.

Example: `1330_get_wall/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1330_get_wall/model.py>`_ 

Attributes:
***********
* | **geometry_object**: GeometryObject
  | Object to which this wall belongs.

* | **wall_index**: int
  | Index of this wall in the object to which this wall belongs.

* | **vertices**: List[Vec3]
  | Vertices of the triangle that represents this wall.

* | **area**: float
  | Area of the wall in um^2.

* | **unit_normal**: Vec3
  | Normal of this wall with unit length of 1 um.
  | There is also a method Model.get_wall_unit_normal that allows to 
  | retrieve just the normal value without the need to prepare this 
  | whole Wall object.

* | **is_movable**: bool = True
  | If True, whis wall can be moved through Model.apply_vertex_moves,
  | if False, wall moves are ignored. 
  | Can be set during simulation.

WallWallHitInfo
===============

This class is used in the return type of Model.apply_vertex_moves.
Contains pair of walls that collided.

Example: `1515_tetrahedron_box_collision_moving_3_w_wall_wall_hit/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/1515_tetrahedron_box_collision_moving_3_w_wall_wall_hit/model.py>`_ 

Attributes:
***********
* | **wall1**: Wall
  | First colliding wall.

* | **wall2**: Wall
  | Second colliding wall.

