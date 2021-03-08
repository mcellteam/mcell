*************
Introspection
*************
WallWallHitInfo
===============

Attributes:
***********
* | **wall1**: Wall

* | **wall2**: Wall

Molecule
========

This is a Python representation of a molecule obtained from Model 
during simulation.

Attributes:
***********
* | **id**: int = ID_INVALID
  | Unique id of this molecule

* | **type**: MoleculeType = MoleculeType.UNSET

* | **species_id**: int = ID_INVALID
  | Species id of this molecule.
  | The id value is only temporary and can be invalidated by simulating an iteration.

* | **pos3d**: Vec3 = None
  | Contains position of a molecule in 3D space.

* | **orientation**: Orientation = Orientation.NOT_SET
  | Contains orientation for surface molecule. Volume molecules 
  | have always orientation set to Orientation.NONE.

* | **pos2d**: Vec2 = None
  | Set only for surface molecules.

* | **geometry_object**: GeometryObject = None
  | Set only for surface molecules.
  | Object on whose surface is the molecule located.

* | **wall_index**: int = -1
  | Set only for surface molecules.
  | Index of wall belonging to the geometry_object where is the 
  | molecule located.


Methods:
*********
* | **remove**


  | Removes this molecule from simulation. Any subsequent modifications
  | of this object won't have any effect.



Wall
====

This is a Python representation of a molecule obtained from Model 
during simulation.

Attributes:
***********
* | **geometry_object**: GeometryObject
  | Object to which this wall belongs.

* | **wall_index**: int
  | Index of this wall in the object to which this wall belongs.

* | **vertices**: List[Vec3]
  | Vertices of the triangle that represents this wall.

* | **area**: float

* | **unit_normal**: Vec3
  | Normal of this wall with unit length of 1 um.
  | To get just the unit vector, not the whole wall, there is also method Model.get_wall_unit_normal.

* | **is_movable**: bool = True
  | If True, whis wall can be moved through Model.apply_vertex_moves,
  | if False, wall moves are ignored.

Introspection
=============

This class is used only as a base class to Model, it is not provided through API. Provides methods to introspect simulation state.


Methods:
*********
* | **get_molecule_ids**

   * | pattern: Complex = None
   * | return type: List[int]


  | Returns a list of ids of molecules.
  | If the arguments pattern is not set, the list of all molecule ids is returned.  
  | If the argument pattern is set, the list of all molecule ids whose species match 
  | the pattern is returned. Matching of patterns with compartments works exactly in the 
  | same was as in observables.


* | **get_molecule**

   * | id: int
   * | return type: Molecule


  | Returns a molecule from the simulated environment, None if the molecule does not exist


* | **get_species_name**

   * | species_id: int
   * | return type: str


  | Returns a string representing canonical species name in the BNGL format.


* | **get_vertex**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in object's walls (wall_list).

   * | return type: Vec3


  | Returns coordinates of a vertex.


* | **get_wall**

   * | object: GeometryObject
   * | wall_index: int
     | This is the index of the wall in object's walls (wall_list).

   * | return type: Wall


  | Returns information about a wall belonging to a given object.


* | **get_vertex_unit_normal**

   * | object: GeometryObject
   * | vertex_index: int
     | This is the index of the vertex in object's vertex_list.

   * | return type: Vec3


  | Returns sum of all wall normals that use this vertex converted to a unit vector of length 1um.
  | This represents the unit vector pointing outwards from the vertex.


* | **get_wall_unit_normal**

   * | object: GeometryObject
   * | wall_index: int
     | This is the index of the vertex in object's walls (wall_list).

   * | return type: Vec3


  | Returns wall normal converted to a unit vector of length 1um.



