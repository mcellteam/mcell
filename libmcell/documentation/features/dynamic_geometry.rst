.. _dynamic_geometry_section:

****************
Dynamic Geometry
****************

Overview
########

Dynamic geometry allows manipulation with individual mesh geometry vertices during simulation.
This section provides an overview of related methods.

NOTE: this documentation is not generated automatically and may get out of date. 
Check the generated API reference when needed. 

Vertex Representation
#####################

In MCell, vertices are represented with class Vec3.
It has attributes: x, y, z, and a method .tolist() that returns the coordinates as a Python list with three elements. 
Class Vec3 also provides usual operators such as +, -, *, /, ==.


Geometry Introspection
######################

This subsection lists several methods used for introspecting geometry state:

* :ref:`Model.get_vertex <Model__get_vertex>` - Returns coordinates of a vertex as Vec3.
* :ref:`Model.get_wall <Model__get_wall>` - Returns information about a wall belonging to a given object. 
* :ref:`Model.get_vertex_unit_normal <Model__get_vertex_unit_normal>` - Returns sum of all wall normals that use this vertex normalized to a unit vector of length 1um. This represents the unit vector pointing outwards from the vertex.
* :ref:`Model.get_wall_unit_normal <Model__get_wall_unit_normal>` - Returns wall normal of length 1um.


Geometry Modifications
######################

* :ref:`Model.add_vertex_move <Model__add_vertex_move>` - Adds a displacement for a given object's vertex, only stored until apply_vertex_moves is called.
* :ref:`Model.apply_vertex_moves <Model__apply_vertex_moves>` - Applies all the vertex moves specified with add_vertex_move call. See API documentation for more details.


Wall Hit Callbacks
##################

* :ref:`Model.register_mol_wall_hit_callback <Model__register_mol_wall_hit_callback>` - 

Special Attributes and Methods
##############################

* :ref:`Wall.is_movable <Wall__is_movable>` - The attribute says whether the wall can be moved through calls to Model.apply_vertex_moves. If False, all changes to vertices are ignored.  
* :ref:`Model.pair_molecules <Model__pair_molecules>` - Allows to pair surface molecules on different objects. When a wall of one molecule is moved, wall of the second paired molecule is also moved. See :ref:`Model.apply_vertex_moves <Model__apply_vertex_moves>` for more details.  



