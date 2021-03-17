*********
Callbacks
*********
MolWallHitInfo
==============

Attributes:
***********
* | **molecule_id**: int

* | **geometry_object**: GeometryObject
  | Object that was hit.

* | **wall_index**: int
  | Index of wall belonging to the geometry_object.

* | **time**: float
  | Time of the hit

* | **pos3d**: Vec3
  | Position of the hit

* | **time_before_hit**: float
  | Time when the molecule started to diffuse towards the hit wall. 
  | It is either the start of the molecule's diffusion or 
  | if a wall was hit later then the time of last wall hit.

* | **pos3d_before_hit**: Vec3
  | Position of the molecule at time_before_hit

ReactionInfo
============

Data structure passed to a reaction callback.

Attributes:
***********
* | **type**: ReactionType
  | Specifies whether the reaction is unimolecular or bimolecular and
  | also provides information in reactant types.

* | **reactant_ids**: List[int]
  | IDs of the reacting molecules, contains 1 ID for a unimolecular reaction, 2 IDs for a bimolecular reaction.
  | For a bimolecular reaction, the first ID is always the molecule that was diffused and the second one 
  | is the molecule that was hit.
  | IDs can be used to obtain location of the molecules. The position of the first molecule obtained through 
  | model.get_molecule() is the position of the diffusing molecule before the collision.
  | All the reactants are removed after return from this callback, unless they are kept by the reaction such as in A + B -> A + C.

* | **product_ids**: List[int]
  | IDs of reaction product molecules. They already exist in the simulated system together with reactants, however reactants 
  | will be removed after return from this callback.

* | **reaction_rule**: ReactionRule
  | Reaction rule of the reaction.

* | **time**: float
  | Time of the reaction

* | **pos3d**: Vec3
  | Specifies where reaction occured in the 3d space, specific meaning depends on the reaction type\:
  | - unimolecular reaction - position of the reacting molecule,
  | - volume-volume or surface-surface reaction - position of the first reactant,
  | - volume-surface reaction - position where the volume molecule hit the wall with the surface molecule.

* | **geometry_object**: GeometryObject = None
  | Set only for surface reactions.
  | Object on whose surface where the reaction occured.

* | **wall_index**: int = -1
  | Set only for surface reactions.
  | Index of wall belonging to the geometry_object where the reaction occured, 
  | i.e. where the volume molecule hit the wall with a surface molecule or
  | wall where the diffusing surface reactant reacted.

* | **pos2d**: Vec2 = None
  | Set only for surface reactions.
  | Specifies where reaction occured in the 2d UV coordinates defined by the wall where the reaction occured, 
  | specific meaning depends on the reaction type\:
  | - unimolecular reaction - position of the reacting molecule,
  | - volume-surface and surface-surface reaction - position of the second reactant.

