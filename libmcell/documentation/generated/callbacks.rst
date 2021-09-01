.. _api-callbacks:

*********
Callbacks
*********
MolWallHitInfo
==============

Data structure passed to a callback function registered through
Model.register_mol_wall_hit_callback.

Example: `1300_wall_hit_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1300_wall_hit_callback/model.py>`_ 

Attributes:
***********
.. _MolWallHitInfo__molecule_id:

molecule_id: int
----------------

  | Id of molecule that hit the wall.


.. _MolWallHitInfo__geometry_object:

geometry_object: GeometryObject
-------------------------------

  | Object that was hit.


.. _MolWallHitInfo__wall_index:

wall_index: int
---------------

  | Index of the wall belonging to the geometry_object.


.. _MolWallHitInfo__time:

time: float
-----------

  | Time of the hit.


.. _MolWallHitInfo__pos3d:

pos3d: Vec3
-----------

  | Position of the hit.


.. _MolWallHitInfo__time_before_hit:

time_before_hit: float
----------------------

  | The time when the molecule started to diffuse towards the hit wall. 
  | It is either the start of the molecule's diffusion or 
  | when the molecule reflected from another wall.


.. _MolWallHitInfo__pos3d_before_hit:

pos3d_before_hit: Vec3
----------------------

  | Position of the molecule at time_before_hit.


ReactionInfo
============

Data structure passed to a reaction callback registered with 
Model.register_reaction_callback.

Example: `1800_vol_rxn_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1800_vol_rxn_callback/model.py>`_ 

Attributes:
***********
.. _ReactionInfo__type:

type: ReactionType
------------------

  | Specifies whether the reaction is unimolecular or bimolecular and
  | also provides information on reactant types.


.. _ReactionInfo__reactant_ids:

reactant_ids: List[int]
-----------------------

  | IDs of the reacting molecules, contains 1 ID for a unimolecular or a molecule+surface class reaction, 
  | 2 IDs for a bimolecular reaction.
  | For a bimolecular reaction, the first ID is always the molecule that diffused and the second one 
  | is the molecule that was hit.
  | IDs can be used to obtain the location of the molecules. The position of the first molecule obtained through 
  | model.get_molecule() is the position of the diffusing molecule before the collision.
  | All the reactants are removed after return from this callback, unless they are kept by the reaction such as A in A + B -> A + C.


.. _ReactionInfo__product_ids:

product_ids: List[int]
----------------------

  | IDs of reaction product molecules. They already exist in the simulated system together with reactants; however reactants 
  | will be removed after return from this callback.


.. _ReactionInfo__reaction_rule:

reaction_rule: ReactionRule
---------------------------

  | Reaction rule of the reaction that occured.


.. _ReactionInfo__time:

time: float
-----------

  | Time of the reaction.


.. _ReactionInfo__pos3d:

pos3d: Vec3
-----------

  | Specifies where reaction occurred in the 3d space, the specific meaning depends on the reaction type\:
  | - unimolecular reaction - position of the reacting molecule,
  | - volume-volume or surface-surface reaction - position of the first reactant,
  | - volume-surface reaction - position where the volume molecule hit the wall with the surface molecule.


.. _ReactionInfo__geometry_object:

geometry_object: GeometryObject
-------------------------------

  | The object on whose surface where the reaction occurred.
  | Set only for surface reactions or reactions with surface classes.
  | - default argument value in constructor: None

.. _ReactionInfo__wall_index:

wall_index: int
---------------

  | Set only for surface reactions or reactions with surface classes.
  | Index of wall belonging to the geometry_object where the reaction occured, 
  | i.e. wall where a volume molecule hit the surface molecule or
  | wall where the diffusing surface reactant reacted.
  | - default argument value in constructor: -1

.. _ReactionInfo__pos2d:

pos2d: Vec2
-----------

  | Set only for surface reactions or reactions with surface classes.
  | Specifies where reaction occurred in the 2d UV coordinates defined by the wall where the reaction occured, 
  | the rspecific meaning depends on the reaction type\:
  | - unimolecular reaction - position of the reacting molecule,
  | - volume-surface and surface-surface reaction - position of the second reactant.
  | - default argument value in constructor: None

