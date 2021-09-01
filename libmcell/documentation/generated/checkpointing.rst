.. _api-checkpointing:

*************
Checkpointing
*************
BaseChkptMol
============

Base class for checkpointed molecules.
Not to be used directly. All times are in seconds.

Attributes:
***********
.. _BaseChkptMol__id:

id: int
-------



.. _BaseChkptMol__species:

species: Species
----------------



.. _BaseChkptMol__diffusion_time:

diffusion_time: float
---------------------



.. _BaseChkptMol__birthday:

birthday: float
---------------



.. _BaseChkptMol__flags:

flags: int
----------



.. _BaseChkptMol__unimol_rxn_time:

unimol_rxn_time: float
----------------------

  | - default argument value in constructor: None

ChkptSurfMol
============

Class representing a checkpointed surface molecule.
Not to be used directly.

Attributes:
***********
.. _ChkptSurfMol__pos:

pos: Vec2
---------



.. _ChkptSurfMol__orientation:

orientation: Orientation
------------------------



.. _ChkptSurfMol__geometry_object:

geometry_object: GeometryObject
-------------------------------



.. _ChkptSurfMol__wall_index:

wall_index: int
---------------



.. _ChkptSurfMol__grid_tile_index:

grid_tile_index: int
--------------------



.. _ChkptSurfMol__id:

id: int
-------



.. _ChkptSurfMol__species:

species: Species
----------------



.. _ChkptSurfMol__diffusion_time:

diffusion_time: float
---------------------



.. _ChkptSurfMol__birthday:

birthday: float
---------------



.. _ChkptSurfMol__flags:

flags: int
----------



.. _ChkptSurfMol__unimol_rxn_time:

unimol_rxn_time: float
----------------------

  | - default argument value in constructor: None

ChkptVolMol
===========

Class representing a checkpointed volume molecule.
Not to be used directly.

Attributes:
***********
.. _ChkptVolMol__pos:

pos: Vec3
---------



.. _ChkptVolMol__id:

id: int
-------



.. _ChkptVolMol__species:

species: Species
----------------



.. _ChkptVolMol__diffusion_time:

diffusion_time: float
---------------------



.. _ChkptVolMol__birthday:

birthday: float
---------------



.. _ChkptVolMol__flags:

flags: int
----------



.. _ChkptVolMol__unimol_rxn_time:

unimol_rxn_time: float
----------------------

  | - default argument value in constructor: None

RngState
========

Internal checkpointing structure holding state of the random number generator.

Attributes:
***********
.. _RngState__randcnt:

randcnt: int
------------



.. _RngState__aa:

aa: int
-------



.. _RngState__bb:

bb: int
-------



.. _RngState__cc:

cc: int
-------



.. _RngState__randslr:

randslr: List[int]
------------------

  | Must contain RNG_SIZE items.


.. _RngState__mm:

mm: List[int]
-------------

  | Must contain RNG_SIZE items.


.. _RngState__rngblocks:

rngblocks: int
--------------

  | Must contain RNG_SIZE items.


