*************
Checkpointing
*************
BaseChkptMol
============

Base class for checkpointed molecules, not to be used directly. All times are in us (microseconds).

Attributes:
***********
* | **id**: int

* | **species**: Species

* | **diffusion_time**: float

* | **birthday**: float

* | **flags**: int

* | **unimol_rx_time**: float = None

ChkptSurfMol
============

Class representing a checkpointed surface molecule.
Not to be used directly.

Attributes:
***********
* | **pos**: Vec2

* | **orientation**: Orientation

* | **geometry_object**: GeometryObject

* | **wall_index**: int

* | **grid_tile_index**: int

* | **id**: int

* | **species**: Species

* | **diffusion_time**: float

* | **birthday**: float

* | **flags**: int

* | **unimol_rx_time**: float = None

ChkptVolMol
===========

Class representing a checkpointed volume molecule.
Not to be used directly.

Attributes:
***********
* | **pos**: Vec3

* | **id**: int

* | **species**: Species

* | **diffusion_time**: float

* | **birthday**: float

* | **flags**: int

* | **unimol_rx_time**: float = None

RngState
========

Internal checkpointing structure holding state of the random number generator.

Attributes:
***********
* | **randcnt**: int

* | **aa**: int

* | **bb**: int

* | **cc**: int

* | **randslr**: List[int]
  | Must contain RNG_SIZE items.

* | **mm**: List[int]
  | Must contain RNG_SIZE items.

* | **rngblocks**: int
  | Must contain RNG_SIZE items.

