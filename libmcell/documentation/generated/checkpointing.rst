*************
Checkpointing
*************
BaseChkptMol
============

All times are in us (microseconds).

Attributes:
***********
* | **id**: int

* | **species**: Species

* | **diffusion_time**: float

* | **birthday**: float

* | **flags**: int

* | **unimol_rx_time**: float = None

ChkptVolMol
===========

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

ChkptSurfMol
============

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

