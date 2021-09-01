.. _api-submodules:

**********
Submodules
**********
bngl_utils
==========

Submodule bngl_utils provides a function to load parameters from a BNGL file.


Methods:
*********
.. _bngl_utils__load_bngl_parameters:

load_bngl_parameters (file_name: str, parameter_overrides: Dict[str, float]=None) -> Dict[str, float]
-----------------------------------------------------------------------------------------------------


  | Load parameters section from a BNGL file and return it as a dictionary name->value.

* | file_name: str
  | Path to the BNGL file to be loaded.

* | parameter_overrides: Dict[str, float] = None
  | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
  | its value is ignored and instead value parameter_overrides[k] is used.

  | Example: `2010_bng_parameter_override/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/2010_bng_parameter_override/model.py>`_ 



data_utils
==========

Submodule data_utils provides data loading and manipulation functions.


Methods:
*********
.. _data_utils__load_dat_file:

load_dat_file (file_name: str) -> List[List[float]]
---------------------------------------------------


  | Loads a two-column file where the first column is usually time and the second is a 
  | floating point value. Returns a two-column list. 
  | Can be used to load a file with variable rate constants.

* | file_name: str
  | Path to the .dat file to be loaded.



geometry_utils
==============

Submodule geometry_utils provides several functions to define 
model geometry. Rather limited for now.


Methods:
*********
.. _geometry_utils__create_box:

create_box (name: str, edge_dimension: float=None, xyz_dimensions: List[float]=None) -> GeometryObject
------------------------------------------------------------------------------------------------------


  | Creates a GeometryObject in the shape of a cube whose center is at (0, 0, 0).

* | name: str
  | Name of the created geometry object.

* | edge_dimension: float = None
  | Specifies length of each edge of the box in um. 
  | None of x/y/z dimensions can be set.

* | xyz_dimensions: List[float] = None
  | Specifies x/y/z sizes of the box in um. Parameter edge_dimension must not be set.

  | Examples: `1400_rel_site_for_each_it/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1400_rel_site_for_each_it/model.py>`_ `1105_point_release_w_create_box_not_cube/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1105_point_release_w_create_box_not_cube/model.py>`_ 


.. _geometry_utils__create_icosphere:

create_icosphere (name: str, radius: float, subdivisions: int) -> GeometryObject
--------------------------------------------------------------------------------


  | Creates a GeometryObject in the shape of an icosphere whose center is at (0, 0, 0).

* | name: str
  | Name of the created geometry object.

* | radius: float
  | Specifies radius of the sphere.

* | subdivisions: int
  | Number of subdivisions from the initial icosphere. 
  | The higher this value will be the smoother the icosphere will be.
  | Allowed range is between 1 and 8.

  | Example: `1110_point_release_w_create_icosphere/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1110_point_release_w_create_icosphere/model.py>`_ 


.. _geometry_utils__validate_volumetric_mesh:

validate_volumetric_mesh (model: Model, geometry_object: GeometryObject)
------------------------------------------------------------------------


  | Checks that the mesh was correctly analyzed, that it has volume and 
  | all edges have neighboring walls.
  | Must be called after model initialization. 
  | Throws exception with detained message if validation did not pass.

* | model: Model
  | Model object after initialization.

* | geometry_object: GeometryObject
  | Geometry object to be checked.



run_utils
=========

Submodule run_utils provides functions used by checkpointing.


Methods:
*********
.. _run_utils__get_last_checkpoint_dir:

get_last_checkpoint_dir (seed: int) -> str
------------------------------------------


  | Searches the directory checkpoints for the last checkpoint for the given 
  | parameters and returns the directory name if such a directory exists. 
  | Returns empty string if no checkpoint directory was found.
  | Currently supports only the seed argument.

* | seed: int

.. _run_utils__remove_cwd:

remove_cwd (paths: List[str]) -> List[str]
------------------------------------------


  | Removes all directory names items pointing to the current working directory from a list and 
  | returns a new list.

* | paths: List[str]


