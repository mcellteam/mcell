**********
Submodules
**********
bngl_utils
==========

Submodule bngl_utils provides a function to load parameters from a BNGL file.


Methods:
*********
* | **load_bngl_parameters**

   * | file_name: str
   * | parameter_overrides: Dict[str, float] = None
   * | return type: Dict[str, float]



geometry_utils
==============

Submodule geometry_utils provides several functions to define 
model geometry. Rather limited for now.


Methods:
*********
* | **create_box**

   * | name: str
     | Name of the created geometry object.

   * | edge_length: float
     | Specifies length of each edge of the box in um.

   * | return type: GeometryObject


  | Creates a GeometryObject whose center is at (0, 0, 0).



run_utils
=========

Submodule run_utils provides functions used by checkpointing.


Methods:
*********
* | **get_last_checkpoint_dir**

   * | seed: int
   * | return type: str


  | Searches the directory checkpoints for the last checkpoint for the given 
  | parameters and returns the directory name if such a directory exists. 
  | Returns empty string if no checkpoint directory was found.
  | Currently supports only the seed argument.


* | **remove_cwd**

   * | paths: List[str]
   * | return type: List[str]


  | Removes all directory names items pointing to the current working directory from a list and 
  | returns a new list.



