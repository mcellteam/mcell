**********
Submodules
**********
bngl_utils
==========


Methods:
*********
* | **load_bngl_parameters**

   * | file_name: str
   * | parameter_overrides: Dict[str, float] = None
   * | return type: Dict[str, float]



run_utils
=========


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



geometry_utils
==============


Methods:
*********
* | **create_box**

   * | name: str
     | Name of the created geometry object

   * | edge_length: float
     | Specifies length of each edge of the box.

   * | return type: GeometryObject


  | Creates a GeometryObject whose center is at (0, 0, 0).



