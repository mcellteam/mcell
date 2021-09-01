.. _api-observables:

***********
Observables
***********
Count
=====

Represents a molecule or reaction count observable.
What is counted is defined through a CounTerm tree and a reference to 
the root of this tree is stored into attribute expression. 
This tree usually contains just one node.

Example: `1500_region_release_must_set_compartment/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1500_region_release_must_set_compartment/model.py>`_ 

Attributes:
***********
.. _Count__name:

name: str
---------

  | Name of a count may be specified when one needs to search for them later. 
  | When the count is created when a BNGL file is loaded, its name is set, for instance
  | when the following BNGL code is loaded\:
  | 
  | begin observables
  |    Molecules Acount A
  | end observables
  | 
  | the name is set to Acount.
  | - default argument value in constructor: None

.. _Count__file_name:

file_name: str
--------------

  | File name where this observable values will be stored.
  | File extension or setting explicit output_format determines the output format.
  | A) When not set, the value is set using seed during model initialization as follows: 
  | file_name = './react_data/seed_' + str(model.config.seed).zfill(5) + '/' + name + '.dat'
  | and the output format is set to CountOutputFormat.DAT in the constructor.
  | B) When the file_name is set explicitly by the user and the extension is .dat such as here:
  | file_name = './react_data/seed_' + str(SEED).zfill(5) + '/' + name + '.dat'
  | and the output format is set to CountOutputFormat.DAT in the constructor.
  | File names for individual Counts must be different.
  | C) When the file_name is set explicitly by the user and the extension is .gdat such as here:
  | file_name = './react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat'
  | and the output format is set to CountOutputFormat.GDAT in the constructor.
  | The file name is usually the same for all counts but one can 
  | create multiple gdat files with different observables.
  | All observables that are stored into a single .gdat file must have the same 
  | periodicity specified by attribute every_n_timesteps.
  | Must be set.
  | - default argument value in constructor: None

.. _Count__expression:

expression: CountTerm
---------------------

  | The expression must be set to a root of an expression tree composed of CountTerms. 
  | In the usual cases, there is just one CountTerm in this expression tree and its 
  | node_type is ExprNodeType.LEAF.
  | The count expression tree defines CountTerm objects that are added or subtracted
  | from each other.
  | - default argument value in constructor: None

.. _Count__multiplier:

multiplier: float
-----------------

  | In some cases it might be useful to multiply the whole count by a constant to get 
  | for instance concentration. The expression tree allows only addition and subtraction 
  | of count terms so such multiplication can be done through this attribute.
  | It can be also used to divide the resulting count by passing an inverse of the divisor (1/d).
  | - default argument value in constructor: 1

.. _Count__every_n_timesteps:

every_n_timesteps: float
------------------------

  | Specifies periodicity of this count's output.
  | Value is truncated (floored) to an integer.
  | If value is set to 0, this Count is used only on-demand through calls to its
  | get_current_value method.
  | - default argument value in constructor: 1

.. _Count__output_format:

output_format: CountOutputFormat
--------------------------------

  | Listed as the last attribute because the automatic default value
  | is sufficient in most cases. 
  | Selection of output format. Default setting uses file extension  
  | from attribute file_name. 
  | When set to CountOutputFormat.AUTOMATIC_FROM_EXTENSION, 
  | this output_format is set automatically only in the Count's constructor.
  | - default argument value in constructor: CountOutputFormat.AUTOMATIC_FROM_EXTENSION


Methods:
*********
.. _Count__get_current_value:

get_current_value () -> float
-----------------------------


  | Returns the current value for this count. Can be used to count both molecules and reactions.
  | Reaction counting starts at the beginning of the simulation.
  | The model must be initialized with this Count present as one of the observables.

  | Examples: `2600_get_current_mol_count/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/2600_get_current_mol_count/model.py>`_ `2650_get_current_rxn_count/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/2650_get_current_rxn_count/model.py>`_ 



CountTerm
=========

A count observable can be defined as an expression composed of addition
or subtraction individual count terms. This class represents one count term
in this expression.

Attributes:
***********
.. _CountTerm__species_pattern:

species_pattern: Complex
------------------------

  | Count the number of molecules that match the given complex instance pattern.
  | This corresponds to the BNGL 'Species' specifier in the BNGL seed species section.
  | Counts each molecule exactly once. 
  | If the pattern has a compartment set, this specifies the counted region.
  | Exactly one of species_pattern, molecules_pattern, and reaction_rule must be set.
  | - default argument value in constructor: None

.. _CountTerm__molecules_pattern:

molecules_pattern: Complex
--------------------------

  | Count the number of matches of the given pattern on molecules.
  | This corresponds to the BNGL 'Molecules' specifier in the BNGL seed species section.
  | The observable will increment the count every time the pattern matches the molecule.
  | For instance, pattern A will match a complex A(a!1).B(a!1,a!2).A(b!2) twice. 
  | When the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a 
  | molecule A(b.a!1).A(a!1,b!2).B(a!2) will be counted twice because the 
  | pattern may match in two different ways. 
  | If the pattern has a compartment set, the compartment is used to filter out the molecules.   
  | Exactly one of species_pattern, molecules_pattern, and reaction_rule must be set.
  | - default argument value in constructor: None

.. _CountTerm__reaction_rule:

reaction_rule: ReactionRule
---------------------------

  | Count the number of applications of this specific reactions that occurred since the
  | start of the simulation.
  | Exactly one of species_pattern, molecules_pattern, and reaction_rule must be set.
  | - default argument value in constructor: None

.. _CountTerm__region:

region: Region
--------------

  | Only a GeometryObject or SurfaceRegion can be passed as the region argument, 
  | compound regions (created with +, -, \*) are not supproted yet.   
  | Can be combined with a compartment specified in the species_pattern or molecules_pattern.
  | If compartment in species_pattern or molecules_pattern is not specified and 
  | region is left unset, counting is done in the whole world.
  | - default argument value in constructor: None

.. _CountTerm__node_type:

node_type: ExprNodeType
-----------------------

  | Internal, used to specify what type of count expression node this object represents.
  | - default argument value in constructor: ExprNodeType.LEAF

.. _CountTerm__left_node:

left_node: CountTerm
--------------------

  | Internal, when node_type is not Leaf, this is the left operand.
  | - default argument value in constructor: None

.. _CountTerm__right_node:

right_node: CountTerm
---------------------

  | Internal, when node_type is not Leaf, this is the right operand.
  | - default argument value in constructor: None

.. _CountTerm__initial_reactions_count:

initial_reactions_count: int
----------------------------

  | Used for checkpointing, allows to set initial count of reactions that occurred.
  | Ignored when molecules are counted.
  | - default argument value in constructor: 0


Methods:
*********
.. _CountTerm____add__:

__add__ (op2: CountTerm) -> CountTerm
-------------------------------------


  | Create a new CountTerm that represents addition of two count terms.
  | Usually used through operator '+' such as in ct1 + ct2.

* | op2: CountTerm

.. _CountTerm____sub__:

__sub__ (op2: CountTerm) -> CountTerm
-------------------------------------


  | Create a new CountTerm that represents subtraction of two count terms.
  | Usually used through operator '-' such as in ct1 - ct2.

* | op2: CountTerm


Observables
===========

Container used to hold observables-related model data. 
Observables are the measured values of the system. 
This class also includes information on visualization of simulation.

Example: `2600_get_current_mol_count/observables.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/2600_get_current_mol_count/observables.py>`_ 

Attributes:
***********
.. _Observables__viz_outputs:

viz_outputs: List[VizOutput]
----------------------------

  | List of visualization outputs to be included in the model.
  | There is usually just one VizOutput object.
  | - default argument value in constructor: None

.. _Observables__counts:

counts: List[Count]
-------------------

  | List of counts to be included in the model.
  | - default argument value in constructor: None


Methods:
*********
.. _Observables__add_viz_output:

add_viz_output (viz_output: VizOutput)
--------------------------------------


  | Adds a reference to the viz_output object to the list of visualization output specifications.

* | viz_output: VizOutput

.. _Observables__add_count:

add_count (count: Count)
------------------------


  | Adds a reference to the count object to the list of count specifications.

* | count: Count

.. _Observables__find_count:

find_count (name: str) -> Count
-------------------------------


  | Finds a count object by its name, returns None if no such count is present.

* | name: str

.. _Observables__load_bngl_observables:

load_bngl_observables (file_name: str, observables_path_or_file: str=None, parameter_overrides: Dict[str, float]=None, observables_output_format: CountOutputFormat=CountOutputFormat.AUTOMATIC_FROM_EXTENSION)
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.

* | file_name: str
  | Path to the BNGL file.

* | observables_path_or_file: str = None
  | Directory prefix or file name where observable values will be stored.
  | If a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty 
  | string/unset is used, each observable gets its own file and the output file format for created Count 
  | objects is CountOutputFormat.DAT.
  | When not set, this path is used: './react_data/seed_' + str(model.config.seed).zfill(5) + '/'.
  | If a file has a .gdat extension such as 
  | './react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat', all observable are stored in this 
  | file and the output file format for created Count objects is CountOutputFormat.GDAT.
  | Must not be empty when observables_output_format is explicitly set to CountOutputFormat.GDAT.

* | parameter_overrides: Dict[str, float] = None
  | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
  | its value is ignored and instead value parameter_overrides[k] is used.

* | observables_output_format: CountOutputFormat = CountOutputFormat.AUTOMATIC_FROM_EXTENSION
  | Selection of output format. Default setting uses automatic detection
  | based on contents of the 'observables_path_or_file' attribute.

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 



VizOutput
=========

Defines a visualization output with locations of molecules 
that can be then loaded by CellBlender.

Example: `1100_point_release/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1100_point_release/model.py>`_ 

Attributes:
***********
.. _VizOutput__output_files_prefix:

output_files_prefix: str
------------------------

  | Prefix for the viz output files.
  | When not set, the default prefix value is computed from the simulation seed
  | when the model is initialized to\: 
  | './viz_data/seed_' + str(seed).zfill(5) + '/Scene'.
  | - default argument value in constructor: None

.. _VizOutput__species_list:

species_list: List[Species]
---------------------------

  | Specifies a list of species to be visualized, when empty, all_species will be generated.
  | - default argument value in constructor: None

.. _VizOutput__mode:

mode: VizMode
-------------

  | Specified the output format of the visualization files. 
  | VizMode.ASCII is a readable representation, VizMode.CELLBLENDER is a binary representation 
  | that cannot be read using a text editor but is faster to generate.
  | - default argument value in constructor: VizMode.ASCII

.. _VizOutput__every_n_timesteps:

every_n_timesteps: float
------------------------

  | Specifies periodicity of visualization output.
  | Value is truncated (floored) to an integer.
  | Value 0 means that the viz output is ran only once at iteration 0.
  | - default argument value in constructor: 1

