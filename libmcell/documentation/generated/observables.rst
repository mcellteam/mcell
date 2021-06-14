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

Example: `1500_region_release_must_set_compartment/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/1500_region_release_must_set_compartment/model.py>`_ 

Attributes:
***********
* | **name**: str = None
  | Name of a count may be specified when one needs to search for them later. 
  | When the count is created when a BNGL file is loaded, its name is set, for instance
  | when the following BNGL code is loaded\:
  | 
  | begin observables
  |    Molecules Acount A
  | end observables
  | 
  | the name is set to Acount.

* | **file_name**: str = None
  | File name where this observable values will be stored.
  | File extension or setting explicit output_format determines the output format.
  | A) When the file_name extension is .dat such as here:
  | file_name = './react_data/seed_' + str(SEED).zfill(5) + '/' + name + '.dat'
  | the output format is set to CountOutputFormat.DAT in the constructor.
  | File names for individual Counts must be different.
  | B) When the file_name extension is .gdat such as here:
  | file_name = './react_data/seed_' + str(SEED).zfill(5) + '/counts.gdat'
  | the output format is set to CountOutputFormat.GDAT in the constructor.
  | The file name is usually the same for all counts but one can 
  | create multiple gdat files with different observables.
  | All observables that are stored into a single .gdat file must have the same 
  | periodicity specified by attribute every_n_timesteps.
  | Must be set.

* | **expression**: CountTerm = None
  | The expression must be set to a root of an expression tree composed of CountTerms. 
  | In the usual cases, there is just one CountTerm in this expression tree and its 
  | node_type is ExprNodeType.LEAF.
  | The count expression tree defines CountTerm objects that are added or subtracted
  | from each other.

* | **multiplier**: float = 1
  | In some cases it might be useful to multiply the whole count by a constant to get 
  | for instance concentration. The expression tree allows only addition and subtraction 
  | of count terms so such multiplication can be done through this attribute.
  | It can be also used to divide the resulting count by passing an inverse of the divisor (1/d).

* | **every_n_timesteps**: float = 1
  | Specifies periodicity of this count's output.
  | Value is truncated (floored) to an integer.
  | If value is set to 0, this Count is used only on-demand through calls to its
  | get_current_value method.

* | **output_format**: CountOutputFormat = CountOutputFormat.AUTOMATIC_FROM_EXTENSION
  | Listed as the last attribute because the automatic default value
  | is sufficient in most cases. 
  | Selection of output format. Default setting uses file extension  
  | from attribute file_name. 
  | When set to CountOutputFormat.AUTOMATIC_FROM_EXTENSION, 
  | this output_format is set automatically only in the Count's constructor.


Methods:
*********
* | **get_current_value**

   * | return type: float


  | Returns the current value for this count. Cannot be used to count reactions.
  | The model must be initialized with this Count present as one of the observables.

  | Example: `2600_get_current_mol_count/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/2600_get_current_mol_count/model.py>`_ 



CountTerm
=========

A count observable can be defined as an expression composed of addition
or subtraction individual count terms. This class represents one count term
in this expression.

Attributes:
***********
* | **species_pattern**: Complex = None
  | Count the number of molecules that match the given complex instance pattern.
  | This corresponds to the BNGL 'Species' specifier in the BNGL seed species section.
  | Counts each molecule exactly once. 
  | If the pattern has a compartment set, this specifies the counted region.
  | Exactly one of species_pattern, molecules_pattern, and reaction_rule must be set.

* | **molecules_pattern**: Complex = None
  | Count the number of matches of the given pattern on molecules.
  | This corresponds to the BNGL 'Molecules' specifier in the BNGL seed species section.
  | The observable will increment the count every time the pattern matches the molecule.
  | For instance, pattern A will match a complex A(a!1).B(a!1,a!2).A(b!2) twice. 
  | When the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a 
  | molecule A(b.a!1).A(a!1,b!2).B(a!2) will be counted twice because the 
  | pattern may match in two different ways. 
  | If the pattern has a compartment set, the compartment is used to filter out the molecules.   
  | Exactly one of species_pattern, molecules_pattern, and reaction_rule must be set.

* | **reaction_rule**: ReactionRule = None
  | Count the number of applications of this specific reactions that occurred since the
  | start of the simulation.
  | Exactly one of species_pattern, molecules_pattern, and reaction_rule must be set.

* | **region**: Region = None
  | Only a GeometryObject or SurfaceRegion can be passed as the region argument, 
  | compound regions (created with +, -, \*) are not supproted yet.   
  | Can be combined with a compartment specified in the species_pattern or molecules_pattern.
  | If compartment in species_pattern or molecules_pattern is not specified and 
  | region is left unset, counting is done in the whole world.

* | **node_type**: ExprNodeType = ExprNodeType.LEAF
  | Internal, used to specify what type of count expression node this object represents.

* | **left_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the left operand.

* | **right_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the right operand.

* | **initial_reactions_count**: int = 0
  | Used for checkpointing, allows to set initial count of reactions that occurred.
  | Ignored when molecules are counted.


Methods:
*********
* | **__add__**

   * | op2: CountTerm
   * | return type: CountTerm


  | Create a new CountTerm that represents addition of two count terms.
  | Usually used through operator '+' such as in ct1 + ct2.


* | **__sub__**

   * | op2: CountTerm
   * | return type: CountTerm


  | Create a new CountTerm that represents subtraction of two count terms.
  | Usually used through operator '-' such as in ct1 - ct2.



Observables
===========

Container used to hold observables-related model data. 
Observables are the measured values of the system. 
This class also includes information on visualization of simulation.

Example: `2600_get_current_mol_count/observables.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4_positive/2600_get_current_mol_count/observables.py>`_ 

Attributes:
***********
* | **viz_outputs**: List[VizOutput] = None
  | List of visualization outputs to be included in the model.
  | There is usually just one VizOutput object.

* | **counts**: List[Count] = None
  | List of counts to be included in the model.


Methods:
*********
* | **add_viz_output**

   * | viz_output: VizOutput

  | Adds a reference to the viz_output object to the list of visualization output specifications.


* | **add_count**

   * | count: Count

  | Adds a reference to the count object to the list of count specifications.


* | **find_count**

   * | name: str
   * | return type: Count


  | Finds a count object by its name, returns None if no such count is present.


* | **load_bngl_observables**

   * | file_name: str
     | Path to the BNGL file.

   * | observables_path_or_file: str = ''
     | Directory prefix or file name where observable values will be stored.
     | If a directory such as './react_data/seed_' + str(SEED).zfill(5) + '/' or an empty 
     | string is used, each observable gets its own file and the output file format for created Count 
     | objects is CountOutputFormat.DAT.
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


  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 



VizOutput
=========

Defines a visualization output with locations of molecules 
that can be then loaded by CellBlender.

Example: `1100_point_release/model.py <https://github.com/mcellteam/mcell_tests/blob/mcell4_dev/tests/pymcell4/1100_point_release/model.py>`_ 

Attributes:
***********
* | **output_files_prefix**: str
  | Prefix for the viz output files, the prefix value is computed from the simulation seed: 
  | output_files_prefix = './viz_data/seed_' + str(SEED).zfill(5) + '/Scene'.

* | **species_list**: List[Species] = None
  | Specifies a list of species to be visualized, when empty, all_species will be generated.

* | **mode**: VizMode = VizMode.ASCII
  | Specified the output format of the visualization files. 
  | VizMode.ASCII is a readable representation, VizMode.CELLBLENDER is a binary representation 
  | that cannot be read using a text editor but is faster to generate.

* | **every_n_timesteps**: float = 1
  | Specifies periodicity of visualization output.
  | Value is truncated (floored) to an integer.
  | Value 0 means that the viz output is ran only once at iteration 0.

