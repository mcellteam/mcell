***********
Observables
***********
Count
=====

Represents a molecule or reaction count observable. 
Inherits all members from class CountTerm. 
Usually just one observable is counted and a member of this class count_expression 
is not set. If an expression is needed, it is stored in the count_expression
attribute.

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
  | The usual file_name, while using this count's name is\:
  | file_name = './react_data/seed_' + str(SEED).zfill(5) + '/' + name + '.dat'.

* | **count_expression**: CountTerm = None
  | If this count uses more than one count term one must define a count term expression 
  | and set a reference to it to this attribute.  
  | The count expression must be composed only from CountTerm objects that are added or 
  | subtracted.

* | **multiplier**: float = 1
  | In some cases it might be useful to multiply the whole count by a constant to get 
  | for instance concentration. The count_expression allows only addition and subtraction 
  | of count terms so such multiplication can be done through this attribute.
  | It can be also used to divide the resulting count by passing an inverse of the divisor (1/d).

* | **every_n_timesteps**: float = 1
  | Specifies periodicity of this count's output.
  | Value is truncated (floored) to an integer.
  | If value is set to 0, this Count is used only on-demand through calls to its
  | get_current_value method.

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
* | **get_current_value**

   * | return type: float


  | Returns the current value for this count. Cannot be used to count reactions.
  | The model must be initialized with this Count present as one of the observables.


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

   * | output_files_prefix: str = ''
     | Prefix to be used when creating files with observable values.
     | The usual value is './react_data/seed_' + str(SEED).zfill(5) + '/'.

   * | parameter_overrides: Dict[str, float] = None
     | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
     | its value is ignored and instead value parameter_overrides[k] is used.


  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.



VizOutput
=========

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

