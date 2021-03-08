***********
Observables
***********
Count
=====

Attributes:
***********
* | **name**: str = None
  | Name of a count may be specified when one needs to search for them later. 
  | Also when the count is created while loading a BNGL file, its name is set.

* | **file_name**: str = None
  | File name with an optional path must be set. It is not deduced automatically.

* | **count_expression**: CountTerm = None
  | The count expression must be composed only from CountTerm objects that are added or 
  | subtracted.

* | **multiplier**: float = 1
  | In some cases it might be useful to multiply the whole count by a constant to get 
  | for instance concentration. The count_expression is not an arbitrary expression
  | and such multiplication can be done through this attribute.

* | **every_n_timesteps**: float = 1
  | Value is truncated (floored) to an integer.
  | If value is set to 0, this Count is used only on-demand through calls to its
  | get_current_value method.

* | **species_pattern**: Complex = None
  | Count the number of molecules that match the given complex instance pattern.
  | Counts each molecule exactly once. 
  | If the pattern has a compartment set, this specifies the counted region.

* | **molecules_pattern**: Complex = None
  | Count the number of matches of the given pattern on molecules.
  | The observable will count a molecule every time it matches the pattern.
  | When the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a 
  | molecule A(a!1).A(a!1,b!2).B(a!2) will be counted twice because the 
  | pattern may match in two different ways. 
  | If the pattern has a compartment set, this specifies the counted region.

* | **reaction_rule**: ReactionRule = None
  | Count the number of reactions that occurred since the start of the simulation.

* | **region**: Region = None
  | Only a GeometryObject or SurfaceRegion can be passed as the region argument, 
  | compound regions (created with +, -, \*) are not supproted yet.   
  | Cannot be set when 'species_pattern' or 'molecules_pattern' has a  
  | compartment specified.
  | If pattern compartment is not specified and 'region' is left 'unset', 
  | counting is done in the whole world.

* | **node_type**: ExprNodeType = ExprNodeType.LEAF
  | Internal, used to represent an expression

* | **left_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the right operand

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


* | **__sub__**

   * | op2: CountTerm
   * | return type: CountTerm



Observables
===========

Neither VizOutput, nor Count have name, therefore there are no find_* methods.

Attributes:
***********
* | **viz_outputs**: List[VizOutput] = None

* | **counts**: List[Count] = None


Methods:
*********
* | **add_viz_output**

   * | viz_output: VizOutput

* | **add_count**

   * | count: Count

* | **find_count**

   * | name: str
   * | return type: Count


* | **load_bngl_observables**

   * | file_name: str
     | BNGL file name.

   * | output_files_prefix: str = ''
     | Prefix to be used when creating files with observable values.

   * | parameter_overrides: Dict[str, float] = None

  | Loads section observables from a BNGL file and creates Count objects according to it.
  | All elementary molecule types used in the seed species section must be defined in subsystem.



VizOutput
=========

Attributes:
***********
* | **output_files_prefix**: str

* | **species_list**: List[Species] = None
  | Specifies a list of species to be visualized, when empty, all_species will be generated.

* | **mode**: VizMode = VizMode.ASCII

* | **every_n_timesteps**: float = 1
  | Value is truncated (floored) to an integer.
  | Value 0 means that the viz output is ran only once at iteration 0.

CountTerm
=========

Attributes:
***********
* | **species_pattern**: Complex = None
  | Count the number of molecules that match the given complex instance pattern.
  | Counts each molecule exactly once. 
  | If the pattern has a compartment set, this specifies the counted region.

* | **molecules_pattern**: Complex = None
  | Count the number of matches of the given pattern on molecules.
  | The observable will count a molecule every time it matches the pattern.
  | When the pattern is symmetric, e.g. as in A(a!1).A(a!1) then a 
  | molecule A(a!1).A(a!1,b!2).B(a!2) will be counted twice because the 
  | pattern may match in two different ways. 
  | If the pattern has a compartment set, this specifies the counted region.

* | **reaction_rule**: ReactionRule = None
  | Count the number of reactions that occurred since the start of the simulation.

* | **region**: Region = None
  | Only a GeometryObject or SurfaceRegion can be passed as the region argument, 
  | compound regions (created with +, -, \*) are not supproted yet.   
  | Cannot be set when 'species_pattern' or 'molecules_pattern' has a  
  | compartment specified.
  | If pattern compartment is not specified and 'region' is left 'unset', 
  | counting is done in the whole world.

* | **node_type**: ExprNodeType = ExprNodeType.LEAF
  | Internal, used to represent an expression

* | **left_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the left operand

* | **right_node**: CountTerm = None
  | Internal, when node_type is not Leaf, this is the right operand

* | **initial_reactions_count**: int = 0
  | Used for checkpointing, allows to set initial count of reactions that occurred.
  | Ignored when molecules are counted.


Methods:
*********
* | **__add__**

   * | op2: CountTerm
   * | return type: CountTerm


* | **__sub__**

   * | op2: CountTerm
   * | return type: CountTerm



