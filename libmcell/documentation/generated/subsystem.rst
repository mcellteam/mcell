.. _api-subsystem:

*********
Subsystem
*********
Complex
=======

This class represents a complex molecule composed of molecule instances.
It is either defined using a BNGL string or using a list of elementary molecule instances.
On top of that, orientation may be defined.
This class is most often by calling its constructor as m.Complex(bngl_string) in cases where a 
fully qualified instance (such as for molecule releases) or a pattern (in observable counts) is needed.  
Comparison operator __eq__ first converts complexes to their canonical representation and 
then does comparison so for instance m.Complex('A(b!1).B(a!1)') == m.Complex('B(a!2).A(b!2)').

Example: `0040_to_bngl_str/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/0040_to_bngl_str/model.py>`_ 

Attributes:
***********
.. _Complex__name:

name: str
---------

  | When set, this complex instance is initialized from a BNGL string passed as this argument, 
  | the string is parsed and elementary_molecules and compartment are initialized.
  | Only one of name or elementary_molecules can be set.
  | - default argument value in constructor: None

.. _Complex__elementary_molecules:

elementary_molecules: List[ElementaryMolecule]
----------------------------------------------

  | Individual molecule instances contained in the complex.
  | Only one of name or elementary_molecules can be set.
  | - default argument value in constructor: None

.. _Complex__orientation:

orientation: Orientation
------------------------

  | Specifies orientation of a molecule. 
  | When Orientation.DEFAULT if kept then during model initialization is
  | 'orientation' set to Orientation.NONE for volume complexes and to 
  | Orientation.UP for surface complexes.
  | Ignored by derived class Species.
  | - default argument value in constructor: Orientation.DEFAULT

.. _Complex__compartment_name:

compartment_name: str
---------------------

  | Specifies compartment name of this Complex.
  | Only one of 'orientation' and 'compartment_name' can be set. 
  | Corresponds to BNGL specification of a compartment for the whole complex '\@COMP\:'.
  | If a 2D/surface compartment is specified, the complex must be a surface complex and 
  | orientation is set to Orientation.UP.
  | If a 3D/volume compartment is specified, the complex must be a volume complex and
  | orientation is set to Orientation.NONE.
  | Sets compartment to all elementary molecules whose compartment is unset. Does not override 
  | specific compartments of elementary molecules that were already set.
  | If this is a volume complex (all elementary molecules have their diffusion_constant_3d set), 
  | all compartments of elementary molecules must be the same volume compartment.
  | If this is a surface complex (at least one elementary molecule has its their diffusion_constant_2d 
  | set), all compartments of surface elementary molecules must be the same, and
  | all compartments of volume elementary molecules must be from the two neighboring 
  | volume compartments.
  | - default argument value in constructor: None


Methods:
*********
.. _Complex__to_bngl_str:

to_bngl_str () -> str
---------------------


  | Creates a string that corresponds to its BNGL representation including compartments.


.. _Complex__as_species:

as_species () -> Species
------------------------


  | Returns a Species object based on this Complex. All species-specific 
  | attributes are set to their default values and 'name' is set to value returned by 
  | 'to_bngl_str()'.



Component
=========

Instance of a component type belonging to a molecule instance.
A component instance must have its state set if there is at least one allowed state.
It is also used to connect molecule instance in a complex instance through bonds.

Example: `0040_to_bngl_str/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/0040_to_bngl_str/model.py>`_ 

Attributes:
***********
.. _Component__component_type:

component_type: ComponentType
-----------------------------

  | Reference to a component type.


.. _Component__state:

state: str
----------

  | Specific state value of this component instance.
  | - default argument value in constructor: STATE_UNSET

.. _Component__bond:

bond: int
---------

  | Specific bond for this component instance.
  | It is either a numberical value such as in A(c!1),
  | or one of special values BOND_UNBOUND in A(c), 
  | BOND_BOUND in A(c!+) or BOND_ANY in A(c!?).
  | - default argument value in constructor: BOND_UNBOUND


Methods:
*********
.. _Component__to_bngl_str:

to_bngl_str () -> str
---------------------


  | Creates a string that corresponds to this component's BNGL representation.



ComponentType
=============

Multiple functional attributes for each molecule type are described using components. And this class defines a type of a component. For example, proteins have multiple functional substructures such as domains, motifs, and binding sites. These components can be unchanging (called stateless) or exist in one of many different internal states For example, certain binding motifs may have different behaviors depending on whether they are unphosphorylated or phosphorylated.

Example: `0040_to_bngl_str/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/0040_to_bngl_str/model.py>`_ 

Attributes:
***********
.. _ComponentType__name:

name: str
---------

  | Name of this component type.


.. _ComponentType__states:

states: List[str]
-----------------

  | List of states allowed by this component.
  | - default argument value in constructor: None


Methods:
*********
.. _ComponentType__inst:

inst (state: str=STATE_UNSET, bond: int=BOND_UNBOUND) -> Component
------------------------------------------------------------------


  | Instantiate a component from this component type.

* | state: str = STATE_UNSET
  | Selected state, must be from the list of the allowed states.

* | bond: int = BOND_UNBOUND
  | Bond information for the created component instance.


.. _ComponentType__inst:

inst (state: int=STATE_UNSET_INT, bond: int=BOND_UNBOUND) -> Component
----------------------------------------------------------------------


  | Instantiate a component from this component type.

* | state: int = STATE_UNSET_INT
  | Selected state, must be from the list of the allowed, converted to string.

* | bond: int = BOND_UNBOUND
  | Bond information for the created component instance.


.. _ComponentType__to_bngl_str:

to_bngl_str () -> str
---------------------


  | Creates a string that corresponds to its BNGL representation.



ElementaryMolecule
==================

Instance of an elementary molecule type. A BNGL complex is composed of elementary molecules.

Example: `0040_to_bngl_str/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/0040_to_bngl_str/model.py>`_ 

Attributes:
***********
.. _ElementaryMolecule__elementary_molecule_type:

elementary_molecule_type: ElementaryMoleculeType
------------------------------------------------

  | Reference to the type of this elementary molecule.


.. _ElementaryMolecule__components:

components: List[Component]
---------------------------

  | List of component instances. Not all components need to be specified 
  | in case when this elementary molecule is used in a pattern.
  | - default argument value in constructor: None

.. _ElementaryMolecule__compartment_name:

compartment_name: str
---------------------

  | Optional BNGL compartment name for this elemenrary molecule. If a 2D/surface compartment is specified, the elementary moelcule must be of surface type. If a 3D/volume compartment is specified, the elementary moelcule must be of volume type.
  | - default argument value in constructor: None


Methods:
*********
.. _ElementaryMolecule__to_bngl_str:

to_bngl_str (with_compartment: bool=True) -> str
------------------------------------------------


  | Creates a string that corresponds to its BNGL representation

* | with_compartment: bool = True
  | Include compartment name in the returned BNGL string.



ElementaryMoleculeType
======================

An elementary molecule type is a base indivisible entity. It is the same as  
a molecule type in BNGL entered in section molecule types. 
The 'elementary' prefix was added to distinguish it clearly from molecules in 
simulation.

Example: `0040_to_bngl_str/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/0040_to_bngl_str/model.py>`_ 

Attributes:
***********
.. _ElementaryMoleculeType__name:

name: str
---------

  | Name of this elementary molecule type.


.. _ElementaryMoleculeType__components:

components: List[ComponentType]
-------------------------------

  | List of components used by this elementary molecule type.
  | - default argument value in constructor: None

.. _ElementaryMoleculeType__diffusion_constant_2d:

diffusion_constant_2d: float
----------------------------

  | Elementary molecule based on this type is constrained to a surface
  | and diffuses with the specified diffusion constant.
  | D can be zero, in which case the molecule doesn’t move. 
  | The units of D are cm^2 /s.
  | - default argument value in constructor: None

.. _ElementaryMoleculeType__diffusion_constant_3d:

diffusion_constant_3d: float
----------------------------

  | Elementary molecule based on this type diffuses in space with the 
  | specified diffusion constant D. 
  | D can be zero, in which case the molecule doesn’t move. 
  | The units of D are cm^2 /s.
  | - default argument value in constructor: None

.. _ElementaryMoleculeType__custom_time_step:

custom_time_step: float
-----------------------

  | This molecule should take timesteps of length custom_time_step (in seconds). 
  | Use either this or custom_time_step, not both.
  | - default argument value in constructor: None

.. _ElementaryMoleculeType__custom_space_step:

custom_space_step: float
------------------------

  | This molecule should take steps of average length given by the custom_space_step value (in microns). 
  | Use either this or custom_time_step, not both.
  | - default argument value in constructor: None

.. _ElementaryMoleculeType__target_only:

target_only: bool
-----------------

  | This molecule will not initiate reactions when it runs into other molecules. This
  | setting can speed up simulations when applied to a molecule at high concentrations 
  | that reacts with a molecule at low concentrations (it is more efficient for
  | the low-concentration molecule to trigger the reactions). This directive does
  | not affect unimolecular reactions.
  | - default argument value in constructor: False


Methods:
*********
.. _ElementaryMoleculeType__inst:

inst (components: List[Component]=None, compartment_name: str=None) -> ElementaryMolecule
-----------------------------------------------------------------------------------------


  | Create an elementary molecule based on this elementary molecule type.

* | components: List[Component] = None
  | Instances of components for the the created elementary molecule.
  | Not all components need to be specified in case when the elementary 
  | molecule is used in a pattern.

* | compartment_name: str = None
  | Optional specification of compartment name for the created elementary molecule.


.. _ElementaryMoleculeType__to_bngl_str:

to_bngl_str () -> str
---------------------


  | Creates a string that corresponds to its BNGL representation.



ReactionRule
============

Represents a BioNetGen Language (BNGL) reaction rule. 
In BNGL, a reaction is simply one or more transformations
applied simultaneously to one or more species. The following
transformations (and their combinations) are allowed\:
  * Forming a bond, e.g. A(b) + B(a) -> A(b!0).B(a!0)
  * Breaking a bond, e.g. A(b!0).B(a!0)-> A(b)+ B(a)
  * Changing of component state, e.g. X(y~0) -> X(y~p)
  * Creating a molecule, e.g. A(b) -> A(b) + C(d)
  * Destroying a molecule, e.g. A(b) + B(a) -> A(b) or A -> Null 
    (Null here means that there is no product)
  * Changing species of a bound molecule when the molecule type has the 
    same components, e.g. A(b!0).B(a!0)-> A(b!0).C(a!0)
  
Also compartments may be specified in reactants (patterns) and for products.
Special compartment classes supported by MCell4 are @IN and @OUT.
They can be used in surface reactions to constrain a reaction with a volume molecule 
hitting a surface molecule from the inside or outside of the compartment, 
e.g. A(s)@IN + S(a) -> S(a!1).A(s!1) and/or to define the location of the 
product, e.g. S(a!1).A(s!1) -> S(a) + A(s)@OUT.

Examples: `0040_to_bngl_str/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/0040_to_bngl_str/model.py>`_ `1840_vol_plus_surf_class_rxn_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1840_vol_plus_surf_class_rxn_callback/model.py>`_ 

Attributes:
***********
.. _ReactionRule__name:

name: str
---------

  | Name of the reaction. If this is a reversible reaction, then it is the name of the 
  | reaction in forward direction.
  | - default argument value in constructor: None

.. _ReactionRule__reactants:

reactants: List[Complex]
------------------------

  | List of reactant patterns. Must contain one or two patterns.
  | - default argument value in constructor: None

.. _ReactionRule__products:

products: List[Complex]
-----------------------

  | List of reactant patterns. May be empty.
  | - default argument value in constructor: None

.. _ReactionRule__fwd_rate:

fwd_rate: float
---------------

  | The units of the reaction rate for uni- and bimolecular reactions are\:
  |   \* [s^-1] for unimolecular reactions,
  |   \* [N^-1\*s^-1] bimolecular reactions between two surface molecules on different objects 
  |     (this is a highly experimental feature and the unit will likely change in the future, 
  |      not sure if probability is computed correctly, it works the way that the surface molecule 
  |      is first diffused and then a potential collisions within the distance of Config.intermembrane_interaction_radius
  |      are evaluated). 
  | Other bimolecular reaction units depend on Model.config.use_bng_units settings.
  | When use_bng_units is False (default), traditional MCell units are used:  
  |   \* [M^-1\*s^-1] for bimolecular reactions between either two volume molecules, a volume molecule 
  |                 and a surface (molecule), 
  |   \* [um^2\*N^-1\*s^-1] bimolecular reactions between two surface molecules on the same surface, and
  | When use_bng_units is True, units compatible with BioNetGen's ODE, SSA, and PLA solvers are used:
  |   \* [um^3\*N^-1\*s^-1] for any bimolecular reactions, surface-surface reaction rate conversion assumes 10nm membrane thickness. 
  | M is the molarity of the solution and N the number of reactants.
  | May be changed after model initialization. 
  | Setting of value is ignored if the rate does not change. 
  | If the new value differs from previous, updates all information related 
  | to the new rate including recomputation of reaction times for molecules if this is a
  | unimolecular reaction.
  | - default argument value in constructor: None

  | Examples: `2500_rxn_rate_change_bimol_box_it_100/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2500_rxn_rate_change_bimol_box_it_100/model.py>`_ `2700_concentration_based_rxn_rate/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2700_concentration_based_rxn_rate/model.py>`_ 


.. _ReactionRule__rev_name:

rev_name: str
-------------

  | Name of the reaction in reverse direction.
  | - default argument value in constructor: None

.. _ReactionRule__rev_rate:

rev_rate: float
---------------

  | Reverse reactions rate, reaction is unidirectional when not specified.
  | May be changed after model initialization, in the case behaves the same was as for 
  | changing the 'fwd_rate'. 
  | Uses the same units as 'fwd_rate'.
  | - default argument value in constructor: None

.. _ReactionRule__variable_rate:

variable_rate: List[List[float]]
--------------------------------

  | The array passed as this argument must have as its items a pair of floats (time in s, rate).
  | Must be sorted by time (this is not checked).      
  | Variable rate is applicable only for irreversible reactions.
  | When simulation starts and the table does not contain value for time 0, the initial fwd_rate is set to 0.
  | When time advances after the last time in this table, the last rate is used for all subsequent iterations.   
  | Members fwd_rate and rev_rate must not be set when setting this attribute through a constructor. 
  | When this attribute is set outside of the class constructor, fwd_rate is automatically reset to an 'unset' value.
  | Cannot be set after model initialization.
  | - default argument value in constructor: None

.. _ReactionRule__is_intermembrane_surface_reaction:

is_intermembrane_surface_reaction: bool
---------------------------------------

  | Experimental, see addintinal explanation in 'fwd' rate.
  | Then set to true, this is a special type of surface-surface reaction that 
  | allows for two surface molecules to react when they are on different geometrical objects. 
  | Note\: This support is limited for now, the reaction rule must be in the form of A + B -> C + D 
  | where all reactants and products must be surface molecules and their orientation must be 'any' (default).
  | - default argument value in constructor: False

  | Example: `3000_intermembrane_rxns/customization.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/3000_intermembrane_rxns/customization.py>`_ 



Methods:
*********
.. _ReactionRule__to_bngl_str:

to_bngl_str () -> str
---------------------


  | Creates a string that corresponds to the reaction rule's BNGL representation, does not contain rates.



Species
=======

There are three ways how to use this class\:
1) definition of simple species - in this case 'name' is 
a single identifier and at least 'diffusion_constant_2d' or 
'diffusion_constant_3d' must be provided.
Example\: m.Species('A', diffusion_constant_3d=1e-6). 
Such a definition must be added to subsystem or model so that  
during model initialization this species is transformed to MCell 
representation and an ElementaryMoleculeType 'A' with a given 
diffusion constant is created as well.
2) full definition of complex species - in this case the 
inherited attribute 'elementary_molecules' from Complex
is used as a definition of the complex and this gives information 
on diffusion constants of the used elementary molecules.
Example\: m.Species(elementary_molecules=[ei1, ei2]). 
Such a definition must be added to subsystem or model.   
3) declaration of species - in this case only 'name' in the form of 
an BNGL string is provided. The complex instance specified by the name 
must be fully qualified (i.e. all components are present and those 
components that have a state have their state set).
No information on diffusion constants and other properties of 
used elementary molecules is provided, it must be provided elsewhere.
Example\: m.Species('A(b!1).B(a!1)').
This is a common form of usage when reaction rules are provided in a BNGL file.
Such declaration does no need to be added to subsystem or model.
This form is used as argument in cases where a fully qualified instance  
must be provided such as in molecule releases.

Example: `0040_to_bngl_str/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/0040_to_bngl_str/model.py>`_ 

Attributes:
***********
.. _Species__name:

name: str
---------

  | Name of the species in the BNGL format. 
  | One must either specify name or elementary_molecules (inherited from Complex). 
  | This argument name is parsed during model initialization.
  | - default argument value in constructor: None

.. _Species__diffusion_constant_2d:

diffusion_constant_2d: float
----------------------------

  | This molecule is constrained to surface  with diffusion constant D. 
  | D can be zero, in which case the molecule doesn’t move. 
  | The units of D are cm^2/s.
  | - default argument value in constructor: None

.. _Species__diffusion_constant_3d:

diffusion_constant_3d: float
----------------------------

  | This molecule diffuses in space with diffusion constant D. 
  | D can be zero, in which case the molecule doesn’t move. 
  | The units of D are cm^2/s.
  | - default argument value in constructor: None

.. _Species__custom_time_step:

custom_time_step: float
-----------------------

  | Optional setting of a custom time step for this specific species. 
  | A molecule of this species should take timesteps of length custom_time_step (in seconds). 
  | Use either this or custom_time_step.
  | - default argument value in constructor: None

.. _Species__custom_space_step:

custom_space_step: float
------------------------

  | Optional setting of a custom space step for this specific species. 
  | A molecule of this species should take steps of average length custom_space_step (in microns). 
  | Use either this or custom_time_step.
  | - default argument value in constructor: None

.. _Species__target_only:

target_only: bool
-----------------

  | A molecule of this species will not initiate reactions when it runs into other molecules. This
  | setting can speed up simulations when applied to a molecule at high concentrations 
  | that reacts with a molecule at low concentrations (it is more efficient for
  | the low-concentration molecule to trigger the reactions). This directive does
  | not affect unimolecular reactions.
  | - default argument value in constructor: False

.. _Species__name:

name: str
---------

  | When set, this complex instance is initialized from a BNGL string passed as this argument, 
  | the string is parsed and elementary_molecules and compartment are initialized.
  | Only one of name or elementary_molecules can be set.
  | - default argument value in constructor: None

.. _Species__elementary_molecules:

elementary_molecules: List[ElementaryMolecule]
----------------------------------------------

  | Individual molecule instances contained in the complex.
  | Only one of name or elementary_molecules can be set.
  | - default argument value in constructor: None

.. _Species__orientation:

orientation: Orientation
------------------------

  | Specifies orientation of a molecule. 
  | When Orientation.DEFAULT if kept then during model initialization is
  | 'orientation' set to Orientation.NONE for volume complexes and to 
  | Orientation.UP for surface complexes.
  | Ignored by derived class Species.
  | - default argument value in constructor: Orientation.DEFAULT

.. _Species__compartment_name:

compartment_name: str
---------------------

  | Specifies compartment name of this Complex.
  | Only one of 'orientation' and 'compartment_name' can be set. 
  | Corresponds to BNGL specification of a compartment for the whole complex '\@COMP\:'.
  | If a 2D/surface compartment is specified, the complex must be a surface complex and 
  | orientation is set to Orientation.UP.
  | If a 3D/volume compartment is specified, the complex must be a volume complex and
  | orientation is set to Orientation.NONE.
  | Sets compartment to all elementary molecules whose compartment is unset. Does not override 
  | specific compartments of elementary molecules that were already set.
  | If this is a volume complex (all elementary molecules have their diffusion_constant_3d set), 
  | all compartments of elementary molecules must be the same volume compartment.
  | If this is a surface complex (at least one elementary molecule has its their diffusion_constant_2d 
  | set), all compartments of surface elementary molecules must be the same, and
  | all compartments of volume elementary molecules must be from the two neighboring 
  | volume compartments.
  | - default argument value in constructor: None


Methods:
*********
.. _Species__inst:

inst (orientation: Orientation=Orientation.DEFAULT, compartment_name: str=None) -> Complex
------------------------------------------------------------------------------------------


  | Creates a copy of a Complex from this Species with specified orientation and compartment name.

* | orientation: Orientation = Orientation.DEFAULT
  | Maximum one of orientation or compartment_name can be set, not both.

* | compartment_name: str = None
  | Maximum one of orientation or compartment_name can be set, not both.


.. _Species__to_bngl_str:

to_bngl_str () -> str
---------------------


  | Creates a string that corresponds to its BNGL representation including compartments.


.. _Species__as_species:

as_species () -> Species
------------------------


  | Returns a Species object based on this Complex. All species-specific 
  | attributes are set to their default values and 'name' is set to value returned by 
  | 'to_bngl_str()'.



Subsystem
=========

Subsystem usually defines a reaction network. It is a collection of 
species and reaction rules that use these species. 
The main motivation for introducing such an object to MCell4 is to have 
a class independent on that particular initial model state and observables that 
only contains reactions. This way, one can define independent reusable subsystems
and possibly merge them together when creating a model that includes multiple reaction 
networks.

Example: `2550_variable_rate_unimol_w_rxn_class_cleanup/sybsystem.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2550_variable_rate_unimol_w_rxn_class_cleanup/sybsystem.py>`_ 

Attributes:
***********
.. _Subsystem__species:

species: List[Species]
----------------------

  | List of species to be included in the model for initialization.
  | Used usually only for simple species (species that are defined using a
  | single molecule type without components such as 'A').
  | Other species may be created inside simulation
  | - default argument value in constructor: None

.. _Subsystem__reaction_rules:

reaction_rules: List[ReactionRule]
----------------------------------

  | - default argument value in constructor: None

.. _Subsystem__surface_classes:

surface_classes: List[SurfaceClass]
-----------------------------------

  | - default argument value in constructor: None

.. _Subsystem__elementary_molecule_types:

elementary_molecule_types: List[ElementaryMoleculeType]
-------------------------------------------------------

  | Contains list of elementary molecule types with their diffusion constants and other information. 
  | Populated when a BNGL file is loaded and also on initialization from Species objects present in 
  | the species list.
  | - default argument value in constructor: None


Methods:
*********
.. _Subsystem__add_species:

add_species (s: Species)
------------------------


  | Add a reference to a Species object to the species list.

* | s: Species

.. _Subsystem__find_species:

find_species (name: str) -> Species
-----------------------------------


  | Find a Species object using name in the species list. 
  | Returns None if no such species is found.

* | name: str

.. _Subsystem__add_reaction_rule:

add_reaction_rule (r: ReactionRule)
-----------------------------------


  | Add a reference to a ReactionRule object to the reaction_rules list.

* | r: ReactionRule

.. _Subsystem__find_reaction_rule:

find_reaction_rule (name: str) -> ReactionRule
----------------------------------------------


  | Find a ReactionRule object using name in the reaction_rules list. 
  | Returns None if no such reaction rule is found.

* | name: str

.. _Subsystem__add_surface_class:

add_surface_class (sc: SurfaceClass)
------------------------------------


  | Add a reference to a SurfaceClass object to the surface_classes list.

* | sc: SurfaceClass

.. _Subsystem__find_surface_class:

find_surface_class (name: str) -> SurfaceClass
----------------------------------------------


  | Find a SurfaceClass object using name in the surface_classes list. 
  | Returns None if no such surface class is found.

* | name: str

.. _Subsystem__add_elementary_molecule_type:

add_elementary_molecule_type (mt: ElementaryMoleculeType)
---------------------------------------------------------


  | Add a reference to an ElementaryMoleculeType object to the elementary_molecule_types list.

* | mt: ElementaryMoleculeType

.. _Subsystem__find_elementary_molecule_type:

find_elementary_molecule_type (name: str) -> ElementaryMoleculeType
-------------------------------------------------------------------


  | Find an ElementaryMoleculeType object using name in the elementary_molecule_types list. 
  | Returns None if no such elementary molecule type is found.

* | name: str

.. _Subsystem__load_bngl_molecule_types_and_reaction_rules:

load_bngl_molecule_types_and_reaction_rules (file_name: str, parameter_overrides: Dict[str, float]=None)
--------------------------------------------------------------------------------------------------------


  | Parses a BNGL file, only reads molecule types and reaction rules sections, 
  | i.e. ignores observables and seed species. 
  | Parameter values are evaluated and the result value is directly used.  
  | Compartments names are stored in rxn rules as strings because compartments belong 
  | to geometry objects and the subsystem is independent on specific geometry.
  | However, the compartments and their objects must be defined before initialization.

* | file_name: str
  | Path to the BNGL file to be loaded.

* | parameter_overrides: Dict[str, float] = None
  | For each key k in the parameter_overrides, if it is defined in the BNGL's parameters section,
  | its value is ignored and instead value parameter_overrides[k] is used.

  | Example: `2100_gradual_bngl_load/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/2100_gradual_bngl_load/model.py>`_ 



SurfaceClass
============

Defining a surface class allows surfaces to behave like species. For instance, one may wish 
to specify that a surface does not block the diffusion of molecules. Each type of surface is defined
by name, and each surface name must be unique in the simulation and should not match any molecule names.
To define a reaction with a surface class, use constructor call m.Complex(name) as one of the reactants.

Examples: `1600_crossing_transparent_compartment_wall/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4/1600_crossing_transparent_compartment_wall/model.py>`_ `1840_vol_plus_surf_class_rxn_callback/model.py <https://github.com/mcellteam/mcell_tests/blob/master/tests/pymcell4_positive/1840_vol_plus_surf_class_rxn_callback/model.py>`_ 

Attributes:
***********
.. _SurfaceClass__name:

name: str
---------

  | Name of the surface class.


.. _SurfaceClass__properties:

properties: List[SurfaceProperty]
---------------------------------

  | A surface class can either have a list of properties or just one property.
  | In the usual case of having one property, one can set the attributes 
  | type, affected_species, etc. inherited from SurfaceProperty directly.
  | - default argument value in constructor: None

.. _SurfaceClass__type:

type: SurfacePropertyType
-------------------------

  | Must be set. See SurfacePropertyType for options.
  | - default argument value in constructor: SurfacePropertyType.UNSET

.. _SurfaceClass__affected_complex_pattern:

affected_complex_pattern: Complex
---------------------------------

  | A complex pattern with optional orientation must be set.
  | Default orientation means that the pattern matches any orientation.
  | For concentration or flux clamp the orientation specifies on which side  
  | will be the concentration held (UP is front or outside, DOWN is back or 
  | inside, and DEFAULT, ANY or NONE is on both sides).
  | The complex pattern must not use compartments.
  | - default argument value in constructor: None

.. _SurfaceClass__concentration:

concentration: float
--------------------

  | Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or 
  | SurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side 
  | of the wall that has this concentration or flux clamped.
  | - default argument value in constructor: None

SurfaceProperty
===============

Single property for a SurfaceClass.

Attributes:
***********
.. _SurfaceProperty__type:

type: SurfacePropertyType
-------------------------

  | Must be set. See SurfacePropertyType for options.
  | - default argument value in constructor: SurfacePropertyType.UNSET

.. _SurfaceProperty__affected_complex_pattern:

affected_complex_pattern: Complex
---------------------------------

  | A complex pattern with optional orientation must be set.
  | Default orientation means that the pattern matches any orientation.
  | For concentration or flux clamp the orientation specifies on which side  
  | will be the concentration held (UP is front or outside, DOWN is back or 
  | inside, and DEFAULT, ANY or NONE is on both sides).
  | The complex pattern must not use compartments.
  | - default argument value in constructor: None

.. _SurfaceProperty__concentration:

concentration: float
--------------------

  | Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or 
  | SurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side 
  | of the wall that has this concentration or flux clamped.
  | - default argument value in constructor: None

