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

Attributes:
***********
* | **name**: str = None
  | When set, this complex instance is initialized from a BNGL string passed as this argument, 
  | the string is parsed and elementary_molecules and compartment are initialized.
  | Only one of name or elementary_molecules can be set.

* | **elementary_molecules**: List[ElementaryMolecule] = None
  | Individual molecule instances contained in the complex.
  | Only one of name or elementary_molecules can be set.

* | **orientation**: Orientation = Orientation.DEFAULT
  | Specifies orientation of a molecule. 
  | When Orientation.DEFAULT if kept then during model initialization is
  | 'orientation' set to Orientation.NONE for volume complexes and to 
  | Orientation.UP for surface complexes.
  | Ignored by derived class Species.

* | **compartment_name**: str = None
  | Specifies compartment name of this Complex.
  | Only one of 'orientation' and 'compartment_name' can be set. 
  | Corresponds to BNGL specification of a compartment for the whole complex '@COMP\:'.
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


Methods:
*********
* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation including compartments.


* | **as_species**

   * | return type: Species


  | Returns a Species object based on this Complex. All species-specific 
  | attributes are set to their default values and 'name' is set to value returned by 
  | 'to_bngl_str()'.



Component
=========

Instance of a component type belonging to a molecule instance.
A component instance must have its state set if there is at least one allowed state.
It is also used to connect molecule instance in a complex instance through bonds.

Attributes:
***********
* | **component_type**: ComponentType
  | Reference to a component type.

* | **state**: str = STATE_UNSET
  | Specific state value of this component instance.

* | **bond**: int = BOND_UNBOUND
  | Specific bond for this component instance.
  | It is either a numberical value such as in A(c!1),
  | or one of special values BOND_UNBOUND in A(c), 
  | BOND_BOUND in A(c!+) or BOND_ANY in A(c!?).


Methods:
*********
* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to this component's BNGL representation.



ComponentType
=============

Multiple functional attributes for each molecule type are described using components. And this class defines a type of a component. For example, proteins have multiple functional substructures such as domains, motifs, and binding sites. These components can be unchanging (called stateless) or exist in one of many different internal states For example, certain binding motifs may have different behaviors depending on whether they are unphosphorylated or phosphorylated.

Attributes:
***********
* | **name**: str
  | Name of this component type.

* | **states**: List[str] = None
  | List of states allowed by this component.


Methods:
*********
* | **inst**

   * | state: str = STATE_UNSET
     | Selected state, must be from the list of the allowed states.

   * | bond: int = BOND_UNBOUND
     | Bond information for the created component instance.

   * | return type: Component


  | Instantiate a component from this component type.


* | **inst**

   * | state: int = STATE_UNSET_INT
     | Selected state, must be from the list of the allowed, converted to string.

   * | bond: int = BOND_UNBOUND
     | Bond information for the created component instance.

   * | return type: Component


  | Instantiate a component from this component type.


* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation.



ElementaryMolecule
==================

Instance of an elementary molecule type. A BNGL complex is composed of elementary molecules.

Attributes:
***********
* | **elementary_molecule_type**: ElementaryMoleculeType
  | Reference to the type of this elementary molecule.

* | **components**: List[Component] = None
  | List of component instances. Not all components need to be specified 
  | in case when this elementary molecule is used in a pattern.

* | **compartment_name**: str = None
  | Optional BNGL compartment name for this elemenrary molecule. If a 2D/surface compartment is specified, the elementary moelcule must be of surface type. If a 3D/volume compartment is specified, the elementary moelcule must be of volume type.


Methods:
*********
* | **to_bngl_str**

   * | with_compartment: bool = True
     | Include compartment name in the returned BNGL string.

   * | return type: str


  | Creates a string that corresponds to its BNGL representation



ElementaryMoleculeType
======================

An elementary molecule type is a base indivisible entity. It is the same as  
a molecule type in BNGL entered in section molecule types. 
The 'elementary' prefix was added to distinguish it clearly from molecules in 
simulation.

Attributes:
***********
* | **name**: str
  | Name of this elementary molecule type.

* | **components**: List[ComponentType] = None
  | List of components used by this elementary molecule type.

* | **diffusion_constant_2d**: float = None
  | Elementary molecule based on this type is constrained to a surface
  | and diffuses with the specified diffusion constant.
  | D can be zero, in which case the molecule doesn’t move. 
  | The units of D are cm^2 /s.

* | **diffusion_constant_3d**: float = None
  | Elementary molecule based on this type diffuses in space with the 
  | specified diffusion constant D. 
  | D can be zero, in which case the molecule doesn’t move. 
  | The units of D are cm^2 /s.

* | **custom_time_step**: float = None
  | This molecule should take timesteps of length custom_time_step (in seconds). 
  | Use either this or custom_time_step, not both.

* | **custom_space_step**: float = None
  | This molecule should take steps of average length given by the custom_space_step value (in microns). 
  | Use either this or custom_time_step, not both.

* | **target_only**: bool = False
  | This molecule will not initiate reactions when it runs into other molecules. This
  | setting can speed up simulations when applied to a molecule at high concentrations 
  | that reacts with a molecule at low concentrations (it is more efficient for
  | the low-concentration molecule to trigger the reactions). This directive does
  | not affect unimolecular reactions.


Methods:
*********
* | **inst**

   * | components: List[Component] = None
     | Instances of components for the the created elementary molecule.
     | Not all components need to be specified in case when the elementary 
     | molecule is used in a pattern.

   * | compartment_name: str = None
     | Optional specification of compartment name for the created elementary molecule.

   * | return type: ElementaryMolecule


  | Create an elementary molecule based on this elementary molecule type.


* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation.



ReactionRule
============

Attributes:
***********
* | **name**: str = None
  | Name of the reaction. If this is a reversible reaction, then it is the name of the 
  | reaction in forward direction.

* | **reactants**: List[Complex] = None

* | **products**: List[Complex] = None

* | **fwd_rate**: float = None
  | Rates have following units\: unimolecular [s^-1], volume bimolecular [M^-1\*s^-1], 
  | The units of the reaction rate for uni- and bimolecular reactions are
  |   \* [s^-1] for unimolecular reactions,
  |   \* [M^-1\*s^-1] for bimolecular reactions between either two volume molecules, a volume molecule 
  |                 and a surface (molecule), 
  |   \* [um^2\*N^-1\*s^-1] bimolecular reactions between two surface molecules on the same surface, and
  |   \* [N^-1\*s^-1] bimolecular reactions between two surface molecules on different objects 
  |     (this is a highly experimental feature and the unit will likely change in the future, 
  |      not sure if probability is computed correctly, it works the way that the surface molecule 
  |      is first diffused and then a potential collisions within the distance of Config.intermembrane_interaction_radius
  |      are evaluated). 
  | Here, M is the molarity of the solution and N the number of reactants.
  | May be changed after model initialization. 
  | Setting of value is ignored if the rate does not change. 
  | If the new value differs from previous, updates all information related 
  | to the new rate including recomputation of reaction times for molecules if this is a
  | unimolecular reaction.

* | **rev_name**: str = None
  | Name of the reaction in reverse direction.

* | **rev_rate**: float = None
  | Reverse reactions rate, reaction is unidirectional when not specified.
  | May be changed after model initialization, in the case behaves the same was as for 
  | changing the 'fwd_rate'.

* | **variable_rate**: List[List[float]] = None
  | The array passed as this argument must have as its items a pair of floats (time, rate).      
  | Variable rate is applicable only for irreversible reactions. 
  | Members fwd_rate and rev_rate must not be set when setting this attribute through a constructor. 
  | When this attribute is set outside of the class constructor, fwd_rate is automatically reset to an 'unset' value.
  | Cannot be set after model initialization.

* | **is_intermembrane_surface_reaction**: bool = False
  | Experimental, see addintinal explanation in 'fwd' rate.
  | Then set to true, this is a special type of surface-surface reaction that 
  | allows for two surface molecules to react when they are on different geometrical objects. 
  | This support is limited for now, the reaction rule must be in the form of A + B -> C + D 
  | where all reactants and products must be surface molecules and 
  | their orientation must be 'any' (default).


Methods:
*********
* | **to_bngl_str**

   * | return type: str


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

Attributes:
***********
* | **name**: str = None
  | Name of the species in the BNGL format. 
  | One must either specify 'name' or 'elementary_molecules' 
  | (inherited from Complex). This argument 'name' is parsed during model 
  | initialization.

* | **diffusion_constant_2d**: float = None
  | This molecule is constrained to a surface and diffuses with diffusion constant D.

* | **diffusion_constant_3d**: float = None
  | This molecule diffuses in space with diffusion constant D. D can be zero, in which case the molecule doesn’t move. The units of D are cm 2 /s.

* | **custom_time_step**: float = None
  | This molecule should take timesteps of length t (in seconds). Use either this or custom_time_step.

* | **custom_space_step**: float = None
  | This molecule should take steps of average length L (in microns). Use either this or custom_time_step.

* | **target_only**: bool = False
  | This molecule will not initiate reactions when it runs into other molecules. This
  | setting can speed up simulations when applied to a molecule at high concentrations 
  | that reacts with a molecule at low concentrations (it is more efficient for
  | the low-concentration molecule to trigger the reactions). This directive does
  | not affect unimolecular reactions.

* | **name**: str = None
  | When set, this complex instance is initialized from a BNGL string passed as this argument, 
  | the string is parsed and elementary_molecules and compartment are initialized.
  | Only one of name or elementary_molecules can be set.

* | **elementary_molecules**: List[ElementaryMolecule] = None
  | Individual molecule instances contained in the complex.
  | Only one of name or elementary_molecules can be set.

* | **orientation**: Orientation = Orientation.DEFAULT
  | Specifies orientation of a molecule. 
  | When Orientation.DEFAULT if kept then during model initialization is
  | 'orientation' set to Orientation.NONE for volume complexes and to 
  | Orientation.UP for surface complexes.
  | Ignored by derived class Species.

* | **compartment_name**: str = None
  | Specifies compartment name of this Complex.
  | Only one of 'orientation' and 'compartment_name' can be set. 
  | Corresponds to BNGL specification of a compartment for the whole complex '@COMP\:'.
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


Methods:
*********
* | **inst**

   * | orientation: Orientation = Orientation.DEFAULT
     | Maximum one of orientation or compartment_name can be set, not both.

   * | compartment_name: str = None
     | Maximum one of orientation or compartment_name can be set, not both.

   * | return type: Complex


  | Creates a copy of a Complex from this Species with specified orientation and compartment name.


* | **to_bngl_str**

   * | return type: str


  | Creates a string that corresponds to its BNGL representation including compartments.


* | **as_species**

   * | return type: Species


  | Returns a Species object based on this Complex. All species-specific 
  | attributes are set to their default values and 'name' is set to value returned by 
  | 'to_bngl_str()'.



Subsystem
=========

Attributes:
***********
* | **species**: List[Species] = None

* | **reaction_rules**: List[ReactionRule] = None

* | **surface_classes**: List[SurfaceClass] = None

* | **elementary_molecule_types**: List[ElementaryMoleculeType] = None
  | Used mainly when a BNGL file is loaded, if BNGL species is defined through 
  | Python API, this array is populated automatically


Methods:
*********
* | **add_species**

   * | s: Species

* | **find_species**

   * | name: str
   * | return type: Species


* | **add_reaction_rule**

   * | r: ReactionRule

* | **find_reaction_rule**

   * | name: str
   * | return type: ReactionRule


* | **add_surface_class**

   * | sc: SurfaceClass

* | **find_surface_class**

   * | name: str
   * | return type: SurfaceClass


* | **add_elementary_molecule_type**

   * | mt: ElementaryMoleculeType

* | **find_elementary_molecule_type**

   * | name: str
   * | return type: ElementaryMoleculeType


* | **load_bngl_molecule_types_and_reaction_rules**

   * | file_name: str
   * | parameter_overrides: Dict[str, float] = None

  | Parses a BNGL file and only reads molecule types and
  | reaction rules sections, e.g. ignores observables. 
  | Parameter values are evaluated and the result value 
  | is directly used.  
  | Compartments names are stored in rxn rules as strings because
  | compartments belong to geometry objects and the subsystem is independent
  | on specific geometry.
  | However they must be defined on initialization.



SurfaceClass
============

Defining a surface class allows surfaces to behave like species. For instance, one may wish 
to specify that a surface does not block the diffusion of molecules. Each type of surface is defined
by name, and each surface name must be unique in the simulation and should not match any molecule names.
To define a reaction with a surface class, use constructor call m.Complex(name) as one of the reactants.

Attributes:
***********
* | **name**: str
  | Name of the surface class

* | **properties**: List[SurfaceProperty] = None
  | A surface class can either have a list of properties or just one property.
  | In the usual case of having one property, one can set the attributes 
  | type, affected_species, etc. inherited from SurfaceProperty directly.

* | **type**: SurfacePropertyType = SurfacePropertyType.UNSET
  | Must be set. See SurfacePropertyType for options.

* | **affected_complex_pattern**: Complex = None
  | A complex pattern with optional orientation must be set.
  | Default orientation means that the pattern matches any orientation.
  | For concentration or flux clamp the orientation specifies on which side  
  | will be the concentration held 
  | (UP is front or outside, DOWN is back or inside, and DEFAULT, ANY or NONE is on both sides).
  | The complex pattern must not have any compartment.

* | **concentration**: float = None
  | Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or 
  | SurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side 
  | of the wall that has this concentration or flux clamped.

SurfaceProperty
===============

Single property for a SurfaceClass.

Attributes:
***********
* | **type**: SurfacePropertyType = SurfacePropertyType.UNSET
  | Must be set. See SurfacePropertyType for options.

* | **affected_complex_pattern**: Complex = None
  | A complex pattern with optional orientation must be set.
  | Default orientation means that the pattern matches any orientation.
  | For concentration or flux clamp the orientation specifies on which side  
  | will be the concentration held 
  | (UP is front or outside, DOWN is back or inside, and DEFAULT, ANY or NONE is on both sides).
  | The complex pattern must not have any compartment.

* | **concentration**: float = None
  | Specifies concentration when type is SurfacePropertyType.CLAMP_CONCENTRATION or 
  | SurfacePropertyType.CLAMP_FLUX. Represents concentration of the imagined opposide side 
  | of the wall that has this concentration or flux clamped.

