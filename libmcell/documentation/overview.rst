************
API Overview
************

Model Execution Phases
######################

There are two phases in model execution. 
The first phase is **build phase**, where all the components 
of the simulated system are defined by creating objects such as ReactionRule, 
ReleaseSite, or Count. All class attributes of these objects 
are accessible and any modifications are allowed in order to define the model.

The second phase is entered by calling Model.initialize() and is called
**simulation phase**. What precisely Model.initialize() does is that it
transforms information in API objects (such as the ReactionRule object) 
into an internal MCell4 representation. Thus, modifying attributes or 
adding new objects e.g. to reaction rules or geometry objects won't have 
any effect and attempting to do modifications raises an error in most cases. 
There are exceptions that are explicitly described in the API documentation,
for instance changing rate of a reaction is allowed.  

Accessing API Object Attributes
###############################

Attributes of API objects can be divided into **scalar**, **object**, 
and **vector** attributes.
Types of scalar attributes are str (string), bool, float, integer or enumerations.
Types of object attributes are API classes and represent references to 
such objects.     
And the vector attributes are those that use ``List`` in their data type.
 
Internally, getter and setter methods are generated for attributes and they 
are called when an attribute is read or set. 

Scalar Attributes
*****************

As an example, let's take a look on forward rate of a ReactionRule object 
that is a scalar attribute.    

.. code-block:: python

      # get a reference to a ReactionRule object
      reaction = model.find_reaction_rule('my_reaction')

      # reading calls C++ method ReactionRule.get_fwd_rate()
      a = reaction.fwd_rate
      
      # writing calls C++ method ReactionRule.set_fwd_rate(float)
      reaction.fwd_rate = 1e5
      
Reads of a scalar attribute are always allowed and a copy of the value is returned, 
therefore no change is propagated back.
When a write to a scalar attribute is executed, the setter method is executed and 
when we are in the build phase, MCell allows any changes. In the simulation phase, 
the majority of writes cause an exception that disallows modification except for 
several cases where the modification is explicitly allowed such as for the fwd_rate 
in our example. They way how the modifications are handled is that the ReactionRule 
contains a link to the internal MCell4 reaction representation and the call to 
set_fwd_rate besides of modifying the attribute value also informs the MCell4 
engine to change the rate wherever needed. 

Object Attributes
*****************

Object attributes behave practically in the same way as scalar attributes. 
The only difference that instead of returing a value, a reference to an object 
is returned.

.. code-block:: python

      # get a reference to a ReleaseSite object
      rel_site = model.find_release_site('rel_a')
      
      # get a reference to a Region object used by the ReleaseSite
      region = rel_site.region
      
      # create a 'box' object that we will use as a new region
      new_region = geometry_utils.create_box('box', 1)
      
      # and replace the original region, allowed only in build phase
      rel_site.region = new_region
       
Reading a reference is allowed at all times. Accessing individual attributes 
of the object follows the same rules as for any other API object.
Writes (replacing the reference) are only allowed in the build phase, 
there is no object type attribute that is allowed to be changed in the simulation phase.

Vector Attributes
*****************

In order to allow propagation of modifications of lists on the Python side to the 
C++ side and vice versa, special vector classes are used on the Python side instead 
of the standard list type.  
These vector attributes try to emulate Python list type, but currently allow only a 
subset of operations provided for lists.

An example for states of a ComponentType is shown here: 

.. code-block:: python

      # create an object of class ComponentType
      # initially only one state with name 's1' is set 
      ct_u = m.ComponentType('u', states = ['s1'])
      
      # accessing ct_u.states causes a getter method to be called that 
      # returns a reference to an object of class VectorStr that  
      # provides method append
      ct_u.states.append('s2')

      # we can also modify a specific item in this vector       
      ct_u.states[1] = '2'
  
Reading a vector attribute is always allowed. The returned reference to the Vector 
object is not guarded against writes and there are no semantic checks. 
In the build phase, such modifications are used when the model is initialized. 
However, in the simulation phase, such modifications are ignored by the MCell4 engine
and no error is reported. 

When writing to a vector attribute, the original vector is replaced by the new one. 
This is allowed in the build phase and when attempting to replace the whole vector in the
simulation phase, an error is reported.  

Object Cloning Support
######################

API objects support shallow and deep copy operations provided through Python methods
*copy.copy(x)* and *copy.deepcopy(x[, memo])*.

Cloning is allowed even if the model was already initialized.
However, all links in the cloned object to the initialized model are lost. E.g. it is not possible 
to clone a *Count* object and then call the clone's method *get_current_value* because the new object 
will be uninitialized and won't know which model's internal count it is referencing.

Due to MCell4 being implemented primarily in C++, there is one significant difference 
in *copy* from Python semantics. All lists are copied 
by value, not by reference as Python's lists since they are internally implemented 
with C++ std::vector. This behavior is shown in the following code snippet:

.. code-block:: python

   ct = m.ComponentType('u', states = ['0', '1'])
   
   ct_copy = copy.copy(ct3)
   ct_copy.states[0] = 'X' # change item in a copied list
   
   assert ct.states == ['0', '1']
   assert ct_copy.states == ['X', '1']  

For *copy.deepcopy(x[, memo])*, the optional *memo* argument is ignored.

Object Debug Printouts
######################

Each of the API objects provides method *__str__* to convert it to 
a string representation that shows the contents of this object. 
This method is used when a method *print* or cast *str(...)* is used. 
By default, not all details are shown for all objects because that would make the 
output too lengthy (especially for the *GeometryObject* and *Complex* classes). 

The method *__str__* has two arguments *all_details* (default False) and 
*ind* (indent, default ""). To obtain access to all details, set *all_details* 
to True.

.. code-block:: python

   cplx = m.Complex('A(x~0)')
   
   print(cplx) # prints only 'A(x~0)'
   
   print(cplx.__str__(True)) # prints a detailed representation 
