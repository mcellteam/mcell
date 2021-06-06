*****************
Model Conversions
*****************

Overview
########

This section desribes how to create an MCell4 model from 
an MDL or Data Model representations.
MDL is a model definition language used by MCell3 
(`MDL Quick Reference Guide <https://mcell.org/documentation/MCellQuickReferenceGuide.pdf>`).
Data Model is an internal JSON representation used primarily in CellBlender.
When converting from MDL to MCell4 Python representation,  

The guide below expect that a system variable *MCELL_PATH* was set as described in 

MDL to Data Model
#################

The fist step is to convert MDL representation to JSON Data Model. 
Navigate to the directory where your MDL files reside and run:

.. code-block:: text

      $MCELL_PATH/mcell Scene.main.mdl -mdl2datamodel4


MCell loads the MDL files, creates internal representation, and 
dumps it into a file called *data_model.json*.

This representation can be then eeither converted to an MCell4 model as described in the 
following section or imported by CellBlender through 
**File** -> **Import** -> **CellBlender Model and Geometry (JSON)**.

A drawback of this process is that all parameter names and 
their expressions will be lost because the source for the data model generation
is an already processed model with all parameters evaluated.
If the parameters and expressions are needed, it is necessary to correct 
the resulting model manually.  


Data Model to MCell4 Python Code
################################

The next step, also used internally by CellBlender is to convvert the data model into
Python 

.. code-block:: text

