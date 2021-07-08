*******************************
Simulation Environment Overview
*******************************

MCell is an agent-based reaction-diffusion program that simulates the movement and interaction of molecules in 3D space.
MCell 4, the latest version of MCell, provides two different user experiences, one through its visual interface as an 
add-on in Blender 2.93, known as `CellBlender <https://mcell.org/tutorials_iframe.html>`__, the other user experience one through a new Python interface.

MCell 4 provides the users with the flexibility to change between both experiences, or to run the simulations using Python and visualize the simulations in Blender. 
Furthermore, support for the file format MDL that older versions of MCell had, is also provided for backwards compatibility. 

In the newly developed MCell 4, reactions can be written following `BioNetGen <https://www.bionetgen.org>`__ syntax, 
and therefore models can be exported from MCell and imported in BioNetGen for simulation and vice versa, 
integrating two different simulation engines.
