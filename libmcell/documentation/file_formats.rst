************
File Formats
************

This section describes some of the output file formats that MCell produces.


Molecule Visualization Data
###########################


MCell4 provides two output formats for molecule visualization. 
A textual and binary format. 
VizMode.ASCII mode selects a readable representation, VizMode.CELLBLENDER mode
selects a binary representation that cannot be read using a text editor but 
is faster to generate and read.
The code below shows an example of the VizOutput object with a selection of output 
file format. 

.. code-block:: python

   viz_output = m.VizOutput(
       mode = m.VizMode.ASCII, # or m.VizMode.CELLBLENDER
       output_files_prefix = './viz_data/seed_' + str(SEED).zfill(5) + '/Scene'
   )
   model.add_viz_output(viz_output)


Text Visualization Data
***********************

The text format is composed of lines where each line contains data on one molecule 
in the following format:

.. code-block:: text

   species id x y z nx ny nz  

Here, *species* is the string representation of the molecule species, 
*id* is a unique integer identifier of the molecule, *x* *y* *z* 
represent the 3D location of the molecule in um (micrometers), and *nx* *ny* *nz* is 
the normal vector. For volume molecules, the normal vector is always 0 0 0. 
For surface molecules, it corresponds to the normal of the wall where the molecule is 
located and can be either identical (when the molecule's orientation is UP), or 
negated (when the molecule's orientation is DOWN).
   
Example:

.. code-block:: text

   va@CP 1 0.0601783518 -0.0796582 -0.0180513811 0 0 0
   sb@PM 2 0.125 -0.115164035 0.066618027 1 0 0
   sc@PM 3 -0.113072306 0.00838499159 0.125 0 0 1

Here we have three molecules, first is a volume molecule with species va\@CP 
(compartment is a part of molecule's species), then two surface molecules with 
species sb\@PM and sc\@PM with their normal vectors set.

Binary Visualization Data
*************************

The binary format uses the following structure.
All floating point values (float) are encoded as 32-bit 
IEEE 745 single-precision floating-point format and the 
location uses the um (micrometer) unit. Values of indented addresses
represent offset from start of each block.  


.. code-block:: text

   0x0: 1 (4 bytes, uint), constant, used to distinguish binary and ascii formats
   0x4: here start species blocks, each block contains:
        0x0:        name_len (1 byte, unsigned char)
        0x1:        name (sequence of 1-byte characters representing the species name without the terminating zero)
        name_len+1: species_type (1 byte, unsigned char), 0 for volume molecules and 1 for surface molecules  
        name_len+2: num_xyz_positions (4 bytes, uint)
        name_len+6: x,y,z triplets, repeated num_xyz_positions times 
            0x0: x (4 bytes, float)
            0x4: y (4 bytes, float)
            0x8: z (4 bytes, float)
        name_len+2+num_xyz_positions*12: if species_type is 1, on this address start nx,ny,nz triplets, otherwise empty (0 bytes)      
            0x0: nx (4 bytes, float)
            0x4: ny (4 bytes, float)
            0x8: nz (4 bytes, float)
  
Example:

.. code-block:: text

   00000000  01 00 00 00 05 76 61 40  43 50 00 06 00 00 00 7f  |.....va@CP......|
   00000010  50 9f 3e 13 b3 14 3e d5  b2 a4 3e 38 54 4c bd 3c  |P.>...>...>8TL.<|
   00000020  b9 d4 be 3a 9f da 3e 05  73 62 40 50 4d 01 06 00  |...:..>.sb@PM...|
   00000030  00 00 00 00 00 3f 40 f7  82 bd e0 d1 2e 3e e8 30  |.....?@......>.0|
   00000040  9d 3e 00 00 00 3f 0a 50  32 be 00 00 80 3f 00 00  |.>...?.P2....?..|
   00000050  00 00 00 00 00 00 00 00  00 80 00 00 80 3f 00 00  |.............?..|
   00000060  00 00                                             |..|
   00000062

This hex dump shows a binary representation of molecule positions equivalent to
ascii output shown here:  

.. code-block:: text

   va@CP 0 0.311161 0.145214365 0.32167688 0 0 0
   va@CP 1 -0.0498850037 -0.415475712 0.426996041 0 0 0
   sb@PM 2 0.5 -0.0639481536 0.170722492 1 0 0
   sb@PM 3 0.307013747 0.5 -0.174133457 -0 1 0 

In a more detail, this is how the data is encoded:

.. code-block:: text

   00000000  01 00 00 00     - constant '1' (in little-endian representation)
   00000004  05              - length of string 'va@CP'
   00000005  76 61 40 43 50  - characters of string 'va@CP'
   0000000A  00              - '0' telling that these are volume molecules
   0000000B  06 00 00 00     - 6 num_xyz_positions
   0000000F  7f 50 9f 3e  13 b3 14 3e  d5 b2 a4 3e - xyz position of the first molecule
   0000001B  38 54 4c bd  3c b9 d4 be  3a 9f da 3e - xyz position of the second molecule
   00000027  05              - length of string 'sb@PM'
   ... (species name, positions, and normals of sb@PM molecules follow) 
   
.. todo - link to tests/pymcell4_positive/0500_cellblender_viz_output
