#!/usr/bin/env python3

import sys
import os

MCELL_PATH = os.environ.get('MCELL_PATH', '')
if MCELL_PATH:
    sys.path.append(os.path.join(MCELL_PATH, 'lib'))
else:
    print("Error: variable MCELL_PATH that is used to find the mcell library was not set.")
    sys.exit(1)
    
import mcell as m

C = m.ComponentType('C', ['0', '1', '2'])
N = m.ComponentType('N', ['0', '1', '2'])


print("** ComponentType:")
print(C)

#C_inst = C.inst(state='1', bond=0)
C_inst = C.inst('1', 0)

C_inst2 = C.inst(1)

print("** ComponentInstance:")
print(C_inst)

# diffusion constant here or with the complex?
# maybe here as well for automatic 
CaM = m.ElementaryMoleculeType('CaM', [C, N], diffusion_constant_3d = 1e-6)

print("** MoleculeType:")
print(CaM)

CaM_inst = CaM.inst([C.inst(0), N.inst(1)])

print("** MoleculeInstance:")
print(CaM_inst)


cplx_inst = m.ComplexInstance(
    [CaM.inst([C.inst(0), N.inst(1, bond=1)]), CaM.inst([C.inst(0), N.inst(1, bond=0)])])

print("** ComplexInstance:")
print(cplx_inst)

# TODO: can I somehow reorder arguments so that molecule_types is second?
# we will probably see according to other examples
CaMC0N1_species = m.Species('CaM(C~0,N~1)', elementary_molecule_instances = [ CaM.inst([C.inst(0), N.inst(1)]) ] )

print("** Species:")
print(CaMC0N1_species)


d = m.ComponentType('d') # no states  
l = m.ComponentType('l')
r = m.ComponentType('r')
Y286 = m.ComponentType('Y286', ['0','P'])
S306 = m.ComponentType('S306', ['0','P'])
cam = m.ComponentType('cam')

CaMKII = m.ElementaryMoleculeType(
    'CaMKII', 
    [d, r, l, Y286, S306],    
    diffusion_constant_3d = 1e-6
)

V = 0.125*1e-15 # um^3 -> liters
NA = 6.022e23/1e6
k_onCaMKII = 50/(NA*V) #1/uM 1/s 
k_offCaMKII = 60 #1/s 

rxn_rule = m.ReactionRule(
    name = "sixth rxn",
    rev_name = "sixth rxn rev",
    reactants=[
        m.ComplexInstance( 
            [ CaMKII.inst( [ l.inst(), r.inst(), Y286.inst('0'), cam.inst(bond=m.BOND_BOUND) ] ) ] 
        ),
        m.ComplexInstance( 
            [ CaMKII.inst( [ l.inst(), r.inst(), cam.inst(bond=m.BOND_BOUND) ] ) ]  
        )
    ], 
    products=[
        m.ComplexInstance( 
            [ CaMKII.inst( [ l.inst(1, bond=1), r.inst(), Y286.inst('0'), cam.inst(bond=m.BOND_BOUND) ] ), 
              CaMKII.inst( [ l.inst(), r.inst(bond=1), cam.inst(bond=m.BOND_BOUND) ] )
            ]
        )
    ],
    fwd_rate = k_onCaMKII,
    rev_rate = k_offCaMKII
)

print("** ReactionRule:")
print(rxn_rule)
