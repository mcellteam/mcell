#!/usr/bin/env python3

import mcell as m

C = m.ComponentType('C', ['0', '1', '2'])
N = m.ComponentType('N', ['0', '1', '2'])


print("ComponentType:")
print(C)

#C_inst = C.inst(state='1', bond=0)
C_inst = C.inst('1', 0)

C_inst2 = C.inst(1)

print("ComponentInstance:")
print(C_inst)


# diffusion constant here or with the complex?
# maybe here as well for automatic 
CaM = m.MoleculeType('CaM', [C, N])

print("MoleculeType:")
print(CaM)

CaM_inst = CaM.inst([C.inst(0), N.inst(1)])

print("MoleculeInstance:")
print(CaM_inst)


cplx_inst = m.ComplexInstance([CaM.inst([C.inst(0), N.inst(1, bond=1)]), CaM.inst([C.inst(0), N.inst(1, bond=0)])])

print("ComplexInstance:")
print(cplx_inst)


CaMC0N1_species = m.Species('CaM(C~0,N~1)', [ CaM.inst([C.inst(0), N.inst(1)]) ] )

print("Species:")
print(CaMC0N1_species)