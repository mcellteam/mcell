#!/usr/bin/env python3

import mcell as m

C = m.ComponentType('C', ['0', '1', '2'])
N = m.ComponentType('N', ['0', '1', '2'])


print("ComponentType:")
print(C)

#C_inst = C.inst(state='1', bond=0)
C_inst = C.inst('1', bond=0)

C_inst2 = C.inst(1)

print("ComponentInstance:")
print(C_inst)


# diffusion constant here or with the complex?
# maybe here as well for automatic 
CaM = m.MoleculeType('CaM', [C, N])

print("MoleculeType:")
print(CaM)


C_inst0 = C.inst(0)
N_inst1 = N.inst(1)

CaM_inst = CaM.inst([C_inst0, N_inst1])

# does not work..
#CaM_inst = CaM.inst([C.inst(0), N.inst(1)])

print("MoleculeInstance:")
print(CaM_inst)
