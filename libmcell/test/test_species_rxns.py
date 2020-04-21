#!/usr/bin/env python3

import mcell as m

C = m.ComponentType('C', ['0', '1', '2'])

#C_inst = C.inst(state='1', bond=0)
C_inst = C.inst('1', bond=0)


print("ComponentType:")
print(C)

print("ComponentInstance:")
print(C_inst)