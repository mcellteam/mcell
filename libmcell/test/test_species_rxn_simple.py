#!/usr/bin/env python3

# variant with simple complexes
import mcell as m

# when no ComplexInstance is specified,
# Species is assumed to be simple species such as in MCell  
A = m.Species('A', diffusion_constant_3d = 1e-6)
B = m.Species('B', diffusion_constant_2d = 1e-3)
C = m.Species('C', diffusion_constant_3d = 1e-6)

# A' + B' -> C, + B'
# we still need to make instances 
rxn_rule = m.ReactionRule(
    name = "a_plus_b_to_c",
    reactants = [A.inst(orientation=m.ORIENTATION_UP), B.inst(orientation=m.ORIENTATION_UP)],
    products = [C.inst(orientation=m.ORIENTATION_DOWN), B.inst(orientation=m.ORIENTATION_UP)],
    fwd_rate = 1e5
)

print("** ReactionRule:")
print(rxn_rule)

# alternative ?
# m.ReactionRule("A' + B' -> C, + B'", 1e5)

   