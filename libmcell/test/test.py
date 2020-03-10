#!/usr/bin/env python3

import mcell                                                                                                                                                                          

s = mcell.Species(name = 'my_species', diffusion_constant_2d = 1e-6)

q = mcell.ReleaseSite(name = 'a', molecule = s, shape= 'sphere', location = mcell.Vec3(1, 2, 3))

print(s)
print(s.check_semantics())
print(q)
print(q.check_semantics())
