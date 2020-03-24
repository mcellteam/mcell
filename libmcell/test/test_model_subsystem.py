#!/usr/bin/env python3

import mcell as m

iterations = 100

s1 = m.Species(name = 'my_species', diffusion_constant_2d = 1e-6)

# define a subsystem with its species and reactions
subs = m.Subsystem()
subs.add_species(s1)

# add subsystem to a model
model = m.Model()
model.add_subsystem(subs)

for i in range(iterations):
    model.run_iterations(1)
