#!/usr/bin/env python3

# based on mcell_tests/tests/mdl/0320_2_mols_react_in_box_it_10

import mcell as m

SEED = 1

a = m.Species('a', diffusion_constant_3d = 1e-6)
b = m.Species('b', diffusion_constant_3d = 0.5e-6)
c = m.Species('c', diffusion_constant_3d = 1e-6)

react_a_and_b = m.ReactionRule(
    name = 'react_a_and_b',
    reactants = [ a.inst(), b.inst()],
    products = [ c.inst() ],
    fwd_rate = 5e8
)

# first we create a subsystem with
subs = m.Subsystem()
subs.add_species(a)
subs.add_species(b)
subs.add_species(c)
#TODO: subs.add_species([b, c]) # same name because species is both singular and plural
subs.add_reaction_rule(react_a_and_b)

# subs.add_reaction_rules([react_a_and_b])
 
# ---- box ----
box_vertex_list = [
    [-0.05, -0.05, 0.05], 
    [-0.05, 0.05, -0.05], 
    [-0.05, -0.05, -0.05], 
    [-0.05, 0.05, 0.05], 
    [0.05, 0.05, -0.05], 
    [0.05, 0.05, 0.05], 
    [0.05, -0.05, -0.05], 
    [0.05, -0.05, 0.05]
] # box_vertex_list

box_element_connections = [
    [0, 1, 2], 
    [3, 4, 1], 
    [5, 6, 4], 
    [7, 2, 6], 
    [4, 2, 1], 
    [3, 7, 5], 
    [0, 3, 1], 
    [3, 5, 4], 
    [5, 7, 6], 
    [7, 0, 2], 
    [4, 6, 2], 
    [3, 0, 7]
] # box_element_connections

box = m.GeometryObject(
    name = 'box',
    vertex_list = box_vertex_list,
    element_connections = box_element_connections,
    surface_regions = []
)
# ^^^^ box ^^^^
 
 
rel_a = m.ReleaseSite(
    name = 'rel_a',
    species = a,
    shape = m.Shape.Spherical, # second option is shape which accepts geometry object or region?, or use names?
    #location = m.Vec3(0, 0, 0),
    #number_to_release = 100,
    #release_probability = 1
)
 
rel_b = m.ReleaseSite(
    name = 'rel_b',
    shape = m.Shape.Spherical, # second option is shape which accepts geometry object or region?, or use names?
    location = m.Vec3(0.005, 0, 0),
    species = a,
    number_to_release = 100,
    release_probability = 1
)
 
inst = m.InstantiationData()

inst.instantiate_geometry_object(box)
inst.instantiate_release_site(rel_a)
inst.instantiate_release_site(rel_b) # TODO: lists will be allowed as well


viz_output = m.VizOutput(
    mode = m.VizMode.Ascii,
    filename = './viz_data/seed_' + str(SEED) + '/Scene',
    species_list = [a, b, c], # how top put 'all'?, maybe let's keep it explicit for now
    # species_patterns = [ '*' ] - do we want string pattern matching wildcards?, what about just BNG patterns?
    every_n_timesteps = 1 # do we need range? 
)


# Count is used both for molecules and reactions
count_react_a_and_b_in_box = m.Count(
    reaction_rule = react_a_and_b,
    region = box, # rxn count in world if enclosed_in_object is not set
    filename = './react_data/seed_' + str(SEED) + '/react_a_and_b.World.dat',
    every_n_timesteps = 1 # default
)

count_a_in_world = m.Count(
    species = a,
    # species_pattern / complex_pattern?
    filename = './react_data/seed_' + str(SEED) + '/a.World.dat', # try to deduce default filename? 
    every_n_timesteps = 1 # default
)

count_term_b_in_world = m.CountTerm(
    species = b
)


count_term_b_in_box = m.CountTerm(
    species = b,
    region = box
)

# will be zero because box is not transparent 
# maybe we can change release of b
# count is derived from CountTerm?
count_b_not_in_box = m.Count(
    count_expression = count_term_b_in_world - count_term_b_in_box, 
    # species_pattern / complex_pattern?
    filename = './react_data/seed_' + str(SEED) + '/b.World.dat', # try to deduce default filename? 
    every_n_timesteps = 1 # default
)


observables = m.Observables()
observables.add_viz_output(viz_output)
observables.add_count(count_react_a_and_b_in_box)
observables.add_count(count_b_not_in_box)


model = m.Model()
model.add_subsystem(subs)
#model.add_subsystems([subs])
model.add_instantiation_data(inst)
#model.add_instantiation_data([inst])
model.add_observables(observables)


model.config.time_step = 1e-6 # default
model.config.seed = SEED

model.initialize()
model.run_iterations(10)

