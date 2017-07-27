import pymcell as m

iterations = 100
time_step = 1e-5
seed = 1
world = m.MCellSim(seed)
world.set_time_step(time_step)
world.set_iterations(iterations)
world.set_output_freq(1)

dm = m.read_json_data_model("./lipidraft.json")
spec_dict = m.create_species_from_dm(dm)
for spec in spec_dict.values():
    world.add_species(spec)
rxn_list = m.create_reactions_from_dm(dm, spec_dict)
# for rxn in rxn_list:
#     world.add_reaction(rxn)
meshobj_dict = m.create_meshobjs_from_dm(dm)
for meshobj_name in meshobj_dict:
    world.add_geometry(meshobj_dict[meshobj_name])
m.create_release_sites_from_dm(dm, world, meshobj_dict, spec_dict)
world.run_sim()
