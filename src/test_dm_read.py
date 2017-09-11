import pymcell as m

seed = 1
world = m.MCellSim(seed)
world.set_output_freq(1)

dm = m.read_json_data_model("./lipidraft.json")
m.create_initializations_from_dm(dm, world)
spec_dict = m.create_species_from_dm(dm)
for spec in spec_dict.values():
    world.add_species(spec)
rxn_list = m.create_reactions_from_dm(dm, spec_dict)
for rxn in rxn_list:
    world.add_reaction(rxn)
meshobj_dict = m.create_meshobjs_from_dm(dm)
for meshobj_name in meshobj_dict:
    world.add_geometry(meshobj_dict[meshobj_name])
sc_dict = m.create_surface_classes_from_dm(dm, world, spec_dict)
m.create_mod_surf_reg_from_dm(dm, world, sc_dict, meshobj_dict)
m.create_release_sites_from_dm(dm, world, meshobj_dict, spec_dict)
m.create_reaction_data_from_dm(dm, world, meshobj_dict, spec_dict)
m.create_viz_data_from_dm(dm, world, spec_dict)
world.run_sim()
