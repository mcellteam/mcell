import pymcell as m

dm = m.read_json_data_model("./lipidraft.json")
spec_dict = m.create_species_from_dm(dm)
rxn_list = m.create_reactions_from_dm(dm, spec_dict)
