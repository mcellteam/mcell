import pymcell as m
import json
from typing import List, Dict


def read_json_data_model(file_name: str):
    with open(file_name, 'r') as f:
        json_model = f.read()
        data_model = json.loads(json_model)
    return data_model


def create_species_from_dm(data_model: Dict) -> List[m.Species]:
    species_dm_list = data_model['mcell']['define_molecules']['molecule_list']
    species_dict = {}
    for species_dm in species_dm_list:
        species_name = species_dm['mol_name']
        dc = float(species_dm['diffusion_constant'])
        species_type = species_dm['mol_type']
        surface = False
        if species_type == '2D':
            surface = True
        species_dm = m.Species(species_name, dc, surface)
        species_dict[species_name] = species_dm
    return species_dict


def make_spec_orient_list(mol_str_list, species):
    spec_orient_list = []
    for r in mol_str_list:
        if r.endswith("'") or r.endswith(","):
            r_str = r[:-1]
            r_orient = r[-1]
            if r_orient == "'":
                orient = m.Orient.up
            elif r_orient == ",":
                orient = m.Orient.down
            else:
                orient = m.Orient.mix
        else:
            r_str = r
            orient = m.Orient.mix
        spec = species[r_str]
        spec_orient = (spec, orient)
        spec_orient_list.append(spec_orient)
    return spec_orient_list


def create_reactions_from_dm(
        data_model: Dict, species: Dict[str, m.Species]) -> List[m.Reaction]:
    rxn_dm_list = data_model['mcell']['define_reactions']['reaction_list']
    rxn_list = []
    for rxn_dm in rxn_dm_list:
        rxn_name = rxn_dm['rxn_name']
        fwd_rate = float(rxn_dm['fwd_rate'])
        try:
            bkwd_rate = float(rxn_dm['bkwd_rate'])
        except ValueError:
            pass
        reactants_str_list = rxn_dm['reactants'].split(" + ")
        products_str_list = rxn_dm['products'].split(" + ")
        r_list = make_spec_orient_list(reactants_str_list, species)
        p_list = make_spec_orient_list(products_str_list, species)
        rxn_type = rxn_dm['rxn_type']
        rxn_dm = m.Reaction(r_list, p_list, fwd_rate)
        rxn_list.append(rxn_dm)
    return rxn_list


def create_meshobjs_from_dm(dm: Dict):
    meshobj_dm_list = dm['mcell']['geometrical_objects']['object_list']
    meshobj_list = []
    for meshobj in meshobj_dm_list:
        name = meshobj['name']
        vert_list = meshobj['vertex_list']
        face_list = meshobj['element_connections']
        meshobj_list.append(m.MeshObj(name, vert_list, face_list))
    return meshobj_list
