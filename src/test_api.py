#!/usr/bin/env python

import pymcell as m


class Vector3(object):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

"""
def create_count(world, where, mol_sym, file_path):
    report_flags = m.REPORT_CONTENTS
    c_list = m.output_column_list()
    # XXX: m.ORIENT_NOT_SET is using -100 instead of SHRT_MIN (used typemap for
    # mcell_create_count in mcell_react_out.i) because limits.h does not work
    # well with swig
    count_list = m.mcell_create_count(
            world, mol_sym, m.ORIENT_NOT_SET, where, report_flags, None,
            c_list)

    os = m.output_set()
    os = m.mcell_create_new_output_set(
        None, 0, count_list.column_head, m.FILE_SUBSTITUTE, file_path)

    out_times = m.output_times_inlist()
    out_times.type = m.OUTPUT_BY_STEP
    out_times.step = 1e-5

    output = m.output_set_list()
    output.set_head = os
    output.set_tail = os

    m.mcell_add_reaction_output_block(world, output, 10000, out_times)

    return (count_list, os, out_times, output)


def create_species(world, name, D, is_2d):
    species_def = m.mcell_species_spec()
    species_def.name = name
    species_def.D = D
    is_2d = 1 if is_2d else 0
    species_def.is_2d = is_2d
    species_def.custom_time_step = 0
    species_def.target_only = 0
    species_def.max_step_length = 0

    species_temp_sym = m.mcell_symbol()
    species_sym = m.mcell_create_species(world, species_def, species_temp_sym)

    return species_sym


def create_reaction(
        world, reactants, products, rate_constant, backward_rate_constant=None,
        surfs=None, name=None):

    if surfs:
        pass  # XXX: need to add this bit
    else:
        surfs = m.mcell_add_to_species_list(None, False, 0, None)

    arrow = m.reaction_arrow()
    # reversible reaction e.g. A<->B
    if backward_rate_constant:
        arrow.flags = m.ARROW_BIDIRECTIONAL
        rate_constant = m.mcell_create_reaction_rates(
                m.RATE_CONSTANT, rate_constant, m.RATE_CONSTANT,
                backward_rate_constant)
    # irreversible reaction e.g. A->B
    else:
        arrow.flags = m.REGULAR_ARROW
        rate_constant = m.mcell_create_reaction_rates(
                m.RATE_CONSTANT, rate_constant, m.RATE_UNSET, 0)
    arrow.catalyst = m.mcell_species()
    arrow.catalyst.next = None
    arrow.catalyst.mol_type = None
    arrow.catalyst.orient_set = 0
    arrow.catalyst.orient = 0

    if (name):
        name_sym = m.mcell_new_rxn_pathname(world, name)
    else:
        name_sym = None
    m.mcell_add_reaction_simplified(
            world, reactants, arrow, surfs, products, rate_constant, name_sym)


def create_instance_object(world, name):
    scene_temp = m.object()
    return m.mcell_create_instance_object(world, name, scene_temp)


def create_surf_class(world, name):
    sc_temp = m.mcell_symbol()
    return m.mcell_create_surf_class(world, name, sc_temp)


def create_release_site(world, scene, pos, diam, shape, number, mol_sym, name):
    position = m.vector3()
    position.x = pos.x
    position.y = pos.y
    position.z = pos.z
    diameter = m.vector3()
    diameter.x = diam.x
    diameter.y = diam.y
    diameter.z = diam.z

    mol_list = m.mcell_add_to_species_list(mol_sym, False, 0, None)
    rel_object = m.object()
    release_site_name = name  # This should derive from the mol_sym name
    release_object = m.mcell_create_geometrical_release_site(
            world, scene, release_site_name, shape, position, diameter,
            mol_list, number, 1, None, rel_object)
    m.mcell_delete_species_list(mol_list)
    return (position, diameter, release_object)


def create_box_verts_elems(half_length):
    hl = half_length
    verts = m.mcell_add_to_vertex_list(hl, hl, -hl, None)
    verts = m.mcell_add_to_vertex_list(hl, -hl, -hl, verts)
    verts = m.mcell_add_to_vertex_list(-hl, -hl, -hl, verts)
    verts = m.mcell_add_to_vertex_list(-hl, hl, -hl, verts)
    verts = m.mcell_add_to_vertex_list(hl, hl, hl, verts)
    verts = m.mcell_add_to_vertex_list(hl, -hl, hl, verts)
    verts = m.mcell_add_to_vertex_list(-hl, -hl, hl, verts)
    verts = m.mcell_add_to_vertex_list(-hl, hl, hl, verts)

    elems = m.mcell_add_to_connection_list(1, 2, 3, None)
    elems = m.mcell_add_to_connection_list(7, 6, 5, elems)
    elems = m.mcell_add_to_connection_list(0, 4, 5, elems)
    elems = m.mcell_add_to_connection_list(1, 5, 6, elems)
    elems = m.mcell_add_to_connection_list(6, 7, 3, elems)
    elems = m.mcell_add_to_connection_list(0, 3, 7, elems)
    elems = m.mcell_add_to_connection_list(0, 1, 3, elems)
    elems = m.mcell_add_to_connection_list(4, 7, 5, elems)
    elems = m.mcell_add_to_connection_list(1, 0, 5, elems)
    elems = m.mcell_add_to_connection_list(2, 1, 6, elems)
    elems = m.mcell_add_to_connection_list(2, 6, 3, elems)
    elems = m.mcell_add_to_connection_list(4, 0, 7, elems)

    return (verts, elems)
"""

def main():
    world = m.mcell_create()
    m.mcell_init_state(world)

    dt = 1e-6
    iterations = 100
    m.mcell_set_time_step(world, dt)
    m.mcell_set_iterations(world, iterations)

    # Define a volume molecule named "x"
    sm1_sym = m.create_species(m, world, "sm1", 1e-6, True)
    vm1_sym = m.create_species(m, world, "vm1", 1e-6, False)
    vm2_sym = m.create_species(m, world, "vm2", 1e-6, False)
    vm3_sym = m.create_species(m, world, "vm3", 1e-6, False)

    # Define reactions
    # vm1 + vm2 -> vm3 [1e8]
    reactants1 = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
    reactants1 = m.mcell_add_to_species_list(vm2_sym, False, 0, reactants1)
    products1 = m.mcell_add_to_species_list(vm3_sym, False, 0, None)
    m.create_reaction(m, world, reactants1, products1, 1e8)
    # vm3 -> NULL [1e5]
    reactants2 = m.mcell_add_to_species_list(vm3_sym, False, 0, None)
    m.create_reaction(m, world, reactants2, None, 0.01, name="rxn")

    scene_name = "Scene"
    scene = m.create_instance_object(m, world, scene_name)

    # Create a spherical release site
    pos_vec3 = Vector3()
    diam_vec3 = Vector3(0.015, 0.015, 0.015)
    # XXX: It seems to be necessary to return some or all of these objects in
    # order to have a functioning release site even though we don't use them
    # anywhere after this call.
    position, diameter, release_object = m.create_release_site(
            m, world, scene, pos_vec3, diam_vec3, m.SHAPE_SPHERICAL, 500,
            vm1_sym, "vm1_rel")
    pos_vec3b = Vector3(0.05, 0.05, 0.00)
    diam_vec3b = Vector3(0.025, 0.025, 0.05)
    position2, diameter2, release_object2 = m.create_release_site(
            m, world, scene, pos_vec3b, diam_vec3b, m.SHAPE_CUBIC, 500,
            vm2_sym, "vm2_rel")

    # Create box object
    verts, elems = m.create_box_verts_elems(m, 0.1)

    pobj = m.poly_object()
    obj_name = "Cube"
    pobj.obj_name = obj_name
    pobj.vertices = verts
    pobj.num_vert = 8
    pobj.connections = elems
    pobj.num_conn = 12

    mesh_temp = m.object()
    mesh = m.mcell_create_poly_object(world, scene, pobj, mesh_temp)

    # Create surface region on the box consisting of two triangles
    # XXX: Creating a region is currently required when creating mesh objects
    test_region = m.mcell_create_region(world, mesh, "reg")
    region_list = m.mcell_add_to_region_list(None, 0)
    region_list = m.mcell_add_to_region_list(region_list, 6)
    m.mcell_set_region_elements(test_region, region_list, 1)

    # create surface class
    sc_sm1 = m.create_surf_class(m, world, "sc_release_y")
    # create releases using a surface class (i.e. not a release object)
    # mdl equivalent: MOLECULE_DENSITY {sm1' = 1000}
    sm1 = m.mcell_add_to_species_list(sm1_sym, True, 1, None)
    smd = m.mcell_add_mol_release_to_surf_class(
      world, sc_sm1, sm1, 10000, 0, None)
    # create surface class that is reflective to sm1
    m.mcell_add_surf_class_properties(world, m.RFLCT, sc_sm1, sm1_sym, 0)
    # m.mcell_add_surf_class_properties(world, m.SINK, sc_sm1, sm1_sym, 0)
    m.mcell_assign_surf_class_to_region(sc_sm1, test_region)

    m.mcell_delete_species_list(sm1)

    
    # Create reaction on sc_sm1 sm1, -> sm1'
    reactantsSurf = m.mcell_add_to_species_list(sm1_sym, True, 1, None)
    productsSurf = m.mcell_add_to_species_list(sm1_sym, True, -1, None)
    productsSurf = m.mcell_add_to_species_list(sm1_sym, True, -1, productsSurf)
    m.create_reaction(m, world, reactantsSurf, productsSurf, .1, surfs=sm1, name="rxnSurf")
    

    # Create viz data
    viz_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
    viz_list = m.mcell_add_to_species_list(vm2_sym, False, 0, viz_list)
    viz_list = m.mcell_add_to_species_list(vm3_sym, False, 0, viz_list)
    viz_list = m.mcell_add_to_species_list(sm1_sym, True, 0, viz_list)
    m.mcell_create_viz_output(
            world, "./viz_data/test", viz_list, 0, iterations, 1)

    # Create reaction data
    mesh_sym = m.mcell_get_obj_sym(mesh)
    count_list1, os1, out_times1, output1 = m.create_count(
            m, world, mesh_sym, vm1_sym, "react_data/vm1_cube.dat")
    reg_sym = m.mcell_get_reg_sym(test_region)
    count_list2, os2, out_times2, output2 = m.create_count(
            m, world, reg_sym, sm1_sym, "react_data/sm1_reg.dat")
    count_list3, os3, out_times3, output3 = m.create_count(
            m, world, None, vm1_sym, "react_data/vm1_world.dat")
    count_list4, os4, out_times4, output4 = m.create_count(
            m, world, mesh_sym, vm3_sym, "react_data/vm3_cube.dat")

    m.mcell_init_simulation(world)
    m.mcell_init_output(world)

    output_freq = 10
    for i in range(iterations):
        vm3_count = m.mcell_get_count(
                "vm3", "%s.%s,ALL" % (scene_name, obj_name), world)
        # When vm3 hits some arbitrary threshold value (i.e. 150), ramp up the
        # rate constant of vm3->NULL. This is just a simple test, but we'll
        # need to do something analagous when interfacing with pyNEURON.
        if (vm3_count > 150): 
            m.mcell_modify_rate_constant(world, "rxn", 1e8)
        m.mcell_run_iteration(world, output_freq, 0)
    m.mcell_flush_data(world)
    # m.mcell_print_final_warnings(world)
    # m.mcell_print_final_statistics(world)

if __name__ == "__main__":
    main()
