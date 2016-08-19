#!/usr/bin/env python

import pymcell as m
import torus


class Vector3(object):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z


def main():
    world = m.mcell_create()
    m.mcell_init_state(world)

    dt = 1e-6
    iterations = 100
    m.mcell_set_time_step(world, dt)
    m.mcell_set_iterations(world, iterations)

    # Define one surface molecule and three volume molecules
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
    position, diameter, sphere_release_object = m.create_release_site(
            m, world, scene, pos_vec3, diam_vec3, m.SHAPE_SPHERICAL, 500,
            vm1_sym, "vm1_rel")
    pos_vec3b = Vector3(0.05, 0.05, 0.00)
    diam_vec3b = Vector3(0.025, 0.025, 0.05)
    position2, diameter2, cube_release_object = m.create_release_site(
            m, world, scene, pos_vec3b, diam_vec3b, m.SHAPE_CUBIC, 500,
            vm2_sym, "vm2_rel")

    # Create box object
    # verts, elems = m.create_box_verts_elems(m, 0.1)

    obj_name = "Torus"
    mesh = m.create_polygon_object(
        m, world, torus.vert_list, torus.face_list, scene, obj_name)

    # Create surface region on half the torus
    # XXX: Creating a region is currently required when creating mesh objects
    region_name = "half_torus"
    test_region = m.create_surface_region(
        m, world, mesh, torus.surf_reg_face_list, region_name)

    region_release_object = m.create_region_release_site(
        m, world, scene, mesh, "vm1_torus_rel", "ALL", 1000, vm1_sym)

    # create surface class
    sc_sm1_sym = m.create_surf_class(m, world, "sc_release_y")
    # create releases using a surface class (i.e. not a release object)
    # mdl equivalent: MOLECULE_DENSITY {sm1' = 1000}
    sm1 = m.mcell_add_to_species_list(sm1_sym, True, 1, None)
    smd = m.mcell_add_mol_release_to_surf_class( world, sc_sm1, sm1, 1000, 0, None)

    # create surface class that is reflective to sm1
    m.mcell_add_surf_class_properties(world, m.RFLCT, sc_sm1_sym, sm1_sym, 0)
    # m.mcell_add_surf_class_properties(world, m.SINK, sc_sm1, sm1_sym, 0)
    m.mcell_assign_surf_class_to_region(sc_sm1_sym, test_region)

    m.mcell_delete_species_list(sm1)

    
    # Create reaction on sc_sm1 sm1, -> sm1'
    sc_surf = m.mcell_add_to_species_list(sc_sm1_sym,True, 1, None)
    reactantsSurf = m.mcell_add_to_species_list(sm1_sym, True, 1, None)
    productsSurf = m.mcell_add_to_species_list(sm1_sym, True, -1, None)
    productsSurf = m.mcell_add_to_species_list(sm1_sym, True, -1, productsSurf)
    m.create_reaction(
        m, world, reactantsSurf, productsSurf, .1, surf_class=sm1,
        name="rxnSurf")
    

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
            m, world, mesh_sym, vm1_sym, "react_data/vm1_%s.dat" % obj_name)
    reg_sym = m.mcell_get_reg_sym(test_region)
    count_list2, os2, out_times2, output2 = m.create_count(
            m, world, reg_sym, sm1_sym, "react_data/sm1_reg.dat")
    count_list3, os3, out_times3, output3 = m.create_count(
            m, world, None, vm1_sym, "react_data/vm1_world.dat")
    count_list4, os4, out_times4, output4 = m.create_count(
            m, world, mesh_sym, vm3_sym, "react_data/vm3_%s.dat" % obj_name)

    m.mcell_init_simulation(world)
    m.mcell_init_output(world)

    output_freq = 10
    for i in range(iterations):
        # vm3_count = m.mcell_get_count(
        #         "vm3", "%s.%s,ALL" % (scene_name, obj_name), world)
        # # When vm3 hits some arbitrary threshold value (i.e. 150), ramp up the
        # # rate constant of vm3->NULL. This is just a simple test, but we'll
        # # need to do something analagous when interfacing with pyNEURON.
        # if (vm3_count > 150): 
        #     m.mcell_modify_rate_constant(world, "rxn", 1e8)
        m.mcell_run_iteration(world, output_freq, 0)
    m.mcell_flush_data(world)
    # m.mcell_print_final_warnings(world)
    # m.mcell_print_final_statistics(world)

if __name__ == "__main__":
    main()
