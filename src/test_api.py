#!/usr/bin/env python

import pymcell as m
import torus
import box


def main():
    world = m.mcell_create()
    m.mcell_init_state(world)

    dt = 1e-5
    iterations = 100
    m.mcell_set_time_step(world, dt)
    m.mcell_set_iterations(world, iterations)

    # Define one surface molecule and three volume molecules
    sm1_sym = m.create_species(world, "sm1", 0, True)
    vm1_sym = m.create_species(world, "vm1", 1e-6, False)
    vm2_sym = m.create_species(world, "vm2", 1e-6, False)
    vm3_sym = m.create_species(world, "vm3", 0, False)

    # Define reactions
    # vm1 + vm2 -> vm3 [1e8]
    reactants1 = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
    reactants1 = m.mcell_add_to_species_list(vm2_sym, False, 0, reactants1)
    products1 = m.mcell_add_to_species_list(vm3_sym, False, 0, None)
    m.create_reaction(world, reactants1, products1, 1e7)
    # vm3 -> NULL [1e5]
    reactants2 = m.mcell_add_to_species_list(vm3_sym, False, 0, None)
    m.create_reaction(world, reactants2, None, 0, name="rxn_vm3_null")

    # This is the world object. We just call it "Scene" here to be consistent
    # with the MDL output from Blender.
    scene_name = "Scene"
    scene = m.create_instance_object(world, scene_name)

    # Create a spherical release site
    # Release 500 vm1 molecules in a spherical shape
    pos_vec3 = m.Vector3()
    diam_vec3 = m.Vector3(0.015, 0.015, 0.015)
    # XXX: It seems to be necessary to return some or all of these objects in
    # order to have a functioning release site even though we don't use them
    # anywhere after this call.
    position, diameter, sphere_release_object = m.create_release_site(
        world, scene, pos_vec3, diam_vec3, m.SHAPE_SPHERICAL, 5000, 0, vm1_sym,
        "vm1_rel")
    # Release 500 vm2 molecules in a cubic shape
    pos_vec3b = m.Vector3(0.05, 0.05, 0.00)
    diam_vec3b = m.Vector3(0.025, 0.025, 0.05)
    position2, diameter2, cube_release_object = m.create_release_site(
        world, scene, pos_vec3b, diam_vec3b, m.SHAPE_CUBIC, 5000, 0, vm2_sym,
        "vm2_rel")

    pos_vec3c = m.Vector3(1.0, 0.0, 0.0)
    diam_vec3c = m.Vector3(0.015, 0.015, 0.015)
    position3, diameter3, sphere_release_object2 = m.create_release_site(
        world, scene, pos_vec3c, diam_vec3c, m.SHAPE_SPHERICAL, 500, 0,
        vm1_sym, "vm3_rel")
    pos_vec3d = m.Vector3(1.05, 0.05, 0.00)
    diam_vec3d = m.Vector3(0.025, 0.025, 0.05)
    position4, diameter4, cube_release_object2 = m.create_release_site(
        world, scene, pos_vec3d, diam_vec3d, m.SHAPE_CUBIC, 500, 0, vm2_sym,
        "vm4_rel")

    box_name = "Box"
    box_mesh = m.create_box(world, scene, 0.3, box_name)
    box_region_name = "side"
    box_region = m.create_surface_region(
        world, box_mesh, box.surf_reg_face_list, box_region_name)

    torus_name = "Torus"
    torus_mesh = m.create_polygon_object(
        world, torus.vert_list, torus.face_list, scene, torus_name)

    # Create surface region on half the torus
    # XXX: Creating a region is currently required when creating mesh objects
    torus_region_name = "half_torus"
    torus_region = m.create_surface_region(
        world, torus_mesh, torus.surf_reg_face_list, torus_region_name)

    # Release 1000 vm1 molecules in the torus
    region_release_object = m.create_region_release_site(
        world, scene, torus_mesh, "vm1_torus_rel", "ALL", 1000, 0, vm1_sym)

    # create surface class
    sc_sm1_sym = m.create_surf_class(world, "sc")
    # create releases using a surface class (i.e. not a release object)
    # mdl equivalent: MOLECULE_DENSITY {sm1' = 1000}
    sm1 = m.mcell_add_to_species_list(sm1_sym, True, 1, None)
    smd = m.mcell_add_mol_release_to_surf_class(
        world, sc_sm1_sym, sm1, 1000, 0, None)

    # make "sc" surface class reflective to sm1 and apply it to half the torus
    m.mcell_add_surf_class_properties(world, m.RFLCT, sc_sm1_sym, sm1_sym, 0)
    # m.mcell_add_surf_class_properties(world, m.SINK, sc_sm1, sm1_sym, 0)
    m.mcell_assign_surf_class_to_region(sc_sm1_sym, torus_region)

    m.mcell_delete_species_list(sm1)

    # Create reaction sm1' @ sc' -> sm1,
    sc_surf = m.mcell_add_to_species_list(sc_sm1_sym, True, 1, None)
    reactants_surf = m.mcell_add_to_species_list(sm1_sym, True, 1, None)
    products_surf = m.mcell_add_to_species_list(sm1_sym, True, -1, None)
    m.create_reaction(
        world, reactants_surf, products_surf, 0, surf_class=sc_surf,
        name="rxn_SC")

    # Create viz data
    viz_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
    viz_list = m.mcell_add_to_species_list(vm2_sym, False, 0, viz_list)
    viz_list = m.mcell_add_to_species_list(vm3_sym, False, 0, viz_list)
    viz_list = m.mcell_add_to_species_list(sm1_sym, True, 0, viz_list)
    m.mcell_create_viz_output(
        world, "./viz_data/test", viz_list, 0, iterations, 1)

    # Create reaction data
    # box_sym = m.mcell_get_obj_sym(box_mesh)
    torus_sym = m.mcell_get_obj_sym(torus_mesh)
    count_list1, os1, out_times1, output1 = m.create_count(
        world, torus_sym, vm3_sym, "react_data/vm3_%s.dat" % torus_name, 1e-5)
    torus_reg_sym = m.mcell_get_reg_sym(torus_region)
    count_list2, os2, out_times2, output2 = m.create_count(
        world, torus_reg_sym, sm1_sym, "react_data/sm1_reg.dat", dt)
    count_list3, os3, out_times3, output3 = m.create_count(
        world, None, vm1_sym, "react_data/vm1_world.dat", dt)
    # count_list4, os4, out_times4, output4 = m.create_count(
    #     world, box_sym, vm3_sym, "react_data/vm3_%s.dat" % box_name, 1e-5)

    # Set partitions
    m.create_partitions(world, 0, -1.3, 1.3, 0.05)
    m.create_partitions(world, 1, -1.3, 1.3, 0.05)
    m.create_partitions(world, 2, -0.275, 0.275, 0.05)

    m.mcell_init_simulation(world)
    m.mcell_init_output(world)

    output_freq = 10
    torus_obj = m.PolyObj("Torus", "half_torus", torus.vert_list,
                          torus.face_list, torus.surf_reg_face_list)
    box_obj = m.PolyObj("Box", "side", box.vert_list,
                        box.face_list, box.surf_reg_face_list)
    obj_list = [torus_obj, box_obj]
    for i in range(iterations+1):
        vm3_count = m.mcell_get_count(
            "vm3", "%s.%s,ALL" % (scene_name, torus_name), world)
        # print(vm3_count)
        # When vm3 hits some arbitrary threshold value (i.e. 400), ramp up the
        # rate constant of vm3->NULL. This is just a simple test, but we'll
        # need to do something analagous when interfacing with pyNEURON.
        if (vm3_count > 1500):
            print("changing rxn_vm3_null rate constant at iteration %d" % i)
            m.mcell_modify_rate_constant(world, "rxn_vm3_null", 1e8)
        # When we hit 50 iterations, crank up the rxn on the torus
        if (i == 50):
            print("changing rxn_SC rate constant at iteration %d" % i)
            m.mcell_modify_rate_constant(world, "rxn_SC", 1e8)
        m.mcell_run_iteration(world, output_freq, 0)
    m.mcell_flush_data(world)
    m.mcell_print_final_warnings(world)
    m.mcell_print_final_statistics(world)

if __name__ == "__main__":
    main()
