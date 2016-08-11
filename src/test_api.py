import mcellSwig as m


def create_species(world, name, D, is_2d=0):
    species_def = m.mcell_species_spec()
    species_def.name = name
    species_def.D = D
    species_def.is_2d = is_2d
    species_def.custom_time_step = 0
    species_def.target_only = 0
    species_def.max_step_length = 0

    species_temp_sym = m.mcell_symbol()
    species_sym = m.mcell_create_species(world, species_def, species_temp_sym)

    return species_sym


def main():
    world = m.mcell_create()
    m.mcell_init_state(world)

    dt = 1e-6
    iterations = 100
    m.mcell_set_time_step(world, dt)
    m.mcell_set_iterations(world, iterations)

    # Define a volume molecule named "x"
    x_mol_sym = create_species(world, "x", 1e-6, 0)

    scene_temp = m.object()
    scene = m.mcell_create_instance_object(world, "Scene", scene_temp)

    # Create a spherical release site
    position = m.vector3()
    pos_value = 0.0
    position.x = pos_value
    position.y = pos_value
    position.z = pos_value
    diameter = m.vector3()
    diam_value = 0.001
    diameter.x = diam_value
    diameter.y = diam_value
    diameter.z = diam_value

    mol_list = m.mcell_add_to_species_list(x_mol_sym, False, 0, None)
    rel_object = m.object()
    release_object = m.mcell_create_geometrical_release_site(
            world, scene, "X_releaser", m.SHAPE_SPHERICAL, position, diameter,
            mol_list, 5000, 1, None, rel_object)
    m.mcell_delete_species_list(mol_list)

    # Create box object
    verts = m.mcell_add_to_vertex_list(0.2, 0.2, -0.2, None);
    verts = m.mcell_add_to_vertex_list(0.2, -0.2, -0.2, verts);
    verts = m.mcell_add_to_vertex_list(-0.2, -0.2, -0.2, verts);
    verts = m.mcell_add_to_vertex_list(-0.2, 0.2, -0.2, verts);
    verts = m.mcell_add_to_vertex_list(0.2, 0.2, 0.2, verts);
    verts = m.mcell_add_to_vertex_list(0.2, -0.2, 0.2, verts);
    verts = m.mcell_add_to_vertex_list(-0.2, -0.2, 0.2, verts);
    verts = m.mcell_add_to_vertex_list(-0.2, 0.2, 0.2, verts);

    elems = m.mcell_add_to_connection_list(1, 2, 3, None);
    elems = m.mcell_add_to_connection_list(7, 6, 5, elems);
    elems = m.mcell_add_to_connection_list(0, 4, 5, elems);
    elems = m.mcell_add_to_connection_list(1, 5, 6, elems);
    elems = m.mcell_add_to_connection_list(6, 7, 3, elems);
    elems = m.mcell_add_to_connection_list(0, 3, 7, elems);
    elems = m.mcell_add_to_connection_list(0, 1, 3, elems);
    elems = m.mcell_add_to_connection_list(4, 7, 5, elems);
    elems = m.mcell_add_to_connection_list(1, 0, 5, elems);
    elems = m.mcell_add_to_connection_list(2, 1, 6, elems);
    elems = m.mcell_add_to_connection_list(2, 6, 3, elems);
    elems = m.mcell_add_to_connection_list(4, 0, 7, elems);

    pobj = m.poly_object()
    pobj.obj_name = "aBox"
    pobj.vertices = verts
    pobj.num_vert = 8
    pobj.connections = elems
    pobj.num_conn = 12
    
    mesh_temp = m.object()
    mesh = m.mcell_create_poly_object(world, scene, pobj, mesh_temp)

    # Create surface region on box
    # XXX: Creating a region is currently required when creating mesh objects
    test_region = m.mcell_create_region(world, mesh, "reg")
    region_list = m.mcell_add_to_region_list(None, 0)
    region_list = m.mcell_add_to_region_list(region_list, 1)
    m.mcell_set_region_elements(test_region, region_list, 1)

    # Create viz data
    viz_list = m.mcell_add_to_species_list(x_mol_sym, False, 0, None)
    m.mcell_create_viz_output(
            world, "./viz_data/test", viz_list, 0, iterations, 1)

    # Create reaction data
    report_flags = m.REPORT_WORLD
    c_list = m.output_column_list()
    # -100 is in the place of ORIENT_NOT_SET (used typemap for
    # mcell_create_count in mcell_react_out.i) because limits.h stuff does not
    # work well with swig.......
    # count_list = m.mcell_create_count(
    #         world, x_mol_sym, -100, None, report_flags, None, c_list)
    # print("count list")
    # print(count_list)

    """
    os = m.output_set()
    os = m.mcell_create_new_output_set(
        None, 0, count_list.column_head, m.FILE_SUBSTITUTE,
        "react_data/foobar.dat")

    outTimes = m.output_times_inlist()
    outTimes.type = OUTPUT_BY_STEP
    outTimes.step = 1e-5

    output = m.output_set_list()
    output.set_head = os
    output.set_tail = os

    m.mcell_add_reaction_output_block(world, output, 10000, outTimes)
    """

    m.mcell_init_simulation(world)
    m.mcell_init_output(world)
    m.mcell_run_simulation(world)

if __name__ == "__main__":
    main()
