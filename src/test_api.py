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
    diam_value = 0.00999
    diameter.x = diam_value
    diameter.y = diam_value
    diameter.z = diam_value

    mol_list = m.mcell_add_to_species_list(x_mol_sym, False, 0, None)
    rel_object = m.object()
    release_object = m.mcell_create_geometrical_release_site(
            world, scene, "X_releaser", m.SHAPE_SPHERICAL, position, diameter,
            mol_list, 5000, 1, None, rel_object)
    m.mcell_delete_species_list(mol_list)

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
