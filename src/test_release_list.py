#!/usr/bin/env python

import pymcell as m

def main():
    world = m.mcell_create()
    m.mcell_init_state(world)

    dt = 1e-5
    iterations = 5
    m.mcell_set_time_step(world, dt)
    m.mcell_set_iterations(world, iterations)

    # Define one volume molecule
    vm1_sym = m.create_species(world, "vm1", 1e-6, False)

    scene_name = "Scene"
    scene = m.create_instance_object(world, scene_name)

    rel_obj = m.create_list_release_site(world,scene,[vm1_sym,vm1_sym],[1.0,0.0],[0.0,0.0],[0.0,0.0],2,"rel_name")

    # Create viz data
    viz_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
    m.mcell_create_viz_output(
        world, "./test_release_list_viz/test", viz_list, 0, iterations, 1)

    m.mcell_init_simulation(world)
    m.mcell_init_output(world)

    output_freq = 1
    for i in range(iterations):
        m.mcell_run_iteration(world, output_freq, 0)
    m.mcell_flush_data(world)
    m.mcell_print_final_warnings(world)
    m.mcell_print_final_statistics(world)

if __name__ == "__main__":
    main()
