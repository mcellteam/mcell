#!/usr/bin/env python

import pymcell as m
import torus


def main():
    world = m.MCellSim()
    world.set_time_step(1e-5)
    world.set_iterations(100)

    # define species
    vm1 = m.Species("vm1", 1e-6)
    vm2 = m.Species("vm2", 1e-6)
    vm3 = m.Species("vm3", 1e-6)

    # define reaction
    rxn = m.Reaction((vm1, vm2), (vm3, ), 1e6)
    world.add_reaction(rxn)

    # create torus object
    torus_obj = m.PolyObj("Torus", "half", torus.vert_list,
                          torus.face_list, torus.surf_reg_face_list)
    world.add_geometry(torus_obj)

    world.release_into_obj(torus_obj, vm1)

    world.add_viz((vm1, vm2, vm3))

    world.run_sim()

    # # Define one surface molecule and three volume molecules
    # sm1_sym = m.create_species(world, "sm1", 0, True)
    # vm1_sym = m.create_species(world, "vm1", 1e-6, False)
    # vm2_sym = m.create_species(world, "vm2", 1e-6, False)
    # vm3_sym = m.create_species(world, "vm3", 0, False)

    # # This is the world object. We just call it "Scene" here to be consistent
    # # with the MDL output from Blender.
    # scene_name = "Scene"
    # scene = m.create_instance_object(world, scene_name)

    # # Set partitions
    # m.create_partitions(world, 0, -1.3, 1.3, 0.05)
    # m.create_partitions(world, 1, -1.3, 1.3, 0.05)
    # m.create_partitions(world, 2, -0.275, 0.275, 0.05)

    # m.mcell_init_simulation(world)
    # m.mcell_init_output(world)

    # output_freq = 10
    # torus_obj = m.PolyObj("Torus", "half_torus", torus.vert_list,
    #                       torus.face_list, torus.surf_reg_face_list)
    # box_obj = m.PolyObj("Box", "side", box.vert_list,
    #                     box.face_list, box.surf_reg_face_list)
    # obj_list = [torus_obj, box_obj]
    # for i in range(iterations+1):
    #     m.mcell_run_iteration(world, output_freq, 0)
    # m.mcell_flush_data(world)
    # m.mcell_print_final_warnings(world)
    # m.mcell_print_final_statistics(world)

if __name__ == "__main__":
    main()
