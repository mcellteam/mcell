#!/usr/bin/env python

import pymcell as m
import torus


def main():
    iterations = 100
    world = m.MCellSim()
    world.set_time_step(1e-5)
    world.set_iterations(iterations)
    world.set_output_freq(1)

    # define species
    vm1 = m.Species("vm1", 1e-6)
    vm2 = m.Species("vm2", 1e-6)
    vm3 = m.Species("vm3", 1e-6)

    # define reaction
    rxn = m.Reaction((vm1, vm2), (vm3, ), 1e8)
    world.add_reaction(rxn)

    # create torus object
    torus_obj = m.PolyObj("Torus", "half", torus.vert_list,
                          torus.face_list, torus.surf_reg_face_list)
    world.add_geometry(torus_obj)

    # release molecules into torus
    world.release_into_obj(torus_obj, vm1, 1000)
    world.release_into_obj(torus_obj, vm2, 1000)

    # viz and reaction data
    world.add_viz((vm1, vm2, vm3))
    world.add_count(vm1, torus_obj)
    world.add_count(vm2, torus_obj)
    world.add_count(vm3, torus_obj)

    # set partitions
    world.add_partitions('x', -1.3, 1.3, 0.05)
    world.add_partitions('y', -1.3, 1.3, 0.05)
    world.add_partitions('z', -0.275, 0.275, 0.05)

    # run the simulation! :)
    # for i in range(iterations+1):
    #     world.run_iteration()
    # world.end_sim()
    world.run_sim()


if __name__ == "__main__":
    main()
