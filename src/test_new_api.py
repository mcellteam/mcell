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

    world.release_into_obj(torus_obj, vm1, 1000)
    world.release_into_obj(torus_obj, vm2, 1000)

    world.add_viz((vm1, vm2, vm3))

    # Set partitions
    world.add_partitions('x', -1.3, 1.3, 0.05)
    world.add_partitions('y', -1.3, 1.3, 0.05)
    world.add_partitions('z', -0.275, 0.275, 0.05)

    world.run_sim()

if __name__ == "__main__":
    main()
