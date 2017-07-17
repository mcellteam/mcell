#!/usr/bin/env python

import pymcell as m
import torus


def main():
    iterations = 100
    time_step = 1e-5
    seed = 1
    world = m.MCellSim(seed)
    world.set_time_step(time_step)
    world.set_iterations(iterations)
    world.set_output_freq(1)

    # define species
    vm1 = m.Species("vm1", 1e-6)
    vm2 = m.Species("vm2", 1e-6)
    vm3 = m.Species("vm3", 1e-6)

    # define reaction
    rxn = m.Reaction((vm1, vm2), (vm3, ), 1e8, "create_vm3")
    world.add_reaction(rxn)

    # create torus object
    torus_obj = m.PolyObj("Torus", torus.vert_list, torus.face_list)
    torus_reg = m.SurfaceRegion(torus_obj, 'half', torus.surf_reg_face_list)
    world.add_geometry(torus_obj)

    # Create surface class to absorb vm1
    # sc = m.SurfaceClass('sc', 'absorptive', vm1)
    # world.assign_surf_class(sc, torus_reg)

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
    #     vm3_count = world.get_molecule_count(vm3, torus_obj)
    #     if (vm3_count > 50):
    #         world.modify_rate_constant(rxn, 1e12)
    #     world.run_iteration()
    world.run_sim()
    world.end_sim()


if __name__ == "__main__":
    main()
