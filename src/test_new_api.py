#!/usr/bin/env python

import pymcell as m
from pymcell import SC
import torus
import logging


def main():
    world = m.MCellSim(seed=1)
    world.set_time_step(time_step=1e-5)
    world.set_iterations(iterations=100)
    world.enable_logging()
    # world.silence_notifications()
    world.set_output_freq(1)

    # define species
    vm1 = m.Species("vm1", 1e-6)
    sm1 = m.Species("sm1", 1e-6, surface=True)
    vm2 = m.Species("vm2", 1e-6)
    world.add_species((vm1, sm1, vm2))

    # define reaction
    rxn = m.Reaction((vm1.down(), sm1.up()), (vm2.down(), ), 1e8, name="create_vm2")
    world.add_reaction(rxn)

    # create torus object
    torus_obj = m.MeshObj("Torus", torus.vert_list, torus.face_list)
    torus_reg = m.SurfaceRegion(torus_obj, 'half', torus.surf_reg_face_list)
    world.add_geometry(torus_obj)

    # Create surface class to absorb vm1
    sc = m.SurfaceClass(SC.absorb, sm1.mix())
    world.assign_surf_class(sc, torus_reg)

    # release molecules into torus
    sm1_torus_rel = m.ObjectRelease(sm1.up(), number=1000, region=torus_reg)
    world.release(sm1_torus_rel)
    vm1_torus_rel = m.ObjectRelease(vm1, number=1000, mesh_obj=torus_obj)
    world.release(vm1_torus_rel)
    # pos_vec3 = m.Vector3()
    # diam_vec3 = m.Vector3(0.015, 0.015, 0.015)
    # world.create_release_site(vm1, 100, "spherical", diam_vec3=diam_vec3)

    # viz and reaction data
    world.add_viz((vm1, vm2, sm1))
    world.add_count(vm1, torus_obj)
    world.add_count(vm2, torus_obj)
    world.add_count(sm1, torus_obj)

    # set partitions
    world.add_partitions('x', -1.3, 1.3, 0.05)
    world.add_partitions('y', -1.3, 1.3, 0.05)
    world.add_partitions('z', -0.275, 0.275, 0.05)

    # run the simulation! :)
    # for i in range(iterations+1):
    #     vm3_count = world.get_species_count(vm3, torus_obj)
    #     if (vm3_count > 50):
    #         world.modify_rate_constant(rxn, 1e12)
    #     world.run_iteration()
    world.run_sim()
    world.end_sim()


if __name__ == "__main__":
    main()
