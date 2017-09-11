#!/usr/bin/env python

import pymcell as m
from pymcell import SC
import torus as t
import box as b
import logging


def main():
    world = m.MCellSim(seed=1)
    world.set_time_step(time_step=1e-5)
    world.set_iterations(iterations=100)
    world.enable_logging()
    world.silence_notifications()
    world.set_output_freq(1)

    slow_dc = 1e-9
    fast_dc = 1e-6
    # define species
    vm1 = m.Species("vm1", fast_dc)
    vm2 = m.Species("vm2", fast_dc)
    vm3 = m.Species("vm3", slow_dc)
    sm_list = []
    for i in range(1, 11):
        sm_list.append(m.Species("sm{}".format(i), slow_dc, surface=True))
    sm1 = sm_list[0]
    spec_orient_list = [sm.up() for sm in sm_list]
    all_mols = sm_list + [vm1, vm2, vm3]
    world.add_species(all_mols)

    # define reaction
    rxn = m.Reaction((vm1.down(), sm1.up()), vm2.down(), 1e8, name="create_vm2")
    world.add_reaction(rxn)

    # create torus object
    torus_obj = m.MeshObj(
        "Torus", t.vert_list, t.face_list, translation=(0, 0, 0))
    torus_reg = m.SurfaceRegion(torus_obj, 'half', t.surf_reg_face_list)
    world.add_geometry(torus_obj)

    # create box object
    box_obj = m.MeshObj(
        "Box", b.vert_list, b.face_list, translation=(0, 0, 0))
    box_reg = m.SurfaceRegion(box_obj, 'top', b.surf_reg_face_list)
    world.add_geometry(box_obj)

    # create surface class to absorb vm1
    sc = m.SurfaceClass(SC.absorb, sm1.mix())
    world.assign_surf_class(sc, torus_reg)

    # release sm1 molecules into torus
    sm1_torus_rel = m.ObjectRelease(sm1.up(), number=1000, region=torus_reg)
    world.release(sm1_torus_rel)
    vm1_torus_rel = m.ObjectRelease(vm1, number=1000, mesh_obj=torus_obj)
    world.release(vm1_torus_rel)

    # release vm3 molecules in a line (list based release)
    vm3_pos = [(x*0.01, 0, 0) for x in range(-5, 5)]
    vm3_list_release = m.ListRelease(vm3, vm3_pos)
    world.release(vm3_list_release)

    # release sm molecules in a line on top of box (list based release)
    sm_pos = [(x*0.01, 0, 0.1) for x in range(-5, 5)]
    sm_list_release = m.ListRelease(spec_orient_list, sm_pos)
    world.release(sm_list_release)
    # pos_vec3 = m.Vector3()
    # diam_vec3 = m.Vector3(0.015, 0.015, 0.015)
    # world.create_release_site(vm1, 100, "spherical", diam_vec3=diam_vec3)

    # viz and reaction data
    world.add_viz(all_mols)
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
