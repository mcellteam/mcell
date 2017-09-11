"""This is a simple example which shows mcells ability to communicate
with neuron through a python interface. This example is simple because
when they communicate with eachother in the simulation it more resembles
two people shouting over vhf to check that their radios work than the
passing of information between two timescales in a biological simulation

"""


import math

import neuron
import pymcell as m
import torus


# Define some functions fo the iteraction of mcell and neuron

# Main
if __name__ == "__main__":

        # Create pyMCell instance
        world = m.mcell_create()
        m.mcell_init_state(world)

        dt = 1e-6
        iterations = 300
        m.mcell_set_time_step(world, dt)
        m.mcell_set_iterations(world, iterations)

        # Define one surface molecule and three volume molecules
        vm1_sym = m.create_species(world, "vm1", 1e-6, False)

        # vm1 -> NULL [1e5]
        reactants2 = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
        m.create_reaction(world, reactants2, None, 0.01, name="rxn")

        # Create Scene for simulation
        scene_name = "Scene"
        scene = m.create_instance_object(world, scene_name)

        # Create a spherical release site
        pos_vec3 = m.Vector3()
        diam_vec3 = m.Vector3(0.015, 0.015, 0.015)

        position, diameter, sphere_release_object = m.create_release_site(
            world, scene, pos_vec3, diam_vec3, m.SHAPE_SPHERICAL, 500, 0, vm1_sym,
            "vm1_rel")

        obj_name = "Torus"
        mesh = m.create_polygon_object(
            world, torus.vert_list, torus.face_list, scene, obj_name)

        region_release_object = m.create_region_release_site(
            world, scene, mesh, "vm1_torus_rel", "ALL", 1000, 0, vm1_sym)

        # Create viz data
        viz_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
        m.mcell_create_viz_output(
            world, "./viz_data/test", viz_list, 0, iterations, 1)

        # Create reaction data
        mesh_sym = m.mcell_get_obj_sym(mesh)
        count_list1, os1, out_times1, output1 = m.create_count(
            world, mesh_sym, vm1_sym, "react_data/vm1_%s.dat" % obj_name)

        # Initialize simulation
        m.mcell_init_simulation(world)
        m.mcell_init_output(world)

        # Define soma compartment with HH mechanisms
        soma = neuron.h.Section(name="soma")
        soma.nseg = 1
        soma.diam = 10
        soma.L = 10
        soma.insert("hh")
        meca = soma(.5).hh

        # Record time from NEURON
        rec_t = neuron.h.Vector()
        rec_t.record(neuron.h._ref_t)
        #record voltage from center of soma
        rec_v = neuron.h.Vector()
        rec_v.record(soma(.5)._ref_v)

        # Set up stimulus
        stim = neuron.h.IClamp(soma(.5))

        stim.delay = 100
        stim.dur = 100
        stim.amp = 20

        # Initialize and Run
        neuron.h.finitialize(-65)
        tstop = 7
        neuron.h.dt = .025

        #neuron.run(tstop)
        t = 0.0
        output_freq = 10
        r_rate = 20

        while t < tstop:

                # Advance neuron
                print("-"*40)
                vm1_count = m.mcell_get_count(
                    "vm1", "%s.%s,ALL" % (scene_name, obj_name), world)
                print("vm1_count: %d" % vm1_count)
                soma.gnabar_hh = vm1_count*100
                print("soma.gnabar_hh: %g" % soma.gnabar_hh)
                neuron.h.fadvance()
                t += neuron.h.dt

                # Advance the MCell main model to keep up
                if t > 3:
                        r_rate = rec_v[math.ceil(t)]*100
                print("r_rate is %g" % r_rate)
                m.mcell_modify_rate_constant(world, "rxn", r_rate)
                m.mcell_run_iteration(world, output_freq, 0)

        # Plot and print the results
        time = rec_t.to_python()
        voltage = rec_v.to_python()

        # print(voltage)

        m.mcell_flush_data(world)
