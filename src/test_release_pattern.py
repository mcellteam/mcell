#!/usr/bin/env python

import pymcell as m

import os

#############
# Helper function to make rectangular boxes
#############

def create_rectangle(world,scene,r,ctr,name):

	xp = ctr[0]+r[0]
	yp = ctr[1]+r[1]
	zp = ctr[2]+r[2]
	xm = ctr[0]-r[0]
	ym = ctr[1]-r[1]
	zm = ctr[2]-r[2]

	verts = m.mcell_add_to_vertex_list(xp, yp, zm, None)
	verts = m.mcell_add_to_vertex_list(xp, ym, zm, verts)
	verts = m.mcell_add_to_vertex_list(xm, ym, zm, verts)
	verts = m.mcell_add_to_vertex_list(xm, yp, zm, verts)
	verts = m.mcell_add_to_vertex_list(xp, yp, zp, verts)
	verts = m.mcell_add_to_vertex_list(xp, ym, zp, verts)
	verts = m.mcell_add_to_vertex_list(xm, ym, zp, verts)
	verts = m.mcell_add_to_vertex_list(xm, yp, zp, verts)

	elems = m.mcell_add_to_connection_list(1, 2, 3, None)
	elems = m.mcell_add_to_connection_list(7, 6, 5, elems)
	elems = m.mcell_add_to_connection_list(0, 4, 5, elems)
	elems = m.mcell_add_to_connection_list(1, 5, 6, elems)
	elems = m.mcell_add_to_connection_list(6, 7, 3, elems)
	elems = m.mcell_add_to_connection_list(0, 3, 7, elems)
	elems = m.mcell_add_to_connection_list(0, 1, 3, elems)
	elems = m.mcell_add_to_connection_list(4, 7, 5, elems)
	elems = m.mcell_add_to_connection_list(1, 0, 5, elems)
	elems = m.mcell_add_to_connection_list(2, 1, 6, elems)
	elems = m.mcell_add_to_connection_list(2, 6, 3, elems)
	elems = m.mcell_add_to_connection_list(4, 0, 7, elems)

	pobj = m.poly_object()
	pobj.obj_name = name
	pobj.vertices = verts
	pobj.num_vert = 8
	pobj.connections = elems
	pobj.num_conn = 12

	mesh_temp = m.object()
	mesh = m.mcell_create_poly_object(world, scene, pobj, mesh_temp)
	
	return mesh


#############
# Main pymcell function
#############



def main():

	# Create the MCell world
	world = m.mcell_create()
	m.mcell_init_state(world)

	# Set iterations, step size
	dt = 1e-5
	iterations = 100
	m.mcell_set_time_step(world, dt)
	m.mcell_set_iterations(world, iterations)

	# Define volume molecule
	vm1_sym = m.create_species(world, "vm1", 1e-6, False)

	# Create a scene
	scene_name = "Scene"
	scene = m.create_instance_object(world, scene_name)

	# Create release pattern with a delay
	rel_pattern = m.create_release_pattern(world,"Release Pattern",delay=5e-4)

	# Create box objects
	box_name = "Box_Union_Outer"
	r = (0.1,0.1,0.1)
	ctr = (-0.3,0.3,0)
	box_mesh = create_rectangle(world, scene, r, ctr, box_name)

	# List of mols to release
	mol_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)

	# Release
	rel_object = m.object()
	box_rel_object = m.mcell_create_region_release(
		world, scene, box_mesh, "Box Release", "ALL", mol_list, 1000, 0, 1, rel_pattern, rel_object)
	m.mcell_delete_species_list(mol_list)

	# Create viz data
	viz_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
	m.mcell_create_viz_output(world, "./viz_data/test", viz_list, 0, iterations, 1)

	m.mcell_init_simulation(world)
	m.mcell_init_output(world)

	output_freq = 10
	for i in range(iterations):
		m.mcell_run_iteration(world, output_freq, 0)
	m.mcell_flush_data(world)
	m.mcell_print_final_warnings(world)
	m.mcell_print_final_statistics(world)

##########
# Main
##########


if __name__ == "__main__":

	# Run pymcell
	main()

