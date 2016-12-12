#!/usr/bin/env python

import pymcell as m
import torus

import os

import ply.lex as lex
import ply.yacc as yacc

##########
# Tokens
##########

# List of tokens
tokens = (
	'UNION',
	'SUBTRACTION',
	'INTERSECTION',
	'LPAREN',
	'RPAREN',
	'LBRACK',
	'RBRACK',
	'VAR'
	)

# Regular expression rules
t_UNION = r'\+'
t_SUBTRACTION = r'\-'
t_INTERSECTION = r'\*'
t_LPAREN  = r'\('
t_RPAREN  = r'\)'
t_LBRACK  = r'\['
t_RBRACK  = r'\]'
t_VAR = r'[a-zA-Z_][a-zA-Z_0-9]*'

###
# Ignored characters
###

# Space
t_ignore = " \t"

# Newline
def t_newline(t):
	r'\n+'
	t.lexer.lineno += t.value.count("\n")

# Error function
def t_error(t):
	print("Illegal character '%s'" % t.value[0])
	t.lexer.skip(1)

##########
# Parser
##########

# Release region expression
# CAUTION:
# this function HAS to go ABOVE p_existing_reg et al. !
def p_release_region_expr(t):
	'''release_region_expr : existing_region
						   | LPAREN release_region_expr RPAREN
						   | release_region_expr UNION release_region_expr
						   | release_region_expr SUBTRACTION release_region_expr
						   | release_region_expr INTERSECTION release_region_expr'''                       
	if len(t) == 2:
		print("Called p_release_region_expr: case 1: Defined release_region_expr from existing_region")
		tmp_sym = m.sym_entry()
		tmp_sym.value = t[1]
		t[0] = m.new_release_region_expr_term(tmp_sym)
		t.lexer.release_evaluator = t[0]
	elif t[1] == '(' and t[3] == ')':
		print("Called p_release_region_expr: case 2: Defined release_region_expr from (release_region_expr)")
		t[0] = t[2]
	elif t[2] == '+':
		print("Called p_release_region_expr: case 3: Defined release_region_expr from union")
		t[0] = m.new_release_region_expr_binary(t[1],t[3],0x02)
		t.lexer.release_evaluator = t[0]
	elif t[2] == '-':
		print("Called p_release_region_expr: case 4: Defined release_region_expr from subtraction")
		t[0] = m.new_release_region_expr_binary(t[1],t[3],0x08)
		t.lexer.release_evaluator = t[0]
	elif t[2] == '*':
		print("Called p_release_region_expr: case 5: Defined release_region_expr from intersection")
		t[0] = m.new_release_region_expr_binary(t[1],t[3],0x04)
		t.lexer.release_evaluator = t[0]

# Existing region
# CAUTION:
# this function HAS to go ABOVE p_var AND ABOVE p_existing_obj !
def p_existing_reg(t):
	'existing_region : existing_object LBRACK var RBRACK'
	print("Called p_existing_reg")
	t[0] = m.existing_region(t.lexer.mcell_world, t[1].sym, t[3]).value

# Existing object
# CAUTION:
# this function HAS to go ABOVE p_var BUT BELOW p_existing_reg !
def p_existing_obj(t):
	'existing_object : var'
	print("Called p_existing_obj")
	t[0] = t.lexer.obj_dict[t[1]]

# Var
def p_var(t):
	'var : VAR'
	print("Called p_var")
	t[0] = t[1]

# Parsing error
def p_error(t):
	print("Syntax error at '%s'" % t.value)


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



def main(lexer, parser):

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

	##########
	# Create each of the three binary operations
	# Each example consists of two box objects
	##########

	# Dictionary of release evaluator objects
	rel_eval_dict = {}

	##########
	# Union
	##########

	# Create box objects
	box_union_outer_name = "Box_Union_Outer"
	r_union_outer = (0.1,0.1,0.1)
	ctr_union_outer = (-0.3,0.3,0)
	box_union_outer_mesh = create_rectangle(world, scene, r_union_outer, ctr_union_outer, box_union_outer_name)

	box_union_inner_name = "Box_Union_Inner"
	r_union_inner = (0.1,0.1,0.12)
	ctr_union_inner = (-0.3+0.1,0.3+0.1,0.0)
	box_union_inner_mesh = create_rectangle(world, scene, r_union_inner, ctr_union_inner, box_union_inner_name)

	# Add to the dictionary of objects for the parser to see
	lexer.obj_dict[box_union_outer_name] = box_union_outer_mesh
	lexer.obj_dict[box_union_inner_name] = box_union_inner_mesh

	# We will use the "ALL" region specification
	# We could also define a surface region using the following code
	'''
	box_union_outer_region_name = "Box_Union_Outer_Reg"
	box_union_outer_face_list = [0,1,2,3,4,5,6,7,8,9,10,11]
	box_union_outer_region = m.create_surface_region(world, box_union_outer_mesh, box_union_outer_face_list, box_union_outer_region_name)
	box_union_inner_region_name = "Box_Union_Inner_Reg"
	box_union_inner_face_list = [0,1,2,3,4,5,6,7,8,9,10,11]
	box_union_inner_region = m.create_surface_region(world, box_union_inner_mesh, box_union_inner_face_list, box_union_inner_region_name)
	'''

	# Parse a string to create a release_evaluator
	s_union = box_union_outer_name+"[ALL]" + " + " + box_union_inner_name+"[ALL]"
	print("> Parsing: " + s_union)
	lexer.mcell_world = world
	parser.parse(s_union)
	rel_eval_dict[s_union] = lexer.release_evaluator
	print("> Finished parsing.")

	##########
	# Subtraction
	##########

	# Create box objects
	box_sub_outer_name = "Box_Sub_Outer"
	r_sub_outer = (0.1,0.1,0.1)
	ctr_sub_outer = (0.3,0.3,0)
	box_sub_outer_mesh = create_rectangle(world, scene, r_sub_outer, ctr_sub_outer, box_sub_outer_name)

	box_sub_inner_name = "Box_Sub_Inner"
	r_sub_inner = (0.1,0.1,0.12)
	ctr_sub_inner = (0.3+0.1,0.3+0.1,0.0)
	box_sub_inner_mesh = create_rectangle(world, scene, r_sub_inner, ctr_sub_inner, box_sub_inner_name)

	# Add to the dictionary of objects for the parser to see
	lexer.obj_dict[box_sub_outer_name] = box_sub_outer_mesh
	lexer.obj_dict[box_sub_inner_name] = box_sub_inner_mesh

	# Parse a string to create a release_evaluator
	s_sub = box_sub_outer_name+"[ALL]" + " - " + box_sub_inner_name+"[ALL]"
	print("> Parsing: " + s_sub)
	lexer.mcell_world = world
	parser.parse(s_sub)
	rel_eval_dict[s_sub] = lexer.release_evaluator
	print("> Finished parsing.")

	##########
	# Intersection
	##########

	# Create box objects
	box_inter_outer_name = "Box_Inter_Outer"
	r_inter_outer = (0.1,0.1,0.1)
	ctr_inter_outer = (0.0,-0.3,0.0)
	box_inter_outer_mesh = create_rectangle(world, scene, r_inter_outer, ctr_inter_outer, box_inter_outer_name)

	box_inter_inner_name = "Box_Inter_Inner"
	r_inter_inner = (0.1,0.1,0.12)
	ctr_inter_inner = (0.1,-0.3+0.1,0.0)
	box_inter_inner_mesh = create_rectangle(world, scene, r_inter_inner, ctr_inter_inner, box_inter_inner_name)

	# Add to the dictionary of objects for the parser to see
	lexer.obj_dict[box_inter_outer_name] = box_inter_outer_mesh
	lexer.obj_dict[box_inter_inner_name] = box_inter_inner_mesh

	# Parse a string to create a release_evaluator
	s_inter = box_inter_outer_name+"[ALL]" + " * " + box_inter_inner_name+"[ALL]"
	print("> Parsing: " + s_inter)
	lexer.mcell_world = world
	parser.parse(s_inter)
	rel_eval_dict[s_inter] = lexer.release_evaluator
	print("> Finished parsing.")

	##########
	# Set up all the release sites for each: union/subtraction/intersection
	##########

	# List of mols to release
	mol_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)

	# Union
	rel_object = m.object()
	box_union_outer_rel_object = m.mcell_create_region_release_boolean(
		world, scene, "Box Union Outer Release", mol_list, 1000, 0, 1, None, rel_eval_dict[s_union], rel_object)

	# Subtraction
	rel_object = m.object()
	box_sub_outer_rel_object = m.mcell_create_region_release_boolean(
		world, scene, "Box Subtraction Outer Release", mol_list, 1000, 0, 1, None, rel_eval_dict[s_sub], rel_object)

	# Union
	mol_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
	rel_object = m.object()
	box_inter_outer_rel_object = m.mcell_create_region_release_boolean(
		world, scene, "Box Intersection Outer Release", mol_list, 1000, 0, 1, None, rel_eval_dict[s_inter], rel_object)

	m.mcell_delete_species_list(mol_list)

	##########
	# Create viz data
	##########

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

	# Clear old parsing files
	os.system("rm ./parser.out")
	os.system("rm ./parsetab.py")

	# Build the lexer
	lexer = lex.lex()
	lexer.release_evaluator = None
	lexer.mcell_world = None
	lexer.obj_dict = {}

	# Build the parser
	parser = yacc.yacc()

	# Run pymcell
	main(lexer, parser)


	# Fragments

	# These are code fragments of the functions that the parser actually calls


	# Start the release
	'''
	tmp_sym = m.sym_entry()
	tmp_sym.name = "Box Outer Release"
	tmp_sym.value = box_outer_mesh
	box_outer_rel_object = m.object()
	m.mcell_start_release_site(world,tmp_sym, box_outer_rel_object)
	'''

	# Pass to mcell_set_release_site_geometry_region
	# m.mcell_set_release_site_geometry_region(world, box_outer_rel_obj, box_outer_mesh, lexer.release_evaluator)

	# Finish the release
	# ...