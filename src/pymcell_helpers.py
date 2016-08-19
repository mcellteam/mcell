def create_count(self, world, where, mol_sym, file_path):
    report_flags = self.REPORT_CONTENTS
    c_list = self.output_column_list()
    # XXX: self.ORIENT_NOT_SET is using -100 instead of SHRT_MIN (used typemap
    # for mcell_create_count in mcell_react_out.i) because limits.h does not
    # work well with swig
    count_list = self.mcell_create_count(
        world, mol_sym, self.ORIENT_NOT_SET, where, report_flags, None, c_list)

    os = self.output_set()
    os = self.mcell_create_new_output_set(
        None, 0, count_list.column_head, self.FILE_SUBSTITUTE, file_path)

    out_times = self.output_times_inlist()
    out_times.type = self.OUTPUT_BY_STEP
    out_times.step = 1e-5

    output = self.output_set_list()
    output.set_head = os
    output.set_tail = os

    self.mcell_add_reaction_output_block(world, output, 10000, out_times)

    return (count_list, os, out_times, output)


def create_species(self, world, name, D, is_2d):
    species_def = self.mcell_species_spec()
    species_def.name = name
    species_def.D = D
    is_2d = 1 if is_2d else 0
    species_def.is_2d = is_2d
    species_def.custom_time_step = 0
    species_def.target_only = 0
    species_def.max_step_length = 0

    species_temp_sym = self.mcell_symbol()
    species_sym = self.mcell_create_species(
        world, species_def, species_temp_sym)

    return species_sym


def create_reaction(
        self, world, reactants, products, rate_constant,
        backward_rate_constant=None, surf_class=None, name=None):

    if surf_class:
        # Do nothing, surf_class has been added and a null object is not needed
        pass
    else:
        surf_class = self.mcell_add_to_species_list(None, False, 0, None)

    arrow = self.reaction_arrow()
    # reversible reaction e.g. A<->B
    if backward_rate_constant:
        arrow.flags = self.ARROW_BIDIRECTIONAL
        rate_constant = self.mcell_create_reaction_rates(
            self.RATE_CONSTANT, rate_constant, self.RATE_CONSTANT,
            backward_rate_constant)
    # irreversible reaction e.g. A->B
    else:
        arrow.flags = self.REGULAR_ARROW
        rate_constant = self.mcell_create_reaction_rates(
            self.RATE_CONSTANT, rate_constant, self.RATE_UNSET, 0)
    arrow.catalyst = self.mcell_species()
    arrow.catalyst.next = None
    arrow.catalyst.mol_type = None
    arrow.catalyst.orient_set = 0
    arrow.catalyst.orient = 0

    if (name):
        name_sym = self.mcell_new_rxn_pathname(world, name)
    else:
        name_sym = None
    self.mcell_add_reaction_simplified(
        world, reactants, arrow, surf_class, products, rate_constant, name_sym)


def create_instance_object(self, world, name):
    scene_temp = self.object()
    return self.mcell_create_instance_object(world, name, scene_temp)


def create_surf_class(self, world, name):
    sc_temp = self.mcell_symbol()
    return self.mcell_create_surf_class(world, name, sc_temp)


def create_release_site(
        self, world, scene, pos, diam, shape, number, mol_sym, name):
    position = self.vector3()
    position.x = pos.x
    position.y = pos.y
    position.z = pos.z
    diameter = self.vector3()
    diameter.x = diam.x
    diameter.y = diam.y
    diameter.z = diam.z

    mol_list = self.mcell_add_to_species_list(mol_sym, False, 0, None)
    rel_object = self.object()
    release_object = self.mcell_create_geometrical_release_site(
        world, scene, name, shape, position, diameter, mol_list, number, 1,
        None, rel_object)
    self.mcell_delete_species_list(mol_list)

    return (position, diameter, release_object)

def create_region_release_site(
        self, world, scene, mesh, release_name, reg_name, number, mol_sym):

    mol_list = self.mcell_add_to_species_list(mol_sym, False, 0, None)
    rel_object = self.object()
    release_object = self.mcell_create_region_release(
        world, scene, mesh, release_name, reg_name, mol_list, number, 1, None,
        rel_object)
    self.mcell_delete_species_list(mol_list)

    return release_object

def create_box_verts_elems(self, half_length):
    hl = half_length
    verts = self.mcell_add_to_vertex_list(hl, hl, -hl, None)
    verts = self.mcell_add_to_vertex_list(hl, -hl, -hl, verts)
    verts = self.mcell_add_to_vertex_list(-hl, -hl, -hl, verts)
    verts = self.mcell_add_to_vertex_list(-hl, hl, -hl, verts)
    verts = self.mcell_add_to_vertex_list(hl, hl, hl, verts)
    verts = self.mcell_add_to_vertex_list(hl, -hl, hl, verts)
    verts = self.mcell_add_to_vertex_list(-hl, -hl, hl, verts)
    verts = self.mcell_add_to_vertex_list(-hl, hl, hl, verts)

    elems = self.mcell_add_to_connection_list(1, 2, 3, None)
    elems = self.mcell_add_to_connection_list(7, 6, 5, elems)
    elems = self.mcell_add_to_connection_list(0, 4, 5, elems)
    elems = self.mcell_add_to_connection_list(1, 5, 6, elems)
    elems = self.mcell_add_to_connection_list(6, 7, 3, elems)
    elems = self.mcell_add_to_connection_list(0, 3, 7, elems)
    elems = self.mcell_add_to_connection_list(0, 1, 3, elems)
    elems = self.mcell_add_to_connection_list(4, 7, 5, elems)
    elems = self.mcell_add_to_connection_list(1, 0, 5, elems)
    elems = self.mcell_add_to_connection_list(2, 1, 6, elems)
    elems = self.mcell_add_to_connection_list(2, 6, 3, elems)
    elems = self.mcell_add_to_connection_list(4, 0, 7, elems)

    return (verts, elems)


def create_polygon_object(self, world, vert_list, face_list, scene, name):
    verts = None
    for x, y, z in vert_list:
        verts = self.mcell_add_to_vertex_list(x, y, z, verts)

    elems = None
    for x, y, z in face_list:
        elems = self.mcell_add_to_connection_list(x, y, z, elems)

    pobj = self.poly_object()
    pobj.obj_name = name
    pobj.vertices = verts
    pobj.num_vert = len(vert_list)
    pobj.connections = elems
    pobj.num_conn = len(face_list)

    mesh_temp = self.object()
    mesh = self.mcell_create_poly_object(world, scene, pobj, mesh_temp)

    return mesh


def create_surface_region(self, world, mesh, surf_reg_face_list, region_name):

    surface_region = self.mcell_create_region(world, mesh, region_name)

    surf_reg_faces = None
    for idx in surf_reg_face_list:
        surf_reg_faces = self.mcell_add_to_region_list(surf_reg_faces, idx)

    self.mcell_set_region_elements(surface_region, surf_reg_faces, 1)

    return surface_region
