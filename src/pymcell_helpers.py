"""pyMCell helper functions.

This defines functions to help pyMCell interface with the user. It
combines calls of the low level swig wrapped c code so that the
resulting python function calls for users using pyMCell are more
intuitive.

"""

import pymcell as m
from typing import List


class PolyObj(object):
    """ An entire polygon object and its associated surface regions. """
    def __init__(
            self, obj_name: str, vert_list: List, face_list: List) -> None:
        self.obj_name = obj_name
        self.vert_list = vert_list
        self.face_list = face_list
        self.regions = [] # type: List[SurfaceRegion]


class SurfaceRegion(object):
    """ Subsets of a surface.
    Examples of uses: molecules can be released on to these and surface classes
    can be assigned to them.
    """
    def __init__(self, obj: PolyObj,
                 reg_name: str,
                 surf_reg_face_list: List) -> None:
        self.reg_name = reg_name
        self.full_reg_name = "%s[%s]" % (obj.obj_name, reg_name)
        self.surf_reg_face_list = surf_reg_face_list
        # Relationship is symmetrical. Every region knows its object. Object
        # knows its regions.
        self.obj = obj
        obj.regions.append(self)


class Vector3(object):
    """ Just a generic 3d  vector to be used for positions and whatnot. """
    def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> None:
        self.x = x
        self.y = y
        self.z = z


class SurfaceClass(object):
    """ These describe how molecules interact with various surfaces/meshes.
    ex: X molecules are absorbed when they hit the front of a surface.
    """
    def __init__(
            self, sc_name: str, sc_type: str, species: Species,
            orient: int = 0) -> None:
        self.sc_name = sc_name
        self.sc_type = sc_type
        self.species = species
        self.orient = orient


class MCellSim(object):
    """ Everything needed to run a pyMCell simulation. """
    def __init__(self, seed: int):
        self._world = m.mcell_create()
        self._started = False
        self._species = {}
        self._surface_classes = {}
        self._objects = {}
        self._regions = {}
        self._releases = {}
        self._counts = {}
        self._iterations = 0
        self._current_iteration = 0
        self._finished = False
        self._output_freq = 10
        self._seed = seed
        m.mcell_set_seed(self._world, seed)
        m.mcell_init_state(self._world)
        # This is the top level instance object. We just call it "Scene" here
        # to be consistent with the MDL output from Blender.
        self._scene = m.create_instance_object(self._world, "Scene")

    def __del__(self):
        self.end_sim()

    def set_output_freq(self, output_freq: int) -> None:
        """ How often do we output reaction data. """
        self._output_freq = output_freq

    def set_time_step(self, dt: float) -> None:
        """ Set time step in seconds. """
        m.mcell_set_time_step(self._world, dt)

    def set_iterations(self, iterations: int) -> None:
        """ Set number of iterations """
        m.mcell_set_iterations(self._world, iterations)
        self._iterations = iterations

    # def set_seed(self, seed):
    #     self._seed = seed
    #     m.mcell_set_seed(self._world, seed)

    def add_reaction(self, rxn: Reaction) -> None:
        """ Add a reaction object. """
        r_spec_list = None
        p_spec_list = None
        for r in rxn.reactants:
            r_sym = m.create_species(
                self._world, r.name, r.diffusion_constant, r.surface)
            if r.name not in self._species:
                self._species[r.name] = r_sym
            r_spec_list = m.mcell_add_to_species_list(
                r_sym, False, 0, r_spec_list)
        for p in rxn.products:
            p_sym = m.create_species(
                self._world, p.name, p.diffusion_constant, p.surface)
            if p.name not in self._species:
                self._species[p.name] = p_sym
            p_spec_list = m.mcell_add_to_species_list(
                p_sym, False, 0, p_spec_list)
        m.create_reaction(
            self._world, r_spec_list, p_spec_list, rxn.rate, name=rxn.name)

    def add_geometry(self, geom: PolyObj) -> None:
        """ Add a polygon object. """
        mesh = m.create_polygon_object(
            self._world,
            geom.vert_list,
            geom.face_list,
            self._scene,
            geom.obj_name)
        if geom.obj_name not in self._objects:
            self._objects[geom.obj_name] = mesh
        for reg in geom.regions:
            region_swig_obj = m.create_surface_region(
                self._world, mesh, reg.surf_reg_face_list, reg.reg_name)
            full_reg_name = "%s[%s]" % (geom.obj_name, reg.reg_name)
            if reg.reg_name not in self._regions:
                self._regions[full_reg_name] = region_swig_obj

    def add_viz(self, species: Species) -> None:
        """ Add this to the list of species to be visualized. """
        viz_list = None
        for spec in species:
            viz_list = m.mcell_add_to_species_list(
                self._species[spec.name], False, 0, viz_list)
        m.mcell_create_viz_output(
            self._world, "./viz_data/seed_%04i/Scene" % self._seed, viz_list,
            0, self._iterations, 1)

    def release_into_obj(
            self, geom: PolyObj, mol: Species, count: int) -> None:
        """ Release the specified amount of molecules into an object. """
        rel_name = "%s_%s_rel" % (mol.name, geom.obj_name)
        release_object = m.create_region_release_site(
            self._world, self._scene, self._objects[geom.obj_name],
            rel_name, "ALL", count, 0, self._species[mol.name])
        if rel_name not in self._releases:
            self._releases[rel_name] = release_object

    def create_release_site(
            self, species: Species, count: int, shape: str,
            pos_vec3: Vector3 = None, diam_vec3=None) -> None:
        """ Create a spherical/cubic release shape. """
        if pos_vec3 is None:
            pos_vec3 = m.Vector3()
        if diam_vec3 is None:
            diam_vec3 = m.Vector3()
        species_sym = self._species[species.name]
        rel_name = "%s_%s_rel" % (species.name, shape)
        if shape == "spherical":
            shape = m.SHAPE_SPHERICAL
        elif shape == "cubic":
            shape = m.SHAPE_CUBIC
        position, diameter, release_object = m.create_release_site(
            self._world, self._scene, pos_vec3, diam_vec3, shape,
            count, 0, species_sym, rel_name)
        if rel_name not in self._releases:
            self._releases[rel_name] = (position, diameter, release_object)

    def add_partitions(
            self, axis: str, start: float, stop: float, step: float) -> None:
        """ Add partitions to speed up the simulation. """
        if axis == "x":
            axis_num = 0
        elif axis == "y":
            axis_num = 1
        elif axis == "z":
            axis_num = 2
        m.create_partitions(self._world, axis_num, start, stop, step)

    def add_count(self, species: Species, geom: PolyObj) -> None:
        """ Add this to the list of species to be counted. """
        species_sym = self._species[species.name]
        mesh = self._objects[geom.obj_name]
        mesh_sym = m.mcell_get_obj_sym(mesh)
        count_str = "react_data/seed_%04d/%s_%s" % (
                self._seed, species.name, geom.obj_name)
        count_list, os, out_times, output = m.create_count(
            self._world, mesh_sym, species_sym, count_str, 1e-5)
        self._counts[count_str] = (count_list, os, out_times, output)

    def assign_surf_class(
            self, sc: SurfaceClass, region: SurfaceRegion) -> None:
        """ Assign a surface class to a region. """

        if sc.sc_type == "reflect":
            sc_type = m.RFLCT
        elif sc.sc_type == "transparent":
            sc_type = m.TRANSP
        elif sc.sc_type == "absorptive":
            sc_type = m.SINK
        else:
            print("Surface class '%s' is not a valid option" % sc_type)
        sc_sym = m.create_surf_class(self._world, sc.sc_name)
        self._surface_classes[sc.sc_name] = sc_sym
        spec_sym = self._species[sc.species.name]
        m.mcell_add_surf_class_properties(
            self._world, sc_type, sc_sym, spec_sym, sc.orient)
        region_swig_obj = self._regions[region.full_reg_name]
        m.mcell_assign_surf_class_to_region(sc_sym, region_swig_obj)

    def get_molecule_count(self, molecule: Species, geom: PolyObj) -> int:
        """ Get the current count of a molecule. """
        return m.mcell_get_count(
            molecule.name, "Scene.%s,ALL" % geom.obj_name, self._world)

    def modify_rate_constant(
            self, rxn: Reaction, new_rate_constant: float) -> None:
        """ Modify the rate constant of the specified reaction. """
        if not rxn.name:
            print("You can only change a named reaction.")
        else:
            m.mcell_modify_rate_constant(
                self._world, rxn.name, new_rate_constant)

    def run_iteration(self) -> None:
        """ Run a single iteration. """
        if self._finished:
            print("The simulation is done running")
            return

        if not self._started:
            m.mcell_init_simulation(self._world)
            m.mcell_init_output(self._world)
            self._started = True
        if self._current_iteration <= self._iterations:
            m.mcell_run_iteration(self._world, self._output_freq, 0)
        # You have to kill it now or it will hang
        if self._current_iteration == self._iterations:
            m.mcell_flush_data(self._world)
            m.mcell_print_final_warnings(self._world)
            m.mcell_print_final_statistics(self._world)
            self._finished = True
        self._current_iteration += 1

    def run_sim(self) -> None:
        """ Run the entire simulation without interruption. """
        if self._finished:
            print("The simulation is done running")
            return

        m.mcell_init_simulation(self._world)
        m.mcell_init_output(self._world)

        for i in range(self._iterations+1):
            m.mcell_run_iteration(self._world, self._output_freq, 0)
        m.mcell_flush_data(self._world)
        m.mcell_print_final_warnings(self._world)
        m.mcell_print_final_statistics(self._world)
        self._finished = True

    def end_sim(self) -> None:
        # Call this if we end the simulation early
        if self._started and not self._finished:
            m.mcell_flush_data(self._world)
            m.mcell_print_final_warnings(self._world)
            m.mcell_print_final_statistics(self._world)
            self._finished = True


class Species(object):
    """ A type of molecule. """
    def __init__(
            self, name: str, diffusion_constant: float, surface: bool = False):
        self.name = name
        self.diffusion_constant = diffusion_constant
        self.surface = surface


class Reaction(object):
    """ A reaction involving one or more molecules.
    ex: vm1 -> vm2 + vm3
    - Can be unimolecular or bimolecular.
    - Involves surface and/or volume molecules.
    """
    def __init__(self, reactants: List[int], products, rate, name=None):
        self.reactants = reactants
        self.products = products
        self.rate = rate
        self.name = name


def create_partitions(world, axis, start, stop, step):
    expr_list = m.num_expr_list_head()
    expr_list.value_head = None
    expr_list.value_tail = None
    expr_list.value_count = 0
    expr_list.shared = 1
    m.mcell_generate_range(expr_list, start, stop, step)
    expr_list.shared = 1
    m.mcell_set_partition(world, axis, expr_list)


def create_release_pattern(
    world, name, delay=0, release_interval=1e20, train_interval=1e20,
        train_duration=1e20, number_of_trains=1):
    """Create a release pattern

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        name -- name of the release pattern
        other arguments -- as listed
    """

    return m.mcell_create_release_pattern(
        world, name, delay, release_interval, train_interval, train_duration,
        number_of_trains)


def create_count(world, where, mol_sym, file_path, step):
    """Creates a count for a specified molecule in a specified region
    and initializes an output block for the count data that will be
    generated.

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        where (sym_entry) -- symbol entry for the location you want to
            record
        mol_sym (sym_entry) -- symbol entry for the molecule
        file_path (dir) -- name of the file path to output the data to
        step -- frequency of output in seconds

    Returns:
        The return values count list, output set, output times and
        output structure

    """
    report_flags = m.REPORT_CONTENTS
    c_list = m.output_column_list()
    # XXX: m.ORIENT_NOT_SET is using -100 instead of SHRT_MIN (used typemap
    # for mcell_create_count in mcell_react_out.i) because limits.h does not
    # work well with swig
    count_list = m.mcell_create_count(
        world, mol_sym, m.ORIENT_NOT_SET, where, report_flags, None, c_list)

    os = m.output_set()
    os = m.mcell_create_new_output_set(
        None, 0, count_list.column_head, m.FILE_SUBSTITUTE, file_path)

    out_times = m.output_times_inlist()
    out_times.type = m.OUTPUT_BY_STEP
    out_times.step = step

    output = m.output_set_list()
    output.set_head = os
    output.set_tail = os

    m.mcell_add_reaction_output_block(world, output, 10000, out_times)

    return (count_list, os, out_times, output)


def create_species(world, name, D, is_2d, custom_time_step=0):
    """Creates a molecule species

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        name (string) -- Name of the molecule species that will be
            generated
        D (double) -- Diffusion Coefficient for the molecule species
            that will be generated.
        is_2d (bool) -- Boolean describing whether new species is a
            surface molecule
        custom_time_step -- Custom time step (< 0.0 for a custom space step,
                       >0.0 for custom timestep, 0.0 for default timestep)
    Returns:
        (mcell_symbol) Returns a species sym_entry

    """
    species_def = m.mcell_species_spec()
    species_def.name = name
    species_def.D = D
    is_2d = 1 if is_2d else 0
    species_def.is_2d = is_2d
    species_def.custom_time_step = custom_time_step
    species_def.target_only = 0
    species_def.max_step_length = 0

    species_temp_sym = m.mcell_symbol()
    species_sym = m.mcell_create_species(
        world, species_def, species_temp_sym)

    return species_sym


def create_reaction(
        world, reactants, products, rate_constant,
        backward_rate_constant=0.0, surf_class=None, name=None):
    """Creates a molecular reaction

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        reactants (mcell_species_list) -- list of species that are the
            reactants for the reaction
        products (mcell_species_list) -- list of species that are the
            products for the reaction
        rate_constant (double) -- the rate constant for the forward
            direction reaction -> product
        backward_rate_constant (double)(optional) -- the rate constant
            for the backward direction reaction <- product
        surf_class (mcell species list surface class)(optional) -- the
            surface class upon which the reaction will happen
        name (string)(optional) -- Name of the reaction

    Returns:
        void -- creates a reaction, by generating reaction_rates
            structure

    """

    if surf_class:
        # Do nothing, surf_class has been added and a null object is not needed
        pass
    else:
        surf_class = m.mcell_add_to_species_list(None, False, 0, None)

    arrow = m.reaction_arrow()
    # reversible reaction e.g. A<->B
    if backward_rate_constant:
        arrow.flags = m.ARROW_BIDIRECTIONAL
        rate_constant = m.mcell_create_reaction_rates(
            m.RATE_CONSTANT, rate_constant, m.RATE_CONSTANT,
            backward_rate_constant)
    # irreversible reaction e.g. A->B
    else:
        arrow.flags = m.REGULAR_ARROW
        rate_constant = m.mcell_create_reaction_rates(
            m.RATE_CONSTANT, rate_constant, m.RATE_UNSET, 0)
    arrow.catalyst = m.mcell_species()
    arrow.catalyst.next = None
    arrow.catalyst.mol_type = None
    arrow.catalyst.orient_set = 0
    arrow.catalyst.orient = 0

    if (name):
        name_sym = m.mcell_new_rxn_pathname(world, name)
    else:
        name_sym = None
    m.mcell_add_reaction_simplified(
        world, reactants, arrow, surf_class, products, rate_constant, name_sym)


def create_instance_object(world, name):
    """Creates an instance object. Simple translation from wrapped code
    to python function. Frees the user from having to initialize the
    scene object and then pass it in and generate the object.

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        name (string) -- name of the instance object

    Returns:
        instance object

    """
    scene_temp = m.object()
    return m.mcell_create_instance_object(world, name, scene_temp)


def create_surf_class(world, name):
    """Creates a surface class. Simple translation from wrapped code to
    python function Frees the user from having to initialize the surface
    class symbol and then pass it in and generate the object.

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        name (string) -- name of the instance object

    Returns:
        mcell_symbol for surface class

    """

    sc_temp = m.mcell_symbol()
    return m.mcell_create_surf_class(world, name, sc_temp)


def create_list_release_site(
    world, scene, mol_list, xpos, ypos, zpos, name, surf_flags=None,
        orientations=None, diameter=1e-4):
    '''
    Creates a list release site
    All is self explanatory except mol_list:
    This is a list of "mol_sym" that you get back when you create the species.
    This is a Python list - it is converted to a species list in this function
    for you. By default, assumes all molecules are volume molecules. Else, need
    to pass surf_flags=[True,True,False,...] and their orientations =
    [1,0,1,...] Diameter is the diameter we search for to place a surface mol
    It can be None (= NULL in C) but then we do a poor job of placing surface
    mols
    '''

    # Check that they're all the same length
    n = len(mol_list)
    if len(xpos) != n or len(ypos) != n or len(zpos) != n:
        raise ValueError("All lists must have the same length.")

    # Check that if there are surface molecules
    if surf_flags is not None:
        # Check that there are enough
        if len(surf_flags) != n or len(orientations) != n:
            raise ValueError(
                "surf_flags and orientations lists must have the same lengths "
                "as the others.")

    # Convert to floats (can't be int)
    xpos = [float(q) for q in xpos]
    ypos = [float(q) for q in ypos]
    zpos = [float(q) for q in zpos]

    # Diameter
    diam = m.vector3()
    diam.x = diameter
    diam.y = diameter
    diam.z = diameter

    # Mols
    mol_list.reverse()
    species_list = None
    # All volume molecules
    if surf_flags is None:
        for mol_sym in mol_list:
            species_list = m.mcell_add_to_species_list(
                mol_sym, False, 0, species_list)
    else:
        for i, mol_sym in enumerate(mol_list):
            species_list = m.mcell_add_to_species_list(
                mol_sym, surf_flags[i], orientations[i], species_list)

    rel_object = m.object()
    ret = m.mcell_create_list_release_site(
        world, scene, name, species_list, xpos, ypos, zpos, n, diam,
        rel_object)
    # Delete the species list
    m.mcell_delete_species_list(species_list)

    # VERY IMPORTANT HERE - MUST RETURN "ret"
    # If we throw this away, all is lost....
    return (rel_object, ret)


def create_release_site(
        world, scene, pos, diam, shape, number, number_type, mol_sym, name):
    """Creates a molecule release site

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        scene (instance object) -- scene for mcell simulation
        pos (vector3) -- position of release site
        diam (vector3) -- diameter of release site
        number (int or float) -- number to be release at release site
        number_type (int) -- 0 for NUMBER, 1 for CONCENTRATION
        mol_sym (mcell_symbol) -- species to be released
        name (string) -- name of the release site

    Returns:
        void -- generates a species release site

    """

    position = m.vector3()
    position.x = pos.x
    position.y = pos.y
    position.z = pos.z
    diameter = m.vector3()
    diameter.x = diam.x
    diameter.y = diam.y
    diameter.z = diam.z

    mol_list = m.mcell_add_to_species_list(mol_sym, False, 0, None)
    rel_object = m.object()
    release_object = m.mcell_create_geometrical_release_site(
        world, scene, name, shape, position, diameter, mol_list, float(number),
        number_type, 1, None, rel_object)
    m.mcell_delete_species_list(mol_list)

    return (position, diameter, release_object)


def create_region_release_site(
        world, scene, mesh, release_name, reg_name, number, number_type,
        mol_sym):
    """Creates a release site on a specific region

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        scene (instance object) -- scene for mcell simulation
        mesh (mesh object) -- scene object where the release will occur
        release_name (string) -- name of the region release site
        reg_name (string) -- name of the region for the release site
        number (int or float) -- number to be released at the region release
            site
        number_type (int) -- 0 for NUMBER, 1 for CONCENTRATION
        mol_sym (mcell_symbol) -- species to be released

    Returns:
        release object (object)

    """

    mol_list = m.mcell_add_to_species_list(mol_sym, False, 0, None)
    rel_object = m.object()
    release_object = m.mcell_create_region_release(
        world, scene, mesh, release_name, reg_name, mol_list, float(number),
        number_type, 1, None, rel_object)
    m.mcell_delete_species_list(mol_list)

    return release_object


def create_box(world, scene, half_length, name):
    """Creates the verteces and lines of a cube object at the origin

    Args:
        half_length (double) -- half length of the cube

    Returns:
        vertex list and element connection list

    """

    hl = half_length
    verts = m.mcell_add_to_vertex_list(hl, hl, -hl, None)
    verts = m.mcell_add_to_vertex_list(hl, -hl, -hl, verts)
    verts = m.mcell_add_to_vertex_list(-hl, -hl, -hl, verts)
    verts = m.mcell_add_to_vertex_list(-hl, hl, -hl, verts)
    verts = m.mcell_add_to_vertex_list(hl, hl, hl, verts)
    verts = m.mcell_add_to_vertex_list(hl, -hl, hl, verts)
    verts = m.mcell_add_to_vertex_list(-hl, -hl, hl, verts)
    verts = m.mcell_add_to_vertex_list(-hl, hl, hl, verts)

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


def change_geometry(world, scene_name, obj_list):

    pobj_list = None
    verts = None
    elems = None

    for p in obj_list:

        verts = None
        for x, y, z in p.vert_list:
            verts = m.mcell_add_to_vertex_list(x, y, z, verts)

        elems = None
        for x, y, z in p.face_list:
            elems = m.mcell_add_to_connection_list(x, y, z, elems)

        pobj = m.poly_object()
        pobj.obj_name = p.obj_name
        pobj.vertices = verts
        pobj.num_vert = len(p.vert_list)
        pobj.connections = elems
        pobj.num_conn = len(p.face_list)

        surf_reg_faces = None
        for idx in p.surf_reg_face_list:
            surf_reg_faces = m.mcell_add_to_region_list(surf_reg_faces, idx)

        pobj_list = m.mcell_add_to_poly_obj_list(
                pobj_list, p.obj_name, verts, len(p.vert_list), elems,
                len(p.face_list), surf_reg_faces, p.reg_name)

    m.mcell_change_geometry(world, pobj_list)


def create_polygon_object(world, vert_list, face_list, scene, name):
    """Creates a polygon object from a vertex and element lest

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        vert_list (vertex list) -- verteces for the polygon
        face_list (element connection list) -- faces for the polygon
        scene (instance object) -- scene for mcell simulation
        name (string) -- name of polygon object that will be created

    Returns:
        polygon object
    """

    verts = None
    vert_list = vert_list[::-1]
    for x, y, z in vert_list:
        verts = m.mcell_add_to_vertex_list(x, y, z, verts)

    elems = None
    face_list = face_list[::-1]
    for x, y, z in face_list:
        elems = m.mcell_add_to_connection_list(x, y, z, elems)

    pobj = m.poly_object()
    pobj.obj_name = name
    pobj.vertices = verts
    pobj.num_vert = len(vert_list)
    pobj.connections = elems
    pobj.num_conn = len(face_list)

    mesh_temp = m.object()
    mesh = m.mcell_create_poly_object(world, scene, pobj, mesh_temp)

    return mesh


def create_surface_region(world, mesh, surf_reg_face_list, region_name):
    """Creates a surface region

    Args:
        world (object) -- the world object which has been generated by
            mcell create_instance_object
        mesh (polygon object) -- object where surface region will reside
        surf_reg_face_list (element connection list) -- list of surface
            region faces
        region_name (string) -- name of the surface region being created

    Returns:
        region object

    """

    surface_region = m.mcell_create_region(world, mesh, region_name)

    surf_reg_faces = None
    for idx in surf_reg_face_list:
        surf_reg_faces = m.mcell_add_to_region_list(surf_reg_faces, idx)

    m.mcell_set_region_elements(surface_region, surf_reg_faces, 1)

    return surface_region
