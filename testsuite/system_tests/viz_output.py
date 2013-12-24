###################################################################################
#                                                                                 #
# Copyright (C) 2006-2013 by                                                      #
# The Salk Institute for Biological Studies and                                   #
# Pittsburgh Supercomputing Center, Carnegie Mellon University                    #
#                                                                                 #
# This program is free software; you can redistribute it and/or                   #
# modify it under the terms of the GNU General Public License                     #
# as published by the Free Software Foundation; either version 2                  #
# of the License, or (at your option) any later version.                          #
#                                                                                 #
# This program is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                  #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   #
# GNU General Public License for more details.                                    #
#                                                                                 #
# You should have received a copy of the GNU General Public License               #
# along with this program; if not, write to the Free Software                     #
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. #
#                                                                                 #
###################################################################################

import os
from testutils import assertFileExists
from testutils import assertFileNotExists
from testutils import assertFileNonempty
from testutils import assertFileSymlink
from testutils import assertFileDir
from testutils import safe_concat


def assertValidVizFileAscii(fname, sstates=None, vstates=None):
    """
    Give an error if the specified ASCII viz file is invalid, according to the
    specified parameters.

    Parameters:
        fname        - the filename to check
        sstates      - allowed surface states, or None to skip id validation
        vstates      - allowed volume states, or None to skip id validation

    Example usage:
        assertValidVizFileAscii("viz_dat/molecules.0.dat", [1, 2, 3], [3, 4, 5])

             This will complain if the ASCII-format viz output file is
             malformed, or if any surface molecules are emitted with states
             other than 1, 2, or 3 (or, alternatively, if any volume molecules
             are emitted with a normal vector). Likewise, it will complain if
             any volume molecules are emitted with states other than 3, 4, or 5,
             or if any surface molecules are emitted with normals of 0, 0, 0.

        assertValidVizFileAscii("viz_dat/molecules.0.dat")

             The above invocation will check only for validity, ignoring the
             molecule state values.

    """

    try:
        got_contents = open(fname).read()
    except:
        assert False, "Expected ASCII viz output file '%s' was not created" % fname

    # Rend file into tiny pieces
    lines = [l for l in got_contents.split('\n') if l != '']

    # Check each line
    idx = 1
    for idx in range(len(lines)):
        l = lines[idx]
        v = l.split()

        # Check number of columns
        assert len(v) == 8, "In ASCII viz output file '%s', line %d is malformed (%d columns instead of 8)" % (fname, idx + 1, len(v))

        # Check format of columns
        try:
            float(v[1])
            float(v[2])
            float(v[3])
        except:
            assert False, "In ASCII viz output file '%s', line %d, position is malformed (was %s %s %s, should be a triple of floating point values)" % (fname, idx + 1, v[1], v[2], v[3])
        try:
            float(v[4])
            float(v[5])
            float(v[6])
        except:
            assert False, "In ASCII viz output file '%s', line %d, normal vector is malformed (was %s %s %s, should be a triple of floating point values)" % (fname, idx + 1, v[4], v[5], v[6])

        # Maybe check states
        if sstates is not None or vstates is not None:
            id = v[0]
            if sstates is not None:
                if type(id) == int:
                    assert id in sstates, "In ASCII viz output file '%s', line %d looks like a surface molecule, but has the state of a volume molecule (state=%d)" % (fname, idx + 1, id)
            if vstates is not None:
                if type(id) == int:
                    assert id in vstates, "In ASCII viz output file '%s', line %d looks like a volume molecule, but has the state of a surface molecule (state=%d)" % (fname, idx + 1, id)


class RequireVizAscii:
    def __init__(self, basename, iters, sstates=None, vstates=None, astates=None):
        self.basename = basename
        self.iters = iters
        self.args = {}
        if sstates is not None:
            self.args["sstates"] = safe_concat(sstates, astates)
        if vstates is not None:
            self.args["vstates"] = safe_concat(vstates, astates)

    def check(self):
        for i in self.iters:
            assertValidVizFileAscii(self.basename + ".ascii.%03d" % i + ".dat", **self.args)


def assertValidVizFilesDx(dir, molfile=None, objprefixes=None, alliters=None, mpositers=None, mstateiters=None, epositers=None, estateiters=None, opositers=None, ostateiters=None):
    """
    Give an error if the specified DX mode viz file is invalid, according to the
    specified parameters.

    Note: Currently, this test does NOT validate that unexpected outputs were
    NOT produced -- only that the expected outputs WERE produced.    Ideally, it
    should check that we only got the output types we expected, and only on the
    iterations where we expected them.

    Parameters:
        dir         - the directory containing the output files
        molfile     - MOLECULE_FILE_PREFIX from viz output block
        objprefixes - OBJECT_FILE_PREFIXES from viz output block
        alliters    - iteration numbers where all types of output were produced
        mpositers   - iteration numbers where vol mol position was produced
        mstateiters - iteration numbers where vol mol state was produced
        epositers   - iteration numbers where grid mol position was produced
        estateiters - iteration numbers where grid mol state was produced
        opositers   - iteration numbers where mesh position was produced
        ostateiters - iteration numbers where mesh state was produced

    """

    if mpositers is None:
        mpositers = []
    if mstateiters is None:
        mstateiters = []
    if epositers is None:
        epositers = []
    if estateiters is None:
        estateiters = []
    if opositers is None:
        opositers = []
    if ostateiters is None:
        ostateiters = []
    if alliters is not None:
        mpositers.extend(alliters)
        mstateiters.extend(alliters)
        epositers.extend(alliters)
        estateiters.extend(alliters)
        opositers.extend(alliters)
        ostateiters.extend(alliters)

    assertFileDir(dir)
    if molfile is not None:
        for i in mpositers:
            assertFileNonempty(os.path.join(dir, "%s.molecule_positions.%d.dx" % (molfile, i)))
        for i in mstateiters:
            assertFileNonempty(os.path.join(dir, "%s.molecule_states.%d.dx" % (molfile, i)))
    if objprefixes is not None:
        for b in objprefixes:
            for i in epositers:
                assertFileNonempty(os.path.join(dir, "%s.effector_site_positions.%d.dx" % (b, i)))
            for i in estateiters:
                assertFileNonempty(os.path.join(dir, "%s.effector_site_states.%d.dx" % (b, i)))
            for i in opositers:
                assertFileNonempty(os.path.join(dir, "%s.mesh_elements.%d.dx" % (b, i)))
            for i in ostateiters:
                assertFileNonempty(os.path.join(dir, "%s.mesh_element_states.%d.dx" % (b, i)))


class RequireVizDX:
    def __init__(self, dir, molfile=None, objprefixes=None, alliters=None, mpositers=None, mstateiters=None, epositers=None, estateiters=None, opositers=None, ostateiters=None):
        self.dir = dir
        self.molfile = molfile
        self.objprefixes = objprefixes
        self.mpositers = mpositers
        self.mstateiters = mstateiters
        self.epositers = epositers
        self.estateiters = estateiters
        self.opositers = opositers
        self.ostateiters = ostateiters
        self.alliters = alliters

    def check(self):
        assertValidVizFilesDx(self.dir,
                              molfile=self.molfile,
                              objprefixes=self.objprefixes,
                              alliters=self.alliters,
                              mpositers=self.mpositers, mstateiters=self.mstateiters,
                              epositers=self.epositers, estateiters=self.estateiters,
                              opositers=self.opositers, ostateiters=self.ostateiters)


def assertValidVizFilesDreammV3(dir, name, n_iters=None, n_times=None):
    """
    Give an error if the specified DREAMM V3 (non-grouped) mode viz output set
    is invalid, according to the specified parameters. This only checks the
    top-level files, not the frame data. The frame data will be checked by
    separate assertions which vary for binary and ascii output formats.

    Note: Currently, this test does NOT validate that unexpected outputs were
    NOT produced -- only that the expected outputs WERE produced. Ideally, it
    should check that we only got the output types we expected, and only on the
    iterations where we expected them.

    Parameters:
        dir     - the directory containing the output files
        name    - name of the output file set
        n_iters - number of output iterations, or None to not check
        n_times - number of distinct output times, or None to not check

    """

    assertFileNonempty(os.path.join(dir, name + ".dx"))
    if n_iters is not None:
        assertFileNonempty(os.path.join(dir, name + ".iteration_numbers.bin"), 12 * n_iters)
    else:
        assertFileNonempty(os.path.join(dir, name + ".iteration_numbers.bin"))
    if n_times is not None:
        assertFileNonempty(os.path.join(dir, name + ".time_values.bin"), 8 * n_times)
    else:
        assertFileNonempty(os.path.join(dir, name + ".time_values.bin"))


class RequireVizDreammV3:
    def __init__(self, dir, name, n_iters=None, n_times=None):
        self.dir = dir
        self.name = name
        self.args = {}
        if n_iters is not None:
            self.args["n_iters"] = n_iters
        if n_times is not None:
            self.args["n_times"] = n_times

    def check(self):
        assertValidVizFilesDreammV3(self.dir, self.name, **self.args)


def assertValidVizFilesDreammV3MolsBin(dir, alliters,
                                       surfpositers=None,
                                       surforientiters=None,
                                       surfstateiters=[],
                                       surfnonempty=True,
                                       volpositers=None,
                                       volorientiters=None,
                                       volstateiters=[],
                                       volnonempty=True):
    """
    Give an error if the specified DREAMM V3 (non-grouped) mode viz output set
    contains invalid (or fails to contain valid) binary molecule data files.

    Parameters:
        dir             - the directory containing the output files
        alliters        - sorted iteration numbers where any output was produced
        surfpositers    - iters producing surface mol positions
        surforientiters - iters producing surface mol orientations
        surfstateiters  - iters producing surface mol states
        surfnonempty    - expect at least one surface molecule
        volpositers     - iters producing volume mol positions
        volorientiters  - iters producing volume mol orientations
        volstateiters   - iters producing volume mol states
        volnonempty     - expect at least one volume molecule

    """

    last_spos = [None]
    last_sorients = [None]
    last_sstate = [None]
    last_vpos = [None]
    last_vorients = [None]
    last_vstate = [None]

    # Reset last iterations if our output types are not totally synchronized
    def check_unset(s):
        if last_spos[0] != s:
            last_spos[0] = None
        if last_sorients[0] != s:
            last_sorients[0] = None
        if last_sstate[0] != s:
            last_sstate[0] = None
        if last_vpos[0] != s:
            last_vpos[0] = None
        if last_vorients[0] != s:
            last_vorients[0] = None
        if last_vstate[0] != s:
            last_vstate[0] = None

    # Make sure iteration list is sorted, create sets for quick searching
    alliters.sort()
    if surfpositers is None:
        surfpositers = alliters
    surfpositers = set(surfpositers)
    if surforientiters is None:
        surforientiters = alliters
    surforientiters = set(surforientiters)
    if surfstateiters is None:
        surfstateiters = alliters
    surfstateiters = set(surfstateiters)
    if volpositers is None:
        volpositers = alliters
    volpositers = set(volpositers)
    if volorientiters is None:
        volorientiters = alliters
    volorientiters = set(volorientiters)
    if volstateiters is None:
        volstateiters = alliters
    volstateiters = set(volstateiters)
    surfiters = surfpositers.union(surforientiters).union(surfstateiters)
    voliters = volpositers.union(volorientiters).union(volstateiters)
    moliters = surfiters.union(voliters)

    # Check each iteration
    for iter in alliters:
        basedir = os.path.join(dir, "frame_data")
        basedir = os.path.join(basedir, "iteration_%d" % iter)
        had_sframe = 0
        had_vframe = 0

        surf_pos = os.path.join(basedir, "surface_molecules_positions.bin")
        surf_orient = os.path.join(basedir, "surface_molecules_orientations.bin")
        surf_states = os.path.join(basedir, "surface_molecules_states.bin")
        surf_header = os.path.join(basedir, "surface_molecules.dx")
        vol_pos = os.path.join(basedir, "volume_molecules_positions.bin")
        vol_orient = os.path.join(basedir, "volume_molecules_orientations.bin")
        vol_states = os.path.join(basedir, "volume_molecules_states.bin")
        vol_header = os.path.join(basedir, "volume_molecules.dx")

        if iter in surfpositers:
            had_sframe = 1
            last_spos[0] = iter
            check_unset(iter)

            if surfnonempty:
                assertFileNonempty(surf_pos)
            else:
                assertFileExists(surf_pos)
        elif last_spos[0] is not None and iter not in moliters:
            assertFileSymlink(surf_pos, "../iteration_%d/surface_molecules_positions.bin" % last_spos[0])
        else:
            assertFileNotExists(surf_pos)

        if iter in surforientiters:
            had_sframe = 1
            last_sorients[0] = iter
            check_unset(iter)

            if surfnonempty:
                assertFileNonempty(surf_orient)
            else:
                assertFileExists(surf_orient)
        elif last_sorients[0] is not None and iter not in moliters:
            assertFileSymlink(surf_orient, "../iteration_%d/surface_molecules_orientations.bin" % last_sorients[0])
        else:
            assertFileNotExists(surf_orient)

        if iter in surfstateiters:
            had_sframe = 1
            last_sstate[0] = iter
            check_unset(iter)

            if surfnonempty:
                assertFileNonempty(surf_states)
            else:
                assertFileExists(surf_states)
        elif last_sstate[0] is not None and iter not in moliters:
            assertFileSymlink(surf_states, "../iteration_%d/surface_molecules_states.bin" % last_sstate[0])
        else:
            assertFileNotExists(surf_states)

        if had_sframe:
            assertFileNonempty(surf_header)
        elif last_spos[0] is not None:
            assertFileSymlink(surf_header, "../iteration_%d/surface_molecules.dx" % last_spos[0])
        elif last_sorients[0] is not None:
            assertFileSymlink(surf_header, "../iteration_%d/surface_molecules.dx" % last_sorients[0])
        elif last_sstate[0] is not None:
            assertFileSymlink(surf_header, "../iteration_%d/surface_molecules.dx" % last_sstate[0])
        else:
            assertFileNotExists(surf_header)

        if iter in volpositers:
            had_vframe = 1
            last_vpos[0] = iter
            check_unset(iter)

            if volnonempty:
                assertFileNonempty(vol_pos)
            else:
                assertFileExists(vol_pos)
        elif last_vpos[0] is not None and iter not in moliters:
            assertFileSymlink(vol_pos, "../iteration_%d/volume_molecules_positions.bin" % last_vpos[0])
        else:
            assertFileNotExists(vol_pos)

        if iter in volorientiters:
            had_vframe = 1
            last_vorients[0] = iter
            check_unset(iter)

            if volnonempty:
                assertFileNonempty(vol_orient)
            else:
                assertFileExists(vol_orient)
        elif last_vorients[0] is not None and iter not in moliters:
            assertFileSymlink(vol_orient, "../iteration_%d/volume_molecules_orientations.bin" % last_vorients[0])
        else:
            assertFileNotExists(vol_orient)

        if iter in volstateiters:
            had_vframe = 1
            last_vstate[0] = iter
            check_unset(iter)

            if volnonempty:
                assertFileNonempty(vol_states)
            else:
                assertFileExists(vol_states)
        elif last_vstate[0] is not None and iter not in moliters:
            assertFileSymlink(vol_states, "../iteration_%d/volume_molecules_states.bin" % last_vstate[0])
        else:
            assertFileNotExists(vol_states)

        if had_vframe:
            assertFileNonempty(vol_header)
        elif last_vpos[0] is not None:
            assertFileSymlink(vol_header, "../iteration_%d/volume_molecules.dx" % last_vpos[0])
        elif last_vorients[0] is not None:
            assertFileSymlink(vol_header, "../iteration_%d/volume_molecules.dx" % last_vorients[0])
        elif last_vstate[0] is not None:
            assertFileSymlink(vol_header, "../iteration_%d/volume_molecules.dx" % last_vstate[0])
        else:
            assertFileNotExists(vol_header)


class RequireVizDreammV3MolsBin:
    def __init__(self, dir, alliters,
                 surfpositers=None,
                 surforientiters=None,
                 surfstateiters=None,
                 surfnonempty=None,
                 volpositers=None,
                 volorientiters=None,
                 volstateiterss=None,
                 volnonempty=None):
        self.dir = dir
        self.alliters = alliters
        self.args = {}
        if surfpositers is not None:
            self.args["surfpositers"] = surfpositers
        if surforientiters is not None:
            self.args["surforientiters"] = surforientiters
        if surfstateiters is not None:
            self.args["surfstateiters"] = surfstateiters
        if surfnonempty is not None:
            self.args["surfnonempty"] = surfnonempty
        if volpositers is not None:
            self.args["volpositers"] = volpositers
        if volorientiters is not None:
            self.args["volorientiters"] = volorientiters
        if volstateiterss is not None:
            self.args["volstateiterss"] = volstateiterss
        if volnonempty is not None:
            self.args["volnonempty"] = volnonempty

    def check(self):
        assertValidVizFilesDreammV3MolsBin(self.dir, self.alliters, **self.args)


def assertValidVizFilesDreammV3MolsAscii(dir, alliters,
                                         molnames,
                                         positers=None,
                                         orientiters=None,
                                         stateiters=[]):

    """
    Give an error if the specified DREAMM V3 (non-grouped) mode viz output set
    contains invalid (or fails to contain valid) ascii molecule data files.

    Parameters:
        dir         - the directory containing the output files
        alliters    - sorted iteration numbers where any output was produced
        molnames    - list of molecule names to be written
        positers    - iters producing surface mol positions
        orientiters - iters producing surface mol orientations
        stateiters  - iters producing surface mol states

    """

    last_pos = [None]
    last_orients = [None]
    last_state = [None]

    # Reset last iterations if our output types are not totally synchronized
    def check_unset(s):
        if last_pos[0] != s:
            last_pos[0] = None
        if last_orients[0] != s:
            last_orients[0] = None
        if last_state[0] != s:
            last_state[0] = None

    # Make sure iteration list is sorted, create sets for quick searching
    alliters.sort()
    if positers is None:
        positers = alliters
    positers = set(positers)
    if orientiters is None:
        orientiters = alliters
    orientiters = set(orientiters)
    if stateiters is None:
        stateiters = alliters
    stateiters = set(stateiters)
    moliters = positers.union(orientiters).union(stateiters)

    # Check each iteration
    for iter in alliters:
        basedir = os.path.join(dir, "frame_data")
        basedir = os.path.join(basedir, "iteration_%d" % iter)
        had_sframe = 0
        #had_vframe = 0

        surf_header = os.path.join(basedir, "surface_molecules.dx")
        vol_header = os.path.join(basedir, "volume_molecules.dx")

        if iter in positers:
            had_sframe = 1
            last_pos[0] = iter
            check_unset(iter)

            for mol in molnames:
                pos = os.path.join(basedir, "%s.positions.dat" % mol)
                assertFileExists(pos)
        elif last_pos[0] is not None and iter not in moliters:
            for mol in molnames:
                pos = os.path.join(basedir, "%s.positions.dat" % mol)
                assertFileSymlink(pos, "../iteration_%d/%s.positions.dat" % (last_pos[0], mol))
        else:
            for mol in molnames:
                pos = os.path.join(basedir, "%s.positions.dat" % mol)
                assertFileNotExists(pos)

        #f_orient = os.path.join(basedir, "orientations.dat")
        #f_states = os.path.join(basedir, "states.dat")

        if iter in orientiters:
            had_sframe = 1
            last_orients[0] = iter
            check_unset(iter)

            for mol in molnames:
                orient = os.path.join(basedir, "%s.orientations.dat" % mol)
                assertFileExists(orient)
        elif last_orients[0] is not None and iter not in moliters:
            for mol in molnames:
                orient = os.path.join(basedir, "%s.orientations.dat" % mol)
                assertFileSymlink(orient, "../iteration_%d/%s.orientations.dat" % (last_orients[0], mol))
        else:
            for mol in molnames:
                orient = os.path.join(basedir, "%s.orientations.dat" % mol)
                assertFileNotExists(orient)

        if iter in stateiters:
            had_sframe = 1
            last_state[0] = iter
            check_unset(iter)

            for mol in molnames:
                states = os.path.join(basedir, "%s.states.dat" % mol)
                assertFileNonempty(states)
        elif last_state[0] is not None and iter not in moliters:
            for mol in molnames:
                states = os.path.join(basedir, "%s.states.dat" % mol)
                assertFileSymlink(states, "../iteration_%d/%s.states.dat" % (last_state[0], mol))
        else:
            for mol in molnames:
                states = os.path.join(basedir, "%s.states.dat" % mol)
                assertFileNotExists(states)

        if had_sframe:
            assertFileExists(surf_header)
            assertFileExists(vol_header)
        elif last_pos[0] is not None:
            assertFileSymlink(surf_header, "../iteration_%d/surface_molecules.dx" % last_pos[0])
            assertFileSymlink(vol_header, "../iteration_%d/volume_molecules.dx" % last_pos[0])
        elif last_orients[0] is not None:
            assertFileSymlink(surf_header, "../iteration_%d/surface_molecules.dx" % last_orients[0])
            assertFileSymlink(vol_header, "../iteration_%d/volume_molecules.dx" % last_orients[0])
        elif last_state[0] is not None:
            assertFileSymlink(surf_header, "../iteration_%d/surface_molecules.dx" % last_state[0])
            assertFileSymlink(vol_header, "../iteration_%d/volume_molecules.dx" % last_state[0])
        else:
            assertFileNotExists(surf_header)
            assertFileNotExists(vol_header)


class RequireVizDreammV3MolsAscii:
    def __init__(self, dir, alliters,
                 molnames,
                 positers=None,
                 orientiters=None,
                 stateiters=None):
        self.dir = dir
        self.alliters = alliters
        self.molnames = molnames
        self.args = {}
        if positers is not None:
            self.args["positers"] = positers
        if orientiters is not None:
            self.args["orientiters"] = orientiters
        if stateiters is not None:
            self.args["stateiters"] = stateiters

    def check(self):
        assertValidVizFilesDreammV3MolsAscii(self.dir, self.alliters, self.molnames, **self.args)


def assertValidVizFilesDreammV3MeshBin(dir, alliters, positers=None, regioniters=None, stateiters=[], meshnonempty=True):
    """
    Give an error if the specified DREAMM V3 (non-grouped) mode viz output set
    contains invalid (or fails to contain valid) binary mesh data files.

    Parameters:
        dir          - the directory containing the output files
        alliters     - iteration numbers where any output was produced
        positers     - iteration numbers where mesh pos output was produced, if different from alliters
        regioniters  - iteration numbers where region data output was produced, if different from alliters
        stateiters   - iteration numbers where mesh state was produced, if any (it is assumed that no mesh state data was produced, otherwise)
        meshnonempty - expect at least one mesh

    """

    last_pos = [None]
    last_rgn = [None]
    last_state = [None]

    # Reset last iterations if our output types are not totally synchronized
    def check_unset(s):
        if last_pos[0] != s:
            last_pos[0] = None
        if last_rgn[0] != s:
            last_rgn[0] = None
        if last_state[0] != s:
            last_state[0] = None

    alliters.sort()
    if positers is None:
        positers = alliters
    positers = set(positers)
    if regioniters is None:
        regioniters = alliters
    regioniters = set(regioniters)
    if stateiters is None:
        stateiters = alliters
    stateiters = set(stateiters)

    for iter in alliters:
        basedir = os.path.join(dir, "frame_data")
        basedir = os.path.join(basedir, "iteration_%d" % iter)
        had_frame = 0

        mesh_pos = os.path.join(basedir, "mesh_positions.bin")
        mesh_rgn = os.path.join(basedir, "region_indices.bin")
        mesh_state = os.path.join(basedir, "mesh_states.bin")
        mesh_header = os.path.join(basedir, "meshes.dx")

        # Check for positions file or symlink
        if iter in positers:
            had_frame = 1
            last_pos[0] = iter
            check_unset(iter)

            if meshnonempty:
                assertFileNonempty(mesh_pos)
            else:
                assertFileExists(mesh_pos)
        elif last_pos[0] is not None and iter not in regioniters and iter not in stateiters:
            assertFileSymlink(mesh_pos, "../iteration_%d/mesh_positions.bin" % last_pos[0])
        else:
            assertFileNotExists(mesh_pos)

        # Check for regions file or symlink
        if iter in regioniters:
            had_frame = 1
            last_rgn[0] = iter
            check_unset(iter)

            if meshnonempty:
                assertFileNonempty(mesh_rgn)
            else:
                assertFileExists(mesh_rgn)
        elif last_rgn[0] is not None and iter not in stateiters:
            assertFileSymlink(mesh_rgn, "../iteration_%d/region_indices.bin" % last_rgn[0])
        else:
            assertFileNotExists(mesh_rgn)

        # Check for states file or symlink
        if iter in stateiters:
            had_frame = 1
            last_state[0] = iter
            check_unset(iter)

            if meshnonempty:
                assertFileNonempty(mesh_state)
            else:
                assertFileExists(mesh_state)
        elif last_state[0] is not None:
            assertFileSymlink(mesh_state, "../iteration_%d/mesh_states.bin" % last_state[0])
        else:
            assertFileNotExists(mesh_state)

        # Check for meshes file or symlink
        if had_frame:
            assertFileNonempty(mesh_header)
        elif last_pos[0] is not None:
            assertFileSymlink(mesh_header, "../iteration_%d/meshes.dx" % last_pos[0])
        elif last_rgn[0] is not None:
            assertFileSymlink(mesh_header, "../iteration_%d/meshes.dx" % last_rgn[0])
        elif last_state[0] is not None:
            assertFileSymlink(mesh_header, "../iteration_%d/meshes.dx" % last_state[0])
        else:
            assertFileNotExists(mesh_header)


class RequireVizDreammV3MeshBin:
    def __init__(self, dir, alliters, positers=None, regioniters=None, stateiters=None, meshnonempty=None):
        self.dir = dir
        self.alliters = alliters
        self.args = {}
        if positers is not None:
            self.args["positers"] = positers
        if regioniters is not None:
            self.args["regioniters"] = regioniters
        if stateiters is not None:
            self.args["stateiters"] = stateiters
        if meshnonempty is not None:
            self.args["meshnonempty"] = meshnonempty

    def check(self):
        assertValidVizFilesDreammV3MeshBin(self.dir, self.alliters, **self.args)


def assertValidVizFilesDreammV3MeshAscii(dir, alliters, objnames, objswithregions=None, positers=None, regioniters=None, stateiters=[], meshnonempty=True):
    """
    Give an error if the specified DREAMM V3 (non-grouped) mode viz output set
    contains invalid (or fails to contain valid) ascii mesh data files.

    Parameters:
        dir             - the directory containing the output files
        alliters        - iteration numbers where any output was produced
        objnames        - names of objects whose files to look for
        objswithregions - names of objects which have regions other than ALL
        positers        - iteration numbers where mesh pos output was produced, if different from alliters
        regioniters     - iteration numbers where region data output was produced, if different from alliters
        stateiters      - iteration numbers where mesh state was produced, if any (it is assumed that no mesh state data was produced, otherwise)
        meshnonempty    - expect at least one mesh

    """

    last_pos = [None]
    last_rgn = [None]
    last_state = [None]

    # Reset last iterations if our output types are not totally synchronized
    def check_unset(s):
        if last_pos[0] != s:
            last_pos[0] = None
        if last_rgn[0] != s:
            last_rgn[0] = None
        if last_state[0] != s:
            last_state[0] = None

    alliters.sort()
    if positers is None:
        positers = alliters
    positers = set(positers)
    if regioniters is None:
        regioniters = alliters
    regioniters = set(regioniters)
    if stateiters is None:
        stateiters = alliters
    stateiters = set(stateiters)

    if objswithregions is None:
        objswithregions = objnames

    for iter in alliters:
        basedir = os.path.join(dir, "frame_data")
        basedir = os.path.join(basedir, "iteration_%d" % iter)
        had_frame = 0

        mesh_header = os.path.join(basedir, "meshes.dx")

        # Check for positions file or symlink
        if iter in positers:
            had_frame = 1
            last_pos[0] = iter
            check_unset(iter)

            for obj in objnames:
                mesh_pos = os.path.join(basedir, obj + ".positions.dat")
                mesh_conn = os.path.join(basedir, obj + ".connections.dat")
                assertFileNonempty(mesh_pos)
                assertFileNonempty(mesh_conn)
        elif last_pos[0] is not None and iter not in regioniters and iter not in stateiters:
            for obj in objnames:
                mesh_pos = os.path.join(basedir, obj + ".positions.dat")
                mesh_conn = os.path.join(basedir, obj + ".connections.dat")
                assertFileSymlink(mesh_pos, "../iteration_%d/%s.positions.dat" % (last_pos[0], obj))
                assertFileSymlink(mesh_conn, "../iteration_%d/%s.connections.dat" % (last_pos[0], obj))
        else:
            for obj in objnames:
                mesh_pos = os.path.join(basedir, obj + ".positions.dat")
                mesh_conn = os.path.join(basedir, obj + ".connections.dat")
                assertFileNotExists(mesh_pos)
                assertFileNotExists(mesh_conn)

        # Check for regions file or symlink
        if iter in regioniters:
            had_frame = 1
            last_rgn[0] = iter
            check_unset(iter)

            for obj in objswithregions:
                mesh_rgn = os.path.join(basedir, obj + ".region_indices.dat")
                assertFileExists(mesh_rgn)
        elif last_rgn[0] is not None and iter not in stateiters:
            for obj in objswithregions:
                mesh_rgn = os.path.join(basedir, obj + ".region_indices.dat")
                assertFileSymlink(mesh_rgn, "../iteration_%d/%s.region_indices.dat" % (last_rgn[0], obj))
        else:
            for obj in objswithregions:
                mesh_rgn = os.path.join(basedir, obj + ".region_indices.dat")
                assertFileNotExists(mesh_rgn)

        # Check for states file or symlink
        if iter in stateiters:
            had_frame = 1
            last_state[0] = iter
            check_unset(iter)

            for obj in objnames:
                mesh_state = os.path.join(basedir, obj + ".states.dat")
                assertFileNonempty(mesh_state)
        elif last_state[0] is not None:
            for obj in objnames:
                mesh_state = os.path.join(basedir, obj + ".states.dat")
                assertFileSymlink(mesh_state, "../iteration_%d/%s.states.dat" % (last_state[0], obj))
        else:
            for obj in objnames:
                mesh_state = os.path.join(basedir, obj + ".states.dat")
                assertFileNotExists(mesh_state)

        # Check for meshes file or symlink
        if had_frame:
            assertFileNonempty(mesh_header)
        elif last_pos[0] is not None:
            assertFileSymlink(mesh_header, "../iteration_%d/meshes.dx" % last_pos[0])
        elif last_rgn[0] is not None:
            assertFileSymlink(mesh_header, "../iteration_%d/meshes.dx" % last_rgn[0])
        elif last_state[0] is not None:
            assertFileSymlink(mesh_header, "../iteration_%d/meshes.dx" % last_state[0])
        else:
            assertFileNotExists(mesh_header)


class RequireVizDreammV3MeshAscii:
    def __init__(self, dir, alliters, objnames, objswithregions=None, positers=None, regioniters=None, stateiters=None, meshnonempty=None):
        self.dir = dir
        self.alliters = alliters
        self.objnames = objnames
        self.args = {}
        if objswithregions is not None:
            self.args["objswithregions"] = objswithregions
        if positers is not None:
            self.args["positers"] = positers
        if regioniters is not None:
            self.args["regioniters"] = regioniters
        if stateiters is not None:
            self.args["stateiters"] = stateiters
        if meshnonempty is not None:
            self.args["meshnonempty"] = meshnonempty

    def check(self):
        assertValidVizFilesDreammV3MeshAscii(self.dir, self.alliters, self.objnames, **self.args)


def assertValidVizFilesDreammV3Grouped(dir, name,
                                       cpno=1,
                                       n_iters=None,
                                       n_times=None,
                                       meshpos=True,
                                       rgnindx=True,
                                       meshstate=False,
                                       meshnonempty=True,
                                       molpos=True,
                                       molorient=True,
                                       molstate=False,
                                       molsnonempty=True):
    """
    Give an error if the specified DREAMM V3 (grouped) mode viz output set is
    invalid.

    Parameters:
        dir          - the directory containing the output files
        name         - name of the output set
        cpno         - checkpoint sequence number (assumed to be 1 otherwise)
        n_iters      - total number of iterations with output (None to skip check)
        n_times      - total number of distinct time points with output (None to skip check)
        meshpos      - true to expect mesh pos file, false to not
        rgnindx      - true to expect region indices file, false to not
        meshstate    - true to expect mesh state file, false to not
        meshnonempty - true to expect at least 1 mesh
        molpos       - true to expect mol pos file, false to not
        molorient    - true to expect mol orient file, false to not
        molstate     - true to expect mol state file, false to not
        molsnonempty - true to expect at least 1 molecule

    """

    fmtargs = name, cpno
    assertFileNonempty(os.path.join(dir, "%s.%d.dx" % fmtargs))

    path = os.path.join(dir, "%s.mesh_positions.%d.bin" % fmtargs)
    if meshpos:
        if meshnonempty:
            assertFileNonempty(path)
        else:
            assertFileExists(path)
    else:
        assertFileNotExists(path)

    path = os.path.join(dir, "%s.region_indices.%d.bin" % fmtargs)
    if rgnindx:
        assertFileExists(path)
    else:
        assertFileNotExists(path)

    path = os.path.join(dir, "%s.mesh_states.%d.bin" % fmtargs)
    if meshstate:
        if meshnonempty:
            assertFileNonempty(path)
        else:
            assertFileExists(path)
    else:
        assertFileNotExists(path)

    path = os.path.join(dir, "%s.molecule_positions.%d.bin" % fmtargs)
    if molpos:
        if molsnonempty:
            assertFileNonempty(path)
        else:
            assertFileExists(path)
    else:
        assertFileNotExists(path)

    path = os.path.join(dir, "%s.molecule_orientations.%d.bin" % fmtargs)
    if molorient:
        if molsnonempty:
            assertFileNonempty(path)
        else:
            assertFileExists(path)
    else:
        assertFileNotExists(path)

    path = os.path.join(dir, "%s.molecule_states.%d.bin" % fmtargs)
    if molstate:
        if molsnonempty:
            assertFileNonempty(path)
        else:
            assertFileExists(path)
    else:
        assertFileNotExists(path)

    ipath = os.path.join(dir, "%s.iteration_numbers.%d.bin" % fmtargs)
    tpath = os.path.join(dir, "%s.time_values.%d.bin" % fmtargs)
    if n_iters is not None:
        assertFileNonempty(ipath, 12 * n_iters)
    else:
        assertFileNonempty(ipath)
    if n_times is not None:
        assertFileNonempty(tpath, 8 * n_times)
    else:
        assertFileNonempty(tpath)


class RequireVizDreammV3Grouped:
    def __init__(self, dir, name,
                 cpno=None,
                 n_iters=None,
                 n_times=None,
                 meshpos=None,
                 rgnindx=None,
                 meshstate=None,
                 meshnonempty=None,
                 molpos=None,
                 molorient=None,
                 molstate=None,
                 molsnonempty=None):
        self.dir = dir
        self.name = name
        self.args = {}
        if cpno is not None:
            self.args["cpno"] = cpno
        if n_iters is not None:
            self.args["n_iters"] = n_iters
        if n_times is not None:
            self.args["n_times"] = n_times
        if meshpos is not None:
            self.args["meshpos"] = meshpos
        if rgnindx is not None:
            self.args["rgnindx"] = rgnindx
        if meshstate is not None:
            self.args["meshstate"] = meshstate
        if meshnonempty is not None:
            self.args["meshnonempty"] = meshnonempty
        if molpos is not None:
            self.args["molpos"] = molpos
        if molorient is not None:
            self.args["molorient"] = molorient
        if molstate is not None:
            self.args["molstate"] = molstate
        if molsnonempty is not None:
            self.args["molsnonempty"] = molsnonempty

    def check(self):
        assertValidVizFilesDreammV3Grouped(self.dir, self.name, **self.args)
        # invalid.
        #
        # Parameters:
        #    dir          - the directory containing the output files
        #    name         - name of the output set
        #    cpno         - checkpoint sequence number (assumed to be 1 otherwise)
        #    n_iters      - total number of iterations with output (None to skip check)
        #    n_times      - total number of distinct time points with output (None to skip check)
        #    meshpos      - true to expect mesh pos file, false to not
        #    rgnindx      - true to expect region indices file, false to not
        #    meshstate    - true to expect mesh state file, false to not
        #    meshnonempty - true to expect at least 1 mesh
        #    molpos       - true to expect mol pos file, false to not
        #    molorient    - true to expect mol orient file, false to not
        #    molstate     - true to expect mol state file, false to not
        #    molsnonempty - true to expect at least 1 molecule
        ###################################################################
        assertValidVizFilesDreammV3Grouped(self.dir, self.name, **self.args)
