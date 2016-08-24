/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

MCELL_STATUS mcell_add_surf_class_properties(
    MCELL_STATE *state,
    int reaction_type,
    mcell_symbol *sc_sym,
    mcell_symbol *reactant_sym,
    short orient);

%typemap(in) mcell_symbol **sc_sym (mcell_symbol *temp) {
  $1 = &temp;
}

%typemap(argout) struct sym_entry **sc_sym {
  %set_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

MCELL_STATUS mcell_create_surf_class(
    MCELL_STATE *state,
    char *surf_class_name,
    mcell_symbol **sc_sym);

struct sm_dat *mcell_add_mol_release_to_surf_class(
    MCELL_STATE *state,
    struct sym_entry *sc_sym,
    struct mcell_species *sm_info,
    double quantity,
    int density_or_num,
    struct sm_dat *smd_list);

MCELL_STATUS mcell_assign_surf_class_to_region(
    struct sym_entry *sc_sym,
    struct region *rgn);
