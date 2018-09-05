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
%include "typemaps.i"

typedef struct sym_entry {
  struct sym_entry *next; /* Chain to next symbol in this bin of the hash */
  int sym_type;           /* Symbol Type */
  char *name;             /* Name of symbol*/
  void *value;            /* Stored value, cast by sym_type */
} mcell_symbol;


struct mcell_species_spec {
  char *name;
  double D;
  int is_2d;               // 3D = 0; 2D = 1
  double custom_time_step; // default is 0.0
  int target_only;         // default is 0
  double max_step_length;  // default is 0.0
  double space_step;
};

struct mcell_species {
  struct mcell_species *next;
  struct sym_entry *mol_type;
  short orient_set;
  short orient;
};

struct mcell_species_list {
  struct mcell_species *mol_type_head;
  struct mcell_species *mol_type_tail;
};

%typemap(in) mcell_symbol **species_ptr (mcell_symbol *temp) {
  $1 = &temp;
}

%typemap(argout) struct sym_entry **species_ptr {
  %set_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

MCELL_STATUS mcell_create_species(MCELL_STATE *state,
                                  struct mcell_species_spec *species,
                                  mcell_symbol **species_ptr);

struct mcell_species *
mcell_add_to_species_list(mcell_symbol *species_ptr, bool is_oriented,
                          int orientation, struct mcell_species *species_list);

void mcell_delete_species_list(struct mcell_species *species);

int new_mol_species(MCELL_STATE *state, char *name, struct sym_entry **sym_ptr);
