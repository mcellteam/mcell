/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
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

#pragma once

MCELL_STATUS mcell_create_viz_output(MCELL_STATE *state, char *filename,
                                     struct mcell_species *mol_viz_list,
                                     long long start, long long end,
                                     long long step);

void mcell_new_viz_output_block(struct viz_output_block *vizblk);

struct frame_data_list *
mcell_create_viz_frame(int time_type, int type,
                       struct num_expr_list *iteration_list);

int mcell_set_molecule_viz_state(struct viz_output_block *vizblk,
                                 struct species *specp, int viz_state);
