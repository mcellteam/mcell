/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

MCELL_STATUS mcell_create_viz_output(MCELL_STATE *state, char *filename,
                                     struct mcell_species *mol_viz_list,
                                     long long start, long long end,
                                     long long step, bool ascii_output);

void mcell_new_viz_output_block(struct viz_output_block *vizblk);

struct frame_data_list *
mcell_create_viz_frame(int time_type, int type,
                       struct num_expr_list *iteration_list);

int mcell_set_molecule_viz_state(struct viz_output_block *vizblk,
                                 struct species *specp, int viz_state);
