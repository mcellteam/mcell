#ifndef VIZ_OUTPUT_H
#define VIZ_OUTPUT_H

#include "mcell_structs.h"

/* Header file for visualization output routines */

int update_frame_data_list(struct viz_output_block *vizblk);
int init_frame_data_list(struct viz_output_block *vizblk);
int finalize_viz_output(struct viz_output_block *vizblk);

#endif
