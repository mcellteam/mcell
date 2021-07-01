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

/* status of libMCell API calls */
typedef int MCELL_STATUS;

#define MCELL_SUCCESS 0
#define MCELL_FAIL 1

/* state of mcell simulation */
typedef struct volume MCELL_STATE;

struct num_expr_list_head {
  struct num_expr_list *value_head;
  struct num_expr_list *value_tail;
  int value_count;
  int shared;
  
  // MCell4 - original values used to create the list
  bool start_end_step_set;
  double start;
  double end;
  double step;  
};

void mcell_set_seed(MCELL_STATE *state, int seed);
void mcell_set_with_checks_flag(MCELL_STATE *state, int value);
void mcell_set_randomize_smol_pos(MCELL_STATE *state, int value);

MCELL_STATE *mcell_create(void);

MCELL_STATUS mcell_init_state(MCELL_STATE *state);

MCELL_STATUS mcell_init_simulation(MCELL_STATE *state);

MCELL_STATUS mcell_init_read_checkpoint(MCELL_STATE *state);

MCELL_STATUS mcell_init_output(MCELL_STATE *state);

MCELL_STATUS mcell_set_partition(MCELL_STATE *state, int dim,
                                 struct num_expr_list_head *head);

MCELL_STATUS mcell_set_time_step(MCELL_STATE *state, double step);

MCELL_STATUS mcell_set_iterations(MCELL_STATE *state, long long iterations);

MCELL_STATUS mcell_silence_notifications(MCELL_STATE *state);
MCELL_STATUS mcell_enable_notifications(MCELL_STATE *state);
MCELL_STATUS mcell_silence_warnings(MCELL_STATE *state);
MCELL_STATUS mcell_enable_warnings(MCELL_STATE *state);
