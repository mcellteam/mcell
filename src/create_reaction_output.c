/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2014 by *
 * The Salk Institute for Biological Studies and *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or *
 * modify it under the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 *
 * of the License, or (at your option) any later version. *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *
 * GNU General Public License for more details. *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License *
 * along with this program; if not, write to the Free Software *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 *USA. *
 *                                                                                 *
 ***********************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "logging.h"
#include "create_reaction_output.h"


/* static function definitions */
static long long pick_buffer_size(MCELL_STATE *state, struct output_block *obp,
  long long n_output);

/**************************************************************************
 new_output_column:
    Create a new output column for an output set.

 In: Nothing
 Out: output column, or NULL if allocation fails
**************************************************************************/
struct output_column *new_output_column() {

  struct output_column *oc;
  oc = CHECKED_MALLOC_STRUCT(struct output_column,
                             "reaction data output column");
  if (oc == NULL)
    return NULL;

  oc->data_type = COUNT_UNSET;
  oc->initial_value = 0.0;
  oc->buffer = NULL;
  oc->expr = NULL;
  oc->next = NULL;

  return oc;
}

/**************************************************************************
 new_output_block:
    Allocate a new reaction data output block, with a specified buffer size.

 In: buffersize: requested buffer size
 Out: output block, or NULL if an error occurs
**************************************************************************/
struct output_block *new_output_block(int buffersize) {

  struct output_block *obp;
  obp =
      CHECKED_MALLOC_STRUCT(struct output_block, "reaction data output block");
  if (obp == NULL)
    return NULL;

  obp->t = 0.0;
  obp->timer_type = OUTPUT_BY_STEP;
  obp->step_time = FOREVER;
  obp->time_list_head = NULL;
  obp->time_now = NULL;
  obp->buffersize = 0;
  obp->trig_bufsize = 0;
  obp->buf_index = 0;
  obp->data_set_head = NULL;

  /* COUNT buffer size might get modified later if there isn't that much to
   * output */
  /* TRIGGER buffer size won't get modified since we can't know what to expect
   */
  obp->buffersize = buffersize;
  obp->trig_bufsize = obp->buffersize;

  obp->time_array = CHECKED_MALLOC_ARRAY(double, obp->buffersize,
                                         "reaction data output times array");
  if (obp->time_array == NULL) {
    free(obp);
    return NULL;
  }

  return obp;
}

/**************************************************************************
 set_reaction_output_timer_step:
    Set the output timer for reaction data output to a time step.

 In: parse_state: parser state
     obp:  output block whose timer is to be set
     step: time step interval
 Out: output timer is updated
**************************************************************************/
void
set_reaction_output_timer_step(MCELL_STATE *state, struct output_block *obp,
  double step) {

  long long output_freq;
  obp->timer_type = OUTPUT_BY_STEP;

  obp->step_time = step;
  output_freq = (long long)(obp->step_time / state->time_unit);

  /* Clip the step time to a good range */
  if (output_freq > state->iterations && output_freq > 1) {
    output_freq = (state->iterations > 1) ? state->iterations : 1;
    obp->step_time = output_freq * state->time_unit;
    if (state->notify->invalid_output_step_time != WARN_COPE)
      mcell_log("Output step time too long\n\tSetting output "
                "step time to %g microseconds\n", obp->step_time * 1.0e6);
  } else if (output_freq < 1) {
    output_freq = 1;
    obp->step_time = output_freq * state->time_unit;
    if (state->notify->invalid_output_step_time != WARN_COPE)
      mcell_log("Output step time too short\n\tSetting output "
                "step time to %g microseconds\n", obp->step_time * 1.0e-6);
  }

  /* Pick a good buffer size */
  long long n_output;
  if (state->chkpt_iterations)
    n_output =
        (long long)(state->chkpt_iterations / output_freq + 1);
  else
    n_output = (long long)(state->iterations / output_freq + 1);
  obp->buffersize = pick_buffer_size(state, obp, n_output);

  no_printf("Default output step time definition:\n");
  no_printf("  output step time = %g\n", obp->step_time);
  no_printf("  output buffersize = %u\n", obp->buffersize);
}


/**************************************************************************
 pick_buffer_size:
    Choose an appropriate output buffer size for our reaction output data,
    based on the total number of outputs expected and the requested buffer
    size.

 In: parse_state: parser state
     obp:  output block whose buffer_size to set
     n_output: maximum number of outputs expected
 Out: 0 on success, 1 on failure
**************************************************************************/
long long pick_buffer_size(MCELL_STATE *state, struct output_block *obp,
  long long n_output) {
  if (state->chkpt_iterations)
    return min3ll(state->chkpt_iterations - state->start_time + 1,
                  n_output, obp->buffersize);
  else
    return min3ll(state->iterations - state->start_time + 1,
                  n_output, obp->buffersize);
}

/**************************************************************************
 set_reaction_output_timer_iterations:
    Set the output timer for reaction data output to a list of iterations.

 In: parse_state: parser state
     obp:  output block whose timer is to be set
     step_values: list of iterations
 Out: 0 on success, 1 on failure; output timer is updated
**************************************************************************/
int set_reaction_output_timer_iterations(MCELL_STATE *state,
  struct output_block *obp, struct num_expr_list_head *step_values) {
  obp->timer_type = OUTPUT_BY_ITERATION_LIST;
  obp->buffersize =
      pick_buffer_size(state, obp, step_values->value_count);
  if (step_values->shared) {
    obp->time_list_head =
        mcell_copysort_numeric_list(step_values->value_head);
    if (obp->time_list_head == NULL)
      return 1;
  } else {
    mcell_sort_numeric_list(step_values->value_head);
    obp->time_list_head = step_values->value_head;
  }
  obp->time_now = NULL;
  return 0;
}

/**************************************************************************
 set_reaction_output_timer_times:
    Set the output timer for reaction data output to a list of times.

 In: parse_state: parser state
     obp:  output block whose timer is to be set
     nstep: number of times
     step_values: list of times
 Out: output timer is updated
**************************************************************************/
int set_reaction_output_timer_times(MCELL_STATE *state, struct output_block *obp,
  struct num_expr_list_head *step_values) {
  obp->timer_type = OUTPUT_BY_TIME_LIST;
  obp->buffersize =
      pick_buffer_size(state, obp, step_values->value_count);
  if (step_values->shared) {
    obp->time_list_head =
        mcell_copysort_numeric_list(step_values->value_head);
    if (obp->time_list_head == NULL)
      return 1;
  } else {
    mcell_sort_numeric_list(step_values->value_head);
    obp->time_list_head = step_values->value_head;
  }
  obp->time_now = NULL;
  return 0;
}


/**************************************************************************
 output_block_finalize:
    Finalizes a reaction data output block, checking for errors, and allocating
    the output buffer.

 In: parse_state: parser state
     obp:  the output block to finalize
 Out: 0 on success, 1 on failure
**************************************************************************/
int output_block_finalize(MCELL_STATE *state, struct output_block *obp) {
  struct output_set *os1;
  for (os1 = obp->data_set_head; os1 != NULL; os1 = os1->next) {
    /* Check for duplicated filenames */
    struct output_set *os2;
    for (os2 = os1->next; os2 != NULL; os2 = os2->next) {
      if (strcmp(os1->outfile_name, os2->outfile_name) == 0) {
        mcell_error("COUNT statements in the same reaction data "
                    "output block should have unique output file "
                    "names (\"%s\" appears more than once)",
                     os1->outfile_name);
        return MCELL_FAIL;
      }
    }

    /* Allocate buffers */
    struct output_column *oc;
    for (oc = os1->column_head; oc != NULL; oc = oc->next) {
      switch (oc->expr->expr_flags & OEXPR_TYPE_MASK) {
      case OEXPR_TYPE_INT:
        oc->data_type = COUNT_INT;
        oc->buffer = CHECKED_MALLOC_ARRAY(int, obp->buffersize,
                                          "reaction data output buffer");
        break;

      case OEXPR_TYPE_DBL:
        oc->data_type = COUNT_DBL;
        oc->buffer = CHECKED_MALLOC_ARRAY(double, obp->buffersize,
                                          "reaction data output buffer");
        break;

      case OEXPR_TYPE_TRIG:
        oc->data_type = COUNT_TRIG_STRUCT;
        oc->buffer =
            CHECKED_MALLOC_ARRAY(struct output_trigger_data, obp->trig_bufsize,
                                 "reaction data output buffer");
        break;

      default:
        mcell_error("Could not figure out what type of count data to store");
        return MCELL_FAIL;
      }
      if (oc->buffer == NULL)
        return MCELL_FAIL;
    }
  }

  return MCELL_SUCCESS;
}
