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

#include "config.h"

#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "sym_table.h"
#include "logging.h"
#include "react_output.h"
#include "mcell_misc.h"
#include "mcell_react_out.h"
#include "mdlparse_util.h"

#include "dyngeom_parse_extras.h"
#include "strfunc.h"
#include "count_util.h"

/* static helper functions */
static struct output_column *new_output_column(void);

static struct output_block *new_output_block(int buffersize);

static void set_reaction_output_timer_step(MCELL_STATE *state,
                                           struct output_block *obp,
                                           double step);

static int
set_reaction_output_timer_iterations(MCELL_STATE *state,
                                     struct output_block *obp,
                                     struct num_expr_list_head *step_values);

static int
set_reaction_output_timer_times(MCELL_STATE *state, struct output_block *obp,
                                struct num_expr_list_head *step_values);

static int output_block_finalize(struct output_block *obp);

static long long pick_buffer_size(MCELL_STATE *state, struct output_block *obp,
                                  long long n_output);

static struct output_column *
get_counter_trigger_column(MCELL_STATE *state, const char *counter_name,
                           int column_id);

/*************************************************************************
 mcell_get_count:
  Get the count of a molecule in a certain region..
 
 In:  mol_name - molecule you want the count of
      reg_name - region where you want the count
      world - the instance world of the object volume.
Out: int mol_count - count of molecule in the region
*************************************************************************/
int mcell_get_count(char *mol_name, char *reg_name, struct volume *world) {
  // Get the hash values for mol and reg---------------// 

  // Get hash value for molecule  

  struct sym_entry *mol_sym = NULL;  
  mol_sym = retrieve_sym(mol_name, world->mol_sym_table);

  // Make sure mol_sym has been initialized. If not, return error code.
  if (mol_sym == NULL)
    return -5;

  // Cast mol_sym (sym_entry, void pointer) to a species pointer
  struct species *mol = mol_sym->value;
  // Get hash value for molecule 
  u_int mol_hashval = mol->hashval;

  // Get hash value for region

  struct sym_entry *reg_sym = NULL;  
  reg_sym = retrieve_sym(reg_name, world->reg_sym_table);

  // Make sure reg_sym has been initialized, if not return garbage

  if (reg_sym == NULL)
    return -6;  
 
  // Cast mol_sym (sym_entry,void pointer) to a species pointer 
  struct region *reg = reg_sym->value;
  // Get hash value for region
  u_int reg_hashval = reg-> hashval;

  //---------------------------------------------------//

  // Use the hash values to get the molecule count in the region from the count
  // hash //

  // Combine hash values for molecule in that region
  int hash_bin = (mol_hashval + reg_hashval) & world->count_hashmask;
 
  // Make sure hash_bin exists
  if (world->count_hash[hash_bin] == NULL)
    return -7;

  // Get the count of molecule in and on the region  
  int mol_count_vol = world->count_hash[hash_bin]->data.move.n_enclosed;
  int mol_count_sur = world->count_hash[hash_bin]->data.move.n_at;

  return mol_count_vol + mol_count_sur;

}


/*************************************************************************
 mcell_new_output_request:
    Create a new output request.

 In:  state: MCell state
      target: what are we counting
      orientation: how is it oriented?
      location: where are we counting?
      report_flags: what type of events are we counting?
 Out: output request item, or NULL if an error occurred
*************************************************************************/
struct output_request *mcell_new_output_request(MCELL_STATE *state,
                                                struct sym_entry *target,
                                                short orientation,
                                                struct sym_entry *location,
                                                struct periodic_image *img,
                                                int report_flags) {
  struct output_request *orq;
  struct output_expression *oe;

  orq = CHECKED_MEM_GET(state->outp_request_mem, "count request");
  if (orq == NULL)
    return NULL;

  oe = new_output_expr(state->oexpr_mem);
  if (oe == NULL) {
    mem_put(state->outp_request_mem, orq);
    mcell_allocfailed("Failed to allocate a count expression.");
    return NULL;
  }
  orq->next = NULL;
  orq->requester = oe;
  orq->count_target = target;
  orq->count_orientation = orientation;
  orq->count_location = location;
  orq->report_type = report_flags;
  orq->periodic_box = img;

  oe->left = orq;
  oe->oper = '#';
  oe->expr_flags = OEXPR_LEFT_REQUEST;
  struct sym_entry *sym = NULL;
  if (location) {
    char *name = location->name;
    // Counting in/on a region. XXX: Using strchr seems inefficient
    if (strchr(location->name, ',')) {
      sym = retrieve_sym(name, state->reg_sym_table);
    }
    // Counting in/on an object
    else {
      sym = retrieve_sym(name, state->obj_sym_table);
    }
  }

  // If the object/region will exist at some point in the future, but not at
  // the beginning of the simulation.
  if (sym && !(is_object_instantiated(sym, state->root_instance))) 
    oe->expr_flags = OEXPR_TYPE_UNDEF;
  else if (orq->report_type & REPORT_TRIGGER)
    oe->expr_flags |= OEXPR_TYPE_TRIG;
  else if ((orq->report_type & REPORT_TYPE_MASK) != REPORT_CONTENTS)
    oe->expr_flags |= OEXPR_TYPE_DBL;
  else
    oe->expr_flags |= OEXPR_TYPE_INT;
  return orq;
}

/******************************************************************************
 *
 * mcell_create_count creates a single count expression and returns it as a
 * output_column_list.
 * Inputs are:
 *    - symbol table entry for target (molecule or reaction)
 *    - orientation for molecule counts
 *    - symbol table entry for count location (NULL implies WORLD)
 *    - report flags (REPORT_WORLD, REPORT_CONTENTS, ...)
 *    - custom header (or NULL if not wanted)
 *    - pointer to empty count list to which count expression will be added
 *
 *****************************************************************************/
MCELL_STATUS
mcell_create_count(MCELL_STATE *state, struct sym_entry *target,
                   short orientation, struct sym_entry *location,
                   int report_flags, char *custom_header,
                   struct output_column_list *count_list) {

  struct output_request *output_A = NULL;
  if ((output_A = mcell_new_output_request(state, target, orientation, location,
    NULL, report_flags)) == NULL) {
    return MCELL_FAIL;
  }
  output_A->next = state->output_request_head;
  state->output_request_head = output_A;

  return mcell_prepare_single_count_expr(count_list, output_A->requester,
                                         custom_header);
}

/*************************************************************************
 mcell_create_new_output_set
    Create a new output set. Here output set refers to a count/trigger
    block which goes to a single data output file.

 In:  comment: textual comment describing the data set or NULL
      exact_time: request exact_time output for trigger statements
      col_head: head of linked list of output columns
      file_flags: file creation flags for output file
      outfile_name: name of output file
 Out: output request item, or NULL if an error occurred
*************************************************************************/
struct output_set *mcell_create_new_output_set(char *comment, int exact_time,
                                               struct output_column *col_head,
                                               int file_flags,
                                               char *outfile_name) {

  struct output_set *os =
      CHECKED_MALLOC_STRUCT(struct output_set, "reaction data output set");
  if (os == NULL) {
    return NULL;
  }

  os->outfile_name = CHECKED_STRDUP(outfile_name, "count outfile_name");
  os->file_flags = file_flags;
  os->exact_time_flag = exact_time;
  os->chunk_count = 0;
  os->block = NULL;
  os->next = NULL;

  struct output_column *oc = col_head;
  os->column_head = oc;

  if (comment == NULL)
    os->header_comment = NULL;
  else if (comment[0] == '\0')
    os->header_comment = "";
  else {
    os->header_comment = strdup(comment);
    if (os->header_comment == NULL) {
      free(os);
      return NULL;
    }
  }

  for (; oc != NULL; oc = oc->next)
    oc->set = os;

  if (check_reaction_output_file(os)) {
    free(os);
    return NULL;
  }

  return os;
}

/*****************************************************************************
 *
 * mcell_prepare_single_count_expression prepares a count expression for
 * inclusion in an output set
 *
 *****************************************************************************/
MCELL_STATUS
mcell_prepare_single_count_expr(struct output_column_list *list,
                                struct output_expression *expr,
                                char *custom_header) {
  list->column_head = NULL;
  list->column_tail = NULL;

  if (custom_header != NULL) {
    expr->title = custom_header;
  }

  /* If we have a list of results, go through to build column stack */
  struct output_expression *oe;
  struct output_column *oc;
  for (oe = first_oexpr_tree(expr); oe != NULL; oe = next_oexpr_tree(oe)) {
    if ((oc = new_output_column()) == NULL)
      return MCELL_FAIL;

    if (!list->column_head)
      list->column_head = list->column_tail = oc;
    else
      list->column_tail = list->column_tail->next = oc;

    oc->expr = oe;
    set_oexpr_column(oe, oc);
  }

  return MCELL_SUCCESS;
}

/*****************************************************************************
 *
 * mcell_add_reaction_output_block creates a new reaction data output block
 * and adds it to the world.
 *
 *****************************************************************************/
MCELL_STATUS
mcell_add_reaction_output_block(MCELL_STATE *state,
                                struct output_set_list *osets, int buffer_size,
                                struct output_times_inlist *otimes) {

  struct output_block *obp;
  struct output_set *os;

  if ((obp = new_output_block(buffer_size)) == NULL)
    return 1;

  if (otimes->type == OUTPUT_BY_STEP)
    set_reaction_output_timer_step(state, obp, otimes->step);
  else if (otimes->type == OUTPUT_BY_ITERATION_LIST) {
    if (set_reaction_output_timer_iterations(state, obp, &otimes->values)) {
      free(obp);
      return MCELL_FAIL;
    }
  } else if (otimes->type == OUTPUT_BY_TIME_LIST) {
    if (set_reaction_output_timer_times(state, obp, &otimes->values)) {
      free(obp);
      return MCELL_FAIL;
    }
  } else {
    mcell_error("Internal error: Invalid output timer def (%d)", otimes->type);
    free(obp);
    return MCELL_FAIL;
  }
  obp->data_set_head = osets->set_head;
  for (os = obp->data_set_head; os != NULL; os = os->next)
    os->block = obp;
  if (output_block_finalize(obp))
    return 1;
  obp->next = state->output_block_head;
  state->output_block_head = obp;
  return MCELL_SUCCESS;
}

/************************************************************************
 *
 * function for retrieving the current value of a given count
 * expression
 *
 * The call expects:
 *
 * - MCELL_STATE
 * - counter_name: a string containing the name of the count statement to
 *   be retrieved. Currently, the name is identical to the full path to which
 *   the corresponding reaction output will be written but this may change
 *   in the future
 * - column: int describing the column to be retrieved
 * - count_data: a *double which will receive the actual value
 * - count_data_type: a *count_type_t which will receive the type of the
 *   data (for casting of count_data)
 *
 * NOTE: This function can be called anytime after the
 *       REACTION_DATA_OUTPUT has been either parsed or
 *       set up with API calls.
 *
 * Returns 1 on error and 0 on success
 *
 ************************************************************************/
MCELL_STATUS
mcell_get_counter_value(MCELL_STATE *state,
                        const char *counter_name,
                        int column_id,
                        double *count_data,
                        enum count_type_t *count_data_type) {
  struct output_column *column = NULL;
  if ((column = get_counter_trigger_column(state, counter_name, column_id)) ==
      NULL) {
    return MCELL_FAIL;
  }

  // if we happen to encounter trigger data we bail
  if (column->buffer[0].data_type == COUNT_TRIG_STRUCT) {
    return MCELL_FAIL;
  }

  // evaluate the expression and retrieve it
  eval_oexpr_tree(column->expr, 1);
  *count_data = (double)column->expr->value;
  *count_data_type = column->buffer[0].data_type;

  return MCELL_SUCCESS;
}

/******************************************************************************
 *
 * static helper functions
 *
 *****************************************************************************/

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
void set_reaction_output_timer_step(MCELL_STATE *state,
                                    struct output_block *obp, double step) {

  long long output_freq;
  obp->timer_type = OUTPUT_BY_STEP;

  obp->step_time = step;
  output_freq = (long long)(obp->step_time / state->time_unit);

  /* Clip the step time to a good range */
  if (output_freq > state->iterations && output_freq > 1) {
    output_freq = (state->iterations > 1) ? state->iterations : 1;
    obp->step_time = output_freq * state->time_unit;
    if (state->notify->invalid_output_step_time != WARN_COPE)
      mcell_warn("output step time too long.\n  Setting output step time to "
                 "%g seconds.", obp->step_time);
  } else if (output_freq < 1) {
    output_freq = 1;
    obp->step_time = output_freq * state->time_unit;
    if (state->notify->invalid_output_step_time != WARN_COPE)
      mcell_warn("output step time too short.\n  Setting output step time to "
                 "%g seconds.", obp->step_time);
  }

  /* Pick a good buffer size */
  long long n_output;
  if (state->chkpt_iterations)
    n_output = (long long)(state->chkpt_iterations / output_freq + 1);
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
    return min3ll(state->chkpt_iterations - state->start_iterations + 1, n_output,
                  obp->buffersize);
  else
    return min3ll(state->iterations - state->start_iterations + 1, n_output,
                  obp->buffersize);
}

/**************************************************************************
 set_reaction_output_timer_iterations:
    Set the output timer for reaction data output to a list of iterations.

 In: parse_state: parser state
     obp:  output block whose timer is to be set
     step_values: list of iterations
 Out: 0 on success, 1 on failure; output timer is updated
**************************************************************************/
int
set_reaction_output_timer_iterations(MCELL_STATE *state,
                                     struct output_block *obp,
                                     struct num_expr_list_head *step_values) {
  obp->timer_type = OUTPUT_BY_ITERATION_LIST;
  obp->buffersize = pick_buffer_size(state, obp, step_values->value_count);
  if (step_values->shared) {
    obp->time_list_head = mcell_copysort_numeric_list(step_values->value_head);
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
int set_reaction_output_timer_times(MCELL_STATE *state,
                                    struct output_block *obp,
                                    struct num_expr_list_head *step_values) {
  obp->timer_type = OUTPUT_BY_TIME_LIST;
  obp->buffersize = pick_buffer_size(state, obp, step_values->value_count);
  if (step_values->shared) {
    obp->time_list_head = mcell_copysort_numeric_list(step_values->value_head);
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

 In: obp:  the output block to finalize
 Out: 0 on success, 1 on failure
**************************************************************************/
int output_block_finalize(struct output_block *obp) {
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
      // Counting on meshes that don't exist at the beginning of the sim.
      case OEXPR_TYPE_UNDEF:
        oc->buffer = CHECKED_MALLOC_ARRAY(struct output_buffer,
                                          obp->buffersize,
                                          "reaction data output buffer");
        for (u_int i = 0; i < obp->buffersize; ++i) {
          oc->buffer[i].data_type = COUNT_UNSET;
          oc->buffer[i].val.cval = 'X';
        }
        break;
      case OEXPR_TYPE_INT:
        oc->buffer = CHECKED_MALLOC_ARRAY(struct output_buffer,
                                          obp->buffersize,
                                          "reaction data output buffer");
        for (u_int i = 0; i < obp->buffersize; ++i) {
          oc->buffer[i].data_type = COUNT_INT;
          oc->buffer[i].val.ival = 0;
        }
        break;

      case OEXPR_TYPE_DBL:
        oc->buffer = CHECKED_MALLOC_ARRAY(struct output_buffer,
                                          obp->buffersize,
                                          "reaction data output buffer");
        for (u_int i = 0; i < obp->buffersize; ++i) {
          oc->buffer[i].data_type = COUNT_DBL;
          oc->buffer[i].val.dval = 0.0;
        }
        break;

      case OEXPR_TYPE_TRIG:
        oc->buffer = CHECKED_MALLOC_ARRAY(struct output_buffer,
                                          obp->trig_bufsize,
                                          "reaction data output buffer");
        for (u_int i = 0; i < obp->trig_bufsize; ++i) {
          oc->buffer[i].data_type = COUNT_TRIG_STRUCT;
          oc->buffer[i].val.tval = CHECKED_MALLOC_STRUCT(
              struct output_trigger_data,
              "reaction data output buffer");
          oc->buffer[i].val.tval->name = NULL;
        }
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

/************************************************************************
 *
 * get_counter_trigger_column retrieves the output_column corresponding
 * to a given count or trigger statement.
 *
 ************************************************************************/
struct output_column *get_counter_trigger_column(MCELL_STATE *state,
                                                 const char *counter_name,
                                                 int column_id) {
  // retrieve the counter for the requested counter_name
  struct sym_entry *counter_sym =
      retrieve_sym(counter_name, state->counter_by_name);
  if (counter_sym == NULL) {
    mcell_log("Failed to retrieve symbol for counter %s.", counter_name);
    return NULL;
  }
  struct output_set *counter = (struct output_set *)(counter_sym->value);

  // retrieve the requested column
  struct output_column *column = counter->column_head;
  int count = 0;
  while (count < column_id && column != NULL) {
    count++;
    column = column->next;
  }
  if (count != column_id || column == NULL) {
    return NULL;
  }

  return column;
}


