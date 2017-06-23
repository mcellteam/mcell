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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>

#include "logging.h"
#include "sched_util.h"
#include "mcell_structs.h"
#include "react_output.h"
#include "mdlparse_util.h"
#include "strfunc.h"

// XXX: This global state should be removed. Currently
// we need it for cleanup via signals.
static struct volume *global_state;

/**************************************************************************
truncate_output_file:
  In: filename string
      value that we will start outputting to the file
  Out: 0 if file preparation is successful, 1 if not.  The file is
       truncated at start of the line containing the first entry
       greater than or equal to the value to be printed out.
**************************************************************************/

int truncate_output_file(char *name, double start_value) {
  FILE *f = NULL;

  /* Check if the file exists */
  struct stat fs;
  int i = stat(name, &fs);
  if (i == -1) {
    mcell_perror(errno, "Failed to stat reaction data output file '%s' in "
                        "preparation for truncation.",
                 name);
  }
  if (fs.st_size == 0)
    return 0; /* File already is empty */

  /* Set the buffer size */
  off_t bsize;
  if (fs.st_size < (1 << 20)) {
    bsize = fs.st_size;
  } else
    bsize = (1 << 20);

  /* Allocate a buffer for the file */
  char *buffer= CHECKED_MALLOC_ARRAY_NODIE(
    char, bsize + 1, "reaction data file scan buffer");
  if (buffer == NULL)
    goto failure;

  /* Open the file */
  f = fopen(name, "r+");
  if (!f) {
    mcell_perror(
        errno, "Failed to open reaction data output file '%s' for truncation.",
        name);
    /*goto failure;*/
  }

  /* Iterate over the entire file */
  int where = 0; /* Byte offset in file */
  int start = 0; /* Byte offset in buffer */
  while (ftell(f) != fs.st_size) {
    /* Refill the buffer */
    long long n = (long long)fread(buffer + start, 1, bsize - start, f);

    /* Until the current buffer runs dry */
    int ran_out = 0;
    i = 0;
    n += start;
    int lf = 0;
    while (!ran_out) {
      /* Skip leading horizontal whitespace */
      while (i < n && (buffer[i] == ' ' || buffer[i] == '\t'))
        i++;

      /* Scan over leading numeric characters */
      int j = i;
      for (;
           j < n && (isdigit(buffer[j]) || strchr("eE-+.", buffer[j]) != NULL);
           j++) {
      }

      /* If we had a leading number... */
      if (j > i && j < (n - 1)) {
        char *done = NULL;

        /* Parse and validate the number */
        buffer[j] = '\0';
        double my_value = strtod(buffer + i, &done) + EPS_C;

        /* If it was a valid number and it was >= our start time */
        if (done != buffer + i && my_value >= start_value) {
          if (fseek(f, 0, SEEK_SET)) {
            mcell_perror(
                errno,
                "Failed to seek to beginning of reaction data output file '%s'",
                name);
            /*goto failure;*/
          }

          if (ftruncate(fileno(f), where + lf)) {
            mcell_perror(errno,
                         "Failed to truncate reaction data output file '%s'",
                         name);
            /*goto failure;*/
          }
          fclose(f);
          free(buffer);
          return 0;
        }
      }

      /* We will keep this line.  Scan until we hit a newline */
      for (i = j; i < n && buffer[i] != '\n' && buffer[i] != '\r'; i++) {
      }

      /* If we're not at the end of the buffer, scan over all crs/lfs. */
      if (i < n) {
        for (j = i; j < n && (buffer[j] == '\n' || buffer[j] == '\r'); j++) {
        }
        lf = j;
        i = j; /* If we have run out, we'll catch it next time through */
      }

      /* If we hit 'n' and the last read was only partial (i.e. EOF), break out
         */
      else if (n < bsize) {
        ran_out = 1;
      }

      /* Last read was a full read and we've scanned the whole buffer */
      else {
        /* If we need to keep some data, resituate it in the buffer */
        if (lf > start) {
          /* discard all but n - lf bytes */
          memmove(buffer, buffer + lf, n - lf);
          where += lf;
          start = n - lf;
        } else {
          /* discard all bytes */
          where += n;
          start = 0;
        }
        lf = 0;
        ran_out = 1;
      }
    }
  }
  fclose(f);
  free(buffer);
  return 0;

failure:
  if (f != NULL)
    fclose(f);
  if (buffer != NULL)
    free(buffer);
  return 1;
}

/**************************************************************************
emergency_output:
  In: No arguments.
  Out: Number of errors encountered while trying to make an emergency
       dump of output buffers (0 indicates emergency save went fine).
       The code assumes a memory error, dumps any memory it can to
       try to make enough room to create a file to save results held
       in buffers.
  Note: The simulation is COMPLETELY TRASHED after this function is
        called.  Do NOT use this in normal operation.  Do NOT try to
    continue running after this function is called.  This function
    will deallocate all molecules, walls, etc. to try to recover
    memory!  You should only print messages and exit after running
    this function.
**************************************************************************/
static int emergency_output(struct volume *world) {
  struct storage_list *mem;

  /* PANIC--delete everything we can get our pointers on! */
  delete_mem(world->coll_mem);
  delete_mem(world->exdv_mem);
  for (mem = world->storage_head; mem != NULL; mem = mem->next) {
    delete_mem(mem->store->list);
    delete_mem(mem->store->mol);
    delete_mem(mem->store->smol);
    delete_mem(mem->store->face);
    delete_mem(mem->store->join);
    delete_mem(mem->store->grids);
    delete_mem(mem->store->regl);
  }
  delete_mem(world->storage_allocator);

  return flush_reaction_output(world);
}

/**************************************************************************
 emergency_output_hook_enabled:
    Flag to disable emergency output hook when the program exits successfully.
**************************************************************************/
int emergency_output_hook_enabled = 1;

/**************************************************************************
 emergency_output_hook:
    This is an atexit hook to flush reaction output to disk in case an error is
    occurred.  Set emergency_output_hook_enabled to 0 to prevent it from being
    called (say, on successful exit).

  In: No arguments.
  Out: None.

**************************************************************************/
static void emergency_output_hook(void) {
  if (emergency_output_hook_enabled) {
    /* Disable the emergency output hook in case a signal is received while
     * producing emergency output. */
    emergency_output_hook_enabled = 0;

    int n_errors = emergency_output(global_state);
    if (n_errors == 0)
      mcell_warn("Reaction output was successfully flushed to disk.");
    else if (n_errors == 1)
      mcell_warn("An error occurred while flushing reaction output to disk.");
    else
      mcell_warn("%d errors occurred while flushing reaction output to disk.",
                 n_errors);
  }
}

/**************************************************************************
 emergency_output_signal_handler:
    This is a signal handler to catch any abnormal termination signals, and
    flush the reaction output data (and anything else which should be flushed).

  In: No arguments.
  Out: None.
**************************************************************************/
static void emergency_output_signal_handler(int signo)
    __attribute__((noreturn));

static void emergency_output_signal_handler(int signo) {

  if (emergency_output_hook_enabled) {
    emergency_output_hook_enabled = 0;

    int n_errors = flush_reaction_output(global_state);
    if (n_errors == 0)
      mcell_error_raw("Reaction output was successfully flushed to disk.\n");
    else if (n_errors == 1)
      mcell_error_raw(
          "An error occurred while flushing reaction output to disk.\n");
    else
      mcell_error_raw(
          "%d errors occurred while flushing reaction output to disk.\n",
          n_errors);
  }
  raise(signo);

  /* We shouldn't get here, but we might if, for instance, SA_NODEFER is
   * unsupported. */
  _exit(128 + signo);
}

/**************************************************************************
 install_emergency_output_signal_handler:
    This installs a handler for a single signal which will print out a sensible
    message, then try to flush as much output data to disk as possible before
    dying.

  In: No arguments.
  Out: None.
**************************************************************************/
static void install_emergency_output_signal_handler(int signo) {
#ifdef _WIN32 /* fixme: for Windows do a better job than just signal(), need   \
                 to find out what other things the *nix version is doing */
  signal(signo, &emergency_output_signal_handler);
#else
  struct sigaction sa, saPrev;
  sa.sa_sigaction = NULL;
  sa.sa_handler = &emergency_output_signal_handler;
  sa.sa_flags = SA_RESTART | SA_RESETHAND | SA_NODEFER;
  sigfillset(&sa.sa_mask);

  if (sigaction(signo, &sa, &saPrev) != 0)
    mcell_warn("Failed to install emergency output signal handler.");
#endif
}

/**************************************************************************
 install_emergency_output_hooks:
    Installs all relevant hooks for catching invalid program termination and
    flushing output to disk, where possible.

  In: No arguments.
  Out: None.
**************************************************************************/
void install_emergency_output_hooks(struct volume *world) {
  global_state = world;

  if (atexit(&emergency_output_hook) != 0)
    mcell_warn("Failed to install emergency output hook.");

  install_emergency_output_signal_handler(
      SIGILL); /* not generated on Windows but can be raised manually */
  install_emergency_output_signal_handler(SIGABRT);
  install_emergency_output_signal_handler(SIGFPE);
  install_emergency_output_signal_handler(SIGSEGV);
#ifdef SIGBUS
  install_emergency_output_signal_handler(SIGBUS);
#endif
}

/*************************************************************************
add_trigger_output:
   In: counter of thing that just happened (trigger of some sort)
       request structure saying who wanted to know that it happened
       number of times that thing happened
   Out: The event is added to the buffer of the requester, and the
        buffer is written and flushed if the buffer is full.
   Note: Hits are assumed to happen with only one molecule; positive
        means the front face was hit, negative means the back was hit.
        This is reported as orientation instead.
*************************************************************************/
void add_trigger_output(struct volume *world, struct counter *c,
                        struct output_request *ear, int n, short flags) {

  struct output_column *first_column;
  first_column = ear->requester->column->set->column_head;

  int idx = (int)first_column->initial_value;

  struct output_trigger_data *otd;
  otd = first_column->buffer[idx].val.tval;

  if (first_column->set->block->timer_type == OUTPUT_BY_ITERATION_LIST)
    otd->t_iteration = world->current_iterations;
  else
    otd->t_iteration = world->current_iterations * world->time_unit;

  otd->event_time = c->data.trig.t_event * world->time_unit;
  otd->loc.x = c->data.trig.loc.x * world->length_unit;
  otd->loc.y = c->data.trig.loc.y * world->length_unit;
  otd->loc.z = c->data.trig.loc.z * world->length_unit;
  if (flags & TRIG_IS_HIT) {
    otd->how_many = 1;
    if (n > 0)
      otd->orient = 1;
    else
      otd->orient = -1;
  } else {
    otd->how_many = n;
    otd->orient = c->data.trig.orient;
  }
  otd->flags = flags;
  otd->name = ear->requester->column->expr->title;

  first_column->initial_value += 1.0;
  idx = (int)first_column->initial_value;
  if (idx >= (int)first_column->set->block->trig_bufsize) {
    if (write_reaction_output(world, first_column->set))
      mcell_error("Failed to write triggered count output to file '%s'.",
                  first_column->set->outfile_name);
    first_column->initial_value = 0;
  }
}

/*************************************************************************
flush_reaction_output:
   In: nothing
   Out: 0 on success, 1 on error (memory allocation or file I/O).
        Writes all remaining trigger events in buffers to disk.
        (Do this before ending the simulation.)
*************************************************************************/
int flush_reaction_output(struct volume *world) {
  struct schedule_helper *sh;
  struct output_block *ob;
  struct output_set *os;
  int i;
  int n_errors = 0;

  for (sh = world->count_scheduler; sh != NULL; sh = sh->next_scale) {
    for (i = 0; i <= sh->buf_len; i++) {
      if (i == sh->buf_len)
        ob = (struct output_block *)sh->current;
      else
        ob = (struct output_block *)sh->circ_buf_head[i];

      for (; ob != NULL; ob = ob->next) {
        for (os = ob->data_set_head; os != NULL; os = os->next) {
          if (write_reaction_output(world, os))
            n_errors++;
        }
      }
    }
  }

  return n_errors;
}

/**************************************************************************
update_reaction_output:
  In: the output_block we want to update
  Out: 0 on success, 1 on failure.
       The counters in this block are updated, and the block is
       rescheduled for the next output time.  The counters are saved
       to an internal buffer, and written out when full.
**************************************************************************/
int update_reaction_output(struct volume *world, struct output_block *block) {
  int report_as_non_trigger = 1;
  int i = block->buf_index;
  if (block->data_set_head != NULL &&
      block->data_set_head->column_head != NULL &&
      block->data_set_head->column_head->buffer[i].data_type == COUNT_TRIG_STRUCT)
    report_as_non_trigger = 0;

  if (report_as_non_trigger) {
    switch (world->notify->reaction_output_report) {
    case NOTIFY_NONE:
      break;

    case NOTIFY_BRIEF:
      mcell_log(
          "Updating reaction output scheduled at time %.15g on iteration %lld.",
          block->t, world->current_iterations);
      break;

    case NOTIFY_FULL:
      mcell_log("Updating reaction output scheduled at time %.15g on iteration"
                " %lld.\n  Buffer fill level is at %u/%u.",
                block->t, world->current_iterations, block->buf_index,
                block->buffersize);
      break;

    default:
      UNHANDLED_CASE(world->notify->reaction_output_report);
    }
  }

  /* update all counters */

  block->t /= (1. + EPS_C);
  if (world->chkpt_seq_num == 1) {
    if (block->timer_type == OUTPUT_BY_ITERATION_LIST)
      block->time_array[i] = block->t;
    else
      block->time_array[i] = block->t * world->time_unit;
  } else {
    if (block->timer_type == OUTPUT_BY_ITERATION_LIST) {
      block->time_array[i] = block->t;
    } else if (block->timer_type == OUTPUT_BY_TIME_LIST) {
      if (block->time_now == NULL) {
        return 0;
      } else {
        block->time_array[i] = block->time_now->value;
      }
    } else {
      /* OUTPUT_BY_STEP */
      block->time_array[i] = convert_iterations_to_seconds(
          world->start_iterations, world->time_unit,
          world->simulation_start_seconds, block->t);
    }
  }

  struct output_set *set;
  struct output_column *column;
  // Each file
  for (set = block->data_set_head; set != NULL; set = set->next) 
  {
    if (report_as_non_trigger) {
      if (world->notify->reaction_output_report == NOTIFY_FULL)
        mcell_log("  Processing reaction output file '%s'.", set->outfile_name);
    }
    // Each column
    for (column = set->column_head; column != NULL; column = column->next) 
    {
      if (column->buffer[i].data_type != COUNT_TRIG_STRUCT) {
        eval_oexpr_tree(column->expr, 1);
        switch (column->buffer[i].data_type) {
        case COUNT_INT:
          column->buffer[i].val.ival = (int)column->expr->value;
          break;

        case COUNT_DBL:
          column->buffer[i].val.dval = (double)column->expr->value;
          break;

        case COUNT_UNSET:
          column->buffer[i].val.cval = 'X';
          break;

        case COUNT_TRIG_STRUCT:
        default:
          UNHANDLED_CASE(column->buffer[i].data_type);
        }
      }
    }
  }
  block->buf_index++;

  int final_chunk_flag = 0; // flag signaling an end to the scheduled
                            // reaction outputs. Takes values {0,1}.
                            // 0 - end not reached yet,
                            // 1 - end reached.
  /* Pick time of next output, if any */
  if (block->timer_type == OUTPUT_BY_STEP)
    block->t += block->step_time / world->time_unit;
  else if (block->time_now != NULL) {
    block->time_now = block->time_now->next;
    if (block->time_now == NULL)
      final_chunk_flag = 1;
    else {
      if (block->timer_type == OUTPUT_BY_ITERATION_LIST)
        block->t = block->time_now->value;
      else {
        /* OUTPUT_BY_TIME_LIST */
        if (world->chkpt_seq_num == 1) {
          block->t = block->time_now->value / world->time_unit;
        } else {
          block->t = world->start_iterations +
                     (block->time_now->value - world->simulation_start_seconds) /
                         world->time_unit;
        }
      }
    }
  } else
    final_chunk_flag = 1;

  /* Schedule next output event--even if we're at the end, since triggers may
   * not yet be written */
  double actual_t;
  if (final_chunk_flag == 1) {
    actual_t = block->t;
    block->t = FOREVER;
  } else
    actual_t = -1;
  block->t *= (1. + EPS_C);
  if (schedule_add(world->count_scheduler, block)) {
    mcell_allocfailed_nodie("Failed to add count to scheduler.");
    return 1;
  }

  if (distinguishable(actual_t, -1, EPS_C))
    block->t = actual_t; /* Fix time for output */

  if (report_as_non_trigger &&
      world->notify->reaction_output_report == NOTIFY_FULL) {
    mcell_log("  Next output for this block scheduled at time %.15g.",
              block->t);
  }

  if (block->t >= world->iterations + 1)
    final_chunk_flag = 1;

  /* write data to outfile */
  if (block->buf_index == block->buffersize || final_chunk_flag) {
    for (set = block->data_set_head; set != NULL; set = set->next) {
      if (set->column_head->buffer[i].data_type == COUNT_TRIG_STRUCT)
        continue;
      if (write_reaction_output(world, set)) {
        mcell_error_nodie("Failed to write reaction output to file '%s'.",
                          set->outfile_name);
        return 1;
      }
    }
    block->buf_index = 0;
    no_printf("Done updating reaction output\n");
  }

  if (distinguishable(actual_t, -1, EPS_C))
    block->t = FOREVER; /* Back to infinity if we're done */

  return 0;
}

/**************************************************************************
 check_reaction_output_file:
    Check that the reaction output file is writable within the policy set by
    the user.  Creates and/or truncates the file to 0 bytes, as appropriate.
    Note that for SUBSTITUTE, the truncation is done later on, during
    initialization.

 In: parse_state: parser state
     os: output set containing file details
 Out: 0 if file preparation is successful, 1 if not.  The file named will be
      created and emptied or truncated as requested.
**************************************************************************/
int check_reaction_output_file(struct output_set *os) {
  FILE *f;
  char *name;
  struct stat fs;
  int i;

  name = os->outfile_name;

  if (make_parent_dir(name)) {
    //    mdlerror_fmt(parse_state,
    //                 "Directory for %s does not exist and could not be
    // created.",
    //                 name);
    return 1;
  }

  switch (os->file_flags) {
  case FILE_OVERWRITE:
    f = fopen(name, "w");
    if (!f) {
      switch (errno) {
      case EACCES:
        //        mdlerror_fmt(parse_state, "Access to %s denied.", name);
        return 1;
      case ENOENT:
        //        mdlerror_fmt(parse_state, "Directory for %s does not exist",
        // name);
        return 1;
      case EISDIR:
        //        mdlerror_fmt(parse_state, "%s already exists and is a
        // directory", name);
        return 1;
      default:
        //        mdlerror_fmt(parse_state, "Unable to open %s for writing",
        // name);
        return 1;
      }
    }
    fclose(f);
    break;
  case FILE_SUBSTITUTE:
    f = fopen(name, "a+");
    if (!f) {
      switch (errno) {
      case EACCES:
        //        mdlerror_fmt(parse_state, "Access to %s denied.", name);
        return 1;
      case ENOENT:
        //        mdlerror_fmt(parse_state, "Directory for %s does not exist",
        // name);
        return 1;
      case EISDIR:
        //        mdlerror_fmt(parse_state, "%s already exists and is a
        // directory", name);
        return 1;
      default:
        //        mdlerror_fmt(parse_state, "Unable to open %s for writing",
        // name);
        return 1;
      }
    }
    i = fstat(fileno(f), &fs);
    if (!i && fs.st_size == 0)
      os->file_flags = FILE_OVERWRITE;
    fclose(f);
    break;
  case FILE_APPEND:
  case FILE_APPEND_HEADER:
    f = fopen(name, "a");
    if (!f) {
      switch (errno) {
      case EACCES:
        //        mdlerror_fmt(parse_state, "Access to %s denied.", name);
        return 1;
      case ENOENT:
        //        mdlerror_fmt(parse_state, "Directory for %s does not exist",
        // name);
        return 1;
      case EISDIR:
        //        mdlerror_fmt(parse_state, "%s already exists and is a
        // directory", name);
        return 1;
      default:
        //        mdlerror_fmt(parse_state, "Unable to open %s for writing",
        // name);
        return 1;
      }
    }
    i = fstat(fileno(f), &fs);
    if (!i && fs.st_size == 0)
      os->file_flags = FILE_APPEND_HEADER;
    fclose(f);
    break;
  case FILE_CREATE:
    i = access(name, F_OK);
    if (!i) {
      i = stat(name, &fs);
      if (!i && fs.st_size > 0) {
        //        mdlerror_fmt(parse_state,
        //                     "Cannot create new file %s: it already exists",
        // name);
        return 1;
      }
    }
    f = fopen(name, "w");
    if (f == NULL) {
      switch (errno) {
      case EEXIST:
        //        mdlerror_fmt(parse_state, "Cannot create %s because it already
        // exists",
        //                     name);
        return 1;
      case EACCES:
        //        mdlerror_fmt(parse_state, "Access to %s denied.", name);
        return 1;
      case ENOENT:
        //        mdlerror_fmt(parse_state, "Directory for %s does not exist",
        // name);
        return 1;
      case EISDIR:
        //        mdlerror_fmt(parse_state, "%s already exists and is a
        // directory", name);
        return 1;
      default:
        //        mdlerror_fmt(parse_state, "Unable to open %s for writing",
        // name);
        return 1;
      }
    }
    fclose(f);
    break;

  default:
    UNHANDLED_CASE(os->file_flags);
  }
  return 0;
}

/**************************************************************************
write_reaction_output:
  In: the output_set we want to write to disk
      the flag that signals an end to the scheduled reaction outputs
  Out: 0 on success, 1 on failure.
       The reaction output buffer is flushed and written to disk.
       Indices are not reset; that's the job of the calling function.
**************************************************************************/
int write_reaction_output(struct volume *world, struct output_set *set) {

  FILE *fp;
  struct output_column *column;
  char *mode;
  u_int n_output;
  u_int i;

  switch (set->file_flags) {
  case FILE_OVERWRITE:
  case FILE_CREATE:
    if (set->chunk_count == 0)
      mode = "w";
    else
      mode = "a";
    break;
  case FILE_SUBSTITUTE:
    if (world->chkpt_seq_num == 1 && set->chunk_count == 0)
      mode = "w";
    else
      mode = "a";
    break;
  case FILE_APPEND:
  case FILE_APPEND_HEADER:
    mode = "a";
    break;
  default:
    mcell_internal_error(
        "Bad file output code %d for reaction data output file '%s'.",
        set->file_flags, set->outfile_name);
  }

  fp = open_file(set->outfile_name, mode);
  if (fp == NULL)
    return 1;

  /*int idx = set->block->buf_index;*/
  if (set->column_head->buffer[0].data_type != COUNT_TRIG_STRUCT) {
    n_output = set->block->buffersize;
    if (set->block->buf_index < set->block->buffersize)
      n_output = set->block->buf_index;

    if (world->notify->file_writes == NOTIFY_FULL)
      mcell_log("Writing %d lines to output file %s.", n_output,
                set->outfile_name);

    /* Write headers */
    if (set->chunk_count == 0 && set->header_comment != NULL &&
        set->file_flags != FILE_APPEND &&
        (world->chkpt_seq_num == 1 || set->file_flags == FILE_APPEND_HEADER ||
         set->file_flags == FILE_CREATE || set->file_flags == FILE_OVERWRITE)) {
      if (set->block->timer_type == OUTPUT_BY_ITERATION_LIST)
        fprintf(fp, "%sIteration_#", set->header_comment);
      else
        fprintf(fp, "%sSeconds", set->header_comment);

      for (column = set->column_head; column != NULL; column = column->next) {
        if (column->expr->title == NULL)
          fprintf(fp, " untitled");
        else
          fprintf(fp, " %s", column->expr->title);
      }
      fprintf(fp, "\n");
    }

    /* Write data */
    for (i = 0; i < n_output; i++) {
      fprintf(fp, "%.15g", set->block->time_array[i]);

      for (column = set->column_head; column != NULL; column = column->next) {
        switch (column->buffer[i].data_type) {
        case COUNT_INT:
          fprintf(fp, " %d", (column->buffer[i].val.ival));
          break;

        case COUNT_DBL:
          fprintf(fp, " %.9g", (column->buffer[i].val.dval));
          break;

        case COUNT_UNSET:
          fprintf(fp, " X");
          break;

        case COUNT_TRIG_STRUCT:
        default:
          if (column->expr->title != NULL)
            mcell_warn(
                "Unexpected data type in column titled '%s' -- skipping.",
                column->expr->title);
          else
            mcell_warn("Unexpected data type in untitled column -- skipping.");
          break;
        }
      }
      fprintf(fp, "\n");
    }
  } else /* Write accumulated trigger data */
  {
    struct output_trigger_data *trig;
    char event_time_string[1024]; /* Wouldn't run out of space even if we
                                     printed out DBL_MAX in non-exponential
                                     notation! */

    n_output = (u_int)set->column_head->initial_value;
    for (i = 0; i < n_output; i++) {
      trig = set->column_head->buffer[i].val.tval;

      if (set->exact_time_flag)
        sprintf(event_time_string, "%.12g ", trig->event_time);
      else
        strcpy(event_time_string, "");

      if (trig->flags & TRIG_IS_RXN) /* Just need time, pos, name */
      {
        fprintf(fp, "%.15g %s%.9g %.9g %.9g %s\n", trig->t_iteration,
                event_time_string, trig->loc.x, trig->loc.y, trig->loc.z,
                (trig->name == NULL) ? "" : trig->name);
      } else if (trig->flags & TRIG_IS_HIT) /* Need orientation also */
      {
        fprintf(fp, "%.15g %s%.9g %.9g %.9g %d %s\n", trig->t_iteration,
                event_time_string, trig->loc.x, trig->loc.y, trig->loc.z,
                trig->orient, (trig->name == NULL) ? "" : trig->name);
      } else /* Molecule count -- need both number and orientation */
      {
        fprintf(fp, "%.15g %s%.9g %.9g %.9g %d %d %s\n", trig->t_iteration,
                event_time_string, trig->loc.x, trig->loc.y, trig->loc.z,
                trig->orient, trig->how_many,
                (trig->name == NULL) ? "" : trig->name);
      }
    }
  }

  set->chunk_count++;

  fclose(fp);
  return 0;
}

/*************************************************************************
new_output_expr:
   In: mem_helper used to allocate output_expressions
   Out: New, initialized output_expression, or NULL if out of memory.
*************************************************************************/
struct output_expression *new_output_expr(struct mem_helper *oexpr_mem) {
  struct output_expression *oe;

  oe = (struct output_expression *)mem_get(oexpr_mem);
  if (oe == NULL)
    return NULL;

  oe->column = NULL;
  oe->expr_flags = 0;
  oe->up = NULL;
  oe->left = NULL;
  oe->right = NULL;
  oe->oper = '\0';
  oe->value = 0;
  oe->title = NULL;

  return oe;
}

/*************************************************************************
set_oexpr_column:
   In: output_expression who needs to have its column set (recursively)
       output_column that owns this expression
   Out: No return value.  Every output_expression in the tree gets its
        column set.
*************************************************************************/
void set_oexpr_column(struct output_expression *oe, struct output_column *oc) {
  for (; oe != NULL;
       oe = ((oe->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_OEXPR)
                ? (struct output_expression *)oe->right
                : NULL) {
    oe->column = oc;
    if ((oe->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_OEXPR)
      set_oexpr_column((struct output_expression *)oe->left, oc);
  }
}

/*************************************************************************
learn_oexpr_oexpr:
   In: output_expression whose flags may not reflect its children properly
   Out: No return value.  Flags are updated for output_expression passed
        in.  (Not recursive.)
*************************************************************************/
void learn_oexpr_flags(struct output_expression *oe) {
  struct output_expression *oel, *oer;

  oel = (struct output_expression *)oe->left;
  oer = (struct output_expression *)oe->right;

  if (oer == NULL) {
    if (oel == NULL)
      oe->expr_flags = OEXPR_TYPE_CONST | OEXPR_TYPE_DBL;
    else {
      oe->expr_flags =
          (oel->expr_flags & (OEXPR_TYPE_MASK | OEXPR_TYPE_CONST)) |
          OEXPR_LEFT_OEXPR;
      if (oel->expr_flags & OEXPR_TYPE_CONST)
        oe->expr_flags |= OEXPR_LEFT_CONST;
    }
  } else {
    oer = (struct output_expression *)oe->right;
    oe->expr_flags = OEXPR_LEFT_OEXPR | OEXPR_RIGHT_OEXPR;
    if (oel->expr_flags & OEXPR_TYPE_CONST)
      oe->expr_flags |= OEXPR_LEFT_CONST;
    if (oer->expr_flags & OEXPR_TYPE_CONST)
      oe->expr_flags |= OEXPR_RIGHT_CONST;
    if (oel->expr_flags & oer->expr_flags & OEXPR_TYPE_CONST)
      oe->expr_flags |= OEXPR_TYPE_CONST;
    if ((oel->expr_flags & OEXPR_TYPE_MASK) ==
        (oer->expr_flags & OEXPR_TYPE_MASK))
      oe->expr_flags |= oel->expr_flags & OEXPR_TYPE_MASK;
    else {
      if ((oel->expr_flags & OEXPR_TYPE_MASK) == OEXPR_TYPE_TRIG)
        oe->expr_flags |= OEXPR_TYPE_TRIG;
      else if ((oer->expr_flags & OEXPR_TYPE_MASK) == OEXPR_TYPE_TRIG)
        oe->expr_flags |= OEXPR_TYPE_TRIG;
      else if ((oel->expr_flags & OEXPR_TYPE_MASK) == OEXPR_TYPE_DBL)
        oe->expr_flags |= OEXPR_TYPE_DBL;
      else if ((oer->expr_flags & OEXPR_TYPE_MASK) == OEXPR_TYPE_DBL)
        oe->expr_flags |= OEXPR_TYPE_DBL;
      else
        oe->expr_flags |= OEXPR_TYPE_INT;
    }
  }
}

/*************************************************************************
first_oexpr_tree:
   In: expression tree
   Out: leftmost stem in that tree (joined by ',' operator)
*************************************************************************/
struct output_expression *first_oexpr_tree(struct output_expression *root) {
  while (root->oper == ',')
    root = (struct output_expression *)root->left;
  return root;
}

/*************************************************************************
next_oexpr_tree:
   In: stem in an expression tree that is joined by ',' operator
   Out: next stem to the right joined by ',' operator, or NULL if the
        current stem is the rightmost
*************************************************************************/
struct output_expression *next_oexpr_tree(struct output_expression *leaf) {
  for (; leaf->up != NULL; leaf = leaf->up) {
    if (leaf->up->left == leaf)
      return first_oexpr_tree((struct output_expression *)leaf->up->right);
  }
  return NULL;
}


/*************************************************************************
dupl_oexpr_tree:
   In: output_expression tree
       place to allocate new output_expressions
   Out: a copy of the tree, or NULL if there was a memory error.
   Note: leaves are not copied, just the expression structure.
*************************************************************************/
struct output_expression *dupl_oexpr_tree(struct output_expression *root,
                                          struct mem_helper *oexpr_mem) {
  struct output_expression *sprout;

  sprout = (struct output_expression *)mem_get(oexpr_mem);
  if (sprout == NULL)
    return NULL;

  memcpy(sprout, root, sizeof(struct output_expression));
  if (root->left != NULL &&
      (root->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_OEXPR) {
    sprout->left = dupl_oexpr_tree(root->left, oexpr_mem);
    if (sprout->left == NULL)
      return NULL;
  }
  if (root->right != NULL &&
      (root->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_OEXPR) {
    sprout->right = dupl_oexpr_tree(root->right, oexpr_mem);
    if (sprout->right == NULL)
      return NULL;
  }

  return sprout;
}

/*************************************************************************
eval_oexpr_tree:
   In: root of an output_expression tree
       flag indicating whether to recalculate values marked CONST
   Out: no return value.  The value member variable of each
        output_expression in the tree is updated to be accurate given
        current leaf values.
*************************************************************************/
void eval_oexpr_tree(struct output_expression *root, int skip_const) {
  double lval = 0.0;
  double rval = 0.0;

  if ((root->expr_flags & OEXPR_TYPE_CONST) && skip_const)
    return;
  if (root->left != NULL) {
    if ((root->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_INT)
      lval = (double)*((int *)root->left);
    else if ((root->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_DBL)
      lval = *((double *)root->left);
    else if ((root->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_OEXPR) {
      eval_oexpr_tree((struct output_expression *)root->left, skip_const);
      lval = ((struct output_expression *)root->left)->value;
    }
  }
  if (root->right != NULL) {
    if ((root->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_INT)
      rval = (double)*((int *)root->right);
    else if ((root->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_DBL)
      rval = *((double *)root->right);
    else if ((root->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_OEXPR) {
      eval_oexpr_tree((struct output_expression *)root->right, skip_const);
      rval = ((struct output_expression *)root->right)->value;
    }
  }
  switch (root->oper) {
  case '=':
    break;
  case '(':
  case '#':
  case '@':
    if (root->right != NULL)
      root->value = lval + rval;
    else
      root->value = lval;
    break;
  case '_':
    root->value = -lval;
    break;
  case '+':
    root->value = lval + rval;
    break;
  case '-':
    root->value = lval - rval;
    break;
  case '*':
    root->value = lval * rval;
    break;
  case '/':
    root->value = (!distinguishable(rval, 0, EPS_C)) ? 0 : lval / rval;
    break;
  default:
    break;
  }
}

/*************************************************************************
oexpr_flood_convert
   In: root of an expression tree
       operator that we don't want any more
       operator that we want to replace it
   Out: no return value.  The old operator is replaced with the new
        one recursively, but only as deep as the old operator remains
        unbroken (so it's like a flood fill starting with the root).
*************************************************************************/
void oexpr_flood_convert(struct output_expression *root, char old_oper,
                         char new_oper) {
  for (; root != NULL;
       root = ((root->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_OEXPR)
                  ? (struct output_expression *)root->right
                  : NULL) {
    if (root->oper != old_oper)
      return;
    root->oper = new_oper;
    if ((root->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_OEXPR)
      oexpr_flood_convert((struct output_expression *)root->left, old_oper,
                          new_oper);
  }
}

/*************************************************************************
oexpr_title:
   In: root of an expression tree
   Out: text string that describes what is in that tree, or NULL if
        there is a memory error.
   Note: the value is recursively generated, but title member variables
         of the expression are not set.  The function keeps track of
         them on the fly.
*************************************************************************/
char *oexpr_title(struct output_expression *root) {
  struct output_request *orq;
  char *lstr, *rstr, *str;

  lstr = rstr = NULL;

  if (root->expr_flags & OEXPR_TYPE_CONST) {
    if ((root->expr_flags & OEXPR_TYPE_MASK) == OEXPR_LEFT_INT)
      return alloc_sprintf("%d", (int)root->value);
    else if ((root->expr_flags & OEXPR_TYPE_MASK) == OEXPR_TYPE_DBL)
      return alloc_sprintf("%.8g", root->value);
    else
      return NULL;
  }

  if (root->left != NULL) {
    if ((root->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_INT)
      lstr = alloc_sprintf("%d", *((int *)root->left));
    else if ((root->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_DBL)
      lstr = alloc_sprintf("%.8g", *((double *)root->left));
    else if ((root->expr_flags & OEXPR_LEFT_MASK) == OEXPR_LEFT_OEXPR)
      lstr = oexpr_title((struct output_expression *)root->left);
  }
  if (root->right != NULL) {
    if ((root->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_INT)
      rstr = alloc_sprintf("%d", *((int *)root->right));
    else if ((root->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_DBL)
      rstr = alloc_sprintf("%.8g", *((double *)root->right));
    else if ((root->expr_flags & OEXPR_RIGHT_MASK) == OEXPR_RIGHT_OEXPR)
      rstr = oexpr_title((struct output_expression *)root->right);
  }

  switch (root->oper) {
  case '=':
    free(rstr);
    return lstr;

  case '@':
    free(lstr);
    free(rstr);
    return CHECKED_STRDUP("(complex)", NULL);

  case '#':
    if ((root->expr_flags & OEXPR_LEFT_MASK) != OEXPR_LEFT_REQUEST ||
        root->left == NULL) {
      free(lstr);
      free(rstr);
      return NULL;
    }
    orq = (struct output_request *)root->left;
    free(lstr);
    free(rstr);
    return strdup(orq->count_target->name);

  case '_':
    if (lstr == NULL) {
      free(rstr);
      return NULL;
    }
    str = alloc_sprintf("-%s", lstr);
    free(lstr);
    free(rstr);
    return str;

  case '(':
    if (lstr == NULL) {
      free(rstr);
      return NULL;
    }
    str = alloc_sprintf("(%s)", lstr);
    free(lstr);
    free(rstr);
    return str;

  case '+':
  case '-':
  case '*':
  case '/':
    if (lstr == NULL || rstr == NULL) {
      free(lstr);
      free(rstr);
      return NULL;
    }
    str = alloc_sprintf("%s%c%s", lstr, root->oper, rstr);
    free(lstr);
    free(rstr);
    return str;

  default:
    free(lstr);
    free(rstr);
    return NULL;
  }
  return NULL;
}
