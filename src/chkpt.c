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

/**************************************************************************\
** File: chkpt.c
**
** Purpose: Writes and reads MCell checkpoint files.
**
*/

#include "config.h"

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/stat.h>
#include <string.h>

#include "mcell_structs.h"
#include "logging.h"
#include "vol_util.h"
#include "chkpt.h"
#include "grid_util.h"
#include "count_util.h"
#include "react.h"
#include "strfunc.h"

/* MCell checkpoint API version */
#define CHECKPOINT_API 1

/* Endian-ness markers */
#define MCELL_BIG_ENDIAN 16
#define MCELL_LITTLE_ENDIAN 17

/* Checkpoint section commands */
#define CURRENT_TIME_CMD 1
#define CURRENT_ITERATION_CMD 2
#define CHKPT_SEQ_NUM_CMD 3
#define RNG_STATE_CMD 4
#define MCELL_VERSION_CMD 5
#define SPECIES_TABLE_CMD 6
#define MOL_SCHEDULER_STATE_CMD 7
#define BYTE_ORDER_CMD 8
#define NUM_CHKPT_CMDS 9
#define CHECKPOINT_API_CMD 10

/* Newbie flags */
#define HAS_ACT_NEWBIE 1
#define HAS_NOT_ACT_NEWBIE 0

#define HAS_ACT_CHANGE 1
#define HAS_NOT_ACT_CHANGE 0

/* these are needed for the chkpt signal handler */
int *chkpt_continue_after_checkpoint;
char **chkpt_initialization_state;
enum checkpoint_request_type_t *chkpt_checkpoint_requested;

/* ============================= */
/* General error-checking macros */

/* If 'op' is true, fail with a message about corrupt checkpoint data.  msg,
 * ... are printf-style. */
#define DATACHECK(op, msg, ...)                                                \
  do {                                                                         \
    if (op) {                                                                  \
      mcell_warn("Corrupted checkpoint data: " msg, ##__VA_ARGS__);            \
      return 1;                                                                \
    }                                                                          \
  } while (0)

/* If 'op' is true, fail with a message about an internal error.  msg, ... are
 * printf-style. */
#define INTERNALCHECK(op, msg, ...)                                            \
  do {                                                                         \
    if (op) {                                                                  \
      mcell_warn("%s internal: " msg, __func__, ##__VA_ARGS__);                \
      return 1;                                                                \
    }                                                                          \
  } while (0)

/* Wrapper for an I/O (write) operation to print out appropriate error messages
 * in case of failure. */
#define WRITECHECK(op, sect)                                                   \
  do {                                                                         \
    if (op) {                                                                  \
      mcell_perror_nodie(errno, "Error while writing '%s' to checkpoint file", \
                         SECTNAME);                                            \
      return 1;                                                                \
    }                                                                          \
  } while (0)

/* Write a raw field to the output stream. */
#define WRITEFIELD(f) WRITECHECK(fwrite(&(f), sizeof(f), 1, fs) != 1, SECTNAME)

/* Write a raw array of fields to the output stream. */
#define WRITEARRAY(f, len)                                                     \
  WRITECHECK(fwrite(f, sizeof(f[0]), (len), fs) != (len), SECTNAME)

/* Write an unsigned integer to the output stream in endian- and
 * size-independent format. */
#define WRITEUINT(f) WRITECHECK(write_varint(fs, (f)), SECTNAME)

/* Write a signed integer to the output stream in endian- and size-independent
 * format. */
#define WRITEINT(f) WRITECHECK(write_svarint(fs, (f)), SECTNAME)

/* Write an unsigned 64-bit niteger to the output stream in endian- and
 * size-independent format. */
#define WRITEUINT64(f) WRITECHECK(write_varintl(fs, (f)), SECTNAME)

/* Write a string to the output stream as a count followed by raw ASCII
 * character data. */
#define WRITESTRING(f)                                                         \
  do {                                                                         \
    uint32_t len = (uint32_t)strlen(f);                                        \
    WRITEUINT(len);                                                            \
    WRITEARRAY(f, len);                                                        \
  } while (0)

/* Wrapper for an I/O (read) operation to print out appropriate error messages
 * in case of failure. */
#define READCHECK(op, sect)                                                    \
  do {                                                                         \
    if (op) {                                                                  \
      mcell_perror_nodie(                                                      \
          errno, "Error while reading '%s' from checkpoint file", SECTNAME);   \
      return 1;                                                                \
    }                                                                          \
  } while (0)

/* Read a raw field from the input stream (i.e. never byteswap). */
#define READFIELDRAW(f) READCHECK(fread(&(f), sizeof(f), 1, fs) != 1, SECTNAME)

/* Maybe byteswap a field. */
#define READBSWAP(f)                                                           \
  do {                                                                         \
    if (state->byte_order_mismatch)                                            \
      byte_swap(&(f), sizeof(f));                                              \
  } while (0)

/* Read a field from the input stream, byteswapping if needed. */
#define READFIELD(f)                                                           \
  do {                                                                         \
    READFIELDRAW(f);                                                           \
    READBSWAP(f);                                                              \
  } while (0)

/* Read a raw array of chars from the input stream. */
#define READSTRING(f, len)                                                     \
  do {                                                                         \
    READCHECK(fread(f, sizeof(f[0]), (len), fs) != (len), SECTNAME);           \
    f[len] = '\0';                                                             \
  } while (0)

/* Read an array of items from the input stream, byteswapping if needed. */
#define READARRAY(f, len)                                                      \
  do {                                                                         \
    int i;                                                                     \
    READCHECK(fread(f, sizeof(f[0]), (len), fs) != (len), SECTNAME);           \
    for (i = 0; i < len; ++i)                                                  \
      READBSWAP(f[i]);                                                         \
  } while (0)

/* Read an unsigned integer from the input stream in endian- and
 * size-independent format. */
#define READUINT(f) READCHECK(read_varint(fs, &(f)), SECTNAME)

/* Read a signed integer from the input stream in endian- and size-independent
 * format. */
#define READINT(f) READCHECK(read_svarint(fs, &(f)), SECTNAME)

/* Read an unsigned 64-bit niteger from the input stream in endian- and
 * size-independent format. */
#define READUINT64(f) READCHECK(read_varintl(fs, &(f)), SECTNAME)

/**
 * Structure to hold any requisite state during the checkpoint read process.
 */
struct chkpt_read_state {
  byte byte_order_mismatch;
};

/* Handlers for individual checkpoint commands */
static int read_current_time_seconds(struct volume *world, FILE *fs,
                                     struct chkpt_read_state *state);
static int read_current_iteration(struct volume *world, FILE *fs,
                                  struct chkpt_read_state *state);
static int read_chkpt_seq_num(struct volume *world, FILE *fs,
                              struct chkpt_read_state *state);
static int read_rng_state(struct volume *world, FILE *fs,
                          struct chkpt_read_state *state);
static int read_byte_order(FILE *fs, struct chkpt_read_state *state);
static int read_mcell_version(FILE *fs, struct chkpt_read_state *state);
static int read_api_version(FILE *fs, struct chkpt_read_state *state,
  uint32_t *api_version);
static int read_species_table(struct volume *world, FILE *fs);
static int read_mol_scheduler_state_real(struct volume *world, FILE *fs,
                                         struct chkpt_read_state *state,
                                         uint32_t api_version);
static int write_mcell_version(FILE *fs, const char *mcell_version);
static int write_current_time_seconds(FILE *fs, double current_time_seconds);
static int write_current_iteration(FILE *fs, long long current_iterations,
                                   double current_time_seconds);
static int write_chkpt_seq_num(FILE *fs, u_int chkpt_seq_num);
static int write_rng_state(FILE *fs, u_int seed_seq, struct rng_state *rng);
static int write_species_table(FILE *fs, int n_species,
                               struct species **species_list);
static int write_mol_scheduler_state_real(FILE *fs,
                                          struct storage_list *storage_head,
                                          double simulation_start_seconds,
                                          double start_iterations,
                                          double time_unit);
static int write_byte_order(FILE *fs);

static int write_api_version(FILE *fs);

static int create_molecule_scheduler(struct storage_list *storage_head,
                                     long long start_iterations);

/********************************************************************
 * this function initializes to global variables
 *
 *     int chkpt_continue_after_checkpoint
 *     char *chkpt_initialization_state
 *     enum checkpoint_request_type_t chkpt_checkpoint_requested
 *
 * needed by the signal handler chkpt_signal_handler
*********************************************************************/
int set_checkpoint_state(struct volume *world) {
  chkpt_continue_after_checkpoint = &world->continue_after_checkpoint;
  chkpt_initialization_state = &world->initialization_state;
  chkpt_checkpoint_requested = &world->checkpoint_requested;

  return 0;
}

/***************************************************************************
 chkpt_signal_handler:
 In:  signo - the signal number that triggered the checkpoint
 Out: Records the checkpoint request

 Note: This function is not to be called during normal program execution.  It is
 registered as a signal handler for SIGUSR1, SIGUSR2, and possibly SIGALRM
 signals.
***************************************************************************/
void chkpt_signal_handler(int signo) {
  if (*chkpt_initialization_state) {
    if (signo != SIGALRM || !*chkpt_continue_after_checkpoint) {
      mcell_warn("Checkpoint requested while %s.  Exiting.",
                 *chkpt_initialization_state);
      exit(EXIT_FAILURE);
    }
  }

#ifndef _WIN32 /* fixme: Windows does not support USR signals */
  if (signo == SIGUSR1)
    *chkpt_checkpoint_requested = CHKPT_SIGNAL_CONT;
  else if (signo == SIGUSR2)
    *chkpt_checkpoint_requested = CHKPT_SIGNAL_EXIT;
  else
#endif
      if (signo == SIGALRM) {
    if (*chkpt_continue_after_checkpoint)
      *chkpt_checkpoint_requested = CHKPT_ALARM_CONT;
    else
      *chkpt_checkpoint_requested = CHKPT_ALARM_EXIT;
  }
}

/***************************************************************************
 create_chkpt:
 In:  filename - the name of the checkpoint file to create
 Out: returns 1 on failure, 0 on success.  On success, checkpoint file is
      written to the appropriate filename.  On failure, the old checkpoint file
      is left unmolested.
***************************************************************************/
int create_chkpt(struct volume *world, char const *filename) {
  FILE *outfs = NULL;

  /* Create temporary filename */
  char *tmpname = alloc_sprintf("%s.tmp", filename);
  if (tmpname == NULL)
    mcell_allocfailed("Out of memory creating temporary checkpoint filename "
                      "for checkpoint '%s'.",
                      filename);

  /* Open the file */
  if ((outfs = fopen(tmpname, "wb")) == NULL)
    mcell_perror(errno, "Failed to write checkpoint file '%s'", tmpname);

  /* Write checkpoint */
  world->current_time_seconds = world->current_time_seconds +
      (world->current_iterations - world->start_iterations) * world->time_unit;
  // These are normally set when reading a checkpoint. They need to be set here
  // in case we checkpoint without exiting (i.e. using NOEXIT). Otherwise,
  // world->current_time_seconds will be set incorrectly upon subsequent calls
  // to create_chkpt
  world->start_iterations = world->current_iterations;
  world->simulation_start_seconds = world->current_time_seconds;
  if (write_chkpt(world, outfs))
    mcell_error("Failed to write checkpoint file %s\n", filename);
  fclose(outfs);

  /* keep previous checkpoint file if requested by appending the current
   * iteration */
  if (world->keep_chkpts) {
    /* check if previous checkpoint file exists - may not exist initially */
    struct stat buf;
    if (stat(filename, &buf) == 0) {
      char *keepName = alloc_sprintf("%s.%lld", filename, world->current_iterations);
      if (keepName == NULL) {
        mcell_allocfailed("Out of memory creating filename for checkpoint");
      }

      if (rename(filename, keepName) != 0) {
        mcell_error("Failed to save previous checkpoint file %s to %s",
                    filename, keepName);
      }
      free(keepName);
    }
  }

  /* Move it into place */
  if (rename(tmpname, filename) != 0)
    mcell_error("Successfully wrote checkpoint to file '%s', but failed to "
                "atomically replace checkpoint file '%s'.\nThe simulation may "
                "be resumed from '%s'.",
                tmpname, filename, tmpname);

  free(tmpname);
  return 0;
}

/***************************************************************************
 write_varintl: Size- and endian-agnostic saving of unsigned long long values.
 In:  fs - file handle to which to write
      val - value to write to file
 Out: returns 1 on failure, 0 on success.  On success, value is written to file.
***************************************************************************/
static int write_varintl(FILE *fs, unsigned long long val) {
  unsigned char buffer[40];
  size_t len = 0;

  buffer[sizeof(buffer) - 1 - len] = val & 0x7f;
  val >>= 7;
  ++len;

  while (val != 0) {
    buffer[sizeof(buffer) - 1 - len] = (val & 0x7f) | 0x80;
    val >>= 7;
    ++len;
  }

  if (fwrite(buffer + sizeof(buffer) - len, 1, len, fs) != len)
    return 1;
  return 0;
}

/***************************************************************************
 write_svarintl: Size- and endian-agnostic saving of signed long long values.
 In:  fs - file handle to which to write
      val - value to write to file
 Out: returns 1 on failure, 0 on success.  On success, value is written to file.
***************************************************************************/
static int write_svarintl(FILE *fs, long long val) {
  if (val < 0)
    return write_varintl(fs, (unsigned long long)(((-val) << 1) | 1));
  else
    return write_varintl(fs, (unsigned long long)(((val) << 1)));
}

/***************************************************************************
 read_varintl: Size- and endian-agnostic loading of unsigned long long values.
 In:  fs - file handle from which to read
      dest - pointer to value to receive data from file
 Out: returns 1 on failure, 0 on success.  On success, value is read from file.
***************************************************************************/
static int read_varintl(FILE *fs, unsigned long long *dest) {
  unsigned long long accum = 0;
  unsigned char ch;
  do {
    if (fread(&ch, 1, 1, fs) != 1)
      return 1;
    accum <<= 7;
    accum |= ch & 0x7f;
  } while (ch & 0x80);

  *dest = accum;
  return 0;
}

/***************************************************************************
 read_svarintl: Size- and endian-agnostic loading of signed long long values.
 In:  fs - file handle from which to read
      dest - pointer to value to receive data from file
 Out: returns 1 on failure, 0 on success.  On success, value is read from file.
***************************************************************************/
static int read_svarintl(FILE *fs, long long *dest) {
  unsigned long long tmp = 0;
  if (read_varintl(fs, &tmp))
    return 1;

  if (tmp & 1)
    *dest = (long long)-(tmp >> 1);
  else
    *dest = (long long)(tmp >> 1);
  return 0;
}

/***************************************************************************
 write_varint: Size- and endian-agnostic saving of unsigned int values.
 In:  fs - file handle to which to write
      val - value to write to file
 Out: returns 1 on failure, 0 on success.  On success, value is written to file.
***************************************************************************/
static int write_varint(FILE *fs, unsigned int val) {
  return write_varintl(fs, (unsigned long long)val);
}

/***************************************************************************
 write_svarint: Size- and endian-agnostic saving of signed int values.
 In:  fs - file handle to which to write
      val - value to write to file
 Out: returns 1 on failure, 0 on success.  On success, value is written to file.
***************************************************************************/
static int write_svarint(FILE *fs, int val) {
  return write_svarintl(fs, (long long)val);
}

/***************************************************************************
 read_varint: Size- and endian-agnostic loading of unsigned int values.
 In:  fs - file handle from which to read
      dest - pointer to value to receive data from file
 Out: returns 1 on failure, 0 on success.  On success, value is read from file.
***************************************************************************/
static int read_varint(FILE *fs, unsigned int *dest) {
  unsigned long long val;
  if (read_varintl(fs, &val))
    return 1;

  *dest = val;
  if ((unsigned long long)*dest != val)
    return 1;

  return 0;
}

/***************************************************************************
 read_svarint: Size- and endian-agnostic loading of signed int values.
 In:  fs - file handle from which to read
      dest - pointer to value to receive data from file
 Out: returns 1 on failure, 0 on success.  On success, value is read from file.
***************************************************************************/
static int read_svarint(FILE *fs, int *dest) {
  long long val;
  if (read_svarintl(fs, &val))
    return 1;

  *dest = val;
  if ((long long)*dest != val)
    return 1;

  return 0;
}

/***************************************************************************
 write_chkpt:
 In:  fs - checkpoint file to write to.
 Out: Writes the checkpoint file with all information needed for the
      simulation to restart.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
int write_chkpt(struct volume *world, FILE *fs) {
  return (write_byte_order(fs) ||
          write_api_version(fs) ||
          write_mcell_version(fs, world->mcell_version) ||
          write_current_time_seconds(fs, world->current_time_seconds) ||
          write_current_iteration(fs, world->current_iterations,
                                  world->current_time_seconds) ||
          write_chkpt_seq_num(fs, world->chkpt_seq_num) ||
          write_rng_state(fs, world->seed_seq, world->rng) ||
          write_species_table(fs, world->n_species, world->species_list) ||
          write_mol_scheduler_state_real(fs, world->storage_head,
              world->simulation_start_seconds, world->start_iterations,
              world->time_unit));
}

/***************************************************************************
 read_preamble:
    Read the required first two sections of an MCell checkpoint file.

 In:  fs:  checkpoint file to read from.
 Out: Reads preamble from checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_preamble(FILE *fs, struct chkpt_read_state *state, uint32_t *api_version) {
  byte cmd;

  /* Read and handle required first command (byte order). */
  size_t count = fread(&cmd, 1, sizeof(cmd), fs);
  DATACHECK(feof(fs), "Checkpoint file is empty.");
  if (count != sizeof(cmd)) {
    return 1; 
  }
  DATACHECK(cmd != BYTE_ORDER_CMD,
            "Checkpoint file does not have the required byte-order command.");
  if (read_byte_order(fs, state))
    return 1;

  /* Read the checkpoint API version if present. If it is missing we assume
   * the lagacy API version 0 */
  count = fread(&cmd, 1, sizeof(cmd), fs);
  DATACHECK(feof(fs), "Checkpoint file is too short (no api or version info).");
  if (count != sizeof(cmd)) {
    return 1; 
  }
  if (cmd != CHECKPOINT_API_CMD) {
    *api_version = 0;
  } else {
    read_api_version(fs, state, api_version);
    count = fread(&cmd, 1, sizeof(cmd), fs);
  }

  /* Read and handle required second command (version). */
  DATACHECK(feof(fs), "Checkpoint file is too short (no version info).");
  if (count != sizeof(cmd)) {
    return 1; 
  }
  DATACHECK(cmd != MCELL_VERSION_CMD,
            "Checkpoint file does not contain required MCell version command.");
  return read_mcell_version(fs, state);
}

/***************************************************************************
 read_chkpt:
 In:  fs - checkpoint file to read from.
 Out: Reads checkpoint file.  Sets the values of multiple parameters
      in the simulation.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
int read_chkpt(struct volume *world, FILE *fs) {
  byte cmd;

  int seen_section[NUM_CHKPT_CMDS];
  memset(seen_section, 0, sizeof(int)*NUM_CHKPT_CMDS);

  struct chkpt_read_state state;
  state.byte_order_mismatch = 0;

  /* Read the required pre-amble sections */
  uint32_t api_version;
  if (read_preamble(fs, &state, &api_version))
    return 1;
  seen_section[BYTE_ORDER_CMD] = 1;
  seen_section[MCELL_VERSION_CMD] = 1;

  /* Handle all other commands */
  while (1) {
    size_t count = fread(&cmd, sizeof(cmd), 1, fs);
    if (feof(fs)) {
      break;
    }
    if (count != sizeof(cmd)) {
      return 1; 
    }

    /* Check that it's a valid command-type */
    DATACHECK(cmd < 1 || cmd >= NUM_CHKPT_CMDS,
              "Unrecognized command-type in checkpoint file.  "
              "Checkpoint file cannot be loaded.");

    /* Check that we haven't seen it already */
    DATACHECK(seen_section[cmd], "Duplicate command-type in checkpoint file.");
    seen_section[cmd] = 1;

    /* Process normal commands */
    switch (cmd) {
    case CURRENT_TIME_CMD:
      if (read_current_time_seconds(world, fs, &state))
        return 1;
      break;

    case CURRENT_ITERATION_CMD:
      if (read_current_iteration(world, fs, &state) ||
          create_molecule_scheduler(world->storage_head, world->start_iterations))
        return 1;
      break;

    case CHKPT_SEQ_NUM_CMD:
      if (read_chkpt_seq_num(world, fs, &state))
        return 1;
      break;

    case RNG_STATE_CMD:
      if (read_rng_state(world, fs, &state))
        return 1;
      break;

    case SPECIES_TABLE_CMD:
      if (read_species_table(world, fs))
        return 1;
      break;

    case MOL_SCHEDULER_STATE_CMD:
      DATACHECK(
          !seen_section[CURRENT_ITERATION_CMD],
          "Current iteration command must precede molecule scheduler command.");
      DATACHECK(
          !seen_section[SPECIES_TABLE_CMD],
          "Species table command must precede molecule scheduler command.");
      if (read_mol_scheduler_state_real(world, fs, &state, api_version))
        return 1;
      break;

    case BYTE_ORDER_CMD:
    case MCELL_VERSION_CMD:
    default:
      /* We should have already filtered out these cases, so if we get here,
       * an internal error has occurred. */
      assert(0);
      break;
    }
  }

  /* Check for required sections */
  DATACHECK(!seen_section[CURRENT_TIME_CMD],
            "Current time command is not present.");
  DATACHECK(!seen_section[CHKPT_SEQ_NUM_CMD],
            "Checkpoint sequence number command is not present.");
  DATACHECK(!seen_section[RNG_STATE_CMD], "RNG state command is not present.");
  DATACHECK(!seen_section[MOL_SCHEDULER_STATE_CMD],
            " Molecule scheduler state command is not present.");

  return 0;
}

/***************************************************************************
 write_byte_order:
 In:  fs - checkpoint file to write to.
 Out: Writes byte order of the machine that creates checkpoint file
         to the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_byte_order(FILE *fs) {
  static const char SECTNAME[] = "byte order";
  static const byte cmd = BYTE_ORDER_CMD;
#ifdef WORDS_BIGENDIAN
  static const unsigned int byte_order = MCELL_BIG_ENDIAN;
#else
  static const unsigned int byte_order = MCELL_LITTLE_ENDIAN;
#endif

  WRITEFIELD(cmd);
  WRITEFIELD(byte_order);
  return 0;
}

/***************************************************************************
 read_byte_order:
 In:  fs - checkpoint file to read from.
 Out: Reads byte order  from the checkpoint file.
      Returns 1 on error, and 0 - on success.
      Reports byte order mismatch between the machines that writes to and
        reads from the checkpoint file,
***************************************************************************/
static int read_byte_order(FILE *fs, struct chkpt_read_state *state) {
  static const char SECTNAME[] = "byte order";
#ifdef WORDS_BIGENDIAN
  static const unsigned int byte_order_present = MCELL_BIG_ENDIAN;
#else
  static const unsigned int byte_order_present = MCELL_LITTLE_ENDIAN;
#endif

  unsigned int byte_order_read;
  READFIELDRAW(byte_order_read);

  /* find whether there is mismatch between two machines */
  state->byte_order_mismatch = (byte_order_read != byte_order_present);

  return 0;
}

/***************************************************************************
 write_api_version:
 In:  fs - checkpoint file to write to.
 Out: Writes the current checkpoint file api version to the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_api_version(FILE *fs) {
  static const char SECTNAME[] = "api version";
  static const byte cmd = CHECKPOINT_API_CMD;
  uint32_t api_version = CHECKPOINT_API;

  WRITEFIELD(cmd);
  WRITEFIELD(api_version);
  return 0;
}

/***************************************************************************
 read_api_version:
 In:  fs - checkpoint file to read from.
 Out: Reads the api version of the current checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_api_version(FILE *fs, struct chkpt_read_state *state,
  uint32_t *api_version) {
  static const char SECTNAME[] = "api version";
  READFIELD(*api_version);
  return 0;
}

/***************************************************************************
 write_mcell_version:
 In:  fs - checkpoint file to write to.
 Out: MCell3 software version is written to the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_mcell_version(FILE *fs, const char *mcell_version) {
  static const char SECTNAME[] = "MCell version";
  static const byte cmd = MCELL_VERSION_CMD;

  uint32_t len = (uint32_t)strlen(mcell_version);

  WRITEFIELD(cmd);
  WRITEFIELD(len);
  WRITEARRAY(mcell_version, len);
  return 0;
}

/***************************************************************************
 read_mcell_version:
 In:  fs - checkpoint file to read from.
 Out: MCell3 software version is read from the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_mcell_version(FILE *fs, struct chkpt_read_state *state) {
  static const char SECTNAME[] = "MCell version";

  /* Read and error-check the MCell version string length */
  unsigned int version_length;
  READFIELD(version_length);
  DATACHECK(version_length >= 100000,
            "Length field for MCell version is too long (%u).", version_length);

  char mcell_version[version_length + 1];
  READSTRING(mcell_version, version_length);

  mcell_log("Checkpoint file was created with MCell Version %s.",
            mcell_version);

  return 0;
}

/***************************************************************************
 write_current_time_seconds:
 In:  fs - checkpoint file to write to.
 Out: Writes current real time (in the terms of sec) in the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_current_time_seconds(FILE *fs, double current_time_seconds) {
  static const char SECTNAME[] = "current real time";
  static const byte cmd = CURRENT_TIME_CMD;

  WRITEFIELD(cmd);
  WRITEFIELD(current_time_seconds);
  return 0;
}

/***************************************************************************
 read_current_time_seconds:
 In:  fs - checkpoint file to read from.
 Out: Reads current real time (in the terms of sec) from the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_current_time_seconds(struct volume *world, FILE *fs,
                                  struct chkpt_read_state *state) {
  static const char SECTNAME[] = "current real time";
  READFIELD(world->simulation_start_seconds);
  return 0;
}

/***************************************************************************
 create_molecule_scheduler:
 In:  none.
 Out: Creates global molecule scheduler using checkpoint file values.
      Returns 0 on success. Error message and exit on failure.
***************************************************************************/
static int create_molecule_scheduler(struct storage_list *storage_head,
                                     long long start_iterations) {
  struct storage_list *stg;
  for (stg = storage_head; stg != NULL; stg = stg->next) {
    if ((stg->store->timer = create_scheduler(1.0, 100.0, 100, start_iterations)) ==
        NULL) {
      mcell_error("Out of memory while creating molecule scheduler.");
    }
    stg->store->current_time = start_iterations;
  }

  return 0;
}

/***************************************************************************
 write_current_iteration:
 In:  fs - checkpoint file to write to.
 Out: Writes current iteration number to the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_current_iteration(FILE *fs, long long current_iterations,
                                   double current_time_seconds) {
  static const char SECTNAME[] = "current iteration";
  static const byte cmd = CURRENT_ITERATION_CMD;

  WRITEFIELD(cmd);
  WRITEFIELD(current_iterations);
  WRITEFIELD(current_time_seconds);
  return 0;
}

/***************************************************************************
 read_current_iteration:
 In:  fs - checkpoint file to read from.
 Out: Reads current iteration number from the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_current_iteration(struct volume *world, FILE *fs,
                                  struct chkpt_read_state *state) {
  static const char SECTNAME[] = "current iteration";
  READFIELD(world->start_iterations);
  READFIELD(world->chkpt_start_time_seconds);
  world->current_time_seconds = world->chkpt_start_time_seconds;
  return 0;
}

/***************************************************************************
 write_chkpt_seq_num:
 In:  fs - checkpoint file to write to.
 Out: Writes checkpoint sequence number to the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_chkpt_seq_num(FILE *fs, u_int chkpt_seq_num) {
  static const char SECTNAME[] = "checkpoint sequence number";
  static const byte cmd = CHKPT_SEQ_NUM_CMD;

  WRITEFIELD(cmd);
  WRITEFIELD(chkpt_seq_num);
  return 0;
}

/***************************************************************************
 read_chkpt_seq_num:
 In:  fs - checkpoint file to read from.
 Out: Reads checkpoint sequence number from the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_chkpt_seq_num(struct volume *world, FILE *fs,
                              struct chkpt_read_state *state) {
  static const char SECTNAME[] = "checkpoint sequence number";
  READFIELD(world->chkpt_seq_num);
  ++world->chkpt_seq_num;
  return 0;
}

/***************************************************************************
 write_an_rng_state:

  In:  fs: checkpoint file to write to.
       rng: rng whose state to write
  Out: Writes random number generator state to the checkpoint file.
       Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_an_rng_state(FILE *fs, struct rng_state *rng) {
  static const char SECTNAME[] = "RNG state";

#ifdef USE_MINIMAL_RNG
  static const char RNG_MINRNG = 'M';
  WRITEFIELD(RNG_MINRNG);
  WRITEFIELD(rng->a);
  WRITEFIELD(rng->b);
  WRITEFIELD(rng->c);
  WRITEFIELD(rng->d);
#else
  static const char RNG_ISAAC = 'I';
  WRITEFIELD(RNG_ISAAC);
  WRITEUINT(rng->randcnt);
  WRITEFIELD(rng->aa);
  WRITEFIELD(rng->bb);
  WRITEFIELD(rng->cc);
  WRITEARRAY(rng->randrsl, RANDSIZ);
  WRITEARRAY(rng->mm, RANDSIZ);
#endif
  return 0;
}

/***************************************************************************
 write_rng_state:
 In:  fs - checkpoint file to write to.
 Out: Writes random number generator state to the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_rng_state(FILE *fs, u_int seed_seq, struct rng_state *rng) {
  static const char SECTNAME[] = "RNG state";
  static const byte cmd = RNG_STATE_CMD;

  WRITEFIELD(cmd);
  WRITEUINT(seed_seq);
  if (write_an_rng_state(fs, rng))
    return 1;
  return 0;
}

/***************************************************************************
 read_an_rng_state:
 In:  fs: checkpoint file to read from.
      state: contextual state for reading checkpoint file
      rng: pointer to RNG state to read
 Out: Reads random number generator state from the checkpoint file.
     Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_an_rng_state(FILE *fs, struct chkpt_read_state *state,
                             struct rng_state *rng) {
  static const char SECTNAME[] = "RNG state";

#ifdef USE_MINIMAL_RNG
  static const char RNG_MINRNG = 'M';
  char rngtype;
  READFIELD(rngtype);
  DATACHECK(rngtype != RNG_MINRNG, "Invalid RNG type stored in checkpoint file "
                                   "(in this version of MCell, only Bob "
                                   "Jenkins' \"small PRNG\" is supported).");
  READFIELD(rng->a);
  READFIELD(rng->b);
  READFIELD(rng->c);
  READFIELD(rng->d);

#else
  static const char RNG_ISAAC = 'I';
  char rngtype;
  READFIELD(rngtype);
  DATACHECK(rngtype != RNG_ISAAC, "Invalid RNG type stored in checkpoint file "
                                  "(in this version of MCell, only ISAAC64 is "
                                  "supported).");

  READUINT(rng->randcnt);
  READFIELD(rng->aa);
  READFIELD(rng->bb);
  READFIELD(rng->cc);
  READARRAY(rng->randrsl, RANDSIZ);
  READARRAY(rng->mm, RANDSIZ);
  rng->rngblocks = 1;
#endif

  return 0;
}

/***************************************************************************
 read_rng_state:
 In:  fs - checkpoint file to read from.
 Out: Reads random number generator state from the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_rng_state(struct volume *world, FILE *fs,
                          struct chkpt_read_state *state) {
  static const char SECTNAME[] = "RNG state";

  /* Load seed_seq from chkpt file to compare with seed_seq from command line.
   * If seeds match, we will continue with the saved RNG state from the chkpt
   * file, otherwise, we'll reinitialize rng to beginning of new seed sequence.
   */
  unsigned int old_seed;
  READUINT(old_seed);

  /* Whether we're going to use it or not, we need to read the RNG state in
   * order to advance to the next cmd in the chkpt file.
   */
  if (read_an_rng_state(fs, state, world->rng))
    return 1;

  /* Reinitialize rngs to beginning of new seed sequence, if necessary. */
  if (world->seed_seq != old_seed)
    rng_init(world->rng, world->seed_seq);

  return 0;
}

/***************************************************************************
 write_species_table:
 In:  fs - checkpoint file to write to.
 Out: Writes species data to the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_species_table(FILE *fs, int n_species,
                               struct species **species_list) {
  static const char SECTNAME[] = "species table";
  static const byte cmd = SPECIES_TABLE_CMD;

  WRITEFIELD(cmd);

  /* Write total number of existing species with at least one representative. */
  unsigned int non_empty_species_count = 0;
  int i;
  for (i = 0; i < n_species; i++) {
    if (species_list[i]->population > 0)
      ++non_empty_species_count;
  }
  WRITEUINT(non_empty_species_count);

  /* Write out all species which have at least one representative. */
  unsigned int external_species_id = 0;
  for (i = 0; i < n_species; i++) {
    if (species_list[i]->population == 0)
      continue;

    /* Write species name and external id. */
    WRITESTRING(species_list[i]->sym->name);
    WRITEUINT(external_species_id);

    /* Allocate and write the external species id for this species.  Stored
     * value is used in wrote_mol_scheduler_state to assign external species
     * ids to mols as they are written out. */
    species_list[i]->chkpt_species_id = external_species_id++;
  }

  return 0;
}

/***************************************************************************
 read_species_table:
 In:  fs - checkpoint file to read from.
 Out: Reads species data from the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int read_species_table(struct volume *world, FILE *fs) {
  static const char SECTNAME[] = "species table";

  /* Read total number of species contained in checkpoint file. */
  unsigned int total_species;
  READUINT(total_species);

  /* Scan over species table, reading in species data */
  for (unsigned int i = 0; i < total_species; i++) {
    unsigned int species_name_length;
    READUINT(species_name_length);

    /* Jed, 2007-07-05: This check guards against the security flaw caused by
     * allocating a variable-length array on the stack.  We should have
     * similar checks whenever we allocate a stack buffer based upon a
     * user-supplied value rather than upon an actual count.
     */
    DATACHECK(
        species_name_length >= 100000,
        "Species table has species name greater than 100000 characters (%u).",
        species_name_length);

    /* Read species fields. */
    char species_name[species_name_length + 1];
    unsigned int external_species_id;

    READSTRING(species_name, species_name_length);
    READUINT(external_species_id);

    /* find this species among world->species */
    int j;
    for (j = 0; j < world->n_species; j++) {
      if ((strcmp(world->species_list[j]->sym->name, species_name) == 0)) {
        world->species_list[j]->chkpt_species_id = external_species_id;
        break;
      }
    }
    DATACHECK(j == world->n_species, "Checkpoint file contains data for "
                                     "species '%s', which does not exist in "
                                     "this simulation.",
              species_name);
  }

  return 0;
}

/***************************************************************************
 count_items_in_scheduler:
 In:  None
 Out: Number of non-defunct molecules in the molecule scheduler
***************************************************************************/
unsigned long long
count_items_in_scheduler(struct storage_list *storage_head) {
  unsigned long long total_items = 0;

  for (struct storage_list *slp = storage_head; slp != NULL; slp = slp->next) {
    for (struct schedule_helper *shp = slp->store->timer; shp != NULL;
         shp = shp->next_scale) {
      for (int i = -1; i < shp->buf_len; i++) {
        for (struct abstract_element *aep = (i < 0) ? shp->current
                                                    : shp->circ_buf_head[i];
             aep != NULL; aep = aep->next) {
          struct abstract_molecule *amp = (struct abstract_molecule *)aep;
          if (amp->properties == NULL)
            continue;

          /* There should never be a surface class in the scheduler... */
          assert(!(amp->properties->flags & IS_SURFACE));
          ++total_items;
        }
      }
    }
  }

  return total_items;
}

/***************************************************************************
 write_mol_scheduler_state_real:
 In:  fs - checkpoint file to write to.
 Out: Writes molecule scheduler data to the checkpoint file.
      Returns 1 on error, and 0 - on success.
***************************************************************************/
static int write_mol_scheduler_state_real(FILE *fs,
                                          struct storage_list *storage_head,
                                          double simulation_start_seconds,
                                          double start_iterations,
                                          double time_unit) {
  static const char SECTNAME[] = "molecule scheduler state";
  static const byte cmd = MOL_SCHEDULER_STATE_CMD;

  WRITEFIELD(cmd);

  /* write total number of items in the scheduler */
  unsigned long long total_items = count_items_in_scheduler(storage_head);
  WRITEUINT64(total_items);

  /* Iterate over all molecules in the scheduler to produce checkpoint */
  for (struct storage_list *slp = storage_head; slp != NULL; slp = slp->next) {
    for (struct schedule_helper *shp = slp->store->timer; shp != NULL;
         shp = shp->next_scale) {
      for (int i = -1; i < shp->buf_len; i++) {
        for (struct abstract_element *aep = (i < 0) ? shp->current
                                                    : shp->circ_buf_head[i];
             aep != NULL; aep = aep->next) {
          struct abstract_molecule *amp = (struct abstract_molecule *)aep;
          if (amp->properties == NULL)
            continue;

          /* Grab the location and orientation for this molecule */
          struct vector3 where;
          short orient = 0;
          byte act_newbie_flag =
              (amp->flags & ACT_NEWBIE) ? HAS_ACT_NEWBIE : HAS_NOT_ACT_NEWBIE;
          byte act_change_flag =
              (amp->flags & ACT_CHANGE) ? HAS_ACT_CHANGE : HAS_NOT_ACT_CHANGE;
          if ((amp->properties->flags & NOT_FREE) == 0) {
            struct volume_molecule *vmp = (struct volume_molecule *)amp;
            INTERNALCHECK(vmp->previous_wall != NULL && vmp->index >= 0,
                          "The value of 'previous_grid' is not NULL.");
            where.x = vmp->pos.x;
            where.y = vmp->pos.y;
            where.z = vmp->pos.z;
            orient = 0;
          } else if ((amp->properties->flags & ON_GRID) != 0) {
            struct surface_molecule *smp = (struct surface_molecule *)amp;
            uv2xyz(&smp->s_pos, smp->grid->surface, &where);
            orient = smp->orient;
          } else
            continue;

          /* Check for valid chkpt_species ID. */
          INTERNALCHECK(amp->properties->chkpt_species_id == UINT_MAX,
                        "Attempted to write out a molecule of species '%s', "
                        "which has not been assigned a checkpoint species id.",
                        amp->properties->sym->name);

          /* write molecule fields */
          WRITEUINT(amp->properties->chkpt_species_id);
          WRITEFIELD(act_newbie_flag);
          WRITEFIELD(act_change_flag);

          // NOTE: we write all times as real times (seconds) *not* as
          // "iterations" (or "scaled times") in order to be able to
          // re-schedule them properly upon restart

          // The scheduling time (t) is essentially iterations, and since time
          // steps can change when checkpointing, we can't directly convert
          // iterations to real time (seconds). We need to correct for this by
          // only converting the iterations of the current simulation
          // [(t-start_iterations)*time_unit] and adding the real time at the
          // start of the simulation (simulation_start_seconds).
          double t = convert_iterations_to_seconds(
              start_iterations, time_unit, simulation_start_seconds, amp->t);
          WRITEFIELD(t);
          // We do a simple conversion for the lifetime t2, since this
          // corresponds to some event in the future and can be directly
          // computed without using an offset.
          double t2 = amp->t2 * time_unit;
          WRITEFIELD(t2);
          // Birthday is now always treated as real time in seconds, not
          // "scaled" time or iterations.
          double bday = amp->birthday;
          WRITEFIELD(bday);
          WRITEFIELD(where);
          WRITEINT(orient);

          static const unsigned char NON_COMPLEX = '\0';
          WRITEFIELD(NON_COMPLEX);
        }
      }
    }
  }

  return 0;
}

/***************************************************************************
 read_mol_scheduler_state_real:
 In:  fs - checkpoint file to read from.
 Out: Reads molecule scheduler data from the checkpoint file.
      Returns 0 on success. Error message and exit on failure.
***************************************************************************/
static int read_mol_scheduler_state_real(struct volume *world, FILE *fs,
                                         struct chkpt_read_state *state,
                                         uint32_t api_version) {
  static const char SECTNAME[] = "molecule scheduler state";

  struct volume_molecule vm;
  struct volume_molecule *vmp = NULL;
  struct abstract_molecule *amp = NULL;
  struct volume_molecule *guess = NULL;

  /* Clear template vol mol structure */
  memset(&vm, 0, sizeof(struct volume_molecule));
  vmp = &vm;
  amp = (struct abstract_molecule *)vmp;

  /* read total number of items in the scheduler. */
  unsigned long long total_items;
  READUINT64(total_items);

  for (unsigned long long n_mol = 0; n_mol < total_items; n_mol++) {
    /* Normal molecule fields */
    unsigned int external_species_id;
    byte act_newbie_flag;
    byte act_change_flag;
    double sched_time;
    double lifetime;
    double birthday;
    double x_coord, y_coord, z_coord;
    int orient;

    /* read molecule fields */
    READUINT(external_species_id);
    READFIELDRAW(act_newbie_flag);
    READFIELDRAW(act_change_flag);
    READFIELD(sched_time);
    READFIELD(lifetime);
    READFIELD(birthday);
    READFIELD(x_coord);
    READFIELD(y_coord);
    READFIELD(z_coord);
    READINT(orient);

    // starting with API version 1, convert the sched_time, lifetime and
    // birthday into scaled time based on the current timestep
    if (api_version >= 1) {
      // This will force lifetimes to be recomputed. This is necessary if
      // unimolecular rate constants change between checkpoints.
      lifetime = 0;
      sched_time = world->start_iterations;
      act_change_flag = HAS_ACT_CHANGE;
    }

    unsigned int complex_no = 0;
    READUINT(complex_no);

    /* Find this species by its external species id */
    struct species *properties = NULL;
    for (int species_idx = 0; species_idx < world->n_species; species_idx++) {
      if (world->species_list[species_idx]->chkpt_species_id ==
          external_species_id) {
        properties = world->species_list[species_idx];
        break;
      }
    }
    DATACHECK(properties == NULL,
              "Found molecule with unknown species id (%d).",
              external_species_id);

    /* Create and add molecule to scheduler */
    struct periodic_image periodic_box = { .x = 0,
                                           .y = 0,
                                           .z = 0
                                         };
    if ((properties->flags & NOT_FREE) == 0) { /* 3D molecule */

      /* set molecule characteristics */
      amp->t = sched_time;
      amp->t2 = lifetime;
      amp->birthday = birthday;
      amp->properties = properties;
      vmp->previous_wall = NULL;
      vmp->index = -1;
      vmp->pos.x = x_coord;
      vmp->pos.y = y_coord;
      vmp->pos.z = z_coord;
      amp->periodic_box = &periodic_box;

      /* Set molecule flags */
      amp->flags = TYPE_VOL | IN_VOLUME;
      if (act_newbie_flag == HAS_ACT_NEWBIE)
        amp->flags |= ACT_NEWBIE;

      if (act_change_flag == HAS_ACT_CHANGE)
        amp->flags |= ACT_CHANGE;

      amp->flags |= IN_SCHEDULE;
      if ((amp->properties->flags & CAN_SURFWALL) != 0 ||
          trigger_unimolecular(world->reaction_hash, world->rx_hashsize,
                               amp->properties->hashval, amp) != NULL)
        amp->flags |= ACT_REACT;
      if (amp->properties->space_step > 0.0)
        amp->flags |= ACT_DIFFUSE;

      /* Insert copy of vm into world */
      guess = insert_volume_molecule(world, vmp, guess);
      if (guess == NULL) {
        mcell_error("Cannot insert copy of molecule of species '%s' into "
                    "world.\nThis may be caused by a shortage of memory.",
                    vmp->properties->sym->name);
      }

    } else { /* surface_molecule */
      struct vector3 where;

      where.x = x_coord;
      where.y = y_coord;
      where.z = z_coord;


      struct surface_molecule *smp = insert_surface_molecule(
          world, properties, &where, orient, CHKPT_GRID_TOLERANCE, sched_time,
          NULL, NULL, NULL, &periodic_box);

      if (smp == NULL) {
        mcell_warn("Could not place molecule %s at (%f,%f,%f).",
                   properties->sym->name, where.x * world->length_unit,
                   where.y * world->length_unit,
                   where.z * world->length_unit);
        continue;
      }

      smp->t2 = lifetime;
      smp->birthday = birthday;
      if (act_newbie_flag == HAS_NOT_ACT_NEWBIE)
        smp->flags &= ~ACT_NEWBIE;

      if (act_change_flag == HAS_ACT_CHANGE) {
        smp->flags |= ACT_CHANGE;
      }
    }
  }

  return 0;
}
