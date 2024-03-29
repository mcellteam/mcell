/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "config.h"

#include "argparse.h"

#include "mcell_structs.h" /* for struct volume */
#include "logging.h"
#include "version_info.h" /* for print_version, print_full_version */

#include <stdarg.h>       /* for va_start, va_end, va_list */
#include <string.h>       /* for strdup */
#include <stdio.h>        /* for *printf functions */
#include <stdlib.h>       /* for strtol, strtoll, strtoul, free */

#ifndef _MSC_VER
#include <getopt.h>
#else
#include "win_getopt/win_getopt.h"
#endif

#include <nfsim_c.h> /* for the nfsim initialization stuff */


/* Display a formatted error message. */
static void argerror(char const *s, ...) PRINTF_FORMAT(1);

/* Command-line arguments structure:
 *     long arg name
 *     has argument
 *     pointer to flag (should always be 0)
 *     short argument letter
 *
 * To add an argument, add a new entry to this array, assigning it a unique
 * short-argument letter.  Then add the short argument letter to the gibberish
 * string passed to getopt_long_only below, followed by a colon if the flag
 * takes an argument.  Finally, add a case to the switch statement to handle
 * the letter you selected for the argument.
 */
static struct option long_options[] = { { "help", 0, 0, 'h' },
                                        { "version", 0, 0, 'v' },
                                        { "fullversion", 0, 0, 'V' },
                                        { "seed", 1, 0, 's' },
                                        { "iterations", 1, 0, 'i' },
                                        { "checkpoint_infile", 1, 0, 'c' },
                                        { "checkpoint_outfile", 1, 0, 'C' },
                                        { "logfile", 1, 0, 'l' },
                                        { "logfreq", 1, 0, 'f' },
                                        { "errfile", 1, 0, 'e' },
                                        { "bond_angle", 1, 0, 'b'},
                                        { "z_options", 1, 0, 'z' },
                                        { "dump", 1, 0, 'd'},
                                        { "quiet", 0, 0, 'q' },
                                        { "with_checks", 1, 0, 'w' },
                                        { "rules", 1, 0, 'r'},
																				{ "mcell4", 0, 0, 'n'},
																				{ "dump_mcell3", 0, 0, 't'},
																				{ "dump_mcell4", 0, 0, 'o'},
																				{ "dump_mcell4_with_geometry", 0, 0, 'g'},
                                        { "mdl2datamodel4", 0, 0, 'u'},
                                        { "mdl2datamodel4viz", 0, 0, 'a'},
                                        { NULL, 0, 0, 0 } };

/* print_usage: Write the usage message for mcell to a file handle.
 *
 *    f:     file handle to which to write message
 *    argv0: command name to use for executable in message
 */
void print_usage(FILE *f, char const *argv0) {
  fprintf(f, "Usage: %s [options] mdl_file_name\n\n", argv0);
  fprintf(
      f,
      "  options:\n"
      "     [-help]                  print this help message\n"
      "     [-version]               print the program version and exit\n"
      "     [-fullversion]           print the detailed program version report and exit\n"
      "     [-seed n]                choose random sequence number (default: 1)\n"
      "     [-iterations n]          override iterations in mdl_file_name\n"
      "     [-logfile log_file_name] send output log to file (default: stdout)\n"
      "     [-logfreq n]             output log frequency\n"
      "     [-errfile err_file_name] send errors log to file (default: stderr)\n"
      "     [-checkpoint_infile checkpoint_file_name]   read checkpoint file\n"
      "     [-checkpoint_outfile checkpoint_file_name]  write checkpoint file\n"
      "     [-z_options opt_int]     additional visualization options (defaults to 0 which is none)\n"
      "     [-bond_angle angle]      bond angle to use for all bonds (defaults to 0)\n"
      "     [-dump level]            print additional information based on level (0 is none, >0 is more)\n"
      "     [-quiet]                 suppress all unrequested output except for errors\n"
      "     [-with_checks ('yes'/'no', default 'yes')]   performs check of the geometry for coincident walls\n"
      "     [-rules rules_file_name] run in MCell-R mode\n"
			"     [-mcell4]                run new experimental MCell 4 version\n"
      "     [-dump_mcell3]           dump initial MCell 3 state for MCell 4 development\n"
			"     [-dump_mcell4]           dump initial MCell 4 state without geometry\n"
      "     [-dump_mcell4_with_geometry] dump initial MCell 4 state with geometry\n"
      "     [-mdl2datamodel4]        convert MDL to datamodel using mcell 4 state, the resulting file will be called 'data_model.json'\n"
      "     [-mdl2datamodel4viz]     convert MDL to datamodel using mcell 4 state, only for visualization purposes the resulting file will be called 'data_model_viz.json'\n"
      "\n");
}

/* argerror: Display a message about an error which occurred during the
 *           argument parsing.  The message goes into the current "err_file",
 *           which defaults to stderr.
 *
 *   vol: the volume into which to imbue the parsed options
 *   fmt: a C "printf"-style format string
 */
static void argerror(char const *fmt, ...) {

  mcell_error_raw("\nMCell: command-line argument syntax error: ");

  va_list args;
  va_start(args, fmt);
  char error_msg[256];
  vsnprintf(error_msg, sizeof(error_msg), fmt, args);
  va_end(args);

  mcell_error_raw("%s", error_msg);
  mcell_error_raw("\n");
}

/* argparse_init: Parse the command-line arguments, imbuing the options into
 *                'vol'.
 *
 *   argc: the number of command-line arguments including the program name
 *   argv: an array of the command-line arguments
 *   vol:  the volume into which to imbue our options
 *
 *   Returns 1 if for any reason the simulation should not be run (i.e. an
 *   error, or '-help'), 0 if the simulation should proceed.  If 1 is returned,
 *   the caller should display a usage message.  (Perhaps the usage message
 *   should actually be displayed from here?)
 */
int argparse_init(int argc, char *const argv[], struct volume *vol) {
  char *endptr = NULL;
  int log_file_specified = 0, err_file_specified = 0;
  FILE *fhandle = NULL;
  char *with_checks_option;
  char *rules_xml_file = NULL; // for nfsim
  //argerror('whats up');
  /* Loop over all arguments */
  while (1) {

    /* get the next argument */
    int c = getopt_long_only(argc, argv, "?h", long_options, NULL);
    if (c == -1)
      break;

    switch (c) {
    case '?':
    case 'h': /* -help */
      return 1;

    case 'v': /* -version */
      print_version(mcell_get_log_file());
      exit(1);

    case 'V': /* -fullversion */
      print_full_version(mcell_get_log_file());
      exit(1);

    case 'q': /* -quiet */
      vol->quiet_flag = 1;
      break;

    case 'd': /* -dump */
      vol->dump_level = strtol(optarg, &endptr, 0);
      if (endptr == optarg || *endptr != '\0') {
        argerror("Dump level must be an integer: %s", optarg);
        return 1;
      }

      if (vol->dump_level < 0) {
        argerror("Dump level %ld is less than 0", vol->dump_level);
        return 1;
      }

      if (vol->dump_level > 0) {
        fprintf(stdout, "Dump level has been set to %ld\n", vol->dump_level);
      }
      break;

    case 'b': /* -bond_angle */
      vol->bond_angle = (double)strtod(optarg, &endptr);
      if (endptr == optarg || *endptr != '\0') {
        argerror("Bond angle must be a double: %s", optarg);
        return 1;
      }
      break;

    case 'z': /* -viz_options */
      vol->viz_options = strtol(optarg, &endptr, 0);
      fprintf ( stdout, "Parsed Visualization Options = %lx\n", vol->viz_options );
      if (endptr == optarg || *endptr != '\0') {
        argerror("viz_options must be an integer: %s", optarg);
        return 1;
      }
      break;

    case 'w': /* walls coincidence check (maybe other checks in future) */
      with_checks_option = strdup(optarg);
      if (with_checks_option == NULL) {
        argerror("File '%s', Line %u: Out of memory while parsing "
                 "command-line arguments: %s\n",
                 __FILE__, __LINE__, optarg);
        return 1;
      }
      if (((strcmp(with_checks_option, "yes") == 0) ||
           (strcmp(with_checks_option, "YES") == 0))) {
        vol->with_checks_flag = 1;
      } else if (((strcmp(with_checks_option, "no") == 0) ||
                  (strcmp(with_checks_option, "NO") == 0))) {
        vol->with_checks_flag = 0;
      } else {
        argerror("-with_checks option should be 'yes' or 'no'.");
        free(with_checks_option);
        return 1;
      }
      free(with_checks_option);
      break;

    case 's': /* -seed */
      vol->seed_seq = (int)strtol(optarg, &endptr, 0); // initialized to '1'
      if (endptr == optarg || *endptr != '\0') {
        argerror("Random seed must be an integer: %s", optarg);
        return 1;
      }
      break;

    case 'i': /* -iterations */
      vol->iterations = strtoll(optarg, &endptr, 0);
      if (endptr == optarg || *endptr != '\0') {
        argerror("Iteration count must be an integer: %s", optarg);
        return 1;
      }

      if (vol->iterations < 0) {
        argerror("Iteration count %lld is less than 0",
                 (long long int)vol->iterations);
        return 1;
      }
      break;

    case 'c': /* -checkpoint_infile */
      vol->chkpt_infile = strdup(optarg);
      if (vol->chkpt_infile == NULL) {
        argerror("File '%s', Line %u: Out of memory while parsing "
                 "command-line arguments: %s\n",
                 __FILE__, __LINE__, optarg);
        return 1;
      }

      if ((fhandle = fopen(vol->chkpt_infile, "rb")) == NULL) {
        argerror("Cannot open input checkpoint file: %s", vol->chkpt_infile);
        free(vol->chkpt_infile);
        vol->chkpt_infile = NULL;
        vol->chkpt_init = 1;
        return 1;
      }

      vol->chkpt_init = 0;
      vol->chkpt_flag = 1;
      fclose(fhandle);
      break;

    case 'C': /* -checkpoint_outfile */
      vol->chkpt_outfile = strdup(optarg);
      if (vol->chkpt_outfile == NULL) {
        argerror("File '%s', Line %u: Out of memory while parsing "
                 "command-line arguments: %s\n",
                 __FILE__, __LINE__, optarg);
        return 1;
      }

      vol->chkpt_flag = 1;
      break;

    case 'r': /* nfsim */
			vol->nfsim_flag = 1;
      rules_xml_file = strdup(optarg);
      break;

    case 'l': /* -logfile */
      if (log_file_specified) {
        argerror("-logfile argument specified more than once: %s", optarg);
        return 1;
      }

      if ((fhandle = fopen(optarg, "w")) == NULL) {
        argerror("Cannot open output log file: %s", optarg);
        return 1;
      } else {
        mcell_set_log_file(fhandle);
        log_file_specified = 1;
      }
      break;

    case 'f': /* -logfreq */
      if (vol->log_freq != ULONG_MAX) {
        argerror("-logfreq specified more than once: %s", optarg);
        return 1;
      }

      vol->log_freq = strtoul(optarg, &endptr, 0);
      if (endptr == optarg || *endptr != '\0') {
        argerror("Iteration report interval must be an integer: %s", optarg);
        return 1;
      }
      if (vol->log_freq == ULONG_MAX) {
        argerror("Iteration report interval must be an integer n such "
                 "that 1 <= n < %lu: %s",
                 ULONG_MAX, optarg);
        return 1;
      }
      if (vol->log_freq < 1) {
        argerror("Iteration report interval must be at least 1 iteration: %s",
                 optarg);
        return 1;
      }
      break;

    case 'e': /* -errfile */
      if (err_file_specified) {
        // XXX: This is almost identical to the -logfile error but doesn't work
        // for some reason.
        argerror("-errfile argument specified more than once: %s", optarg);
        return 1;
      }

      if ((fhandle = fopen(optarg, "w")) == NULL) {
        argerror("Cannot open output error file: %s", optarg);
        return 1;
      } else {
        mcell_set_error_file(fhandle);
        err_file_specified = 1;
      }
      break;

    case 'n':
      vol->use_mcell4 = 1;
      break;

    case 't':
      vol->dump_mcell3 = 1;
      break;

    case 'o':
      vol->dump_mcell4 = 1;
      break;

    case 'g':
      vol->dump_mcell4_with_geometry = 1;
      break;

    case 'u':
      vol->mdl2datamodel4 = 1;
      vol->mdl2datamodel4_only_viz = 0;
      break;

    case 'a':
      vol->mdl2datamodel4 = 1;
      vol->mdl2datamodel4_only_viz = 1;
      break;

    default:
      argerror("Internal error: getopt returned character code 0x%02x",
               (unsigned int)c);
      return 1;
    }
  }

  /* Handle any left-over arguments, which we assume to be MDL files. */
  if (optind < argc) {
    FILE *f;
    if (argc - optind > 1) {
      argerror("%d MDL file names specified: %s, %s, ...", argc - optind,
               argv[optind], argv[optind + 1]);
      return 1;
    }

    vol->mdl_infile_name = strdup(argv[optind]);
    if (vol->mdl_infile_name == NULL) {
      argerror("File '%s', Line %ld: Out of memory while parsing command "
               "line arguments: %s",
               __FILE__, (long)__LINE__, argv[optind]);
      return 1;
    }

    if ((f = fopen(argv[optind], "r")) == NULL) {
      /* XXX: Should probably use perror to explain why... */
      argerror("Cannot read MDL file: %s", argv[optind]);
      return 1;
    }
    fclose(f);
  } else {
    argerror("No MDL file name specified");
    return 1;
  }

  /* Initialize NFSim if requested */
  if (vol->nfsim_flag) {
    int nfsimStatus = setupNFSim_c(rules_xml_file, vol->seed_seq, vol->dump_level > 0);
    free(rules_xml_file);
    if (nfsimStatus != 0){
      argerror("nfsim model could not be properly initialized: %s", optarg);
      return 1;
    }
  }

  return 0;
}
