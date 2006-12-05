#include "argparse.h"

#include "mcell_structs.h"  /* for struct volume */
#include <stdarg.h>         /* for va_start, va_end, va_list */
#include <string.h>         /* for strdup */
#include <getopt.h>         /* for getopt_long_only, struct option, ... */
#include <stdio.h>          /* for *printf functions */
#include <stdlib.h>         /* for strtol, strtoll, strtoul, free */

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
static struct option long_options[] = {
  {"help",              0, 0, 'h'},
  {"seed",              1, 0, 's'},
  {"iterations",        1, 0, 'i'},
  {"checkpoint_infile", 1, 0, 'c'},
  {"logfile",           1, 0, 'l'},
  {"logfreq",           1, 0, 'f'},
  {"errfile",           1, 0, 'e'}
};

/* argerror: Display a message about an error which occurred during the
 *           argument parsing.  The message goes into the current "err_file",
 *           which defaults to stderr.
 *
 *   vol: the volume into which to imbue the parsed options
 *   fmt: a C "printf"-style format string
 */
void argerror(struct volume *vol, char const *fmt, ...)
{
  va_list args;
  fprintf(vol->err_file,"\nMCell: command-line argument syntax error: ");

  va_start(args, fmt);
  vfprintf(vol->err_file, fmt, args);
  va_end(args);

  fprintf(vol->err_file, "\n");
  fflush(vol->err_file);
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
int argparse_init(int argc, char * const argv[], struct volume *vol)
{
  char *endptr = NULL;

  /* Set up default values */
  vol->log_freq = -1;
  vol->log_file_name = NULL;
  vol->err_file_name = NULL;
  vol->log_file = stdout;
  vol->err_file = stderr;
  vol->seed_seq = 1;
  vol->mdl_infile_name = NULL;

  /* Loop over all arguments */
  while (1)
  {

    /* get the next argument */
    int c = getopt_long_only(argc, argv, "?hs:i:c:l:f:e:", long_options, NULL);
    if (c == -1)
        break;

    switch (c)
    {
      case '?':
      case 'h':  /* -help */
        return 1;

      case 's':  /* -seed */
        vol->seed_seq = (int) strtol(optarg, &endptr, 0);
        if (endptr == optarg || *endptr != '\0') {
          argerror(vol, "Random seed must be an integer: %s", optarg);
          return 1;
        }

#ifdef USE_RAN4
        if (vol->seed_seq < 1 || vol->seed_seq > 3000) {
          argerror(vol, "Random sequence number %d not in range 1 to 3000", vol->seed_seq);
          return 1;
        }
#endif
        break;

      case 'i':  /* -iterations */
        vol->iterations = strtoll(optarg, &endptr, 0);
        if (endptr == optarg || *endptr != '\0') {
          argerror(vol, "Iteration count must be an integer: %s", optarg);
          return 1;
        }

        if (vol->iterations < 0) {
          argerror(vol, "Iteration count %lld is less than 0", (long long int) vol->iterations);
          return 1;
        }
        break;

      case 'c':  /* -checkpoint_infile */
        vol->chkpt_infile = strdup(optarg);
        if (vol->chkpt_infile == NULL) {
          argerror(vol, "File '%s', Line %u: Out of memory while parsing command-line arguments: %s\n", __FILE__, __LINE__, optarg);
          return 1;
        }

        if ((vol->chkpt_infs = fopen(vol->chkpt_infile, "rb")) == NULL) {
          argerror(vol, "Cannot open input checkpoint file: %s", vol->chkpt_infile);
          free(vol->chkpt_infile);
          vol->chkpt_infile = NULL;
          vol->chkpt_init = 1;
          return 1;
        }

        vol->chkpt_init=0;
        vol->chkpt_flag = 1;
        fclose(vol->chkpt_infs);
        vol->chkpt_infs=NULL;
        break;

      case 'l':  /* -logfile */
        if (vol->log_file_name != NULL) {
          argerror(vol, "-l or --logfile argument specified more than once: %s", optarg);
          return 1;
        }

        vol->log_file_name = strdup(optarg);
        if (vol->log_file_name == NULL) {
          argerror(vol, "File '%s', Line %u: Out of memory while parsing command-line arguments: %s\n", __FILE__, __LINE__, optarg);
          return 1;
        }

        if ((vol->log_file=fopen(vol->log_file_name, "w")) == NULL) {
          vol->log_file = stdout;
          argerror(vol, "Cannot open output log file: %s", vol->log_file_name);
          free(vol->log_file_name);
          vol->log_file_name = NULL;
          return 1;
        }
        break;

      case 'f':  /* -logfreq */
        if (vol->log_freq != (unsigned long) -1) {
          argerror(vol, "-f or --logfreq specified more than once: %s", optarg);
          return 1;
        }

        vol->log_freq = (int) strtoul(optarg, &endptr, 0);
        if (endptr == optarg || *endptr != '\0') {
          argerror(vol, "Logging frequency must be an integer: %s", optarg);
          return 1;
        }

        if (vol->log_freq < 1) {
          argerror(vol, "Iteration report update interval must be at least 1 iteration");
          return 1;
        }
        break;

      case 'e':  /* -errfile */
        if (vol->err_file_name != NULL) {
          argerror(vol, "-e or --errfile argument specified more than once");
          return 1;
        }

        vol->err_file_name = strdup(optarg);
        if (vol->err_file_name == NULL) {
          argerror(vol, "File '%s', Line %ld: Out of memory while parsing command line arguments: %s", __FILE__, (long)__LINE__, optarg);
          return 1;
        }

        if ((vol->err_file = fopen(vol->err_file_name, "w")) == NULL) {

          /* what should happen in this case is debatable.  stderr, stdout,
           * or log_file...
           */
          vol->err_file = stderr;
          argerror(vol, "Cannot open output error file: %s", vol->err_file_name);
          free(vol->err_file_name);
          vol->err_file_name = NULL;
          return 1;
        }
        break;

      default:
        argerror(vol, "Internal error: getopt returned character code 0x%02x", (unsigned int) c);
        return 1;
    }
  }

  /* The old code would set err_file to log_file if err_file was NULL.  As
   * far as I can tell, this could never occur.  Perhaps this behavior is
   * what was intended?
   */
  /*
     if (vol->err_file == stderr  &&  vol->log_file != stdout)
         vol->err_file = vol->log_file;
   */

  /* Handle any left-over arguments, which we assume to be MDL files. */
  if (optind < argc) {
    if (argc - optind > 1) {
      argerror(vol, "%d MDL file names specified: %s, %s, ...", argc - optind, argv[optind], argv[optind+1]);
      return 1;
    }

    vol->mdl_infile_name = strdup(argv[optind]);
    if (vol->mdl_infile_name == NULL) {
      argerror(vol, "File '%s', Line %ld: Out of memory while parsing command line arguments: %s", __FILE__, (long)__LINE__, argv[optind]);
      return 1;
    }
  }
  else {
    argerror(vol, "No MDL file name specified");
    return 1;
  }

  return 0;
}
