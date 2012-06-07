#include "logging.h"
#include "mem_util.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#include <string.h>

/* Our log file */
static FILE *mcell_log_file = NULL;

/* Our warning/error file */
static FILE *mcell_error_file = NULL;

/* Get the log file. */
FILE *mcell_get_log_file(void)
{
  if (mcell_log_file == NULL)
    mcell_log_file = stdout;
  return mcell_log_file;
}

/* Get the error file. */
FILE *mcell_get_error_file(void)
{
  if (mcell_error_file == NULL)
    mcell_error_file = stderr;
  return mcell_error_file;
}

/* Set the log file. */
void mcell_set_log_file(FILE *f)
{
  if (mcell_log_file != NULL    &&
      mcell_log_file != stdout  &&
      mcell_log_file != stderr)
    fclose(mcell_log_file);

  mcell_log_file = f;
  setvbuf(mcell_log_file, NULL, _IOLBF, 128);
}

/* Set the error file. */
void mcell_set_error_file(FILE *f)
{
  if (mcell_error_file != NULL    &&
      mcell_error_file != stdout  &&
      mcell_error_file != stderr)
    fclose(mcell_error_file);
  mcell_error_file = f;
  setvbuf(mcell_error_file, NULL, _IOLBF, 128);
}

/* Log a message. */
void mcell_log_raw(char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_logv_raw(fmt, args);
  va_end(args);
}

/* Log a message (va_list version). */
void mcell_logv_raw(char const *fmt, va_list args)
{
  vfprintf(mcell_get_log_file(), fmt, args);
}

/* Log a message. */
void mcell_error_raw(char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_errorv_raw(fmt, args);
  va_end(args);
}

/* Log a message (va_list version). */
void mcell_errorv_raw(char const *fmt, va_list args)
{
  vfprintf(mcell_get_error_file(), fmt, args);
}

/* Log a message. */
void mcell_log(char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_logv(fmt, args);
  va_end(args);
}

/* Log a message (va_list version). */
void mcell_logv(char const *fmt, va_list args)
{
  mcell_logv_raw(fmt, args);
  fprintf(mcell_get_log_file(), "\n");
}

/* Log a warning. */
void mcell_warn(char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_warnv(fmt, args);
  va_end(args);
}

/* Log a warning (va_list version). */
void mcell_warnv(char const *fmt, va_list args)
{
  fprintf(mcell_get_error_file(), "Warning: ");
  mcell_errorv_raw(fmt, args);
  fprintf(mcell_get_error_file(), "\n");
}

/* Log an error and carry on. */
void mcell_error_nodie(char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_errorv_nodie(fmt, args);
  va_end(args);
}

/* Log an error and carry on (va_list version). */
void mcell_errorv_nodie(char const *fmt, va_list args)
{
  fprintf(mcell_get_error_file(), "Fatal error: ");
  mcell_errorv_raw(fmt, args);
  fprintf(mcell_get_error_file(), "\n");
}

/* Log an error and exit. */
void mcell_error(char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_errorv(fmt, args);
  va_end(args);
}

/* Log an error and exit (va_list version). */
void mcell_errorv(char const *fmt, va_list args)
{
  mcell_errorv_nodie(fmt, args);
  mcell_die();
}

/* Log an internal error and exit. */
void mcell_internal_error_(char const *file,
                           unsigned int line,
                           char const *func,
                           char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_internal_errorv_(file, line, func, fmt, args);
  va_end(args);
}

/* Log an error and exit (va_list version). */
void mcell_internal_errorv_(char const *file,
                           unsigned int line,
                           char const *func,
                           char const *fmt,
                           va_list args)
{
  fprintf(mcell_get_error_file(), "****************\n");
  fprintf(mcell_get_error_file(), "INTERNAL ERROR at %s:%u [%s]: ", file, line, func);
  mcell_errorv_raw(fmt, args);
  fprintf(mcell_get_error_file(), "\n");
  fprintf(mcell_get_error_file(), "MCell has detected an internal program error.\n");
  fprintf(mcell_get_error_file(), "Please report this error to the MCell developers at <%s>.\n", PACKAGE_BUGREPORT);
  fprintf(mcell_get_error_file(), "****************\n");
  mcell_die();
}

/* Get a copy of a string giving an error message. */
char *mcell_strerror(int err)
{
  char buffer[2048];
#ifdef STRERROR_R_CHAR_P
  char *pbuf = strerror_r(err, buffer, sizeof(buffer));
  if (pbuf != NULL)
    return CHECKED_STRDUP(pbuf, "error description");
  else
    return CHECKED_SPRINTF("UNIX error code %d.", err);
#else
  if (strerror_r(err, buffer, sizeof(buffer)) == 0)
    return CHECKED_STRDUP(buffer, "error description");
  else
    return CHECKED_SPRINTF("UNIX error code %d.", err);
#endif
}

/* Log an error due to a failed standard library call, and exit. */
void mcell_perror_nodie(int err, char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_perrorv_nodie(err, fmt, args);
  va_end(args);
}

/* Log an error due to a failed standard library call, and exit (va_list version). */
void mcell_perrorv_nodie(int err, char const *fmt, va_list args)
{
  char buffer[2048];
  fprintf(mcell_get_error_file(), "Fatal error: ");
  mcell_errorv_raw(fmt, args);
#ifdef STRERROR_R_CHAR_P
  fprintf(mcell_get_error_file(), ": %s\n", strerror_r(err, buffer, sizeof(buffer)));
#else
  if (strerror_r(err, buffer, sizeof(buffer)) == 0)
    fprintf(mcell_get_error_file(), ": %s\n", buffer);
  else
    fprintf(mcell_get_error_file(), "\n");
#endif
}

/* Log an error due to a failed standard library call, and exit. */
void mcell_perror(int err, char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_perrorv(err, fmt, args);
  va_end(args);
}

/* Log an error due to a failed standard library call, and exit (va_list version). */
void mcell_perrorv(int err, char const *fmt, va_list args)
{
  mcell_perrorv_nodie(err, fmt, args);
  mcell_die();
}

/* Log an error due to failed memory allocation, but do not exit. */
void mcell_allocfailed_nodie(char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_allocfailedv_nodie(fmt, args);
  va_end(args);
}

/* Log an error due to failed memory allocation, but do not exit (va_list version). */
void mcell_allocfailedv_nodie(char const *fmt, va_list args)
{
  fprintf(mcell_get_error_file(), "Fatal error: ");
  mcell_errorv_raw(fmt, args);
  fprintf(mcell_get_error_file(), "\n");
  fprintf(mcell_get_error_file(), "Fatal error: Out of memory\n\n");
  mem_dump_stats(mcell_get_error_file());
}

/* Log an error due to failed memory allocation, and exit. */
void mcell_allocfailed(char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  mcell_allocfailedv_nodie(fmt, args);
  va_end(args);
  mcell_die();
}

/* Log an error due to failed memory allocation, and exit (va_list version). */
void mcell_allocfailedv(char const *fmt, va_list args)
{
  mcell_allocfailedv_nodie(fmt, args);
  mcell_die();
}

/* Terminate program execution due to an error. */
void mcell_die(void)
{
  exit(EXIT_FAILURE);
}
