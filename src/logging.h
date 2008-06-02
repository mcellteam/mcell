#ifndef INCLUDED_LOGGING_H
#define INCLUDED_LOGGING_H

#include <stdarg.h>
#include <stdio.h>

#ifndef __GNUC__
#ifndef __attribute__
#define __attribute__(x) /* empty */
#endif
#endif

#if __GNUC__ < 3
#ifndef __attribute__
#define __attribute__(x) /* empty */
#endif
#endif

#ifdef DEBUG
#define no_printf(fmt, ...) printf(fmt, ## __VA_ARGS__)
#else
#define no_printf(fmt, ...) do { /* do nothing */ } while(0)
#define printf(fmt, ...) DO_NOT_USE_PRINTF(fmt, ## __VA_ARGS__)
#endif

#define UNHANDLED_CASE(v) mcell_internal_error("Unhandled case in switch statement (%d).", v)

/* Get the log file. */
FILE *mcell_get_log_file(void);

/* Get the error file. */
FILE *mcell_get_error_file(void);

/* Set the log file. */
void mcell_set_log_file(FILE *f);

/* Set the error file. */
void mcell_set_error_file(FILE *f);

/********************************************************
 * Raw I/O to log and error streams
 ********************************************************/

/* Log a message. */
void mcell_log_raw(char const *fmt, ...)
  __attribute__((format (printf, 1, 2)));

/* Log a message (va_list version). */
void mcell_logv_raw(char const *fmt, va_list args);

/* Log a message. */
void mcell_error_raw(char const *fmt, ...)
  __attribute__((format (printf, 1, 2)));

/* Log a message (va_list version). */
void mcell_errorv_raw(char const *fmt, va_list args);

/********************************************************
 * Standard I/O to log and error streams -- may add standard prefix/suffix to
 * lines, may add newlines, etc.  May terminate program execution.
 ********************************************************/

/* Log a message. */
void mcell_log(char const *fmt, ...)
  __attribute__((format (printf, 1, 2)));

/* Log a message (va_list version). */
void mcell_logv(char const *fmt, va_list args);

/* Log a warning. */
void mcell_warn(char const *fmt, ...)
  __attribute__((format (printf, 1, 2)));

/* Log a warning (va_list version). */
void mcell_warnv(char const *fmt, va_list args);

/* Log an error and carry on. */
void mcell_error_nodie(char const *fmt, ...)
  __attribute__((format (printf, 1, 2)));

/* Log an error and carry on (va_list version). */
void mcell_errorv_nodie(char const *fmt, va_list args);

/* Log an error and exit. */
void mcell_error(char const *fmt, ...)
  __attribute__((format (printf, 1, 2), noreturn));

/* Log an error and exit (va_list version). */
void mcell_errorv(char const *fmt, va_list args)
  __attribute__((noreturn));

/* Log an internal error and exit. */
#define mcell_internal_error(fmt, ...)                                      \
  mcell_internal_error_(__FILE__, __LINE__, __func__, fmt, ## __VA_ARGS__)
void mcell_internal_error_(char const *file,
                           unsigned int line,
                           char const *func,
                           char const *fmt, ...)
  __attribute__((format (printf, 4, 5), noreturn));

/* Log an internal error and exit (va_list version). */
#define mcell_internal_errorv(fmt, args)                                    \
  mcell_internal_errorv_(__FILE__, __LINE__, __func__, fmt, args)
void mcell_internal_errorv_(char const *file,
                            unsigned int line,
                            char const *func,
                            char const *fmt,
                            va_list args)
  __attribute__((noreturn));

/* Get a copy of a string giving an error message. */
char *mcell_strerror(int err);

/* Log an error due to a failed standard library call, but do not exit. */
void mcell_perror_nodie(int err, char const *fmt, ...)
  __attribute__((format (printf, 2, 3)));

/* Log an error due to a failed standard library call, but do not exit (va_list version). */
void mcell_perrorv_nodie(int err, char const *fmt, va_list args);

/* Log an error due to a failed standard library call, and exit. */
void mcell_perror(int err, char const *fmt, ...)
  __attribute__((format (printf, 2, 3), noreturn));

/* Log an error due to a failed standard library call, and exit (va_list version). */
void mcell_perrorv(int err, char const *fmt, va_list args)
  __attribute__((noreturn));

/* Log an error due to failed memory allocation, but do not exit. */
void mcell_allocfailed_nodie(char const *fmt, ...)
  __attribute__((format (printf, 1, 2)));

/* Log an error due to failed memory allocation, but do not exit (va_list version). */
void mcell_allocfailedv_nodie(char const *fmt, va_list args);

/* Log an error due to failed memory allocation, and exit. */
void mcell_allocfailed(char const *fmt, ...)
  __attribute__((format (printf, 1, 2), noreturn));

/* Log an error due to failed memory allocation, and exit (va_list version). */
void mcell_allocfailedv(char const *fmt, va_list args)
  __attribute__((noreturn));

/* Terminate program execution due to an error. */
void mcell_die(void)
  __attribute__((noreturn));

#endif
