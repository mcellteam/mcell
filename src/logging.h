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

#pragma once

#include <stdarg.h>
#include <stdio.h>

#ifdef DEBUG
#define no_printf(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define no_printf(fmt, ...)                                                    \
  do { /* do nothing */                                                        \
  } while (0)
#ifdef printf
#undef printf
#endif
#define printf(fmt, ...) DO_NOT_USE_PRINTF(fmt, ##__VA_ARGS__)
#endif

#define UNHANDLED_CASE(v)                                                      \
  mcell_internal_error("Unhandled case in switch statement (%d).", v)

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
void mcell_log_raw(char const *fmt, ...) PRINTF_FORMAT(1);

/* Log a message (va_list version). */
void mcell_logv_raw(char const *fmt, va_list args) PRINTF_FORMAT_V(1);

/* Log a message. */
void mcell_error_raw(char const *fmt, ...) PRINTF_FORMAT(1);

/* Log a message (va_list version). */
void mcell_errorv_raw(char const *fmt, va_list args) PRINTF_FORMAT_V(1);

/********************************************************
 * Standard I/O to log and error streams -- may add standard prefix/suffix to
 * lines, may add newlines, etc.  May terminate program execution.
 ********************************************************/

/* Log a message. */
void mcell_log(char const *fmt, ...) PRINTF_FORMAT(1);

/* Log a message (va_list version). */
void mcell_logv(char const *fmt, va_list args) PRINTF_FORMAT_V(1);

/* Log a warning. */
void mcell_warn(char const *fmt, ...) PRINTF_FORMAT(1);

/* Log a warning (va_list version). */
void mcell_warnv(char const *fmt, va_list args) PRINTF_FORMAT_V(1);

/* Log an error and carry on. */
void mcell_error_nodie(char const *fmt, ...) PRINTF_FORMAT(1);

/* Log an error and carry on (va_list version). */
void mcell_errorv_nodie(char const *fmt, va_list args) PRINTF_FORMAT_V(1);

/* Log an error and exit. */
void mcell_error(char const *fmt, ...) PRINTF_FORMAT(1)
    __attribute__((noreturn));

/* Log an error and exit (va_list version). */
void mcell_errorv(char const *fmt, va_list args) PRINTF_FORMAT_V(1)
    __attribute__((noreturn));

/* Log an internal error and exit. */
#define mcell_internal_error(fmt, ...)                                         \
  mcell_internal_error_(__FILE__, __LINE__, __func__, fmt, ##__VA_ARGS__)
void mcell_internal_error_(char const *file, unsigned int line,
                           char const *func, char const *fmt, ...)
    PRINTF_FORMAT(4) __attribute__((noreturn));

/* Log an internal error and exit (va_list version). */
#define mcell_internal_errorv(fmt, args)                                       \
  mcell_internal_errorv_(__FILE__, __LINE__, __func__, fmt, args)
void mcell_internal_errorv_(char const *file, unsigned int line,
                            char const *func, char const *fmt, va_list args)
    PRINTF_FORMAT_V(4) __attribute__((noreturn));

/* Get a copy of a string giving an error message. */
char *mcell_strerror(int err);

/* Log an error due to a failed standard library call, but do not exit. */
void mcell_perror_nodie(int err, char const *fmt, ...) PRINTF_FORMAT(2);

/* Log an error due to a failed standard library call, but do not exit (va_list
 * version). */
void mcell_perrorv_nodie(int err, char const *fmt, va_list args)
    PRINTF_FORMAT_V(2);

/* Log an error due to a failed standard library call, and exit. */
void mcell_perror(int err, char const *fmt, ...) PRINTF_FORMAT(2)
    __attribute__((noreturn));

/* Log an error due to a failed standard library call, and exit (va_list
 * version). */
void mcell_perrorv(int err, char const *fmt, va_list args) PRINTF_FORMAT_V(2)
    __attribute__((noreturn));

/* Log an error due to failed memory allocation, but do not exit. */
void mcell_allocfailed_nodie(char const *fmt, ...) PRINTF_FORMAT(1);

/* Log an error due to failed memory allocation, but do not exit (va_list
 * version). */
void mcell_allocfailedv_nodie(char const *fmt, va_list args) PRINTF_FORMAT_V(1);

/* Log an error due to failed memory allocation, and exit. */
void mcell_allocfailed(char const *fmt, ...) PRINTF_FORMAT(1)
    __attribute__((noreturn));

/* Log an error due to failed memory allocation, and exit (va_list version). */
void mcell_allocfailedv(char const *fmt, va_list args) PRINTF_FORMAT_V(1)
    __attribute__((noreturn));

/* Terminate program execution due to an error. */
void mcell_die(void) __attribute__((noreturn));
