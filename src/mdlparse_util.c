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
 *USA.
 *
 ******************************************************************************/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <ctype.h>

#include "logging.h"
#include "strfunc.h"
#include "sym_table.h"
#include "mdlparse_util.h"
#include "react_output.h"
#include "diffuse_util.h"

#include "mcell_misc.h"
#include "mcell_structs.h"
#include "mcell_viz.h"
#include "mcell_release.h"
#include "dyngeom_parse_extras.h"
#include "mem_util.h"

extern void chkpt_signal_handler(int sn);

/* Free a variable value, leaving the symbol free for reassignment to another
 * type. */
static int mdl_free_variable_value(struct mdlparse_vars *parse_state,
                                   struct sym_entry *sym);

/*************************************************************************
 mdl_strip_quotes:
    Strips the quotes from a quoted string.

 In:  in: string to strip
 Out: stripped string, or NULL on error
*************************************************************************/
char *mdl_strip_quotes(char *in) {

  char *q = strip_quotes(in);
  free(in);
  if (q != NULL)
    return q;
  else {
    mcell_allocfailed("Failed to strip quotes from a quoted string.");
    return NULL;
  }
}

/*************************************************************************
 mdl_strcat:
    Concatenates two strings, freeing the original strings.

 In:  s1: first string to concatenate
      s2: second string to concatenate
 Out: conjoined string, or NULL if an error occurred
*************************************************************************/
char *mdl_strcat(char *s1, char *s2) {

  char *cat = my_strcat(s1, s2);
  free(s1);
  free(s2);
  if (cat != NULL)
    return cat;
  else {
    mcell_allocfailed("Failed to concatenate two strings.");
    return NULL;
  }
}

/*************************************************************************
 mdl_strdup:
    Duplicates a string.

 In:  s1: string to duplicate
 Out: allocated string, or NULL if an error occurred
*************************************************************************/
char *mdl_strdup(char const *s1) {

  char *s = strdup(s1);
  if (s1 == NULL)
    mcell_allocfailed("Failed to duplicate a string.");
  return s;
}

/*************************************************************************
 mdl_warning:
    Display a warning message about something encountered during the parse
    process.

 In:  parse_state: parser state
      fmt: Format string (as for printf)
      ...: arguments
 Out: Warning message is printed in the log_file
*************************************************************************/
void mdl_warning(struct mdlparse_vars *parse_state, char const *fmt, ...) {
  va_list arglist;
  if (parse_state->vol->procnum != 0)
    return;

  mcell_log_raw("MCell: Warning on line: %d of file: %s  ",
                parse_state->line_num[parse_state->include_stack_ptr - 1],
                parse_state->vol->curr_file);

  va_start(arglist, fmt);
  mcell_logv_raw(fmt, arglist);
  va_end(arglist);

  mcell_log_raw("\n");
}

/**************************************************************************
 mdl_valid_file_mode:
    Check that the speficied file mode string is valid for an fopen statement.

 In: parse_state: parser state
     mode: mode string to check
 Out: 1 if
**************************************************************************/
int mdl_valid_file_mode(struct mdlparse_vars *parse_state, char *mode) {
  char c = mode[0];
  if (c != 'r' && c != 'w' && c != 'a') {
    mdlerror_fmt(parse_state, "Invalid file mode: %s", mode);
    return 1;
  }
  if (c == 'r') {
    mdlerror(parse_state,
             "MCell models do not currently support opening files for reading");
    return 1;
  }

  return 0;
}

/**************************************************************************
 mdl_expr_log:
    Compute a natural log, reporting any range or domain errors.

 In:  parse_state: parser state
      in:   value whose log to compute
      out:  destination for computed value
 Out: 0 on success, 1 on failure.  *out is updated on success
**************************************************************************/
int mdl_expr_log(struct mdlparse_vars *parse_state, double in, double *out) {
  if (in <= 0.0) {
    mdlerror_fmt(parse_state,
                 "Attempt to compute LOG(%.15g), which is not defined.", in);
    return 1;
  }

  *out = log(in);
  return 0;
}

/**************************************************************************
 mdl_expr_log10:
    Compute a base-10 log, reporting any range or domain errors.

 In:  parse_state: parser state
      in:   value whose log to compute
      out:  destination for computed value
 Out: 0 on success, 1 on failure.  *out is updated on success
**************************************************************************/
int mdl_expr_log10(struct mdlparse_vars *parse_state, double in, double *out) {
  if (in <= 0.0) {
    mdlerror_fmt(parse_state,
                 "Attempt to compute LOG10(%.15g), which is not defined.", in);
    return 1;
  }

  *out = log10(in);
  return 0;
}

/**************************************************************************
 mdl_expr_mod:
    Compute a floating point modulo operator, reporting any domain errors.

 In:  parse_state: parser state
      in:   value
      divisor: divisor for modulo operator
      out:  destination for computed value
 Out: 0 on success, 1 on failure.  *out is updated on success
**************************************************************************/
int mdl_expr_mod(struct mdlparse_vars *parse_state, double in, double divisor,
                 double *out) {
  if (divisor == 0.0) {
    mdlerror_fmt(parse_state,
                 "Attempt to compute MOD(%.15g, 0.0), which is not defined.",
                 in);
    return 1;
  }
  *out = fmod(in, divisor);
  return 0;
}

/**************************************************************************
 mdl_expr_div:
    Compute a division, reporting any domain or range errors.

 In:  parse_state: parser state
      in:   value
      divisor: divisor
      out:  destination for computed value
 Out: 0 on success, 1 on failure.  *out is updated on success
**************************************************************************/
int mdl_expr_div(struct mdlparse_vars *parse_state, double in, double divisor,
                 double *out) {
  if (divisor == 0.0) {
    mdlerror_fmt(parse_state, "Attempt to divide %.15g by zero", in);
    return 1;
  }

  *out = in / divisor;
  if (isinf(*out)) {
    mdlerror_fmt(parse_state,
                 "Cannot compute %.15g / %.15g: result is too large", in,
                 divisor);
    return 1;
  }
  return 0;
}

/**************************************************************************
 mdl_expr_pow:
    Compute an exponentiation, reporting any domain or range errors.

 In:  parse_state: parser state
      in:   value
      exponent: exponent for exponentiation
      out:  destination for computed value
 Out: 0 on success, 1 on failure.  *out is updated on success
**************************************************************************/
int mdl_expr_pow(struct mdlparse_vars *parse_state, double in, double exponent,
                 double *out) {
  if (in < 0.0) {
    if (exponent != (int)exponent) {
      mdlerror_fmt(parse_state, "Cannot compute %.15g^%.15g: negative value "
                                "raised to non-integral power would be complex",
                   in, exponent);
      return 1;
    }
  }

  *out = pow(in, exponent);
  if (isinf(*out)) {
    mdlerror_fmt(parse_state, "Cannot compute %.15g^%.15g: result is too large",
                 in, exponent);
    return 1;
  }
  return 0;
}

/**************************************************************************
 mdl_expr_rng_uniform:
    Compute a uniform random variate.

 In:  parse_state: parser state
 Out: random variate, uniform on [0, 1)
**************************************************************************/
double mdl_expr_rng_uniform(struct mdlparse_vars *parse_state) {
  return rng_dbl(parse_state->vol->rng);
}

/**************************************************************************
 mdl_expr_roundoff:
    Round the input value off to n significant figures.

 In:  in:   value to round off
      ndigits: number of digits
 Out: value rounded to n digits
**************************************************************************/
double mdl_expr_roundoff(double in, int ndigits) {
  char fmt_string[1024];
  fmt_string[0] = '\0';
  snprintf(fmt_string, 1024, "%.*g", ndigits, in);
  return strtod(fmt_string, (char **)NULL);
}

/**************************************************************************
 mdl_expr_string_to_double:
    Convert a string value to a double, freeing the string value.

 In:  parse_state: parser state
      str:  string form of value
      out:  location to receive parsed value
 Out: 0 on success, 1 on failure.  *out is updated.
**************************************************************************/
int mdl_expr_string_to_double(struct mdlparse_vars *parse_state, char *str,
                              double *out) {
  *out = strtod(str, (char **)NULL);
  if (errno == ERANGE) {
    mdlerror_fmt(parse_state, "Error converting string to number: %s", str);
    free(str);
    return 1;
  }
  free(str);
  return 0;
}

/**************************************************************************
 mdl_new_filehandle:
    Create a new filehandle in the global symbol table.

 In:  parse_state: parser state
      name: name for file symbol
 Out: symbol or NULL on error
**************************************************************************/
struct sym_entry *mdl_new_filehandle(struct mdlparse_vars *parse_state,
                                     char *name) {
  struct sym_entry *sym;
  sym = retrieve_sym(name, parse_state->vol->fstream_sym_table);

  /* If this file is already open, close it. */
  if (sym != NULL) {
    struct file_stream *filep = (struct file_stream *)sym->value;
    if (filep->stream != NULL) {
      fclose(filep->stream);
      free(filep->name);
      filep->stream = NULL;
      filep->name = NULL;
    }
  }

  /* Otherwise, create it */
  else if ((sym = store_sym(name, FSTRM, parse_state->vol->fstream_sym_table,
                            NULL)) == NULL) {
    mdlerror_fmt(parse_state, "Out of memory while creating file stream: %s",
                 name);
    free(name);
    return NULL;
  }

  free(name);
  return sym;
}

/**************************************************************************
 mdl_fopen:
    Process an fopen statement, opening a new file handle.

 In:  parse_state: parser state
      filesym: symbol for the file
      name: filename
      mode: file open mode
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_fopen(struct mdlparse_vars *parse_state, struct sym_entry *filesym,
              char *name, char *mode) {
  struct file_stream *filep = (struct file_stream *)filesym->value;
  filep->name = name;
  if ((filep->stream = fopen(filep->name, mode)) == NULL) {
    mdlerror_fmt(parse_state, "Cannot open file: %s", filep->name);
    free(mode);
    /* XXX: Free filep? */
    return 1;
  }
  free(mode);

  return 0;
}

/**************************************************************************
 mdl_fclose:
    Process an fclose statement, closing an existing file handle.

 In:  parse_state: parser state
      filesym: symbol for the file
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_fclose(struct mdlparse_vars *parse_state, struct sym_entry *filesym) {
  struct file_stream *filep = (struct file_stream *)filesym->value;
  if (filep->stream == NULL)
    return 0;

  if (fclose(filep->stream) != 0) {
    filep->stream = NULL;
    mdlerror_fmt(parse_state, "Error closing file: %s", filep->name);
    free(filep->name);
    filep->name = NULL;
    return 1;
  }
  filep->stream = NULL;
  free(filep->name);
  filep->name = NULL;
  return 0;
}

/**************************************************************************
 double_dup:
 In:  parse_state: parser state
      value: value for newly allocated double
 Out: new double value that is copy of the input double value
**************************************************************************/
static double *double_dup(double value) {

  double *dup_value;
  dup_value = CHECKED_MALLOC_STRUCT(double, "numeric variable");
  if (dup_value != NULL)
    *dup_value = value;
  return dup_value;
}

/*************************************************************************
 mdl_new_printf_arg_double:
    Create a new double argument for a printf argument list.

 In: d: value for argument
 Out: The new argument, or NULL if an error occurred.
*************************************************************************/
struct arg *mdl_new_printf_arg_double(double d) {
  struct arg *a = CHECKED_MALLOC_STRUCT(struct arg, "format argument");
  if (a == NULL)
    return NULL;

  if ((a->arg_value = double_dup(d)) == NULL) {
    free(a);
    return NULL;
  }

  a->arg_type = DBL;
  a->next = NULL;
  return a;
}

/*************************************************************************
 mdl_new_printf_arg_string:
    Create a new string argument for a printf argument list.

 In: str: value for argument
 Out: The new argument, or NULL if an error occurred.
*************************************************************************/
struct arg *mdl_new_printf_arg_string(char const *str) {

  struct arg *a = CHECKED_MALLOC_STRUCT(struct arg, "format argument");
  if (a == NULL)
    return NULL;

  a->arg_value = (void *)str;
  a->arg_type = STR;
  a->next = NULL;
  return a;
}

/*************************************************************************
 mdl_free_printf_arg_list:
    Frees a printf argument list.

 In: args: list to free
 Out: all arguments are freed
*************************************************************************/
static void mdl_free_printf_arg_list(struct arg *args) {
  while (args != NULL) {
    struct arg *arg_next = args->next;
    if (args->arg_type == DBL)
      free(args->arg_value);
    free(args);
    args = arg_next;
  }
}

/*************************************************************************
 mdl_expand_string_escapes:
    Expands C-style escape sequences in the string.

 In:  in: string to expand
 Out: string with expanded escape sequences
*************************************************************************/
char *mdl_expand_string_escapes(char *in) {

  char *out;
  char *in_start = in;

  /* Allocate buffer for new string */
  char *expanded = CHECKED_MALLOC(strlen(in) + 1, "printf format string");
  if (expanded == NULL)
    return NULL;

  /* Copy expanded string to new buffer */
  out = expanded;
  while (*in != '\0') {
    if (*in != '\\')
      *out++ = *in++;
    else {
      ++in;
      if (*in == 'n')
        *out++ = '\n';
      else if (*in == 't')
        *out++ = '\t';
      else if (*in == 'e')
        *out++ = '\x1b';
      else if (*in == 'a')
        *out++ = '\a';
      else if (*in == 'f')
        *out++ = '\f';
      else if (*in == 'v')
        *out++ = '\v';
      else if (*in == 'b')
        *out++ = '\b';
      else if (*in == 'x') {
        if (isxdigit(in[1]) && isxdigit(in[2])) {
          char buffer[3];
          buffer[0] = in[1];
          buffer[1] = in[2];
          buffer[2] = '\0';
          *out++ = (char)strtol(buffer, NULL, 16);
          if (out[-1] == '\0')
            out[-1] = ' ';
          in += 2;
        } else
          *out++ = 'x';
      } else if (*in == '0') {
        if ('0' <= in[1] && in[1] <= '7' && '0' <= in[2] && in[2] <= '7' &&
            '0' <= in[3] && in[3] <= '7') {
          char buffer[4];
          buffer[0] = in[1];
          buffer[1] = in[2];
          buffer[2] = in[3];
          buffer[3] = '\0';
          *out++ = (char)strtol(buffer, NULL, 8);
          if (out[-1] == '\0')
            out[-1] = ' ';
          in += 3;
        } else
          *out++ = '0';
      } else if (*in == '\\')
        *out++ = '\\';
      else if (*in == '"')
        *out++ = '"';
      else if (*in == '\0') {
        *out++ = '\\';
        break;
      } else
        *out++ = *in;
      ++in;
    }
  }
  *out++ = '\0';
  free(in_start);
  return expanded;
}

/*************************************************************************
 Internal code used to designate the type required to fulfill a particular
 format specifier within a printf-style format string.
*************************************************************************/
enum {
  PRINTF_INVALID = -1,
  PRINTF_DOUBLE,
  PRINTF_SIGNED_CHAR,
  PRINTF_UNSIGNED_CHAR,
  PRINTF_SHORT_INT,
  PRINTF_U_SHORT_INT,
  PRINTF_INT,
  PRINTF_U_INT,
  PRINTF_LONG_INT,
  PRINTF_U_LONG_INT,
  PRINTF_LONG_LONG_INT,
  PRINTF_U_LONG_LONG_INT,
  PRINTF_STRING
};

/*************************************************************************
 get_printf_conversion_specifier:
    Find the details of the conversion specifier found in the given segment of
    format string.

 In:  parse_state: parse arguments
      fmt: format string to scan
      num_asterisks: counter for number of asterisks in conversion specifier
 Out: code indicating the type required for this conversion specifier, or
      PRINTF_INVALID if it's invalid or disallowed.
*************************************************************************/
static int get_printf_conversion_specifier(struct mdlparse_vars *parse_state,
                                           char const *fmt,
                                           int *num_asterisks) {
  static char CONVERSION_SPECIFIERS[] = "diouxXeEfFgGaAcCsSpnm";

  char const *const fmt_orig = fmt;
  char length = '\0';
  *num_asterisks = 0;
  for (++fmt; *fmt != '\0'; ++fmt) {
    if (isalpha(*fmt) && strchr(CONVERSION_SPECIFIERS, *fmt) != NULL)
      break;

    /* Scan any options in the string, bombing out if we find something
     * invalid, and keeping track of the length specifiers and asterisks we
     * find. */
    switch (*fmt) {
    case 'l':
      if (length == '\0')
        length = 'l';
      else if (length == 'l')
        length = 'L';
      else {
        mdlerror_fmt(parse_state, "Format string segment '%s' has an invalid "
                                  "length modifier inside a conversion "
                                  "specification.",
                     fmt_orig);
        return PRINTF_INVALID;
      }
      break;

    case 'h':
      if (length == '\0')
        length = 'h';
      else if (length == 'h')
        length = 'H';
      else {
        mdlerror_fmt(parse_state, "Format string segment '%s' has an invalid "
                                  "length modifier inside a conversion "
                                  "specification.",
                     fmt_orig);
        return PRINTF_INVALID;
      }
      break;

    case 'L':
    case 'q':
    case 'j':
    case 'z':
    case 't':
      mdlerror_fmt(parse_state, "Format string segment '%s' has an unsupported "
                                "length modifier (%c) inside a conversion "
                                "specification.",
                   fmt_orig, *fmt);
      return PRINTF_INVALID;

    case '*':
      ++*num_asterisks;
      if (*num_asterisks > 2) {
        mdlerror_fmt(parse_state, "Format string segment '%s' has more than "
                                  "two asterisks inside a conversion "
                                  "specification, which is unsupported by "
                                  "MCell.",
                     fmt_orig);
        return PRINTF_INVALID;
      }
      break;

    case '$':
      mdlerror_fmt(parse_state, "Format string segment '%s' has a '$' inside a "
                                "conversion specification, which is "
                                "unsupported by MCell.",
                   fmt_orig);
      return PRINTF_INVALID;

    default:
      /* Skip this char. */
      break;
    }
  }

  /* Filter invalid types */
  if (*fmt == '\0') {
    mdlerror_fmt(parse_state,
                 "Format string segment '%s' contains no conversion specifier.",
                 fmt_orig);
    return PRINTF_INVALID;
  } else if (*fmt == 'n') {
    mdlerror_fmt(parse_state, "Format string segment '%s' attempts to use "
                              "dangerous conversion specifier 'n'.",
                 fmt_orig);
    return PRINTF_INVALID;
  }

  /* Now, handle the format specifier itself */
  switch (*fmt) {
  case 'd':
  case 'i':
    if (length == '\0')
      return PRINTF_INT;
    else if (length == 'l')
      return PRINTF_LONG_INT;
    else if (length == 'L')
      return PRINTF_LONG_LONG_INT;
    else if (length == 'h')
      return PRINTF_SHORT_INT;
    else if (length == 'H')
      return PRINTF_SIGNED_CHAR;
    else {
      mdlerror_fmt(parse_state, "Format string segment '%s' uses unsupported "
                                "length specifier '%c' for integer.",
                   fmt_orig, length);
      return PRINTF_INVALID;
    }

  case 'o':
  case 'u':
  case 'x':
  case 'X':
    if (length == '\0')
      return PRINTF_U_INT;
    else if (length == 'l')
      return PRINTF_U_LONG_INT;
    else if (length == 'L')
      return PRINTF_U_LONG_LONG_INT;
    else if (length == 'h')
      return PRINTF_U_SHORT_INT;
    else if (length == 'H')
      return PRINTF_UNSIGNED_CHAR;
    else {
      mdlerror_fmt(parse_state, "Format string segment '%s' uses unsupported "
                                "length specifier '%c' for unsigned integer.",
                   fmt_orig, length);
      return PRINTF_INVALID;
    }

  case 'e':
  case 'E':
  case 'f':
  case 'F':
  case 'g':
  case 'G':
  case 'a':
  case 'A':
    if (length == '\0')
      return PRINTF_DOUBLE;
    else if (length == 'l')
      return PRINTF_DOUBLE;
    else {
      mdlerror_fmt(parse_state, "Format string segment '%s' uses unsupported "
                                "length specifier '%c' for double-precision "
                                "floating point.",
                   fmt_orig, length);
      return PRINTF_INVALID;
    }

  case 'c':
    if (length == '\0')
      return PRINTF_SIGNED_CHAR;
    else {
      mdlerror_fmt(parse_state, "Format string segment '%s' uses unsupported "
                                "length specifier '%c' for char.",
                   fmt_orig, length);
      return PRINTF_INVALID;
    }

  case 's':
    if (length == '\0')
      return PRINTF_STRING;
    else {
      mdlerror_fmt(parse_state, "Format string segment '%s' uses unsupported "
                                "length specifier '%c' for string.",
                   fmt_orig, length);
      return PRINTF_INVALID;
    }

  case 'C':
  case 'S':
  case 'p':
  case 'm':
  default:
    mdlerror_fmt(parse_state, "Format string segment '%s' uses unsupported "
                              "conversion specifier '%c'.",
                 fmt_orig, *fmt);
    return PRINTF_INVALID;
  }

  return PRINTF_INVALID;
}

/*************************************************************************
 my_sprintf_segment:
    fprintf a segment of a format string, accounting for data type conversions.

 In:  parse_state: parse arguments
      fmt_seg: format string segment
      argpp: pointer to args list (will be advanced if we consume arguments)
 Out: formatted segment on success, NULL on failure; *argpp is updated.
*************************************************************************/
static char *my_sprintf_segment(struct mdlparse_vars *parse_state,
                                char *fmt_seg, struct arg **argpp) {
  int num_asterisks = 0;
  struct arg *argp = *argpp;

  int spec_type =
      get_printf_conversion_specifier(parse_state, fmt_seg, &num_asterisks);
  if (spec_type == PRINTF_INVALID)
    return NULL;

  /* Pull out arguments for asterisks */
  int ast1 = 0, ast2 = 0;
  switch (num_asterisks) {
  case 0:
    break;

  case 2:
    if (argp->arg_type == DBL)
      ast2 = (int)(*(double *)(argp->arg_value) + EPS_C);
    else {
      mdlerror_fmt(parse_state, "Argument to be consumed by '*' in conversion "
                                "specification segment '%s' is not numeric",
                   fmt_seg);
      return NULL;
    }
    if ((argp = argp->next) == NULL) {
      mdlerror_fmt(
          parse_state,
          "Too few arguments for conversion specification segment '%s'",
          fmt_seg);
      return NULL;
    }
  /* FALL THROUGH */

  case 1:
    if (argp->arg_type == DBL)
      ast1 = (int)(*(double *)(argp->arg_value) + EPS_C);
    else {
      mdlerror_fmt(parse_state, "Argument to be consumed by '*' in conversion "
                                "specification segment '%s' is not numeric",
                   fmt_seg);
      return NULL;
    }
    if ((argp = argp->next) == NULL) {
      mdlerror_fmt(
          parse_state,
          "Too few arguments for conversion specification segment '%s'",
          fmt_seg);
      return NULL;
    }
    break;

  default:
    mdlerror_fmt(
        parse_state,
        "Invalid asterisk count in conversion specification segment '%s'",
        fmt_seg);
    return NULL;
  }

  /* Type check */
  if (argp->arg_type == STR) {
    if (spec_type != PRINTF_STRING && spec_type != PRINTF_SIGNED_CHAR &&
        spec_type != PRINTF_UNSIGNED_CHAR) {
      mdlerror_fmt(parse_state, "Argument in conversion specification segment "
                                "'%s' is a string, but a numeric value is "
                                "required",
                   fmt_seg);
      return NULL;
    }
  } else {
    if (spec_type == PRINTF_STRING) {
      mdlerror_fmt(parse_state, "Argument in conversion specification segment "
                                "'%s' is numeric, but a string value is "
                                "required",
                   fmt_seg);
      return NULL;
    }
  }

  /* Now, convert! */
  int arg_value;
  char *formatted = NULL;
  switch (spec_type) {

  case PRINTF_STRING:
    switch (num_asterisks) {
    case 0:
      formatted = alloc_sprintf(fmt_seg, (char *)argp->arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(fmt_seg, ast1, (char *)argp->arg_value);
      break;
    case 2:
      formatted = alloc_sprintf(fmt_seg, ast2, ast1, (char *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_SIGNED_CHAR:
  case PRINTF_UNSIGNED_CHAR:
    if (argp->arg_type == STR) {
      if (strlen((char *)argp->arg_value) > 1) {
        mdlerror_fmt(parse_state, "Argument in conversion specification "
                                  "segment '%s' has too many characters",
                     fmt_seg);
        return NULL;
      }

      if (spec_type == PRINTF_SIGNED_CHAR)
        arg_value = *(signed char *)argp->arg_value;
      else
        arg_value = *(unsigned char *)argp->arg_value;
    } else {
      if (spec_type == PRINTF_SIGNED_CHAR)
        arg_value = (signed char)*(double *)argp->arg_value;
      else
        arg_value = (unsigned char)*(double *)argp->arg_value;
    }

    switch (num_asterisks) {
    case 0:
      formatted = alloc_sprintf(fmt_seg, arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(fmt_seg, ast1, arg_value);
      break;
    case 2:
      formatted = alloc_sprintf(fmt_seg, ast2, ast1, arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_DOUBLE:
    switch (num_asterisks) {
    case 0:
      formatted = alloc_sprintf(fmt_seg, *(double *)argp->arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(fmt_seg, ast1, *(double *)argp->arg_value);
      break;
    case 2:
      formatted =
          alloc_sprintf(fmt_seg, ast2, ast1, *(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_INT:
    switch (num_asterisks) {
    case 0:
      formatted = alloc_sprintf(fmt_seg, (int)*(double *)argp->arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(fmt_seg, ast1, (int)*(double *)argp->arg_value);
      break;
    case 2:
      formatted =
          alloc_sprintf(fmt_seg, ast2, ast1, (int)*(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_LONG_INT:
    switch (num_asterisks) {
    case 0:
      formatted = alloc_sprintf(fmt_seg, (long int)*(double *)argp->arg_value);
      break;
    case 1:
      formatted =
          alloc_sprintf(fmt_seg, ast1, (long int)*(double *)argp->arg_value);
      break;
    case 2:
      formatted = alloc_sprintf(fmt_seg, ast2, ast1,
                                (long int)*(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_LONG_LONG_INT:
    switch (num_asterisks) {
    case 0:
      formatted =
          alloc_sprintf(fmt_seg, (long long int)*(double *)argp->arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(fmt_seg, ast1,
                                (long long int)*(double *)argp->arg_value);
      break;
    case 2:
      formatted = alloc_sprintf(fmt_seg, ast2, ast1,
                                (long long int)*(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_SHORT_INT:
    switch (num_asterisks) {
    case 0:
      formatted = alloc_sprintf(fmt_seg, (short int)*(double *)argp->arg_value);
      break;
    case 1:
      formatted =
          alloc_sprintf(fmt_seg, ast1, (short int)*(double *)argp->arg_value);
      break;
    case 2:
      formatted = alloc_sprintf(fmt_seg, ast2, ast1,
                                (short int)*(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_U_INT:
    switch (num_asterisks) {
    case 0:
      formatted =
          alloc_sprintf(fmt_seg, (unsigned int)*(double *)argp->arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(fmt_seg, ast1,
                                (unsigned int)*(double *)argp->arg_value);
      break;
    case 2:
      formatted = alloc_sprintf(fmt_seg, ast2, ast1,
                                (unsigned int)*(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_U_LONG_INT:
    switch (num_asterisks) {
    case 0:
      formatted =
          alloc_sprintf(fmt_seg, (unsigned long int)*(double *)argp->arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(fmt_seg, ast1,
                                (unsigned long int)*(double *)argp->arg_value);
      break;
    case 2:
      formatted = alloc_sprintf(fmt_seg, ast2, ast1,
                                (unsigned long int)*(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_U_LONG_LONG_INT:
    switch (num_asterisks) {
    case 0:
      formatted = alloc_sprintf(
          fmt_seg, (unsigned long long int)*(double *)argp->arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(
          fmt_seg, ast1, (unsigned long long int)*(double *)argp->arg_value);
      break;
    case 2:
      formatted =
          alloc_sprintf(fmt_seg, ast2, ast1,
                        (unsigned long long int)*(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  case PRINTF_U_SHORT_INT:
    switch (num_asterisks) {
    case 0:
      formatted = alloc_sprintf(fmt_seg,
                                (unsigned short int)*(double *)argp->arg_value);
      break;
    case 1:
      formatted = alloc_sprintf(fmt_seg, ast1,
                                (unsigned short int)*(double *)argp->arg_value);
      break;
    case 2:
      formatted = alloc_sprintf(fmt_seg, ast2, ast1,
                                (unsigned short int)*(double *)argp->arg_value);
      break;
    default:
      UNHANDLED_CASE(num_asterisks);
    }
    break;

  default:
    UNHANDLED_CASE(spec_type);
  }

  if (formatted == NULL)
    mcell_allocfailed("Failed to format a string for an MDL "
                      "printf/fprintf/spritnf statement.");
  else
    *argpp = argp->next;
  return formatted;
}

/*************************************************************************
 Temporary structure used by sprintf implementation.
*************************************************************************/
struct sprintf_output_list {
  struct sprintf_output_list *next; /* Next item in list */
  size_t len;                       /* Length of this item */
  char *segment;                    /* Data for this item */
};

/*************************************************************************
 my_sprintf:
    sprintf-like formatting of MDL arguments.

 In:  parse_state: parser state
      format: string to expand
      argp: argument list
 Out: 0 on success, 1 on failure
*************************************************************************/
static char *my_sprintf(struct mdlparse_vars *parse_state, char *format,
                        struct arg *argp) {
  char *this_start, *this_end;
  char const *const format_orig = format;

  struct sprintf_output_list head, *tail;
  head.segment = NULL;
  head.next = NULL;
  head.len = 0u;
  tail = &head;

  /* Find the start of the first format code */
  this_start = strchr(format, '%');
  while (this_start != NULL && this_start[1] == '%')
    this_start = strchr(this_start + 2, '%');

  /* If we have text before the first conversion specification, copy it. */
  if (this_start != NULL) {
    tail = tail->next = CHECKED_MALLOC_STRUCT(struct sprintf_output_list,
                                              "SPRINTF output list segment");
    if (tail == NULL)
      goto failure;
    memset(tail, 0, sizeof(*tail));
    *this_start = '\0';
    tail->segment = mdl_strdup(format);
    *this_start = '%';
    if (tail->segment == NULL)
      goto failure;
    tail->len = strlen(tail->segment);
    head.len += tail->len;
    format = this_start;
  }

  /* Process each format code */
  while (this_start != NULL) {
    /* If we've run out of arguments... */
    if (argp == NULL) {
      mdlerror_fmt(parse_state,
                   "Insufficient arguments for printf-style format string '%s'",
                   format_orig);
      goto failure;
    }

    /* Find the start of the first format code */
    this_end = strchr(this_start + 1, '%');
    while (this_end != NULL && this_end[1] == '%')
      this_end = strchr(this_end + 2, '%');

    /* If we have only a single format code, do the rest of the string */
    if (this_end == NULL) {
      tail = tail->next = CHECKED_MALLOC_STRUCT(struct sprintf_output_list,
                                                "SPRINTF output list segment");
      if (tail == NULL)
        goto failure;
      memset(tail, 0, sizeof(*tail));

      tail->segment = my_sprintf_segment(parse_state, format, &argp);
      if (tail->segment == NULL)
        goto failure;
      tail->len = strlen(tail->segment);
      head.len += tail->len;
      format = NULL;
      break;
    }

    /* Otherwise, print the entire string up to this point. */
    else {
      /* Print this segment */
      *this_end = '\0';
      tail = tail->next = CHECKED_MALLOC_STRUCT(struct sprintf_output_list,
                                                "SPRINTF output list segment");
      if (tail == NULL)
        goto failure;
      memset(tail, 0, sizeof(*tail));

      tail->segment = my_sprintf_segment(parse_state, format, &argp);
      if (tail->segment == NULL)
        goto failure;
      tail->len = strlen(tail->segment);
      head.len += tail->len;
      *this_end = '%';

      /* Advance to the next segment of the format string */
      format = this_end;
      this_start = strchr(format, '%');
      while (this_start != NULL && this_start[1] == '%')
        this_start = strchr(this_start + 2, '%');
    }
  }

  /* Now, build the string! */
  int total_len = head.len + 1;
  if (format != NULL)
    total_len += strlen(format);
  char *buffer = CHECKED_MALLOC_ARRAY(char, total_len, "SPRINTF data");
  char *pos = buffer;
  if (buffer == NULL)
    goto failure;
  while (head.next != NULL) {
    struct sprintf_output_list *l = head.next;
    memcpy(pos, l->segment, l->len);
    pos += l->len;
    free(l->segment);
    head.next = l->next;
    free(l);
  }
  if (format != NULL)
    strcpy(pos, format);
  else
    *pos = '\0';

  return buffer;

failure:
  while (head.next != NULL) {
    struct sprintf_output_list *l = head.next;
    if (l->segment != NULL)
      free(l->segment);
    head.next = l->next;
    free(l);
  }
  return NULL;
}

/*************************************************************************
 my_fprintf:
    fprintf-like formatting of MDL arguments.

 In:  parse_state: parser state
      outfile: file stream to which to write
      format: string to expand
      argp: argument list
 Out: 0 on success, 1 on failure
*************************************************************************/
static int my_fprintf(struct mdlparse_vars *parse_state, FILE *outfile,
                      char *format, struct arg *argp) {
  char *str = my_sprintf(parse_state, format, argp);
  if (str == NULL)
    return 1;

  if (fputs(str, outfile) < 0) {
    free(str);
    return 1;
  }

  free(str);
  return 0;
}

/*************************************************************************
 defang_format_string:
    Remove any potentially dangerous characters from a format string so it can
    be safely used in error reporting.

 In:  fmt: string to defang
 Out: None.  fmt string is updated
*************************************************************************/
static void defang_format_string(char *fmt) {
  while (*fmt != '\0') {
    if (iscntrl(*fmt))
      *fmt = '.';
    else if (!isascii(*fmt))
      *fmt = '.';
    ++fmt;
  }
}

/*************************************************************************
 mdl_printf:
    printf-like formatting of MDL arguments.  Prints to the defined log_file.

 In:  parse_state: parser state
      fmt: string to expand
      arg_head: argument list
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_printf(struct mdlparse_vars *parse_state, char *fmt,
               struct arg *arg_head) {
  if (my_fprintf(parse_state, mcell_get_log_file(), fmt, arg_head)) {
    mdl_free_printf_arg_list(arg_head);
    defang_format_string(fmt);
    mdlerror_fmt(parse_state, "Could not print to logfile: %s", fmt);
    free(fmt);
    return 1;
  }
  mdl_free_printf_arg_list(arg_head);
  free(fmt);
  return 0;
}

/*************************************************************************
 mdl_fprintf:
    fprintf-like formatting of MDL arguments.

 In:  parse_state: parser state
      filep: file stream to receive output
      fmt: string to expand
      arg_head: argument list
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_fprintf(struct mdlparse_vars *parse_state, struct file_stream *filep,
                char *fmt, struct arg *arg_head) {
  if (my_fprintf(parse_state, filep->stream, fmt, arg_head)) {
    mdl_free_printf_arg_list(arg_head);
    free(fmt);
    mdlerror_fmt(parse_state, "Could not print to file: %s", filep->name);
    return 1;
  }
  mdl_free_printf_arg_list(arg_head);
  free(fmt);
  return 0;
}

/*************************************************************************
 mdl_string_format:
    Expression-friendly sprintf-like formatting of MDL arguments.

 In:  parse_state: parser state
      fmt: string to expand
      arg_head: argument list
 Out: string on success, NULL on failure
*************************************************************************/
char *mdl_string_format(struct mdlparse_vars *parse_state, char *fmt,
                        struct arg *arg_head) {
  char *str = my_sprintf(parse_state, fmt, arg_head);
  if (str == NULL) {
    mdl_free_printf_arg_list(arg_head);
    free(fmt);
    mdlerror_fmt(parse_state, "Could not format string\n");
    return NULL;
  }
  free(fmt);
  mdl_free_printf_arg_list(arg_head);

  return str;
}

/*************************************************************************
 mdl_sprintf:
    sprintf-like formatting of MDL arguments.

 In:  parse_state: parser state
      assign_var: variable to receive formatted string
      fmt: string to expand
      arg_head: argument list
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_sprintf(struct mdlparse_vars *parse_state, struct sym_entry *assign_var,
                char *fmt, struct arg *arg_head) {
  char *str = my_sprintf(parse_state, fmt, arg_head);
  if (str == NULL) {
    mdl_free_printf_arg_list(arg_head);
    free(fmt);
    mdlerror_fmt(parse_state, "Could not sprintf to variable: %s",
                 assign_var->name);
    return 1;
  }
  free(fmt);
  mdl_free_printf_arg_list(arg_head);

  /* If the symbol had a value, try to free it */
  if (assign_var->value && mdl_free_variable_value(parse_state, assign_var)) {
    free(str);
    return 1;
  }

  /* Store new value */
  assign_var->sym_type = STR;
  assign_var->value = str;

  /* Log the updated value */
  no_printf("\n%s is equal to: %s\n", assign_var->name, str);
  return 0;
}

/*************************************************************************
 mdl_fprint_time:
    strtime-like formatting of current time.

 In:  parse_state: parser state
      filep_sym: file stream to receive output
      fmt: string to expand
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_fprint_time(struct mdlparse_vars *parse_state,
                    struct sym_entry *filep_sym, char *fmt) {
  char time_str[128];
  time_t the_time;
  struct file_stream *filep = (struct file_stream *)filep_sym->value;

  the_time = time(NULL);
  strftime(time_str, 128, fmt, localtime(&the_time));
  free(fmt);
  if (fprintf(filep->stream, "%s", time_str) == EOF) {
    mdlerror_fmt(parse_state, "Could not print to file: %s", filep_sym->name);
    return 1;
  }

  return 0;
}

/*************************************************************************
 mdl_print_time:
    strtime-like formatting of current time.  Prints to logfile.

 In:  parse_state: parser state
      filep: file stream to receive output
      fmt: string to expand
 Out: 0 on success, 1 on failure
*************************************************************************/
void mdl_print_time(struct mdlparse_vars *parse_state, char *fmt) {
  char time_str[128];
  time_t the_time = time(NULL);
  strftime(time_str, 128, fmt, localtime(&the_time));
  free(fmt);
  if (parse_state->vol->procnum == 0)
    fprintf(mcell_get_log_file(), "%s", time_str);
}

/*************************************************************************
 mdl_generate_range:
    Generate a num_expr_list containing the numeric values from start to end,
    incrementing by step.

 In:  parse_state:  parser state
      list:  destination to receive list of values
      start: start of range
      end:   end of range
      step:  increment
 Out: 0 on success, 1 on failure.  On success, list is filled in.
*************************************************************************/
int mdl_generate_range(struct mdlparse_vars *parse_state,
                       struct num_expr_list_head *list, double start,
                       double end, double step) {
  if (step == 0.0) {
    mdlerror(parse_state, "A step size of 0 was requested, which would "
                          "generate an infinite list.");
    return 1;
  }

  /* Above a certain point, we probably don't want this.  If we need arrays
   * this large, we need to reengineer this. */
  if (fabs((end - start) / step) > 100000000.) {
    mdlerror(parse_state, "A range was requested that encompasses too many "
                          "values (maximum 100,000,000)");
    return 1;
  }

  mcell_generate_range(list, start, end, step);
  return 0;
}

/*************************************************************************
 mdl_add_range_value:
    Add a value to a numeric list.

 In:  lh:   list to receive value
      value: value for list
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_add_range_value(struct num_expr_list_head *lh, double value) {
  if (lh->value_head == NULL)
    return mcell_generate_range_singleton(lh, value);

  struct num_expr_list *nel =
      CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric array");
  if (nel == NULL)
    return 1;
  nel->next = NULL;
  nel->value = value;
  lh->value_tail = lh->value_tail->next = nel;
  ++lh->value_count;
  return 0;
}

#ifdef DEBUG
/*************************************************************************
 mdl_debug_dump_array:
    Display a human-readable representation of the array to stdout.

 In:  nel:  the list to free
 Out: displayed to stdout
*************************************************************************/
void mdl_debug_dump_array(struct num_expr_list *nel) {
  struct num_expr_list *elp;
  no_printf("\nArray expression: [");
  for (elp = nel; elp != NULL; elp = elp->next) {
    no_printf("%lf", elp->value);
    if (elp->next != NULL)
      no_printf(",");
  }
  no_printf("]\n");
}
#endif

/*************************************************************************
 num_expr_list_to_array:
    Convert a numeric list into an array.  The list count will be stored into
    the count pointer, if it is not NULL.

 In:  lh: list head
      count: location to receive count of items in list
 Out: an array of all doubles in the list, or NULL if allocation fails or if
      the list is empty.  If count is non-zero and NULL is returned, allocation
      has failed.
*************************************************************************/
static double *num_expr_list_to_array(struct num_expr_list_head *lh,
                                      int *count) {

  double *arr = NULL, *ptr = NULL;
  struct num_expr_list *i;

  if (count != NULL)
    *count = lh->value_count;

  if (lh->value_count == 0)
    return NULL;

  arr = CHECKED_MALLOC_ARRAY(double, lh->value_count, "numeric array");
  if (arr != NULL) {
    for (i = (struct num_expr_list *)lh->value_head, ptr = arr; i != NULL;
         i = i->next)
      *ptr++ = i->value;
  }
  return arr;
}

/*************************************************************************
 mdl_point:
    Create a 3-D vector from a numeric array.

 In:  parse_state: parser state
      vals: values to put into vector
 Out: a 3-D vector, or NULL if an error occurs.
*************************************************************************/
struct vector3 *mdl_point(struct mdlparse_vars *parse_state,
                          struct num_expr_list_head *vals) {
  if (vals->value_count != 3) {
    mdlerror(parse_state, "Three dimensional value required");
    return NULL;
  }

  struct vector3 *vec = CHECKED_MALLOC_STRUCT(struct vector3, "3-D vector");
  if (!vec)
    return NULL;

  vec->x = vals->value_head->value;
  vec->y = vals->value_head->next->value;
  vec->z = vals->value_tail->value;
  if (!vals->shared)
    mcell_free_numeric_list(vals->value_head);
  return vec;
}

/*************************************************************************
 mdl_point_scalar:
    Create a 3-D vector equal to s*[1, 1, 1] for some scalar s.

 In:  val: scalar
 Out: a 3-D vector, or NULL if an error occurs.
*************************************************************************/
struct vector3 *mdl_point_scalar(double val) {

  struct vector3 *vec = CHECKED_MALLOC_STRUCT(struct vector3, "3-D vector");
  if (!vec)
    return NULL;

  vec->x = val;
  vec->y = val;
  vec->z = val;
  return vec;
}

/**************************************************************************
 mdl_free_variable_value:
    Free the value in a given variable entry.  Currently, arrays are not freed,
    as they are also not copied when variable assignments occur.  We would need
    to either reference count arrays or directly copy them on assignment in
    order to be able to free them.

 In: parse_state: the parse variables structure
     sym: the symbol whose value to free
    Out: 0 on success, 1 if the symbol is not a double, string, or array
**************************************************************************/
static int mdl_free_variable_value(struct mdlparse_vars *parse_state,
                                   struct sym_entry *sym) {
  switch (sym->sym_type) {
  case DBL:
  case STR:
    free(sym->value);
    break;

  case ARRAY:
    /* XXX: For the moment, we can't free these since they might be shared by
     * two different variables... */
    break;

  default:
    mdlerror_fmt(
        parse_state,
        "Internal error: Attempt to free variable value of incorrect type: %s",
        sym->name);
    return 1;
  }

  return 0;
}

/**************************************************************************
 mdl_get_or_create_variable:
    Get a named variable if it exists, or create it if it doesn't.  Caller is
    no longer responsible for deallocation of 'name'.

 In: parse_state: the parse variables structure
     name: the name of the variable
 Out: a symbol table entry, or NULL if we ran out of memory
**************************************************************************/
struct sym_entry *mdl_get_or_create_variable(struct mdlparse_vars *parse_state,
                                             char *name) {
  /* Attempt to fetch existing variable */
  struct sym_entry *st = NULL;
  if ((st = retrieve_sym(name, parse_state->vol->var_sym_table)) != NULL) {
    free(name);
    return st;
  }

  /* Create the variable */
  if ((st = store_sym(name, TMP, parse_state->vol->var_sym_table, NULL)) ==
      NULL)
    mcell_allocfailed("Failed to store a variable symbol in the symbol table.");

  free(name);
  return st;
}

/**************************************************************************
 mdl_assign_variable_double:
    Assign a "double" value to a variable, freeing any previous value.

 In: parse_state: the parse variables structure
     sym: the symbol whose value to free
     value: the value to assign
 Out: 0 on success, 1 if the symbol is not a double, string, or array, or if
      memory allocation fails
**************************************************************************/
int mdl_assign_variable_double(struct mdlparse_vars *parse_state,
                               struct sym_entry *sym, double value) {
  /* If the symbol had a value, try to free it */
  if (sym->value && mdl_free_variable_value(parse_state, sym))
    return 1;

  /* Allocate space for the new value */
  if ((sym->value = CHECKED_MALLOC_STRUCT(double, "numeric variable")) == NULL)
    return 1;

  /* Update the value */
  *(double *)sym->value = value;
  sym->sym_type = DBL;
  no_printf("\n%s is equal to: %f\n", sym->name, value);
  return 0;
}

/**************************************************************************
 mdl_assign_variable_string:
    Assign a string value to a variable, freeing any previous value.

 In: parse_state: the parse variables structure
     sym: the symbol whose value to free
     value: the value to assign
 Out: 0 on success, 1 if the symbol is not a double, string, or array, or if
      memory allocation fails.
**************************************************************************/
int mdl_assign_variable_string(struct mdlparse_vars *parse_state,
                               struct sym_entry *sym, char *value) {
  /* If the symbol had a value, try to free it */
  if (sym->value && mdl_free_variable_value(parse_state, sym))
    return 1;

  /* Allocate space for the new value */
  if ((sym->value = mdl_strdup(value)) == NULL)
    return 1;

  free(value);

  /* Update the value */
  sym->sym_type = STR;
  no_printf("\n%s is equal to: %s\n", sym->name, value);
  return 0;
}

/**************************************************************************
 mdl_assign_variable_array:
    Assign an array value to a variable, freeing any previous value.

 In: parse_state: the parse variables structure
     sym: the symbol whose value to free
     value: the value to assign
 Out: 0 on success, 1 if the symbol is not a double, string, or array
**************************************************************************/
int mdl_assign_variable_array(struct mdlparse_vars *parse_state,
                              struct sym_entry *sym,
                              struct num_expr_list *value) {
  /* If the symbol had a value, try to free it */
  if (sym->value && mdl_free_variable_value(parse_state, sym))
    return 1;

  sym->sym_type = ARRAY;
  sym->value = value;
#ifdef DEBUG
  struct num_expr_list *elp = (struct num_expr_list *)sym->value;
  no_printf("\n%s is equal to: [", sym->name);
  while (elp != NULL) {
    no_printf("%lf", elp->value);
    elp = elp->next;
    if (elp != NULL) {
      no_printf(",");
    }
  }
  no_printf("]\n");
#endif
  return 0;
}

/**************************************************************************
 mdl_assign_variable:
    Assign one variable value to another variable, freeing any previous value.

 In: parse_state: the parse variables structure
     sym: the symbol whose value to free
     value: the value to assign
 Out: 0 on success, 1 if the symbol is not a double, string, or array
**************************************************************************/
int mdl_assign_variable(struct mdlparse_vars *parse_state,
                        struct sym_entry *sym, struct sym_entry *value) {
  switch (value->sym_type) {
  case DBL:
    if (mdl_assign_variable_double(parse_state, sym, *(double *)value->value))
      return 1;
    break;

  case STR:
    if (mdl_assign_variable_string(parse_state, sym,
                                   (char *)value->value))
      return 1;
    break;

  case ARRAY:
    if (mdl_assign_variable_array(parse_state, sym,
                                  (struct num_expr_list *)value->value))
      return 1;
    break;

  default:
    UNHANDLED_CASE(value->sym_type);
  }

  return 0;
}

/*************************************************************************
 mdl_set_all_notifications:
    Set all notification levels to a particular value.  Corresponds to
    ALL_NOTIFICATIONS rule in grammar.

 In:  vol: the world
      notify_value: level to which to set notifications
 Out: notification levels are changed
*************************************************************************/
void mdl_set_all_notifications(struct volume *vol, byte notify_value) {
  vol->notify->progress_report = notify_value;
  vol->notify->diffusion_constants = notify_value;
  vol->notify->reaction_probabilities = notify_value;
  vol->notify->time_varying_reactions = notify_value;
  vol->notify->partition_location = notify_value;
  vol->notify->box_triangulation = notify_value;
  if (vol->log_freq == ULONG_MAX)
    vol->notify->iteration_report = notify_value;
  vol->notify->throughput_report = notify_value;
  vol->notify->checkpoint_report = notify_value;
  vol->notify->release_events = notify_value;
  vol->notify->file_writes = notify_value;
  vol->notify->final_summary = notify_value;
  vol->notify->molecule_collision_report = notify_value;
}

/*************************************************************************
 mdl_set_iteration_report_freq:
    Set the iteration report frequency.

 In:  parse_state: parser state
      interval: report frequency to request
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_iteration_report_freq(struct mdlparse_vars *parse_state,
                                  long long interval) {
  /* Only if not set on command line */
  if (parse_state->vol->log_freq == ULONG_MAX) {
    if (interval < 1) {
      mdlerror(parse_state,
               "Invalid iteration number reporting interval: use value >= 1");
      return 1;
    } else {
      parse_state->vol->notify->custom_iteration_value = interval;
      parse_state->vol->notify->iteration_report = NOTIFY_FULL;
    }
  }
  return 0;
}

/*************************************************************************
 mdl_set_all_warnings:
    Set all warning levels to a particular value.  Corresponds to
    ALL_WARNINGS rule in grammar.

 In:  vol: the world
      warning_value: level to which to set warnings
 Out: warning levels are changed
*************************************************************************/
void mdl_set_all_warnings(struct volume *vol, byte warning_level) {
  vol->notify->neg_diffusion = warning_level;
  vol->notify->neg_reaction = warning_level;
  vol->notify->high_reaction_prob = warning_level;
  vol->notify->close_partitions = warning_level;
  vol->notify->degenerate_polys = warning_level;
  vol->notify->overwritten_file = warning_level;
  vol->notify->mol_placement_failure = warning_level;

  if (warning_level == WARN_ERROR)
    warning_level = WARN_WARN;
  vol->notify->short_lifetime = warning_level;
  vol->notify->missed_reactions = warning_level;
  vol->notify->missed_surf_orient = warning_level;
  vol->notify->useless_vol_orient = warning_level;
  vol->notify->invalid_output_step_time = warning_level;
  vol->notify->large_molecular_displacement = warning_level;
  vol->notify->add_remove_mesh_warning = warning_level;
}

/*************************************************************************
 mdl_set_lifetime_warning_threshold:
    Set the lifetime warning threshold.

 In:  parse_state: parser state
      lifetime: threshold to set
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_lifetime_warning_threshold(struct mdlparse_vars *parse_state,
                                       long long lifetime) {
  if (lifetime < 0) {
    mdlerror(
        parse_state,
        "Molecule lifetimes are measured in iterations and cannot be negative");
    return 1;
  }
  parse_state->vol->notify->short_lifetime_value = lifetime;
  return 0;
}

/*************************************************************************
 mdl_set_missed_reaction_warning_threshold:
    Set the missed reaction warning threshold.

 In:  parse_state: parser state
      rxfrac: threshold to set
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_missed_reaction_warning_threshold(struct mdlparse_vars *parse_state,
                                              double rxfrac) {
  if (rxfrac < 0.0 || rxfrac > 1.0) {
    mdlerror(
        parse_state,
        "Values for fraction of reactions missed should be between 0 and 1");
    return 1;
  }
  parse_state->vol->notify->missed_reaction_value = rxfrac;
  return 0;
}

/*************************************************************************
 mdl_set_time_step:
    Set the global timestep for the simulation.

 In:  parse_state: parser state
      step: timestep to set
 Out: 0 on success, 1 on failure; global timestep is updated
*************************************************************************/
int mdl_set_time_step(struct mdlparse_vars *parse_state, double step) {

  int error_code = mcell_set_time_step(parse_state->vol, step);
  if (error_code == 2) {
    mdlerror_fmt(
        parse_state,
        "Time step of %.15g requested; time step must be a positive value",
        step);
    return 1;
  } else if (error_code == 3) {
    mdlerror_fmt(parse_state, "Time step of %.15g requested, but the time step "
                              "was already set to %.15g",
                 step, parse_state->vol->time_unit);
    return 1;
  }
  no_printf("Time unit = %g\n", parse_state->vol->time_unit);
  return 0;
}

/*************************************************************************
 mdl_set_max_time_step:
    Set the maximum timestep for the simulation.

 In:  parse_state: parser state
      step: maximum timestep to set
 Out: 0 on success, 1 on failure; on success, maximum timestep is updated
*************************************************************************/
int mdl_set_max_time_step(struct mdlparse_vars *parse_state, double step) {
  if (step <= 0) {
    mdlerror_fmt(parse_state, "Maximum time step of %.15g requested; maximum "
                              "time step must be a positive value",
                 step);
    return 1;
  }

  if (distinguishable(parse_state->vol->time_step_max, 0, EPS_C)) {
    mdlerror_fmt(parse_state, "Maximum time step of %.15g requested, but the "
                              "maximum time step was already set to %.15g",
                 step, parse_state->vol->time_step_max);
    return 1;
  }

  parse_state->vol->time_step_max = step;
  no_printf("Maximum time step = %g\n", parse_state->vol->time_step_max);
  return 0;
}

/*************************************************************************
 mdl_set_space_step:
    Set the global space step for the simulation.

 In:  parse_state: parser state
      step: global space step to set
 Out: 0 on success, 1 on failure; on success, global space step is updated
*************************************************************************/
int mdl_set_space_step(struct mdlparse_vars *parse_state, double step) {
  if (step <= 0) {
    mdlerror_fmt(
        parse_state,
        "Space step of %.15g requested; space step must be a positive value",
        step);
    return 1;
  }

  if (distinguishable(parse_state->vol->space_step, 0, EPS_C)) {
    mdlerror_fmt(parse_state, "Space step of %.15g requested, but the space "
                              "step was already set to %.15g",
                 step, parse_state->vol->space_step *
                           parse_state->vol->r_length_unit);
    return 1;
  }

  parse_state->vol->space_step = step;
  no_printf("Space step = %g\n", parse_state->vol->space_step);

  /* Use internal units */
  parse_state->vol->space_step *= parse_state->vol->r_length_unit;
  return 0;
}

/*************************************************************************
 mdl_set_num_iterations:
    Set the number of iterations for the simulation.

 In:  parse_state: parser state
      numiters: number of iterations to run
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_num_iterations(struct mdlparse_vars *parse_state,
                           long long numiters) {
  /* If the iteration count was not overridden on the command-line... */
  if (parse_state->vol->iterations == INT_MIN) {
    parse_state->vol->iterations = numiters;
    if (parse_state->vol->iterations < 0) {
      mdlerror(parse_state, "ITERATIONS value is negative");
      return 1;
    }
  }
  no_printf("Iterations = %lld\n", parse_state->vol->iterations);
  return 0;
}

/*************************************************************************
 mdl_set_num_radial_directions:
    Set the number of radial directions.

 In:  parse_state: parser state
      numdirs: number of radial directions to use
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_num_radial_directions(struct mdlparse_vars *parse_state,
                                  int numdirs) {
  parse_state->vol->radial_directions = numdirs;
  if (parse_state->vol->radial_directions <= 0) {
    mdlerror(parse_state, "RADIAL_DIRECTIONS must be a positive number.");
    return 1;
  }

  parse_state->vol->num_directions = 0;
  if (parse_state->vol->d_step != NULL)
    free(parse_state->vol->d_step);
  if ((parse_state->vol->d_step =
           init_d_step(parse_state->vol->radial_directions,
                       &parse_state->vol->num_directions)) == NULL) {
    mcell_allocfailed("Failed to allocate the diffusion directions table.");
    /*return 1;*/
  }

  /* Mask must contain at least every direction */
  /* Num directions, rounded up to the nearest 2^n - 1 */
  parse_state->vol->directions_mask = parse_state->vol->num_directions;
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 1);
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 2);
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 4);
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 8);
  parse_state->vol->directions_mask |=
      (parse_state->vol->directions_mask >> 16);
  if (parse_state->vol->directions_mask > (1 << 18)) {
    mdlerror(parse_state, "Too many RADIAL_DIRECTIONS requested (max 131072).");
    return 1;
  }

  no_printf("desired radial directions = %d\n",
            parse_state->vol->radial_directions);
  no_printf("actual radial directions = %d\n",
            parse_state->vol->num_directions);
  return 0;
}

/*************************************************************************
 mdl_set_num_radial_subdivisions:
    Set the number of radial subdivisions.

 In:  parse_state: parser state
      numdivs: number of radial subdivisions to use
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_num_radial_subdivisions(struct mdlparse_vars *parse_state,
                                    int numdivs) {
  parse_state->vol->radial_subdivisions = numdivs;
  if (parse_state->vol->radial_subdivisions <= 0) {
    mdlerror(parse_state, "RADIAL_SUBDIVISIONS must be a positive number.");
    return 1;
  }

  /* 'x & (x-1)' clears the lowest-order bit of x.  thus, only for 0 or a  */
  /* power of 2 will the expression be zero.  we've eliminated the case of */
  /* zero in the previous test.  (this condition is required so that we    */
  /* use 'x & (RSD - 1)' as an optimized form of 'x % RSD'.)               */
  if ((parse_state->vol->radial_subdivisions &
       (parse_state->vol->radial_subdivisions - 1)) != 0) {
    mdlerror(parse_state, "Radial subdivisions must be a power of two");
    return 1;
  }

  if (parse_state->vol->r_step != NULL)
    free(parse_state->vol->r_step);
  if (parse_state->vol->r_step_surface != NULL)
    free(parse_state->vol->r_step_surface);
  if (parse_state->vol->r_step_release != NULL)
    free(parse_state->vol->r_step_release);

  parse_state->vol->r_step = init_r_step(parse_state->vol->radial_subdivisions);
  parse_state->vol->r_step_surface =
      init_r_step_surface(parse_state->vol->radial_subdivisions);
  parse_state->vol->r_step_release =
      init_r_step_3d_release(parse_state->vol->radial_subdivisions);

  if (parse_state->vol->r_step == NULL ||
      parse_state->vol->r_step_surface == NULL ||
      parse_state->vol->r_step_release == NULL) {
    mcell_allocfailed(
        "Failed to allocate the diffusion radial subdivisions table.");
    /*return 1;*/
  }

  no_printf("radial subdivisions = %d\n",
            parse_state->vol->radial_subdivisions);
  return 0;
}

/*************************************************************************
 mdl_set_interaction_radius:
    Set the interaction radius.

 In:  parse_state: parser state
      interaction_radius
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_interaction_radius(struct mdlparse_vars *parse_state,
                               double interaction_radius) {
  parse_state->vol->rx_radius_3d = interaction_radius;
  if (parse_state->vol->rx_radius_3d <= 0) {
    mdlerror(parse_state, "INTERACTION_RADIUS must be a positive number.");
    return 1;
  }

  no_printf("interaction radius = %f\n", parse_state->vol->rx_radius_3d);
  return 0;
}

/*************************************************************************
 mdl_set_grid_density:
    Set the surface grid density.

 In:  parse_state: parser state
      density: global surface grid density for simulation
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_grid_density(struct mdlparse_vars *parse_state, double density) {
  if (density <= 0) {
    mdlerror_fmt(parse_state, "EFFECTOR_GRID_DENSITY must be greater than 0.0 "
                              "(value provided was %lg)",
                 density);
    return 1;
  }
  parse_state->vol->grid_density = density;
  no_printf("Max density = %f\n", parse_state->vol->grid_density);

  parse_state->vol->space_step *= parse_state->vol->length_unit;
  parse_state->vol->r_length_unit = sqrt(parse_state->vol->grid_density);
  parse_state->vol->length_unit = 1.0 / parse_state->vol->r_length_unit;
  parse_state->vol->space_step *= parse_state->vol->r_length_unit;

  no_printf("Length unit = %f\n", parse_state->vol->length_unit);
  return 0;
}

/*************************************************************************
 schedule_async_checkpoint:
    Schedule an asynchronous checkpoint.

 In:  parse_state: parser state
      dur: length of real (wall-clock) time until we should checkpoint
 Out: 0 on success, 1 on failure; signal handler is installed and an alarm is
      scheduled.
*************************************************************************/
static int schedule_async_checkpoint(struct mdlparse_vars *parse_state,
                                     unsigned int dur) {
#ifdef _WIN32
  set_alarm_handler(&chkpt_signal_handler);
#else
  struct sigaction sa, saPrev;
  sa.sa_sigaction = NULL;
  sa.sa_handler = &chkpt_signal_handler;
  sa.sa_flags = SA_RESTART;
  sigfillset(&sa.sa_mask);

  if (sigaction(SIGALRM, &sa, &saPrev) != 0) {
    mdlerror(parse_state, "Failed to install ALRM signal handler");
    return 1;
  }
#endif
  alarm(dur);

  return 0;
}

/*************************************************************************
 mdl_set_realtime_checkpoint:
    Schedule an asynchronous checkpoint.

 In:  parse_state: parser state
      duration: length of real (wall-clock) time until we should checkpoint
      cont_after_cp: 1 if we should continue after checkpoint, 0 if we should
                     exit
 Out: 0 on success, 1 on failure; signal handler is installed and an alarm is
      scheduled.
*************************************************************************/
int mdl_set_realtime_checkpoint(struct mdlparse_vars *parse_state,
                                long duration, int cont_after_cp) {
  time_t t;
  time(&t);

  if (duration <= 0) {
    mdlerror_fmt(parse_state, "Realtime checkpoint requested at %ld seconds, "
                              "but the checkpoint interval must be positive.",
                 duration);
    return 1;
  } else if (parse_state->vol->checkpoint_alarm_time != 0) {
    mdlerror_fmt(parse_state, "Realtime checkpoint requested at %ld seconds, "
                              "but a checkpoint is already scheduled for %u "
                              "seconds.",
                 duration, parse_state->vol->checkpoint_alarm_time);
    return 1;
  }

  parse_state->vol->checkpoint_alarm_time = (unsigned int)duration;
  parse_state->vol->continue_after_checkpoint = cont_after_cp;
  parse_state->vol->chkpt_flag = 1;
  if (t - parse_state->vol->begin_timestamp > 0) {
    if (duration <= t - parse_state->vol->begin_timestamp) {
      mdlerror_fmt(parse_state, "Checkpoint scheduled for %ld seconds, but %ld "
                                "seconds have already elapsed during parsing.  "
                                "Exiting.",
                   duration, (long)(t - parse_state->vol->begin_timestamp));
      return 1;
    }

    duration -= (t - parse_state->vol->begin_timestamp);
  }

  if (schedule_async_checkpoint(parse_state, (unsigned int)duration)) {
    mdlerror_fmt(parse_state, "Failed to schedule checkpoint for %ld seconds",
                 duration);
    return 1;
  }

  return 0;
}

/*************************************************************************
 mdl_set_checkpoint_infile:
    Set the input checkpoint file to use.

 In:  parse_state: parser state
      name: name of checkpoint file
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_checkpoint_infile(struct mdlparse_vars *parse_state, char *name) {
  FILE *file;
  if (parse_state->vol->chkpt_infile == NULL) {
    parse_state->vol->chkpt_infile = name;
    if ((file = fopen(parse_state->vol->chkpt_infile, "r")) == NULL) {
      parse_state->vol->chkpt_init = 1;
    } else {
      parse_state->vol->chkpt_init = 0;
      fclose(file);
    }
    parse_state->vol->chkpt_flag = 1;
  }
  return 0;
}

/*************************************************************************
 mdl_set_checkpoint_outfile:
    Set the output checkpoint file to use.

 In:  parse_state: parser state
      name: name of checkpoint file
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_checkpoint_outfile(struct mdlparse_vars *parse_state, char *name) {
  if (parse_state->vol->chkpt_outfile == NULL) {
    parse_state->vol->chkpt_outfile = name;
    parse_state->vol->chkpt_flag = 1;
  }
  return 0;
}

/*************************************************************************
 mdl_set_checkpoint_iterations:
    Set the number of iterations between checkpoints.

 In:  parse_state: parser state
      iters: number of iterations between checkpoints
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_set_checkpoint_interval(struct mdlparse_vars *parse_state,
                                long long iters, int continueAfterChkpt) {
  parse_state->vol->chkpt_iterations = iters;
  if (parse_state->vol->chkpt_iterations <= 0) {
    mdlerror(parse_state, "CHECKPOINT_ITERATIONS must be a positive integer");
    return 1;
  }
  parse_state->vol->chkpt_flag = 1;
  parse_state->vol->continue_after_checkpoint = continueAfterChkpt;
  return 0;
}

/*************************************************************************
 mdl_keep_checkpoint_files:
    Select if previous checkpoint files should be kept rather than be
    overwritten each time a new checkpoint starts. This option only
    affects checkpointed simulations run with the NOEXIT option.

 In:  parse_state: parser state
      keepFiles: boolean variables selecting if intermediate checkpointing
                 files should be kept or not
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_keep_checkpoint_files(struct mdlparse_vars *parse_state,
                              int keepFiles) {

  parse_state->vol->keep_chkpts = keepFiles;
  return 0;
}

/*************************************************************************
 mdl_make_new_object:
    Create a new object, adding it to the global symbol table.  the object must
    not be defined yet.

 In:  parse_state: parser state
      obj_name: fully qualified object name
 Out: the newly created object
*************************************************************************/
static struct object *mdl_make_new_object(struct mdlparse_vars *parse_state,
                                          char *obj_name) {

  int error_code = 0;
  struct object *obj_ptr = make_new_object(
      parse_state->vol->dg_parse,
      parse_state->vol->obj_sym_table,
      obj_name,
      &error_code);

  return obj_ptr;
}

/*************************************************************************
 mdl_start_object:
    Create a new object, adding it to the global symbol table.  the object must
    not be defined yet.  The qualified name of the object will be built by
    adding to the object_name_list, and the object is made the "current_object"
    in the mdl parser state.  Because of these side effects, it is vital to
    call mdl_finish_object at the end of the scope of the object created here.

 In:  parse_state: parser state
      obj_name: unqualified object name
 Out: the newly created object
*************************************************************************/
struct sym_entry *mdl_start_object(struct mdlparse_vars *parse_state,
                                   char *name) {
  // Create new fully qualified name.
  char *new_name;
  struct object_creation obj_creation;
  obj_creation.object_name_list = parse_state->object_name_list;
  obj_creation.object_name_list_end = parse_state->object_name_list_end;
  if ((new_name = push_object_name(&obj_creation, name)) == NULL) {
    mdlerror_fmt(parse_state, "Out of memory while creating object: %s", name);
    free(name);
    return NULL;
  }
  parse_state->object_name_list = obj_creation.object_name_list;
  parse_state->object_name_list_end = obj_creation.object_name_list_end;

  // Create the symbol, if it doesn't exist yet.
  int error_code = 0;
  struct object *obj_ptr = make_new_object(
      parse_state->vol->dg_parse,
      parse_state->vol->obj_sym_table,
      new_name,
      &error_code);
  if (error_code == 1) {
    mdlerror_fmt(parse_state,"Object '%s' is already defined", new_name);
  }
  else if (error_code == 2) {
    mdlerror_fmt(parse_state, "Out of memory while creating object: %s",
                 new_name);
  }


  struct sym_entry *sym_ptr = obj_ptr->sym;
  obj_ptr->last_name = name;
  no_printf("Creating new object: %s\n", new_name);

  // Set parent object, make this object "current".
  obj_ptr->parent = parse_state->current_object;
  parse_state->current_object = obj_ptr;

  return sym_ptr;
}

/*************************************************************************
 mdl_finish_object:
    "Finishes" a new object, undoing the state changes that occurred when the
    object was "started".  (This means popping the name off of the object name
    stack, and resetting current_object to its value before this object was
    defined.

 In:  parse_state: parser state
      obj_name: unqualified object name
 Out: the newly created object
*************************************************************************/
void mdl_finish_object(struct mdlparse_vars *parse_state) {
  struct object_creation obj_creation;
  obj_creation.object_name_list_end = parse_state->object_name_list_end;

  pop_object_name(&obj_creation);
  parse_state->object_name_list_end = obj_creation.object_name_list_end;
  parse_state->current_object = parse_state->current_object->parent;
}

/*************************************************************************
 object_list_singleton:
    Adds the first element to an empty object list.

 In:  head: object list head
      objp: object to add
 Out: none
*************************************************************************/
void mdl_object_list_singleton(struct object_list *head, struct object *objp) {
  objp->next = NULL;
  head->obj_tail = head->obj_head = objp;
}

/*************************************************************************
 add_object_to_list:
    Adds an element to an object list.

 In:  head: object list head
      objp: object to add
 Out: none
*************************************************************************/
void mdl_add_object_to_list(struct object_list *head, struct object *objp) {
  objp->next = NULL;
  head->obj_tail = head->obj_tail->next = objp;
}

/*************************************************************************
 SYMBOL_TYPE_DESCRIPTIONS:
    Human-readable symbol type descriptions, indexed by symbol type (MOL, OBJ,
    STR, etc.)
*************************************************************************/
static const char *SYMBOL_TYPE_DESCRIPTIONS[] = {
  "reaction",        "reaction pathway", "molecule",          "object",
  "release pattern", "region",           "numeric variable",  "string variable",
  "array variable",  "file stream",      "placeholder value", "trigger",
};

/*************************************************************************
 SYMBOL_TYPE_ARTICLES:
    Indefinite articles to go with human-readable symbol type descriptions,
    indexed by symbol type (MOL, OBJ, STR, etc.)
*************************************************************************/
static const char *SYMBOL_TYPE_ARTICLES[] = { "a", "a", "a",  "an", "a", "a",
                                              "a", "a", "an", "a",  "a", "a", };

/*************************************************************************
 mdl_symbol_type_name:
    Get a human-readable description of a symbol type, given its numerical
    code.

 In:  type: numeric type code
 Out: the symbol type name
*************************************************************************/
static char const *mdl_symbol_type_name(enum symbol_type_t type) {

  if (type <= 0 || type >= (int)COUNT_OF(SYMBOL_TYPE_DESCRIPTIONS)) {
    mcell_internal_error("Invalid symbol type '%d'", type);
    /*return "(unknown symbol type)";*/
  }

  return SYMBOL_TYPE_DESCRIPTIONS[type];
}

/*************************************************************************
 mdl_symbol_type_name_article:
    Get an indefinite article to precede the human-readable description of a
    symbol type, given its numerical code.

 In:  type: numeric type code
 Out: the article
*************************************************************************/
static char const *mdl_symbol_type_name_article(enum symbol_type_t type) {

  if (type <= 0 || type >= (int)COUNT_OF(SYMBOL_TYPE_ARTICLES)) {
    mcell_internal_error("Invalid symbol type '%d'", type);
    /*return "an";*/
  }

  return SYMBOL_TYPE_ARTICLES[type];
}

/*************************************************************************
 mdl_existing_symbol:
    Find an existing symbol or print an error message.  Also print an error
    message if the symbol is not of the right type.

 In:  parse_state: parser state
      name: name of symbol to find
      tab:  table to search
      type: type of symbol to find
 Out: returns the symbol, or NULL if none found
*************************************************************************/
static struct sym_entry *mdl_existing_symbol(struct mdlparse_vars *parse_state,
                                             char *name,
                                             struct sym_table_head *tab,
                                             int type) {
  struct sym_entry *symp = retrieve_sym(name, tab);
  if (symp == NULL)
    mdlerror_fmt(parse_state, "Undefined %s: %s", mdl_symbol_type_name(type),
                 name);
  else if (symp->sym_type != type) {
    mdlerror_fmt(
        parse_state, "Invalid type for symbol %s: expected %s, but found %s",
        name, mdl_symbol_type_name(type), mdl_symbol_type_name(symp->sym_type));
    symp = NULL;
  } else if (strcmp(symp->name, "GENERIC_MOLECULE") == 0) {
    mdlerror_fmt(parse_state, "The keyword 'GENERIC_MOLECULE' is obsolete. "
                              "Please use instead 'ALL_VOLUME_MOLECULES', "
                              "'ALL_SURFACE_MOLECULES' or 'ALL_MOLECULES'.");
    symp = NULL;
  }
  free(name);

  return symp;
}

/*************************************************************************
 mdl_existing_symbol_2types:
    Find an existing symbol or print an error message.

 In:  parse_state: parser state
      name: name of symbol to find
      type1: 1st type of symbol to find
      tab1:  1st table
      type2: 2nd type of symbol to find
      tab2:  2nd table
 Out: returns the symbol, or NULL if none found
*************************************************************************/
static struct sym_entry *
mdl_existing_symbol_2types(struct mdlparse_vars *parse_state, char *name,
                           struct sym_table_head *tab1, int type1,
                           struct sym_table_head *tab2, int type2) {
  struct sym_entry *symp;
  symp = retrieve_sym(name, tab1);
  if (symp == NULL) {
    symp = retrieve_sym(name, tab2);
    if (symp == NULL)
      mdlerror_fmt(parse_state, "Undefined %s or %s: %s",
                   mdl_symbol_type_name(type1), mdl_symbol_type_name(type2),
                   name);
  } else {
    if (retrieve_sym(name, tab2) != NULL) {
      mdlerror_fmt(parse_state, "Named object '%s' could refer to %s %s or %s "
                                "%s.  Please rename one of them.",
                   name, mdl_symbol_type_name_article(type1),
                   mdl_symbol_type_name(type1),
                   mdl_symbol_type_name_article(type2),
                   mdl_symbol_type_name(type2));
      symp = NULL;
    }
  }
  free(name);
  return symp;
}

/*************************************************************************
 mdl_find_symbols_by_wildcard:
    Find all symbols of a particular type matching a particular wildcard.

 In:  parse_state: parser state
      wildcard: wildcard to match
      tab: table to search for symbols
      type: type of symbol to match
 Out: linked list of matching symbols
*************************************************************************/
static struct sym_table_list *
mdl_find_symbols_by_wildcard(struct mdlparse_vars *parse_state,
                             char const *wildcard, struct sym_table_head *tab,
                             int type) {
  struct sym_table_list *symbols = NULL, *stl;
  for (int i = 0; i < tab->n_bins; i++) {
    for (struct sym_entry *sym_t = tab->entries[i]; sym_t != NULL;
         sym_t = sym_t->next) {
      if (sym_t->sym_type != type)
        continue;

      if (is_wildcard_match((char *)wildcard, sym_t->name)) {
        stl =
            CHECKED_MEM_GET(parse_state->sym_list_mem, "wildcard symbol list");
        if (stl == NULL) {
          if (symbols)
            mem_put_list(parse_state->sym_list_mem, symbols);
          return NULL;
        }

        stl->node = sym_t;
        stl->next = symbols;
        symbols = stl;
      }
    }
  }

  if (symbols == NULL) {
    switch (type) {
    case OBJ:
      mdlerror_fmt(parse_state, "No objects found matching wildcard \"%s\"",
                   wildcard);
      break;
    case MOL:
      mdlerror_fmt(parse_state, "No molecules found matching wildcard \"%s\"",
                   wildcard);
      break;
    default:
      mdlerror_fmt(parse_state, "No items found matching wildcard \"%s\"",
                   wildcard);
      break;
    }
  }

  return symbols;
}

/*************************************************************************
 compare_sym_names:
    Comparator for symbol list sorting.

 In:  a: first symbol for comparison
      b: second symbol for comparison
 Out: 0 if a > b, 1 if a <= b
*************************************************************************/
static int compare_sym_names(void *a, void *b) {
  return strcmp(((struct sym_entry *)a)->name, ((struct sym_entry *)b)->name) <=
         0;
}

/*************************************************************************
 sort_sym_list_by_name:
    Sort a symbol list in collation order by name.

 In:  unsorted: unsorted list
 Out: a sorted list of symbols
*************************************************************************/
static struct sym_table_list *
sort_sym_list_by_name(struct sym_table_list *unsorted) {
  return (struct sym_table_list *)void_list_sort_by(
      (struct void_list *)unsorted, compare_sym_names);
}

/*************************************************************************
 mdl_existing_object:
    Find an existing object.  Print an error message if the object isn't found.

 In:  parse_state: parser state
      name: fully qualified object name
 Out: the object, or NULL if not found
*************************************************************************/
struct sym_entry *mdl_existing_object(struct mdlparse_vars *parse_state,
                                      char *name) {
  return mdl_existing_symbol(
      parse_state, name, parse_state->vol->obj_sym_table, OBJ);
}

/*************************************************************************
 mdl_existing_objects_wildcard:
    Find a list of existing objects matching a particular wildcard.

 In:  parse_state: parser state
      wildcard: fully qualified object name
 Out: the meshes, or NULL if none were found or allocation failed
*************************************************************************/
struct sym_table_list *
mdl_existing_objects_wildcard(struct mdlparse_vars *parse_state,
                              char *wildcard) {
  char *stripped = mdl_strip_quotes(wildcard);
  return mdl_meshes_by_wildcard(parse_state, stripped);
}

/*************************************************************************
 mdl_existing_region:
    Find an existing region.  Print an error message if it isn't found.

 In:  parse_state: parser state
      obj_symp: object on which to find the region
      name: region name
 Out: the region, or NULL if not found
*************************************************************************/
struct sym_entry *mdl_existing_region(struct mdlparse_vars *parse_state,
                                      struct sym_entry *obj_symp, char *name) {
  char *region_name = CHECKED_SPRINTF("%s,%s", obj_symp->name, name);
  if (region_name == NULL) {
    free(name);
    return NULL;
  }

  free(name);
  return mdl_existing_symbol(
      parse_state, region_name, parse_state->vol->reg_sym_table, REG);
}

/*************************************************************************
 mdl_existing_molecule:
    Find an existing molecule species.  Print an error message if it isn't
    found.

 In:  parse_state: parser state
      name: species name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_entry *mdl_existing_molecule(struct mdlparse_vars *parse_state,
                                        char *name) {
  return mdl_existing_symbol(parse_state, name, parse_state->vol->mol_sym_table,
                             MOL);
}

/**************************************************************************
 mdl_singleton_symbol_list:
    Turn a single symbol into a singleton symbol list.

 In: parse_state: parser state
     sym:  the symbol
 Out: the symbol list, or NULL
**************************************************************************/
struct sym_table_list *
mdl_singleton_symbol_list(struct mdlparse_vars *parse_state,
                          struct sym_entry *sym) {
  struct sym_table_list *stl =
      CHECKED_MEM_GET(parse_state->sym_list_mem, "symbol list item");
  if (stl != NULL) {
    stl->next = NULL;
    stl->node = sym;
  }
  return stl;
}

/*************************************************************************
 mdl_existing_molecule_list:
    Find an existing molecule species, and return it in a singleton list.
    Print an error message if it isn't found.

 In:  parse_state: parser state
      name: species name
 Out: a list containing only the symbol, or NULL if not found or allocation
      failed
*************************************************************************/
struct sym_table_list *
mdl_existing_molecule_list(struct mdlparse_vars *parse_state, char *name) {
  struct sym_entry *symp = mdl_existing_molecule(parse_state, name);
  if (symp == NULL)
    return NULL;

  return mdl_singleton_symbol_list(parse_state, symp);
}

/*************************************************************************
 mdl_existing_molecules_wildcard:
    Find a list of all molecule species matching the specified wildcard.  Print
    an error message if it doesn't match any.

 In:  parse_state: parser state
      wildcard: species wildcard (will be freed by this function)
 Out: a list containing the symbols, or NULL if an error occurred
*************************************************************************/
struct sym_table_list *
mdl_existing_molecules_wildcard(struct mdlparse_vars *parse_state,
                                char *wildcard) {
  struct sym_table_list *stl;
  char *wildcard_string;
  if (!(wildcard_string = mdl_strip_quotes(wildcard)))
    return NULL;

  stl = mdl_find_symbols_by_wildcard(parse_state, wildcard_string,
                                     parse_state->vol->mol_sym_table, MOL);
  if (stl == NULL) {
    free(wildcard_string);
    return NULL;
  }
  free(wildcard_string);

  return sort_sym_list_by_name(stl);
}

/*************************************************************************
 mdl_existing_surface_molecule:
    Find an existing surface molecule species.  Print an error message if it
    isn't found, or isn't a surface molecule.

 In:  parse_state: parser state
      name: species name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_entry *
mdl_existing_surface_molecule(struct mdlparse_vars *parse_state, char *name) {
  struct sym_entry *symp = mdl_existing_molecule(parse_state, name);
  if (symp == NULL)
    return NULL;

  struct species *sp = (struct species *)symp->value;
  if (!(sp->flags & ON_GRID)) {
    mdlerror_fmt(parse_state, "Molecule '%s' is not a surface molecule",
                 symp->name);
    return NULL;
  }

  return symp;
}

/*************************************************************************
 mdl_existing_surface_class:
    Find an existing surface class species.  Print an error message if it isn't
    found, or isn't a surface class.

 In:  parse_state: parser state
      name: species name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_entry *mdl_existing_surface_class(struct mdlparse_vars *parse_state,
                                             char *name) {
  struct sym_entry *symp = mdl_existing_molecule(parse_state, name);
  if (symp == NULL)
    return NULL;

  struct species *sp = (struct species *)symp->value;
  if (!(sp->flags & IS_SURFACE)) {
    mdlerror_fmt(parse_state, "Species '%s' is not a surface class",
                 sp->sym->name);
    return NULL;
  }

  return symp;
}

/**************************************************************************
 mdl_existing_variable:
    Find a named variable if it exists, or print an error if it does not.

 In: parse_state: parser state
     name: the name of the variable
 Out: a symbol table entry, or NULL if we couldn't find the variable
**************************************************************************/
struct sym_entry *mdl_existing_variable(struct mdlparse_vars *parse_state,
                                        char *name) {
  struct sym_entry *st = NULL;

  /* Attempt to fetch existing variable */
  if ((st = retrieve_sym(name, parse_state->vol->var_sym_table)) != NULL) {
    free(name);
    return st;
  }

  mdlerror_fmt(parse_state, "Undefined variable: %s", name);
  free(name);
  return NULL;
}

/*************************************************************************
 mdl_existing_array:
    Find an existing array symbol.  Print an error message if it isn't found.

 In:  parse_state: parser state
      name: symbol name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_entry *mdl_existing_array(struct mdlparse_vars *parse_state,
                                     char *name) {
  return mdl_existing_symbol(parse_state, name, parse_state->vol->var_sym_table,
                             ARRAY);
}

/**************************************************************************
 mdl_existing_double:
    Find a named numeric variable if it exists.  Print an error message if it
    isn't found.

 In: parse_state: the parse variables structure
     name: the name of the variable
 Out: a symbol table entry, or NULL if we couldn't find the variable
**************************************************************************/
struct sym_entry *mdl_existing_double(struct mdlparse_vars *parse_state,
                                      char *name) {
  return mdl_existing_symbol(parse_state, name, parse_state->vol->var_sym_table,
                             DBL);
}

/**************************************************************************
 mdl_existing_string:
    Find a named string variable if it exists.  Print an error message if it
    isn't found.

 In: parse_state: the parse variables structure
     name: the name of the variable
 Out: a symbol table entry, or NULL if we couldn't find the variable
**************************************************************************/
struct sym_entry *mdl_existing_string(struct mdlparse_vars *parse_state,
                                      char *name) {
  return mdl_existing_symbol(parse_state, name, parse_state->vol->var_sym_table,
                             STR);
}

/**************************************************************************
 mdl_existing_num_or_array:
    Find a named numeric or array variable if it exists.  Print an error message
    if it isn't found.

 In: parse_state: the parse variables structure
     name: the name of the variable
 Out: a symbol table entry, or NULL if we couldn't find the variable
**************************************************************************/
struct sym_entry *mdl_existing_num_or_array(struct mdlparse_vars *parse_state,
                                            char *name) {
  struct sym_entry *st = NULL;

  /* Attempt to fetch existing variable */
  if ((st = retrieve_sym(name, parse_state->vol->var_sym_table)) != NULL) {
    if (st->sym_type == STR) {
      mdlerror_fmt(parse_state,
                   "Incorrect type (got string, expected number or array): %s",
                   name);
      return NULL;
    }
    return st;
  }

  mdlerror_fmt(parse_state, "Undefined variable: %s", name);
  return NULL;
}

/*************************************************************************
 mdl_existing_rxn_pathname_or_molecule:
    Find an existing named reaction pathway or molecule.  Print an error
    message if it isn't found.

 In:  parse_state: parser state
      name: symbol name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_entry *
mdl_existing_rxn_pathname_or_molecule(struct mdlparse_vars *parse_state,
                                      char *name) {
  return mdl_existing_symbol_2types(parse_state, name,
                                    parse_state->vol->rxpn_sym_table, RXPN,
                                    parse_state->vol->mol_sym_table, MOL);
}

/*************************************************************************
 mdl_existing_release_pattern_or_rxn_pathname:
    Find an existing reaction pathway or release pattern.  Print an error
    message if it isn't found, or if the name could refer to either a release
    pattern or a reaction pathway.

 In:  parse_state: parser state
      name: symbol name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_entry *
mdl_existing_release_pattern_or_rxn_pathname(struct mdlparse_vars *parse_state,
                                             char *name) {
  return mdl_existing_symbol_2types(parse_state, name,
                                    parse_state->vol->rpat_sym_table, RPAT,
                                    parse_state->vol->rxpn_sym_table, RXPN);
}

/*************************************************************************
 mdl_existing_file_stream:
    Find an existing file stream.  Print an error message if the stream isn't
    found.

 In:  parse_state: parser state
      name: stream name
 Out: the stream, or NULL if not found
*************************************************************************/
struct sym_entry *mdl_existing_file_stream(struct mdlparse_vars *parse_state,
                                           char *name) {
  struct sym_entry *sym = mdl_existing_symbol(
      parse_state, name, parse_state->vol->fstream_sym_table, FSTRM);
  if (sym == NULL)
    return sym;

  struct file_stream *filep = (struct file_stream *)sym->value;
  if (filep->stream == NULL) {
    mdlerror_fmt(parse_state, "File stream '%s' has already been closed",
                 sym->name);
    return NULL;
  }

  return sym;
}

/*************************************************************************
 mdl_meshes_by_wildcard:
    Find all mesh objects (polygons and boxes) that match a given wildcard.

    XXX: Should we include meta-objects in this list?

 In:  parse_state: parser state
      wildcard: the wildcard to match
 Out: a list of all matching symbols
*************************************************************************/
struct sym_table_list *mdl_meshes_by_wildcard(struct mdlparse_vars *parse_state,
                                              char *wildcard) {
  /* Scan for objects matching the wildcard */
  struct sym_table_list *matches = mdl_find_symbols_by_wildcard(
      parse_state, wildcard, parse_state->vol->obj_sym_table, OBJ);
  if (matches == NULL) {
    free(wildcard);
    return NULL;
  }

  /* Scan through the list, discarding inappropriate objects */
  struct sym_table_list *cur_match, *next_match, **prev = &matches;
  for (cur_match = matches; cur_match != NULL; cur_match = next_match) {
    next_match = cur_match->next;

    struct object *objp = (struct object *)cur_match->node->value;
    if (objp->object_type != POLY_OBJ && objp->object_type != BOX_OBJ) {
      *prev = cur_match->next;
      mem_put(parse_state->sym_list_mem, cur_match);
    } else
      prev = &cur_match->next;
  }

  if (matches == NULL) {
    mdlerror_fmt(parse_state,
                 "No matches for the wildcard '%s' were mesh objects",
                 wildcard);
    free(wildcard);
    return NULL;
  }

  free(wildcard);
  return sort_sym_list_by_name(matches);
}

/*************************************************************************
 mdl_transform_rotate:
    Apply a rotation to the given transformation matrix.

 In:  parse_state: parser state
      mat: transformation matrix
      axis: axis of rotation
      angle: angle of rotation (degrees!)
 Out: 0 on success, 1 on failure; rotation is right-multiplied into the
      transformation matrix
*************************************************************************/
int mdl_transform_rotate(struct mdlparse_vars *parse_state, double (*mat)[4],
                         struct vector3 *axis, double angle) {
  if (transform_rotate(mat, axis, angle)) {
    mdlerror(parse_state, "Rotation vector has zero length.");
    return 1;
  }
  free(axis);
  return 0;
}

/*************************************************************************
 mdl_make_new_region:
    Create a new region, adding it to the global symbol table.  The region must
    not be defined yet.  The region is not added to the object's list of
    regions.

    full region names of REG type symbols stored in main symbol table have the
    form:
         metaobj.metaobj.poly,region_last_name

 In:  parse_state: parse vars for error reporting
      obj_name: fully qualified object name
      region_last_name: name of the region to define
 Out: The newly created region
*************************************************************************/
static struct region *mdl_make_new_region(struct mdlparse_vars *parse_state,
                                          char *obj_name,
                                          char *region_last_name) {
  char *region_name;
  region_name = CHECKED_SPRINTF("%s,%s", obj_name, region_last_name);
  if (region_name == NULL)
    return NULL;

  struct sym_entry *gp;
  if ((gp = retrieve_sym(region_name, parse_state->vol->reg_sym_table)) != NULL) {
    if (gp->count == 0) {
      gp->count = 1;
      free(region_name);
      return (struct region *)gp->value;
    }
    else {
      mdlerror_fmt(parse_state, "Region already defined: %s", region_name);
    }
  }

  if ((gp = store_sym(region_name, REG, parse_state->vol->reg_sym_table,
                      NULL)) == NULL) {
    free(region_name);
    mcell_allocfailed("Failed to store a region in the region symbol table.");
  }

  free(region_name);
  return (struct region *)gp->value;
}

/*************************************************************************
 mdl_copy_object_regions:
    Duplicate src_obj's regions on dst_obj.

 In:  parse_state: parser state
      dst_obj: destination object
      src_obj: object from which to copy
 Out: 0 on success, 1 on failure
*************************************************************************/
static int mdl_copy_object_regions(struct mdlparse_vars *parse_state,
                                   struct object *dst_obj,
                                   struct object *src_obj) {
  struct region_list *src_rlp;
  struct region *dst_reg, *src_reg;
  struct sm_dat *dst_sm, *src_sm;

  /* Copy each region */
  for (src_rlp = src_obj->regions; src_rlp != NULL; src_rlp = src_rlp->next) {
    src_reg = src_rlp->reg;
    if ((dst_reg = mdl_create_region(parse_state, dst_obj,
                                     src_reg->region_last_name)) == NULL)
      return 1;

    /* Copy over simple region attributes */
    dst_reg->surf_class = src_reg->surf_class;
    dst_reg->flags = src_reg->flags;
    dst_reg->area = src_reg->area;
    dst_reg->bbox = src_reg->bbox;
    dst_reg->manifold_flag = src_reg->manifold_flag;

    /* Copy membership data */
    if (src_reg->membership != NULL) {
      dst_reg->membership = duplicate_bit_array(src_reg->membership);
      if (dst_reg->membership == NULL) {
        mcell_allocfailed("Failed allocation for membership array in %s",
                          dst_obj->sym->name);
        /*return 1;*/
      }
    } else
      mdl_warning(parse_state, "No membership data for %s\n",
                  src_reg->sym->name);

    /* Copy surface molecule data list */
    struct sm_dat *smdp_tail = NULL;
    for (src_sm = src_reg->sm_dat_head; src_sm != NULL; src_sm = src_sm->next) {
      dst_sm = CHECKED_MALLOC_STRUCT(struct sm_dat, "surface molecule info");
      if (dst_sm == NULL)
        return 1;

      if (smdp_tail != NULL)
        smdp_tail->next = dst_sm;
      else
        dst_reg->sm_dat_head = dst_sm;
      smdp_tail = dst_sm;

      dst_sm->sm = src_sm->sm;
      dst_sm->quantity_type = src_sm->quantity_type;
      dst_sm->quantity = src_sm->quantity;
      dst_sm->orientation = src_sm->orientation;
      dst_sm->next = NULL;
    }
  }
  return 0;
}

/*************************************************************************
 find_corresponding_region:
    When an object is instantiated or included in another object, we may need
    to fix up region references for region releases.  This function performs
    this fixup.

 In: old_r: a region
     old_ob: the object that referred to that region
     new_ob: a new object that should refer to its corresponding region
     instance: the root object that begins the instance tree
     symhash: the main symbol hash table for the world
 Out: a pointer to the region that has the same relationship to the new object
      as the given region has to the old object, or NULL if there is no
      corresponding region.
 Note: correspondence is computed by name mangling; if an object A.B.C refers
       to region A.B.D,R and we have a new object E.F.G.H, then the
       corresponding region is considered to be E.F.G.D,R (i.e. go back one to
       find common ancestor, then go forward into the thing labeled "D,R").
*************************************************************************/
static struct region *
find_corresponding_region(struct region *old_r, struct object *old_ob,
                          struct object *new_ob, struct object *instance,
                          struct sym_table_head *symhash) {
  struct object *ancestor;
  struct object *ob;
  struct sym_entry *gp;

  ancestor = common_ancestor(old_ob, old_r->parent);

  if (ancestor == NULL) {
    for (ob = old_r->parent; ob != NULL; ob = ob->parent) {
      if (ob == instance)
        break;
    }

    if (ob == NULL)
      return NULL;
    else
      return old_r; /* Point to same already-instanced object */
  } else {
    /* Length of old prefix and name */
    const size_t old_prefix_idx = strlen(ancestor->sym->name);
    const size_t old_name_len = strlen(old_ob->sym->name);

    /* Suffix is everything in the old name after the old prefix. */
    const size_t suffix_len = old_name_len - old_prefix_idx;

    /* New prefix length is new name, less suffix length. */
    const size_t new_name_len = strlen(new_ob->sym->name);
    /* If we get here, we must have a common ancestor, so this had better work.
     */
    assert(new_name_len > suffix_len);
    const size_t new_prefix_idx = new_name_len - suffix_len;

    /* Length of the "region" spec from the old region name is the length of
     * the region symbol, less the length of the old object symbol. */
    const size_t region_name_len = strlen(old_r->sym->name);
    assert(region_name_len > old_prefix_idx);
    const size_t just_region_name_len = region_name_len - old_prefix_idx;

    /* Buffer size needed for new name is new object name length + region name
     * length + 1 (for NUL-termination). */
    const size_t max_name_len = new_name_len + just_region_name_len + 1;

    /* Build new name. */
    char new_name[max_name_len];
    strncpy(new_name, new_ob->sym->name, new_prefix_idx);
    strncpy(new_name + new_prefix_idx,         /* new prefix */
            old_r->sym->name + old_prefix_idx, /* old suffix + region name */
            max_name_len - new_prefix_idx);

    /* Finally, retrieve symbol from newly-constructed name. */
    gp = retrieve_sym(new_name, symhash);
    if (gp == NULL)
      return NULL;
    else
      return (struct region *)gp->value;
  }
}

/*************************************************************************
 duplicate_rel_region_expr:
    Create a deep copy of a region release expression.

 In: parse_state: parser state
     expr: a region expression tree
     old_self: the object containing that tree
     new_self: a new object for which we want to build a corresponding tree
     instance: the root object that begins the instance tree
 Out: the newly constructed expression tree for the new object, or
      NULL if no such tree can be built
*************************************************************************/
static struct release_evaluator *duplicate_rel_region_expr(
    struct mdlparse_vars *parse_state, struct release_evaluator *expr,
    struct object *old_self, struct object *new_self, struct object *instance) {
  struct release_evaluator *nexp = CHECKED_MALLOC_STRUCT(
      struct release_evaluator, "region release expression");
  if (nexp == NULL)
    return NULL;

  nexp->op = expr->op;

  struct region *r;
  if (expr->left != NULL) {
    if (expr->op & REXP_LEFT_REGION) {
      r = find_corresponding_region(expr->left, old_self, new_self, instance,
                                    parse_state->vol->reg_sym_table);

      if (r == NULL) {
        mdlerror_fmt(
            parse_state,
            "Can't find new region corresponding to %s for %s (copy of %s)",
            ((struct region *)expr->left)->sym->name, new_self->sym->name,
            old_self->sym->name);
        free(nexp);
        return NULL;
      }

      nexp->left = r;
    } else
      nexp->left = duplicate_rel_region_expr(parse_state, expr->left, old_self,
                                             new_self, instance);
  } else
    nexp->left = NULL;

  if (expr->right != NULL) {
    if (expr->op & REXP_RIGHT_REGION) {
      r = find_corresponding_region(expr->right, old_self, new_self, instance,
                                    parse_state->vol->reg_sym_table);

      if (r == NULL) {
        mdlerror_fmt(
            parse_state,
            "Can't find new region corresponding to %s for %s (copy of %s)",
            ((struct region *)expr->right)->sym->name, new_self->sym->name,
            old_self->sym->name);
        free(nexp);
        return NULL;
      }

      nexp->right = r;
    } else
      nexp->right = duplicate_rel_region_expr(parse_state, expr->right,
                                              old_self, new_self, instance);
  } else
    nexp->right = NULL;

  return nexp;
}

/*************************************************************************
 duplicate_release_site:
    Create a deep copy of a release-site.

 In: parse_state: parser state
     old: an existing release site object
     new_self: the object that is to contain a duplicate release site object
     instance: the root object that begins the instance tree
 Out: a duplicated release site object, or NULL if the release site cannot be
      duplicated.

 N.B.: In order to give a proper report, we always need to duplicate the
       release site, since the copy of the release site needs to have a
       different name.  Otherwise, the user will see confusing reports of
       several releases from the same release site.
*************************************************************************/
static struct release_site_obj *
duplicate_release_site(struct mdlparse_vars *parse_state,
                       struct release_site_obj *old, struct object *new_self,
                       struct object *instance) {
  struct release_site_obj *rel_site_obj =
      CHECKED_MALLOC_STRUCT(struct release_site_obj, "release site");
  if (rel_site_obj == NULL) {
    return NULL;
  }

  if (old->location != NULL) {
    if ((rel_site_obj->location = CHECKED_MALLOC_STRUCT(
             struct vector3, "release site location")) == NULL) {
      free(rel_site_obj);
      return NULL;
    }
    *(rel_site_obj->location) = *(old->location);
  } else {
    rel_site_obj->location = NULL;
  }
  rel_site_obj->mol_type = old->mol_type;
  rel_site_obj->release_number_method = old->release_number_method;
  rel_site_obj->release_shape = old->release_shape;
  rel_site_obj->orientation = old->orientation;
  rel_site_obj->release_number = old->release_number;
  rel_site_obj->mean_diameter = old->mean_diameter;
  rel_site_obj->concentration = old->concentration;
  rel_site_obj->standard_deviation = old->standard_deviation;
  rel_site_obj->diameter = old->diameter;
  rel_site_obj->release_prob = old->release_prob;
  rel_site_obj->pattern = old->pattern;
  rel_site_obj->mol_list = old->mol_list;
  rel_site_obj->periodic_box = old->periodic_box;
  rel_site_obj->name = NULL;

  if (old->region_data != NULL) {
    struct release_region_data *rel_reg_data = CHECKED_MALLOC_STRUCT(
        struct release_region_data, "release region data");
    if (rel_reg_data == NULL) {
      free(rel_site_obj);
      return NULL;
    }

    memcpy(&(rel_reg_data->llf), &(old->region_data->llf),
           sizeof(struct vector3));
    memcpy(&(rel_reg_data->urb), &(old->region_data->urb),
           sizeof(struct vector3));
    rel_reg_data->n_walls_included = -1;
    rel_reg_data->cum_area_list = NULL;
    rel_reg_data->wall_index = NULL;
    rel_reg_data->obj_index = NULL;
    rel_reg_data->n_objects = -1;
    rel_reg_data->owners = NULL;
    rel_reg_data->in_release = NULL;
    rel_reg_data->self = new_self;

    rel_reg_data->expression =
        duplicate_rel_region_expr(parse_state, old->region_data->expression,
                                  old->region_data->self, new_self, instance);
    if (rel_reg_data->expression == NULL) {
      free(rel_site_obj);
      free(rel_reg_data);
      return NULL;
    }

    rel_site_obj->region_data = rel_reg_data;
  } else {
    rel_site_obj->region_data = NULL;
  }

  return rel_site_obj;
}

/*************************************************************************
 mdl_deep_copy_object:
    Deep copy an object.  The destination object should already be added to the
    symbol table, but should be otherwise unpopulated, as no effort is made to
    free any existing data contained in the object.

 In:  parse_state: parse vars for error reporting
      dst_obj: object into which to copy
      src_obj: object from which to copy
 Out: The newly created region
*************************************************************************/
int mdl_deep_copy_object(struct mdlparse_vars *parse_state,
                         struct object *dst_obj, struct object *src_obj) {

  /* Copy over simple object attributes */
  dst_obj->object_type = src_obj->object_type;
  dst_obj->n_walls = src_obj->n_walls;
  dst_obj->n_walls_actual = src_obj->n_walls_actual;
  dst_obj->walls = src_obj->walls;
  dst_obj->wall_p = src_obj->wall_p;
  dst_obj->n_verts = src_obj->n_verts;
  dst_obj->vertices = src_obj->vertices;

  /* Copy over regions */
  if (mdl_copy_object_regions(parse_state, dst_obj, src_obj))
    return 1;

  /* Inherit object coordinate transformations */
  mult_matrix(dst_obj->t_matrix, src_obj->t_matrix, dst_obj->t_matrix, 4, 4, 4);

  switch (dst_obj->object_type) {
  case META_OBJ:
    /* Copy children */
    for (struct object *src_child = src_obj->first_child; src_child != NULL;
         src_child = src_child->next) {
      struct object *dst_child;
      char *child_obj_name =
          CHECKED_SPRINTF("%s.%s", dst_obj->sym->name, src_child->last_name);
      if (child_obj_name == NULL)
        return 1;

      /* Create child object */
      if ((dst_child = mdl_make_new_object(parse_state, child_obj_name)) ==
          NULL) {
        free(child_obj_name);
        return 1;
      }
      free(child_obj_name);

      /* Copy in last name */
      dst_child->last_name = mdl_strdup(src_child->last_name);
      if (dst_child->last_name == NULL)
        return 1;

      /* Recursively copy object and its children */
      if (mdl_deep_copy_object(parse_state, dst_child, src_child))
        return 1;
      dst_child->parent = dst_obj;
      dst_child->next = NULL;
      add_child_objects(dst_obj, dst_child, dst_child);
    }
    break;

  case REL_SITE_OBJ:
    dst_obj->contents =
        duplicate_release_site(parse_state, src_obj->contents, dst_obj,
                               parse_state->vol->root_instance);
    if (dst_obj->contents == NULL)
      return 1;
    struct release_site_obj *rso = (struct release_site_obj *)dst_obj->contents;
    rso->name = mdl_strdup(dst_obj->sym->name);
    if (rso->name == NULL)
      return 1;
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (parse_state->vol->disable_polygon_objects) {
      mdlerror(
          parse_state,
          "When using dynamic geometries, polygon objects should only be "
          "defined/instantiated through the dynamic geometry file.");
    }
    dst_obj->contents = src_obj->contents;
    struct polygon_object *poly_obj = \
        (struct polygon_object *)src_obj->contents;
    // Effectively, this tracks the instances of this object, which we need for
    // cleaning up after dynamic geometry events.
    poly_obj->references++;
    dst_obj->periodic_x = src_obj->periodic_x;
    dst_obj->periodic_y = src_obj->periodic_y;
    dst_obj->periodic_z = src_obj->periodic_z;
    break;

  case VOXEL_OBJ:
  default:
    mdlerror_fmt(parse_state, "Error: bad object type %d",
                 dst_obj->object_type);
    return 1;
  }

  return 0;
}

/*************************************************************************
 * Cuboid processing
 ************************************************************************/

/*************************************************************************
 init_cuboid:

 In: parse_state: parser state
     p1: llf corner of a cube
     p2: urb corner of a cube
 Out: returns a subdivided_box struct, with no subdivisions and corners as
      specified.  NULL is returned if there is no memory or the urb corner is
      not up from, to the right of, and behind the llf corner
*************************************************************************/
static struct subdivided_box *init_cuboid(struct mdlparse_vars *parse_state,
                                          struct vector3 *p1,
                                          struct vector3 *p2) {

  if (p2->x - p1->x < EPS_C || p2->y - p1->y < EPS_C || p2->z - p1->z < EPS_C) {
    mdlerror(parse_state, "Box vertices out of order or box is degenerate.");
    return NULL;
  }

  struct subdivided_box *b = CHECKED_MALLOC_STRUCT(
      struct subdivided_box, "subdivided box");
  if (b == NULL)
    return NULL;

  b->nx = b->ny = b->nz = 2;
  if ((b->x = CHECKED_MALLOC_ARRAY(
      double, b->nx, "subdivided box X partitions")) == NULL) {
    free(b);
    return NULL;
  }
  if ((b->y = CHECKED_MALLOC_ARRAY(
      double, b->ny, "subdivided box Y partitions")) == NULL) {
    free(b);
    return NULL;
  }
  if ((b->z = CHECKED_MALLOC_ARRAY(
      double, b->nz, "subdivided box Z partitions")) == NULL) {
    free(b);
    return NULL;
  }

  b->x[0] = p1->x;
  b->x[1] = p2->x;
  b->y[0] = p1->y;
  b->y[1] = p2->y;
  b->z[0] = p1->z;
  b->z[1] = p2->z;

  return b;
}

/*************************************************************************
 refine_cuboid:
 In: parse_state: parser state
     p1: 3D vector that is one corner of the patch
     p2: 3D vector that is the other corner of the patch
     b: a subdivided box upon which the patch will be placed
     egd: the surface molecule grid density, which limits how fine divisions
          can be
 Out: returns 1 on failure, 0 on success.  The box has additional subdivisions
      added that fall at the edges of the patch so that the patch can be
      specified in terms of subdivisions (i.e. can be constructed by triangles
      that tile each piece of the subdivided surface).
*************************************************************************/
static int refine_cuboid(struct mdlparse_vars *parse_state, struct vector3 *p1,
                         struct vector3 *p2, struct subdivided_box *b,
                         double egd) {
  int j, k;
  double *new_list;
  int new_n;

  int i = check_patch(b, p1, p2, egd);
  if (i == 0) {
    mdlerror(parse_state,
             "Could not refine box to include patch: Invalid patch specified");
    return 1;
  }

  if (i & BRANCH_X) {
    new_n = b->nx + 2;
    for (j = 0; j < b->nx; j++) {
      if (!distinguishable(p1->x, b->x[j], EPS_C))
        new_n--;
      if (!distinguishable(p2->x, b->x[j], EPS_C))
        new_n--;
    }
    if (new_n > b->nx) {
      new_list = CHECKED_MALLOC_ARRAY(double, new_n,
                                      "refined subdivided box X partitions");
      if (new_list == NULL)
        return 1;

      for (j = k = 0; b->x[j] < p1->x; j++)
        new_list[k++] = b->x[j];
      if (distinguishable(b->x[j], p1->x, EPS_C))
        new_list[k++] = p1->x;
      for (; b->x[j] < p2->x; j++)
        new_list[k++] = b->x[j];
      if ((distinguishable(p1->x, p2->x, EPS_C)) &&
          (distinguishable(b->x[j], p2->x, EPS_C))) {
        new_list[k++] = p2->x;
      }
      for (; j < b->nx; j++)
        new_list[k++] = b->x[j];

      free(b->x);
      b->x = new_list;
      b->nx = new_n;
    }
  }
  if (i & BRANCH_Y) /* Same as above with x->y */
  {
    new_n = b->ny + 2;
    for (j = 0; j < b->ny; j++) {
      if (!distinguishable(p1->y, b->y[j], EPS_C))
        new_n--;
      if (!distinguishable(p2->y, b->y[j], EPS_C))
        new_n--;
    }
    if (new_n > b->ny) {
      new_list = CHECKED_MALLOC_ARRAY(double, new_n,
                                      "refined subdivided box Y partitions");
      if (new_list == NULL)
        return 1;

      for (j = k = 0; b->y[j] < p1->y; j++)
        new_list[k++] = b->y[j];
      if (distinguishable(b->y[j], p1->y, EPS_C))
        new_list[k++] = p1->y;
      for (; b->y[j] < p2->y; j++)
        new_list[k++] = b->y[j];
      if ((distinguishable(p1->y, p2->y, EPS_C)) &&
          (distinguishable(b->y[j], p2->y, EPS_C))) {
        new_list[k++] = p2->y;
      }
      for (; j < b->ny; j++)
        new_list[k++] = b->y[j];

      free(b->y);
      b->y = new_list;
      b->ny = new_n;
    }
  }
  if (i & BRANCH_Z) /* Same again, x->z */
  {
    new_n = b->nz + 2;
    for (j = 0; j < b->nz; j++) {
      if (!distinguishable(p1->z, b->z[j], EPS_C))
        new_n--;
      if (!distinguishable(p2->z, b->z[j], EPS_C))
        new_n--;
    }
    if (new_n > b->nz) {
      new_list = CHECKED_MALLOC_ARRAY(double, new_n,
                                      "refined subdivided box Z partitions");
      if (new_list == NULL)
        return 1;

      for (j = k = 0; b->z[j] < p1->z; j++)
        new_list[k++] = b->z[j];
      if (distinguishable(b->z[j], p1->z, EPS_C))
        new_list[k++] = p1->z;
      for (; b->z[j] < p2->z; j++)
        new_list[k++] = b->z[j];
      if ((distinguishable(p1->z, p2->z, EPS_C) &&
          (distinguishable(b->z[j], p2->z, EPS_C)))) {
        new_list[k++] = p2->z;
      }
      for (; j < b->nz; j++)
        new_list[k++] = b->z[j];

      free(b->z);
      b->z = new_list;
      b->nz = new_n;
    }
  }

  return 0;
}

/*************************************************************************
 divide_cuboid:

 In: parse_state: parser state
     b: a subdivided box to further subdivide
     axis: which axis to divide
     idx: which of the existing divisions should be subdivided
     ndiv: the number of subdivisions to make
 Out: returns 1 on failure, 0 on success.  The requested subdivision(s) are
      added.
*************************************************************************/
static int divide_cuboid(struct mdlparse_vars *parse_state,
                         struct subdivided_box *b, int axis, int idx,
                         int ndiv) {
  double *old_list;
  double *new_list;
  int old_n;
  int new_n;
  int i, j, k;

  if (ndiv < 2)
    ndiv = 2;
  switch (axis) {
  case BRANCH_X:
    old_list = b->x;
    old_n = b->nx;
    break;
  case BRANCH_Y:
    old_list = b->y;
    old_n = b->ny;
    break;
  case BRANCH_Z:
    old_list = b->z;
    old_n = b->nz;
    break;
  default:
    mdlerror_fmt(parse_state, "File '%s', Line %ld: Unknown flag is used.",
                 __FILE__, (long)__LINE__);
    return 1;
  }

  new_n = old_n + ndiv - 1;
  new_list =
      CHECKED_MALLOC_ARRAY(double, new_n, "refined subdivided box partitions");
  if (new_list == NULL)
    return 1;

  for (i = j = 0; i <= idx; i++, j++)
    new_list[j] = old_list[i];
  for (k = 1; k < ndiv; k++)
    new_list[j++] =
        (((double)k / (double)ndiv)) * (old_list[i] - old_list[i - 1]) +
        old_list[i - 1];
  for (; i < old_n; i++, j++)
    new_list[j] = old_list[i];

  switch (axis) {
  case BRANCH_X:
    b->x = new_list;
    b->nx = new_n;
    break;
  case BRANCH_Y:
    b->y = new_list;
    b->ny = new_n;
    break;
  case BRANCH_Z:
    b->z = new_list;
    b->nz = new_n;
    break;
  default:
    return 1;
    /*break;*/
  }

  free(old_list);
  return 0;
}

/*************************************************************************
 reaspect_cuboid:
 In: parse_state: parser state
     b: a subdivided box whose surface is a bunch of rectangles
     ratio: the maximum allowed aspect ratio (long side / short side) of a
rectangle
 Out: returns 1 on failure, 0 on success.  The subdivided box is further
      divided to ensure that all surface subdivisions meet the aspect ratio
      criterion.
 Note: This is not, in general, possible if max_ratio is less than sqrt(2), and
       it is diffucult if max_ratio is less than 2.
*************************************************************************/
static int reaspect_cuboid(struct mdlparse_vars *parse_state,
                           struct subdivided_box *b, double max_ratio) {
  double min_x, min_y, min_z, max_x, max_y, max_z;
  int jx, jy, jz;
  int i, j;
  int changed;

  do {
    changed = 0;

    max_x = min_x = b->x[1] - b->x[0];
    jx = 0;
    for (i = 1; i < b->nx - 1; i++) {
      if (min_x > b->x[i + 1] - b->x[i]) {
        min_x = b->x[i + 1] - b->x[i];
      } else if (max_x < b->x[i + 1] - b->x[i]) {
        max_x = b->x[i + 1] - b->x[i];
        jx = i;
      }
    }

    max_y = min_y = b->y[1] - b->y[0];
    jy = 0;
    for (i = 1; i < b->ny - 1; i++) {
      if (min_y > b->y[i + 1] - b->y[i]) {
        min_y = b->y[i + 1] - b->y[i];
      } else if (max_y < b->y[i + 1] - b->y[i]) {
        max_y = b->y[i + 1] - b->y[i];
        jy = i;
      }
    }

    max_z = min_z = b->z[1] - b->z[0];
    jz = 0;
    for (i = 1; i < b->nz - 1; i++) {
      if (min_z > b->z[i + 1] - b->z[i]) {
        min_z = b->z[i + 1] - b->z[i];
      } else if (max_z < b->z[i + 1] - b->z[i]) {
        max_z = b->z[i + 1] - b->z[i];
        jz = i;
      }
    }

    if (max_y / min_x > max_ratio) {
      j = divide_cuboid(parse_state, b, BRANCH_Y, jy,
                        (int)ceil(max_y / (max_ratio * min_x)));
      if (j)
        return 1;
      changed |= BRANCH_Y;
    } else if (max_x / min_y > max_ratio) {
      j = divide_cuboid(parse_state, b, BRANCH_X, jx,
                        (int)ceil(max_x / (max_ratio * min_y)));
      if (j)
        return 1;
      changed |= BRANCH_X;
    }

    if ((changed & BRANCH_X) == 0 && max_z / min_x > max_ratio) {
      j = divide_cuboid(parse_state, b, BRANCH_Z, jz,
                        (int)ceil(max_z / (max_ratio * min_x)));
      if (j)
        return 1;
      changed |= BRANCH_Z;
    } else if ((changed & BRANCH_X) == 0 && max_x / min_z > max_ratio) {
      j = divide_cuboid(parse_state, b, BRANCH_X, jx,
                        (int)ceil(max_x / (max_ratio * min_z)));
      if (j)
        return 1;
      changed |= BRANCH_X;
    }

    if ((changed & (BRANCH_Y | BRANCH_Z)) == 0 && max_z / min_y > max_ratio) {
      j = divide_cuboid(parse_state, b, BRANCH_Z, jz,
                        (int)ceil(max_z / (max_ratio * min_y)));
      if (j)
        return 1;
      changed |= BRANCH_Z;
    } else if ((changed & (BRANCH_Y | BRANCH_Z)) == 0 &&
               max_y / min_z > max_ratio) {
      j = divide_cuboid(parse_state, b, BRANCH_Y, jy,
                        (int)ceil(max_y / (max_ratio * min_z)));
      if (j)
        return 1;
      changed |= BRANCH_Y;
    }
  } while (changed);

  return 0;
}

/*************************************************************************
 count_cuboid_vertices:
    Trivial utility function that counts # vertices in a box

 In: sb: a subdivided box
 Out: the number of vertices in the box
*************************************************************************/
static int count_cuboid_vertices(struct subdivided_box *sb) {
  return 2 * sb->ny * sb->nz + 2 * (sb->nx - 2) * sb->nz +
         2 * (sb->nx - 2) * (sb->ny - 2);
}

/*************************************************************************
 mdl_normalize_elements:
    Prepare a region description for use in the simulation by creating a
    membership bitmask on the region object.

 In: parse_state: parser state
     reg: a region
     existing: a flag indicating whether the region is being modified or
               created
 Out: returns 1 on failure, 0 on success.  Lists of element specifiers
      are converted into bitmasks that specify whether a wall is in or
      not in that region.  This also handles combinations of regions.
*************************************************************************/
int mdl_normalize_elements(struct mdlparse_vars *parse_state,
                           struct region *reg, int existing) {
  struct bit_array *temp = NULL;
  char op;
  unsigned int num_elems;

  if (reg->element_list_head == NULL) {
    return 0;
  }

  struct polygon_object *poly_obj = NULL;
  if (reg->parent->object_type == BOX_OBJ) {
    poly_obj = (struct polygon_object *)reg->parent->contents;
    num_elems = count_cuboid_elements(poly_obj->sb);
  } else {
    num_elems = reg->parent->n_walls;
  }

  struct bit_array *elem_array;
  if (reg->membership == NULL) {
    elem_array = new_bit_array(num_elems);
    if (elem_array == NULL) {
      mcell_allocfailed("Failed to allocate a region membership bitmask.");
      /*return 1;*/
    }
    reg->membership = elem_array;
  } else {
    elem_array = reg->membership;
  }

  if (reg->element_list_head->special == NULL) {
    set_all_bits(elem_array, 0);
  }
  // Special flag for exclusion
  else if ((void *)reg->element_list_head->special ==
           (void *)reg->element_list_head) {
    set_all_bits(elem_array, 1);
  } else {
    if (reg->element_list_head->special->exclude) {
      set_all_bits(elem_array, 1);
    } else {
      set_all_bits(elem_array, 0);
    }
  }

  int i = 0;
  struct element_list *elem_list;
  for (elem_list = reg->element_list_head; elem_list != NULL;
       elem_list = elem_list->next) {
    if (reg->parent->object_type == BOX_OBJ) {
      assert(poly_obj != NULL);
      i = elem_list->begin;
      switch (i) {
      case X_NEG:
        elem_list->begin = 0;
        elem_list->end =
            2 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) - 1;
        break;
      case X_POS:
        elem_list->begin = 2 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1);
        elem_list->end =
            4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) - 1;
        break;
      case Y_NEG:
        elem_list->begin = 4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1);
        elem_list->end = elem_list->begin +
                         2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1) -
                         1;
        break;
      case Y_POS:
        elem_list->begin = 4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) +
                           2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1);
        elem_list->end = elem_list->begin +
                         2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1) -
                         1;
        break;
      case Z_NEG:
        elem_list->begin = 4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) +
                           4 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1);
        elem_list->end = elem_list->begin +
                         2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->ny - 1) -
                         1;
        break;
      case Z_POS:
        elem_list->end = num_elems - 1;
        elem_list->begin = elem_list->end + 1 -
                           2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->ny - 1);
        break;
      case ALL_SIDES:
        elem_list->begin = 0;
        elem_list->end = num_elems - 1;
        break;
      default:
        UNHANDLED_CASE(i);
        /*return 1;*/
      }
    } else if (elem_list->begin >= (u_int)num_elems ||
               elem_list->end >= (u_int)num_elems) {
      mdlerror_fmt(parse_state, "Region element specifier for region %s[%s] refers to sides "
                                "%u...%u, but polygon has only %u sides.",
                   reg->parent->sym->name, reg->region_last_name, elem_list->begin, elem_list->end, num_elems);
      return 1;
    }

    if (elem_list->special == NULL) {
      set_bit_range(elem_array, elem_list->begin, elem_list->end, 1);
    } else if ((void *)elem_list->special == (void *)elem_list) {
      set_bit_range(elem_array, elem_list->begin, elem_list->end, 0);
    } else {
      if (elem_list->special->referent != NULL) {
        if (elem_list->special->referent->membership == NULL) {
          if (elem_list->special->referent->element_list_head != NULL) {
            i = mdl_normalize_elements(parse_state,
                                       elem_list->special->referent, existing);
            if (i) {
              free_bit_array(temp);
              return i;
            }
          }
        }
        if (elem_list->special->referent->membership != NULL) {
          // What does it mean for the membership array to have length zero?
          if (elem_list->special->referent->membership->nbits == 0) {
            if (elem_list->special->exclude) {
              set_all_bits(elem_array, 0);
            } else {
              set_all_bits(elem_array, 1);
            }
          } else {
            if (elem_list->special->exclude) {
              op = '-';
            } else {
              op = '+';
            }
            bit_operation(elem_array, elem_list->special->referent->membership,
                          op);
          }
        }
      } else {
        int ii;
        if (temp == NULL) {
          temp = new_bit_array(num_elems);
          if (temp == NULL) {
            mcell_allocfailed(
                "Failed to allocate a region membership bitmask.");
            /*return 1;*/
          }
        }
        if (poly_obj == NULL) {
          mcell_internal_error("Attempt to create a PATCH on a POLYGON_LIST.");
          /*return 1;*/
        }
        if (existing) {
          mcell_internal_error(
              "Attempt to create a PATCH on an already triangulated BOX.");
          /*return 1;*/
        }
        if (elem_list->special->exclude) {
          op = '-';
        } else {
          op = '+';
        }

        ii = cuboid_patch_to_bits(poly_obj->sb, &(elem_list->special->corner1),
                                  &(elem_list->special->corner2), temp);
        if (ii) {
          // Something wrong with patch.
          free_bit_array(temp);
          return 1;
        }
        bit_operation(elem_array, temp, op);
      }
    }
  }

  if (temp != NULL) {
    free_bit_array(temp);
  }

  if (existing) {
    bit_operation(
        elem_array,
        ((struct polygon_object *)reg->parent->contents)->side_removed, '-');
  }

  while (reg->element_list_head) {
    struct element_list *next = reg->element_list_head->next;
    if (reg->element_list_head->special) {
      if (reg->element_list_head->special !=
          (struct element_special *)reg->element_list_head) {
        free(reg->element_list_head->special);
      }
    }
    free(reg->element_list_head);
    reg->element_list_head = next;
  }

#ifdef DEBUG
  printf("Normalized membership of %s: ", reg->sym->name);
  for (i = 0; i < reg->membership->nbits; i++) {
    if (get_bit(reg->membership, i)) {
      printf("X");
    } else {
      printf("_");
    }
  }
  printf("\n");
#endif

  return 0;
}

/*************************************************************************
 vertex_at_index:

 In: sb: a subdivided box from which we want to retrieve one surface patch
     ix: the x index of the patch
     iy: the y index of the patch
     iz: the z index of the patch
 Out: returns the index into the array of walls for the first wall in the
      patch; add 1 to get the second triangle.  If an invalid coordinate is
      given, -1 is returned.
 Note: since the patch must be on the surface, at least one of ix, iy, iz must
       be either 0 or at its maximum.
*************************************************************************/
static int vertex_at_index(struct subdivided_box *sb, int ix, int iy, int iz) {

  if (ix == 0 || ix == sb->nx - 1) {
    int i = sb->ny * iz + iy;
    if (ix == 0)
      return i;
    else
      return i + sb->ny * sb->nz;
  } else if (iy == 0 || iy == sb->ny - 1) {
    int i = 2 * sb->ny * sb->nz + (sb->nx - 2) * iz + (ix - 1);
    if (iy == 0)
      return i;
    else
      return i + (sb->nx - 2) * sb->nz;
  } else if (iz == 0 || iz == sb->nz - 1) {
    int i = 2 * sb->ny * sb->nz + 2 * (sb->nx - 2) * sb->nz +
            (sb->nx - 2) * (iy - 1) + (ix - 1);
    if (iz == 0)
      return i;
    else
      return i + (sb->nx - 2) * (sb->ny - 2);
  } else {
    mcell_internal_error(
        "Asking for point %d %d %d but limits are [0 0 0] to [%d %d %d].", ix,
        iy, iz, sb->nx - 1, sb->ny - 1, sb->nz - 1);
    /*return -1;*/
  }
}

/*************************************************************************
 polygonalize_cuboid:
 In: opp: an ordered polygon object that we will create
     sb: a subdivided box
 Out: returns 1 on failure, 0 on success.  The partitions along each axis of
      the subdivided box are considered to be grid lines along which we
      subdivide the surface of the box.  Walls corresponding to these surface
      elements are created and placed into the polygon_object.
*************************************************************************/
static int polygonalize_cuboid(struct polygon_object *pop,
                               struct subdivided_box *sb) {
  struct vector3 *v;
  struct element_data *e;
  struct vertex_list *head = NULL;

  pop->n_verts = count_cuboid_vertices(sb);

  struct vector3 *vert_array =
      CHECKED_MALLOC_ARRAY(struct vector3, pop->n_verts, "cuboid vertices");
  if (vert_array == NULL)
    return 1;

  pop->n_walls = count_cuboid_elements(sb);
  pop->element =
      CHECKED_MALLOC_ARRAY(struct element_data, pop->n_walls, "cuboid walls");
  if (pop->element == NULL) {
    free(vert_array);
    vert_array = NULL;
    return 1;
  }

  /*  for (a=0;a<2;a++) for (b=0;b<2;b++) for (c=0;c<2;c++)
   * printf("%d,%d,%d->%d\n",a,b,c,vertex_at_index(sb,a,b,c)); */

  /* Set vertices and elements on X faces */
  int ii = 0;
  int bb = 0;
  int cc = 2 * (sb->nz - 1) * (sb->ny - 1);
  int b = 0;
  int c = sb->nz * sb->ny;
  for (int j = 0; j < sb->nz; j++) {
    int a = sb->ny;
    for (int i = 0; i < sb->ny; i++) {
      /*printf("Setting indices %d %d\n",b+j*a+i,c+j*a+i);*/
      v = &(vert_array[b + j * a + i]);
      v->x = sb->x[0];
      v->y = sb->y[i];
      v->z = sb->z[j];
      v = &(vert_array[c + j * a + i]);
      v->x = sb->x[sb->nx - 1];
      v->y = sb->y[i];
      v->z = sb->z[j];

      if (i > 0 && j > 0) {
        e = &(pop->element[bb + ii]);
        e->vertex_index[0] = vertex_at_index(sb, 0, i - 1, j - 1);
        e->vertex_index[2] = vertex_at_index(sb, 0, i, j - 1);
        e->vertex_index[1] = vertex_at_index(sb, 0, i - 1, j);
        e = &(pop->element[bb + ii + 1]);
        e->vertex_index[0] = vertex_at_index(sb, 0, i, j);
        e->vertex_index[1] = vertex_at_index(sb, 0, i, j - 1);
        e->vertex_index[2] = vertex_at_index(sb, 0, i - 1, j);
        e = &(pop->element[cc + ii]);
        e->vertex_index[0] = vertex_at_index(sb, sb->nx - 1, i - 1, j - 1);
        e->vertex_index[1] = vertex_at_index(sb, sb->nx - 1, i, j - 1);
        e->vertex_index[2] = vertex_at_index(sb, sb->nx - 1, i - 1, j);
        e = &(pop->element[cc + ii + 1]);
        e->vertex_index[0] = vertex_at_index(sb, sb->nx - 1, i, j);
        e->vertex_index[2] = vertex_at_index(sb, sb->nx - 1, i, j - 1);
        e->vertex_index[1] = vertex_at_index(sb, sb->nx - 1, i - 1, j);
        /*printf("Setting elements %d %d %d %d of
         * %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,pop->n_walls);*/

        ii += 2;
      }
    }
  }

  /* Set vertices and elements on Y faces */
  bb = ii;
  cc = bb + 2 * (sb->nx - 1) * (sb->nz - 1);
  b = 2 * sb->nz * sb->ny;
  c = b + sb->nz * (sb->nx - 2);
  for (int j = 0; j < sb->nz; j++) {
    int a = sb->nx - 2;
    for (int i = 1; i < sb->nx; i++) {
      if (i < sb->nx - 1) {
        /*printf("Setting indices %d %d of
         * %d\n",b+j*a+(i-1),c+j*a+(i-1),pop->n_verts);*/
        v = &(vert_array[b + j * a + (i - 1)]);
        v->x = sb->x[i];
        v->y = sb->y[0];
        v->z = sb->z[j];
        v = &(vert_array[c + j * a + (i - 1)]);
        v->x = sb->x[i];
        v->y = sb->y[sb->ny - 1];
        v->z = sb->z[j];
      }

      if (j > 0) {
        e = &(pop->element[bb + ii]);
        e->vertex_index[0] = vertex_at_index(sb, i - 1, 0, j - 1);
        e->vertex_index[1] = vertex_at_index(sb, i, 0, j - 1);
        e->vertex_index[2] = vertex_at_index(sb, i - 1, 0, j);
        e = &(pop->element[bb + ii + 1]);
        e->vertex_index[0] = vertex_at_index(sb, i, 0, j);
        e->vertex_index[2] = vertex_at_index(sb, i, 0, j - 1);
        e->vertex_index[1] = vertex_at_index(sb, i - 1, 0, j);
        e = &(pop->element[cc + ii]);
        e->vertex_index[0] = vertex_at_index(sb, i - 1, sb->ny - 1, j - 1);
        e->vertex_index[2] = vertex_at_index(sb, i, sb->ny - 1, j - 1);
        e->vertex_index[1] = vertex_at_index(sb, i - 1, sb->ny - 1, j);
        e = &(pop->element[cc + ii + 1]);
        e->vertex_index[0] = vertex_at_index(sb, i, sb->ny - 1, j);
        e->vertex_index[1] = vertex_at_index(sb, i, sb->ny - 1, j - 1);
        e->vertex_index[2] = vertex_at_index(sb, i - 1, sb->ny - 1, j);
        /*printf("Setting elements %d %d %d %d of
         * %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,pop->n_walls);*/

        ii += 2;
      }
    }
  }

  /* Set vertices and elements on Z faces */
  bb = ii;
  cc = bb + 2 * (sb->nx - 1) * (sb->ny - 1);
  b = 2 * sb->nz * sb->ny + 2 * (sb->nx - 2) * sb->nz;
  c = b + (sb->nx - 2) * (sb->ny - 2);
  for (int j = 1; j < sb->ny; j++) {
    int a = sb->nx - 2;
    for (int i = 1; i < sb->nx; i++) {
      if (i < sb->nx - 1 && j < sb->ny - 1) {
        /*printf("Setting indices %d %d of
         * %d\n",b+(j-1)*a+(i-1),c+(j-1)*a+(i-1),pop->n_verts);*/
        v = &(vert_array[b + (j - 1) * a + (i - 1)]);
        v->x = sb->x[i];
        v->y = sb->y[j];
        v->z = sb->z[0];
        v = &(vert_array[c + (j - 1) * a + (i - 1)]);
        v->x = sb->x[i];
        v->y = sb->y[j];
        v->z = sb->z[sb->nz - 1];
      }

      e = &(pop->element[bb + ii]);
      e->vertex_index[0] = vertex_at_index(sb, i - 1, j - 1, 0);
      e->vertex_index[2] = vertex_at_index(sb, i, j - 1, 0);
      e->vertex_index[1] = vertex_at_index(sb, i - 1, j, 0);
      e = &(pop->element[bb + ii + 1]);
      e->vertex_index[0] = vertex_at_index(sb, i, j, 0);
      e->vertex_index[1] = vertex_at_index(sb, i, j - 1, 0);
      e->vertex_index[2] = vertex_at_index(sb, i - 1, j, 0);
      e = &(pop->element[cc + ii]);
      e->vertex_index[0] = vertex_at_index(sb, i - 1, j - 1, sb->nz - 1);
      e->vertex_index[1] = vertex_at_index(sb, i, j - 1, sb->nz - 1);
      e->vertex_index[2] = vertex_at_index(sb, i - 1, j, sb->nz - 1);
      e = &(pop->element[cc + ii + 1]);
      e->vertex_index[0] = vertex_at_index(sb, i, j, sb->nz - 1);
      e->vertex_index[2] = vertex_at_index(sb, i, j - 1, sb->nz - 1);
      e->vertex_index[1] = vertex_at_index(sb, i - 1, j, sb->nz - 1);

      /*printf("Setting elements %d %d %d %d of
       * %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,pop->n_walls);*/

      ii += 2;
    }
  }

  /* build the head node of the linked list "pop->parsed_vertices" */

  struct vertex_list *vlp = CHECKED_MALLOC_STRUCT(
      struct vertex_list, "vertex_list");
  if (vlp == NULL) {
    free(vert_array);
    return 1;
  }
  vlp->vertex = CHECKED_MALLOC_STRUCT(struct vector3, "vertex");
  if (vlp->vertex == NULL) {
    free(vert_array);
    free(vlp);
    return 1;
  }
  memcpy(vlp->vertex, &vert_array[0], sizeof(struct vector3));
  vlp->next = head;
  head = vlp;
  struct vertex_list *tail = head;

  /* build other nodes of the linked list "pop->parsed_vertices" */
  int error_code = 0;
  for (int i = 1; i < pop->n_verts; i++) {
    vlp = CHECKED_MALLOC_STRUCT(struct vertex_list, "vertex_list");
    if (vlp == NULL) {
      error_code = 1;
      break;
    }
    vlp->vertex = CHECKED_MALLOC_STRUCT(struct vector3, "vertex");
    if (vlp->vertex == NULL) {
      free(vlp);
      error_code = 1;
      break;
    }
    memcpy(vlp->vertex, &vert_array[i], sizeof(struct vector3));
    vlp->next = tail->next;
    tail->next = vlp;
    tail = tail->next;
  }
  if (error_code == 1) {
    free_vertex_list(head);
    if (vert_array != NULL)
      free(vert_array);
    vert_array = NULL;
    return 1;
  }
  pop->parsed_vertices = head;

#ifdef DEBUG
  printf("BOX has vertices:\n");
  for (int i = 0; i < pop->n_verts; i++)
    printf("  %.5e %.5e %.5e\n", vert_array[i].x, vert_array[i].y,
           vert_array[i].z);
  printf("BOX has walls:\n");
  for (int i = 0; i < pop->n_walls; i++)
    printf("  %d %d %d\n", pop->element[i].vertex_index[0],
           pop->element[i].vertex_index[1], pop->element[i].vertex_index[2]);
  printf("\n");
#endif

  if (vert_array != NULL)
    free(vert_array);
  vert_array = NULL;

  return 0;
}

/*************************************************************************
 mdl_triangulate_box_object:
    Finalizes the polygonal structure of the box, normalizing all regions.

 In:  parse_state: parser state
      box_sym: symbol for the box object
      pop: polygon object for the box
      box_aspect_ratio: aspect ratio for the box
 Out: 0 on success, 1 on failure.  Box is polygonalized and regions normalized.
*************************************************************************/
int mdl_triangulate_box_object(struct mdlparse_vars *parse_state,
                               struct sym_entry *box_sym,
                               struct polygon_object *pop,
                               double box_aspect_ratio) {
  struct object *objp = (struct object *)box_sym->value;

  if (box_aspect_ratio >= 2.0) {
    if (reaspect_cuboid(parse_state, pop->sb, box_aspect_ratio)) {
      mdlerror(parse_state, "Error setting up box geometry");
      return 1;
    }
  }
  for (struct region_list *rlp = objp->regions; rlp != NULL; rlp = rlp->next) {
    if (mdl_normalize_elements(parse_state, rlp->reg, 0))
      return 1;
  }
  if (polygonalize_cuboid(pop, pop->sb)) {
    mdlerror(parse_state, "Could not turn box object into polygons");
    return 1;
  } else if (parse_state->vol->notify->box_triangulation == NOTIFY_FULL) {
    mcell_log("Box object %s converted into %d polygons.", box_sym->name,
              pop->n_walls);
  }

  const unsigned int n_walls = pop->n_walls;
  pop->side_removed = new_bit_array(n_walls);
  if (pop->side_removed == NULL) {
    mcell_allocfailed("Failed to allocate a box object removed side bitmask.");
    /*return 1;*/
  }
  set_all_bits(pop->side_removed, 0);
  return 0;
}

/**************************************************************************
 mdl_check_diffusion_constant:
    Check that the specified diffusion constant is valid, correcting it if
    appropriate.

 In: parse_state: parser state
     d: pointer to the diffusion constant
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_check_diffusion_constant(struct mdlparse_vars *parse_state, double *d) {
  if (parse_state->vol->notify->neg_diffusion == WARN_COPE) {
    if (*d < 0)
      *d = 0.0;
  } else if (parse_state->vol->notify->neg_diffusion == WARN_WARN) {
    if (*d < 0.0) {
      mcell_warn(
          "negative diffusion constant found, setting to zero and continuing.");
      *d = 0.0;
    }
  } else {
    if (*d < 0.0) {
      mdlerror(parse_state, "diffusion constants should be zero or positive.");
      return 1;
    }
  }
  return 0;
}

/*************************************************************************
 report_diffusion_distances:
    Helper function to print average diffusion distances per species.

 In:  spec: the species
      time_unit:
      length_unit:
      lvl:
 Out: Nothing.
*************************************************************************/
static void report_diffusion_distances(struct mcell_species_spec *spec,
                                       double time_unit, double length_unit,
                                       int lvl) {

  double l_perp_bar = 0;
  double l_perp_rms = 0;
  double l_r_bar = 0;
  double l_r_rms = 0;

  if (spec->custom_time_step == 1.0) {
    /* Theoretical average diffusion distances for the molecule need to
     * distinguish between 2D and 3D molecules for computing l_r_bar and
     * friends */

    // Volume molecule
    if ((spec->is_2d) == 0) {
      l_perp_bar = sqrt(4 * 1.0e8 * spec->D * time_unit / MY_PI);
      l_perp_rms = sqrt(2 * 1.0e8 * spec->D * time_unit);
      l_r_bar = 2 * l_perp_bar;
      l_r_rms = sqrt(6 * 1.0e8 * spec->D * time_unit);
    }
    // Surface molecule
    else {
      l_r_bar = sqrt(MY_PI * 1.0e8 * spec->D * time_unit);
    }

    if (lvl == NOTIFY_FULL) {
      mcell_log(
          "MCell: Theoretical average diffusion distances for molecule %s:\n"
          "\tl_r_bar = %.9g microns\n"
          "\tl_r_rms = %.9g microns\n"
          "\tl_perp_bar = %.9g microns\n"
          "\tl_perp_rms = %.9g microns",
          spec->name, l_r_bar, l_r_rms, l_perp_bar, l_perp_rms);
    } else if (lvl == NOTIFY_BRIEF) {
      mcell_log("  l_r_bar=%.9g um for %s", l_r_bar, spec->name);
    }
  } else {
    if (lvl == NOTIFY_FULL) {
      /* the size of the length unit depends on if the molecule is
       * 2D or 3D; the values for step length simply follow from
       * converting space_step = sqrt(4Dt) into the 2D/3D expression
       * for l_r_bar */
      double step_length = 0.0;
      if ((spec->is_2d) == 0) {
        step_length = length_unit * spec->space_step * 2.0 / sqrt(MY_PI);
      } else {
        step_length = length_unit * spec->space_step * sqrt(MY_PI) / 2.0;
      }

      mcell_log("MCell: Theoretical average diffusion time for molecule %s:\n"
                "\tl_r_bar fixed at %.9g microns\n"
                "\tPosition update every %.3e seconds (%.3g timesteps)",
                spec->name, step_length, spec->custom_time_step * time_unit,
                spec->custom_time_step);
    } else if (lvl == NOTIFY_BRIEF) {
      mcell_log("  delta t=%.3g timesteps for %s", spec->custom_time_step,
                spec->name);
    }
  }
}

// TODO: Remove by merging with mdl_print_species_summaries or at least
// eliminate redundancies
void mdl_print_species_summary(MCELL_STATE *state,
                               struct mcell_species_spec *species) {
  if (state->procnum == 0) {
    if (state->notify->diffusion_constants == NOTIFY_BRIEF) {
      mcell_log("Defining molecule with the following diffusion constant:");
    }
    report_diffusion_distances(species, state->time_unit, state->length_unit,
                               state->notify->diffusion_constants);
    no_printf("Molecule %s defined with D = %g\n", species->name, species->D);
    free(species->name);
    free(species);
  }
}

/*************************************************************************
 mdl_print_species_summaries:
    Finish the creation of a series of molecules, undoing any state changes we
    made during the creation of the molecules.  Presently, this just means
    "print the diffusion distances report".

 In: state: the simulation state
 Out: A report is printed to the file handle.
*************************************************************************/
void
mdl_print_species_summaries(struct volume *state,
                            struct parse_mcell_species_list_item *spec_items) {
  if (state->procnum == 0) {
    if (state->notify->diffusion_constants == NOTIFY_BRIEF) {
      mcell_log("Defining molecules with the following theoretical average "
                "diffusion distances:");
    }
    struct parse_mcell_species_list_item *spec_item;
    for (spec_item = spec_items; spec_item != NULL;
         spec_item = spec_item->next) {
      struct mcell_species_spec *spec = spec_item->spec;
      report_diffusion_distances(spec, state->time_unit, state->length_unit,
                                 state->notify->diffusion_constants);
      no_printf("Molecule %s defined with D = %g\n", spec->name, spec->D);
    }
    struct parse_mcell_species_list_item *next;
    for (spec_item = spec_items; NULL != spec_item; spec_item = next) {
      next = spec_item->next;
      free(spec_item->spec->name);
      free(spec_item->spec);
      free(spec_item);
    }
    if (state->notify->diffusion_constants == NOTIFY_BRIEF) {
      mcell_log_raw("\n");
    }
  }
}

/*************************************************************************
 mdl_add_to_species_list:
    Helper function to add a species to a species list.

 In:  species_list_mem:
      list: the list of species
      spec: the species to be added
 Out: 0 on success
*************************************************************************/
int mdl_add_to_species_list(struct parse_mcell_species_list *list,
                            struct mcell_species_spec *spec) {
  struct parse_mcell_species_list_item *spec_item =
      CHECKED_MALLOC_STRUCT(struct parse_mcell_species_list_item,
                            "struct parse_mcell_species_list_item");

  spec_item->spec = spec;
  spec_item->next = NULL;
  if (list->species_count == 0) {
    list->species_tail = list->species_head = spec_item;
    list->species_count = 1;
  } else {
    list->species_tail = list->species_tail->next = spec_item;
    ++list->species_count;
  }

  return 0;
}

/**************************************************************************
 mdl_start_release_site:
    Start parsing the innards of a release site.

 In: parse_state: parser state
     symp: symbol for the release site
     shape: shape for the release site
 Out: 0 on success, 1 on failure
 NOTE: This is just a thin wrapper around start_release_site
**************************************************************************/
int mdl_start_release_site(struct mdlparse_vars *parse_state,
                           struct sym_entry *symp, int shape) {
  struct object *obj_ptr = NULL;
  if (mcell_start_release_site(parse_state->vol, symp, &obj_ptr)) {
    return 1;
  }

  parse_state->current_release_site = obj_ptr->contents;
  if (obj_ptr->contents == NULL) {
    return 1;
  }

  parse_state->current_release_site->release_shape = (int8_t)shape;
  return 0;
}

/**************************************************************************
 mdl_finish_release_site:
    Finish parsing the innards of a release site.

 In: parse_state: parser state
     symp: symbol for the release site
 Out: the object, on success, or NULL on failure
 NOTE: This is just a thin wrapper around finish_release_site
**************************************************************************/
struct object *mdl_finish_release_site(struct mdlparse_vars *parse_state,
                                       struct sym_entry *symp) {
  struct object *objp_new = NULL;
  if (mcell_finish_release_site(symp, &objp_new)) {
    mcell_error_nodie("Failed to create release site %s", symp->name);
    return NULL;
  }

  if (objp_new == NULL) {
    return NULL;
  }
  parse_state->current_release_site = NULL;
  return objp_new;
}

/**************************************************************************
 mdl_is_release_site_valid:
    Validate a release site.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
 Out: 0 if it is valid, 1 if not
 NOTE: This is just a thin wrapper around is_release_site_valid
**************************************************************************/
/*
//XXX: Remove this but port over error messages first
int mdl_is_release_site_valid(struct mdlparse_vars *parse_state,
                              struct release_site_obj *rel_site_obj_ptr) {
  switch (is_release_site_valid(rel_site_obj_ptr)) {
  case 2:
    mdlerror(parse_state,
             "Must specify molecule to release using MOLECULE=molecule_name.");
    return 1;
  case 3:
    mdlerror_fmt(parse_state,
                 "Cannot release surface class '%s' from release site",
                 rel_site_obj_ptr->mol_type->sym->name);
    return 1;
  case 4:
    mdlerror(parse_state, "CONCENTRATION may only be used with molecules that "
                          "can diffuse in 3D.\n"
                          "  Use DENSITY for molecules diffusing in 2D.");
    return 1;
  case 5:
    mdlerror(parse_state,
             "DENSITY may only be used with molecules that can diffuse in 2D.\n"
             "  Use CONCENTRATION for molecules diffusing in 3D.");
    return 1;
  case 6:
    mdlerror(parse_state, "Release site is missing location.");
    return 1;
  }
  return 0;
}
*/

/*************************************************************************
 mdl_check_release_regions:

  In: parse_state: parser state
      rel_eval: an release evaluator (set operations applied to regions)
      parent: the object that owns this release evaluator
      instance: the root object that begins the instance tree
  Out: 0 if all regions refer to instanced objects or to a common ancestor of
       the object with the evaluator, meaning that the object can be found.  1
       if any referred-to region cannot be found.
 NOTE: This is just a thin wrapper around check_release_regions
*************************************************************************/
static int mdl_check_release_regions(struct mdlparse_vars *parse_state,
                                     struct release_evaluator *rel_eval,
                                     struct object *parent,
                                     struct object *instance) {
  switch (check_release_regions(rel_eval, parent, instance)) {
  case 1:
    return 1;
  case 2:
    mdlerror(parse_state,
             "Region neither instanced nor grouped with release site.");
    return 1;
  case 3:
    mdlerror(parse_state, "Region not grouped with release site.");
    return 1;
  }

  return 0;
}

/**************************************************************************
 mdl_set_release_site_geometry_region:
    Set the geometry for a particular release site to be a region expression.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
     obj_ptr: the object representing this release site
     rel_eval:   the release evaluator representing the region of release
 Out: 0 on success, 1 on failure
**************************************************************************/
int
mdl_set_release_site_geometry_region(struct mdlparse_vars *parse_state,
                                     struct release_site_obj *rel_site_obj_ptr,
                                     struct object *obj_ptr,
                                     struct release_evaluator *rel_eval) {
  switch (mcell_set_release_site_geometry_region(
      parse_state->vol, rel_site_obj_ptr, obj_ptr, rel_eval)) {
  case 1:
    return 1;
  case 2:
    mdlerror(
        parse_state,
        "Trying to release on a region that the release site cannot see!\n  "
        "Try grouping the release site and the corresponding geometry with an "
        "OBJECT.");
    return 1;
  }

  return 0;
}

/**************************************************************************
 mdl_set_release_site_geometry_object:
    Set the geometry for a particular release site to be an entire object.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
     obj_ptr: the object upon which to release
 Out: 0 on success, 1 on failure
**************************************************************************/
int
mdl_set_release_site_geometry_object(struct mdlparse_vars *parse_state,
                                     struct release_site_obj *rel_site_obj_ptr,
                                     struct object *obj_ptr) {
  if ((obj_ptr->object_type == META_OBJ) ||
      (obj_ptr->object_type == REL_SITE_OBJ)) {
    mdlerror(
        parse_state,
        "only BOX or POLYGON_LIST objects may be assigned to the SHAPE keyword "
        "in the RELEASE_SITE definition. Metaobjects or release objects are "
        "not allowed here.");
    return 1;
  }

  char *obj_name = obj_ptr->sym->name;
  char *region_name = CHECKED_SPRINTF("%s,ALL", obj_name);
  if (region_name == NULL) {
    return 1;
  }
  struct sym_entry *sym_ptr;
  if (((sym_ptr = retrieve_sym(
      region_name, parse_state->vol->reg_sym_table)) == NULL) ||
      sym_ptr->count == 0) {
    mdlerror_fmt(parse_state, "Undefined region: %s", region_name);
    free(region_name);
    return 1;
  }
  free(region_name);

  struct release_evaluator *rel_eval =
      CHECKED_MALLOC_STRUCT(struct release_evaluator, "release site on region");
  if (rel_eval == NULL) {
    return 1;
  }

  rel_eval->op = REXP_NO_OP | REXP_LEFT_REGION;
  rel_eval->left = sym_ptr->value;
  rel_eval->right = NULL;

  ((struct region *)rel_eval->left)->flags |= COUNT_CONTENTS;

  rel_site_obj_ptr->release_shape = SHAPE_REGION;
  parse_state->vol->place_waypoints_flag = 1;

  if (mdl_check_release_regions(parse_state, rel_eval, obj_ptr,
                                parse_state->vol->root_instance)) {
    mdlerror(
        parse_state,
        "Trying to release on a region that the release site cannot see!\n  "
        "Try grouping the release site and the corresponding geometry with an "
        "OBJECT.");
    free(rel_eval);
    return 1;
  }

  struct release_region_data *rel_reg_data = CHECKED_MALLOC_STRUCT(
      struct release_region_data, "release site on region");
  if (rel_reg_data == NULL) {
    mdlerror(parse_state,
             "Out of memory while trying to create release site on region");
    free(rel_eval);
    return 1;
  }

  rel_reg_data->n_walls_included = -1; /* Indicates uninitialized state */
  rel_reg_data->cum_area_list = NULL;
  rel_reg_data->wall_index = NULL;
  rel_reg_data->obj_index = NULL;
  rel_reg_data->n_objects = -1;
  rel_reg_data->owners = NULL;
  rel_reg_data->in_release = NULL;
  rel_reg_data->self = parse_state->current_object;
  rel_reg_data->expression = rel_eval;
  rel_site_obj_ptr->region_data = rel_reg_data;

  return 0;
}

/**************************************************************************
 mdl_check_valid_molecule_release:
    Check that a particular molecule type is valid for inclusion in a release
    site.  Checks that orientations are present if required, and absent if
    forbidden, and that we aren't trying to release a surface class.

 In: parse_state: parser state
     mol_type: molecule species and (optional) orientation for release
 Out: 0 on success, 1 on failure
**************************************************************************/
static int mdl_check_valid_molecule_release(struct mdlparse_vars *parse_state,
                                            struct mcell_species *mol_type) {
  static const char *EXTRA_ORIENT_MSG =
      "surface orientation not specified for released surface molecule\n"
      "(use ; or ', or ,' for random orientation)";
  static const char *MISSING_ORIENT_MSG =
      "orientation not used for released volume molecule";

  struct species *mol = (struct species *)mol_type->mol_type->value;
  if (mol->flags & ON_GRID) {
    if (!mol_type->orient_set) {
      if (parse_state->vol->notify->missed_surf_orient == WARN_ERROR) {
        mdlerror_fmt(parse_state, "Error: %s", EXTRA_ORIENT_MSG);
        return 1;
      } else if (parse_state->vol->notify->missed_surf_orient == WARN_WARN) {
        mdlerror_fmt(parse_state, "Warning: %s", EXTRA_ORIENT_MSG);
      }
    }
  } else if ((mol->flags & NOT_FREE) == 0) {
    if (mol_type->orient_set) {
      if (parse_state->vol->notify->useless_vol_orient == WARN_ERROR) {
        mdlerror_fmt(parse_state, "Error: %s", MISSING_ORIENT_MSG);
        return 1;
      } else if (parse_state->vol->notify->useless_vol_orient == WARN_WARN) {
        mdlerror_fmt(parse_state, "Warning: %s", MISSING_ORIENT_MSG);
      }
    }
  } else {
    mdlerror(parse_state,
             "cannot release a surface class instead of a molecule.");
    return 1;
  }

  return 0;
}

/**************************************************************************
 mdl_set_release_site_molecule:
    Set the molecule to be released from this release site.

 In: parse_state: parser state
     rel_site_obj_ptr: release site object
     mol_type: molecule species and (optional) orientation for release
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_release_site_molecule(struct mdlparse_vars *parse_state,
                                  struct release_site_obj *rel_site_obj_ptr,
                                  struct mcell_species *mol_type) {
  // Store molecule information
  rel_site_obj_ptr->mol_type = (struct species *)mol_type->mol_type->value;
  if ((rel_site_obj_ptr->mol_type->flags & NOT_FREE) == 0) {
    if (rel_site_obj_ptr->release_shape == SHAPE_REGION)
      parse_state->vol->place_waypoints_flag = 1;
  } else {
    if (rel_site_obj_ptr->release_shape != SHAPE_REGION &&
        rel_site_obj_ptr->release_shape != SHAPE_LIST) {
      mdlerror_fmt(parse_state, "The release site '%s' is a geometric release "
                                "site, and may not be used to \n"
                                "  release the surface molecule '%s'.  Surface "
                                "molecule release sites must \n"
                                "  be either LIST or region release sites.",
                   rel_site_obj_ptr->name, mol_type->mol_type->name);
      return 1;
    }
  }
  rel_site_obj_ptr->orientation = mol_type->orient;

  /* Now, validate molecule information */
  return mdl_check_valid_molecule_release(parse_state, mol_type);
}

/**************************************************************************
 mdl_set_release_site_diameter:
    Set the diameter of a release site.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
     diam: the desired diameter of this release site
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_release_site_diameter(struct mdlparse_vars *parse_state,
                                  struct release_site_obj *rel_site_obj_ptr,
                                  double diam) {
  diam *= parse_state->vol->r_length_unit;

  rel_site_obj_ptr->diameter =
      CHECKED_MALLOC_STRUCT(struct vector3, "release site diameter");
  if (rel_site_obj_ptr->diameter == NULL) {
    return 1;
  }
  rel_site_obj_ptr->diameter->x = diam;
  rel_site_obj_ptr->diameter->y = diam;
  rel_site_obj_ptr->diameter->z = diam;
  return 0;
}

/**************************************************************************
 mdl_set_release_site_diameter_array:
    Set the diameter of the release site along the X, Y, and Z axes.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
     n_diams: dimensionality of the diameters array (should be 3)
     diams: list containing X, Y, and Z diameters for release site
     factor: factor to scale diameter -- 2.0 if diameters are actually radii,
             1.0 for actual diameters
 Out: 0 on success, 1 on failure
**************************************************************************/
int
mdl_set_release_site_diameter_array(struct mdlparse_vars *parse_state,
                                    struct release_site_obj *rel_site_obj_ptr,
                                    int n_diams, struct num_expr_list *diams,
                                    double factor) {
  factor *= parse_state->vol->r_length_unit;

  if (rel_site_obj_ptr->release_shape == SHAPE_LIST) {
    mdlerror(parse_state, "Release list diameters must be single valued.");
    return 1;
  }

  if (n_diams != 3) {
    mdlerror(parse_state, "Three dimensional value required");
    return 1;
  }

  rel_site_obj_ptr->diameter =
      CHECKED_MALLOC_STRUCT(struct vector3, "release site diameter");
  if (rel_site_obj_ptr->diameter == NULL) {
    return 1;
  }
  rel_site_obj_ptr->diameter->x = diams->value * factor;
  rel_site_obj_ptr->diameter->y = diams->next->value * factor;
  rel_site_obj_ptr->diameter->z = diams->next->next->value * factor;
  return 0;
}

/**************************************************************************
 mdl_set_release_site_diameter_var:
    Set the diameters of the release site along the X, Y, and Z axes from a
    variable, either scalar or vector.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
     factor: factor to scale diameter -- 2.0 if diameters are actually radii,
             1.0 for actual diameters
     symp: the variable from which to set

 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_release_site_diameter_var(struct mdlparse_vars *parse_state,
                                      struct release_site_obj *rel_site_obj_ptr,
                                      double factor, struct sym_entry *symp) {
  struct num_expr_list *expr_list_ptr;
  int count = 0;
  rel_site_obj_ptr->diameter =
      CHECKED_MALLOC_STRUCT(struct vector3, "release site diameter");
  if (rel_site_obj_ptr->diameter == NULL) {
    return 1;
  }

  switch (symp->sym_type) {
  case DBL:
    if (mdl_set_release_site_diameter(parse_state, rel_site_obj_ptr,
                                      *(double *)symp->value * factor)) {
      return 1;
    }
    break;

  case ARRAY:
    // Count up to 4 elements -- that's all we need to count to know if it's
    // valid
    for (expr_list_ptr = (struct num_expr_list *)symp->value;
         expr_list_ptr != NULL && count < 4;
         ++count, expr_list_ptr = expr_list_ptr->next)
      ;
    expr_list_ptr = (struct num_expr_list *)symp->value;
    if (mdl_set_release_site_diameter_array(parse_state, rel_site_obj_ptr,
                                            count, expr_list_ptr, factor)) {
      return 1;
    }
    break;

  default:
    mdlerror(parse_state,
             "Diameter must either be a number or a 3-valued vector.");
    return 1;
  }

  return 0;
}

/**************************************************************************
 mdl_set_release_site_periodic_box:

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object
     periodic_box: the periodic box that we want to release molecules into

 Out: 0 on success
**************************************************************************/
int mdl_set_release_site_periodic_box(struct mdlparse_vars *parse_state,
                                      struct release_site_obj *rel_site_obj_ptr,
                                      struct vector3 *periodic_box) {
  rel_site_obj_ptr->periodic_box->x = (int16_t)periodic_box->x;                                    
  rel_site_obj_ptr->periodic_box->y = (int16_t)periodic_box->y;                                    
  rel_site_obj_ptr->periodic_box->z = (int16_t)periodic_box->z;                                    
  return 0;
}

/**************************************************************************
 mdl_set_release_site_probability:
    Set the release probability for a release site.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
     prob: the release probability
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_release_site_probability(struct mdlparse_vars *parse_state,
                                     struct release_site_obj *rel_site_obj_ptr,
                                     double prob) {
  if (!distinguishable(rel_site_obj_ptr->release_prob,
                       MAGIC_PATTERN_PROBABILITY, EPS_C)) {
    mdlerror(parse_state,
             "Ignoring release probability for reaction-triggered releases.");
  } else {
    rel_site_obj_ptr->release_prob = prob;
    if (rel_site_obj_ptr->release_prob < 0) {
      mdlerror(parse_state, "Release probability cannot be less than 0.");
      return 1;
    }
    if (rel_site_obj_ptr->release_prob > 1) {
      mdlerror(parse_state, "Release probability cannot be greater than 1.");
      return 1;
    }
  }

  return 0;
}

/**************************************************************************
 mdl_set_release_site_pattern:
    Set the release pattern to be used by a particular release site.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
     pattern: the release pattern
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_release_site_pattern(struct mdlparse_vars *parse_state,
                                 struct release_site_obj *rel_site_obj_ptr,
                                 struct sym_entry *pattern) {
  rel_site_obj_ptr->pattern = (struct release_pattern *)pattern->value;

  // Careful!  We've put a rxn_pathname into the "pattern" pointer!
  if (pattern->sym_type == RXPN) {
    if (rel_site_obj_ptr->release_prob != 1.0) {
      mdlerror(parse_state,
               "Ignoring release probability for reaction-triggered releases.");
    }
    // Magic number indicating a reaction-triggered release
    rel_site_obj_ptr->release_prob = MAGIC_PATTERN_PROBABILITY;
  }

  return 0;
}

/**************************************************************************
 mdl_set_release_site_molecule_positions:
    Set the molecule positions for a LIST release.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site object to validate
     list: list of release_single_molecule structs
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_release_site_molecule_positions(
    struct mdlparse_vars *parse_state,
    struct release_site_obj *rel_site_obj_ptr,
    struct release_single_molecule_list *list) {
  if (rel_site_obj_ptr->release_shape != SHAPE_LIST) {
    mdlerror(parse_state,
             "You must use the LIST shape to specify molecule positions in a "
             "release.");
    return 1;
  }

  struct release_single_molecule *rsm;
  if (rel_site_obj_ptr->mol_list == NULL) {
    rel_site_obj_ptr->mol_list = list->rsm_head;
  } else {
    for (rsm = rel_site_obj_ptr->mol_list; rsm->next != NULL; rsm = rsm->next)
      ;
    rsm->next = list->rsm_head;
  }
  rel_site_obj_ptr->release_number += list->rsm_count;
  return 0;
}

/**************************************************************************
 mdl_new_release_single_molecule:
    Create a mew single molecule release position for a LIST release site.

 In: parse_state: parser state
     mol_type: molecule type and optional orientation
     pos: 3D position in the world
 Out: molecule release description, or NULL if an error occurred
**************************************************************************/
struct release_single_molecule *
mdl_new_release_single_molecule(struct mdlparse_vars *parse_state,
                                struct mcell_species *mol_type,
                                struct vector3 *pos) {
  struct vector3 temp_v3;
  memcpy(&temp_v3, pos, sizeof(struct vector3));
  free(pos);

  struct release_single_molecule *rsm = CHECKED_MALLOC_STRUCT(
      struct release_single_molecule, "release site molecule position");
  if (rsm == NULL) {
    mdlerror(parse_state, "Out of memory reading molecule positions");
    return NULL;
  }

  rsm->orient = mol_type->orient;
  rsm->loc.x = temp_v3.x * parse_state->vol->r_length_unit;
  rsm->loc.y = temp_v3.y * parse_state->vol->r_length_unit;
  rsm->loc.z = temp_v3.z * parse_state->vol->r_length_unit;
  rsm->mol_type = (struct species *)(mol_type->mol_type->value);
  rsm->next = NULL;

  if (mdl_check_valid_molecule_release(parse_state, mol_type)) {
    free(rsm);
    return NULL;
  }

  return rsm;
}

/**************************************************************************
 mdl_set_release_site_concentration:
    Set a release quantity from this release site based on a fixed
    concentration within the release-site's area.

 In: parse_state: parser state
     rel_site_obj_ptr: the release site
     conc: concentration for release
 Out: 0 on success, 1 on failure.  release site object is updated
 NOTE: This is just a thin wrapper around set_release_site_concentration
**************************************************************************/
int
mdl_set_release_site_concentration(struct mdlparse_vars *parse_state,
                                   struct release_site_obj *rel_site_obj_ptr,
                                   double conc) {
  if (set_release_site_concentration(rel_site_obj_ptr, conc)) {
    mdlerror_fmt(parse_state,
                 "Release site '%s' is a spherical shell; concentration-based "
                 "release is not supported on a spherical shell",
                 rel_site_obj_ptr->name);
    return 1;
  }
  return 0;
}

/**************************************************************************
 mdl_vertex_list_singleton:
    Set an item to be the sole element of a vertex list.

 In: head: the list
     item: ite item
 Out: none.  list is updated
**************************************************************************/
void mdl_vertex_list_singleton(struct vertex_list_head *head,
                               struct vertex_list *item) {
  item->next = NULL;
  head->vertex_tail = head->vertex_head = item;
  head->vertex_count = 1;
}

/**************************************************************************
 mdl_add_vertex_to_list:
    Append a vertex to a list.

 In: head: the list
     item: ite item
 Out: none.  list is updated
**************************************************************************/
void mdl_add_vertex_to_list(struct vertex_list_head *head,
                            struct vertex_list *item) {
  item->next = NULL;
  head->vertex_tail = head->vertex_tail->next = item;
  ++head->vertex_count;
}

/**************************************************************************
 mdl_new_vertex_list_item:
    Allocate an item for a vertex list.

 In: vertex: this vertex
     normal: surface normal at this vertex, or NULL
 Out: the vertex list item, or NULL if an error occurred
**************************************************************************/
struct vertex_list *mdl_new_vertex_list_item(struct vector3 *vertex) {

  struct vertex_list *vlp =
      CHECKED_MALLOC_STRUCT(struct vertex_list, "vertices");
  if (vlp == NULL)
    return NULL;
  vlp->vertex = vertex;
  vlp->next = NULL;
  return vlp;
}

/**************************************************************************
 mdl_element_connection_list_singleton:
    Set an item to be the sole element of an element connection list.

 In: head: the list
     item: ite item
 Out: none.  list is updated
**************************************************************************/
void
mdl_element_connection_list_singleton(struct element_connection_list_head *head,
                                      struct element_connection_list *item) {
  item->next = NULL;
  head->connection_tail = head->connection_head = item;
  head->connection_count = 1;
}

/**************************************************************************
 mdl_add_element_connection_to_list:
    Append an element connection to a list.

 In: head: the list
     item: ite item
 Out: none.  list is updated
**************************************************************************/
void
mdl_add_element_connection_to_list(struct element_connection_list_head *head,
                                   struct element_connection_list *item) {
  item->next = NULL;
  head->connection_tail = head->connection_tail->next = item;
  ++head->connection_count;
}

/**************************************************************************
 mdl_new_element_connection:
    Create an element connection (essentially a triplet of vertex indices).

 In: parse_state: parser state
     indices: the element connections
 Out: the list, or NULL if an error occurred
**************************************************************************/
struct element_connection_list *
mdl_new_element_connection(struct mdlparse_vars *parse_state,
                           struct num_expr_list_head *indices) {
  if (indices->value_count != 3) {
    mdlerror(parse_state,
             "Non-triangular element found in polygon list object");
    return NULL;
  }

  struct element_connection_list *eclp = CHECKED_MALLOC_STRUCT(
      struct element_connection_list, "polygon element commections");
  if (eclp == NULL)
    return NULL;

  eclp->indices = CHECKED_MALLOC_ARRAY(int, 3, "polygon element connections");
  if (eclp->indices == NULL) {
    free(eclp);
    return NULL;
  }
  eclp->indices[0] = (int)indices->value_head->value;
  eclp->indices[1] = (int)indices->value_head->next->value;
  eclp->indices[2] = (int)indices->value_tail->value;
  eclp->n_verts = indices->value_count;
  eclp->next = NULL;

  if (!indices->shared)
    mcell_free_numeric_list(indices->value_head);
  return eclp;
}

/**************************************************************************
 mdl_new_tet_element_connection:
    Create a tetrahedral element connection (essentially a quadruplet of vertex
    indices).

 In: parse_state: parser state
     indices: the element connections
 Out: the list, or NULL if an error occurred
**************************************************************************/
struct element_connection_list *
mdl_new_tet_element_connection(struct mdlparse_vars *parse_state,
                               struct num_expr_list_head *indices) {
  if (indices->value_count != 4) {
    mdlerror(parse_state, "Non-tetrahedron element found in voxel list object");
    return NULL;
  }

  struct element_connection_list *eclp = CHECKED_MALLOC_STRUCT(
      struct element_connection_list, "polygon element commections");
  if (eclp == NULL)
    return NULL;

  eclp->indices = CHECKED_MALLOC_ARRAY(int, 4, "polygon element connections");
  if (eclp->indices == NULL) {
    free(eclp);
    return NULL;
  }
  eclp->indices[0] = (int)indices->value_head->value;
  eclp->indices[1] = (int)indices->value_head->next->value;
  eclp->indices[2] = (int)indices->value_head->next->next->value;
  eclp->indices[3] = (int)indices->value_tail->value;
  eclp->n_verts = indices->value_count;
  eclp->next = NULL;

  if (!indices->shared)
    mcell_free_numeric_list(indices->value_head);
  return eclp;
}

/**************************************************************************
 mdl_new_polygon_list:
    Create a new polygon list object.

 In: parse_state: parser state
     sym: symbol for this polygon list
     n_vertices: count of vertices
     vertices: list of vertices
     n_connections: count of walls
     connections: list of walls
 Out: polygon object, or NULL if there was an error
**************************************************************************/
struct object *
mdl_new_polygon_list(struct mdlparse_vars *parse_state, char *obj_name,
                     int n_vertices, struct vertex_list *vertices,
                     int n_connections,
                     struct element_connection_list *connections) {
  struct object_creation obj_creation;
  obj_creation.object_name_list = parse_state->object_name_list;
  obj_creation.object_name_list_end = parse_state->object_name_list_end;
  obj_creation.current_object = parse_state->current_object;

  if (parse_state->vol->disable_polygon_objects) {
    mdlerror(
        parse_state,
        "When using dynamic geometries, polygon objects should only be "
        "defined/instantiated through the dynamic geometry file.");
  }
  int error_code = 0;
  struct object *obj_ptr =
      start_object(parse_state->vol, &obj_creation, obj_name, &error_code);
  if (error_code == 1) {
    mdlerror_fmt(parse_state,"Object '%s' is already defined", obj_name);
  }
  else if (error_code == 2) {
    mdlerror_fmt(parse_state, "Out of memory while creating object: %s",
                 obj_name);
  }

  struct polygon_object *poly_obj_ptr =
      new_polygon_list(parse_state->vol, obj_ptr, n_vertices, vertices,
                       n_connections, connections);

  parse_state->object_name_list = obj_creation.object_name_list;
  parse_state->object_name_list_end = obj_creation.object_name_list_end;
  parse_state->current_object = obj_ptr;

  parse_state->allow_patches = 0;
  parse_state->current_polygon = poly_obj_ptr;

  return obj_ptr;
}

/**************************************************************************
 mdl_finish_polygon_list:
    Finalize the polygon list, cleaning up any state updates that were made
    when we started creating the polygon.

 In: parse_state: parser state
     symp: symbol for the completed polygon
 Out: 1 on failure, 0 on success
**************************************************************************/
int mdl_finish_polygon_list(struct mdlparse_vars *parse_state,
                            struct object *obj_ptr) {
  struct object_creation obj_creation;
  obj_creation.object_name_list_end = parse_state->object_name_list_end;

  int error_code = 0;
  if (finish_polygon_list(obj_ptr, &obj_creation)) {
    error_code = 1;
  }
  parse_state->object_name_list_end = obj_creation.object_name_list_end;
  parse_state->current_object = parse_state->current_object->parent;
  parse_state->current_polygon = NULL;

  return error_code;
}

/**************************************************************************
 allocate_voxel_object:
    Create a new voxel object.

 In:  Nothing
 Out: voxel object, or NULL if allocation fails
**************************************************************************/
static struct voxel_object *allocate_voxel_object() {

  struct voxel_object *vop;
  if ((vop = CHECKED_MALLOC_STRUCT(struct voxel_object, "voxel list object")) ==
      NULL)
    return NULL;
  vop->vertex = NULL;
  vop->element = NULL;
  vop->neighbor = NULL;
  vop->n_verts = 0;
  vop->n_voxels = 0;

  return vop;
}

/**************************************************************************
 mdl_new_voxel_list:
    Create a new voxel list object.

 In: parse_state: parser state
     sym:  the symbol for this voxel list
     n_vertices: count of vertices in this object
     vertices: list of vertices for this object
     n_connections: count of tetrahedra in this object
     connections: list of tetrahedra
 Out: voxel object, or NULL if there is an error
**************************************************************************/
struct voxel_object *
mdl_new_voxel_list(struct mdlparse_vars *parse_state, struct sym_entry *sym,
                   int n_vertices, struct vertex_list *vertices,
                   int n_connections,
                   struct element_connection_list *connections) {
  struct tet_element_data *tedp;

  struct object *objp = (struct object *)sym->value;
  struct voxel_object *vop = allocate_voxel_object();
  if (vop == NULL)
    goto failure;

  objp->object_type = VOXEL_OBJ;
  objp->contents = vop;

  vop->n_voxels = n_connections;
  vop->n_verts = n_vertices;

  /* Allocate vertices */
  if ((vop->vertex = CHECKED_MALLOC_ARRAY(
           struct vector3, vop->n_verts, "voxel list object vertices")) == NULL)
    goto failure;

  /* Populate vertices */
  for (int i = 0; i < vop->n_verts; i++) {
    struct vertex_list *vlp_temp = vertices;
    vop->vertex[i].x = vertices->vertex->x;
    vop->vertex[i].y = vertices->vertex->y;
    vop->vertex[i].z = vertices->vertex->z;
    free(vertices->vertex);
    vertices = vertices->next;
    free(vlp_temp);
  }

  /* Allocate tetrahedra */
  if ((tedp = CHECKED_MALLOC_ARRAY(struct tet_element_data, vop->n_voxels,
                                   "voxel list object tetrahedra")) == NULL)
    goto failure;
  vop->element = tedp;

  /* Copy in tetrahedra */
  for (int i = 0; i < vop->n_voxels; i++) {
    if (connections->n_verts != 4) {
      mdlerror(parse_state, "All voxels must have four vertices.");
      goto failure;
    }

    struct element_connection_list *eclp_temp = connections;
    memcpy(tedp[i].vertex_index, connections->indices, 4 * sizeof(int));
    connections = connections->next;
    free(eclp_temp);
  }
  return vop;

failure:
  free_vertex_list(vertices);
  free_connection_list(connections);
  if (vop) {
    if (vop->element)
      free(vop->element);
    if (vop->vertex)
      free(vop->vertex);
    free(vop);
  }
  return NULL;
}

struct polygon_object *mdl_create_periodic_box(
    struct mdlparse_vars *parse_state,
    struct vector3 *llf,
    struct vector3 *urb,
    bool isPeriodicX,
    bool isPeriodicY,
    bool isPeriodicZ) {

  struct polygon_object *pop;
  struct region *rp;

  char *name_tmp = "PERIODIC_BOX_OBJ";
  int name_len = strlen(name_tmp) + 1;
  char *name = (char*)malloc(name_len * sizeof(char));
  strcpy(name, name_tmp);

  struct sym_entry *sym = mdl_start_object(parse_state, name);
  struct object *objp = (struct object *)sym->value;

  /* Allocate polygon object */
  pop = allocate_polygon_object("box object");
  if (pop == NULL) {
    free(llf);
    free(urb);
    return NULL;
  }
  objp->object_type = BOX_OBJ;
  objp->contents = pop;

  /* Create object default region on box object: */
  if ((rp = mdl_create_region(parse_state, objp, "ALL")) == NULL) {
    free(pop);
    free(llf);
    free(urb);
    return NULL;
  }
  if ((rp->element_list_head = new_element_list(ALL_SIDES, ALL_SIDES)) ==
      NULL) {
    free(pop);
    free(llf);
    free(urb);
    return NULL;
  }

  /* Scale corners to internal units */
  llf->x *= parse_state->vol->r_length_unit;
  llf->y *= parse_state->vol->r_length_unit;
  llf->z *= parse_state->vol->r_length_unit;
  urb->x *= parse_state->vol->r_length_unit;
  urb->y *= parse_state->vol->r_length_unit;
  urb->z *= parse_state->vol->r_length_unit;

  /* Initialize our subdivided box */
  pop->sb = init_cuboid(parse_state, llf, urb);
  free(llf);
  free(urb);
  if (pop->sb == NULL) {
    free(pop);
    return NULL;
  }

  // mark box as periodic or not
  objp->periodic_x = isPeriodicX;
  objp->periodic_y = isPeriodicY;
  objp->periodic_z = isPeriodicZ;

  parse_state->allow_patches = 1;
  parse_state->current_polygon = pop;

  mdl_triangulate_box_object(parse_state, sym, parse_state->current_polygon, 0.0);

  return pop;
}

int mdl_finish_periodic_box(struct mdlparse_vars *parse_state) {
  struct sym_entry *symp = retrieve_sym("PERIODIC_BOX_OBJ", parse_state->vol->obj_sym_table);
  struct object *objp = (struct object *)symp->value;
  remove_gaps_from_regions(objp);
  objp->n_walls = parse_state->current_polygon->n_walls;
  objp->n_verts = parse_state->current_polygon->n_verts;
  if (check_degenerate_polygon_list(objp)) {
    parse_state->current_polygon = NULL;
    return 1;
  }

  parse_state->current_polygon = NULL;

  // This next bit is a little strange. We are essentially, creating a meta
  // object that contains an instance of the periodic box object. The meta
  // object will be added to the root instance, like a user would do with their
  // "Scene" or "World" objects.
  parse_state->current_object = parse_state->vol->root_instance;

  // Create meta object
  char *meta_name_tmp = "PERIODIC_BOX_META";
  int meta_name_len = strlen(meta_name_tmp) + 1;
  char *meta_name = (char*)malloc(meta_name_len * sizeof(char));
  strcpy(meta_name, meta_name_tmp);

  struct sym_entry *meta_sym = mdl_start_object(parse_state, meta_name);
  struct object *meta_objp = (struct object *)meta_sym->value;

  meta_objp->object_type = META_OBJ;

  // Create instance of PERIODIC_BOX_OBJECT
  char *inst_name_tmp = "PERIODIC_BOX_INSTANT";
  int inst_name_len = strlen(inst_name_tmp) + 1;
  char *inst_name = (char*)malloc(inst_name_len * sizeof(char));
  strcpy(inst_name, inst_name_tmp);

  struct sym_entry *inst_sym = mdl_start_object(parse_state, inst_name);
  struct object *inst_objp = (struct object *)inst_sym->value;

  mdl_deep_copy_object(parse_state, inst_objp, objp);

  // Finish instance object
  mdl_finish_object(parse_state);
  add_child_objects(meta_objp, inst_objp, inst_objp);
  parse_state->vol->periodic_box_obj = inst_objp;
  // Finish meta object
  mdl_finish_object(parse_state);
  add_child_objects(parse_state->vol->root_instance, meta_objp, meta_objp);
  parse_state->current_object = parse_state->vol->root_object;
 
  return 0;
}

/**************************************************************************
 mdl_new_box_object:
    Create a new box object, with particular corners.

 In: parse_state: parser state
     sym:  symbol for this box object
     llf:  lower left front corner
     urb:  upper right back corner
 Out: polygon object for this box, or NULL if there's an error
**************************************************************************/
struct polygon_object *mdl_new_box_object(struct mdlparse_vars *parse_state,
                                          struct sym_entry *sym,
                                          struct vector3 *llf,
                                          struct vector3 *urb) {
  struct polygon_object *pop;
  struct region *rp;
  struct object *objp = (struct object *)sym->value;


  /* Allocate polygon object */
  pop = allocate_polygon_object("box object");
  if (pop == NULL) {
    free(llf);
    free(urb);
    return NULL;
  }
  objp->object_type = BOX_OBJ;
  objp->contents = pop;

  /* Create object default region on box object: */
  if ((rp = mdl_create_region(parse_state, objp, "ALL")) == NULL) {
    free(pop);
    free(llf);
    free(urb);
    return NULL;
  }
  if ((rp->element_list_head = new_element_list(ALL_SIDES, ALL_SIDES)) ==
      NULL) {
    free(pop);
    free(llf);
    free(urb);
    return NULL;
  }

  /* Scale corners to internal units */
  llf->x *= parse_state->vol->r_length_unit;
  llf->y *= parse_state->vol->r_length_unit;
  llf->z *= parse_state->vol->r_length_unit;
  urb->x *= parse_state->vol->r_length_unit;
  urb->y *= parse_state->vol->r_length_unit;
  urb->z *= parse_state->vol->r_length_unit;

  /* Initialize our subdivided box */
  pop->sb = init_cuboid(parse_state, llf, urb);
  free(llf);
  free(urb);
  if (pop->sb == NULL) {
    free(pop);
    return NULL;
  }

  parse_state->allow_patches = 1;
  parse_state->current_polygon = pop;
  return pop;
}

/**************************************************************************
 mdl_finish_box_object:
    Finalize the box object, cleaning up any state updates that were made when
    we started creating the box.

 In: parse_state: parser state
     symp: symbol for the completed box
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_finish_box_object(struct mdlparse_vars *parse_state,
                          struct sym_entry *symp) {
  struct object *objp = (struct object *)symp->value;
  remove_gaps_from_regions(objp);
  objp->n_walls = parse_state->current_polygon->n_walls;
  objp->n_verts = parse_state->current_polygon->n_verts;
  if (check_degenerate_polygon_list(objp)) {
    parse_state->current_polygon = NULL;
    return 1;
  }

  parse_state->current_polygon = NULL;
  return 0;
}

/**************************************************************************
 mdl_create_region:
    Create a named region on an object.

 In: parse_state: parser state
     objp: object upon which to create a region
     name: region name to create
 Out: region object, or NULL if there was an error (region already exists or
      allocation failed)
**************************************************************************/
struct region *mdl_create_region(struct mdlparse_vars *parse_state,
                                 struct object *objp, char *name) {
  struct region *rp;
  struct region_list *rlp;
  no_printf("Creating new region: %s\n", name);
  if ((rp = mdl_make_new_region(parse_state, objp->sym->name, name)) == NULL)
    return NULL;
  if ((rlp = CHECKED_MALLOC_STRUCT(struct region_list, "region list")) ==
      NULL) {
    mdlerror_fmt(parse_state, "Out of memory while creating object region '%s'",
                 rp->sym->name);
    return NULL;
  }
  rp->region_last_name = name;
  rp->parent = objp;
  char *region_name = CHECKED_SPRINTF("%s,%s", objp->sym->name, name);
  if (!mcell_check_for_region(region_name, objp)) {
    rlp->reg = rp;
    rlp->next = objp->regions;
    objp->regions = rlp;
    objp->num_regions++;
  }
  else {
    free(rlp);
  }
  free(region_name);
  return rp;
}

/**************************************************************************
 mdl_get_region:
    Get a region on an object, creating it if it does not exist yet.

 In: parse_state: parser state
     objp: object upon which to create
     name: region to get
 Out: region, or NULL if allocation fails
**************************************************************************/
struct region *mdl_get_region(struct mdlparse_vars *parse_state,
                              struct object *objp, char *name) {
  struct sym_entry *reg_sym;
  char *region_name;
  struct region *rp;

  region_name = CHECKED_SPRINTF("%s,%s", objp->sym->name, name);
  if (region_name == NULL)
    return NULL;

  reg_sym = retrieve_sym(region_name, parse_state->vol->reg_sym_table);
  free(region_name);

  if (reg_sym == NULL)
    rp = mdl_create_region(parse_state, objp, name);
  else
    rp = (struct region *)reg_sym->value;

  return rp;
}

/**************************************************************************
 mdl_start_existing_obj_region_def:
    Begin construction of a region on an existing object.

 In: parse_state: parser state
     obj_symp: symbol of object upon which to create region
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_start_existing_obj_region_def(struct mdlparse_vars *parse_state,
                                      struct sym_entry *obj_symp) {
  struct object *objp = (struct object *)obj_symp->value;
  if (objp->object_type != BOX_OBJ && objp->object_type != POLY_OBJ) {
    mdlerror_fmt(parse_state, "Cannot define region on non-surface object: %s",
                 obj_symp->name);
    return 1;
  }
  parse_state->current_polygon = objp->contents;
  parse_state->current_object = objp;
  parse_state->allow_patches = 0;
  return 0;
}

/**************************************************************************
 mdl_add_elements_to_list:
    Append an element to an element list.

 In: list: list for element
     head: first element to add
     tail: last element to add
 Out: none.  list is updated
**************************************************************************/
void mdl_add_elements_to_list(struct element_list_head *list,
                              struct element_list *head,
                              struct element_list *tail) {
  tail->next = NULL;
  list->elml_tail->next = head;
  list->elml_tail = tail;
}

/**************************************************************************
 mdl_set_elements_to_exclude:
    Marks elements as being excluded, rather than included.  This is done by
    setting the "special" pointer on each element to point to the element
    itself.  The resultant "special" pointer, then, points to the wrong type of
    object, but the pointer itself is used as a signal to the rest of the code.

 In: els:  elements to set to exclude
 Out: none.  list is updated
**************************************************************************/
void mdl_set_elements_to_exclude(struct element_list *els) {
  /* HACK: els->special actually points to an element_list.  Don't dereference
   * it now. */
  for (; els != NULL; els = els->next)
    els->special = (struct element_special *)els;
}

/**************************************************************************
 mdl_new_element_side:
    Create a new element list for a region description based on a side name.

 In: parse_state: parser state
     side: side name constant (ALL_SIDES, etc.)
 Out: element list, or NULL if allocation fails
**************************************************************************/
struct element_list *mdl_new_element_side(struct mdlparse_vars *parse_state,
                                          unsigned int side) {
  unsigned int begin, end;
  if (side == ALL_SIDES &&
      parse_state->current_object->object_type == POLY_OBJ) {
    begin = 0;
    end = parse_state->current_polygon->n_walls - 1;
  } else if (parse_state->current_object->object_type == POLY_OBJ) {
    mdlerror(parse_state,
             "Illegal reference to polygon list element by side-name");
    return NULL;
  } else {
    begin = side;
    end = side;
  }
  return new_element_list(begin, end);
}

/**************************************************************************
 mdl_new_element_previous_region:
    Create a new element list for a "previous region" include/exclude
    statement.

 In: parse_state: parser state
     objp: object containing referent region
     rp_container: region for whom we're creating this element list
     name_region_referent: name of referent region
     exclude: 1 if we're excluding, 0 if including
 Out: element list, or NULL if an error occrs
**************************************************************************/
struct element_list *mdl_new_element_previous_region(
    struct mdlparse_vars *parse_state, struct object *objp,
    struct region *rp_container, char *name_region_referent, int exclude) {
  struct sym_entry *stp;
  char *full_reg_name = NULL;
  struct element_list *elmlp = NULL;

  /* Create element list */
  elmlp = new_element_list(0, 0);
  if (elmlp == NULL)
    goto failure;

  /* Create "special" element description */
  elmlp->special =
      CHECKED_MALLOC_STRUCT(struct element_special, "region element");
  if (elmlp->special == NULL)
    goto failure;
  elmlp->special->exclude = (byte)exclude;

  /* Create referent region full name */
  full_reg_name =
      CHECKED_SPRINTF("%s,%s", objp->sym->name, name_region_referent);
  if (full_reg_name == NULL)
    goto failure;

  /* Look up region or die */
  stp = retrieve_sym(full_reg_name, parse_state->vol->reg_sym_table);
  if (stp == NULL) {
    mdlerror_fmt(parse_state, "Undefined region: %s", full_reg_name);
    goto failure;
  }
  free(full_reg_name);
  full_reg_name = NULL;

  /* Store referent region */
  elmlp->special->referent = (struct region *)stp->value;
  if (elmlp->special->referent == rp_container) {
    mdlerror_fmt(parse_state,
                 "Self-referential region include.  No paradoxes, please.");
    goto failure;
  }

  free(name_region_referent);
  return elmlp;

failure:
  free(name_region_referent);
  if (full_reg_name)
    free(full_reg_name);
  if (elmlp) {
    if (elmlp->special)
      free(elmlp->special);
    free(elmlp);
  }
  return NULL;
}

/**************************************************************************
 mdl_new_element_patch:
    Allocate a new region element list item for an include/exclude PATCH
    statement.

 In: parse_state: parser state
     polygon: polygon upon which we're making a patch
     llf: first corner of patch
     urb: second corner of patch
     exclude: 1 if we're excluding, 0 if including
 Out: element list, or NULL if an error occrs
**************************************************************************/
struct element_list *mdl_new_element_patch(struct mdlparse_vars *parse_state,
                                           struct polygon_object *poly,
                                           struct vector3 *llf,
                                           struct vector3 *urb, int exclude) {
  if (parse_state->current_object->object_type != BOX_OBJ) {
    mdlerror(
        parse_state,
        "INCLUDE_PATCH and EXCLUDE_PATCH may only be used on a BOX object.");
    return NULL;
  }

  if (!parse_state->allow_patches) {
    mdlerror(
        parse_state,
        "Cannot create PATCH on a BOX outside of the original declaration.");
    return NULL;
  }

  struct element_list *elmlp = new_element_list(0, 0);
  if (elmlp == NULL)
    goto failure;

  /* Allocate special element description */
  elmlp->special =
      CHECKED_MALLOC_STRUCT(struct element_special, "region element");
  if (elmlp->special == NULL)
    goto failure;
  elmlp->special->referent = NULL;
  elmlp->special->exclude = (byte)exclude;

  /* Convert to internal units */
  llf->x *= parse_state->vol->r_length_unit;
  llf->y *= parse_state->vol->r_length_unit;
  llf->z *= parse_state->vol->r_length_unit;
  urb->x *= parse_state->vol->r_length_unit;
  urb->y *= parse_state->vol->r_length_unit;
  urb->z *= parse_state->vol->r_length_unit;
  memcpy(&(elmlp->special->corner1), llf, sizeof(struct vector3));
  memcpy(&(elmlp->special->corner2), urb, sizeof(struct vector3));

  /* Refine the cuboid's mesh to accomodate the new patch */
  if (refine_cuboid(parse_state, llf, urb, poly->sb,
                    parse_state->vol->grid_density))
    goto failure;

  free(llf);
  free(urb);
  return elmlp;

failure:
  free(llf);
  free(urb);
  if (elmlp) {
    free(elmlp->special);
    free(elmlp);
  }
  return NULL;
}

/**************************************************************************
 mdl_set_region_elements:
    Set the elements for a region, normalizing the region if it's on a polygon
    list object.

 In: parse_state: parser state
     rgn:  region to receive elements
     elements: elements comprising region
     normalize_now: flag indicating whether to normalize right now
 Out: symbol for new pathway, or NULL if an error occurred
**************************************************************************/
int mdl_set_region_elements(struct mdlparse_vars *parse_state,
                            struct region *rgn, struct element_list *elements,
                            int normalize_now) {
  rgn->element_list_head = elements;
  if (normalize_now)
    return mdl_normalize_elements(parse_state, rgn, 0);
  else
    return 0;
}

/**************************************************************************
 mdl_new_rxn_pathname:
    Create a new named reaction pathway name structure.

 In: parse_state: parser state
     name: name for new named pathway
 Out: symbol for new pathway, or NULL if an error occurred
**************************************************************************/
struct sym_entry *mdl_new_rxn_pathname(struct mdlparse_vars *parse_state,
                                       char *name) {
  if ((retrieve_sym(name, parse_state->vol->rxpn_sym_table)) != NULL) {
    mdlerror_fmt(parse_state, "Named reaction pathway already defined: %s",
                 name);
    free(name);
    return NULL;
  } else if ((retrieve_sym(name, parse_state->vol->mol_sym_table)) != NULL) {
    mdlerror_fmt(parse_state,
                 "Named reaction pathway already defined as a molecule: %s",
                 name);
    free(name);
    return NULL;
  }

  struct sym_entry *symp =
      store_sym(name, RXPN, parse_state->vol->rxpn_sym_table, NULL);
  if (symp == NULL) {
    mdlerror_fmt(parse_state, "Out of memory while creating reaction name: %s",
                 name);
    free(name);
    return NULL;
  }
  free(name);
  return symp;
}

/**************************************************************************
 mdl_add_surf_mol_to_region:
    Adds an surface molecule (or list of surface molecules) to a region.  These
    surface molecules will be placed on the surface at initialization time.

 In: rgn:  the region
     lst:  a list of surface molecules to place
 Out: none.  list is merged into region
**************************************************************************/
void mdl_add_surf_mol_to_region(struct region *rgn, struct sm_dat_list *lst) {
  lst->sm_tail->next = rgn->sm_dat_head;
  rgn->sm_dat_head = lst->sm_head;
}

/**************************************************************************
 mdl_set_region_surface_class:
    Set the surface class of this region, possibly inheriting the viz_value.

 In: parse_state: parser state
     rgn:  the region
     scsymp: symbol for the surface class
 Out: none.  region is updated.
**************************************************************************/
void mdl_set_region_surface_class(struct mdlparse_vars *parse_state,
                                  struct region *rgn,
                                  struct sym_entry *scsymp) {
  if (rgn->surf_class != NULL) {
    mdlerror(parse_state, "ATTENTION: region definition allows only one "
                          "SURFACE_CLASS statement.");
  }
  rgn->surf_class = (struct species *)scsymp->value;
}

/*************************************************************************
 * Reaction output
 *************************************************************************/

/**************************************************************************
 mdl_new_output_set:
    Populate an output set.

 In: parse_state: parser state
     os:   output set
     col_head: head of linked list of output columns
     file_flags: file creation disposition
     outfile_name: file output name
 Out: output set, or NULL if an error occurs
**************************************************************************/
struct output_set *mdl_populate_output_set(struct mdlparse_vars *parse_state,
                                           char *comment, int exact_time,
                                           struct output_column *col_head,
                                           int file_flags, char *outfile_name) {
  if ((parse_state->count_flags & (TRIGGER_PRESENT | COUNT_PRESENT)) ==
      (TRIGGER_PRESENT | COUNT_PRESENT)) {
    mdlerror(parse_state,
             "Cannot mix TRIGGER and COUNT statements.  Use separate files.");
    return NULL;
  }

  struct output_set *os =
      mcell_create_new_output_set(comment, exact_time,
                                  col_head, file_flags, outfile_name);

  return os;
}

/**************************************************************************
 mdl_add_reaction_output_block_to_world:
    Construct and add an output block to the world.

 In: parse_state: parser state
     buffer_size: size of output buffer for this block
     otimes: output timing information for this block
     osets: output sets for this block
 Out: 0 on success, 1 on failure; world is updated with new output block
**************************************************************************/
int mdl_add_reaction_output_block_to_world(struct mdlparse_vars *parse_state,
                                           int buffer_size,
                                           struct output_times_inlist *otimes,
                                           struct output_set_list *osets) {

  return mcell_add_reaction_output_block(parse_state->vol, osets, buffer_size,
                                         otimes);
}

/**************************************************************************
 mdl_join_oexpr_tree:
    Joins two subtrees into a reaction data output expression tree, with a
    specified operation.

 In: parse_state: parser state
     left: left subtree
     right: right subtree
     oper: specified operation
 Out: joined output expression, or NULL if an error occurs
**************************************************************************/
struct output_expression *mdl_join_oexpr_tree(struct mdlparse_vars *parse_state,
                                              struct output_expression *left,
                                              struct output_expression *right,
                                              char oper) {
  struct output_expression *joined;
  struct output_expression *leaf, *new_oe, *up;
  int first_leaf = 1;

  joined = NULL;
  if (left->oper == ',' && (right == NULL || right->oper == ',')) {
    mdlerror(parse_state, "Can't do math on multiple wildcard expressions");
    return NULL;
  }

  if (left->oper != ',' && (right == NULL || right->oper != ',')) {
    joined = new_output_expr(parse_state->vol->oexpr_mem);
    if (joined == NULL)
      return NULL;

    joined->left = (void *)left;
    joined->right = (void *)right;
    joined->oper = oper;
    left->up = joined;
    if (right != NULL)
      right->up = joined;

    learn_oexpr_flags(joined);
    if (joined->expr_flags & OEXPR_TYPE_CONST)
      eval_oexpr_tree(joined, 0);

    return joined;
  } else if (left->oper == ',') {
    for (leaf = first_oexpr_tree(left); leaf != NULL;
         leaf = next_oexpr_tree(leaf)) {
      if (first_leaf) {
        new_oe = right;
        first_leaf = 0;
      } else if (right != NULL) {
        new_oe = dupl_oexpr_tree(right, parse_state->vol->oexpr_mem);
        if (new_oe == NULL)
          return NULL;
      } else
        new_oe = NULL;

      up = leaf->up;
      joined = mdl_join_oexpr_tree(parse_state, leaf, new_oe, oper);
      if (joined == NULL)
        return NULL;
      joined->up = up;
      if (leaf == up->left)
        up->left = joined;
      else
        up->right = joined;
      if (joined->expr_flags & OEXPR_TYPE_CONST)
        eval_oexpr_tree(joined, 0);
      learn_oexpr_flags(up);
      leaf = joined;
    }
    return left;
  } else /* right->oper==',' */
  {
    for (leaf = first_oexpr_tree(right); leaf != NULL;
         leaf = next_oexpr_tree(leaf)) {
      if (first_leaf) {
        new_oe = left;
        first_leaf = 0;
      } else {
        new_oe = dupl_oexpr_tree(left, parse_state->vol->oexpr_mem);
        if (new_oe == NULL)
          return NULL;
      }
      up = leaf->up;
      joined = mdl_join_oexpr_tree(parse_state, new_oe, leaf, oper);
      if (joined == NULL)
        return NULL;
      joined->up = up;
      if (leaf == up->left)
        up->left = joined;
      else
        up->right = joined;
      if (joined->expr_flags & OEXPR_TYPE_CONST)
        eval_oexpr_tree(joined, 0);
      learn_oexpr_flags(up);
      leaf = joined;
    }

    return right;
  }

  return NULL; /* Should never get here */
}

/**************************************************************************
 mdl_sum_expression:
    Convert an output expression tree into a summation.

 In: expr: expression to convert
 Out: modified expression
**************************************************************************/
struct output_expression *mdl_sum_oexpr(struct output_expression *expr) {
  oexpr_flood_convert(expr, ',', '+');
  eval_oexpr_tree(expr, 0);
  return expr;
}

/**************************************************************************
 mdl_new_oexpr_constant:
    Creates a constant output expression for reaction data output.

 In: parse_state: parser state
     value: the value of the constantj
 Out: the output expression, or NULL if allocation fails
**************************************************************************/
struct output_expression *
mdl_new_oexpr_constant(struct mdlparse_vars *parse_state, double value) {
  struct output_expression *oe = new_output_expr(parse_state->vol->oexpr_mem);
  if (oe == NULL) {
    mdlerror(parse_state, "Out of memory creating output expression");
    return NULL;
  }
  oe->expr_flags = OEXPR_TYPE_DBL | OEXPR_TYPE_CONST;
  oe->value = value;
  oe->oper = '=';
  return oe;
}



/**************************************************************************
 mdl_count_syntax_periodic_1:

    Generates a reaction data output expression from the first count syntax
    form (simple molecule, unquoted, no orientation) within a certain
    periodic box.

    example:

      COUNT[foo, region, [periodic_box_x, periodic_box_y, periodic_box_z]]

 In: parse_state: parser state
     what: symbol representing the molecule type
     where: symbol representing the count location (or NULL for WORLD)
     periodicBox: what box are we counting in?
     hit_spec: what are we counting?
     count_flags: is this a count or a trigger?
 Out: 0 on success, 1 on failure
**************************************************************************/
struct output_expression *mdl_count_syntax_periodic_1(
  struct mdlparse_vars *parse_state,
  struct sym_entry *what,
  struct sym_entry *where,
  struct vector3 *periodicBox,
  int hit_spec,
  int count_flags) {

  if (parse_state->vol->periodic_traditional) {
    mdlerror(parse_state, "Counting in virtual periodic boxes is invalid if PERIODIC_TRADITIONAL is TRUE");
  }
  // cannot combine world counting with periodic box since world means everything
  if (where == NULL) {
    mdlerror(parse_state, "Invalid combination of WORLD with periodic box counting");
  }

  byte report_flags = 0;
  if (count_flags & TRIGGER_PRESENT)
    report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED)
    report_flags |= REPORT_ENCLOSED;

  if (what->sym_type == MOL) {
    if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING)
      report_flags |= REPORT_CONTENTS;
    else
      report_flags |= (hit_spec & REPORT_TYPE_MASK);
  } else {
    report_flags |= REPORT_RXNS;
    if ((hit_spec & REPORT_TYPE_MASK) != REPORT_NOTHING) {
      mdlerror_fmt(parse_state,
                   "Invalid counting options used with reaction pathway %s",
                   what->name);
      return NULL;
    }
  }

  /* extract image for periodic image we would like to count */
  struct periodic_image *img = CHECKED_MALLOC_STRUCT(struct periodic_image,
    "periodic image descriptor");
  img->x = (int16_t)periodicBox->x;
  img->y = (int16_t)periodicBox->y;
  img->z = (int16_t)periodicBox->z;
  free(periodicBox);

  struct output_request *orq;
  if ((orq = mcell_new_output_request(parse_state->vol, what, ORIENT_NOT_SET,
                                      where, img, report_flags)) == NULL)
    return NULL;
  orq->next = parse_state->vol->output_request_head;
  parse_state->vol->output_request_head = orq;
  return orq->requester;
}



/**************************************************************************
 mdl_count_syntax_1:
    Generates a reaction data output expression from the first count syntax
    form (simple molecule, unquoted, no orientation).

    example:

      COUNT(foo,WORLD)

 In: parse_state: parser state
     what: symbol representing the molecule type
     where: symbol representing the count location (or NULL for WORLD)
     hit_spec: what are we counting?
     count_flags: is this a count or a trigger?
 Out: 0 on success, 1 on failure
**************************************************************************/
struct output_expression *mdl_count_syntax_1(struct mdlparse_vars *parse_state,
                                             struct sym_entry *what,
                                             struct sym_entry *where,
                                             int hit_spec, int count_flags) {
  byte report_flags = 0;
  struct output_request *orq;
  if (where != NULL && parse_state->vol->periodic_box_obj && !(parse_state->vol->periodic_traditional)) {
    mdlerror(parse_state,
             "If PERIODIC_TRADITIONAL is FALSE, then you must specify virtual counting box.\n"
             "(e.g. COUNT,vm,Scene.box,[1,0,0]).");
  }
  if (where == NULL) {
    report_flags = REPORT_WORLD;
    if (hit_spec != REPORT_NOTHING) {
      mdlerror(parse_state,
               "Invalid combination of WORLD with other counting options");
      return NULL;
    } else if (count_flags & TRIGGER_PRESENT) {
      mdlerror(parse_state, "Invalid combination of WORLD with TRIGGER option");
      return NULL;
    }
  } else
    report_flags = 0;

  if (count_flags & TRIGGER_PRESENT)
    report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED)
    report_flags |= REPORT_ENCLOSED;

  if (what->sym_type == MOL) {
    if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING)
      report_flags |= REPORT_CONTENTS;
    else
      report_flags |= (hit_spec & REPORT_TYPE_MASK);
  } else {
    report_flags |= REPORT_RXNS;
    if ((hit_spec & REPORT_TYPE_MASK) != REPORT_NOTHING) {
      mdlerror_fmt(parse_state,
                   "Invalid counting options used with reaction pathway %s",
                   what->name);
      return NULL;
    }
  }

  if ((orq = mcell_new_output_request(parse_state->vol, what, ORIENT_NOT_SET,
                                      where, NULL, report_flags)) == NULL)
    return NULL;
  orq->next = parse_state->vol->output_request_head;
  parse_state->vol->output_request_head = orq;
  return orq->requester;
}

/**************************************************************************
 mdl_count_syntax_periodic_2:
    Generates a reaction data output expression from the second count syntax
    form (simple molecule, unquoted, orientation in braces) within a certain
    periodic box

    example:

      COUNT[foo{1}, region, [periodic_box_x, periodic_box_y, periodic_box_z]]

 In: parse_state: parser state
     mol_type: symbol representing the molecule type
     orient: orientation specified for the molecule
     where: symbol representing the count location (or NULL for WORLD)
     periodicBox: what box are we counting in?
     hit_spec: what are we counting?
     count_flags: is this a count or a trigger?
 Out: 0 on success, 1 on failure
**************************************************************************/
struct output_expression *mdl_count_syntax_periodic_2(
    struct mdlparse_vars *parse_state,
    struct sym_entry *mol_type,
    short orient,
    struct sym_entry *where,
    struct vector3 *periodicBox,
    int hit_spec,
    int count_flags) {

  if (parse_state->vol->periodic_traditional) {
    mdlerror(
      parse_state,
      "Counting in virtual periodic boxes is invalid if PERIODIC_TRADITIONAL "
      "is TRUE");
  }
  // cant combine world counting with periodic box since world means everything
  if (where == NULL) {
    mdlerror(
      parse_state,
      "Invalid combination of WORLD with periodic box counting");
  }

  byte report_flags = 0;
  if (count_flags & TRIGGER_PRESENT)
    report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED)
    report_flags |= REPORT_ENCLOSED;

  if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING)
    report_flags |= REPORT_CONTENTS;
  else
    report_flags |= (hit_spec & REPORT_TYPE_MASK);

  /* extract image for periodic image we would like to count */
  struct periodic_image *img = CHECKED_MALLOC_STRUCT(struct periodic_image,
    "periodic image descriptor");
  img->x = (int16_t)periodicBox->x;
  img->y = (int16_t)periodicBox->y;
  img->z = (int16_t)periodicBox->z;
  free(periodicBox);

  /* Grab orientation and reset orientation state in parser */
  short orientation;
  if (orient < 0)
    orientation = -1;
  else if (orient > 0)
    orientation = 1;
  else
    orientation = 0;

  struct output_request *orq;
  if ((orq = mcell_new_output_request(parse_state->vol, mol_type, orientation,
                                      where, img, report_flags)) == NULL)
    return NULL;
  orq->next = parse_state->vol->output_request_head;
  parse_state->vol->output_request_head = orq;
  return orq->requester;
}

/**************************************************************************
 mdl_count_syntax_2:
    Generates a reaction data output expression from the second count syntax
    form (simple molecule, unquoted, orientation in braces)

    example:

      COUNT(foo{1},WORLD)

 In: parse_state: parser state
     mol_type: symbol representing the molecule type
     orient: orientation specified for the molecule
     where: symbol representing the count location (or NULL for WORLD)
     hit_spec: what are we counting?
     count_flags: is this a count or a trigger?
 Out: 0 on success, 1 on failure
**************************************************************************/
struct output_expression *mdl_count_syntax_2(struct mdlparse_vars *parse_state,
                                             struct sym_entry *mol_type,
                                             short orient,
                                             struct sym_entry *where,
                                             int hit_spec, int count_flags) {
  byte report_flags = 0;
  struct output_request *orq;
  short orientation;
  if (where != NULL && parse_state->vol->periodic_box_obj && !(parse_state->vol->periodic_traditional)) {
    mdlerror(parse_state,
             "If PERIODIC_TRADITIONAL is FALSE, then you must specify virtual counting box.\n"
             "(e.g. COUNT,vm,Scene.box,[1,0,0]).");
  }
  if (where == NULL) {
    mdlerror(parse_state, "Counting of an oriented molecule in the WORLD is "
                          "not implemented.\nAn oriented molecule may only be "
                          "counted in a regions.");
    return NULL;
  } else
    report_flags = 0;

  if (count_flags & TRIGGER_PRESENT)
    report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED)
    report_flags |= REPORT_ENCLOSED;

  if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING)
    report_flags |= REPORT_CONTENTS;
  else
    report_flags |= (hit_spec & REPORT_TYPE_MASK);

  /* Grab orientation and reset orientation state in parser */
  if (orient < 0)
    orientation = -1;
  else if (orient > 0)
    orientation = 1;
  else
    orientation = 0;

  if ((orq = mcell_new_output_request(parse_state->vol, mol_type, orientation,
                                      where, NULL, report_flags)) == NULL)
    return NULL;
  orq->next = parse_state->vol->output_request_head;
  parse_state->vol->output_request_head = orq;
  return orq->requester;
}

/**************************************************************************
 mdl_get_orientation_from_string:
    Get the orientation from a molecule name+orientation string, removing the
    orientation marks from the string.  This only supports the "tick" notation.

 In: mol_string: string to parse
 Out: the orientation.  mol_string has been modified to remove the orientation
      part
**************************************************************************/
static int mdl_get_orientation_from_string(char *mol_string) {
  short orientation = 0;
  int pos = strlen(mol_string) - 1;

  if (mol_string[pos] == ';') {
    mol_string[pos] = '\0';
    return 0;
  }

  /* Peel off any apostrophes or commas at the end of the string */
  while (pos >= 0) {
    switch (mol_string[pos]) {
    case '\'':
      ++orientation;
      break;
    case ',':
      --orientation;
      break;
    default:
      mol_string[pos + 1] = '\0';
      pos = 0;
      break;
    }
    --pos;
  }

  if (orientation < 0)
    return -1;
  else if (orientation > 0)
    return 1;
  else
    return 0;
}

/**************************************************************************
 mdl_string_has_orientation:
    Check if the string looks like a molecule name+orientation string.
    Otherwise, it will be considered a wildcard.  This only supports the "tick"
    notation.

 In: mol_string: string to parse
 Out: 1 if the string has orientation, 0 if it does not
**************************************************************************/
static int mdl_string_has_orientation(char const *mol_string) {
  char last_char = mol_string[strlen(mol_string) - 1];
  return (last_char == '\'' || last_char == ',' || last_char == ';');
}

/*************************************************************************
 mdl_new_output_requests_from_list:
    Create an output expression from a list that resulted from a wildcard
    match.  The generated output requests are added to the global linked list
    of count output requests.

 In:  parse_state: parser state
      targets: linked list of what to count
      location: where are we counting?
      report_flags: what do we report
      hit_spec: how do we report
 Out: an output expression representing the requested targets
*************************************************************************/
static struct output_expression *mdl_new_output_requests_from_list(
    struct mdlparse_vars *parse_state,
    struct sym_table_list *targets,
    struct sym_entry *location,
    int report_flags,
    int hit_spec,
    struct periodic_image *img) {
  struct output_expression *oe_head = NULL, *oe_tail = NULL;
  struct output_request *or_head = NULL, *or_tail = NULL;
  int report_type;

  assert(targets != NULL);
  for (; targets != NULL; targets = targets->next) {
    if (targets->node->sym_type == MOL) {
      if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING)
        report_type = REPORT_CONTENTS;
      else
        report_type = (hit_spec & REPORT_TYPE_MASK);
    } else {
      report_type = REPORT_RXNS;
      if ((hit_spec & REPORT_TYPE_MASK) != REPORT_NOTHING) {
        mdlerror_fmt(parse_state,
                     "Invalid counting options used with reaction pathway %s",
                     targets->node->name);
        return NULL;
      }
    }

    struct output_request *orq = mcell_new_output_request(
        parse_state->vol, targets->node, ORIENT_NOT_SET, location, img,
        report_type | report_flags);
    if (orq == NULL)
      return NULL;
    struct output_expression *oe = orq->requester;

    if (oe_tail == NULL) {
      oe_head = oe_tail = oe;
      or_head = or_tail = orq;
    } else {
      or_tail->next = orq;
      or_tail = orq;

      struct output_expression *oet =
          new_output_expr(parse_state->vol->oexpr_mem);
      if (oet == NULL) {
        mdlerror(parse_state, "Out of memory storing requested counts");
        return NULL;
      }
      if (oe_tail->up == NULL) {
        oe_head = oet;
      } else {
        oet->up = oe_tail->up;
        if (oet->up->left == oe_tail)
          oet->up->left = oet;
        else
          oet->up->right = oet;
      }
      oet->left = oe_tail;
      oe_tail->up = oet;
      oet->right = oe;
      oe->up = oet;
      oet->oper = ',';
      oet->expr_flags = OEXPR_LEFT_OEXPR | OEXPR_RIGHT_OEXPR |
                        (oe->expr_flags & OEXPR_TYPE_MASK);
      oe_tail = oe;
    }
  }

  or_tail->next = parse_state->vol->output_request_head;
  parse_state->vol->output_request_head = or_head;
  return oe_head;
}

/**************************************************************************
 mdl_find_rxpns_and_mols_by_wildcard:
    Find all molecules and named reaction pathways matching a particular
    wildcard, and return them in a list.

 In: parse_state: parser state
     wildcard: the wildcard to match
 Out: the symbol list, or NULL if an error occurs.
**************************************************************************/
static struct sym_table_list *
mdl_find_rxpns_and_mols_by_wildcard(struct mdlparse_vars *parse_state,
                                    char const *wildcard) {
  struct sym_table_list *symbols = NULL, *stl;
  for (int i = 0; i < parse_state->vol->mol_sym_table->n_bins; i++) {
    for (struct sym_entry *sym_t = parse_state->vol->mol_sym_table->entries[i];
         sym_t != NULL; sym_t = sym_t->next) {
      if (is_wildcard_match((char *)wildcard, sym_t->name)) {
        stl = (struct sym_table_list *)CHECKED_MEM_GET(
            parse_state->sym_list_mem,
            "list of named reactions and molecules for counting");
        if (stl == NULL) {
          if (symbols)
            mem_put_list(parse_state->sym_list_mem, symbols);
          return NULL;
        }

        stl->node = sym_t;
        stl->next = symbols;
        symbols = stl;
      }
    }
  }
  for (int i = 0; i < parse_state->vol->rxpn_sym_table->n_bins; i++) {
    for (struct sym_entry *sym_t = parse_state->vol->rxpn_sym_table->entries[i];
         sym_t != NULL; sym_t = sym_t->next) {
      if (is_wildcard_match((char *)wildcard, sym_t->name)) {
        stl = (struct sym_table_list *)CHECKED_MEM_GET(
            parse_state->sym_list_mem,
            "list of named reactions and molecules for counting");
        if (stl == NULL) {
          if (symbols)
            mem_put_list(parse_state->sym_list_mem, symbols);
          return NULL;
        }

        stl->node = sym_t;
        stl->next = symbols;
        symbols = stl;
      }
    }
  }

  if (symbols == NULL)
    mdlerror_fmt(
        parse_state,
        "No molecules or named reactions found matching wildcard \"%s\"",
        wildcard);

  return symbols;
}

/**************************************************************************
 mdl_count_syntax_periodic_3:
    Generates a reaction data output expression from the third count syntax
    form (quoted string, possibly a wildcard, possibly an oriented molecule)
    within a certain periodic box.

    examples:

      COUNT["foo'", region, [periodic_box_x, periodic_box_y, periodic_box_z]]
      COUNT["foo*", region, [periodic_box_x, periodic_box_y, periodic_box_z]]

 In: parse_state: parser state
     what: string representing the target of this count
     orient: orientation specified for the molecule
     where: symbol representing the count location (or NULL for WORLD)
     hit_spec: what are we counting?
     count_flags: is this a count or a trigger?
 Out: 0 on success, 1 on failure
**************************************************************************/
struct output_expression *mdl_count_syntax_periodic_3(
    struct mdlparse_vars *parse_state,
    char *what,
    struct sym_entry *where,
    struct vector3 *periodicBox,
    int hit_spec,
    int count_flags) {

  if (parse_state->vol->periodic_traditional) {
    mdlerror(
      parse_state,
      "Counting in virtual periodic boxes is invalid if PERIODIC_TRADITIONAL "
      "is TRUE");
  }
  // cant combine world counting with periodic box since world means everything
  if (where == NULL) {
    mdlerror(
      parse_state,
      "Invalid combination of WORLD with periodic box counting");
  }

  byte report_flags = 0;
  if (count_flags & TRIGGER_PRESENT)
    report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED)
    report_flags |= REPORT_ENCLOSED;

  /* extract image for periodic image we would like to count */
  struct periodic_image *img = CHECKED_MALLOC_STRUCT(struct periodic_image,
    "periodic image descriptor");
  img->x = (int16_t)periodicBox->x;
  img->y = (int16_t)periodicBox->y;
  img->z = (int16_t)periodicBox->z;
  free(periodicBox);

  struct output_expression *oe;
  char *what_to_count;
  if ((what_to_count = mdl_strip_quotes(what)) == NULL)
    return NULL;

  /* Oriented molecule specified inside a string */
  if (mdl_string_has_orientation(what_to_count)) {
    struct output_request *orq;
    struct sym_entry *sp;
    short orientation;

    if (where == NULL) {
      mdlerror(parse_state, "Counting of an oriented molecule in the WORLD is "
                            "not implemented.\nAn oriented molecule may only "
                            "be counted in a region.");
      free(img);
      free(what_to_count);
      return NULL;
    }

    orientation = mdl_get_orientation_from_string(what_to_count);
    if ((sp = mdl_existing_molecule(parse_state, what_to_count)) == NULL)
      return NULL;

    if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING)
      report_flags |= REPORT_CONTENTS;
    else
      report_flags |= (hit_spec & REPORT_TYPE_MASK);

    if ((orq = mcell_new_output_request(parse_state->vol, sp, orientation,
                                        where, img, report_flags)) == NULL)
      return NULL;
    orq->next = parse_state->vol->output_request_head;
    parse_state->vol->output_request_head = orq;
    oe = orq->requester;
  }

  /* Wildcard specified inside a string */
  else {
    struct sym_table_list *stl =
        mdl_find_rxpns_and_mols_by_wildcard(parse_state, what_to_count);

    if (stl == NULL) {
      mdlerror(parse_state,
               "Wildcard matching found no matches for count output.");
      free(img);
      free(what_to_count);
      return NULL;
    }

    if (where == NULL) {
      report_flags |= REPORT_WORLD;
      if (hit_spec != REPORT_NOTHING) {
        mdlerror(parse_state,
                 "Invalid combination of WORLD with other counting options");
        free(img);
        free(what_to_count);
        return NULL;
      } else if (count_flags & TRIGGER_PRESENT) {
        mdlerror(parse_state,
                 "Invalid combination of WORLD with TRIGGER option");
        free(img);
        free(what_to_count);
        return NULL;
      }
    }

    if ((oe = mdl_new_output_requests_from_list(
        parse_state, stl, where, report_flags, hit_spec, img)) == NULL) {
      free(img);
      free(what_to_count);
      return NULL;
    }

    /* free allocated memory */
    mem_put_list(parse_state->sym_list_mem, stl);
  }

  return oe;
}

/**************************************************************************
 mdl_count_syntax_3:
    Generates a reaction data output expression from the third count syntax
    form (quoted string, possibly a wildcard, possibly an oriented molecule).

    examples:

      COUNT("Ca_*",WORLD)
      COUNT("AChR'",WORLD)

 In: parse_state: parser state
     what: string representing the target of this count
     orient: orientation specified for the molecule
     where: symbol representing the count location (or NULL for WORLD)
     hit_spec: what are we counting?
     count_flags: is this a count or a trigger?
 Out: 0 on success, 1 on failure
**************************************************************************/
struct output_expression *mdl_count_syntax_3(struct mdlparse_vars *parse_state,
                                             char *what,
                                             struct sym_entry *where,
                                             int hit_spec, int count_flags) {
  struct output_expression *oe;
  char *what_to_count;
  byte report_flags = 0;
  if (where != NULL && parse_state->vol->periodic_box_obj && !(parse_state->vol->periodic_traditional)) {
    mdlerror(parse_state,
             "If PERIODIC_TRADITIONAL is FALSE, then you must specify virtual counting box.\n"
             "(e.g. COUNT,vm,Scene.box,[1,0,0]).");
  }
  if ((what_to_count = mdl_strip_quotes(what)) == NULL)
    return NULL;

  if (count_flags & TRIGGER_PRESENT)
    report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED)
    report_flags |= REPORT_ENCLOSED;

  /* Oriented molecule specified inside a string */
  if (mdl_string_has_orientation(what_to_count)) {
    struct output_request *orq;
    struct sym_entry *sp;
    short orientation;

    if (where == NULL) {
      mdlerror(parse_state, "Counting of an oriented molecule in the WORLD is "
                            "not implemented.\nAn oriented molecule may only "
                            "be counted in a regions.");
      free(what_to_count);
      return NULL;
    }

    orientation = mdl_get_orientation_from_string(what_to_count);
    if ((sp = mdl_existing_molecule(parse_state, what_to_count)) == NULL)
      return NULL;

    if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING)
      report_flags |= REPORT_CONTENTS;
    else
      report_flags |= (hit_spec & REPORT_TYPE_MASK);

    if ((orq = mcell_new_output_request(parse_state->vol, sp, orientation,
                                        where, NULL, report_flags)) == NULL)
      return NULL;
    orq->next = parse_state->vol->output_request_head;
    parse_state->vol->output_request_head = orq;
    oe = orq->requester;
  }

  /* Wildcard specified inside a string */
  else {
    struct sym_table_list *stl =
        mdl_find_rxpns_and_mols_by_wildcard(parse_state, what_to_count);

    if (stl == NULL) {
      mdlerror(parse_state,
               "Wildcard matching found no matches for count output.");
      free(what_to_count);
      return NULL;
    }

    if (where == NULL) {
      report_flags |= REPORT_WORLD;
      if (hit_spec != REPORT_NOTHING) {
        mdlerror(parse_state,
                 "Invalid combination of WORLD with other counting options");
        free(what_to_count);
        return NULL;
      } else if (count_flags & TRIGGER_PRESENT) {
        mdlerror(parse_state,
                 "Invalid combination of WORLD with TRIGGER option");
        free(what_to_count);
        return NULL;
      }
    }

    if ((oe = mdl_new_output_requests_from_list(
        parse_state, stl, where, report_flags, hit_spec, NULL)) == NULL) {
      free(what_to_count);
      return NULL;
    }

    free(what_to_count);
    /* free allocated memory */
    mem_put_list(parse_state->sym_list_mem, stl);
  }

  return oe;
}

/**************************************************************************
 mdl_single_count_expr:
    Prepare a single count expression for inclusion in an output set.

 In: parse_state: parser state
     list: list to receive output columns
     expr: the expression whose columns to add
     custom_header: custom header for this column
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_single_count_expr(struct mdlparse_vars *parse_state,
                          struct output_column_list *list,
                          struct output_expression *expr, char *custom_header) {
  if (expr->oper == ',' && custom_header != NULL) {
    mdlerror(parse_state,
             "Cannot use custom column headers with wildcard expansion");
    return 1;
  }

  return mcell_prepare_single_count_expr(list, expr, custom_header);
}

/*************************************************************************
 * VIZ output
 *************************************************************************/

/**************************************************************************
 mdl_new_viz_output_block:
    Build a new VIZ output block, containing parameters for an output set for
    visualization.

 In: parse_state: parser state
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_new_viz_output_block(struct mdlparse_vars *parse_state) {
  struct viz_output_block *vizblk = CHECKED_MALLOC_STRUCT(
      struct viz_output_block, "visualization data output parameters");
  if (vizblk == NULL)
    return 1;

  mcell_new_viz_output_block(vizblk);

  vizblk->next = parse_state->vol->viz_blocks;
  parse_state->vol->viz_blocks = vizblk;
  return 0;
}

/**************************************************************************
 mdl_set_viz_mode:
    Set the mode for a new VIZ output block.

 In: vizblk: the viz block to check
     mode: the mode to set
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_viz_mode(struct viz_output_block *vizblk, int mode) {

  vizblk->viz_mode = mode;
  return 0;
}

/**************************************************************************
 mdl_set_viz_filename_prefix:
    Set the filename prefix for a new VIZ output block.

 In: parse_state: parser state
     vizblk: the viz block to check
     filename: the filename
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_viz_filename_prefix(struct mdlparse_vars *parse_state,
                                struct viz_output_block *vizblk,
                                char *filename) {
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->file_prefix_name != NULL) {
    mdlerror_fmt(parse_state,
                 "FILENAME may only appear once per output block.");
    free(filename);
    return 1;
  }

  vizblk->file_prefix_name = filename;

  return 0;
}

/**************************************************************************
 mdl_viz_state:
    Sets a flag on all of the listed objects, requesting that they be
    visualized.

 In: parse_state: parser state
     target: destination for the viz state
     value: the raw (floating point!) value
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_viz_state(struct mdlparse_vars *parse_state, int *target,
                  double value) {
  /* Disallow out-of-range values. */
  if (value + 0.001 < (double)EXCLUDE_OBJ ||
      value - 0.001 > (double)INCLUDE_OBJ) {
    mdlerror(parse_state, "Visualization state is out of range "
                          "(must be between -(2^31 - 1) and (2^31 - 2).");
    return 1;
  }

  /* Round to integer. */
  int int_value;
  if (value < 0.0)
    int_value = (int)(value - 0.001);
  else
    int_value = (int)(value + 0.001);

  /* Disallow "special" values. */
  if (int_value == EXCLUDE_OBJ || int_value == INCLUDE_OBJ) {
    mdlerror(parse_state, "Visualization states of -(2^31) and (2^31) - 1 "
                          "are reserved for MCell's internal bookkeeping.");
    return 1;
  }

  *target = int_value;
  return 0;
}

/**************************************************************************
 mdl_set_viz_include_molecules:
    Sets a flag on all of the listed species, requesting that they be
    visualized.

 In: parse_state: parser state
     vizblk: the viz block to check
     list: the list of symbols
     viz_state: the desired viz state
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_viz_include_molecules(struct mdlparse_vars *parse_state,
                                  struct viz_output_block *vizblk,
                                  struct sym_table_list *list, int viz_state) {
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  /* Mark all specified molecules */
  struct sym_table_list *stl;
  for (stl = list; stl != NULL; stl = stl->next) {
    struct species *specp = (struct species *)stl->node->value;
    if (mcell_set_molecule_viz_state(vizblk, specp, viz_state))
      return 1;
  }

  /* free allocated memory  */
  mem_put_list(parse_state->sym_list_mem, list);
  return 0;
}

/**************************************************************************
 mdl_set_viz_include_all_molecules:
    Sets a flag on a viz block, requesting that all species be visualized.

 In: vizblk: the viz block to check
     viz_state: the desired viz state
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_viz_include_all_molecules(struct viz_output_block *vizblk,
                                      int viz_state) {
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (viz_state == INCLUDE_OBJ &&
      (vizblk->viz_output_flag & VIZ_ALL_MOLECULES)) {
    /* Do nothing - we will not override the old value if we have no specific
     * state value.
     */
  } else
    vizblk->default_mol_state = viz_state;
  vizblk->viz_output_flag |= VIZ_ALL_MOLECULES;
  return 0;
}

/**************************************************************************
 mdl_create_viz_mol_frames:
    Create one or more molecule frames for output in the visualization.

 In: parse_state: parser state
     time_type: either OUTPUT_BY_TIME_LIST or OUTPUT_BY_ITERATION_LIST
     type: the type (MOL_POS, MOL_ORIENT, etc.)
     viz_mode: visualization mode
     times: list of iterations/times at which to output
 Out: the frame_data_list object, if successful, or NULL if we ran out of memory
**************************************************************************/
static struct frame_data_list *
mdl_create_viz_mol_frames(struct mdlparse_vars *parse_state, int time_type,
                          int type, int viz_mode,
                          struct num_expr_list_head *times) {
  struct frame_data_list *frames = NULL;
  struct frame_data_list *new_frame;
  struct num_expr_list *times_sorted;
  if (times->shared) {
    times_sorted = mcell_copysort_numeric_list(times->value_head);
    if (times_sorted == NULL)
      return NULL;
  } else {
    mcell_sort_numeric_list(times->value_head);
    times_sorted = times->value_head;
  }

  if (type == MOL_POS || type == ALL_MOL_DATA) {
    if ((new_frame = mcell_create_viz_frame(time_type, type, times_sorted)) ==
        NULL)
      return NULL;
    new_frame->next = frames;
    frames = new_frame;
  } else if (viz_mode == NO_VIZ_MODE) {
    /* Create viz frames consistent with other visualization modes */
    if ((new_frame = mcell_create_viz_frame(time_type, type, times_sorted)) ==
        NULL)
      return NULL;
    new_frame->next = frames;
    frames = new_frame;
  } else {
    mdlerror_fmt(parse_state,
                 "This type of molecule output data (%d) is not valid "
                 "for the selected VIZ output mode (%d).",
                 type, viz_mode);
    return NULL;
  }

  return frames;
}

/**************************************************************************
 mdl_new_viz_mol_frames:
    Adds some new molecule output frames to a list.

 In: parse_state: parser state
     vizblk: the viz block to check
     frames: list to receive frames
     time_type: timing type (OUTPUT_BY_TIME_LIST or ...ITERATION_LIST)
     mol_item_type: MOLECULE_POSITIONS, etc.
     timelist: list of times in appropriate units (as per time_type)
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_new_viz_mol_frames(struct mdlparse_vars *parse_state,
                           struct viz_output_block *vizblk,
                           struct frame_data_list_head *frames, int time_type,
                           int mol_item_type,
                           struct num_expr_list_head *timelist) {
  frames->frame_head = frames->frame_tail = NULL;
  // if (vizblk->viz_mode == NO_VIZ_MODE)
  //   return 0;

  struct frame_data_list *fdlp;
  fdlp = mdl_create_viz_mol_frames(parse_state, time_type, mol_item_type,
                                   vizblk->viz_mode, timelist);
  if (!fdlp)
    return 1;

  frames->frame_head = fdlp;
  while (fdlp->next != NULL)
    fdlp = fdlp->next;
  frames->frame_tail = fdlp;
  return 0;
}

/**************************************************************************
 mdl_new_viz_all_times:
    Build a list of times for VIZ output, one timepoint per iteration in the
    simulation.

 In: parse_state: parser state
     list: location to receive timepoint list
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_new_viz_all_times(struct mdlparse_vars *parse_state,
                          struct num_expr_list_head *list) {
  long long step;
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  for (step = 0; step <= parse_state->vol->iterations; step++) {
    struct num_expr_list *nel;
    nel = CHECKED_MALLOC_STRUCT(struct num_expr_list, "VIZ_OUTPUT time point");
    if (nel == NULL)
      return 1;

    ++list->value_count;
    if (list->value_tail)
      list->value_tail = list->value_tail->next = nel;
    else
      list->value_head = list->value_tail = nel;
    list->value_tail->value = step * parse_state->vol->time_unit;
    list->value_tail->next = NULL;
  }
  return 0;
}

/**************************************************************************
 mdl_new_viz_all_iterations:
    Build a list of iterations for VIZ output, one for each iteration in the
    simulation.

 In: parse_state: parser state
     list: location to receive iteration list
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_new_viz_all_iterations(struct mdlparse_vars *parse_state,
                               struct num_expr_list_head *list) {
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  for (long long step = 0; step <= parse_state->vol->iterations; step++) {
    struct num_expr_list *nel;
    nel = CHECKED_MALLOC_STRUCT(struct num_expr_list, "VIZ_OUTPUT iteration");
    if (nel == NULL)
      return 1;

    ++list->value_count;
    if (list->value_tail)
      list->value_tail = list->value_tail->next = nel;
    else
      list->value_head = list->value_tail = nel;
    list->value_tail->value = step;
    list->value_tail->next = NULL;
  }
  return 0;
}

/*************************************************************************
 * Volume output
 *************************************************************************/

/*************************************************************************
 species_list_to_array:
    Convert a pointer list into an array.  The list count will be stored into
    the count pointer, if it is not NULL.

 In:  lh: list head
      count: location to receive list item count
 Out: array of pointers to species in list, or NULL if allocation fails, or if
      the list is empty.  If NULL is returned and the count field has a
      non-zero value, no error has occurred.
*************************************************************************/
static struct species **species_list_to_array(struct species_list *lh,
                                              int *count) {

  struct species **arr = NULL, **ptr = NULL;
  struct species_list_item *i;

  if (count != NULL)
    *count = lh->species_count;

  if (lh->species_count == 0)
    return NULL;

  arr = CHECKED_MALLOC_ARRAY(struct species *, lh->species_count,
                             "species array");
  if (arr) {
    for (i = (struct species_list_item *)lh->species_head, ptr = arr; i != NULL;
         i = i->next)
      *ptr++ = i->spec;
  }
  return arr;
}

/**************************************************************************
 mdl_new_volume_output_item:
    Create a new volume output request.

 In: parse_state: parser state
     filename_prefix: filename prefix for all files produced from this request
     molecules: list of molecules to output
     location: lower, left, front corner of output box
     voxel_size: dimensions of each voxel
     voxel_count: counts of voxels in x, y, and z directions
     ot: output times for this request
 Out: volume output item, or NULL if an error occurred
**************************************************************************/
struct volume_output_item *mdl_new_volume_output_item(
    struct mdlparse_vars *parse_state, char *filename_prefix,
    struct species_list *molecules, struct vector3 *location,
    struct vector3 *voxel_size, struct vector3 *voxel_count,
    struct output_times *ot) {
  struct volume_output_item *vo =
      CHECKED_MALLOC_STRUCT(struct volume_output_item, "volume output request");
  if (vo == NULL) {
    free(filename_prefix);
    mem_put_list(parse_state->species_list_mem, molecules->species_head);
    free(location);
    free(voxel_size);
    free(voxel_count);
    if (ot->times != NULL)
      free(ot->times);
    mem_put(parse_state->output_times_mem, ot);
    return NULL;
  }
  memset(vo, 0, sizeof(struct volume_output_item));

  vo->filename_prefix = filename_prefix;
  vo->molecules = species_list_to_array(molecules, &vo->num_molecules);
  if (vo->num_molecules != 0 && vo->molecules == NULL) {
    free(filename_prefix);
    mem_put_list(parse_state->species_list_mem, molecules->species_head);
    free(location);
    free(voxel_size);
    free(voxel_count);
    if (ot->times != NULL)
      free(ot->times);
    mem_put(parse_state->output_times_mem, ot);
    free(vo);
    return NULL;
  }

  qsort(vo->molecules, vo->num_molecules, sizeof(void *), &void_ptr_compare);
  mem_put_list(parse_state->species_list_mem, molecules->species_head);

  memcpy(&vo->location, location, sizeof(struct vector3));
  free(location);
  vo->location.x *= parse_state->vol->r_length_unit;
  vo->location.y *= parse_state->vol->r_length_unit;
  vo->location.z *= parse_state->vol->r_length_unit;

  memcpy(&vo->voxel_size, voxel_size, sizeof(struct vector3));
  free(voxel_size);
  vo->voxel_size.x *= parse_state->vol->r_length_unit;
  vo->voxel_size.y *= parse_state->vol->r_length_unit;
  vo->voxel_size.z *= parse_state->vol->r_length_unit;

  vo->nvoxels_x = (int)(voxel_count->x + 0.5);
  vo->nvoxels_y = (int)(voxel_count->y + 0.5);
  vo->nvoxels_z = (int)(voxel_count->z + 0.5);
  free(voxel_count);

  vo->timer_type = ot->timer_type;
  vo->step_time = ot->step_time;
  vo->num_times = ot->num_times;
  vo->times = ot->times;
  mem_put(parse_state->output_times_mem, ot);

  if (vo->timer_type == OUTPUT_BY_ITERATION_LIST ||
      vo->timer_type == OUTPUT_BY_TIME_LIST)
    vo->next_time = vo->times;
  else
    vo->next_time = NULL;

  switch (parse_state->vol->notify->volume_output_report) {
  case NOTIFY_NONE:
    break;

  case NOTIFY_BRIEF:
    mcell_log("Added volume output item '%s', counting in the region "
              "[%.15g,%.15g,%.15g]-[%.15g,%.15g,%.15g]",
              vo->filename_prefix,
              vo->location.x * parse_state->vol->length_unit,
              vo->location.y * parse_state->vol->length_unit,
              vo->location.z * parse_state->vol->length_unit,
              (vo->location.x + vo->voxel_size.x * vo->nvoxels_x) *
                  parse_state->vol->length_unit,
              (vo->location.y + vo->voxel_size.x * vo->nvoxels_x) *
                  parse_state->vol->length_unit,
              (vo->location.z + vo->voxel_size.x * vo->nvoxels_x) *
                  parse_state->vol->length_unit);
    break;

  case NOTIFY_FULL:
    mcell_log("Added volume output item '%s', counting the following molecules "
              "in the region [%.15g,%.15g,%.15g]-[%.15g,%.15g,%.15g]:\n",
              vo->filename_prefix,
              vo->location.x * parse_state->vol->length_unit,
              vo->location.y * parse_state->vol->length_unit,
              vo->location.z * parse_state->vol->length_unit,
              (vo->location.x + vo->voxel_size.x * vo->nvoxels_x) *
                  parse_state->vol->length_unit,
              (vo->location.y + vo->voxel_size.x * vo->nvoxels_x) *
                  parse_state->vol->length_unit,
              (vo->location.z + vo->voxel_size.x * vo->nvoxels_x) *
                  parse_state->vol->length_unit);
    for (int i = 0; i < vo->num_molecules; ++i)
      mcell_log("  %s", vo->molecules[i]->sym->name);
    break;

  default:
    UNHANDLED_CASE(parse_state->vol->notify->volume_output_report);
  }

  return vo;
}

/**************************************************************************
 mdl_new_output_times_default:
    Create new default output timing for volume output.

 In: parse_state: parser state
 Out: output times structure, or NULL if allocation fails
**************************************************************************/
struct output_times *
mdl_new_output_times_default(struct mdlparse_vars *parse_state) {
  struct output_times *ot = CHECKED_MEM_GET(parse_state->output_times_mem,
                                            "output times for volume output");
  if (ot == NULL)
    return NULL;
  memset(ot, 0, sizeof(struct output_times));
  ot->timer_type = OUTPUT_BY_STEP;
  ot->step_time = parse_state->vol->time_unit;
  return ot;
}

/**************************************************************************
 mdl_new_output_times_step:
    Create new "step" output timing for volume output.

 In: parse_state: parser state
     step: time step for volume output
 Out: output times structure, or NULL if allocation fails
 XXX: This is really similar to set_reaction_output_timer_step in
 mcell_react_out.c. Consolidate these.
**************************************************************************/
struct output_times *
mdl_new_output_times_step(struct mdlparse_vars *parse_state, double step) {
  long long output_freq;
  struct output_times *ot = CHECKED_MEM_GET(parse_state->output_times_mem,
                                            "output times for volume output");
  if (ot == NULL)
    return NULL;
  memset(ot, 0, sizeof(struct output_times));
  ot->timer_type = OUTPUT_BY_STEP;
  ot->step_time = step;

  /* Clip step_time to a reasonable range */
  output_freq = ot->step_time / parse_state->vol->time_unit;
  if (output_freq > parse_state->vol->iterations && output_freq > 1) {
    output_freq =
        (parse_state->vol->iterations > 1) ? parse_state->vol->iterations : 1;
    ot->step_time = output_freq * parse_state->vol->time_unit;
    if (parse_state->vol->notify->invalid_output_step_time != WARN_COPE)
      mdl_warning(parse_state, "Output step time too long\n\tSetting output "
                               "step time to %g microseconds\n",
                  ot->step_time * 1.0e6);
  } else if (output_freq < 1) {
    ot->step_time = parse_state->vol->time_unit;
    if (parse_state->vol->notify->invalid_output_step_time != WARN_COPE)
      mdl_warning(parse_state, "Output step time too short\n\tSetting output "
                               "step time to %g microseconds\n",
                  ot->step_time * 1.0e6);
  }
  return ot;
}

/**************************************************************************
 mdl_new_output_times_iterations:
    Create new "iteration list" output timing for volume output.

 In: parse_state: parser state
     iters: iterations on which to give volume output
 Out: output times structure, or NULL if allocation fails
**************************************************************************/
struct output_times *
mdl_new_output_times_iterations(struct mdlparse_vars *parse_state,
                                struct num_expr_list_head *iters) {
  struct output_times *ot = CHECKED_MEM_GET(parse_state->output_times_mem,
                                            "output times for volume output");
  if (ot == NULL) {
    if (!iters->shared)
      mcell_free_numeric_list(iters->value_head);
    return NULL;
  }
  memset(ot, 0, sizeof(struct output_times));

  ot->timer_type = OUTPUT_BY_ITERATION_LIST;
  ot->times = num_expr_list_to_array(iters, &ot->num_times);
  if (ot->times != NULL)
    qsort(ot->times, ot->num_times, sizeof(double), &double_cmp);
  if (!iters->shared)
    mcell_free_numeric_list(iters->value_head);

  if (ot->num_times != 0 && ot->times == NULL) {
    mem_put(parse_state->output_times_mem, ot);
    return NULL;
  }

  return ot;
}

/**************************************************************************
 mdl_new_output_times_time:
    Create new "time list" output timing for volume output.

 In: parse_state: parser state
     times: simulation times at which to give volume output
 Out: output times structure, or NULL if allocation fails
**************************************************************************/
struct output_times *
mdl_new_output_times_time(struct mdlparse_vars *parse_state,
                          struct num_expr_list_head *times) {
  struct output_times *ot = CHECKED_MEM_GET(parse_state->output_times_mem,
                                            "output times for volume output");
  if (ot == NULL) {
    if (!times->shared)
      mcell_free_numeric_list(times->value_head);
    return NULL;
  }
  memset(ot, 0, sizeof(struct output_times));

  ot->timer_type = OUTPUT_BY_TIME_LIST;
  ot->times = num_expr_list_to_array(times, &ot->num_times);
  if (ot->times != NULL)
    qsort(ot->times, ot->num_times, sizeof(double), &double_cmp);
  if (!times->shared)
    mcell_free_numeric_list(times->value_head);

  if (ot->num_times != 0 && ot->times == NULL) {
    mem_put(parse_state->output_times_mem, ot);
    return NULL;
  }

  return ot;
}

/****************************************************************
 * Release patterns
 ***************************************************************/

/**************************************************************************
 mdl_new_release_pattern:
    Create a new release pattern.  There must not yet be a release pattern with
    the given name.

 In: parse_state: parser state
     name: name for the new release pattern
 Out: symbol of the release pattern, or NULL if an error occurred
**************************************************************************/
struct sym_entry *mdl_new_release_pattern(struct mdlparse_vars *parse_state,
                                          char *name) {
  struct sym_entry *st;
  if (retrieve_sym(name, parse_state->vol->rpat_sym_table) != NULL) {
    mdlerror_fmt(parse_state, "Release pattern already defined: %s", name);
    free(name);
    return NULL;
  } else if ((st = store_sym(name, RPAT, parse_state->vol->rpat_sym_table,
                             NULL)) == NULL) {
    mdlerror_fmt(parse_state,
                 "Out of memory while creating release pattern: %s", name);
    free(name);
    return NULL;
  }

  free(name);
  return st;
}

/**************************************************************************
 mdl_set_release_pattern:
    Fill in the details of a release pattern.

 In: parse_state: parser state
     rpat_sym: symbol for the release pattern
     rpat_data: data to fill in for release pattern
 Out: 0 on succes, 1 on failure.
**************************************************************************/
int mdl_set_release_pattern(struct mdlparse_vars *parse_state,
                            struct sym_entry *rpat_sym,
                            struct release_pattern *rpat_data) {
  struct release_pattern *rpatp = (struct release_pattern *)rpat_sym->value;
  if (rpat_data->release_interval <= 0) {
    mdlerror(parse_state, "Release interval must be set to a positive number.");
    return 1;
  }
  if (rpat_data->train_interval <= 0) {
    mdlerror(parse_state, "Train interval must be set to a positive number.");
    return 1;
  }
  if (rpat_data->train_duration > rpat_data->train_interval) {
    mdlerror(parse_state,
             "Train duration must not be longer than the train interval.");
    return 1;
  }
  if (rpat_data->train_duration <= 0) {
    mdlerror(parse_state, "Train duration must be set to a positive number.");
    return 1;
  }

  /* Copy in release pattern */
  rpatp->delay = rpat_data->delay;
  rpatp->release_interval = rpat_data->release_interval;
  rpatp->train_interval = rpat_data->train_interval;
  rpatp->train_duration = rpat_data->train_duration;
  rpatp->number_of_trains = rpat_data->number_of_trains;

  no_printf("Release pattern %s defined:\n", rpat_sym->name);
  no_printf("\tdelay = %f\n", rpatp->delay);
  no_printf("\trelease_interval = %f\n", rpatp->release_interval);
  no_printf("\ttrain_interval = %f\n", rpatp->train_interval);
  no_printf("\ttrain_duration = %f\n", rpatp->train_duration);
  no_printf("\tnumber_of_trains = %d\n", rpatp->number_of_trains);
  return 0;
}

/****************************************************************
 * Molecules
 ***************************************************************/

/**************************************************************************
 mdl_new_mol_species:
    Create a new species. There must not yet be a molecule or named reaction
    pathway with the supplied name.

 In: parse_state: parser state
     name: name for the new species
 Out: symbol for the species, or NULL if an error occurred
**************************************************************************/
struct sym_entry *mdl_new_mol_species(struct mdlparse_vars *parse_state,
                                      char *name) {
  struct sym_entry *sym = NULL;
  if (retrieve_sym(name, parse_state->vol->mol_sym_table) != NULL) {
    mdlerror_fmt(parse_state, "Molecule already defined: %s", name);
    free(name);
    return NULL;
  } else if (retrieve_sym(name, parse_state->vol->rxpn_sym_table) != NULL) {
    mdlerror_fmt(parse_state,
                 "Molecule already defined as a named reaction pathway: %s",
                 name);
    free(name);
    return NULL;
  } else if ((sym = store_sym(name, MOL, parse_state->vol->mol_sym_table,
                              NULL)) == NULL) {
    mdlerror_fmt(parse_state, "Out of memory while creating molecule: %s",
                 name);
    free(name);
    return NULL;
  }

  free(name);
  return sym;
}

/**************************************************************************
 mdl_create_species:
    Assemble a molecule species from its component pieces.

 In: parse_state:      parser state
     name:             name of the molecule
     D:                diffusion constant
     is_2d:            1 if the species is a 2D molecule, 0 if 3D
     custom_time_step: time_step for the molec (< 0.0 for a custom space step,
                       >0.0 for custom timestep, 0.0 for default timestep)
     target_only:      1 if the molecule cannot initiate reactions
     max_step_length:
 Out: Nothing. The molecule is created.
**************************************************************************/
struct mcell_species_spec *mdl_create_species(struct mdlparse_vars *parse_state,
                                              char *name, double D, int is_2d,
                                              double custom_time_step,
                                              int target_only,
                                              double max_step_length) {
  // Can't define molecule before we have a time step.
  // Move this to mcell_create_species?
  double global_time_unit = parse_state->vol->time_unit;
  if (!distinguishable(global_time_unit, 0, EPS_C)) {
    mdlerror_fmt(parse_state,
                 "TIME_STEP not yet specified.  Cannot define molecule: %s",
                 name);
  }

  struct mcell_species_spec *species =
      CHECKED_MALLOC_STRUCT(struct mcell_species_spec, "struct mcell_species");
  species->name = name;
  species->D = D;
  species->is_2d = is_2d;
  species->custom_time_step = custom_time_step;
  species->target_only = target_only;
  species->max_step_length = max_step_length;
  int error_code = mcell_create_species(parse_state->vol, species, NULL);

  switch (error_code) {
  case 2:
    mdlerror_fmt(parse_state, "Molecule already defined: %s", name);
    break;
  case 3:
    mdlerror_fmt(parse_state,
                 "Molecule already defined as a named reaction pathway: %s",
                 name);
    break;
  case 4:
    mdlerror_fmt(parse_state, "Out of memory while creating molecule: %s",
                 name);
    break;
  case 5:
    mdlerror(parse_state,
             "Out of memory while creating r_step data for molecule");
    break;
  case 6:
    mdlerror(parse_state, "Cannot store r_step_surface data.");
    break;
  case 7:
    mdlerror(parse_state,
             "Out of memory while creating d_step data for molecule");
    break;
  case 8:
    mdlerror(parse_state, "Internal error: bad number of default "
                          "RADIAL_DIRECTIONS (max 131072).");
    break;
  }

  return species;
}

/****************************************************************
 * Reactions, surface classes
 ***************************************************************/

/**************************************************************************
 mdl_valid_rate:
    Check whether the reaction rate is valid.

 In: parse_state: parser state
     rate: the unidirectional (either forward or reverse) reaction rate
 Out: 0 if valid, 1 if not
**************************************************************************/
int mdl_valid_rate(struct mdlparse_vars *parse_state,
                   struct reaction_rate *rate) {
  if (rate->rate_type == RATE_UNSET) {
    mdlerror_fmt(parse_state,
                 "File %s, Line %d: Internal error: Rate is not set", __FILE__,
                 __LINE__);
    return 1;
  } else if (rate->rate_type == RATE_CONSTANT) {
    if (rate->v.rate_constant < 0.0) {
      if (parse_state->vol->notify->neg_reaction == WARN_ERROR) {
        mdlerror(parse_state,
                 "reaction rate constants should be zero or positive.");
        return 1;
      } else if (parse_state->vol->notify->neg_reaction == WARN_WARN) {
        mcell_warn("negative reaction rate constant %f; setting to zero and "
                   "continuing.", rate->v.rate_constant);
        rate->v.rate_constant = 0.0;
      }
    }
  } else if (rate->rate_type == RATE_FILE) {
    if (rate->v.rate_file == NULL) {
      mdlerror_fmt(parse_state,
                   "File %s, Line %d: Internal error: Rate filename is not set",
                   __FILE__, __LINE__);
      return 1;
    }
  }

  return 0;
}

/**************************************************************************
 mdl_new_reaction_player:
    Create a new reaction player from a species with optional orientation.

 In: parse_state: parser state
     spec: species with optional orientation
 Out: reaction player, or NULL if allocation failed
**************************************************************************/
static struct mcell_species *
mdl_new_reaction_player(struct mdlparse_vars *parse_state,
                        struct mcell_species *spec) {
  struct mcell_species *new_spec;
  if ((new_spec = (struct mcell_species *)CHECKED_MEM_GET(
           parse_state->mol_data_list_mem, "molecule type")) == NULL)
    return NULL;

  *new_spec = *spec;
  new_spec->next = NULL;
  return new_spec;
}

/**************************************************************************
 mdl_reaction_player_singleton:
    Set the reactant/product list to contain a single item.

 In: parse_state: parser state
     list: list to receive player
     spec: species with optional orientation
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_reaction_player_singleton(struct mdlparse_vars *parse_state,
                                  struct mcell_species_list *list,
                                  struct mcell_species *spec) {
  struct mcell_species *player = mdl_new_reaction_player(parse_state, spec);
  if (player == NULL)
    return 1;
  list->mol_type_head = list->mol_type_tail = player;
  return 0;
}

/**************************************************************************
 mdl_add_reaction_player:
    Add a single item to a reactant/product player list.

 In: parse_state: parser state
     list: list to receive player
     spec: species with optional orientation
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_add_reaction_player(struct mdlparse_vars *parse_state,
                            struct mcell_species_list *list,
                            struct mcell_species *spec) {
  struct mcell_species *player = mdl_new_reaction_player(parse_state, spec);
  if (player == NULL)
    return 1;
  list->mol_type_tail->next = player;
  list->mol_type_tail = player;
  return 0;
}

/**************************************************************************
 mdl_reaction_rate_from_var:
    Set a reaction rate from a variable.

 In: parse_state: parser state
     rate: reaction rate struct
     symp: pointer to symbol for reaction rate
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_reaction_rate_from_var(struct mdlparse_vars *parse_state,
                               struct reaction_rate *rate,
                               struct sym_entry *symp) {
  switch (symp->sym_type) {
  case DBL:
    rate->rate_type = RATE_CONSTANT;
    rate->v.rate_constant = *(double *)symp->value;
    break;

  case STR:
    rate->rate_type = RATE_FILE;
    if ((rate->v.rate_file = mdl_strdup((char *)symp->value)) == NULL)
      return 1;
    break;

  default:
    mdlerror(parse_state, "Invalid variable used for reaction rate: must be a "
                          "number or a filename");
    return 1;
  }
  return 0;
}

/**************************************************************************
 mdl_assemble_reaction:
    Assemble a standard reaction from its component parts.

 In: parse_state: parser state
     reactants: list of reactants
     surface_class: optional surface class
     react_arrow: reaction arrow (uni, bi, catalytic, etc)
     products: list of products
     rates: forward/reverse reaction rates
     pathname: name reaction pathway name, or NULL
 Out: the reaction, or NULL if an error occurred
**************************************************************************/
struct mdlparse_vars *mdl_assemble_reaction(struct mdlparse_vars *parse_state,
                                            struct mcell_species *reactants,
                                            struct mcell_species *surface_class,
                                            struct reaction_arrow *react_arrow,
                                            struct mcell_species *products,
                                            struct reaction_rates *rate,
                                            struct sym_entry *pathname) {
  char *forward_rate_filename = NULL;
  char *backward_rate_filename = NULL;
  if (rate->forward_rate.rate_type == RATE_FILE) {
    forward_rate_filename =
        mcell_find_include_file(rate->forward_rate.v.rate_file,
                                parse_state->vol->curr_file);
  }
  if (rate->backward_rate.rate_type == RATE_FILE) {
    backward_rate_filename =
        mcell_find_include_file(rate->backward_rate.v.rate_file,
                                parse_state->vol->curr_file);
  }

  struct volume *state = parse_state->vol;
  if (mcell_add_reaction(state->notify,
                         &state->r_step_release,
                         state->rxn_sym_table,
                         state->radial_subdivisions,
                         state->vacancy_search_dist2,
                         reactants, react_arrow, surface_class, products,
                         pathname, rate, forward_rate_filename,
                         backward_rate_filename)) {
    return NULL;
  }

  return parse_state;
}

/**************************************************************************
 mdl_assemble_surface_reaction:
    Assemble a surface reaction from its component parts.

 In: parse_state: parser state
     reaction_type: RFLCT, TRANSP, or SINK
     surface_class: surface class
     reactant_sym: symbol for reactant molecule
     orient: orientation
 Out: the reaction, or NULL if an error occurred
**************************************************************************/
struct mdlparse_vars *
mdl_assemble_surface_reaction(struct mdlparse_vars *parse_state,
                              int reaction_type, struct species *surface_class,
                              struct sym_entry *reactant_sym, short orient) {
  if (mcell_add_surface_reaction(parse_state->vol->rxn_sym_table, reaction_type,
                                 surface_class, reactant_sym, orient)) {
    return NULL;
  }

  return parse_state;
}

/**************************************************************************
 mdl_assemble_concentration_clamp_reaction:
    Assemble a concentration clamp reaction from its component parts.

 In: parse_state: parser state
     surface_class: surface class
     mol_sym: symbol for molecule being clamped
     orient: orientation
     conc: desired concentration
 Out: the reaction, or NULL if an error occurred
**************************************************************************/
struct mdlparse_vars *mdl_assemble_concentration_clamp_reaction(
    struct mdlparse_vars *parse_state, struct species *surface_class,
    struct sym_entry *mol_sym, short orient, double conc) {
  if (mcell_add_concentration_clamp(parse_state->vol->rxn_sym_table,
                                    surface_class, mol_sym, orient, conc)) {
    return NULL;
  }

  return parse_state;
}

/**************************************************************************
 mdl_start_surface_class:
    Start a surface class declaration.

    SIDE EFFECT: sets current_surface_class in the parser state

 In: parse_state: parser state
     symp: symbol for the surface class
 Out: 0 on success, 1 on failure
**************************************************************************/
void mdl_start_surface_class(struct mdlparse_vars *parse_state,
                             struct sym_entry *symp) {
  struct species *specp = (struct species *)symp->value;
  specp->flags = IS_SURFACE;
  specp->refl_mols = NULL;
  specp->transp_mols = NULL;
  specp->absorb_mols = NULL;
  specp->clamp_conc_mols = NULL;
  parse_state->current_surface_class = specp;
}

/**************************************************************************
 mdl_finish_surface_class:
    Finish a surface class declaration.  Undoes side effects from
    mdl_start_surface_class.

 In: parse_state: parser state
 Out: 0 on success, 1 on failure
**************************************************************************/
void mdl_finish_surface_class(struct mdlparse_vars *parse_state) {
  parse_state->current_surface_class = NULL;
}

/**************************************************************************
 mdl_new_surf_mol_data:
    Create a new surface molecule data for surface molecule initialization.

 In: parse_state: parser state
     sm_info:
     quant: the amount of surface molecules to release
 Out: 0 on success, 1 on failure
**************************************************************************/
struct sm_dat *mdl_new_surf_mol_data(struct mdlparse_vars *parse_state,
                                     struct mcell_species *sm_info,
                                     double quant) {
  struct sm_dat *smdp;
  struct species *specp = (struct species *)sm_info->mol_type->value;
  if (!(specp->flags & ON_GRID)) {
    mdlerror_fmt(parse_state,
                 "Cannot initialize surface with non-surface molecule '%s'",
                 specp->sym->name);
    return NULL;
  } else if (mdl_check_valid_molecule_release(parse_state, sm_info)) {
    return NULL;
  }

  if ((smdp = CHECKED_MALLOC_STRUCT(struct sm_dat, "surface molecule data")) ==
      NULL)
    return NULL;

  smdp->next = NULL;
  smdp->sm = specp;
  smdp->quantity_type = 0;
  smdp->quantity = quant;
  smdp->orientation = sm_info->orient;
  return smdp;
}

/**********************************************************************
 ***  helper functions for release sites creation
 **********************************************************************/

/**************************************************************************
 add_release_single_molecule_to_list:
    Adds a release molecule descriptor to a list.

 In: list: the list
     mol:  the descriptor
 Out: none.  list is updated
**************************************************************************/
void
add_release_single_molecule_to_list(struct release_single_molecule_list *list,
                                    struct release_single_molecule *mol) {
  list->rsm_tail = list->rsm_tail->next = mol;
  ++list->rsm_count;
}

/**************************************************************************
 release_single_molecule_singleton:
    Populates a list with a single LIST release molecule descriptor.

 In: list: the list
     mol:  the descriptor
 Out: none.  list is updated
**************************************************************************/
void
release_single_molecule_singleton(struct release_single_molecule_list *list,
                                  struct release_single_molecule *mol) {
  list->rsm_tail = list->rsm_head = mol;
  list->rsm_count = 1;
}

/**************************************************************************
 set_release_site_density:
    Set a release quantity from this release site based on a fixed
    density within the release-site's area.  (Hopefully we're talking about a
    surface release here.)

 In: rel_site_obj_ptr: the release site
     dens: density for release
 Out: 0 on success, 1 on failure.  release site object is updated
**************************************************************************/
int set_release_site_density(struct release_site_obj *rel_site_obj_ptr,
                             double dens) {

  rel_site_obj_ptr->release_number_method = DENSITYNUM;
  rel_site_obj_ptr->concentration = dens;
  return 0;
}

/**************************************************************************
 set_release_site_volume_dependent_number:
    Set a release quantity from this release site based on a fixed
    concentration in a sphere of a gaussian-distributed diameter with a
    particular mean and std. deviation.

 In: rel_site_obj_ptr: the release site
     mean: mean value of distribution of diameters
     stdev: std. dev. of distribution of diameters
     conc: concentration for release
 Out: none.  release site object is updated
**************************************************************************/
void set_release_site_volume_dependent_number(
    struct release_site_obj *rel_site_obj_ptr, double mean, double stdev,
    double conc) {
  rel_site_obj_ptr->release_number_method = VOLNUM;
  rel_site_obj_ptr->mean_diameter = mean;
  rel_site_obj_ptr->standard_deviation = stdev;
  rel_site_obj_ptr->concentration = conc;
}

/*************************************************************************
 transform_translate:
    Apply a translation to the given transformation matrix.

 In:  state: system state
      mat: transformation matrix
      xlat: translation vector
 Out: translation is right-multiplied into the transformation matrix
*************************************************************************/
void transform_translate(MCELL_STATE *state, double (*mat)[4],
                         struct vector3 *xlat) {
  double tm[4][4];
  struct vector3 scaled_xlat = *xlat;
  scaled_xlat.x *= state->r_length_unit;
  scaled_xlat.y *= state->r_length_unit;
  scaled_xlat.z *= state->r_length_unit;
  init_matrix(tm);
  translate_matrix(tm, tm, &scaled_xlat);
  mult_matrix(mat, tm, mat, 4, 4, 4);
  free(xlat);
}

/*************************************************************************
 transform_scale:
    Apply a scale to the given transformation matrix.

 In:  mat: transformation matrix
      scale: scale vector
 Out: scale is right-multiplied into the transformation matrix
*************************************************************************/
void transform_scale(double (*mat)[4], struct vector3 *scale) {
  double tm[4][4];
  init_matrix(tm);
  scale_matrix(tm, tm, scale);
  mult_matrix(mat, tm, mat, 4, 4, 4);
  free(scale);
}

/*************************************************************************
 transform_rotate:
    Apply a rotation to the given transformation matrix.

 In:  mat: transformation matrix
      axis: axis of rotation
      angle: angle of rotation (degrees!)
 Out: 0 on success, 1 on failure; rotation is right-multiplied into the
      transformation matrix
*************************************************************************/
int transform_rotate(double (*mat)[4], struct vector3 *axis, double angle) {
  double tm[4][4];
  if (!distinguishable(vect_length(axis), 0.0, EPS_C)) {
    return 1;
  }
  init_matrix(tm);
  rotate_matrix(tm, tm, axis, angle);
  mult_matrix(mat, tm, mat, 4, 4, 4);
  return 0;
}

/*******************************************************************************
 *
 * check_release_regions makes sure that all regions used in release objects
 * correspond to properly instantiated objects.
 *
 *******************************************************************************/
void check_regions(struct object *rootInstance, struct object *child) {

  while (child != NULL) {
    for (struct object *fc = child->first_child; fc != NULL; fc = fc->next) {
      if (fc->object_type == REL_SITE_OBJ) {
        struct release_site_obj *rel =
            (struct release_site_obj *)(fc->contents);
        if (rel->region_data != NULL) {
          struct release_evaluator *eval = rel->region_data->expression;
          if (check_release_regions(eval, fc, rootInstance)) {
            mcell_error("Release object %s contains at least one uninstantiated"
                        "region.",
                        fc->sym->name);
          }
        }
      }
    }
    child = child->next;
  }
}

/**************************************************************************
 finish_polygon_list:
    Finalize the polygon list, cleaning up any state updates that were made
    when we started creating the polygon.

 In: obj_ptr: contains information about the object (name, etc)
     obj_creation: information about object being created
 Out: 1 on failure, 0 on success
 NOTE: This function call might be too low-level for what we want from the API,
       but it is needed to create polygon objects for now.
**************************************************************************/
int finish_polygon_list(struct object *obj_ptr,
                        struct object_creation *obj_creation) {
  pop_object_name(obj_creation);
  remove_gaps_from_regions(obj_ptr);
  // no_printf(" n_verts = %d\n", mpvp->current_polygon->n_verts);
  // no_printf(" n_walls = %d\n", mpvp->current_polygon->n_walls);
  if (check_degenerate_polygon_list(obj_ptr)) {
    return 1;
  }
  return 0;
}

/*************************************************************************
 start_object:
    Create a new object, adding it to the global symbol table. The object must
    not be defined yet. The qualified name of the object will be built by
    adding to the object_name_list, and the object is made the "current_object"
    in the mdl parser state. Because of these side effects, it is vital to call
    finish_object at the end of the scope of the object created here.

 In:  state: the simulation state
      obj_creation: information about object being created
      name: unqualified object name
 Out: the newly created object
 NOTE: This is very similar to mdl_start_object, but there is no parse state.
*************************************************************************/
struct object *start_object(MCELL_STATE *state,
                            struct object_creation *obj_creation,
                            char *name,
                            int *error_code) {
  // Create new fully qualified name.
  char *new_name;
  if ((new_name = push_object_name(obj_creation, name)) == NULL) {
    free(name);
    return NULL;
  }

  struct dyngeom_parse_vars *dg_parse = state->dg_parse;
  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = make_new_object(
      dg_parse,
      state->obj_sym_table,
      new_name,
      error_code);
  if (*error_code == 1) {
    return NULL;
  }

  obj_ptr->last_name = name;
  no_printf("Creating new object: %s\n", new_name);

  // Set parent object, make this object "current".
  obj_ptr->parent = obj_creation->current_object;

  return obj_ptr;
}
