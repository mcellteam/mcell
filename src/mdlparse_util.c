/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by                                                      *
 * The Salk Institute for Biological Studies and                                   *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University                    *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

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
#include "vector.h"
#include "strfunc.h"
#include "sym_table.h"
#include "util.h"
#include "vol_util.h"
#include "mcell_structs.h"
#include "mdlparse_util.h"
#include "mdlparse_aux.h"
#include "react_output.h"
#include "react_util.h"
#include "macromolecule.h"
#include "diffuse_util.h"
#include "create_species.h"
#include "create_release_site.h"
#include "create_geometry.h"
#include "create_object.h"
#include "libmcell.h"

extern void chkpt_signal_handler(int sn);

/* Free a variable value, leaving the symbol free for reassignment to another
 * type. */
static int mdl_free_variable_value(struct mdlparse_vars *parse_state,
                                   struct sym_table *sym);

/* Free a numeric expression list, deallocating all items. */
static void mdl_free_numeric_list(struct num_expr_list *nel);

/*************************************************************************
 double_cmp:
    Comparison function for doubles, to be passed to qsort.

 In:  i1: first item for comparison
      i2: second item for comparison
 Out: -1, 0, or 1 if the first item is less than, equal to, or greater than the
      second, resp.
*************************************************************************/
static int
double_cmp(void const *i1, void const *i2)
{
  double const *d1 = (double const *) i1;
  double const *d2 = (double const *) i2;
  if (*d1 < *d2)
    return -1;
  else if (*d1 > *d2)
    return 1;
  else
    return 0;
}

/*************************************************************************
 mdl_strip_quotes:
    Strips the quotes from a quoted string.

 In:  in: string to strip
 Out: stripped string, or NULL on error
*************************************************************************/
char *
mdl_strip_quotes(char *in)
{

  char *q = strip_quotes(in);
  free(in);
  if (q != NULL)
    return q;
  else
  {
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
char *
mdl_strcat(char *s1, char *s2)
{

  char *cat = my_strcat(s1, s2);
  free(s1);
  free(s2);
  if (cat != NULL)
    return cat;
  else
  {
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
char *
mdl_strdup(char const *s1)
{

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
void mdl_warning(struct mdlparse_vars *parse_state, char const *fmt, ...)
{
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

/************************************************************************
 mdl_find_include_file:
      Find the path for an include file.  For an absolute include path, the
      file path is unmodified, but for a relative path, the resultant path will
      be relative to the currently parsed file.

 In:  path: path from include statement
      cur_path: path of current include file
 Out: allocated buffer containing path of the include file, or NULL if the file
      path couldn't be allocated.  If we ever use a more complex mechanism for
      locating include files, we may also return NULL if no file could be
      located.
 ***********************************************************************/
char *
mdl_find_include_file(struct mdlparse_vars *parse_state,
                      char const *path,
                      char const *cur_path)
{
  char *candidate = NULL;
  if (path[0] == '/')
    candidate = mdl_strdup(path);
  else
  {
    char *last_slash = strrchr(cur_path, '/');
#ifdef _WIN32
    char *last_bslash = strrchr(cur_path, '\\');
    if (last_bslash > last_slash)
      last_slash = last_bslash;
#endif
    if (last_slash == NULL)
      candidate = mdl_strdup(path);
    else
      candidate = CHECKED_SPRINTF("%.*s/%s",
                                  (int) (last_slash - cur_path),
                                  cur_path,
                                  path);
  }

  return candidate;
}

/*************************************************************************
 load_rate_file:
    Read in a time-varying reaction rates file.

 In:  parse_state:  parser state
      rx:    Reaction structure that we'll load the rates into.
      fname: Filename to read the rates from.
      path:  Index of the pathway that these rates apply to.
 Out: Returns 1 on error, 0 on success.
      Rates are added to the prob_t linked list.  If there is a rate given for
      time <= 0, then this rate is stuck into cum_probs and the (time <= 0)
      entries are not added to the list.  If no initial rate is given in the
      file, it is assumed to be zero.
 Note: The file format is assumed to be two columns of numbers; the first
      column is time (in seconds) and the other is rate (in appropriate
      units) that starts at that time.  Lines that are not numbers are
      ignored.
*************************************************************************/
#define RATE_SEPARATORS "\f\n\r\t\v ,;"
#define FIRST_DIGIT "+-0123456789"
static int
load_rate_file(struct mdlparse_vars *parse_state,
               struct rxn *rx,
               char *fname,
               int path)
{
  int i;
  FILE *f = fopen(fname,"r");

  if (!f) return 1;
  else
  {
    struct t_func *tp,*tp2;
    double t,rate;
    char buf[2048];
    char *cp;
    int linecount = 0;
#ifdef DEBUG
    int valid_linecount = 0;
#endif

    tp2 = NULL;
    while (fgets(buf,2048,f))
    {
      linecount++;
      for (i=0;i<2048;i++) { if (!strchr(RATE_SEPARATORS,buf[i])) break; }

      if (i<2048 && strchr(FIRST_DIGIT,buf[i]))
      {
        t = strtod((buf+i) , &cp);
        if (cp == (buf+i)) continue;  /* Conversion error. */

        for (i=cp-buf ; i<2048 ; i++) { if (!strchr(RATE_SEPARATORS,buf[i])) break; }
        rate = strtod((buf+i) , &cp);
        if (cp == (buf+i)) continue;  /* Conversion error */

        /* at this point we need to handle negative reaction rates */
        if (rate < 0.0)
        {
          if (parse_state->vol->notify->neg_reaction==WARN_ERROR)
          {
            mdlerror(parse_state, "Error: reaction rates should be zero or positive.");
            return 1;
          }
          else if (parse_state->vol->notify->neg_reaction == WARN_WARN) {
            mcell_warn("Warning: negative reaction rate %f; setting to zero and continuing.", rate);
            rate = 0.0;
          }
        }


        tp = CHECKED_MEM_GET(parse_state->vol->tv_rxn_mem, "time-varying reaction rate");
        if (tp == NULL)
          return 1;
        tp->next = NULL;
        tp->path = path;
        tp->time = t / parse_state->vol->time_unit;
        tp->value = rate;
#ifdef DEBUG
        valid_linecount++;
#endif

        if (rx->prob_t == NULL)
        {
          rx->prob_t = tp;
          tp2 = tp;
        }
        else
        {
          if (tp2==NULL)
          {
            tp2 = tp;
            tp->next = rx->prob_t;
            rx->prob_t = tp;
          }
          else
          {
            if (tp->time < tp2->time)
              mcell_warn("In rate file '%s', line %d is out of sequence.  Resorting.", fname, linecount);
            tp->next = tp2->next;
            tp2->next = tp;
            tp2 = tp;
          }
        }
      }
    }

#ifdef DEBUG
    mcell_log("Read %d rates from file %s.", valid_linecount, fname);
#endif

    fclose(f);
  }
  return 0;
}
#undef FIRST_DIGIT
#undef RATE_SEPARATORS

/**************************************************************************
 mdl_valid_file_mode:
    Check that the speficied file mode string is valid for an fopen statement.

 In: parse_state: parser state
     mode: mode string to check
 Out: 1 if
**************************************************************************/
int
mdl_valid_file_mode(struct mdlparse_vars *parse_state, char *mode)
{
  char c = mode[0];
  if (c != 'r' && c != 'w' && c != 'a')
  {
    mdlerror_fmt(parse_state, "Invalid file mode: %s", mode);
    return 1;
  }
  if (c == 'r')
  {
    mdlerror(parse_state, "MCell models do not currently support opening files for reading");
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
int
mdl_expr_log(struct mdlparse_vars *parse_state, double in, double *out)
{
  if (in <= 0.0)
  {
    mdlerror_fmt(parse_state, "Attempt to compute LOG(%.15g), which is not defined.", in);
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
int
mdl_expr_log10(struct mdlparse_vars *parse_state, double in, double *out)
{
  if (in <= 0.0)
  {
    mdlerror_fmt(parse_state, "Attempt to compute LOG10(%.15g), which is not defined.", in);
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
int
mdl_expr_mod(struct mdlparse_vars *parse_state,
             double in,
             double divisor,
             double *out)
{
  if (divisor == 0.0)
  {
    mdlerror_fmt(parse_state, "Attempt to compute MOD(%.15g, 0.0), which is not defined.", in);
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
int
mdl_expr_div(struct mdlparse_vars *parse_state,
             double in,
             double divisor,
             double *out)
{
  if (divisor == 0.0)
  {
    mdlerror_fmt(parse_state,
                 "Attempt to divide %.15g by zero",
                 in);
    return 1;
  }

  *out = in / divisor;
  if (isinf(*out))
  {
    mdlerror_fmt(parse_state,
                 "Cannot compute %.15g / %.15g: result is too large",
                 in,
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
int
mdl_expr_pow(struct mdlparse_vars *parse_state,
             double in,
             double exponent,
             double *out)
{
  if (in < 0.0)
  {
    if (exponent != (int) exponent)
    {
      mdlerror_fmt(parse_state,
                   "Cannot compute %.15g^%.15g: negative value raised to non-integral power would be complex",
                   in,
                   exponent);
      return 1;
    }
  }

  *out = pow(in, exponent);
  if (isinf(*out))
  {
    mdlerror_fmt(parse_state,
                 "Cannot compute %.15g^%.15g: result is too large",
                 in,
                 exponent);
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
double
mdl_expr_rng_uniform(struct mdlparse_vars *parse_state)
{
  return rng_dbl(parse_state->vol->rng);
}

/**************************************************************************
 mdl_expr_roundoff:
    Round the input value off to n significant figures.

 In:  in:   value to round off
      ndigits: number of digits
 Out: value rounded to n digits
**************************************************************************/
double
mdl_expr_roundoff(double in, int ndigits)
{
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
int
mdl_expr_string_to_double(struct mdlparse_vars *parse_state,
                          char *str,
                          double *out)
{
  *out = strtod(str, (char **)NULL);
  if (errno == ERANGE)
  {
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
struct sym_table *
mdl_new_filehandle(struct mdlparse_vars *parse_state,
                   char *name)
{
  struct sym_table *sym;
  sym = retrieve_sym(name, parse_state->vol->fstream_sym_table);

  /* If this file is already open, close it. */
  if (sym != NULL)
  {
    struct file_stream *filep = (struct file_stream *) sym->value;
    if (filep->stream != NULL)
    {
      fclose(filep->stream);
      free(filep->name);
      filep->stream = NULL;
      filep->name = NULL;
    }
  }

  /* Otherwise, create it */
  else if ((sym = store_sym(name, FSTRM, parse_state->vol->fstream_sym_table, NULL)) == NULL)
  {
    mdlerror_fmt(parse_state, "Out of memory while creating file stream: %s", name);
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
int
mdl_fopen(struct mdlparse_vars *parse_state,
          struct sym_table *filesym,
          char *name,
          char *mode)
{
  struct file_stream *filep = (struct file_stream *) filesym->value;
  filep->name = name;
  if ((filep->stream = fopen(filep->name, mode)) == NULL)
  {
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
int
mdl_fclose(struct mdlparse_vars *parse_state, struct sym_table *filesym)
{
  struct file_stream *filep = (struct file_stream *) filesym->value;
  if (filep->stream == NULL)
    return 0;

  if (fclose(filep->stream) != 0)
  {
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
static double *
double_dup(double value)
{

  double *dup_value;
  dup_value = CHECKED_MALLOC_STRUCT(double, "numeric variable");
  if (dup_value != NULL)
    *dup_value=value;
  return dup_value;
}

/*************************************************************************
 mdl_new_printf_arg_double:
    Create a new double argument for a printf argument list.

 In: parse_state: parser state
     d: value for argument
 Out: The new argument, or NULL if an error occurred.
*************************************************************************/
struct arg *
mdl_new_printf_arg_double(struct mdlparse_vars *parse_state, double d)
{
  struct arg *a = CHECKED_MALLOC_STRUCT(struct arg, "format argument");
  if (a == NULL)
    return NULL;

  if ((a->arg_value = double_dup(d)) == NULL)
  {
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
struct arg *
mdl_new_printf_arg_string(char const *str)
{

  struct arg *a = CHECKED_MALLOC_STRUCT(struct arg, "format argument");
  if (a == NULL)
    return NULL;

  a->arg_value = (void *) str;
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
static void
mdl_free_printf_arg_list(struct arg *args)
{
  while (args != NULL)
  {
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
char *
mdl_expand_string_escapes(char const *in)
{

  char *out;

  /* Allocate buffer for new string */
  char *expanded = CHECKED_MALLOC(strlen(in)+1, "printf format string");
  if (expanded == NULL)
    return NULL;

  /* Copy expanded string to new buffer */
  out = expanded;
  while (*in != '\0')
  {
    if (*in != '\\')
      *out ++ = *in ++;
    else
    {
      ++ in;
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
      else if (*in == 'x')
      {
        if (isxdigit(in[1]) && isxdigit(in[2]))
        {
          char buffer[3];
          buffer[0] = in[1];
          buffer[1] = in[2];
          buffer[2] = '\0';
          *out++ = (char) strtol(buffer, NULL, 16);
          if (out[-1] == '\0')
            out[-1] = ' ';
          in += 2;
        }
        else
          *out++ = 'x';
      }
      else if (*in == '0')
      {
        if ('0' <= in[1] && in[1] <= '7'  &&
            '0' <= in[2] && in[2] <= '7'  &&
            '0' <= in[3] && in[3] <= '7')
        {
          char buffer[4];
          buffer[0] = in[1];
          buffer[1] = in[2];
          buffer[2] = in[3];
          buffer[3] = '\0';
          *out++ = (char) strtol(buffer, NULL, 8);
          if (out[-1] == '\0')
            out[-1] = ' ';
          in += 3;
        }
        else
          *out++ = '0';
      }
      else if (*in == '\\')
        *out++ = '\\';
      else if (*in == '"')
        *out++ = '"';
      else if (*in == '\0')
      {
        *out++ = '\\';
        break;
      }
      else
        *out++ = *in;
      ++ in;
    }
  }
  *out++ = '\0';
  return expanded;
}

/*************************************************************************
 Internal code used to designate the type required to fulfill a particular
 format specifier within a printf-style format string.
*************************************************************************/
enum
{
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
static int
get_printf_conversion_specifier(struct mdlparse_vars *parse_state,
                                char const *fmt,
                                int *num_asterisks)
{
  static char CONVERSION_SPECIFIERS[] = "diouxXeEfFgGaAcCsSpnm";

  char const * const fmt_orig = fmt;
  char length = '\0';
  *num_asterisks = 0;
  for (++ fmt; *fmt != '\0'; ++ fmt)
  {
    if (isalpha(*fmt)  &&  strchr(CONVERSION_SPECIFIERS, *fmt) != NULL)
      break;

    /* Scan any options in the string, bombing out if we find something
     * invalid, and keeping track of the length specifiers and asterisks we
     * find. */
    switch (*fmt)
    {
      case 'l':
        if (length == '\0')
          length = 'l';
        else if (length == 'l')
          length = 'L';
        else
        {
          mdlerror_fmt(parse_state,
                       "Format string segment '%s' has an invalid length modifier inside a conversion specification.",
                       fmt_orig);
          return PRINTF_INVALID;
        }
        break;

      case 'h':
        if (length == '\0')
          length = 'h';
        else if (length == 'h')
          length = 'H';
        else
        {
          mdlerror_fmt(parse_state,
                       "Format string segment '%s' has an invalid length modifier inside a conversion specification.",
                       fmt_orig);
          return PRINTF_INVALID;
        }
        break;

      case 'L':
      case 'q':
      case 'j':
      case 'z':
      case 't':
        mdlerror_fmt(parse_state,
                     "Format string segment '%s' has an unsupported length modifier (%c) inside a conversion specification.",
                     fmt_orig,
                     *fmt);
        return PRINTF_INVALID;

      case '*':
        ++ *num_asterisks;
        if (*num_asterisks > 2)
        {
          mdlerror_fmt(parse_state,
                       "Format string segment '%s' has more than two asterisks inside a conversion specification, which is unsupported by MCell.",
                       fmt_orig);
          return PRINTF_INVALID;
        }
        break;

      case '$':
        mdlerror_fmt(parse_state,
                     "Format string segment '%s' has a '$' inside a conversion specification, which is unsupported by MCell.",
                     fmt_orig);
        return PRINTF_INVALID;

      default:
        /* Skip this char. */
        break;
    }
  }

  /* Filter invalid types */
  if (*fmt == '\0')
  {
    mdlerror_fmt(parse_state,
                 "Format string segment '%s' contains no conversion specifier.",
                 fmt_orig);
    return PRINTF_INVALID;
  }
  else if (*fmt == 'n')
  {
    mdlerror_fmt(parse_state,
                 "Format string segment '%s' attempts to use dangerous conversion specifier 'n'.",
                 fmt_orig);
    return PRINTF_INVALID;
  }

  /* Now, handle the format specifier itself */
  switch (*fmt)
  {
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
      else
      {
        mdlerror_fmt(parse_state,
                     "Format string segment '%s' uses unsupported length specifier '%c' for integer.",
                     fmt_orig,
                     length);
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
      else
      {
        mdlerror_fmt(parse_state,
                     "Format string segment '%s' uses unsupported length specifier '%c' for unsigned integer.",
                     fmt_orig,
                     length);
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
      else
      {
        mdlerror_fmt(parse_state,
                     "Format string segment '%s' uses unsupported length specifier '%c' for double-precision floating point.",
                     fmt_orig,
                     length);
        return PRINTF_INVALID;
      }

    case 'c':
      if (length == '\0')
        return PRINTF_SIGNED_CHAR;
      else
      {
        mdlerror_fmt(parse_state,
                     "Format string segment '%s' uses unsupported length specifier '%c' for char.",
                     fmt_orig,
                     length);
        return PRINTF_INVALID;
      }

    case 's':
      if (length == '\0')
        return PRINTF_STRING;
      else
      {
        mdlerror_fmt(parse_state,
                     "Format string segment '%s' uses unsupported length specifier '%c' for string.",
                     fmt_orig,
                     length);
        return PRINTF_INVALID;
      }

    case 'C':
    case 'S':
    case 'p':
    case 'm':
    default:
      mdlerror_fmt(parse_state,
                   "Format string segment '%s' uses unsupported conversion specifier '%c'.",
                   fmt_orig,
                   *fmt);
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
static char *
my_sprintf_segment(struct mdlparse_vars *parse_state,
                   char *fmt_seg,
                   struct arg **argpp)
{
  int num_asterisks = 0;
  struct arg *argp = *argpp;

  int spec_type = get_printf_conversion_specifier(parse_state, fmt_seg, &num_asterisks);
  if (spec_type == PRINTF_INVALID)
    return NULL;

  /* Pull out arguments for asterisks */
  int ast1=0, ast2=0;
  switch (num_asterisks)
  {
    case 0:
      break;

    case 2:
      if (argp->arg_type == DBL)
        ast2 = (int) (*(double *) (argp->arg_value) + EPS_C);
      else
      {
        mdlerror_fmt(parse_state,
                     "Argument to be consumed by '*' in conversion specification segment '%s' is not numeric",
                     fmt_seg);
        return NULL;
      }
      if ((argp = argp->next) == NULL)
      {
        mdlerror_fmt(parse_state,
                     "Too few arguments for conversion specification segment '%s'",
                     fmt_seg);
        return NULL;
      }
      /* FALL THROUGH */

    case 1:
      if (argp->arg_type == DBL)
        ast1 = (int) (*(double *) (argp->arg_value) + EPS_C);
      else
      {
        mdlerror_fmt(parse_state,
                     "Argument to be consumed by '*' in conversion specification segment '%s' is not numeric",
                     fmt_seg);
        return NULL;
      }
      if ((argp = argp->next) == NULL)
      {
        mdlerror_fmt(parse_state,
                     "Too few arguments for conversion specification segment '%s'",
                     fmt_seg);
        return NULL;
      }
      break;

    default:
      mdlerror_fmt(parse_state,
                   "Invalid asterisk count in conversion specification segment '%s'",
                   fmt_seg);
      return NULL;
  }

  /* Type check */
  if (argp->arg_type == STR)
  {
    if (spec_type != PRINTF_STRING       &&
        spec_type != PRINTF_SIGNED_CHAR  &&
        spec_type != PRINTF_UNSIGNED_CHAR)
    {
      mdlerror_fmt(parse_state,
                   "Argument in conversion specification segment '%s' is a string, but a numeric value is required",
                   fmt_seg);
      return NULL;
    }
  }
  else
  {
    if (spec_type == PRINTF_STRING)
    {
      mdlerror_fmt(parse_state,
                   "Argument in conversion specification segment '%s' is numeric, but a string value is required",
                   fmt_seg);
      return NULL;
    }
  }

  /* Now, convert! */
  int arg_value;
  char *formatted = NULL;
  switch (spec_type)
  {
    case PRINTF_INVALID:
      return NULL;

    case PRINTF_STRING:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (char *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (char *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (char *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_SIGNED_CHAR:
    case PRINTF_UNSIGNED_CHAR:
      if (argp->arg_type == STR)
      {
        if (strlen((char *) argp->arg_value) > 1)
        {
          mdlerror_fmt(parse_state,
                       "Argument in conversion specification segment '%s' has too many characters",
                       fmt_seg);
          return NULL;
        }

        if (spec_type == PRINTF_SIGNED_CHAR)
          arg_value = * (signed char *) argp->arg_value;
        else
          arg_value = * (unsigned char *) argp->arg_value;
      }
      else
      {
        if (spec_type == PRINTF_SIGNED_CHAR)
          arg_value = (signed char) *(double *) argp->arg_value;
        else
          arg_value = (unsigned char) *(double *) argp->arg_value;
      }

      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_DOUBLE:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_INT:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (int) *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (int) *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (int) *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_LONG_INT:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (long int) *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (long int) *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (long int) *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_LONG_LONG_INT:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (long long int) *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (long long int) *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (long long int) *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_SHORT_INT:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (short int) *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (short int) *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (short int) *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_U_INT:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (unsigned int) *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (unsigned int) *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (unsigned int) *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_U_LONG_INT:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (unsigned long int) *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (unsigned long int) *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (unsigned long int) *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_U_LONG_LONG_INT:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (unsigned long long int) *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (unsigned long long int) *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (unsigned long long int) *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    case PRINTF_U_SHORT_INT:
      switch (num_asterisks)
      {
        case 0: formatted = alloc_sprintf(fmt_seg, (unsigned short int) *(double *) argp->arg_value); break;
        case 1: formatted = alloc_sprintf(fmt_seg, ast1, (unsigned short int) *(double *) argp->arg_value); break;
        case 2: formatted = alloc_sprintf(fmt_seg, ast2, ast1, (unsigned short int) *(double *) argp->arg_value); break;
        default: UNHANDLED_CASE(num_asterisks);
      }
      break;

    default: UNHANDLED_CASE(spec_type);
  }

  if (formatted == NULL)
    mcell_allocfailed("Failed to format a string for an MDL printf/fprintf/spritnf statement.");
  else
    *argpp = argp->next;
  return formatted;
}

/*************************************************************************
 Temporary structure used by sprintf implementation.
*************************************************************************/
struct sprintf_output_list
{
  struct sprintf_output_list *next;     /* Next item in list */
  size_t len;                           /* Length of this item */
  char *segment;                        /* Data for this item */
};

/*************************************************************************
 my_sprintf:
    sprintf-like formatting of MDL arguments.

 In:  parse_state: parser state
      format: string to expand
      argp: argument list
 Out: 0 on success, 1 on failure
*************************************************************************/
static char *
my_sprintf(struct mdlparse_vars *parse_state,
           char *format,
           struct arg *argp)
{
  char *this_start, *this_end;
  char const * const format_orig = format;

  struct sprintf_output_list head, *tail;
  head.segment = NULL;
  head.next = NULL;
  head.len = 0u;
  tail = &head;

  /* Find the start of the first format code */
  this_start = strchr(format, '%');
  while (this_start != NULL  &&  this_start[1] == '%')
    this_start = strchr(this_start+2, '%');

  /* If we have text before the first conversion specification, copy it. */
  if (this_start != NULL)
  {
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
  while (this_start != NULL)
  {
    /* If we've run out of arguments... */
    if (argp == NULL)
    {
      mdlerror_fmt(parse_state,
                   "Insufficient arguments for printf-style format string '%s'",
                   format_orig);
      goto failure;
    }

    /* Find the start of the first format code */
    this_end = strchr(this_start + 1, '%');
    while (this_end != NULL  &&  this_end[1] == '%')
      this_end = strchr(this_end+2, '%');

    /* If we have only a single format code, do the rest of the string */
    if (this_end == NULL)
    {
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
    else
    {
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
      while (this_start != NULL  &&  this_start[1] == '%')
        this_start = strchr(this_start+2, '%');
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
  while (head.next != NULL)
  {
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
  while (head.next != NULL)
  {
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
static int
my_fprintf(struct mdlparse_vars *parse_state,
           FILE *outfile,
           char *format,
           struct arg *argp)
{
  char *str = my_sprintf(parse_state, format, argp);
  if (str == NULL)
    return 1;

  if (fputs(str, outfile) < 0)
  {
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
static void
defang_format_string(char *fmt)
{
  while (*fmt != '\0')
  {
    if (iscntrl(*fmt))
      *fmt = '.';
    else if (! isascii(*fmt))
      *fmt = '.';
    ++ fmt;
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
int
mdl_printf(struct mdlparse_vars *parse_state,
           char *fmt,
           struct arg *arg_head)
{
  if (my_fprintf(parse_state, mcell_get_log_file(), fmt, arg_head))
  {
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
int
mdl_fprintf(struct mdlparse_vars *parse_state,
            struct file_stream *filep,
            char *fmt,
            struct arg *arg_head)
{
  if (my_fprintf(parse_state, filep->stream, fmt, arg_head))
  {
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
char *
mdl_string_format(struct mdlparse_vars *parse_state, char *fmt, struct arg *arg_head)
{
  char *str = my_sprintf(parse_state, fmt, arg_head);
  if (str == NULL)
  {
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
int
mdl_sprintf(struct mdlparse_vars *parse_state,
            struct sym_table *assign_var,
            char *fmt,
            struct arg *arg_head)
{
  char *str = my_sprintf(parse_state, fmt, arg_head);
  if (str == NULL)
  {
    mdl_free_printf_arg_list(arg_head);
    free(fmt);
    mdlerror_fmt(parse_state, "Could not sprintf to variable: %s", assign_var->name);
    return 1;
  }
  free(fmt);
  mdl_free_printf_arg_list(arg_head);

  /* If the symbol had a value, try to free it */
  if (assign_var->value  &&  mdl_free_variable_value(parse_state, assign_var))
  {
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
int
mdl_fprint_time(struct mdlparse_vars *parse_state,
                struct sym_table *filep_sym,
                char *fmt)
{
  char time_str[128];
  time_t the_time;
  struct file_stream *filep = (struct file_stream *) filep_sym->value;

  the_time = time(NULL);
  strftime(time_str, 128, fmt, localtime(&the_time));
  free(fmt);
  if (fprintf(filep->stream, "%s", time_str) == EOF)
  {
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
void
mdl_print_time(struct mdlparse_vars *parse_state, char *fmt)
{
  char time_str[128];
  time_t the_time = time(NULL);
  strftime(time_str, 128, fmt, localtime(&the_time));
  free(fmt);
  if (parse_state->vol->procnum == 0) fprintf(mcell_get_log_file(), "%s", time_str);
}

/************************************************************************
 swap_double:
 In:  x, y: doubles to swap
 Out: Swaps references to two double values.
 ***********************************************************************/
static void
swap_double(double *x, double *y)
{
  double temp;

  temp=*x;
  *x=*y;
  *y=temp;
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
int
mdl_generate_range(struct mdlparse_vars *parse_state,
                   struct num_expr_list_head *list,
                   double start,
                   double end,
                   double step)
{
  list->value_head  = NULL;
  list->value_tail  = NULL;
  list->value_count = 0;
  list->shared = 0;

  if (step == 0.0)
  {
    mdlerror(parse_state, "A step size of 0 was requested, which would generate an infinite list.");
    return 1;
  }

  /* Above a certain point, we probably don't want this.  If we need arrays
   * this large, we need to reengineer this. */
  if (fabs((end - start) / step) > 100000000.)
  {
    mdlerror(parse_state, "A range was requested that encompasses too many values (maximum 100,000,000)");
    return 1;
  }

  if (step > 0)
  {
    /* JW 2008-03-31: In the guard on the loop below, it seems to me that
     * the third condition is redundant with the second.
     */
    for (double tmp_dbl = start;
         tmp_dbl < end                           ||
         ! distinguishable(tmp_dbl, end, EPS_C)  ||
         fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step)
    {
      struct num_expr_list *nel;
      nel = CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric list");
      if (nel == NULL)
      {
        mdl_free_numeric_list(list->value_head);
        list->value_head = list->value_tail = NULL;
        return 1;
      }
      nel->value = tmp_dbl;
      nel->next = NULL;

      ++ list->value_count;
      if (list->value_tail != NULL)
        list->value_tail->next = nel;
      else
        list->value_head = nel;
      list->value_tail = nel;
    }
  }
  else /* if (step < 0) */
  {
    /* JW 2008-03-31: In the guard on the loop below, it seems to me that
     * the third condition is redundant with the second.
     */
    for (double tmp_dbl = start;
         tmp_dbl > end                           ||
         ! distinguishable(tmp_dbl, end, EPS_C)  ||
         fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step)
    {
      struct num_expr_list *nel;
      nel = CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric list");
      if (nel == NULL)
      {
        mdl_free_numeric_list(list->value_head);
        list->value_head = list->value_tail = NULL;
        return 1;
      }
      nel->value = tmp_dbl;
      nel->next = NULL;

      ++ list->value_count;
      if (list->value_tail != NULL)
        list->value_tail->next = nel;
      else
        list->value_head = nel;
      list->value_tail = nel;
    }
  }

  return 0;
}

/*************************************************************************
 mdl_add_range_value:
    Add a value to a numeric list.

 In:  parse_state: parser state
      lh:   list to receive value
      value: value for list
 Out: 0 on success, 1 on failure
*************************************************************************/
int
mdl_add_range_value(struct mdlparse_vars *parse_state,
                    struct num_expr_list_head *lh,
                    double value)
{
  if (lh->value_head == NULL)
    return mdl_generate_range_singleton(lh, value);

  struct num_expr_list *nel = CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric array");
  if (nel == NULL)
    return 1;
  nel->next = NULL;
  nel->value = value;
  lh->value_tail = lh->value_tail->next = nel;
  ++ lh->value_count;
  return 0;
}

/*************************************************************************
 mdl_generate_range_singleton:
    Generate a numeric list containing a single value.

 In:  lh:   list to receive value
      value: value for list
 Out: 0 on success, 1 on failure
*************************************************************************/
int
mdl_generate_range_singleton(struct num_expr_list_head *lh,
                             double value)
{

  struct num_expr_list *nel = CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric array");
  if (nel == NULL)
    return 1;
  lh->value_head = lh->value_tail = nel;
  lh->value_count = 1;
  lh->shared = 0;
  lh->value_head->value = value;
  lh->value_head->next = NULL;
  return 0;
}

/*************************************************************************
 mdl_copysort_numeric_list:
    Copy and sort a num_expr_list in ascending numeric order.

 In:  parse_state:  parser state
      head:  the list to sort
 Out: list is sorted
*************************************************************************/
static struct num_expr_list *
mdl_copysort_numeric_list(struct mdlparse_vars *parse_state,
                          struct num_expr_list *head)
{
  struct num_expr_list_head new_head;
  if (mdl_generate_range_singleton(&new_head, head->value))
    return NULL;

  head = head->next;
  while (head != NULL)
  {
    struct num_expr_list *insert_pt, **prev;
    for (insert_pt = new_head.value_head, prev = &new_head.value_head;
         insert_pt != NULL;
         prev = &insert_pt->next, insert_pt = insert_pt->next)
    {
      if (insert_pt->value >= head->value)
        break;
    }

    struct num_expr_list *new_item = CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric array");
    if (new_item == NULL)
    {
      mdl_free_numeric_list(new_head.value_head);
      return NULL;
    }

    new_item->next = insert_pt;
    new_item->value = head->value;
    *prev = new_item;
    if (insert_pt == NULL)
      new_head.value_tail = new_item;
    head = head->next;
  }

  return new_head.value_head;
}

/*************************************************************************
 mdl_sort_numeric_list:
    Sort a num_expr_list in ascending numeric order.  N.B. This uses bubble
    sort, which is O(n^2).  Don't use it if you expect your list to be very
    long.  The list is sorted in-place.

 In:  head:  the list to sort
 Out: list is sorted
*************************************************************************/
static void
mdl_sort_numeric_list(struct num_expr_list *head)
{
  struct num_expr_list *curr,*next;
  int done = 0;
  while (! done)
  {
    done = 1;
    curr = head;
    while (curr != NULL)
    {
      next = curr->next;
      if (next != NULL)
      {
        if (curr->value > next->value)
        {
          done = 0;
          swap_double(&curr->value, &next->value);
        }
      }
      curr = next;
    }
  }
}

/*************************************************************************
 mdl_free_numeric_list:
    Free a num_expr_list.

 In:  nel:  the list to free
 Out: all elements are freed
*************************************************************************/
static void
mdl_free_numeric_list(struct num_expr_list *nel)
{
  while (nel != NULL)
  {
    struct num_expr_list *n = nel;
    nel = nel->next;
    free(n);
  }
}

#ifdef DEBUG
/*************************************************************************
 mdl_debug_dump_array:
    Display a human-readable representation of the array to stdout.

 In:  nel:  the list to free
 Out: displayed to stdout
*************************************************************************/
void
mdl_debug_dump_array(struct num_expr_list *nel)
{
  struct num_expr_list *elp;
  no_printf("\nArray expression: [");
  for (elp = nel; elp != NULL; elp = elp->next)
  {
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
static double *
num_expr_list_to_array(struct num_expr_list_head *lh,
                       int *count)
{

  double *arr = NULL, *ptr = NULL;
  struct num_expr_list *i;

  if (count != NULL)
    *count = lh->value_count;

  if (lh->value_count == 0)
    return NULL;

  arr = CHECKED_MALLOC_ARRAY(double, lh->value_count, "numeric array");
  if (arr != NULL)
  {
    for (i = (struct num_expr_list *) lh->value_head, ptr = arr; i != NULL; i = i->next)
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
struct vector3 *
mdl_point(struct mdlparse_vars *parse_state,
          struct num_expr_list_head *vals)
{
  struct vector3 *vec;
  if (vals->value_count != 3)
  {
    mdlerror(parse_state, "Three dimensional value required");
    return NULL;
  }

  vec = CHECKED_MALLOC_STRUCT(struct vector3, "3-D vector");
  if (! vec)
    return NULL;

  vec->x = vals->value_head->value;
  vec->y = vals->value_head->next->value;
  vec->z = vals->value_tail->value;
  if (! vals->shared)
    mdl_free_numeric_list(vals->value_head);
  return vec;
}

/*************************************************************************
 mdl_point_scalar:
    Create a 3-D vector equal to s*[1, 1, 1] for some scalar s.

 In:  val: scalar
 Out: a 3-D vector, or NULL if an error occurs.
*************************************************************************/
struct vector3 *
mdl_point_scalar(double val)
{

  struct vector3 *vec;
  vec = CHECKED_MALLOC_STRUCT(struct vector3, "3-D vector");
  if (! vec)
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
static int
mdl_free_variable_value(struct mdlparse_vars *parse_state,
                        struct sym_table *sym)
{
  switch (sym->sym_type)
  {
    case DBL:
    case STR:
      free(sym->value);
      break;

    case ARRAY:
      /* XXX: For the moment, we can't free these since they might be shared by
       * two different variables... */
      break;

    default:
      mdlerror_fmt(parse_state, "Internal error: Attempt to free variable value of incorrect type: %s", sym->name);
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
struct sym_table *
mdl_get_or_create_variable(struct mdlparse_vars *parse_state,
                           char *name)
{
  /* Attempt to fetch existing variable */
  struct sym_table *st = NULL;
  if ((st=retrieve_sym(name, parse_state->vol->var_sym_table)) != NULL)
  {
    free(name);
    return st;
  }

  /* Create the variable */
  if ((st = store_sym(name, TMP, parse_state->vol->var_sym_table, NULL)) == NULL)
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
int
mdl_assign_variable_double(struct mdlparse_vars *parse_state,
                           struct sym_table *sym,
                           double value)
{
  /* If the symbol had a value, try to free it */
  if (sym->value  &&  mdl_free_variable_value(parse_state, sym))
    return 1;

  /* Allocate space for the new value */
  if ((sym->value = CHECKED_MALLOC_STRUCT(double, "numeric variable")) == NULL)
    return 1;

  /* Update the value */
  *(double *) sym->value = value;
  sym->sym_type = DBL;
  no_printf("\n%s is equal to: %f\n",sym->name,value);
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
int
mdl_assign_variable_string(struct mdlparse_vars *parse_state,
                           struct sym_table *sym,
                           char const *value)
{
  /* If the symbol had a value, try to free it */
  if (sym->value  &&  mdl_free_variable_value(parse_state, sym))
    return 1;

  /* Allocate space for the new value */
  if ((sym->value = mdl_strdup(value)) == NULL)
    return 1;

  /* Update the value */
  sym->sym_type = STR;
  no_printf("\n%s is equal to: %s\n",sym->name,value);
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
int
mdl_assign_variable_array(struct mdlparse_vars *parse_state,
                          struct sym_table *sym,
                          struct num_expr_list *value)
{
  /* If the symbol had a value, try to free it */
  if (sym->value  &&  mdl_free_variable_value(parse_state, sym))
    return 1;

  sym->sym_type = ARRAY;
  sym->value = value;
#ifdef DEBUG
  struct num_expr_list *elp=(struct num_expr_list *) sym->value;
  no_printf("\n%s is equal to: [", sym->name);
  while (elp!=NULL) 
  {
    no_printf("%lf",elp->value);
    elp=elp->next;
    if (elp!=NULL)
    {
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
int
mdl_assign_variable(struct mdlparse_vars *parse_state,
                    struct sym_table *sym,
                    struct sym_table *value)
{
  switch (value->sym_type)
  {
    case DBL:
      if (mdl_assign_variable_double(parse_state, sym, *(double *) value->value))
        return 1;
      break;

    case STR:
      if (mdl_assign_variable_string(parse_state, sym, (char const *) value->value))
        return 1;
      break;

    case ARRAY:
      if (mdl_assign_variable_array(parse_state, sym, (struct num_expr_list *) value->value))
        return 1;
      break;

    default: UNHANDLED_CASE(value->sym_type);
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
void
mdl_set_all_notifications(struct volume *vol,
                          byte notify_value)
{
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
int
mdl_set_iteration_report_freq(struct mdlparse_vars *parse_state,
                              long long interval)
{
  /* Only if not set on command line */
  if (parse_state->vol->log_freq == ULONG_MAX)
  {
    if (interval < 1)
    {
      mdlerror(parse_state, "Invalid iteration number reporting interval: use value >= 1");
      return 1;
    }
    else
    {
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
void
mdl_set_all_warnings(struct volume *vol, byte warning_level)
{
  vol->notify->neg_diffusion = warning_level;
  vol->notify->neg_reaction = warning_level;
  vol->notify->high_reaction_prob = warning_level;
  vol->notify->close_partitions = warning_level;
  vol->notify->degenerate_polys = warning_level;
  vol->notify->overwritten_file = warning_level;
  vol->notify->complex_placement_failure = warning_level;
  vol->notify->mol_placement_failure = warning_level;

  if (warning_level==WARN_ERROR) warning_level=WARN_WARN;
  vol->notify->short_lifetime = warning_level;
  vol->notify->missed_reactions = warning_level;
  vol->notify->missed_surf_orient = warning_level;
  vol->notify->useless_vol_orient = warning_level;
  vol->notify->invalid_output_step_time = warning_level;
}

/*************************************************************************
 mdl_set_lifetime_warning_threshold:
    Set the lifetime warning threshold.

 In:  parse_state: parser state
      lifetime: threshold to set
 Out: 0 on success, 1 on failure
*************************************************************************/
int
mdl_set_lifetime_warning_threshold(struct mdlparse_vars *parse_state,
                                   long long lifetime)
{
  if (lifetime < 0)
  {
    mdlerror(parse_state, "Molecule lifetimes are measured in iterations and cannot be negative");
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
int
mdl_set_missed_reaction_warning_threshold(struct mdlparse_vars *parse_state,
                                          double rxfrac)
{
  if (rxfrac < 0.0 || rxfrac > 1.0)
  {
    mdlerror(parse_state, "Values for fraction of reactions missed should be between 0 and 1");
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
int
mdl_set_time_step(struct mdlparse_vars *parse_state, double step)
{

  int error_code = mcell_set_time_step(parse_state->vol, step);
  if (error_code == 2) {
    mdlerror_fmt(
      parse_state, "Time step of %.15g requested; time step must be a positive value",
      step);
    return 1;
  }
  else if (error_code == 3) {
    mdlerror_fmt(
      parse_state,
      "Time step of %.15g requested, but the time step was already set to %.15g",
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
int
mdl_set_max_time_step(struct mdlparse_vars *parse_state, double step)
{
  if (step <= 0)
  {
    mdlerror_fmt(parse_state, "Maximum time step of %.15g requested; maximum time step must be a positive value", step);
    return 1;
  }

  if (parse_state->vol->time_step_max != 0)
  {
    mdlerror_fmt(parse_state, "Maximum time step of %.15g requested, but the maximum time step was already set to %.15g", step, parse_state->vol->time_step_max);
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
int
mdl_set_space_step(struct mdlparse_vars *parse_state, double step)
{
  if (step <= 0)
  {
    mdlerror_fmt(parse_state, "Space step of %.15g requested; space step must be a positive value", step);
    return 1;
  }

  if (parse_state->vol->space_step != 0)
  {
    mdlerror_fmt(parse_state,
                 "Space step of %.15g requested, but the space step was already set to %.15g",
                 step,
                 parse_state->vol->space_step * parse_state->vol->r_length_unit);
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
int
mdl_set_num_iterations(struct mdlparse_vars *parse_state,
                       long long numiters)
{
  /* If the iteration count was not overridden on the command-line... */
  if (parse_state->vol->iterations == INT_MIN)
  {
    parse_state->vol->iterations = numiters;
    if (parse_state->vol->iterations < 0)
    {
      mdlerror(parse_state, "Error: ITERATIONS value is negative");
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
int
mdl_set_num_radial_directions(struct mdlparse_vars *parse_state,
                              int numdirs)
{
  parse_state->vol->radial_directions = numdirs;
  if (parse_state->vol->radial_directions <= 0)
  {
    mdlerror(parse_state, "RADIAL_DIRECTIONS must be a positive number.");
    return 1;
  }

  parse_state->vol->num_directions = 0;
  if (parse_state->vol->d_step != NULL)
    free(parse_state->vol->d_step);
  if ((parse_state->vol->d_step = init_d_step(parse_state->vol->radial_directions, &parse_state->vol->num_directions)) == NULL)
  {
    mcell_allocfailed("Failed to allocate the diffusion directions table.");
    return 1;
  }

  /* Mask must contain at least every direction */
  /* Num directions, rounded up to the nearest 2^n - 1 */
  parse_state->vol->directions_mask = parse_state->vol->num_directions;
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 1);
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 2);
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 4);
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 8);
  parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 16);
  if (parse_state->vol->directions_mask > (1<<18))
  {
    mdlerror(parse_state, "Too many RADIAL_DIRECTIONS requested (max 131072).");
    return 1;
  }

  no_printf("desired radial directions = %d\n",parse_state->vol->radial_directions);
  no_printf("actual radial directions = %d\n",parse_state->vol->num_directions);
  return 0;
}

/*************************************************************************
 mdl_set_num_radial_subdivisions:
    Set the number of radial subdivisions.

 In:  parse_state: parser state
      numdivs: number of radial subdivisions to use
 Out: 0 on success, 1 on failure
*************************************************************************/
int
mdl_set_num_radial_subdivisions(struct mdlparse_vars *parse_state,
                                int numdivs)
{
  parse_state->vol->radial_subdivisions = numdivs;
  if (parse_state->vol->radial_subdivisions <= 0)
  {
    mdlerror(parse_state, "RADIAL_SUBDIVISIONS must be a positive number.");
    return 1;
  }

  /* 'x & (x-1)' clears the lowest-order bit of x.  thus, only for 0 or a  */
  /* power of 2 will the expression be zero.  we've eliminated the case of */
  /* zero in the previous test.  (this condition is required so that we    */
  /* use 'x & (RSD - 1)' as an optimized form of 'x % RSD'.)               */
  if ((parse_state->vol->radial_subdivisions & (parse_state->vol->radial_subdivisions - 1)) != 0)
  {
    mdlerror(parse_state, "Radial subdivisions must be a power of two");
    return 1;
  }

  if (parse_state->vol->r_step!=NULL) free(parse_state->vol->r_step);
  if (parse_state->vol->r_step_surface!=NULL) free(parse_state->vol->r_step_surface);
  if (parse_state->vol->r_step_release!=NULL) free(parse_state->vol->r_step_release);

  parse_state->vol->r_step = init_r_step(parse_state->vol->radial_subdivisions);
  parse_state->vol->r_step_surface = init_r_step_surface(parse_state->vol->radial_subdivisions);
  parse_state->vol->r_step_release = init_r_step_3d_release(parse_state->vol->radial_subdivisions);

  if (parse_state->vol->r_step==NULL || parse_state->vol->r_step_surface==NULL || parse_state->vol->r_step_release==NULL)
  {
    mcell_allocfailed("Failed to allocate the diffusion radial subdivisions table.");
    return 1;
  }

  no_printf("radial subdivisions = %d\n",parse_state->vol->radial_subdivisions);
  return 0;
}

/*************************************************************************
 mdl_set_interaction_radius:
    Set the interaction radius.

 In:  parse_state: parser state
      interaction_radius
 Out: 0 on success, 1 on failure
*************************************************************************/
int
mdl_set_interaction_radius(struct mdlparse_vars *parse_state,
                           double interaction_radius)
{
  parse_state->vol->rx_radius_3d = interaction_radius;
  if (parse_state->vol->rx_radius_3d <= 0)
  {
    mdlerror(parse_state, "INTERACTION_RADIUS must be a positive number.");
    return 1;
  }

  no_printf("interaction radius = %f\n",parse_state->vol->rx_radius_3d);
  return 0;
}

/*************************************************************************
 mdl_set_grid_density:
    Set the effector grid density.

 In:  parse_state: parser state
      density: global effector grid density for simulation
 Out: 0 on success, 1 on failure
*************************************************************************/
int
mdl_set_grid_density(struct mdlparse_vars *parse_state, double density)
{
  if (density <= 0)
  {
    mdlerror_fmt(parse_state, "EFFECTOR_GRID_DENSITY must be greater than 0.0 (value provided was %lg)", density);
    return 1;
  }
  parse_state->vol->grid_density = density;
  no_printf("Max density = %f\n", parse_state->vol->grid_density);

  parse_state->vol->space_step *= parse_state->vol->length_unit;
  parse_state->vol->r_length_unit = sqrt(parse_state->vol->grid_density);
  parse_state->vol->length_unit = 1.0/parse_state->vol->r_length_unit;
  parse_state->vol->space_step *= parse_state->vol->r_length_unit;

  no_printf("Length unit = %f\n", parse_state->vol->length_unit);
  return 0;
}

/*************************************************************************
 mdl_set_complex_placement_attempts:
    Set the number of times to place any particular macromolecule.

 In:  parse_state: parser state
      attempts: number of attempts
 Out: 0 on success, 1 on failure
*************************************************************************/
int
mdl_set_complex_placement_attempts(struct mdlparse_vars *parse_state, double attempts)
{
  if (attempts < 1.0 || attempts > (double) INT_MAX)
  {
    mdlerror_fmt(parse_state,
                 "COMPLEX_PLACEMENT_ATTEMPTS must be an integer between 1 and %d (value provided was %ld)",
                 INT_MAX,
                 (long int) attempts);
    return 1;
  }

  parse_state->vol->complex_placement_attempts = (int) attempts;
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
static int
schedule_async_checkpoint(struct mdlparse_vars *parse_state, unsigned int dur)
{
#ifdef _WIN32
  set_alarm_handler(&chkpt_signal_handler);
#else
  struct sigaction sa, saPrev;
  sa.sa_sigaction = NULL;
  sa.sa_handler = &chkpt_signal_handler;
  sa.sa_flags = SA_RESTART;
  sigfillset(&sa.sa_mask);

  if (sigaction(SIGALRM, &sa, &saPrev) != 0)
  {
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
int
mdl_set_realtime_checkpoint(struct mdlparse_vars *parse_state,
                            long duration,
                            int cont_after_cp)
{
  time_t t;
  time(&t);

  if (duration <= 0)
  {
    mdlerror_fmt(parse_state,
                 "Realtime checkpoint requested at %ld seconds, but the checkpoint interval must be positive.",
                 duration);
    return 1;
  }
  else if (parse_state->vol->checkpoint_alarm_time != 0)
  {
    mdlerror_fmt(parse_state,
                 "Realtime checkpoint requested at %ld seconds, but a checkpoint is already scheduled for %u seconds.",
                 duration,
                 parse_state->vol->checkpoint_alarm_time);
    return 1;
  }

  parse_state->vol->checkpoint_alarm_time = (unsigned int) duration;
  parse_state->vol->continue_after_checkpoint = cont_after_cp;
  parse_state->vol->chkpt_flag = 1;
  if (t - parse_state->vol->begin_timestamp > 0)
  {
    if (duration <= t - parse_state->vol->begin_timestamp)
    {
      mdlerror_fmt(parse_state, "Checkpoint scheduled for %ld seconds, but %ld seconds have already elapsed during parsing.  Exiting.", duration, (long)(t - parse_state->vol->begin_timestamp));
      return 1;
    }

    duration -= (t - parse_state->vol->begin_timestamp);
  }

  if (schedule_async_checkpoint(parse_state, (unsigned int) duration))
  {
    mdlerror_fmt(parse_state, "Failed to schedule checkpoint for %ld seconds", duration);
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
int
mdl_set_checkpoint_infile(struct mdlparse_vars *parse_state, char *name)
{
  FILE *file;
  if (parse_state->vol->chkpt_infile == NULL)
  {
    parse_state->vol->chkpt_infile = name;
    if ((file = fopen(parse_state->vol->chkpt_infile, "r")) == NULL)
    {
      parse_state->vol->chkpt_init = 1;
    }
    else
    {
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
int
mdl_set_checkpoint_outfile(struct mdlparse_vars *parse_state, char *name)
{
  if (parse_state->vol->chkpt_outfile == NULL)
  {
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
int
mdl_set_checkpoint_interval(struct mdlparse_vars *parse_state, long long iters)
{
  parse_state->vol->chkpt_iterations = iters;
  if (parse_state->vol->chkpt_iterations <= 0)
  {
    mdlerror(parse_state, "Error: CHECKPOINT_ITERATIONS must be a positive integer");
    return 1;
  }
  parse_state->vol->chkpt_flag = 1;
  return 0;
}

/*************************************************************************
 mdl_set_partition:
    Set the partitioning in a particular dimension.

 In:  parse_state: parser state
      dim: the dimension whose partitions we'll set
      head: the partitioning
 Out: 0 on success, 1 on failure
*************************************************************************/
int
mdl_set_partition(struct mdlparse_vars *parse_state,
                  int dim,
                  struct num_expr_list_head *head)
{
  unsigned int num_values = 0;
  double *dblp;
  struct num_expr_list *nel;

  /* Allocate array for partitions */
  dblp = CHECKED_MALLOC_ARRAY(double, (head->value_count + 2), "volume partitions");
  if (dblp == NULL)
    return 1;

  /* Copy partitions in sorted order to the array */
  dblp[num_values ++] = -GIGANTIC;
  for (nel = head->value_head; nel != NULL; nel = nel->next)
    dblp[num_values ++] = nel->value * parse_state->vol->r_length_unit;
  dblp[num_values ++] =  GIGANTIC;
  qsort(dblp, num_values, sizeof(double), & double_cmp);

  /* Copy the partitions into the model */
  switch (dim)
  {
    case X_PARTS:
      if (parse_state->vol->x_partitions != NULL)
        free(parse_state->vol->x_partitions);
      parse_state->vol->nx_parts = num_values;
      parse_state->vol->x_partitions = dblp;
      break;

    case Y_PARTS:
      if (parse_state->vol->y_partitions != NULL)
        free(parse_state->vol->y_partitions);
      parse_state->vol->ny_parts = num_values;
      parse_state->vol->y_partitions = dblp;
      break;

    case Z_PARTS:
      if (parse_state->vol->z_partitions != NULL)
        free(parse_state->vol->z_partitions);
      parse_state->vol->nz_parts = num_values;
      parse_state->vol->z_partitions = dblp;
      break;

    default: UNHANDLED_CASE(dim);
  }

  if (! head->shared)
    mdl_free_numeric_list(head->value_head);

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
static struct object *
mdl_make_new_object(struct mdlparse_vars *parse_state, char *obj_name)
{
  if ((retrieve_sym(obj_name, parse_state->vol->obj_sym_table)) != NULL)
  {
    mdlerror_fmt(parse_state,"Object '%s' is already defined", obj_name);
    return NULL;
  }

  struct object *obj_ptr = make_new_object(parse_state->vol, obj_name);

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
struct sym_table *
mdl_start_object(struct mdlparse_vars *parse_state, char *name)
{
  // Create new fully qualified name.
  char *new_name;
  struct object_creation obj_creation;
  obj_creation.object_name_list = parse_state->object_name_list;
  obj_creation.object_name_list_end = parse_state->object_name_list_end;
  if ((new_name = push_object_name(&obj_creation, name)) == NULL)
  {
    mdlerror_fmt(parse_state, "Out of memory while creating object: %s", name);
    free(name);
    return NULL;
  }
  parse_state->object_name_list = obj_creation.object_name_list;
  parse_state->object_name_list_end = obj_creation.object_name_list_end;

  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = mdl_make_new_object(parse_state, new_name);
  if (obj_ptr == NULL)
  {
    free(name);
    free(new_name);
    return NULL;
  }

  struct sym_table *sym_ptr = obj_ptr->sym;
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
void
mdl_finish_object(struct mdlparse_vars *parse_state)
{
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
void
mdl_object_list_singleton(struct object_list *head, struct object *objp)
{
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
void
mdl_add_object_to_list(struct object_list *head, struct object *objp)
{
  objp->next = NULL;
  head->obj_tail = head->obj_tail->next = objp;
}



/*************************************************************************
 SYMBOL_TYPE_DESCRIPTIONS:
    Human-readable symbol type descriptions, indexed by symbol type (MOL, OBJ,
    STR, etc.)
*************************************************************************/
static const char *SYMBOL_TYPE_DESCRIPTIONS[] =
{
  "reaction",
  "reaction pathway",
  "molecule",
  "object",
  "release pattern",
  "region",
  "numeric variable",
  "string variable",
  "array variable",
  "file stream",
  "placeholder value",
  "trigger",
};

/*************************************************************************
 SYMBOL_TYPE_ARTICLES:
    Indefinite articles to go with human-readable symbol type descriptions,
    indexed by symbol type (MOL, OBJ, STR, etc.)
*************************************************************************/
static const char *SYMBOL_TYPE_ARTICLES[] =
{
  "a",
  "a",
  "a",
  "an",
  "a",
  "a",
  "a",
  "a",
  "an",
  "a",
  "a",
  "a",
};

/*************************************************************************
 mdl_symbol_type_name:
    Get a human-readable description of a symbol type, given its numerical
    code.

 In:  type: numeric type code
 Out: the symbol type name
*************************************************************************/
static char const *
mdl_symbol_type_name(enum symbol_type_t type)
{

  if (type <= 0 || type >= (int) COUNT_OF(SYMBOL_TYPE_DESCRIPTIONS))
  {
    mcell_internal_error("Invalid symbol type '%d'", type);
    return "(unknown symbol type)";
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
static char const *
mdl_symbol_type_name_article(enum symbol_type_t type)
{

  if (type <= 0 || type >= (int) COUNT_OF(SYMBOL_TYPE_ARTICLES))
  {
    mcell_internal_error("Invalid symbol type '%d'", type);
    return "an";
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
static struct sym_table *
mdl_existing_symbol(struct mdlparse_vars *parse_state,
                    char *name,
                    struct sym_table_head *tab,
                    int type)
{
  struct sym_table *symp = retrieve_sym(name, tab);
  if (symp == NULL)
    mdlerror_fmt(parse_state, "Undefined %s: %s", mdl_symbol_type_name(type), name);
  else if (symp->sym_type != type)
  {
    mdlerror_fmt(parse_state, "Invalid type for symbol %s: expected %s, but found %s",
                 name,
                 mdl_symbol_type_name(type),
                 mdl_symbol_type_name(symp->sym_type));
    symp = NULL;
  }
  else if (strcmp(symp->name, "GENERIC_MOLECULE") == 0)
  {
    mdlerror_fmt(parse_state, "The keyword 'GENERIC_MOLECULE' is obsolete. Please use instead 'ALL_VOLUME_MOLECULES', 'ALL_SURFACE_MOLECULES' or 'ALL_MOLECULES'.");
    symp = NULL;
  }
  else
  {
#ifdef KELP
    ++ symp->ref_count;
    no_printf("ref_count: %d\n", symp->ref_count);
#endif
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
static struct sym_table *
mdl_existing_symbol_2types(struct mdlparse_vars *parse_state,
                           char *name,
                           struct sym_table_head *tab1,
                           int type1,
                           struct sym_table_head *tab2,
                           int type2)
{
  struct sym_table *symp;
  symp = retrieve_sym(name, tab1);
  if (symp == NULL)
  {
    symp = retrieve_sym(name, tab2);
    if (symp == NULL)
      mdlerror_fmt(parse_state, "Undefined %s or %s: %s",
                   mdl_symbol_type_name(type1),
                   mdl_symbol_type_name(type2),
                   name);
  }
  else
  {
    if (retrieve_sym(name, tab2) != NULL)
    {
      mdlerror_fmt(parse_state, "Named object '%s' could refer to %s %s or %s %s.  Please rename one of them.",
                   name,
                   mdl_symbol_type_name_article(type1),
                   mdl_symbol_type_name(type1),
                   mdl_symbol_type_name_article(type2),
                   mdl_symbol_type_name(type2));
      symp = NULL;
    }
  }

  if (symp)
  {
#ifdef KELP
    ++ symp->ref_count;
    no_printf("ref_count: %d\n", symp->ref_count);
#endif
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
                             char const *wildcard,
                             struct sym_table_head *tab,
                             int type)
{
  struct sym_table_list *symbols = NULL, *stl;
  for (int i = 0; i < tab->n_bins; i++)
  {
    for (struct sym_table *sym_t = tab->entries[i];
         sym_t != NULL;
         sym_t = sym_t->next)
    {
      if (sym_t->sym_type != type) continue;

      if (is_wildcard_match((char *) wildcard, sym_t->name))
      {
        stl = CHECKED_MEM_GET(parse_state->sym_list_mem, "wildcard symbol list");
        if (stl == NULL)
        {
          if (symbols) mem_put_list(parse_state->sym_list_mem, symbols);
          return NULL;
        }

        stl->node = sym_t;
        stl->next = symbols;
        symbols = stl;
      }
    }
  }

  if (symbols == NULL)
  {
    switch (type)
    {
      case OBJ: mdlerror_fmt(parse_state, "No objects found matching wildcard \"%s\"", wildcard); break;
      case MOL: mdlerror_fmt(parse_state, "No molecules found matching wildcard \"%s\"", wildcard); break;
      default:  mdlerror_fmt(parse_state, "No items found matching wildcard \"%s\"", wildcard); break;
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
static int
compare_sym_names(void *a,void *b)
{
  return strcmp(((struct sym_table*)a)->name , ((struct sym_table*)b)->name) <= 0;
}

/*************************************************************************
 sort_sym_list_by_name:
    Sort a symbol list in collation order by name.

 In:  unsorted: unsorted list
 Out: a sorted list of symbols
*************************************************************************/
static struct sym_table_list*
sort_sym_list_by_name(struct sym_table_list *unsorted)
{
  return (struct sym_table_list*)void_list_sort_by((struct void_list*)unsorted, compare_sym_names);
}

/*************************************************************************
 mdl_existing_object:
    Find an existing object.  Print an error message if the object isn't found.

 In:  parse_state: parser state
      name: fully qualified object name
 Out: the object, or NULL if not found
*************************************************************************/
struct sym_table *
mdl_existing_object(struct mdlparse_vars *parse_state, char *name)
{
  return mdl_existing_symbol(parse_state, name, parse_state->vol->obj_sym_table, OBJ);
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
                              char *wildcard)
{
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
struct sym_table *
mdl_existing_region(struct mdlparse_vars *parse_state,
                    struct sym_table *obj_symp,
                    char *name)
{
  struct sym_table *symp;
  char *region_name = CHECKED_SPRINTF("%s,%s",
                                      obj_symp->name,
                                      name);
  if (region_name == NULL)
  {
    free(name);
    return NULL;
  }

  symp = mdl_existing_symbol(parse_state, region_name, parse_state->vol->reg_sym_table, REG);
  free(name);
  return symp;
}

/*************************************************************************
 mdl_existing_molecule:
    Find an existing molecule species.  Print an error message if it isn't
    found.

 In:  parse_state: parser state
      name: species name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_table *
mdl_existing_molecule(struct mdlparse_vars *parse_state, char *name)
{
  return mdl_existing_symbol(parse_state, name, parse_state->vol->mol_sym_table, MOL);
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
                          struct sym_table *sym)
{
  struct sym_table_list *stl = CHECKED_MEM_GET(parse_state->sym_list_mem, "symbol list item");
  if (stl != NULL)
  {
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
mdl_existing_molecule_list(struct mdlparse_vars *parse_state,
                           char *name)
{
  struct sym_table *symp = mdl_existing_molecule(parse_state, name);
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
                                char *wildcard)
{
  struct sym_table_list *stl;
  char *wildcard_string;
  if (! (wildcard_string = mdl_strip_quotes(wildcard)))
    return NULL;

  stl = mdl_find_symbols_by_wildcard(parse_state, wildcard_string, parse_state->vol->mol_sym_table, MOL);
  if (stl == NULL)
  {
    free(wildcard_string);
    return NULL;
  }
  free(wildcard_string);

  return sort_sym_list_by_name(stl);
}

/*************************************************************************
 mdl_existing_macromolecule:
    Find an existing macromolecule species.  Print an error message if it isn't
    found, or isn't a macromolecule.

 In:  parse_state: parser state
      name: species name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_table *
mdl_existing_macromolecule(struct mdlparse_vars *parse_state,
                           char *name)
{
  struct sym_table *symp = mdl_existing_molecule(parse_state, name);
  if (symp == NULL)
    return NULL;

  struct species *sp = (struct species *) symp->value;
  if (! (sp->flags & IS_COMPLEX))
  {
    mdlerror_fmt(parse_state, "Molecule '%s' is not a macromolecule", symp->name);
    return NULL;
  }

  return symp;
}

/*************************************************************************
 mdl_existing_surface_molecule:
    Find an existing surface molecule species.  Print an error message if it
    isn't found, or isn't a surface molecule.

 In:  parse_state: parser state
      name: species name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_table *
mdl_existing_surface_molecule(struct mdlparse_vars *parse_state,
                              char *name)
{
  struct sym_table *symp = mdl_existing_molecule(parse_state, name);
  if (symp == NULL)
    return NULL;

  struct species *sp = (struct species *) symp->value;
  if (! (sp->flags & ON_GRID))
  {
    mdlerror_fmt(parse_state, "Molecule '%s' is not a surface molecule", symp->name);
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
struct sym_table *
mdl_existing_surface_class(struct mdlparse_vars *parse_state,
                           char *name)
{
  struct sym_table *symp = mdl_existing_molecule(parse_state, name);
  if (symp == NULL)
    return NULL;

  struct species *sp = (struct species *) symp->value;
  if (! (sp->flags & IS_SURFACE))
  {
    mdlerror_fmt(parse_state, "Species '%s' is not a surface class", sp->sym->name);
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
struct sym_table *
mdl_existing_variable(struct mdlparse_vars *parse_state, char *name)
{
  struct sym_table *st = NULL;

  /* Attempt to fetch existing variable */
  if ((st=retrieve_sym(name, parse_state->vol->var_sym_table)) != NULL)
  {
    free(name);
#ifdef KELP
    st->ref_count++;
    no_printf("ref_count: %d\n",st->ref_count);
#endif
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
struct sym_table *
mdl_existing_array(struct mdlparse_vars *parse_state, char *name)
{
  return mdl_existing_symbol(parse_state, name, parse_state->vol->var_sym_table, ARRAY);
}

/**************************************************************************
 mdl_existing_double:
    Find a named numeric variable if it exists.  Print an error message if it
    isn't found.

 In: parse_state: the parse variables structure
     name: the name of the variable
 Out: a symbol table entry, or NULL if we couldn't find the variable
**************************************************************************/
struct sym_table *mdl_existing_double(struct mdlparse_vars *parse_state, char *name)
{
  return mdl_existing_symbol(parse_state, name, parse_state->vol->var_sym_table, DBL);
}

/**************************************************************************
 mdl_existing_string:
    Find a named string variable if it exists.  Print an error message if it
    isn't found.

 In: parse_state: the parse variables structure
     name: the name of the variable
 Out: a symbol table entry, or NULL if we couldn't find the variable
**************************************************************************/
struct sym_table *
mdl_existing_string(struct mdlparse_vars *parse_state, char *name)
{
  return mdl_existing_symbol(parse_state, name, parse_state->vol->var_sym_table, STR);
}

/**************************************************************************
 mdl_existing_num_or_array:
    Find a named numeric or array variable if it exists.  Print an error message
    if it isn't found.

 In: parse_state: the parse variables structure
     name: the name of the variable
 Out: a symbol table entry, or NULL if we couldn't find the variable
**************************************************************************/
struct sym_table *
mdl_existing_num_or_array(struct mdlparse_vars *parse_state,
                          char *name)
{
  struct sym_table *st = NULL;

  /* Attempt to fetch existing variable */
  if ((st = retrieve_sym(name, parse_state->vol->var_sym_table)) != NULL)
  {
    if (st->sym_type == STR)
    {
      mdlerror_fmt(parse_state, "Incorrect type (got string, expected number or array): %s", name);
      return NULL;
    }

#ifdef KELP
    st->ref_count++;
    no_printf("ref_count: %d\n",st->ref_count);
#endif
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
struct sym_table *
mdl_existing_rxn_pathname_or_molecule(struct mdlparse_vars *parse_state,
                                      char *name)
{
  return mdl_existing_symbol_2types(parse_state,
                                    name,
                                    parse_state->vol->rxpn_sym_table,
                                    RXPN,
                                    parse_state->vol->mol_sym_table,
                                    MOL);
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
struct sym_table *
mdl_existing_release_pattern_or_rxn_pathname(struct mdlparse_vars *parse_state,
                                             char *name)
{
  return mdl_existing_symbol_2types(parse_state,
                                    name,
                                    parse_state->vol->rpat_sym_table,
                                    RPAT,
                                    parse_state->vol->rxpn_sym_table,
                                    RXPN);
}

/*************************************************************************
 mdl_existing_molecule_or_object:
    Find an existing molecule or object.  Print an error message if it isn't
    found, or if the name could refer to either type of object.  This is only
    used for the NAME_LIST in old-style VIZ_DATA_OUTPUT blocks.

 In:  parse_state: parser state
      name: symbol name
 Out: the symbol, or NULL if not found
*************************************************************************/
struct sym_table *
mdl_existing_molecule_or_object(struct mdlparse_vars *parse_state, char *name)
{
  return mdl_existing_symbol_2types(parse_state,
                                    name,
                                    parse_state->vol->mol_sym_table,
                                    MOL,
                                    parse_state->vol->obj_sym_table,
                                    OBJ);
}

/*************************************************************************
 mdl_existing_file_stream:
    Find an existing file stream.  Print an error message if the stream isn't
    found.

 In:  parse_state: parser state
      name: stream name
 Out: the stream, or NULL if not found
*************************************************************************/
struct sym_table *
mdl_existing_file_stream(struct mdlparse_vars *parse_state,
                         char *name)
{
  struct sym_table *sym = mdl_existing_symbol(parse_state, name, parse_state->vol->fstream_sym_table, FSTRM);
  if (sym == NULL)
    return sym;

  struct file_stream *filep = (struct file_stream *) sym->value;
  if (filep->stream == NULL)
  {
    mdlerror_fmt(parse_state, "File stream '%s' has already been closed", sym->name);
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
struct sym_table_list *
mdl_meshes_by_wildcard(struct mdlparse_vars *parse_state, char *wildcard)
{
  /* Scan for objects matching the wildcard */
  struct sym_table_list *matches = mdl_find_symbols_by_wildcard(parse_state,
                                                                wildcard,
                                                                parse_state->vol->obj_sym_table,
                                                                OBJ);
  if (matches == NULL)
  {
    free(wildcard);
    return NULL;
  }

  /* Scan through the list, discarding inappropriate objects */
  struct sym_table_list *cur_match, *next_match, **prev = &matches;
  for (cur_match = matches; cur_match != NULL; cur_match = next_match)
  {
    next_match = cur_match->next;

    struct object *objp = (struct object *) cur_match->node->value;
    if (objp->object_type != POLY_OBJ && objp->object_type != BOX_OBJ)
    {
      *prev = cur_match->next;
      mem_put(parse_state->sym_list_mem, cur_match);
    }
    else
      prev = &cur_match->next;
  }

  if (matches == NULL)
  {
    mdlerror_fmt(parse_state, "No matches for the wildcard '%s' were mesh objects", wildcard);
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
int
mdl_transform_rotate(struct mdlparse_vars *parse_state,
                     double (*mat)[4],
                     struct vector3 *axis,
                     double angle)
{
  if (transform_rotate(mat, axis, angle))
  {
    mdlerror(parse_state, "Rotation vector has zero length.");
    return 1;
  }
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
static struct region *
mdl_make_new_region(struct mdlparse_vars *parse_state,
                    char *obj_name,
                    char *region_last_name)
{
  struct sym_table *gp;
  char *region_name;

  region_name = CHECKED_SPRINTF("%s,%s",
                                obj_name,
                                region_last_name);
  if (region_name == NULL)
    return NULL;

  if ((retrieve_sym(region_name, parse_state->vol->reg_sym_table)) != NULL)
  {
    mdlerror_fmt(parse_state, "Region already defined: %s", region_name);
    free(region_name);
    return NULL;
  }

  if ((gp = store_sym(region_name, REG, parse_state->vol->reg_sym_table, NULL)) == NULL)
  {
    mcell_allocfailed("Failed to store a region in the region symbol table.");
    free(region_name);
    return NULL;
  }

  free(region_name);
  return (struct region *) gp->value;
}

/*************************************************************************
 mdl_copy_object_regions:
    Duplicate src_obj's regions on dst_obj.

 In:  parse_state: parser state
      dst_obj: destination object
      src_obj: object from which to copy
 Out: 0 on success, 1 on failure
*************************************************************************/
static int
mdl_copy_object_regions(struct mdlparse_vars *parse_state,
                        struct object *dst_obj,
                        struct object *src_obj)
{
  struct region_list *src_rlp;
  struct region *dst_reg, *src_reg;
  struct eff_dat *dst_eff, *src_eff;

  /* Copy each region */
  for (src_rlp = src_obj->regions;
       src_rlp != NULL;
       src_rlp = src_rlp->next)
  {
    src_reg = src_rlp->reg;
    if ((dst_reg = mdl_create_region(parse_state, dst_obj, src_reg->region_last_name)) == NULL)
      return 1;

    /* Copy over simple region attributes */
    dst_reg->surf_class       = src_reg->surf_class;
    dst_reg->flags            = src_reg->flags;
    dst_reg->area             = src_reg->area;
    dst_reg->bbox             = src_reg->bbox;
    dst_reg->manifold_flag    = src_reg->manifold_flag;
    dst_reg->region_viz_value = src_reg->region_viz_value;

    /* Copy membership data */
    if (src_reg->membership != NULL)
    {
      dst_reg->membership = duplicate_bit_array(src_reg->membership);
      if (dst_reg->membership==NULL)
      {
        mcell_allocfailed("Failed allocation for membership array in %s",
                     dst_obj->sym->name);
        return 1;
      }
    }
    else mdl_warning(parse_state, "No membership data for %s\n", src_reg->sym->name);

    /* Copy effector data list */
    struct eff_dat *effdp_tail = NULL;
    for (src_eff = src_reg->eff_dat_head;
         src_eff != NULL;
         src_eff = src_eff->next)
    {
      dst_eff = CHECKED_MALLOC_STRUCT(struct eff_dat, "surface molecule info");
      if (dst_eff == NULL)
        return 1;

      if (effdp_tail != NULL)
        effdp_tail->next = dst_eff;
      else
        dst_reg->eff_dat_head = dst_eff;
      effdp_tail = dst_eff;

      dst_eff->eff           = src_eff->eff;
      dst_eff->quantity_type = src_eff->quantity_type;
      dst_eff->quantity      = src_eff->quantity;
      dst_eff->orientation   = src_eff->orientation;
      dst_eff->next = NULL;
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
static struct region*
find_corresponding_region(struct region *old_r,
                          struct object *old_ob,
                          struct object *new_ob,
                          struct object *instance,
                          struct sym_table_head *symhash)
{
  struct object *ancestor;
  struct object *ob;
  struct sym_table *gp;

  ancestor = common_ancestor(old_ob,old_r->parent);

  if (ancestor==NULL)
  {
    for (ob=old_r->parent ; ob!=NULL ; ob=ob->parent)
    {
      if (ob==instance) break;
    }

    if (ob==NULL) return NULL;
    else return old_r;  /* Point to same already-instanced object */
  }
  else
  {
    /* Length of old prefix and name */
    const size_t old_prefix_idx = strlen(ancestor->sym->name);
    const size_t old_name_len   = strlen(old_ob->sym->name);

    /* Suffix is everything in the old name after the old prefix. */
    const size_t suffix_len     = old_name_len - old_prefix_idx;

    /* New prefix length is new name, less suffix length. */
    const size_t new_name_len   = strlen(new_ob->sym->name);
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
    char new_name[ max_name_len ];
    strncpy(new_name, new_ob->sym->name, new_prefix_idx);
    strncpy(new_name + new_prefix_idx,            /* new prefix */
            old_r->sym->name + old_prefix_idx,    /* old suffix + region name */
            max_name_len - new_prefix_idx);

    /* Finally, retrieve symbol from newly-constructed name. */
    gp = retrieve_sym(new_name, symhash);
    if (gp == NULL) return NULL;
    else return (struct region*) gp->value;
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
static struct release_evaluator*
duplicate_rel_region_expr(struct mdlparse_vars *parse_state,
                          struct release_evaluator *expr,
                          struct object *old_self,
                          struct object *new_self,
                          struct object *instance)
{
  struct region *r;
  struct release_evaluator *nexp;

  nexp = CHECKED_MALLOC_STRUCT(struct release_evaluator, "region release expression");
  if (nexp == NULL)
    return NULL;

  nexp->op = expr->op;

  if (expr->left!=NULL)
  {
    if (expr->op&REXP_LEFT_REGION)
    {
      r = find_corresponding_region(expr->left,old_self,new_self,instance,parse_state->vol->reg_sym_table);

      if (r==NULL)
      {
        mdlerror_fmt(parse_state,
                     "Can't find new region corresponding to %s for %s (copy of %s)",
                     ((struct region *) expr->left)->sym->name,
                     new_self->sym->name,
                     old_self->sym->name);
        return NULL;
      }

      nexp->left = r;
    }
    else nexp->left = duplicate_rel_region_expr(parse_state, expr->left, old_self, new_self, instance);
  }
  else nexp->left = NULL;

  if (expr->right!=NULL)
  {
    if (expr->op&REXP_RIGHT_REGION)
    {
      r = find_corresponding_region(expr->right, old_self, new_self, instance, parse_state->vol->reg_sym_table);

      if (r==NULL)
      {
        mdlerror_fmt(parse_state,
                     "Can't find new region corresponding to %s for %s (copy of %s)",
                     ((struct region *) expr->right)->sym->name,
                     new_self->sym->name,
                     old_self->sym->name);
        return NULL;
      }

      nexp->right = r;
    }
    else nexp->right = duplicate_rel_region_expr(parse_state, expr->right, old_self, new_self, instance);
  }
  else nexp->right = NULL;

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
static struct release_site_obj*
duplicate_release_site(struct mdlparse_vars *parse_state,
                       struct release_site_obj *old,
                       struct object *new_self,
                       struct object *instance)
{
  struct release_site_obj *rel_site_obj = CHECKED_MALLOC_STRUCT(
    struct release_site_obj, "release site");
  if (rel_site_obj == NULL) {
    return NULL;
  }

  if (old->location != NULL)
  {
    if ((rel_site_obj->location = CHECKED_MALLOC_STRUCT(
      struct vector3, "release site location")) == NULL)
    {
      return NULL;
    }
    *(rel_site_obj->location) = *(old->location);
  }
  else {
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
  rel_site_obj->name = NULL;

  if (old->region_data != NULL)
  {
    struct release_region_data *rel_reg_data = CHECKED_MALLOC_STRUCT(
      struct release_region_data, "release region data");
    if (rel_reg_data == NULL) {
      return NULL;
    }

    memcpy(&(rel_reg_data->llf), &(old->region_data->llf), sizeof(struct vector3));
    memcpy(&(rel_reg_data->urb), &(old->region_data->urb), sizeof(struct vector3));
    rel_reg_data->n_walls_included = -1;
    rel_reg_data->cum_area_list = NULL;
    rel_reg_data->wall_index = NULL;
    rel_reg_data->obj_index = NULL;
    rel_reg_data->n_objects = -1;
    rel_reg_data->owners = NULL;
    rel_reg_data->in_release = NULL;
    rel_reg_data->self = new_self;

    rel_reg_data->expression = duplicate_rel_region_expr(
      parse_state,
      old->region_data->expression,
      old->region_data->self,
      new_self,
      instance);
    if (rel_reg_data->expression == NULL) {
      return NULL;
    }

    rel_site_obj->region_data = rel_reg_data;
  }
  else {
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
int
mdl_deep_copy_object(struct mdlparse_vars *parse_state,
                     struct object *dst_obj,
                     struct object *src_obj)
{
  struct object *src_child;

  /* Copy over simple object attributes */
  dst_obj->object_type    = src_obj->object_type;
  dst_obj->n_walls        = src_obj->n_walls;
  dst_obj->n_walls_actual = src_obj->n_walls_actual;
  dst_obj->walls          = src_obj->walls;
  dst_obj->wall_p         = src_obj->wall_p;
  dst_obj->n_verts        = src_obj->n_verts;
  dst_obj->vertices       = src_obj->vertices;

  /* Copy over regions */
  if (mdl_copy_object_regions(parse_state, dst_obj, src_obj))
    return 1;

  /* Inherit object coordinate transformations */
  mult_matrix(dst_obj->t_matrix, src_obj->t_matrix, dst_obj->t_matrix, 4, 4, 4);

  switch (dst_obj->object_type)
  {
    case META_OBJ:
      /* Copy children */
      for (src_child = src_obj->first_child;
           src_child !=NULL;
           src_child = src_child->next)
      {
        struct object *dst_child;
        char *child_obj_name = CHECKED_SPRINTF("%s.%s",
                                               dst_obj->sym->name,
                                               src_child->last_name);
        if (child_obj_name == NULL)
          return 1;

        /* Create child object */
        if ((dst_child = mdl_make_new_object(parse_state, child_obj_name)) == NULL)
        {
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
      dst_obj->contents = duplicate_release_site(parse_state,
                                                 src_obj->contents,
                                                 dst_obj,
                                                 parse_state->vol->root_instance);
      if (dst_obj->contents == NULL) return 1;
      struct release_site_obj *rso = (struct release_site_obj *) dst_obj->contents;
      rso->name = mdl_strdup(dst_obj->sym->name);
      if (rso->name == NULL)
        return 1;
      break;

    case BOX_OBJ:
    case POLY_OBJ:
      dst_obj->contents = src_obj->contents;
      break;

    default:
      mdlerror_fmt(parse_state, "Error: bad object type %d", dst_obj->object_type);
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
static struct subdivided_box*
init_cuboid(struct mdlparse_vars *parse_state,
            struct vector3 *p1,
            struct vector3 *p2)
{
  struct subdivided_box *b;

  if (p2->x-p1->x < EPS_C || p2->y-p1->y < EPS_C || p2->z-p1->z < EPS_C)
  {
    mdlerror(parse_state, "Box vertices out of order or box is degenerate.");
    return NULL;
  }

  b = CHECKED_MALLOC_STRUCT(struct subdivided_box, "subdivided box");
  if (b == NULL)
    return NULL;

  b->nx = b->ny = b->nz = 2;
  if ((b->x = CHECKED_MALLOC_ARRAY(double, b->nx, "subdivided box X partitions")) == NULL) return NULL;
  if ((b->y = CHECKED_MALLOC_ARRAY(double, b->ny, "subdivided box Y partitions")) == NULL) return NULL;
  if ((b->z = CHECKED_MALLOC_ARRAY(double, b->nz, "subdivided box Z partitions")) == NULL) return NULL;

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
     egd: the effector grid density, which limits how fine divisions can be
 Out: returns 1 on failure, 0 on success.  The box has additional subdivisions
      added that fall at the edges of the patch so that the patch can be
      specified in terms of subdivisions (i.e. can be constructed by triangles
      that tile each piece of the subdivided surface).
*************************************************************************/
static int
refine_cuboid(struct mdlparse_vars *parse_state,
              struct vector3 *p1,
              struct vector3 *p2,
              struct subdivided_box *b,
              double egd)
{
  int i,j,k;
  double *new_list;
  int new_n;

  i = check_patch(b,p1,p2,egd);
  if (i==0)
  {
    mdlerror(parse_state, "Could not refine box to include patch: Invalid patch specified");
    return 1;
  }

  if (i&BRANCH_X)
  {
    new_n = b->nx + 2;
    for (j=0;j<b->nx;j++)
    {
      if (p1->x == b->x[j]) new_n--;
      if (p2->x == b->x[j]) new_n--;
    }
    if (new_n > b->nx)
    {
      new_list = CHECKED_MALLOC_ARRAY(double, new_n, "refined subdivided box X partitions");
      if (new_list == NULL)
          return 1;

      for (j=k=0 ; b->x[j]<p1->x ; j++) new_list[k++]=b->x[j];
      if (b->x[j]!=p1->x) new_list[k++]=p1->x;
      for (; b->x[j]<p2->x ; j++) new_list[k++]=b->x[j];
      if (p1->x!=p2->x && b->x[j]!=p2->x) new_list[k++]=p2->x;
      for (; j<b->nx ; j++) new_list[k++]=b->x[j];

      free(b->x);
      b->x = new_list;
      b->nx = new_n;
    }
  }
  if (i&BRANCH_Y) /* Same as above with x->y */
  {
    new_n = b->ny + 2;
    for (j=0;j<b->ny;j++)
    {
      if (p1->y == b->y[j]) new_n--;
      if (p2->y == b->y[j]) new_n--;
    }
    if (new_n > b->ny)
    {
      new_list = CHECKED_MALLOC_ARRAY(double, new_n, "refined subdivided box Y partitions");
      if (new_list==NULL)
        return 1;

      for (j=k=0 ; b->y[j]<p1->y ; j++) new_list[k++]=b->y[j];
      if (b->y[j]!=p1->y) new_list[k++]=p1->y;
      for (; b->y[j]<p2->y ; j++) new_list[k++]=b->y[j];
      if (p1->y!=p2->y && b->y[j]!=p2->y) new_list[k++]=p2->y;
      for (; j<b->ny ; j++) new_list[k++]=b->y[j];

      free(b->y);
      b->y = new_list;
      b->ny = new_n;
    }
  }
  if (i&BRANCH_Z)  /* Same again, x->z */
  {
    new_n = b->nz + 2;
    for (j=0;j<b->nz;j++)
    {
      if (p1->z == b->z[j]) new_n--;
      if (p2->z == b->z[j]) new_n--;
    }
    if (new_n > b->nz)
    {
      new_list = CHECKED_MALLOC_ARRAY(double, new_n, "refined subdivided box Z partitions");
      if (new_list == NULL)
        return 1;

      for (j=k=0 ; b->z[j]<p1->z ; j++) new_list[k++]=b->z[j];
      if (b->z[j]!=p1->z) new_list[k++]=p1->z;
      for (; b->z[j]<p2->z ; j++) new_list[k++]=b->z[j];
      if (p1->z!=p2->z && b->z[j]!=p2->z) new_list[k++]=p2->z;
      for (; j<b->nz ; j++) new_list[k++]=b->z[j];

      free(b->z);
      b->z = new_list;
      b->nz = new_n;
    }
  }

  return 0;
}

#ifdef DEBUG
/*************************************************************************
 print_cuboid:
 In: b: subdivided box to print
 Out: coordinates of subdivided box are printed to stdout
*************************************************************************/
void
print_cuboid(struct subdivided_box *b)
{
  int i;
  printf("X coordinate:\n");
  for (i=0;i<b->nx;i++) printf("  %.8e\n",b->x[i]);
  printf("Y coordinate:\n");
  for (i=0;i<b->ny;i++) printf("  %.8e\n",b->y[i]);
  printf("Z coordinate:\n");
  for (i=0;i<b->nz;i++) printf("  %.8e\n",b->z[i]);
}
#endif

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
static int
divide_cuboid(struct mdlparse_vars *parse_state,
              struct subdivided_box *b,
              int axis,
              int idx,
              int ndiv)
{
  double *old_list;
  double *new_list;
  int old_n;
  int new_n;
  int i,j,k;

  if (ndiv<2) ndiv=2;
  switch(axis)
  {
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
      mdlerror_fmt(parse_state, "File '%s', Line %ld: Unknown flag is used.", __FILE__, (long)__LINE__);
      return 1;
  }

  new_n = old_n + ndiv - 1;
  new_list = CHECKED_MALLOC_ARRAY(double, new_n, "refined subdivided box partitions");
  if (new_list==NULL)
     return 1;

  for (i=j=0 ; i<=idx ; i++,j++) new_list[j] = old_list[i];
  for (k=1;k<ndiv;k++) new_list[j++] = (((double)k/(double)ndiv))*(old_list[i]-old_list[i-1]) + old_list[i-1];
  for (; i<old_n ; i++,j++) new_list[j] = old_list[i];

  switch(axis)
  {
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
      break;
  }

  free(old_list);
  return 0;
}


/*************************************************************************
 reaspect_cuboid:
 In: parse_state: parser state
     b: a subdivided box whose surface is a bunch of rectangles
     ratio: the maximum allowed aspect ratio (long side / short side) of a rectangle
 Out: returns 1 on failure, 0 on success.  The subdivided box is further
      divided to ensure that all surface subdivisions meet the aspect ratio
      criterion.
 Note: This is not, in general, possible if max_ratio is less than sqrt(2), and
       it is diffucult if max_ratio is less than 2.
*************************************************************************/
static int
reaspect_cuboid(struct mdlparse_vars *parse_state,
                struct subdivided_box *b,
                double max_ratio)
{
  double min_x,min_y,min_z,max_x,max_y,max_z;
  int jx,jy,jz;
  int i,j;
  int changed;

  do
  {
    changed = 0;

    max_x = min_x = b->x[1] - b->x[0];
    jx = 0;
    for (i=1;i<b->nx-1;i++)
    {
      if (min_x > b->x[i+1] - b->x[i])
      {
    min_x = b->x[i+1] - b->x[i];
      }
      else if (max_x < b->x[i+1] - b->x[i])
      {
    max_x = b->x[i+1] - b->x[i];
    jx = i;
      }
    }

    max_y = min_y = b->y[1] - b->y[0];
    jy = 0;
    for (i=1;i<b->ny-1;i++)
    {
      if (min_y > b->y[i+1] - b->y[i])
      {
    min_y = b->y[i+1] - b->y[i];
      }
      else if (max_y < b->y[i+1] - b->y[i])
      {
    max_y = b->y[i+1] - b->y[i];
    jy = i;
      }
    }

    max_z = min_z = b->z[1] - b->z[0];
    jz = 0;
    for (i=1;i<b->nz-1;i++)
    {
      if (min_z > b->z[i+1] - b->z[i])
      {
    min_z = b->z[i+1] - b->z[i];
      }
      else if (max_z < b->z[i+1] - b->z[i])
      {
    max_z = b->z[i+1] - b->z[i];
    jz = i;
      }
    }

    if (max_y/min_x > max_ratio)
    {
      j = divide_cuboid(parse_state, b, BRANCH_Y, jy, (int)ceil(max_y/(max_ratio*min_x)));
      if (j) return 1;
      changed |= BRANCH_Y;
    }
    else if (max_x/min_y > max_ratio)
    {
      j = divide_cuboid(parse_state, b, BRANCH_X, jx, (int)ceil(max_x/(max_ratio*min_y)));
      if (j) return 1;
      changed |= BRANCH_X;
    }

    if ((changed&BRANCH_X)==0 && max_z/min_x > max_ratio)
    {
      j = divide_cuboid(parse_state, b, BRANCH_Z, jz, (int)ceil(max_z/(max_ratio*min_x)));
      if (j) return 1;
      changed |= BRANCH_Z;
    }
    else if ((changed&BRANCH_X)==0 && max_x/min_z > max_ratio)
    {
      j = divide_cuboid(parse_state, b, BRANCH_X, jx, (int)ceil(max_x/(max_ratio*min_z)));
      if (j) return 1;
      changed |= BRANCH_X;
    }

    if ((changed&(BRANCH_Y|BRANCH_Z))==0 && max_z/min_y > max_ratio)
    {
      j = divide_cuboid(parse_state, b, BRANCH_Z, jz, (int)ceil(max_z/(max_ratio*min_y)));
      if (j) return 1;
      changed |= BRANCH_Z;
    }
    else if ((changed&(BRANCH_Y|BRANCH_Z))==0 && max_y/min_z > max_ratio)
    {
      j = divide_cuboid(parse_state, b, BRANCH_Y, jy, (int)ceil(max_y/(max_ratio*min_z)));
      if (j) return 1;
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
static int
count_cuboid_vertices(struct subdivided_box *sb)
{
  return 2 * sb->ny * sb->nz + 2 * (sb->nx - 2)*sb->nz + 2 * (sb->nx - 2)*(sb->ny - 2);
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
int
mdl_normalize_elements(struct mdlparse_vars *parse_state,
                       struct region *reg,
                       int existing)
{
  struct bit_array *temp = NULL;
  char op;
  unsigned int num_elems;

  if (reg->element_list_head == NULL) {
    return 0;
  }

  struct polygon_object *poly_obj = NULL;
  if (reg->parent->object_type == BOX_OBJ)
  {
    poly_obj = (struct polygon_object*) reg->parent->contents;
    num_elems = count_cuboid_elements(poly_obj->sb);
  }
  else {
    num_elems = reg->parent->n_walls;
  }

  struct bit_array *elem_array;
  if (reg->membership == NULL)
  {
    elem_array = new_bit_array(num_elems);
    if (elem_array==NULL)
    {
      mcell_allocfailed("Failed to allocate a region membership bitmask.");
      return 1;
    }
    reg->membership = elem_array;
  }
  else {
    elem_array = reg->membership;
  }

  if (reg->element_list_head->special == NULL) {
    set_all_bits(elem_array, 0);
  }
  // Special flag for exclusion
  else if ((void*) reg->element_list_head->special == (void*) reg->element_list_head) {
    set_all_bits(elem_array, 1);
  }
  else
  {
    if (reg->element_list_head->special->exclude) {
      set_all_bits(elem_array, 1);
    }
    else {
      set_all_bits(elem_array, 0);
    }
  }

  int i = 0;
  struct element_list *elem_list;
  for (elem_list = reg->element_list_head; elem_list != NULL; elem_list = elem_list->next)
  {
    if (reg->parent->object_type == BOX_OBJ)
    {
      assert(poly_obj != NULL);
      i = elem_list->begin;
      switch(i)
      {
        case X_NEG:
          elem_list->begin = 0;
          elem_list->end = 2 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) - 1;
          break;
        case X_POS:
          elem_list->begin = 2 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1);
          elem_list->end = 4 * (poly_obj->sb->ny-1) * (poly_obj->sb->nz - 1) - 1;
          break;
        case Y_NEG:
          elem_list->begin = 4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1);
          elem_list->end = elem_list->begin + 2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1) - 1;
          break;
        case Y_POS:
          elem_list->begin = 4 * (poly_obj->sb->ny - 1)*(poly_obj->sb->nz-1) + 2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1);
          elem_list->end = elem_list->begin + 2 * (poly_obj->sb->nx-1) * (poly_obj->sb->nz-1) - 1;
          break;
        case Z_NEG:
          elem_list->begin = 4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) + 4 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1);
          elem_list->end = elem_list->begin + 2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->ny-1) - 1;
          break;
        case Z_POS:
          elem_list->end = num_elems - 1;
          elem_list->begin = elem_list->end + 1 - 2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->ny - 1);
          break;
        case ALL_SIDES:
          elem_list->begin = 0;
          elem_list->end = num_elems - 1;
          break;
        default:
          UNHANDLED_CASE(i);
          return 1;
      }
    }
    else if (elem_list->begin >= (u_int) num_elems || elem_list->end >= (u_int) num_elems)
    {
      mdlerror_fmt(parse_state,
                   "Region element specifier refers to sides %u...%u, but polygon has only %u sides.",
                   elem_list->begin,
                   elem_list->end,
                   num_elems);
      return 1;
    }

    if (elem_list->special == NULL) {
      set_bit_range(elem_array, elem_list->begin, elem_list->end, 1);
    }
    else if ((void*) elem_list->special == (void*) elem_list) {
      set_bit_range(elem_array, elem_list->begin, elem_list->end, 0);
    }
    else
    {
      if (elem_list->special->referent != NULL)
      {
        if (elem_list->special->referent->membership == NULL)
        {
          if (elem_list->special->referent->element_list_head != NULL)
          {
            i = mdl_normalize_elements(parse_state, elem_list->special->referent, existing);
            if (i) {
              return i;
            }
          }
        }
        if (elem_list->special->referent->membership != NULL)
        {
          // What does it mean for the membership array to have length zero?
          if (elem_list->special->referent->membership->nbits == 0)
          {
            if (elem_list->special->exclude) {
              set_all_bits(elem_array, 0);
            }
            else {
              set_all_bits(elem_array, 1);
            }
          }
          else
          {
            if (elem_list->special->exclude) {
              op = '-';
            }
            else {
              op = '+';
            }
            bit_operation(elem_array, elem_list->special->referent->membership, op);
          }
        }
      }
      else
      {
        int ii;
        if (temp == NULL)
        {
          temp = new_bit_array(num_elems);
          if (temp == NULL)
          {
            mcell_allocfailed("Failed to allocate a region membership bitmask.");
            return 1;
          }
        }
        if (poly_obj == NULL)
        {
          mcell_internal_error("Attempt to create a PATCH on a POLYGON_LIST.");
          return 1;
        }
        if (existing)
        {
          mcell_internal_error("Attempt to create a PATCH on an already triangulated BOX.");
          return 1;
        }
        if (elem_list->special->exclude) {
          op = '-';
        }
        else {
          op = '+';
        }

        ii = cuboid_patch_to_bits(poly_obj->sb, &(elem_list->special->corner1), &(elem_list->special->corner2), temp);
        if (ii) {
          // Something wrong with patch. 
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
    bit_operation(elem_array, ((struct polygon_object*) reg->parent->contents)->side_removed, '-');
  }

  while (reg->element_list_head)
  {
    struct element_list *next = reg->element_list_head->next;
    if (reg->element_list_head->special)
    {
      if (reg->element_list_head->special != (struct element_special *) reg->element_list_head) {
        free(reg->element_list_head->special);
      }
    }
    free(reg->element_list_head);
    reg->element_list_head = next;
  }

#ifdef DEBUG
  printf("Normalized membership of %s: ", reg->sym->name);
  for (i=0; i<reg->membership->nbits; i++)
  {
    if (get_bit(reg->membership, i)) {
      printf("X");
    }
    else {
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
static int
vertex_at_index(struct subdivided_box *sb,
                int ix, int iy, int iz)
{

  if (ix==0 || ix==sb->nx-1)
  {
    int i = sb->ny * iz + iy;
    if (ix==0) return i;
    else return i + sb->ny*sb->nz;
  }
  else if (iy==0 || iy==sb->ny-1)
  {
    int i = 2*sb->ny*sb->nz + (sb->nx-2)*iz + (ix-1);
    if (iy==0) return i;
    else return i + (sb->nx-2)*sb->nz;
  }
  else if (iz==0 || iz==sb->nz-1)
  {
    int i = 2*sb->ny*sb->nz + 2*(sb->nx-2)*sb->nz + (sb->nx-2)*(iy-1) + (ix-1);
    if (iz==0) return i;
    else return i + (sb->nx-2)*(sb->ny-2);
  }
  else
  {
    mcell_internal_error("Asking for point %d %d %d but limits are [0 0 0] to [%d %d %d].",ix,iy,iz,sb->nx-1,sb->ny-1,sb->nz-1);
    return -1;
  }
}

/*************************************************************************
 polygonalize_cuboid:
 In: parse_state: parser state
     opp: an ordered polygon object that we will create
     sb: a subdivided box
 Out: returns 1 on failure, 0 on success.  The partitions along each axis of
      the subdivided box are considered to be grid lines along which we
      subdivide the surface of the box.  Walls corresponding to these surface
      elements are created and placed into the polygon_object.
*************************************************************************/
static int
polygonalize_cuboid(struct mdlparse_vars *parse_state,
                    struct polygon_object *pop,
                    struct subdivided_box *sb)
{
  struct vector3 *v;
  struct element_data *e;
  int i,j,a,b,c;
  int ii,bb,cc;
  struct vertex_list *head = NULL, *tail, *vlp;
  struct vector3 *vert_array;

  pop->n_verts = count_cuboid_vertices(sb);

  vert_array = CHECKED_MALLOC_ARRAY(struct vector3,
                                      pop->n_verts,
                                      "cuboid vertices");
  if (vert_array == NULL) return 1;

  pop->n_walls = count_cuboid_elements(sb);
  pop->element = CHECKED_MALLOC_ARRAY(struct element_data,
                                       pop->n_walls,
                                       "cuboid walls");
  if (pop->element == NULL)
  {
    free(vert_array);
    vert_array = NULL;
    return 1;
  }

/*  for (a=0;a<2;a++) for (b=0;b<2;b++) for (c=0;c<2;c++) printf("%d,%d,%d->%d\n",a,b,c,vertex_at_index(sb,a,b,c)); */

  /* Set vertices and elements on X faces */
  ii = 0;
  bb = 0;
  cc = 2*(sb->nz-1)*(sb->ny-1);
  b = 0;
  c = sb->nz*sb->ny;
  for (j=0 ; j<sb->nz ; j++)
  {
    a = sb->ny;
    for (i=0 ; i<sb->ny ; i++)
    {
      /*printf("Setting indices %d %d\n",b+j*a+i,c+j*a+i);*/
      v = &(vert_array[b+j*a+i]);
      v->x = sb->x[0];
      v->y = sb->y[i];
      v->z = sb->z[j];
      v = &(vert_array[c+j*a+i]);
      v->x = sb->x[sb->nx-1];
      v->y = sb->y[i];
      v->z = sb->z[j];

      if (i>0 && j>0)
      {
        e = &(pop->element[bb+ii]);
    e->vertex_index[0] = vertex_at_index(sb,0,i-1,j-1);
    e->vertex_index[2] = vertex_at_index(sb,0,i,j-1);
    e->vertex_index[1] = vertex_at_index(sb,0,i-1,j);
        e = &(pop->element[bb+ii+1]);
    e->vertex_index[0] = vertex_at_index(sb,0,i,j);
    e->vertex_index[1] = vertex_at_index(sb,0,i,j-1);
    e->vertex_index[2] = vertex_at_index(sb,0,i-1,j);
        e = &(pop->element[cc+ii]);
    e->vertex_index[0] = vertex_at_index(sb,sb->nx-1,i-1,j-1);
    e->vertex_index[1] = vertex_at_index(sb,sb->nx-1,i,j-1);
    e->vertex_index[2] = vertex_at_index(sb,sb->nx-1,i-1,j);
        e = &(pop->element[cc+ii+1]);
    e->vertex_index[0] = vertex_at_index(sb,sb->nx-1,i,j);
    e->vertex_index[2] = vertex_at_index(sb,sb->nx-1,i,j-1);
    e->vertex_index[1] = vertex_at_index(sb,sb->nx-1,i-1,j);
        /*printf("Setting elements %d %d %d %d of %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,pop->n_walls);*/

    ii+=2;
      }
    }
  }

  /* Set vertices and elements on Y faces */
  bb = ii;
  cc = bb + 2*(sb->nx-1)*(sb->nz-1);
  b = 2*sb->nz*sb->ny;
  c = b + sb->nz*(sb->nx-2);
  for (j=0 ; j<sb->nz ; j++)
  {
    a = sb->nx-2;
    for (i=1 ; i<sb->nx ; i++)
    {
      if (i<sb->nx-1)
      {
        /*printf("Setting indices %d %d of %d\n",b+j*a+(i-1),c+j*a+(i-1),pop->n_verts);*/
        v = &(vert_array[b+j*a+(i-1)]);
    v->x = sb->x[i];
    v->y = sb->y[0];
    v->z = sb->z[j];
        v = &(vert_array[c+j*a+(i-1)]);
    v->x = sb->x[i];
    v->y = sb->y[sb->ny-1];
    v->z = sb->z[j];
      }

      if (j>0)
      {
        e = &(pop->element[bb+ii]);
    e->vertex_index[0] = vertex_at_index(sb,i-1,0,j-1);
    e->vertex_index[1] = vertex_at_index(sb,i,0,j-1);
    e->vertex_index[2] = vertex_at_index(sb,i-1,0,j);
        e = &(pop->element[bb+ii+1]);
    e->vertex_index[0] = vertex_at_index(sb,i,0,j);
    e->vertex_index[2] = vertex_at_index(sb,i,0,j-1);
    e->vertex_index[1] = vertex_at_index(sb,i-1,0,j);
        e = &(pop->element[cc+ii]);
    e->vertex_index[0] = vertex_at_index(sb,i-1,sb->ny-1,j-1);
    e->vertex_index[2] = vertex_at_index(sb,i,sb->ny-1,j-1);
    e->vertex_index[1] = vertex_at_index(sb,i-1,sb->ny-1,j);
        e = &(pop->element[cc+ii+1]);
    e->vertex_index[0] = vertex_at_index(sb,i,sb->ny-1,j);
    e->vertex_index[1] = vertex_at_index(sb,i,sb->ny-1,j-1);
    e->vertex_index[2] = vertex_at_index(sb,i-1,sb->ny-1,j);
        /*printf("Setting elements %d %d %d %d of %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,pop->n_walls);*/

    ii+=2;
      }
    }
  }

  /* Set vertices and elements on Z faces */
  bb = ii;
  cc = bb + 2*(sb->nx-1)*(sb->ny-1);
  b = 2*sb->nz*sb->ny + 2*(sb->nx-2)*sb->nz;
  c = b + (sb->nx-2)*(sb->ny-2);
  for (j=1 ; j<sb->ny ; j++)
  {
    a = sb->nx-2;
    for (i=1 ; i<sb->nx ; i++)
    {
      if (i<sb->nx-1 && j<sb->ny-1)
      {
        /*printf("Setting indices %d %d of %d\n",b+(j-1)*a+(i-1),c+(j-1)*a+(i-1),pop->n_verts);*/
        v = &(vert_array[b+(j-1)*a+(i-1)]);
    v->x = sb->x[i];
    v->y = sb->y[j];
    v->z = sb->z[0];
        v = &(vert_array[c+(j-1)*a+(i-1)]);
    v->x = sb->x[i];
    v->y = sb->y[j];
    v->z = sb->z[sb->nz-1];
      }

      e = &(pop->element[bb+ii]);
      e->vertex_index[0] = vertex_at_index(sb,i-1,j-1,0);
      e->vertex_index[2] = vertex_at_index(sb,i,j-1,0);
      e->vertex_index[1] = vertex_at_index(sb,i-1,j,0);
      e = &(pop->element[bb+ii+1]);
      e->vertex_index[0] = vertex_at_index(sb,i,j,0);
      e->vertex_index[1] = vertex_at_index(sb,i,j-1,0);
      e->vertex_index[2] = vertex_at_index(sb,i-1,j,0);
      e = &(pop->element[cc+ii]);
      e->vertex_index[0] = vertex_at_index(sb,i-1,j-1,sb->nz-1);
      e->vertex_index[1] = vertex_at_index(sb,i,j-1,sb->nz-1);
      e->vertex_index[2] = vertex_at_index(sb,i-1,j,sb->nz-1);
      e = &(pop->element[cc+ii+1]);
      e->vertex_index[0] = vertex_at_index(sb,i,j,sb->nz-1);
      e->vertex_index[2] = vertex_at_index(sb,i,j-1,sb->nz-1);
      e->vertex_index[1] = vertex_at_index(sb,i-1,j,sb->nz-1);

      /*printf("Setting elements %d %d %d %d of %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,pop->n_walls);*/

      ii+=2;
    }
  }


   /* build the head node of the linked list "pop->parsed_vertices" */

  vlp = CHECKED_MALLOC_STRUCT(struct vertex_list, "vertex_list");
  if (vlp == NULL) return 1;
  vlp->vertex = CHECKED_MALLOC_STRUCT(struct vector3, "vertex");
  if (vlp->vertex == NULL) return 1;
  memcpy(vlp->vertex, &vert_array[0], sizeof(struct vector3));
  vlp->next = head;
  head = vlp;
  tail = head;

  /* build other nodes of the linked list "pop->parsed_vertices" */
  for (i = 1; i < pop->n_verts; i++)
  {
    vlp = CHECKED_MALLOC_STRUCT(struct vertex_list, "vertex_list");
    if (vlp == NULL) return 1;
    vlp->vertex = CHECKED_MALLOC_STRUCT(struct vector3, "vertex");
    if (vlp->vertex == NULL) return 1;
    memcpy(vlp->vertex, &vert_array[i], sizeof(struct vector3));
    vlp->next = tail->next;
    tail->next = vlp;
    tail = tail->next;

  }
  pop->parsed_vertices = head;

#ifdef DEBUG
  printf("BOX has vertices:\n");
  for (i=0;i<pop->n_verts;i++) printf("  %.5e %.5e %.5e\n",vert_array[i].x,vert_array[i].y,vert_array[i].z);
  printf("BOX has walls:\n");
  for (i=0;i<pop->n_walls;i++) printf("  %d %d %d\n",pop->element[i].vertex_index[0],pop->element[i].vertex_index[1],pop->element[i].vertex_index[2]);
  printf("\n");
#endif

  if (vert_array != NULL) free(vert_array);
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
int
mdl_triangulate_box_object(struct mdlparse_vars *parse_state,
                           struct sym_table *box_sym,
                           struct polygon_object *pop,
                           double box_aspect_ratio)
{
  struct region_list *rlp;
  struct object *objp = (struct object *) box_sym->value;

  if (box_aspect_ratio >= 2.0)
  {
    if (reaspect_cuboid(parse_state, pop->sb, box_aspect_ratio))
    {
      mdlerror(parse_state, "Error setting up box geometry");
      return 1;
    }
  }
  for (rlp = objp->regions; rlp != NULL; rlp = rlp->next)
  {
    if (mdl_normalize_elements(parse_state, rlp->reg, 0))
      return 1;
  }
  if (polygonalize_cuboid(parse_state, pop, pop->sb))
  {
    mdlerror(parse_state, "Could not turn box object into polygons");
    return 1;
  }
  else if (parse_state->vol->notify->box_triangulation==NOTIFY_FULL)
  {
    mcell_log("Box object %s converted into %d polygons.", box_sym->name, pop->n_walls);
  }

  const unsigned int n_walls = pop->n_walls;
  pop->side_removed = new_bit_array(n_walls);
  if (pop->side_removed==NULL)
  {
    mcell_allocfailed("Failed to allocate a box object removed side bitmask.");
    return 1;
  }
  set_all_bits(pop->side_removed,0);
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
int
mdl_check_diffusion_constant(struct mdlparse_vars *parse_state, double *d)
{
  if (parse_state->vol->notify->neg_diffusion==WARN_COPE)
  {
    if (*d < 0) *d = 0.0;
  }
  else if (parse_state->vol->notify->neg_diffusion==WARN_WARN)
  {
    if (*d < 0.0)
    {
      mdlerror_fmt(parse_state, "Negative diffusion constant found, setting to zero and continuing.");
      *d = 0.0;
    }
  }
  else
  {
    if (*d < 0.0)
    {
      mdlerror(parse_state, "Error: diffusion constants should be zero or positive.");
      return 1;
    }
  }
  return 0;
}

/*************************************************************************
 mdl_print_species_summary:
    Finish the creation of a molecule, undoing any state changes we made during
    the creation of the molecule.  Presently, this just means "print the
    diffusion distances report".

 In:  parse_state: parser state
      mol:  species finished
 Out: A report is printed to the file handle.
 NOTE: This is just a thin wrapper around print_species_summary
*************************************************************************/
void
mdl_print_species_summary(struct mdlparse_vars *parse_state, struct species *mol)
{
  print_species_summary(parse_state->vol, mol);
}

/*************************************************************************
 mdl_print_species_summaries:
    Finish the creation of a series of molecules, undoing any state changes we
    made during the creation of the molecules.  Presently, this just means
    "print the diffusion distances report".

 In:  parse_state: parser state
      mols: species finished
 Out: A report is printed to the file handle.
 NOTE: This is just a thin wrapper around print_species_summaries
*************************************************************************/
void
mdl_print_species_summaries(struct mdlparse_vars *parse_state,
                            struct species_list_item *mols)
{
  print_species_summaries(parse_state->vol, mols);
  mem_put_list(parse_state->species_list_mem, mols);
}

/*************************************************************************
 mdl_add_to_species_list:
    Add a single species to a species list.

 In:  parse_state: parser state
      list: the list
      spec: the species
 Out: 0 on success, 1 on failure
 NOTE: This is just a thin wrapper around add_to_species_list
*************************************************************************/
int
mdl_add_to_species_list(struct mdlparse_vars *parse_state,
                        struct species_list *list,
                        struct species *spec)
{
  return add_to_species_list(parse_state->species_list_mem, list, spec);
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
int
mdl_start_release_site(struct mdlparse_vars *parse_state,
                       struct sym_table *symp,
                       int shape)
{
  struct object *obj_ptr = start_release_site(parse_state->vol, symp);
  parse_state->current_release_site = obj_ptr->contents;
  if (obj_ptr->contents == NULL) {
    return 1;
  }

  parse_state->current_release_site->release_shape = (int8_t) shape;
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
struct object *
mdl_finish_release_site(struct mdlparse_vars *parse_state,
                        struct sym_table *symp)
{
  struct object *objp_new = finish_release_site(symp);
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
int
mdl_is_release_site_valid(struct mdlparse_vars *parse_state,
                          struct release_site_obj *rel_site_obj_ptr)
{
  switch (is_release_site_valid(rel_site_obj_ptr))
  {
    case 2:
      mdlerror(
        parse_state,
        "Must specify molecule to release using MOLECULE=molecule_name.");
      return 1;
    case 3:
      mdlerror_fmt(
        parse_state,
        "Cannot release surface class '%s' from release site",
        rel_site_obj_ptr->mol_type->sym->name);
      return 1;
    case 4:
      mdlerror(
        parse_state,
        "CONCENTRATION may only be used with molecules that can diffuse in 3D.\n"
        "  Use DENSITY for molecules diffusing in 2D.");
      return 1;
    case 5:
      mdlerror(
        parse_state,
        "DENSITY may only be used with molecules that can diffuse in 2D.\n"
        "  Use CONCENTRATION for molecules diffusing in 3D.");
      return 1;
    case 6:
      mdlerror(parse_state, "Release site is missing location.");
      return 1;
  }
  return 0;
}



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
static int
mdl_check_release_regions(struct mdlparse_vars *parse_state,
                          struct release_evaluator *rel_eval,
                          struct object *parent,
                          struct object *instance)
{
  switch (check_release_regions(rel_eval, parent, instance))
  {
    case 1:
      return 1;
    case 2:
      mdlerror(parse_state, "Region neither instanced nor grouped with release site.");
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
                                     struct release_evaluator *rel_eval)
{
  switch (set_release_site_geometry_region(
    parse_state->vol, rel_site_obj_ptr, obj_ptr, rel_eval))
  {
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
                                     struct object *obj_ptr)
{
  if ((obj_ptr->object_type == META_OBJ) || (obj_ptr->object_type == REL_SITE_OBJ))
  {
    mdlerror(
      parse_state,
      "Error: only BOX or POLYGON_LIST objects may be assigned to the SHAPE "
      "keyword in the RELEASE_SITE definition.  Metaobjects or release objects"
      " are not allowed here.");
    return 1;
  }

  char *obj_name = obj_ptr->sym->name;
  char *region_name = CHECKED_SPRINTF("%s,ALL", obj_name);
  if (region_name == NULL) {
    return 1;
  }
  struct sym_table *sym_ptr;
  if ((sym_ptr = retrieve_sym(
    region_name, parse_state->vol->reg_sym_table)) == NULL)
  {
    mdlerror_fmt(parse_state, "Undefined region: %s", region_name);
    free(region_name);
    return 1;
  }
  free(region_name);

  struct release_evaluator *rel_eval = CHECKED_MALLOC_STRUCT(
    struct release_evaluator, "release site on region");
  if (rel_eval==NULL) {
    return 1;
  }

  rel_eval->op = REXP_NO_OP | REXP_LEFT_REGION;
  rel_eval->left = sym_ptr->value;
  rel_eval->right = NULL;

  ((struct region*)rel_eval->left)->flags |= COUNT_CONTENTS;

  rel_site_obj_ptr->release_shape = SHAPE_REGION;
  parse_state->vol->place_waypoints_flag = 1;

  if (mdl_check_release_regions(
    parse_state, rel_eval, obj_ptr, parse_state->vol->root_instance))
  {
    mdlerror(
      parse_state,
      "Trying to release on a region that the release site cannot see!\n  "
      "Try grouping the release site and the corresponding geometry with an "
      "OBJECT.");
    return 1;
  }

  struct release_region_data *rel_reg_data = CHECKED_MALLOC_STRUCT(
    struct release_region_data, "release site on region");
  if (rel_reg_data==NULL)
  {
    mdlerror(
      parse_state,
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
static int
mdl_check_valid_molecule_release(struct mdlparse_vars *parse_state,
                                 struct species_opt_orient *mol_type)
{
  static const char *EXTRA_ORIENT_MSG =
        "surface orientation not specified for released surface molecule\n"
        "(use ; or ', or ,' for random orientation)";
  static const char *MISSING_ORIENT_MSG =
        "orientation not used for released volume molecule";

  struct species *mol = (struct species *) mol_type->mol_type->value;
  if (mol->flags & ON_GRID)
  {
    if (!mol_type->orient_set)
    {
      if (parse_state->vol->notify->missed_surf_orient==WARN_ERROR)
      {
        mdlerror_fmt(parse_state, "Error: %s", EXTRA_ORIENT_MSG);
        return 1;
      }
      else if (parse_state->vol->notify->missed_surf_orient==WARN_WARN) {
        mdlerror_fmt(parse_state, "Warning: %s", EXTRA_ORIENT_MSG);
      }
    }
  }
  else if ((mol->flags & NOT_FREE) == 0)
  {
    if (mol_type->orient_set)
    {
      if (parse_state->vol->notify->useless_vol_orient==WARN_ERROR)
      {
        mdlerror_fmt(parse_state, "Error: %s", MISSING_ORIENT_MSG);
        return 1;
      }
      else if (parse_state->vol->notify->useless_vol_orient==WARN_WARN) {
        mdlerror_fmt(parse_state, "Warning: %s", MISSING_ORIENT_MSG);
      }
    }
  }
  else
  {
    mdlerror(parse_state, "Error: cannot release a surface class instead of a molecule.");
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
int
mdl_set_release_site_molecule(struct mdlparse_vars *parse_state,
                              struct release_site_obj *rel_site_obj_ptr,
                              struct species_opt_orient *mol_type)
{
  // Store molecule information
  rel_site_obj_ptr->mol_type = (struct species *) mol_type->mol_type->value;
  if ((rel_site_obj_ptr->mol_type->flags & NOT_FREE) == 0)
  {
    if (rel_site_obj_ptr->release_shape == SHAPE_REGION)
      parse_state->vol->place_waypoints_flag = 1;
  }
  else
  {
    if (rel_site_obj_ptr->release_shape != SHAPE_REGION  &&
        rel_site_obj_ptr->release_shape != SHAPE_LIST)
    {
      mdlerror_fmt(
        parse_state,
        "The release site '%s' is a geometric release site, and may not be used to \n"
        "  release the surface molecule '%s'.  Surface molecule release sites must \n"
        "  be either LIST or region release sites.",
        rel_site_obj_ptr->name,
        mol_type->mol_type->name);
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
int
mdl_set_release_site_diameter(struct mdlparse_vars *parse_state,
                              struct release_site_obj *rel_site_obj_ptr,
                              double diam)
{
  diam *= parse_state->vol->r_length_unit;

  rel_site_obj_ptr->diameter = CHECKED_MALLOC_STRUCT(struct vector3,
                                                     "release site diameter");
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
                                    int n_diams,
                                    struct num_expr_list *diams,
                                    double factor)
{
  factor *= parse_state->vol->r_length_unit;

  if (rel_site_obj_ptr->release_shape == SHAPE_LIST)
  {
    mdlerror(parse_state, "Release list diameters must be single valued.");
    return 1;
  }

  if (n_diams != 3)
  {
    mdlerror(parse_state, "Three dimensional value required");
    return 1;
  }

  rel_site_obj_ptr->diameter = CHECKED_MALLOC_STRUCT(struct vector3,
                                                     "release site diameter");
  if (rel_site_obj_ptr->diameter == NULL) {
    return 1;
  }
  rel_site_obj_ptr->diameter->x = diams->value             * factor;
  rel_site_obj_ptr->diameter->y = diams->next->value       * factor;
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
int
mdl_set_release_site_diameter_var(struct mdlparse_vars *parse_state,
                                  struct release_site_obj *rel_site_obj_ptr,
                                  double factor,
                                  struct sym_table *symp)
{
  struct num_expr_list *expr_list_ptr;
  int count = 0;
  rel_site_obj_ptr->diameter = CHECKED_MALLOC_STRUCT(
    struct vector3, "release site diameter");
  if (rel_site_obj_ptr->diameter == NULL) {
    return 1;
  }

  switch (symp->sym_type)
  {
    case DBL:
      if (mdl_set_release_site_diameter(
          parse_state, rel_site_obj_ptr, *(double *) symp->value * factor))
      {
        return 1;
      }
      break;

    case ARRAY:
      // Count up to 4 elements -- that's all we need to count to know if it's valid
      for (expr_list_ptr = (struct num_expr_list *) symp->value; expr_list_ptr != NULL && count < 4; ++ count, expr_list_ptr = expr_list_ptr->next)
        ;
      expr_list_ptr = (struct num_expr_list *) symp->value;
      if (mdl_set_release_site_diameter_array(
          parse_state, rel_site_obj_ptr, count, expr_list_ptr, factor))
      {
        return 1;
      }
      break;

    default:
      mdlerror(
        parse_state,
        "Diameter must either be a number or a 3-valued vector.");
      return 1;
  }

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
int
mdl_set_release_site_probability(struct mdlparse_vars *parse_state,
                                 struct release_site_obj *rel_site_obj_ptr,
                                 double prob)
{
  if (rel_site_obj_ptr->release_prob == MAGIC_PATTERN_PROBABILITY)
  {
    mdlerror(parse_state,
             "Ignoring release probability for reaction-triggered releases.");
  }
  else
  {
    rel_site_obj_ptr->release_prob = prob;
    if (rel_site_obj_ptr->release_prob < 0)
    {
      mdlerror(parse_state, "Release probability cannot be less than 0.");
      return 1;
    }
    if (rel_site_obj_ptr->release_prob>1)
    {
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
int
mdl_set_release_site_pattern(struct mdlparse_vars *parse_state,
                             struct release_site_obj *rel_site_obj_ptr,
                             struct sym_table *pattern)
{
  rel_site_obj_ptr->pattern = (struct release_pattern *) pattern->value;

  // Careful!  We've put a rxn_pathname into the "pattern" pointer! 
  if (pattern->sym_type == RXPN) 
  {
    if (rel_site_obj_ptr->release_prob != 1.0)
    {
      mdlerror(
        parse_state,
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
int
mdl_set_release_site_molecule_positions(struct mdlparse_vars *parse_state,
                                        struct release_site_obj *rel_site_obj_ptr,
                                        struct release_single_molecule_list *list)
{
  if (rel_site_obj_ptr->release_shape != SHAPE_LIST)
  {
    mdlerror(parse_state,
            "You must use the LIST shape to specify molecule positions in a "
            "release.");
    return 1;
  }

  struct release_single_molecule *rsm;
  if (rel_site_obj_ptr->mol_list == NULL) {
    rel_site_obj_ptr->mol_list = list->rsm_head;
  }
  else
  {
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
                                struct species_opt_orient *mol_type,
                                struct vector3 *pos)
{
  struct release_single_molecule *rsm;
  struct vector3 temp_v3;

  memcpy(&temp_v3,pos,sizeof(struct vector3));
  free(pos);

  rsm = CHECKED_MALLOC_STRUCT(struct release_single_molecule, "release site molecule position");
  if (rsm==NULL)
  {
    mdlerror(parse_state, "Out of memory reading molecule positions");
    return NULL;
  }

  rsm->orient = mol_type->orient;
  rsm->loc.x = temp_v3.x * parse_state->vol->r_length_unit;
  rsm->loc.y = temp_v3.y * parse_state->vol->r_length_unit;
  rsm->loc.z = temp_v3.z * parse_state->vol->r_length_unit;
  rsm->mol_type = (struct species*)(mol_type->mol_type->value);
  rsm->next = NULL;

  if (mdl_check_valid_molecule_release(parse_state, mol_type))
  {
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
                                   double conc)
{
  if (set_release_site_concentration(rel_site_obj_ptr, conc))
  {
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
void 
mdl_vertex_list_singleton(struct vertex_list_head *head,
                          struct vertex_list *item)
{
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
void 
mdl_add_vertex_to_list(struct vertex_list_head *head,
                       struct vertex_list *item)
{
  item->next = NULL;
  head->vertex_tail = head->vertex_tail->next = item;
  ++ head->vertex_count;
}

/**************************************************************************
 mdl_new_vertex_list_item:
    Allocate an item for a vertex list.

 In: vertex: this vertex
     normal: surface normal at this vertex, or NULL
 Out: the vertex list item, or NULL if an error occurred
**************************************************************************/
struct vertex_list *
mdl_new_vertex_list_item(struct vector3 *vertex)
{

  struct vertex_list *vlp = CHECKED_MALLOC_STRUCT(struct vertex_list,
                                                   "vertices");
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
                                      struct element_connection_list *item)
{
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
                                   struct element_connection_list *item)
{
  item->next = NULL;
  head->connection_tail = head->connection_tail->next = item;
  ++ head->connection_count;
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
                           struct num_expr_list_head *indices)
{
  if (indices->value_count != 3)
  {
    mdlerror(parse_state, "Non-triangular element found in polygon list object");
    return NULL;
  }

  struct element_connection_list *eclp = CHECKED_MALLOC_STRUCT(struct element_connection_list,
                                                                "polygon element commections");
  if (eclp == NULL)
    return NULL;

  eclp->indices = CHECKED_MALLOC_ARRAY(int, 3, "polygon element connections");
  if (eclp->indices == NULL)
  {
    free(eclp);
    return NULL;
  }
  eclp->indices[0] = (int) indices->value_head->value;
  eclp->indices[1] = (int) indices->value_head->next->value;
  eclp->indices[2] = (int) indices->value_tail->value;
  eclp->n_verts = indices->value_count;
  eclp->next = NULL;

  if (! indices->shared)
    mdl_free_numeric_list(indices->value_head);
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
                               struct num_expr_list_head *indices)
{
  if (indices->value_count != 4)
  {
    mdlerror(parse_state, "Non-tetrahedron element found in voxel list object");
    return NULL;
  }

  struct element_connection_list *eclp = CHECKED_MALLOC_STRUCT(struct element_connection_list,
                                                                "polygon element commections");
  if (eclp == NULL)
    return NULL;

  eclp->indices = CHECKED_MALLOC_ARRAY(int, 4, "polygon element connections");
  if (eclp->indices == NULL)
  {
    free(eclp);
    return NULL;
  }
  eclp->indices[0] = (int) indices->value_head->value;
  eclp->indices[1] = (int) indices->value_head->next->value;
  eclp->indices[2] = (int) indices->value_head->next->next->value;
  eclp->indices[3] = (int) indices->value_tail->value;
  eclp->n_verts = indices->value_count;
  eclp->next = NULL;

  if (! indices->shared)
    mdl_free_numeric_list(indices->value_head);
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
struct polygon_object *
mdl_new_polygon_list(struct mdlparse_vars *parse_state,
                     struct sym_table *sym,
                     int n_vertices,
                     struct vertex_list *vertices,
                     int n_connections,
                     struct element_connection_list *connections)
{

  struct polygon_object *poly_obj_ptr = allocate_polygon_object("polygon list object");
  if (poly_obj_ptr == NULL) {
    goto failure;
  }

  struct object *obj_ptr = (struct object *) sym->value;
  obj_ptr->object_type = POLY_OBJ;
  obj_ptr->contents = poly_obj_ptr;

  poly_obj_ptr->n_walls = n_connections;
  poly_obj_ptr->n_verts = n_vertices;

  // Allocate and initialize removed sides bitmask
  poly_obj_ptr->side_removed = new_bit_array(poly_obj_ptr->n_walls);
  if (poly_obj_ptr->side_removed==NULL)
  {
    mcell_allocfailed("Failed to allocate a polygon list object removed side bitmask.");
    goto failure;
  }
  set_all_bits(poly_obj_ptr->side_removed, 0);

  /* Keep temporarily information about vertices in the form of
     "parsed_vertices" */
  poly_obj_ptr->parsed_vertices = vertices;

  // Copy in vertices and normals
  struct vertex_list *vert_list = poly_obj_ptr->parsed_vertices;
  for (int i = 0; i < poly_obj_ptr->n_verts; i++)
  {
    // Rescale vertices coordinates
    vert_list->vertex->x *= parse_state->vol->r_length_unit;
    vert_list->vertex->y *= parse_state->vol->r_length_unit;
    vert_list->vertex->z *= parse_state->vol->r_length_unit;

    vert_list = vert_list->next;
  }

  // Allocate wall elements
  struct element_data *elem_data_ptr = NULL;
  if ((elem_data_ptr = CHECKED_MALLOC_ARRAY(
        struct element_data,
        poly_obj_ptr->n_walls,
        "polygon list object walls")) == NULL) {
    goto failure;
  }
  poly_obj_ptr->element = elem_data_ptr;

  // Copy in wall elements 
  for (int i = 0; i<poly_obj_ptr->n_walls; i++)
  {
    if (connections->n_verts != 3)
    {
      mdlerror(parse_state, "All polygons must have three vertices.");
      goto failure;
    }

    struct element_connection_list *elem_conn_list_temp = connections;
    memcpy(elem_data_ptr[i].vertex_index, connections->indices, 3*sizeof(int));
    connections = connections->next;
    free(elem_conn_list_temp->indices);
    free(elem_conn_list_temp);
  }

  // Create object default region on polygon list object: 
  struct region *reg_ptr = NULL;
  if ((reg_ptr = mdl_create_region(parse_state, obj_ptr, "ALL")) == NULL)
    goto failure;
  if ((reg_ptr->element_list_head = new_element_list(0, poly_obj_ptr->n_walls - 1)) == NULL)
    goto failure;

  obj_ptr->n_walls = poly_obj_ptr->n_walls;
  obj_ptr->n_verts = poly_obj_ptr->n_verts;
  if (mdl_normalize_elements(parse_state, reg_ptr, 0))
  {
    mdlerror_fmt(parse_state,
                 "Error setting up elements in default 'ALL' region in the polygon object '%s'.",
                 sym->name);
    goto failure;
  }

  parse_state->allow_patches = 0;
  return poly_obj_ptr;

failure:
  free_connection_list(connections);
  free_vertex_list(vertices);
  if (poly_obj_ptr)
  {
    if (poly_obj_ptr->element)
      free(poly_obj_ptr->element);
    if (poly_obj_ptr->side_removed)
      free_bit_array(poly_obj_ptr->side_removed);
    free(poly_obj_ptr);
  }
  return NULL;
}



/**************************************************************************
 mdl_finish_polygon_list:
    Finalize the polygon list, cleaning up any state updates that were made
    when we started creating the polygon.

 In: parse_state: parser state
     symp: symbol for the completed polygon
 Out: 1 on failure, 0 on success
**************************************************************************/
int 
mdl_finish_polygon_list(struct mdlparse_vars *parse_state, struct sym_table *sym_ptr)
{
  if (finish_polygon_list(sym_ptr))
  {
    parse_state->current_polygon = NULL;
    return 1;
  }
  parse_state->current_polygon = NULL;
  return 0;
}

/**************************************************************************
 allocate_voxel_object:
    Create a new voxel object.

 In:  Nothing 
 Out: voxel object, or NULL if allocation fails
**************************************************************************/
static struct voxel_object *
allocate_voxel_object()
{

  struct voxel_object *vop;
  if ((vop = CHECKED_MALLOC_STRUCT(struct voxel_object,
                                    "voxel list object")) == NULL)
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
mdl_new_voxel_list(struct mdlparse_vars *parse_state,
                   struct sym_table *sym,
                   int n_vertices,
                   struct vertex_list *vertices,
                   int n_connections,
                   struct element_connection_list *connections)
{
  struct tet_element_data *tedp;

  struct object *objp = (struct object *) sym->value;
  struct voxel_object *vop = allocate_voxel_object();
  if (vop == NULL)
    goto failure;

  objp->object_type = VOXEL_OBJ;
  objp->contents = vop;

  vop->n_voxels = n_connections;
  vop->n_verts  = n_vertices;

  /* Allocate vertices */
  if ((vop->vertex = CHECKED_MALLOC_ARRAY(struct vector3,
                                           vop->n_verts,
                                           "voxel list object vertices")) == NULL)
    goto failure;

  /* Populate vertices */
  for (int i=0;i<vop->n_verts;i++)
  {
    struct vertex_list *vlp_temp = vertices;
    vop->vertex[i].x = vertices->vertex->x;
    vop->vertex[i].y = vertices->vertex->y;
    vop->vertex[i].z = vertices->vertex->z;
    free(vertices->vertex);
    vertices = vertices->next;
    free(vlp_temp);
  }

  /* Allocate tetrahedra */
  if ((tedp = CHECKED_MALLOC_ARRAY(struct tet_element_data,
                                    vop->n_voxels,
                                    "voxel list object tetrahedra")) == NULL)
    goto failure;
  vop->element = tedp;

  /* Copy in tetrahedra */
  for (int i=0; i<vop->n_voxels; i++)
  {
    if (connections->n_verts != 4)
    {
      mdlerror(parse_state, "All voxels must have four vertices.");
      goto failure;
    }

    struct element_connection_list *eclp_temp = connections;
    memcpy(tedp[i].vertex_index, connections->indices, 4*sizeof(int));
    connections = connections->next;
    free(eclp_temp);
  }
  return vop;

failure:
  free_vertex_list(vertices);
  free_connection_list(connections);
  if (vop)
  {
    if (vop->element)
      free(vop->element);
    if (vop->vertex)
      free(vop->vertex);
    free(vop);
  }
  return NULL;
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
struct polygon_object *
mdl_new_box_object(struct mdlparse_vars *parse_state,
                   struct sym_table *sym,
                   struct vector3 *llf,
                   struct vector3 *urb)
{
  struct polygon_object *pop;
  struct region *rp;
  struct object *objp = (struct object *) sym->value;

  /* Allocate polygon object */
  pop = allocate_polygon_object("box object");
  if (pop == NULL)
  {
    free(llf);
    free(urb);
    return NULL;
  }
  objp->object_type = BOX_OBJ;
  objp->contents = pop;

  /* Create object default region on box object: */
  if ((rp = mdl_create_region(parse_state, objp, "ALL")) == NULL)
  {
    free(pop);
    free(llf);
    free(urb);
    return NULL;
  }
  if ((rp->element_list_head = new_element_list(ALL_SIDES, ALL_SIDES)) == NULL)
  {
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
  if (pop->sb == NULL)
  {
    free(pop);
    free(llf);
    free(urb);
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
int 
mdl_finish_box_object(struct mdlparse_vars *parse_state,
                      struct sym_table *symp)
{
  struct object *objp = (struct object *) symp->value;
  remove_gaps_from_regions(objp);
  objp->n_walls = parse_state->current_polygon->n_walls;
  objp->n_verts = parse_state->current_polygon->n_verts;
  if (check_degenerate_polygon_list(objp))
  {
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
struct region *
mdl_create_region(struct mdlparse_vars *parse_state, struct object *objp, char *name)
{
  struct region *rp;
  struct region_list *rlp;
  no_printf("Creating new region: %s\n", name);
  if ((rp = mdl_make_new_region(parse_state, objp->sym->name, name)) == NULL)
    return NULL;
  if ((rlp = CHECKED_MALLOC_STRUCT(struct region_list, "region list")) == NULL)
  {
    mdlerror_fmt(parse_state,
                 "Out of memory while creating object region '%s'",
                 rp->sym->name);
    return NULL;
  }
  rp->region_last_name = name;
  rp->parent = objp;
  rlp->reg = rp;
  rlp->next = objp->regions;
  objp->regions = rlp;
  objp->num_regions++;
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
struct region *
mdl_get_region(struct mdlparse_vars *parse_state,
               struct object *objp,
               char *name)
{
  struct sym_table *reg_sym;
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
    rp = (struct region *) reg_sym->value;

  return rp;
}

/**************************************************************************
 mdl_start_existing_obj_region_def:
    Begin construction of a region on an existing object.

 In: parse_state: parser state
     obj_symp: symbol of object upon which to create region
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_start_existing_obj_region_def(struct mdlparse_vars *parse_state,
                                  struct sym_table *obj_symp)
{
  struct object *objp = (struct object *) obj_symp->value;
  if (objp->object_type != BOX_OBJ && objp->object_type != POLY_OBJ)
  {
    mdlerror_fmt(parse_state, "Cannot define region on non-surface object: %s", obj_symp->name);
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
void 
mdl_add_elements_to_list(struct element_list_head *list,
                         struct element_list *head,
                         struct element_list *tail)
{
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
void 
mdl_set_elements_to_exclude(struct element_list *els)
{
  /* HACK: els->special actually points to an element_list.  Don't dereference
   * it now. */
  for (; els != NULL; els = els->next)
    els->special = (struct element_special *) els;
}

/**************************************************************************
 mdl_new_element_side:
    Create a new element list for a region description based on a side name.

 In: parse_state: parser state
     side: side name constant (ALL_SIDES, etc.)
 Out: element list, or NULL if allocation fails
**************************************************************************/
struct element_list *
mdl_new_element_side(struct mdlparse_vars *parse_state,
                     unsigned int side)
{
  unsigned int begin, end;
  if (side == ALL_SIDES  &&  parse_state->current_object->object_type == POLY_OBJ)
  {
    begin = 0;
    end = parse_state->current_polygon->n_walls-1;
  }
  else if (parse_state->current_object->object_type==POLY_OBJ)
  {
    mdlerror(parse_state, "Illegal reference to polygon list element by side-name");
    return NULL;
  }
  else 
  {
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
struct element_list *
mdl_new_element_previous_region(struct mdlparse_vars *parse_state,
                                struct object *objp,
                                struct region *rp_container,
                                char *name_region_referent,
                                int exclude)
{
  struct sym_table *stp;
  char *full_reg_name = NULL;
  struct element_list *elmlp = NULL;

  /* Create element list */
  elmlp = new_element_list(0, 0);
  if (elmlp == NULL)
    goto failure;

  /* Create "special" element description */
  elmlp->special = CHECKED_MALLOC_STRUCT(struct element_special, "region element");
  if (elmlp->special==NULL)
    goto failure;
  elmlp->special->exclude = (byte) exclude;

  /* Create referent region full name */
  full_reg_name = CHECKED_SPRINTF("%s,%s",
                                  objp->sym->name,
                                  name_region_referent);
  if (full_reg_name==NULL)
    goto failure;

  /* Look up region or die */
  stp = retrieve_sym(full_reg_name, parse_state->vol->reg_sym_table);
  if (stp == NULL)
  {
    mdlerror_fmt(parse_state, "Undefined region: %s", full_reg_name);
    goto failure;
  }
  free(full_reg_name);
  full_reg_name = NULL;

  /* Store referent region */
  elmlp->special->referent = (struct region*) stp->value;
  if (elmlp->special->referent == rp_container)
  {
    mdlerror_fmt(parse_state, "Self-referential region include.  No paradoxes, please.");
    goto failure;
  }

  free(name_region_referent);
  return elmlp;

failure:
  free(name_region_referent);
  if (full_reg_name)
    free(full_reg_name);
  if (elmlp)
  {
    if (elmlp->special) free(elmlp->special);
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
struct element_list *
mdl_new_element_patch(struct mdlparse_vars *parse_state,
                      struct polygon_object *poly,
                      struct vector3 *llf,
                      struct vector3 *urb,
                      int exclude)
{
  if (parse_state->current_object->object_type != BOX_OBJ)
  {
    mdlerror(parse_state, "INCLUDE_PATCH and EXCLUDE_PATCH may only be used on a BOX object.");
    return NULL;
  }

  if (! parse_state->allow_patches)
  {
    mdlerror(parse_state, "Cannot create PATCH on a BOX outside of the original declaration.");
    return NULL;
  }

  struct element_list *elmlp = new_element_list(0, 0);
  if (elmlp == NULL)
    goto failure;

  /* Allocate special element description */
  elmlp->special = CHECKED_MALLOC_STRUCT(struct element_special, "region element");
  if (elmlp->special == NULL)
    goto failure;
  elmlp->special->referent = NULL;
  elmlp->special->exclude = (byte) exclude;

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
  if (refine_cuboid(parse_state, llf, urb, poly->sb, parse_state->vol->grid_density))
    return NULL;

  free(llf);
  free(urb);
  return elmlp;

failure:
  free(llf);
  free(urb);
  if (elmlp)
  {
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
int 
mdl_set_region_elements(struct mdlparse_vars *parse_state,
                        struct region *rgn,
                        struct element_list *elements,
                        int normalize_now)
{
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
struct sym_table *
mdl_new_rxn_pathname(struct mdlparse_vars *parse_state, char *name)
{
  if ((retrieve_sym(name, parse_state->vol->rxpn_sym_table)) != NULL)
  {
    mdlerror_fmt(parse_state, "Named reaction pathway already defined: %s", name);
    free(name);
    return NULL;
  }
  else if ((retrieve_sym(name, parse_state->vol->mol_sym_table)) != NULL)
  {
    mdlerror_fmt(parse_state, "Named reaction pathway already defined as a molecule: %s", name);
    free(name);
    return NULL;
  }

  struct sym_table *symp = store_sym(name,RXPN,parse_state->vol->rxpn_sym_table, NULL);
  if (symp == NULL)
  {
    mdlerror_fmt(parse_state, "Out of memory while creating reaction name: %s", name);
    free(name);
    return NULL;
  }
  free(name);
  return symp;
}



/**************************************************************************
 mdl_add_effector_to_region:
    Adds an effector (or list of effectors) to a region.  These effectors will
    be placed on the surface at initialization time.

 In: rgn:  the region
     lst:  a list of effectors to place
 Out: none.  list is merged into region
**************************************************************************/
void 
mdl_add_effector_to_region(struct region *rgn,
                           struct eff_dat_list *lst)
{
  lst->eff_tail->next = rgn->eff_dat_head;
  rgn->eff_dat_head = lst->eff_head;
}

/**************************************************************************
 mdl_set_region_surface_class:
    Set the surface class of this region, possibly inheriting the viz_value.

 In: parse_state: parser state
     rgn:  the region
     scsymp: symbol for the surface class
 Out: none.  region is updated.
**************************************************************************/
void 
mdl_set_region_surface_class(struct mdlparse_vars *parse_state,
                             struct region *rgn,
                             struct sym_table *scsymp)
{
  if (rgn->surf_class != NULL)
  {
      mdlerror(parse_state, "ATTENTION: region definition allows only one SURFACE_CLASS  statement.");
  }
  rgn->surf_class = (struct species *) scsymp->value;
  if (rgn->surf_class->region_viz_value > 0)
  {
    if (rgn->region_viz_value > 0)
      mdlerror(parse_state, "ATTENTION: region_viz_value defined both through SURFACE_CLASS and VIZ_VALUE statements; ignoring value from SURFACE_CLASS in this region.");
    else
      rgn->region_viz_value = rgn->surf_class->region_viz_value;
  }
}

/**************************************************************************
 mdl_set_region_region_viz_value:
    Set the VIZ_VALUE for this region.  (NOTE: This is different from the
    VIZ_STATE which can be set on the region.  It's unclear why we need both of
    them, but it seems that one is for newer output (DREAMM) and the other is
    for older output (DX).

 In: parse_state: parser state
     rgn:  the region
     viz_value: viz_value to set for this region
 Out: none.  region is updated.
**************************************************************************/
void 
mdl_set_region_region_viz_value(struct mdlparse_vars *parse_state,
                                struct region *rgn,
                                int viz_value)
{
  if (rgn->surf_class != NULL  && rgn->surf_class->region_viz_value > 0)
    mdlerror(parse_state, "ATTENTION: region_viz_value defined both through SURFACE_CLASS and VIZ_VALUE statements; ignoring value from SURFAE_CLASS in this region.");

  rgn->region_viz_value = viz_value;
}

/*************************************************************************
 * Reaction output
 *************************************************************************/

/**************************************************************************
 mdl_check_reaction_output_file:
    Check that the reaction output file is writable within the policy set by
    the user.  Creates and/or truncates the file to 0 bytes, as appropriate.
    Note that for SUBSTITUTE, the truncation is done later on, during
    initialization.

 In: parse_state: parser state
     os: output set containing file details
 Out: 0 if file preparation is successful, 1 if not.  The file named will be
      created and emptied or truncated as requested.
**************************************************************************/
static int 
mdl_check_reaction_output_file(struct mdlparse_vars *parse_state,
                               struct output_set *os)
{
  FILE *f;
  char *name;
  struct stat fs;
  int i;

  name = os->outfile_name;

  if (make_parent_dir(name))
  {
    mdlerror_fmt(parse_state, "Directory for %s does not exist and could not be created.",name);
    return 1;
  }

  switch (os->file_flags)
  {
    case FILE_OVERWRITE:
      f = fopen(name,"w");
      if (!f)
      {
    switch (errno)
    {
      case EACCES:
        mdlerror_fmt(parse_state,"Access to %s denied.",name);
        return 1;
      case ENOENT:
        mdlerror_fmt(parse_state,"Directory for %s does not exist",name);
        return 1;
      case EISDIR:
        mdlerror_fmt(parse_state,"%s already exists and is a directory",name);
        return 1;
      default:
        mdlerror_fmt(parse_state,"Unable to open %s for writing",name);
        return 1;
    }
      }
      fclose(f);
      break;
    case FILE_SUBSTITUTE:
      f = fopen(name,"a+");
      if (!f)
      {
    switch (errno)
    {
      case EACCES:
        mdlerror_fmt(parse_state,"Access to %s denied.",name);
        return 1;
      case ENOENT:
        mdlerror_fmt(parse_state,"Directory for %s does not exist",name);
        return 1;
      case EISDIR:
        mdlerror_fmt(parse_state,"%s already exists and is a directory",name);
        return 1;
      default:
        mdlerror_fmt(parse_state,"Unable to open %s for writing",name);
        return 1;
    }
      }
      i = fstat(fileno(f),&fs);
      if (!i && fs.st_size==0) os->file_flags = FILE_OVERWRITE;
      fclose(f);
      break;
    case FILE_APPEND:
    case FILE_APPEND_HEADER:
      f = fopen(name,"a");
      if (!f)
      {
    switch (errno)
    {
      case EACCES:
        mdlerror_fmt(parse_state,"Access to %s denied.",name);
        return 1;
      case ENOENT:
        mdlerror_fmt(parse_state,"Directory for %s does not exist",name);
        return 1;
      case EISDIR:
        mdlerror_fmt(parse_state,"%s already exists and is a directory",name);
        return 1;
      default:
        mdlerror_fmt(parse_state,"Unable to open %s for writing",name);
        return 1;
    }
      }
      i = fstat(fileno(f),&fs);
      if (!i && fs.st_size==0) os->file_flags = FILE_APPEND_HEADER;
      fclose(f);
      break;
    case FILE_CREATE:
      i = access(name,F_OK);
      if (!i)
      {
    i = stat(name,&fs);
    if (!i && fs.st_size>0)
    {
      mdlerror_fmt(parse_state,"Cannot create new file %s: it already exists",name);
      return 1;
    }
      }
      f = fopen(name,"w");
      if (f==NULL)
      {
    switch (errno)
    {
      case EEXIST:
        mdlerror_fmt(parse_state,"Cannot create %s because it already exists",name);
        return 1;
      case EACCES:
        mdlerror_fmt(parse_state,"Access to %s denied.",name);
        return 1;
      case ENOENT:
        mdlerror_fmt(parse_state,"Directory for %s does not exist",name);
        return 1;
      case EISDIR:
        mdlerror_fmt(parse_state,"%s already exists and is a directory",name);
        return 1;
      default:
        mdlerror_fmt(parse_state,"Unable to open %s for writing",name);
        return 1;
    }
      }
      fclose(f);
      break;

    default:
      UNHANDLED_CASE(os->file_flags);
      return 1;
  }
  return 0;
}

/**************************************************************************
 mdl_new_output_set:
    Create a new output set for reaction output.

 In: parse_state: parser state
     comment: header comment for output set
     exact_time: exact time column flag
 Out: new output set, or NULL if an error occurs
**************************************************************************/
struct output_set* 
mdl_new_output_set(struct mdlparse_vars *parse_state,
                   char *comment,
                   int exact_time)
{
  struct output_set *os;
  os = CHECKED_MALLOC_STRUCT(struct output_set, "reaction data output set");
  if (os == NULL)
    return NULL;

  os->outfile_name=NULL;
  os->file_flags=FILE_UNDEFINED;
  os->chunk_count=0;
  os->column_head=NULL;
  os->next = NULL;

  if (comment == NULL) os->header_comment = NULL;
  else if (comment[0] == '\0') os->header_comment = "";
  else
  {
    os->header_comment = mdl_strdup(comment);
    if (os->header_comment == NULL)
    {
      free(os);
      return NULL;
    }
  }

  os->exact_time_flag = exact_time;
  parse_state->count_flags = 0;
  return os;
}

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
struct output_set *
mdl_populate_output_set(struct mdlparse_vars *parse_state,
                        struct output_set *os,
                        struct output_column *col_head,
                        int file_flags,
                        char *outfile_name)
{
  struct output_column *oc = col_head;
  os->column_head = oc;
  os->file_flags = file_flags;
  os->outfile_name = outfile_name;
  no_printf("Counter output file set to %s\n", os->outfile_name);

  if ((parse_state->count_flags & (TRIGGER_PRESENT|COUNT_PRESENT)) == (TRIGGER_PRESENT|COUNT_PRESENT))
  {
    mdlerror(parse_state, "Cannot mix TRIGGER and COUNT statements.  Use separate files.");
    return NULL;
  }

  for (; oc != NULL; oc = oc->next)
    oc->set = os;

  if (mdl_check_reaction_output_file(parse_state, os))
    return NULL;

  return os;
}

/**************************************************************************
 mdl_new_output_block:
    Allocate a new reaction data output block, with a specified buffer size.

 In: buffersize: requested buffer size
 Out: output block, or NULL if an error occurs
**************************************************************************/
static struct output_block *
mdl_new_output_block(int buffersize)
{

  struct output_block *obp;
  obp = CHECKED_MALLOC_STRUCT(struct output_block,
                               "reaction data output block");
  if (obp==NULL) return NULL;

  obp->t = 0.0;
  obp->timer_type = OUTPUT_BY_STEP;
  obp->step_time = FOREVER;
  obp->time_list_head = NULL;
  obp->time_now = NULL;
  obp->buffersize = 0;
  obp->trig_bufsize = 0;
  obp->buf_index = 0;
  obp->data_set_head = NULL;

  /* COUNT buffer size might get modified later if there isn't that much to output */
  /* TRIGGER buffer size won't get modified since we can't know what to expect */
  obp->buffersize = buffersize;
  obp->trig_bufsize = obp->buffersize;

  obp->time_array = CHECKED_MALLOC_ARRAY(double,
                                          obp->buffersize,
                                          "reaction data output times array");
  if (obp->time_array == NULL)
  {
    free(obp);
    return NULL;
  }

  return obp;
}

/**************************************************************************
 mdl_output_block_finalize:
    Finalizes a reaction data output block, checking for errors, and allocating
    the output buffer.

 In: parse_state: parser state
     obp:  the output block to finalize
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_output_block_finalize(struct mdlparse_vars *parse_state, struct output_block *obp)
{
  struct output_set *os1;
  for (os1 = obp->data_set_head; os1 != NULL; os1 = os1->next)
  {
    /* Check for duplicated filenames */
    struct output_set *os2;
    for (os2 = os1->next; os2 != NULL; os2 = os2->next)
    {
      if (strcmp(os1->outfile_name, os2->outfile_name) == 0)
      {
        mdlerror_fmt(parse_state, "COUNT statements in the same reaction data output block should have unique output file names (\"%s\" appears more than once)", os1->outfile_name);
        return 1;
      }
    }

    /* Allocate buffers */
    struct output_column *oc;
    for (oc = os1->column_head; oc != NULL; oc = oc->next)
    {
      switch (oc->expr->expr_flags&OEXPR_TYPE_MASK)
      {
        case OEXPR_TYPE_INT:
          oc->data_type = COUNT_INT;
          oc->buffer = CHECKED_MALLOC_ARRAY(int,
                                             obp->buffersize,
                                             "reaction data output buffer");
          break;

        case OEXPR_TYPE_DBL:
          oc->data_type = COUNT_DBL;
          oc->buffer = CHECKED_MALLOC_ARRAY(double,
                                             obp->buffersize,
                                             "reaction data output buffer");
          break;

        case OEXPR_TYPE_TRIG:
          oc->data_type = COUNT_TRIG_STRUCT;
          oc->buffer = CHECKED_MALLOC_ARRAY(struct output_trigger_data,
                                             obp->trig_bufsize,
                                             "reaction data output buffer");
          break;

        default:
          mdlerror(parse_state, "Could not figure out what type of count data to store");
          return 1;
      }
      if (oc->buffer==NULL)
        return 1;
    }
  }

  return 0;
}

/**************************************************************************
 mdl_pick_buffer_size:
    Choose an appropriate output buffer size for our reaction output data,
    based on the total number of outputs expected and the requested buffer
    size.

 In: parse_state: parser state
     obp:  output block whose buffer_size to set
     n_output: maximum number of outputs expected
 Out: 0 on success, 1 on failure
**************************************************************************/
static long long 
mdl_pick_buffer_size(struct mdlparse_vars *parse_state,
                     struct output_block *obp,
                     long long n_output)
{
  if (parse_state->vol->chkpt_iterations)
    return min3ll(parse_state->vol->chkpt_iterations-parse_state->vol->start_time+1, n_output, obp->buffersize);
  else
    return min3ll(parse_state->vol->iterations-parse_state->vol->start_time+1, n_output, obp->buffersize);
}

/**************************************************************************
 mdl_set_reaction_output_timer_step:
    Set the output timer for reaction data output to a time step.

 In: parse_state: parser state
     obp:  output block whose timer is to be set
     step: time step interval
 Out: output timer is updated
**************************************************************************/
static void 
mdl_set_reaction_output_timer_step(struct mdlparse_vars *parse_state,
                                   struct output_block *obp,
                                   double step)
{
  long long output_freq;
  obp->timer_type = OUTPUT_BY_STEP;

  obp->step_time = step;
  output_freq = (long long) (obp->step_time / parse_state->vol->time_unit);

  /* Clip the step time to a good range */
  if (output_freq > parse_state->vol->iterations  &&  output_freq > 1)
  {
    output_freq = (parse_state->vol->iterations > 1) ? parse_state->vol->iterations : 1;
    obp->step_time = output_freq*parse_state->vol->time_unit;
    if (parse_state->vol->notify->invalid_output_step_time != WARN_COPE)
      mdl_warning(parse_state, "Output step time too long\n\tSetting output step time to %g microseconds\n", obp->step_time*1.0e6);
  }
  else if (output_freq < 1)
  {
    output_freq = 1;
    obp->step_time = output_freq*parse_state->vol->time_unit;
    if (parse_state->vol->notify->invalid_output_step_time != WARN_COPE)
      mdl_warning(parse_state, "Output step time too short\n\tSetting output step time to %g microseconds\n",obp->step_time*1.0e-6);
  }

  /* Pick a good buffer size */
  long long n_output;
  if (parse_state->vol->chkpt_iterations)
    n_output = (long long)(parse_state->vol->chkpt_iterations / output_freq + 1);
  else
    n_output = (long long)(parse_state->vol->iterations / output_freq + 1);
  obp->buffersize = mdl_pick_buffer_size(parse_state, obp, n_output);

  no_printf("Default output step time definition:\n");
  no_printf("  output step time = %g\n", obp->step_time);
  no_printf("  output buffersize = %u\n", obp->buffersize);
}

/**************************************************************************
 mdl_set_reaction_output_timer_iterations:
    Set the output timer for reaction data output to a list of iterations.

 In: parse_state: parser state
     obp:  output block whose timer is to be set
     step_values: list of iterations
 Out: 0 on success, 1 on failure; output timer is updated
**************************************************************************/
static int 
mdl_set_reaction_output_timer_iterations(struct mdlparse_vars *parse_state,
                                         struct output_block *obp,
                                         struct num_expr_list_head *step_values)
{
  obp->timer_type = OUTPUT_BY_ITERATION_LIST;
  obp->buffersize = mdl_pick_buffer_size(parse_state, obp, step_values->value_count);
  if (step_values->shared)
  {
    obp->time_list_head = mdl_copysort_numeric_list(parse_state, step_values->value_head);
    if (obp->time_list_head == NULL)
      return 1;
  }
  else
  {
    mdl_sort_numeric_list(step_values->value_head);
    obp->time_list_head = step_values->value_head;
  }
  obp->time_now = NULL;
  return 0;
}

/**************************************************************************
 mdl_set_reaction_output_timer_times:
    Set the output timer for reaction data output to a list of times.

 In: parse_state: parser state
     obp:  output block whose timer is to be set
     nstep: number of times
     step_values: list of times
 Out: output timer is updated
**************************************************************************/
static int 
mdl_set_reaction_output_timer_times(struct mdlparse_vars *parse_state,
                                    struct output_block *obp,
                                    struct num_expr_list_head *step_values)
{
  obp->timer_type = OUTPUT_BY_TIME_LIST;
  obp->buffersize = mdl_pick_buffer_size(parse_state, obp, step_values->value_count);
  if (step_values->shared)
  {
    obp->time_list_head = mdl_copysort_numeric_list(parse_state, step_values->value_head);
    if (obp->time_list_head == NULL)
      return 1;
  }
  else
  {
    mdl_sort_numeric_list(step_values->value_head);
    obp->time_list_head = step_values->value_head;
  }
  obp->time_now = NULL;
  return 0;
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
int 
mdl_add_reaction_output_block_to_world(struct mdlparse_vars *parse_state,
                                       int buffer_size,
                                       struct output_times_inlist *otimes,
                                       struct output_set_list *osets)
{
  struct output_block *obp;
  struct output_set *os;
  if ((obp = mdl_new_output_block(buffer_size)) == NULL)
    return 1;

  if (otimes->type == OUTPUT_BY_STEP)
    mdl_set_reaction_output_timer_step(parse_state, obp, otimes->step);
  else if (otimes->type == OUTPUT_BY_ITERATION_LIST)
  {
    if (mdl_set_reaction_output_timer_iterations(parse_state, obp, & otimes->values))
      return 1;
  }
  else if (otimes->type == OUTPUT_BY_TIME_LIST)
  {
    if (mdl_set_reaction_output_timer_times(parse_state, obp, & otimes->values))
      return 1;
  }
  else
  {
    mdlerror_fmt(parse_state, "Internal error: Invalid output timer def (%d)", otimes->type);
    return 1;
  }
  obp->data_set_head = osets->set_head;
  for (os = obp->data_set_head; os != NULL; os = os->next)
    os->block = obp;
  if (mdl_output_block_finalize(parse_state, obp))
    return 1;
  obp->next = parse_state->vol->output_block_head;
  parse_state->vol->output_block_head = obp;
  return 0;
}

/**************************************************************************
 mdl_new_output_column:
    Create a new output column for an output set.

 In: Nothing
 Out: output column, or NULL if allocation fails
**************************************************************************/
static struct output_column*
mdl_new_output_column()
{

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
 mdl_join_oexpr_tree:
    Joins two subtrees into a reaction data output expression tree, with a
    specified operation.

 In: parse_state: parser state
     left: left subtree
     right: right subtree
     oper: specified operation
 Out: joined output expression, or NULL if an error occurs
**************************************************************************/
struct output_expression*
mdl_join_oexpr_tree(struct mdlparse_vars *parse_state,
                    struct output_expression *left,
                    struct output_expression *right,
                    char oper)
{
  struct output_expression *joined;
  struct output_expression *leaf,*new_oe,*up;
  int first_leaf=1;

  joined=NULL;
  if (left->oper==',' && (right==NULL || right->oper==','))
  {
    mdlerror(parse_state, "Can't do math on multiple wildcard expressions");
    return NULL;
  }

  if (left->oper!=',' && (right==NULL||right->oper!=','))
  {
    joined = new_output_expr(parse_state->vol->oexpr_mem);
    if (joined==NULL) return NULL;

    joined->left=(void*)left;
    joined->right=(void*)right;
    joined->oper=oper;
    left->up=joined;
    if (right!=NULL) right->up=joined;

    learn_oexpr_flags(joined);
    if (joined->expr_flags&OEXPR_TYPE_CONST) eval_oexpr_tree(joined,0);

    return joined;
  }
  else if (left->oper==',')
  {
    for (leaf=first_oexpr_tree(left) ; leaf!=NULL ; leaf=next_oexpr_tree(leaf))
    {
      if (first_leaf)
      {
        new_oe=right;
        first_leaf=0;
      }
      else if (right!=NULL)
      {
        new_oe=dupl_oexpr_tree(right,parse_state->vol->oexpr_mem);
        if (new_oe==NULL) return NULL;
      }
      else new_oe=NULL;

      up=leaf->up;
      joined=mdl_join_oexpr_tree(parse_state, leaf,new_oe,oper);
      joined->up=up;
      if (joined==NULL) return NULL;
      if (leaf==up->left) up->left=joined;
      else up->right=joined;
      if (joined->expr_flags&OEXPR_TYPE_CONST) eval_oexpr_tree(joined,0);
      learn_oexpr_flags(up);
      leaf=joined;
    }
    return left;
  }
  else /* right->oper==',' */
  {
    for (leaf=first_oexpr_tree(right) ; leaf!=NULL ; leaf=next_oexpr_tree(leaf))
    {
      if (first_leaf)
      {
        new_oe=left;
        first_leaf=0;
      }
      else
      {
        new_oe=dupl_oexpr_tree(left,parse_state->vol->oexpr_mem);
        if (new_oe==NULL) return NULL;
      }
      up=leaf->up;
      joined=mdl_join_oexpr_tree(parse_state, new_oe,leaf,oper);
      joined->up=up;
      if (joined==NULL) return NULL;
      if (leaf==up->left) up->left=joined;
      else up->right=joined;
      if (joined->expr_flags&OEXPR_TYPE_CONST) eval_oexpr_tree(joined,0);
      learn_oexpr_flags(up);
      leaf=joined;
    }

    return right;
  }

  return NULL;  /* Should never get here */
}

/**************************************************************************
 mdl_sum_expression:
    Convert an output expression tree into a summation.

 In: expr: expression to convert
 Out: modified expression
**************************************************************************/
struct output_expression *
mdl_sum_oexpr(struct output_expression *expr)
{
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
mdl_new_oexpr_constant(struct mdlparse_vars *parse_state,
                       double value)
{
  struct output_expression *oe = new_output_expr(parse_state->vol->oexpr_mem);
  if (oe == NULL)
  {
    mdlerror(parse_state, "Out of memory creating output expression");
    return NULL;
  }
  oe->expr_flags = OEXPR_TYPE_DBL | OEXPR_TYPE_CONST;
  oe->value = value;
  oe->oper = '=';
  return oe;
}

/*************************************************************************
 mdl_new_output_request:
    Create a new output request.

 In:  parse_state: parser state
      target: what are we counting
      orientation: how is it oriented?
      location: where are we counting?
      report_flags: what type of events are we counting?
 Out: output request item, or NULL if an error occurred
*************************************************************************/
static struct output_request *
mdl_new_output_request(struct mdlparse_vars *parse_state,
                       struct sym_table *target,
                       short orientation,
                       struct sym_table *location,
                       int report_flags)
{
  struct output_request *orq;
  struct output_expression *oe;

  orq = CHECKED_MEM_GET(parse_state->vol->outp_request_mem, "count request");
  if (orq == NULL)
    return NULL;

  oe = new_output_expr(parse_state->vol->oexpr_mem);
  if (oe == NULL)
  {
    mem_put(parse_state->vol->outp_request_mem, orq);
    mcell_allocfailed("Failed to allocate a count expression.");
    return NULL;
  }
  orq->next=NULL;
  orq->requester=oe;
  orq->count_target = target;
  orq->count_orientation = orientation;
  orq->count_location = location;
  orq->report_type = report_flags;

  oe->left = orq;
  oe->oper = '#';
  oe->expr_flags = OEXPR_LEFT_REQUEST;
  if (orq->report_type&REPORT_TRIGGER)
    oe->expr_flags|=OEXPR_TYPE_TRIG;
  else if ((orq->report_type&REPORT_TYPE_MASK)!=REPORT_CONTENTS)
    oe->expr_flags|=OEXPR_TYPE_DBL;
  else
    oe->expr_flags|=OEXPR_TYPE_INT;
  return orq;
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
struct output_expression *
mdl_count_syntax_1(struct mdlparse_vars *parse_state,
                   struct sym_table *what,
                   struct sym_table *where,
                   int hit_spec,
                   int count_flags)
{
  byte report_flags = 0;
  struct output_request *orq;
  if (where == NULL)
  {
    report_flags = REPORT_WORLD;
    if (hit_spec != REPORT_NOTHING)
    {
      mdlerror(parse_state, "Invalid combination of WORLD with other counting options");
      return NULL;
    }
    else if (count_flags & TRIGGER_PRESENT)
    {
      mdlerror(parse_state, "Invalid combination of WORLD with TRIGGER option");
      return NULL;
    }
  }
  else report_flags = 0;

  if (count_flags & TRIGGER_PRESENT) report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED) report_flags |= REPORT_ENCLOSED;

  if (what->sym_type == MOL)
  {
    if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING) report_flags |= REPORT_CONTENTS;
    else report_flags |= (hit_spec & REPORT_TYPE_MASK);
  }
  else
  {
    report_flags |= REPORT_RXNS;
    if ((hit_spec & REPORT_TYPE_MASK) != REPORT_NOTHING)
    {
      mdlerror_fmt(parse_state, "Invalid counting options used with reaction pathway %s", what->name);
      return NULL;
    }
  }

  if ((orq = mdl_new_output_request(parse_state, what, ORIENT_NOT_SET, where, report_flags)) == NULL)
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
struct output_expression *
mdl_count_syntax_2(struct mdlparse_vars *parse_state,
                   struct sym_table *mol_type,
                   short orient,
                   struct sym_table *where,
                   int hit_spec,
                   int count_flags)
{
  byte report_flags = 0;
  struct output_request *orq;
  short orientation;
  if (where == NULL)
  {
    mdlerror(parse_state, "Counting of an oriented molecule in the WORLD is not implemented.\nAn oriented molecule may only be counted in a regions.");
    return NULL;
  }
  else report_flags = 0;

  if (count_flags & TRIGGER_PRESENT) report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED) report_flags |= REPORT_ENCLOSED;

  if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING) report_flags |= REPORT_CONTENTS;
  else report_flags |= (hit_spec & REPORT_TYPE_MASK);

  /* Grab orientation and reset orientation state in parser */
  if (orient < 0)
    orientation = -1;
  else if (orient > 0)
    orientation = 1;
  else
    orientation = 0;

  if ((orq = mdl_new_output_request(parse_state, mol_type, orientation, where, report_flags)) == NULL)
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
static int 
mdl_get_orientation_from_string(char *mol_string)
{
  short orientation = 0;
  int pos = strlen(mol_string) - 1;

  if (mol_string[pos] == ';')
  {
    mol_string[pos] = '\0';
    return 0;
  }

  /* Peel off any apostrophes or commas at the end of the string */
  while (pos >= 0)
  {
    switch (mol_string[pos])
    {
      case '\'': ++ orientation; break;
      case ',':  -- orientation; break;
      default:
        mol_string[pos+1] = '\0';
        pos = 0;
        break;
    }
    -- pos;
  }

  if (orientation < 0) return -1;
  else if (orientation > 0) return 1;
  else return 0;
}

/**************************************************************************
 mdl_string_has_orientation:
    Check if the string looks like a molecule name+orientation string.
    Otherwise, it will be considered a wildcard.  This only supports the "tick"
    notation.

 In: mol_string: string to parse
 Out: 1 if the string has orientation, 0 if it does not
**************************************************************************/
static int 
mdl_string_has_orientation(char const *mol_string)
{
  char last_char = mol_string[ strlen(mol_string) - 1 ];
  return (last_char == '\''  ||  last_char == ','  ||  last_char == ';');
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
static struct output_expression *
mdl_new_output_requests_from_list(struct mdlparse_vars *parse_state,
                                  struct sym_table_list *targets,
                                  struct sym_table *location,
                                  int report_flags,
                                  int hit_spec)
{
  struct output_expression *oe_head=NULL, *oe_tail=NULL;
  struct output_request *or_head=NULL, *or_tail=NULL;
  int report_type;

  assert(targets != NULL);
  for (; targets != NULL; targets = targets->next)
  {
    if (targets->node->sym_type==MOL)
    {
      if ((hit_spec&REPORT_TYPE_MASK)==REPORT_NOTHING) report_type = REPORT_CONTENTS;
      else report_type=(hit_spec&REPORT_TYPE_MASK);
    }
    else
    {
      report_type = REPORT_RXNS;
      if ((hit_spec&REPORT_TYPE_MASK)!=REPORT_NOTHING)
      {
        mdlerror_fmt(parse_state, "Invalid counting options used with reaction pathway %s",targets->node->name);
        return NULL;
      }
    }

    struct output_request *orq = mdl_new_output_request(parse_state,
                                                        targets->node,
                                                        ORIENT_NOT_SET,
                                                        location,
                                                        report_type|report_flags);
    if (orq == NULL)
      return NULL;
    struct output_expression *oe = orq->requester;

    if (oe_tail==NULL)
    {
      oe_head=oe_tail=oe;
      or_head=or_tail=orq;
    }
    else
    {
      or_tail->next=orq;
      or_tail=orq;

      struct output_expression *oet = new_output_expr(parse_state->vol->oexpr_mem);
      if (oet==NULL)
      {
        mdlerror(parse_state, "Out of memory storing requested counts");
        return NULL;
      }
      if (oe_tail->up==NULL)
      {
        oe_head=oet;
      }
      else
      {
        oet->up = oe_tail->up;
        if (oet->up->left==oe_tail) oet->up->left=oet;
        else oet->up->right=oet;
      }
      oet->left=oe_tail;
      oe_tail->up=oet;
      oet->right=oe;
      oe->up=oet;
      oet->oper=',';
      oet->expr_flags = OEXPR_LEFT_OEXPR|OEXPR_RIGHT_OEXPR|(oe->expr_flags&OEXPR_TYPE_MASK);
      oe_tail=oe;
    }
  }

  or_tail->next=parse_state->vol->output_request_head;
  parse_state->vol->output_request_head=or_head;
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
                                    char const *wildcard)
{
  struct sym_table_list *symbols = NULL, *stl;
  for (int i = 0; i < parse_state->vol->mol_sym_table->n_bins; i++)
  {
    for (struct sym_table *sym_t = parse_state->vol->mol_sym_table->entries[i];
        sym_t != NULL;
        sym_t = sym_t->next)
    {
      if (is_wildcard_match((char *)wildcard, sym_t->name))
      {
        stl = (struct sym_table_list *) CHECKED_MEM_GET(parse_state->sym_list_mem,
                                                         "list of named reactions and molecules for counting");
        if (stl == NULL)
        {
          if (symbols) mem_put_list(parse_state->sym_list_mem, symbols);
          return NULL;
        }

        stl->node = sym_t;
        stl->next = symbols;
        symbols = stl;
      }
    }
  }
  for (int i = 0; i < parse_state->vol->rxpn_sym_table->n_bins; i++)
  {
    for (struct sym_table *sym_t = parse_state->vol->rxpn_sym_table->entries[i];
        sym_t != NULL;
        sym_t = sym_t->next)
    {
      if (is_wildcard_match((char *)wildcard, sym_t->name))
      {
        stl = (struct sym_table_list *) CHECKED_MEM_GET(parse_state->sym_list_mem,
                                                         "list of named reactions and molecules for counting");
        if (stl == NULL)
        {
          if (symbols) mem_put_list(parse_state->sym_list_mem, symbols);
          return NULL;
        }

        stl->node = sym_t;
        stl->next = symbols;
        symbols = stl;
      }
    }
  }

  if (symbols == NULL)
    mdlerror_fmt(parse_state, "No molecules or named reactions found matching wildcard \"%s\"", wildcard);

  return symbols;
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
struct output_expression *
mdl_count_syntax_3(struct mdlparse_vars *parse_state,
                   char *what,
                   struct sym_table *where,
                   int hit_spec,
                   int count_flags)
{
  struct output_expression *oe;
  char *what_to_count;
  byte report_flags = 0;
  if ((what_to_count = mdl_strip_quotes(what)) == NULL)
    return NULL;

  if (count_flags & TRIGGER_PRESENT) report_flags |= REPORT_TRIGGER;
  if (hit_spec & REPORT_ENCLOSED) report_flags |= REPORT_ENCLOSED;

  /* Oriented molecule specified inside a string */
  if (mdl_string_has_orientation(what_to_count))
  {
    struct output_request *orq;
    struct sym_table *sp;
    short orientation;

    if (where == NULL)
    {
      mdlerror(parse_state, "Counting of an oriented molecule in the WORLD is not implemented.\nAn oriented molecule may only be counted in a regions.");
      return NULL;
    }

    orientation = mdl_get_orientation_from_string(what_to_count);
    if ((sp = mdl_existing_molecule(parse_state, what_to_count)) == NULL)
      return NULL;

    if ((hit_spec & REPORT_TYPE_MASK) == REPORT_NOTHING) report_flags |= REPORT_CONTENTS;
    else report_flags |= (hit_spec & REPORT_TYPE_MASK);

    if ((orq = mdl_new_output_request(parse_state, sp, orientation, where, report_flags)) == NULL)
      return NULL;
    orq->next = parse_state->vol->output_request_head;
    parse_state->vol->output_request_head = orq;
    oe = orq->requester;
  }

  /* Wildcard specified inside a string */
  else
  {
    struct sym_table_list *stl = mdl_find_rxpns_and_mols_by_wildcard(parse_state, what_to_count);

    if (stl == NULL)
    {
      mdlerror(parse_state, "Wildcard matching found no matches for count output.");
      return NULL;
    }

    if (where == NULL)
    {
      report_flags |= REPORT_WORLD;
      if (hit_spec != REPORT_NOTHING)
      {
        mdlerror(parse_state, "Invalid combination of WORLD with other counting options");
        return NULL;
      }
      else if (count_flags & TRIGGER_PRESENT)
      {
        mdlerror(parse_state, "Invalid combination of WORLD with TRIGGER option");
        return NULL;
      }
    }

    if ((oe = mdl_new_output_requests_from_list(parse_state, stl, where, report_flags, hit_spec)) == NULL)
      return NULL;

    /* free allocated memory */
    mem_put_list(parse_state->sym_list_mem, stl);
  }

  return oe;
}

/*************************************************************************
 macro_new_complex_count:
    Allocate a new complex count structure for subunit counting and an output
    expression to reference it.

 In:  parse_state - the parser state for error logging
      macromol - the macromolecule species
      master_orientation - the relevant orientation of the complex itself
      subunit - the type of the subunit of interest
      subunit_orientation - the orientation of the subunit of interest
      relation_states - a linked list of conditions which must be met in order
                        for the count to match
      location - the location for which to do the counting, or NULL for
                        "WORLD".  This location may be a region or an object.
 Out: A freshly allocated relation state struct, or NULL if allocation fails.
*************************************************************************/
static struct output_expression *
macro_new_complex_count(struct mdlparse_vars *parse_state,
                        struct complex_species *macromol,
                        short master_orientation,
                        struct species *subunit,
                        short subunit_orientation,
                        struct macro_relation_state *relation_states,
                        struct sym_table *location)
{
  struct macro_count_request *mcr;
  mcr = CHECKED_MALLOC_STRUCT(struct macro_count_request,
                               "macromolecule count request");
  if (mcr == NULL)
    return NULL;

  mcr->next = NULL;
  mcr->paired_expression = NULL;
  mcr->the_complex = macromol;
  mcr->master_orientation = master_orientation;
  mcr->subunit_state = subunit;
  mcr->subunit_orientation = subunit_orientation;
  mcr->relation_states = relation_states;
  mcr->location = location;

  mcr->next = parse_state->vol->macro_count_request_head;
  parse_state->vol->macro_count_request_head = mcr;

  mcr->paired_expression = new_output_expr(parse_state->vol->oexpr_mem);
  mcr->paired_expression->left = mcr;
  mcr->paired_expression->oper = '@';
  mcr->paired_expression->expr_flags = OEXPR_LEFT_MACROREQUEST | OEXPR_TYPE_INT;

  return mcr->paired_expression;
}

/**************************************************************************
 mdl_count_syntax_macromol_subunit:
    Generate a reaction data output expression from the macromolecule "subunit"
    syntax variant.

      Example:
        COUNT(SUBUNIT { CamKII' : CamSub_U' [ dimer_partner == CamSub_U' ] }, WORLD)

 In: parse_state: parser state
     macromol: the species of the macromolecule whose subunits we're counting
     master_orientation: optional orientation of the complex as a whole
     subunit: type and orientation of the reference subunit to count
     relation_states: rules to restrict matches by related subunit states
     location: region or object where we are counting, or NULL for WORLD
 Out: output exprsession, or NULL if an error occurs
**************************************************************************/
struct output_expression *
mdl_count_syntax_macromol_subunit(struct mdlparse_vars *parse_state,
                                  struct complex_species *macromol,
                                  struct species_opt_orient *master_orientation,
                                  struct species_opt_orient *subunit,
                                  struct macro_relation_state *relation_states,
                                  struct sym_table *location)
{
  if ((macromol->base.flags & NOT_FREE) == 0)
  {
    if (master_orientation->orient_set)
    {
      if (parse_state->vol->notify->useless_vol_orient==WARN_ERROR)
      {
        mdlerror_fmt(parse_state, "Error: orientation specified for volume complex '%s' in count statement", macromol->base.sym->name);
        return NULL;
      }
      else if (parse_state->vol->notify->useless_vol_orient==WARN_WARN)
        mdlerror_fmt(parse_state, "Warning: orientation specified for volume complex '%s' in count statement", macromol->base.sym->name);
    }

    if (subunit->orient_set)
    {
      if (parse_state->vol->notify->useless_vol_orient==WARN_ERROR)
      {
        mdlerror_fmt(parse_state, "Error: orientation specified for subunit of volume complex '%s' in count statement", macromol->base.sym->name);
        return NULL;
      }
      else if (parse_state->vol->notify->useless_vol_orient==WARN_WARN)
        mdlerror_fmt(parse_state, "Warning: orientation specified for subunit of volume complex '%s' in count statement", macromol->base.sym->name);
    }

    master_orientation->orient = 0;
    subunit->orient = 0;
  }
  else
  {
    if (! master_orientation->orient_set)
      master_orientation->orient = 0;
    if (! subunit->orient_set)
      subunit->orient = 0;
  }

  /* Check that no relation features more than once in the relation states */
  struct macro_relation_state *states1;
  for (states1 = relation_states; states1 != NULL; states1 = states1->next)
  {
    struct macro_relation_state *states2;
    for (states2 = states1->next; states2 != NULL; states2 = states2->next)
    {
      if (states1->relation == states2->relation)
      {
        mdlerror_fmt(parse_state,
                     "Error: The SUBUNIT count statement for the complex '%s' includes multiple references to the relation '%s'",
                     macromol->base.sym->name,
                     macromol->relations[states1->relation].name);
        /* XXX: Free relation states */
        return NULL;
      }
    }
  }

  /* Create macro count request and associated expression */
  parse_state->vol->place_waypoints_flag = 1;
  return macro_new_complex_count(parse_state, macromol, master_orientation->orient, (struct species *) subunit->mol_type->value, subunit->orient, relation_states, location);
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
int 
mdl_single_count_expr(struct mdlparse_vars *parse_state,
                      struct output_column_list *list,
                      struct output_expression *expr,
                      char *custom_header)
{
  struct output_expression *oer,*oe;
  struct output_column *oc;

  list->column_head = NULL;
  list->column_tail = NULL;

  oer = expr; /* Root of count expression */

  if (oer->oper == ',' && custom_header != NULL)
  {
    mdlerror(parse_state, "Cannot use custom column headers with wildcard expansion");
    return 1;
  }

  if (custom_header != NULL) oer->title = custom_header;

  /* If we have a list of results, go through to build column stack */
  for (oe = first_oexpr_tree(oer); oe != NULL; oe = next_oexpr_tree(oe))
  {
    if ((oc = mdl_new_output_column()) == NULL)
      return 1;

    if (! list->column_head)
      list->column_head = list->column_tail = oc;
    else
      list->column_tail = list->column_tail->next = oc;

    oc->expr = oe;
    set_oexpr_column(oe, oc);
  }

  return 0;
}

/*************************************************************************
 * VIZ output
 *************************************************************************/

static struct viz_child *mdl_get_viz_child(struct viz_output_block *vizblk,
                                           struct object *objp);

/**************************************************************************
 mdl_new_viz_output_block:
    Build a new VIZ output block, containing parameters for an output set for
    visualization.

 In: parse_state: parser state
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_new_viz_output_block(struct mdlparse_vars *parse_state)
{
  struct viz_output_block *vizblk = CHECKED_MALLOC_STRUCT(struct viz_output_block,
                                                          "visualization data output parameters");
  if (vizblk == NULL)
    return 1;

  vizblk->frame_data_head = NULL;
  memset(&vizblk->viz_state_info, 0, sizeof(vizblk->viz_state_info));
  vizblk->viz_mode = -1;
  vizblk->molecule_prefix_name = NULL;
  vizblk->file_prefix_name = NULL;
  vizblk->viz_output_flag = 0;
  vizblk->species_viz_states = NULL;

  vizblk->dreamm_object_info = NULL;
  vizblk->dreamm_objects = NULL;
  vizblk->n_dreamm_objects = 0;

  vizblk->dx_obj_head = NULL;
  vizblk->viz_children = init_symtab(1024);
  if (pointer_hash_init(&vizblk->parser_species_viz_states, 32))
    mcell_allocfailed("Failed to initialize viz species states table.");

  vizblk->next = parse_state->vol->viz_blocks;
  parse_state->vol->viz_blocks = vizblk;
  return 0;
}

/**************************************************************************
 mdl_finish_viz_output_block:
    Finalize a new VIZ output block, ensuring all required parameters were set,
    and doing any postprocessing necessary for runtime.

 In: parse_state: parser state
     vizblk: the viz block to check
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_finish_viz_output_block(struct mdlparse_vars *parse_state,
                            struct viz_output_block *vizblk)
{
  if (vizblk->viz_mode == DREAMM_V3_MODE  ||
      vizblk->viz_mode == DREAMM_V3_GROUPED_MODE)
  {
    if (vizblk->file_prefix_name == NULL)
    {
      mdlerror(parse_state, "DREAMM output mode requested, but the required keyword FILENAME has not been supplied.");
      return 1;
    }
  }

  return 0;
}

/**************************************************************************
 mdl_require_old_style_viz:
    Require the mode of the specified VIZ output block to be pre-MCell 3 (i.e.
    neither of the DREAMM_V3 modes.)  It's impossible to create a valid DREAMM
    output using the old-style notation, so this makes the impossibility more
    explicit for user-friendliness.

 In: parse_state: parser state
     mode: the requested mode
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_require_old_style_viz(struct mdlparse_vars *parse_state, int mode)
{
  if (mode == DREAMM_V3_MODE  ||
      mode == DREAMM_V3_GROUPED_MODE)
  {
    mdlerror(parse_state, "DREAMM output modes must use the new VIZ_OUTPUT notation.\n  Please change this VIZ_DATA_OUTPUT block to a VIZ_OUTPUT block.");
    return 1;
  }

  return 0;
}

/**************************************************************************
 mdl_set_viz_mode:
    Set the mode for a new VIZ output block.

 In: vizblk: the viz block to check
     mode: the mode to set
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_mode(struct viz_output_block *vizblk, int mode)
{

  vizblk->viz_mode = mode;
  return 0;
}

/**************************************************************************
 mdl_set_viz_mesh_format:
    Set the mesh format for a new VIZ output block.

 In: parse_state: parser state
     vizblk: the viz block to check
     format: the format to set
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_mesh_format(struct mdlparse_vars *parse_state,
                        struct viz_output_block *vizblk,
                        int format)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->viz_mode != DREAMM_V3_MODE)
  {
    mdlerror_fmt(parse_state, "VIZ_MESH_FORMAT command is allowed only in DREAMM_V3 mode.");
    return 1;
  }
  vizblk->viz_output_flag |= format;
  return 0;
}

/**************************************************************************
 mdl_set_viz_molecule_format:
    Set the molecule format for a new VIZ output block.

 In: parse_state: parser state
     vizblk: the viz block to check
     format: the format to set
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_molecule_format(struct mdlparse_vars *parse_state,
                            struct viz_output_block *vizblk,
                            int format)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->viz_mode != DREAMM_V3_MODE)
  {
    mdlerror_fmt(parse_state, "VIZ_MOLECULE_FORMAT command is allowed only in DREAMM_V3 mode.");
    return 1;
  }
  vizblk->viz_output_flag |= format;
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
int 
mdl_set_viz_filename_prefix(struct mdlparse_vars *parse_state,
                            struct viz_output_block *vizblk,
                            char *filename)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->file_prefix_name != NULL)
  {
    mdlerror_fmt(parse_state, "FILENAME may only appear once per output block.");
    free(filename);
    return 1;
  }

  if (vizblk->viz_mode == ASCII_MODE)
  {
    if (vizblk->molecule_prefix_name != NULL)
    {
      mdlerror_fmt(parse_state, "In non-DREAMM/DX modes, MOLECULE_FILE_PREFIX and FILENAME are aliases, and may not both be specified.");
      free(filename);
      return 1;
    }
  }

  vizblk->file_prefix_name = filename;
  if (vizblk->molecule_prefix_name == NULL)
    vizblk->molecule_prefix_name = filename;

  return 0;
}

/**************************************************************************
 mdl_set_viz_molecule_filename_prefix:
    Set the molecule filename prefix for an old-style VIZ output block.

 In: parse_state: parser state
     vizblk: the viz block to check
     filename: the filename
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_molecule_filename_prefix(struct mdlparse_vars *parse_state,
                                     struct viz_output_block *vizblk,
                                     char *filename)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->molecule_prefix_name != NULL  &&
      vizblk->molecule_prefix_name != vizblk->file_prefix_name)
  {
    mdlerror_fmt(parse_state, "MOLECULE_FILE_PREFIX may only appear once per output block.");
    free(filename);
    return 1;
  }

  if (vizblk->viz_mode == DREAMM_V3_MODE  ||
      vizblk->viz_mode == DREAMM_V3_GROUPED_MODE)
  {
    mdlerror(parse_state, "MOLECULE_FILE_PREFIX canot be used with the DREAMM viz output modes.");
    return 1;
  }

  if (vizblk->viz_mode != DX_MODE)
  {
    if (vizblk->file_prefix_name != NULL)
    {
      mdlerror_fmt(parse_state, "In non-DREAMM/DX output modes, MOLECULE_FILE_PREFIX and FILENAME are aliases, and may not both be specified.");
      free(filename);
      return 1;
    }
    vizblk->file_prefix_name = filename;
  }

  vizblk->molecule_prefix_name = filename;

  return 0;
}

static struct viz_child *get_viz_children(struct mdlparse_vars *parse_state,
                                          struct viz_output_block *vizblk,
                                          struct viz_child *vcp_parent,
                                          struct object *obj)
{
  struct viz_child *vcp = mdl_get_viz_child(vizblk, obj);
  if (vcp == NULL)
    return NULL;

  if (vcp_parent != NULL  &&  vcp->parent == NULL)
  {
    vcp->next = vcp_parent->children;
    vcp_parent->children = vcp;
    vcp->parent = vcp_parent;
  }

  if (obj->object_type == META_OBJ)
  {
    for (struct object *child_objp = obj->first_child;
         child_objp != NULL;
         child_objp = child_objp->next)
      get_viz_children(parse_state, vizblk, vcp, child_objp);
  }

  return vcp;
}

static struct viz_dx_obj *mdl_get_viz_obj(struct mdlparse_vars *parse_state,
                                          struct viz_output_block *vizblk,
                                          struct sym_table *sym)
{
  for (struct viz_dx_obj *viz = vizblk->dx_obj_head;
       viz != NULL;
       viz = viz->next)
  {
    if (viz->obj == sym->value)
      return viz;
  }

  struct viz_dx_obj *viz = CHECKED_MALLOC_STRUCT(struct viz_dx_obj, "viz object");
  if (viz == NULL)
    return NULL;

  viz->name = NULL;
  viz->full_name = mdl_strdup(sym->name);
  if (viz->full_name == NULL)
  {
    free(viz);
    return NULL;
  }
  viz->obj = (struct object *) sym->value;
  viz->viz_child_head = get_viz_children(parse_state, vizblk, NULL, viz->obj);
  viz->parent = vizblk;
  viz->actual_objects = NULL;
  viz->n_actual_objects = 0;
  viz->next = vizblk->dx_obj_head;
  vizblk->dx_obj_head = viz;
  return viz;
}

/**************************************************************************
 mdl_set_viz_object_filename_prefix:
    Set the object filename prefix for an old-style VIZ output block.

 In: parse_state: parser state
     vizblk: the viz block to check
     sym: the object
     filename: the filename
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_object_filename_prefix(struct mdlparse_vars *parse_state,
                                   struct viz_output_block *vizblk,
                                   struct sym_table *sym,
                                   char *filename)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->viz_mode != DX_MODE)
  {
    mdlerror(parse_state, "OBJECT_FILE_PREFIXES is only meaningful with the DX VIZ output mode.");
    return 1;
  }

  struct viz_dx_obj *vizp = mdl_get_viz_obj(parse_state, vizblk, sym);
  if (vizp == NULL)
    return 1;
  free(vizp->name);
  vizp->name = filename;
  return 0;
}

static struct viz_child *mdl_get_viz_child(struct viz_output_block *vizblk,
                                           struct object *objp)
{

  struct sym_table *st = retrieve_sym(objp->sym->name, vizblk->viz_children);
  if (st == NULL)
  {
    struct viz_child *vcp = CHECKED_MALLOC_STRUCT(struct viz_child,
                                                  "visualization child object");
    if (vcp == NULL)
      return NULL;

    vcp->obj = objp;
    vcp->viz_state = NULL;
    vcp->next = NULL;
    vcp->parent = NULL;
    vcp->children = NULL;
    if (store_sym(objp->sym->name, VIZ_CHILD, vizblk->viz_children, vcp) == NULL)
    {
      mcell_allocfailed("Failed to store VIZ child object in symbol table.");
      free(vcp);
      return NULL;
    }

    return vcp;
  }
  else
    return (struct viz_child *) st->value;
}

/**************************************************************************
 set_viz_state_value:
    Set the viz_state value for an object and all of its children.

 In: parse_state: parser state
     objp: the object for which to set the state
     viz_state: state to set
 Out: 0 on success, 1 on failure
**************************************************************************/
static int 
set_viz_state_value(struct mdlparse_vars *parse_state,
                    struct viz_output_block *vizblk,
                    struct object *objp,
                    struct viz_child *vcp_parent,
                    int viz_state)
{
  struct viz_child *vcp = mdl_get_viz_child(vizblk, objp);
  if (vcp == NULL)
    return 1;

  if (vcp_parent != NULL && vcp->parent == NULL)
  {
    vcp->next = vcp_parent->children;
    vcp_parent->children = vcp;
    vcp->parent = vcp_parent;
  }

  switch (objp->object_type)
  {
    case META_OBJ:
      for (struct object *child_objp = objp->first_child;
           child_objp != NULL;
           child_objp=child_objp->next)
      {
        if (set_viz_state_value(parse_state, vizblk, child_objp, vcp, viz_state))
          return 1;
      }
      break;

    case BOX_OBJ:
    case POLY_OBJ:
      if (vcp->viz_state==NULL)
      {
        if ((vcp->viz_state = CHECKED_MALLOC_ARRAY(int, objp->n_walls, "viz_state array for geometry"))==NULL)
          return 1;
        for (int i = 0; i < objp->n_walls; i++)
          vcp->viz_state[i] = viz_state;
      }
      else if (viz_state == INCLUDE_OBJ)
      {
        /* Don't override specific with generic */
        for (int i = 0; i < objp->n_walls; i++)
          if (vcp->viz_state[i] == EXCLUDE_OBJ)
            vcp->viz_state[i] = viz_state;
      }
      else
        for (int i = 0; i < objp->n_walls; i++)
          vcp->viz_state[i] = viz_state;
      break;

    case REL_SITE_OBJ:
      /* just do nothing */
      break;

    default:
      mcell_internal_error("Attempt to set viz_state_value for object '%s', which is of invalid type '%d'.",
                           objp->sym->name,
                           objp->object_type);
      break;
  }

  return 0;
}

/**************************************************************************
 mdl_set_object_viz_state:
    Set the viz_state value for an object and all of its children.

 In: parse_state: parser state
     obj_sym: symbol for the object
     viz_state: state to set
 Out: 0 on success, 1 on failure
**************************************************************************/
static int 
mdl_set_object_viz_state(struct mdlparse_vars *parse_state,
                         struct viz_output_block *vizblk,
                         struct sym_table *obj_sym,
                         int viz_state)
{
  struct object *objp = (struct object *) obj_sym->value;

  /* set viz_state value for the object */
  switch (objp->object_type)
  {
    case META_OBJ:
    case BOX_OBJ:
    case POLY_OBJ:
      return set_viz_state_value(parse_state, vizblk, objp, NULL, viz_state);

    case REL_SITE_OBJ:
      mdlerror(parse_state, "Cannot set viz state value of this type of object");
      return 1;

    default:
      mcell_internal_error("Attempt to set viz_state_value for object '%s', which is of invalid type '%d'.",
                           objp->sym->name,
                           objp->object_type);
      break;
  }

  return 0;
}

/**************************************************************************
 mdl_add_viz_object:
    Adds a viz_child for a particular object to the current viz block.

 In: parse_state: parser state
     vizblk: VIZ_OUTPUT block for this frame list
     obj_sym: the symbol for the object to add
     viz_state: the state for this object
 Out: 0 on success, 1 on failure
**************************************************************************/
static int 
mdl_add_viz_object(struct mdlparse_vars *parse_state,
                   struct viz_output_block *vizblk,
                   struct sym_table *obj_sym,
                   int viz_state)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (parse_state->vol->viz_blocks->file_prefix_name == NULL)
  {
    mdlerror(parse_state, "The keyword FILENAME should be specified before any MESHES are included.");
    return 1;
  }

  if (vizblk->viz_mode == DX_MODE)
  {
    struct viz_dx_obj *vizp = mdl_get_viz_obj(parse_state, vizblk, obj_sym);
    if (vizp == NULL)
      return 1;

    if (vizp->name == NULL)
      vizp->name = mdl_strdup(parse_state->vol->viz_blocks->file_prefix_name);
    if (vizp->name == NULL)
    {
      free(vizp->full_name);
      free(vizp);
      return 1;
    }
  }

  if (mdl_set_object_viz_state(parse_state, vizblk, obj_sym, viz_state))
    return 1;

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
int 
mdl_viz_state(struct mdlparse_vars *parse_state,
              int *target,
              double value)
{
  /* Disallow out-of-range values. */
  if (value + 0.001 < (double) EXCLUDE_OBJ  ||
      value - 0.001 > (double) INCLUDE_OBJ)
  {
    mdlerror(parse_state, "Visualization state is out of range "
                   "(must be between -(2^31 - 1) and (2^31 - 2).");
    return 1;
  }

  /* Round to integer. */
  int int_value;
  if (value < 0.0)
    int_value = (int) (value - 0.001);
  else
    int_value = (int) (value + 0.001);

  /* Disallow "special" values. */
  if (int_value == EXCLUDE_OBJ  ||  int_value == INCLUDE_OBJ)
  {
    mdlerror(parse_state, "Visualization states of -(2^31) and (2^31) - 1 "
                   "are reserved for MCell's internal bookkeeping.");
    return 1;
  }

  *target = int_value;
  return 0;
}

/**************************************************************************
  is_instantiated:
     Check if a given object is instantiated.

  In:  parse_state: the parser state
       objp: the object
  Out: 1 if the object is instantiated, 0 otherwise
**************************************************************************/
static int 
is_instantiated(struct mdlparse_vars *parse_state,
                struct object *objp)
{
  while (objp->parent != NULL)
    objp = objp->parent;

  return (objp == parse_state->vol->root_instance) ? 1 : 0;
}

/**************************************************************************
 mdl_set_viz_include_meshes:
    Sets a flag on all of the listed objects, requesting that they be
    visualized.

 In: parse_state: parser state
     vizblk: the viz block to check
     list: the list of symbols
     state: the state to set
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_include_meshes(struct mdlparse_vars *parse_state,
                           struct viz_output_block *vizblk,
                           struct sym_table_list *list,
                           int viz_state)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  struct sym_table_list *stl;
  if (vizblk->viz_mode == DX_MODE  &&  viz_state == INCLUDE_OBJ)
  {
    mdlerror(parse_state, "In DX MODE the state value for the object must be specified.");
    return 1;
  }

  for (stl = list; stl != NULL; stl = stl->next)
  {
    /* Add this object to the visualization. */
    struct object *objp = (struct object *) stl->node->value;

    /* It is an error to try to include an uninstantiated object */
    if (! is_instantiated(parse_state, objp))
    {
      mdlerror_fmt(parse_state,
                   "Cannot produce visualization for the uninstantiated object '%s'",
                   objp->sym->name);
    }

    if (objp->object_type == REL_SITE_OBJ)
      continue;
    if (mdl_add_viz_object(parse_state, vizblk, stl->node, viz_state))
      return 1;
  }
  mem_put_list(parse_state->sym_list_mem, list);
  return 0;
}

/**************************************************************************
 mdl_set_viz_include_mesh_state:
    Sets the viz state of a particular mesh object, indicating whether it
    should be visualized.

 In: parse_state: parser state
     vizblk: the viz block to check
     obj: the list of symbols
     state: the viz state to set
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_include_mesh_state(struct mdlparse_vars *parse_state,
                               struct viz_output_block *vizblk,
                               struct sym_table *obj,
                               int viz_state)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  vizblk->viz_output_flag |= VIZ_SURFACE_STATES;
  return mdl_add_viz_object(parse_state, vizblk, obj, viz_state);
}

/**************************************************************************
 mdl_set_viz_include_all_meshes:
    Sets a flag on a viz block, requesting that all meshes be visualized.

 In: parse_state: parser state
     vizblk: the viz block to check
     viz_state: the desired viz state
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_include_all_meshes(struct mdlparse_vars *parse_state,
                               struct viz_output_block *vizblk,
                               int viz_state)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->viz_mode == DX_MODE)
  {
    mdlerror(parse_state, "In DX MODE, the ALL_MESHES keyword cannot be used.");
    return 1;
  }
  if (viz_state == INCLUDE_OBJ  &&  (vizblk->viz_output_flag & VIZ_ALL_MESHES))
  {
    /* Do nothing - we will not override the old value if we have no specific
     * state value.
     */
  }
  else
    vizblk->default_mesh_state = viz_state;
  vizblk->viz_output_flag |= VIZ_ALL_MESHES;
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
int 
mdl_set_viz_include_molecules(struct mdlparse_vars *parse_state,
                              struct viz_output_block *vizblk,
                              struct sym_table_list *list,
                              int viz_state)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->viz_mode == DX_MODE  &&  viz_state == INCLUDE_OBJ)
  {
    mem_put_list(parse_state->sym_list_mem, list);
    mdlerror(parse_state, "In DX MODE, the state value for the molecule must be specified.");
    return 1;
  }

  /* Mark all specified molecules */
  struct sym_table_list *stl;
  for (stl = list; stl != NULL; stl = stl->next)
  {
    struct species *specp = (struct species *) stl->node->value;
    if (mdl_set_molecule_viz_state(vizblk, specp, viz_state))
      return 1;
  }

  /* free allocated memory  */
  mem_put_list(parse_state->sym_list_mem, list);
  return 0;
}

/**************************************************************************
 mdl_set_viz_include_all_molecules:
    Sets a flag on a viz block, requesting that all species be visualized.

 In: parse_state: parser state
     vizblk: the viz block to check
     viz_state: the desired viz state
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_set_viz_include_all_molecules(struct mdlparse_vars *parse_state,
                                  struct viz_output_block *vizblk,
                                  int viz_state)
{
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  if (vizblk->viz_mode == DX_MODE)
  {
    mdlerror(parse_state, "In DX MODE, the ALL_MOLECULES keyword cannot be used.");
    return 1;
  }
  if (viz_state == INCLUDE_OBJ  &&  (vizblk->viz_output_flag & VIZ_ALL_MOLECULES))
  {
    /* Do nothing - we will not override the old value if we have no specific
     * state value.
     */
  }
  else
    vizblk->default_mol_state = viz_state;
  vizblk->viz_output_flag |= VIZ_ALL_MOLECULES;
  return 0;
}

/**************************************************************************
 mdl_create_viz_frame:
    Create a frame for output in the visualization.

 In: time_type: either OUTPUT_BY_TIME_LIST or OUTPUT_BY_ITERATION_LIST
     type: the type (MESH_GEOMETRY, MOL_POS, etc.)
     iteration_list: list of iterations/times at which to output
 Out: the frame_data_list object, if successful, or NULL if we ran out of memory
**************************************************************************/
static struct frame_data_list *
mdl_create_viz_frame(int time_type,
                     int type,
                     struct num_expr_list *iteration_list)
{

  struct frame_data_list *fdlp;
  fdlp = CHECKED_MALLOC_STRUCT(struct frame_data_list,
                                "VIZ_OUTPUT frame data");
  if (fdlp == NULL)
    return NULL;

  fdlp->list_type = time_type;
  fdlp->type = type;
  fdlp->viz_iteration = -1;
  fdlp->n_viz_iterations = 0;
  fdlp->iteration_list = iteration_list;
  fdlp->curr_viz_iteration = iteration_list;
  return fdlp;
}

/**************************************************************************
 mdl_create_viz_mesh_frames:
    Create one or more mesh frames for output in the visualization.

 In: parse_state: parser state
     time_type: either OUTPUT_BY_TIME_LIST or OUTPUT_BY_ITERATION_LIST
     type: the type (MESH_GEOMETRY, REGION_DATA, etc.)
     viz_mode: visualization mode
     times: list of iterations/times at which to output
 Out: the frame_data_list object, if successful, or NULL if we ran out of memory
**************************************************************************/
static struct frame_data_list *
mdl_create_viz_mesh_frames(struct mdlparse_vars *parse_state,
                           int time_type,
                           int type,
                           int viz_mode,
                           struct num_expr_list_head *times)
{
  struct frame_data_list *frames = NULL;
  struct frame_data_list *new_frame;
  struct num_expr_list *times_sorted;
  if (times->shared)
  {
    times_sorted = mdl_copysort_numeric_list(parse_state, times->value_head);
    if (times_sorted == NULL)
      return NULL;
  }
  else
  {
    mdl_sort_numeric_list(times->value_head);
    times_sorted = times->value_head;
  }

  if ((viz_mode == DREAMM_V3_GROUPED_MODE) || (viz_mode == DREAMM_V3_MODE))
  {
    if ((new_frame = mdl_create_viz_frame(time_type, type, times_sorted)) == NULL)
      return NULL;
    new_frame->next = frames;
    frames = new_frame;
  }
  else if ((viz_mode == DX_MODE) && (type == REG_DATA))
  {
    /* do nothing */
    mdlerror(parse_state, "REGION_DATA cannot be displayed in DX_MODE; please use DREAMM_V3_GROUPED (or DREAMM_V3) mode.");
  }
  else if (viz_mode == DX_MODE)
  {
    if ((type == MESH_GEOMETRY) || (type == ALL_MESH_DATA))
    {
      /* create two frames - SURF_POS and SURF_STATES */
      if ((new_frame = mdl_create_viz_frame(time_type, SURF_POS, times_sorted)) == NULL)
        return NULL;
      new_frame->next = frames;
      frames = new_frame;

      if ((new_frame = mdl_create_viz_frame(time_type, SURF_STATES, times_sorted)) == NULL)
      {
        free(frames);
        return NULL;
      }
      new_frame->next = frames;
      frames = new_frame;
    }
  }
  else if (viz_mode == NO_VIZ_MODE)
  {
    if ((new_frame = mdl_create_viz_frame(time_type, type, times_sorted)) == NULL)
      return NULL;
    new_frame->next = frames;
    frames = new_frame;
  }
  else
    mdlerror(parse_state, "Meshes may not be included in visualization outputs in the selected output mode.");

  return frames;
}

/**************************************************************************
 mdl_new_viz_mesh_frames:
    Adds some new mesh output frames to a list.

 In: parse_state: parser state
     vizblk: the viz block to check
     frames: list to receive frames
     time_type: timing type (OUTPUT_BY_TIME_LIST or ...ITERATION_LIST)
     mesh_item_type: REGION_DATA, etc.
     timelist: list of times in appropriate units (as per time_type)
 Out: 0 on success, 1 on failure
**************************************************************************/
int 
mdl_new_viz_mesh_frames(struct mdlparse_vars *parse_state,
                        struct viz_output_block *vizblk,
                        struct frame_data_list_head *frames,
                        int time_type,
                        int mesh_item_type,
                        struct num_expr_list_head *timelist)
{
  frames->frame_head = frames->frame_tail = NULL;
 // if (vizblk->viz_mode == NO_VIZ_MODE)
 //   return 0;

  struct frame_data_list *fdlp;
  fdlp = mdl_create_viz_mesh_frames(parse_state,
                                    time_type,
                                    mesh_item_type,
                                    vizblk->viz_mode,
                                    timelist);
  if (! fdlp)
    return 1;

  frames->frame_head = fdlp;
  while (fdlp->next != NULL)
    fdlp = fdlp->next;
  frames->frame_tail = fdlp;
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
 N.B. If this function fails in DX_MODE, it may leak memory.
**************************************************************************/
static struct frame_data_list *
mdl_create_viz_mol_frames(struct mdlparse_vars *parse_state,
                          int time_type,
                          int type,
                          int viz_mode,
                          struct num_expr_list_head *times)
{
  struct frame_data_list *frames = NULL;
  struct frame_data_list *new_frame;
  struct num_expr_list *times_sorted;
  if (times->shared)
  {
    times_sorted = mdl_copysort_numeric_list(parse_state, times->value_head);
    if (times_sorted == NULL)
      return NULL;
  }
  else
  {
    mdl_sort_numeric_list(times->value_head);
    times_sorted = times->value_head;
  }

  if ((viz_mode == DREAMM_V3_GROUPED_MODE) || (viz_mode == DREAMM_V3_MODE))
  {
    if ((new_frame = mdl_create_viz_frame(time_type, type, times_sorted)) == NULL)
      return NULL;
    new_frame->next = frames;
    frames = new_frame;
  }
  else if (viz_mode == DX_MODE)
  {
    if ((type == MOL_POS) || (type == ALL_MOL_DATA))
    {
      /* create four frames */
      if ((new_frame = mdl_create_viz_frame(time_type, EFF_POS, times_sorted)) == NULL)
        return NULL;
      new_frame->next = frames;
      frames = new_frame;

      if ((new_frame = mdl_create_viz_frame(time_type, EFF_STATES, times_sorted)) == NULL)
        return NULL;
      new_frame->next = frames;
      frames = new_frame;

      if ((new_frame = mdl_create_viz_frame(time_type, MOL_POS, times_sorted)) == NULL)
        return NULL;
      new_frame->next = frames;
      frames = new_frame;

      if ((new_frame = mdl_create_viz_frame(time_type, MOL_STATES, times_sorted)) == NULL)
        return NULL;
      new_frame->next = frames;
      frames = new_frame;
    }
    else if (type == MOL_ORIENT)
    {
      /* do nothing */
      mdlerror(parse_state, "MOL_ORIENT cannot be displayed in DX_MODE; please use DREAMM_V3_GROUPED (or DREAMM_V3) mode.");
    }
  }
  else if (type == MOL_POS  ||  type == ALL_MOL_DATA  ||  type == MOL_STATES)
  {
    if ((new_frame = mdl_create_viz_frame(time_type, type, times_sorted)) == NULL)
      return NULL;
    new_frame->next = frames;
    frames = new_frame;
  }
  else if (viz_mode == NO_VIZ_MODE)
  {
    /* Create viz frames consistent with other visualization modes */
    if ((new_frame = mdl_create_viz_frame(time_type, type, times_sorted)) == NULL)
      return NULL;
    new_frame->next = frames;
    frames = new_frame;
  }
  else
  {
    mdlerror_fmt(parse_state, "This type of molecule output data (%d) is not valid "
        "for the selected VIZ output mode (%d).", type, viz_mode);
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
int 
mdl_new_viz_mol_frames(struct mdlparse_vars *parse_state,
                       struct viz_output_block *vizblk,
                       struct frame_data_list_head *frames,
                       int time_type,
                       int mol_item_type,
                       struct num_expr_list_head *timelist)
{
  frames->frame_head = frames->frame_tail = NULL;
 // if (vizblk->viz_mode == NO_VIZ_MODE)
 //   return 0;

  struct frame_data_list *fdlp;
  fdlp = mdl_create_viz_mol_frames(parse_state,
                                   time_type,
                                   mol_item_type,
                                   vizblk->viz_mode,
                                   timelist);
  if (! fdlp)
    return 1;

  frames->frame_head = fdlp;
  while (fdlp->next != NULL)
    fdlp = fdlp->next;
  frames->frame_tail = fdlp;
  return 0;
}

/**************************************************************************
 mdl_new_viz_frames:
    Adds some new output frames to a list.

 In: parse_state: parser state
     vizblk: the viz block to check
     frames: list to receive frames
     time_type: timing type (OUTPUT_BY_TIME_LIST or ...ITERATION_LIST)
     type: the type (ALL_FRAME_DATA, etc.)
     times: list of iterations/times at which to output
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_new_viz_frames(struct mdlparse_vars *parse_state,
                       struct viz_output_block *vizblk,
                       struct frame_data_list_head *frames,
                       int time_type,
                       int type,
                       struct num_expr_list_head *times)
{
  frames->frame_head = frames->frame_tail = NULL;
  if (vizblk->viz_mode == NO_VIZ_MODE)
    return 0;

  struct num_expr_list *times_sorted;
  if (times->shared)
  {
    times_sorted = mdl_copysort_numeric_list(parse_state, times->value_head);
    if (times_sorted == NULL)
      return 1;
  }
  else
  {
    mdl_sort_numeric_list(times->value_head);
    times_sorted = times->value_head;
  }

  struct frame_data_list *fdlp;
  fdlp = mdl_create_viz_frame(time_type,
                              type,
                              times_sorted);
  if (! fdlp)
    return 1;
  frames->frame_tail = frames->frame_head = fdlp;
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
                          struct num_expr_list_head *list)
{
  long long step;
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  for (step = 0; step <= parse_state->vol->iterations; step ++)
  {
    struct num_expr_list *nel;
    nel = CHECKED_MALLOC_STRUCT(struct num_expr_list,
                                 "VIZ_OUTPUT time point");
    if (nel == NULL)
      return 1;

    ++ list->value_count;
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
int mdl_new_viz_all_iterations(struct mdlparse_vars *parse_state, struct num_expr_list_head *list)
{
  long long step;
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  for (step = 0; step <= parse_state->vol->iterations; step ++)
  {
    struct num_expr_list *nel;
    nel = CHECKED_MALLOC_STRUCT(struct num_expr_list,
                                 "VIZ_OUTPUT iteration");
    if (nel == NULL)
      return 1;

    ++ list->value_count;
    if (list->value_tail)
      list->value_tail = list->value_tail->next = nel;
    else
      list->value_head = list->value_tail = nel;
    list->value_tail->value = step;
    list->value_tail->next = NULL;
  }
  return 0;
}

/**************************************************************************
 mdl_set_object_viz_state_by_name:
    Set the viz_state value for an object and all of its children.

 In: parse_state: parser state
     vizblk: visualization block currently being built
     obj_sym: symbol for the object
     viz_state: state to set
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_object_viz_state_by_name(struct mdlparse_vars *parse_state,
                                     struct viz_output_block *vizblk,
                                     struct sym_table *obj_symp,
                                     int viz_state)
{
  return mdl_set_object_viz_state(parse_state, vizblk, obj_symp, viz_state);
}

int mdl_set_molecule_viz_state(struct viz_output_block *vizblk,
                               struct species *specp,
                               int viz_state)
{

  /* Make sure not to override a specific state with a generic state. */
  if (viz_state == INCLUDE_OBJ)
  {
    void * const exclude = (void *) (intptr_t) EXCLUDE_OBJ;

    void *oldval = pointer_hash_lookup_ext(&vizblk->parser_species_viz_states,
                                           specp,
                                           specp->hashval,
                                           exclude);
    if (oldval != exclude)
      return 0;
  }
  else
    vizblk->viz_output_flag |= VIZ_MOLECULES_STATES;

  /* Store new value in the hashtable or die trying. */
  void *val = (void *) (intptr_t) viz_state;
  assert(viz_state == (int) (intptr_t) val);
  if (pointer_hash_add(&vizblk->parser_species_viz_states,
                       specp,
                       specp->hashval,
                       val))
  {
    mcell_allocfailed("Failed to store viz state for molecules of species '%s'.",
                      specp->sym->name);
    return 1;
  }
  return 0;
}

/**************************************************************************
 mdl_set_region_viz_state:
    Set the viz_state for a particular region.

 In: parse_state: parser state
     vizblk: the viz block currently being defined
     rp:   region for which to set the vizstate
     viz_state: state to set
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_set_region_viz_state(struct mdlparse_vars *parse_state,
                             struct viz_output_block *vizblk,
                             struct region *rp,
                             int viz_state)
{
  struct object *objp = rp->parent;
  struct polygon_object *pop = (struct polygon_object *) objp->contents;

  /* Only allow referencing instantiated objects. */
  /* XXX: Should we only make this restriction for DREAMM?  (i.e. not for old
   * output modes?) */
  if (! is_instantiated(parse_state, objp))
  {
    mdlerror_fmt(parse_state,
                 "Cannot produce visualization for the uninstantiated region '%s'",
                 rp->sym->name);
  }

  struct viz_child *vcp = mdl_get_viz_child(vizblk, objp);
  if (vcp == NULL)
    return 1;

  const int n_walls = pop->n_walls;
  if (vcp->viz_state==NULL)
  {
    if ((vcp->viz_state = CHECKED_MALLOC_ARRAY(int, n_walls, "viz state array")) == NULL)
      return 1;
    for (int i=0; i<n_walls; i++)
      vcp->viz_state[i] = EXCLUDE_OBJ;
  }

  for (int i=0; i<rp->membership->nbits; ++ i)
  {
    if (get_bit(rp->membership, i))
      vcp->viz_state[i] = viz_state;
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
                                              int *count)
{

  struct species **arr = NULL, **ptr = NULL;
  struct species_list_item *i;

  if (count != NULL)
    *count = lh->species_count;

  if (lh->species_count == 0)
    return NULL;

  arr = CHECKED_MALLOC_ARRAY(struct species *, lh->species_count, "species array");
  if (arr)
  {
    for (i = (struct species_list_item *) lh->species_head, ptr = arr; i != NULL; i = i->next)
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
struct volume_output_item *mdl_new_volume_output_item(struct mdlparse_vars *parse_state,
                                                      char *filename_prefix,
                                                      struct species_list *molecules,
                                                      struct vector3 *location,
                                                      struct vector3 *voxel_size,
                                                      struct vector3 *voxel_count,
                                                      struct output_times *ot)
{
  struct volume_output_item *vo = CHECKED_MALLOC_STRUCT(struct volume_output_item,
                                                        "volume output request");
  if (vo == NULL)
  {
    free(filename_prefix);
    mem_put_list(parse_state->species_list_mem, molecules->species_head);
    free(location);
    free(voxel_size);
    free(voxel_count);
    if (ot->times != NULL) free(ot->times);
    mem_put(parse_state->output_times_mem, ot);
    return NULL;
  }
  memset(vo, 0, sizeof(struct volume_output_item));

  vo->filename_prefix = filename_prefix;
  vo->molecules = species_list_to_array(molecules, &vo->num_molecules);
  if (vo->num_molecules != 0  &&  vo->molecules == NULL)
  {
    free(filename_prefix);
    mem_put_list(parse_state->species_list_mem, molecules->species_head);
    free(location);
    free(voxel_size);
    free(voxel_count);
    if (ot->times != NULL) free(ot->times);
    mem_put(parse_state->output_times_mem, ot);
    free(vo);
    return NULL;
  }

  qsort(vo->molecules, vo->num_molecules, sizeof(void *), & void_ptr_compare);
  mem_put_list(parse_state->species_list_mem, molecules->species_head);

  memcpy(& vo->location, location, sizeof(struct vector3));
  free(location);
  vo->location.x *= parse_state->vol->r_length_unit;
  vo->location.y *= parse_state->vol->r_length_unit;
  vo->location.z *= parse_state->vol->r_length_unit;

  memcpy(& vo->voxel_size, voxel_size, sizeof(struct vector3));
  free(voxel_size);
  vo->voxel_size.x *= parse_state->vol->r_length_unit;
  vo->voxel_size.y *= parse_state->vol->r_length_unit;
  vo->voxel_size.z *= parse_state->vol->r_length_unit;

  vo->nvoxels_x = (int) (voxel_count->x + 0.5);
  vo->nvoxels_y = (int) (voxel_count->y + 0.5);
  vo->nvoxels_z = (int) (voxel_count->z + 0.5);
  free(voxel_count);

  vo->timer_type = ot->timer_type;
  vo->step_time = ot->step_time;
  vo->num_times = ot->num_times;
  vo->times = ot->times;
  mem_put(parse_state->output_times_mem, ot);

  if (vo->timer_type == OUTPUT_BY_ITERATION_LIST  ||
      vo->timer_type == OUTPUT_BY_TIME_LIST)
    vo->next_time = vo->times;
  else
    vo->next_time = NULL;

  switch (parse_state->vol->notify->volume_output_report)
  {
    case NOTIFY_NONE:
      break;

    case NOTIFY_BRIEF:
      mcell_log("Added volume output item '%s', counting in the region [%.15g,%.15g,%.15g]-[%.15g,%.15g,%.15g]",
                vo->filename_prefix,
                vo->location.x * parse_state->vol->length_unit,
                vo->location.y * parse_state->vol->length_unit,
                vo->location.z * parse_state->vol->length_unit,
                (vo->location.x + vo->voxel_size.x * vo->nvoxels_x) * parse_state->vol->length_unit,
                (vo->location.y + vo->voxel_size.x * vo->nvoxels_x) * parse_state->vol->length_unit,
                (vo->location.z + vo->voxel_size.x * vo->nvoxels_x) * parse_state->vol->length_unit);
      break;

    case NOTIFY_FULL:
      mcell_log("Added volume output item '%s', counting the following molecules in the region [%.15g,%.15g,%.15g]-[%.15g,%.15g,%.15g]:\n",
                vo->filename_prefix,
                vo->location.x * parse_state->vol->length_unit,
                vo->location.y * parse_state->vol->length_unit,
                vo->location.z * parse_state->vol->length_unit,
                (vo->location.x + vo->voxel_size.x * vo->nvoxels_x) * parse_state->vol->length_unit,
                (vo->location.y + vo->voxel_size.x * vo->nvoxels_x) * parse_state->vol->length_unit,
                (vo->location.z + vo->voxel_size.x * vo->nvoxels_x) * parse_state->vol->length_unit);
      for (int i=0; i<vo->num_molecules; ++i)
        mcell_log("  %s", vo->molecules[i]->sym->name);
      break;

    default: UNHANDLED_CASE(parse_state->vol->notify->volume_output_report);
  }

  return vo;
}

/**************************************************************************
 mdl_new_output_times_default:
    Create new default output timing for volume output.

 In: parse_state: parser state
 Out: output times structure, or NULL if allocation fails
**************************************************************************/
struct output_times *mdl_new_output_times_default(struct mdlparse_vars *parse_state)
{
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
**************************************************************************/
struct output_times *mdl_new_output_times_step(struct mdlparse_vars *parse_state,
                                               double step)
{
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
  if (output_freq > parse_state->vol->iterations  &&  output_freq > 1)
  {
    output_freq = (parse_state->vol->iterations > 1) ? parse_state->vol->iterations : 1;
    ot->step_time = output_freq * parse_state->vol->time_unit;
    if (parse_state->vol->notify->invalid_output_step_time != WARN_COPE)
      mdl_warning(parse_state, "Output step time too long\n\tSetting output step time to %g microseconds\n", ot->step_time * 1.0e6);
  }
  else if (output_freq < 1)
  {
    ot->step_time = parse_state->vol->time_unit;
    if (parse_state->vol->notify->invalid_output_step_time != WARN_COPE)
      mdl_warning(parse_state, "Output step time too short\n\tSetting output step time to %g microseconds\n", ot->step_time * 1.0e6);
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
struct output_times *mdl_new_output_times_iterations(struct mdlparse_vars *parse_state,
                                                     struct num_expr_list_head *iters)
{
  struct output_times *ot = CHECKED_MEM_GET(parse_state->output_times_mem,
                                             "output times for volume output");
  if (ot == NULL)
  {
    if (! iters->shared)
      mdl_free_numeric_list(iters->value_head);
    return NULL;
  }
  memset(ot, 0, sizeof(struct output_times));

  ot->timer_type = OUTPUT_BY_ITERATION_LIST;
  ot->times = num_expr_list_to_array(iters, &ot->num_times);
  if (ot->times != NULL)
    qsort(ot->times, ot->num_times, sizeof(double), & double_cmp);
  if (! iters->shared)
    mdl_free_numeric_list(iters->value_head);

  if (ot->num_times != 0  &&  ot->times == NULL)
  {
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
struct output_times *mdl_new_output_times_time(struct mdlparse_vars *parse_state,
                                               struct num_expr_list_head *times)
{
  struct output_times *ot = CHECKED_MEM_GET(parse_state->output_times_mem,
                                             "output times for volume output");
  if (ot == NULL)
  {
    if (! times->shared)
      mdl_free_numeric_list(times->value_head);
    return NULL;
  }
  memset(ot, 0, sizeof(struct output_times));

  ot->timer_type = OUTPUT_BY_TIME_LIST;
  ot->times = num_expr_list_to_array(times, &ot->num_times);
  if (ot->times != NULL)
    qsort(ot->times, ot->num_times, sizeof(double), & double_cmp);
  if (! times->shared)
    mdl_free_numeric_list(times->value_head);

  if (ot->num_times != 0  &&  ot->times == NULL)
  {
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
struct sym_table *mdl_new_release_pattern(struct mdlparse_vars *parse_state,
                                          char *name)
{
  struct sym_table *st;
  if (retrieve_sym(name, parse_state->vol->rpat_sym_table) != NULL)
  {
    mdlerror_fmt(parse_state, "Release pattern already defined: %s", name);
    free(name);
    return NULL;
  }
  else if ((st = store_sym(name, RPAT, parse_state->vol->rpat_sym_table, NULL)) == NULL)
  {
    mdlerror_fmt(parse_state, "Out of memory while creating release pattern: %s", name);
    free(name);
    return NULL;
  }

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
                            struct sym_table *rpat_sym,
                            struct release_pattern *rpat_data)
{
  struct release_pattern *rpatp = (struct release_pattern *) rpat_sym->value;
  if (rpat_data->release_interval <= 0)
  {
    mdlerror(parse_state, "Release interval must be set to a positive number.");
    return 1;
  }
  if (rpat_data->train_interval <= 0)
  {
    mdlerror(parse_state, "Train interval must be set to a positive number.");
    return 1;
  }
  if (rpat_data->train_duration > rpat_data->train_interval)
  {
    mdlerror(parse_state, "Train duration must not be longer than the train interval.");
    return 1;
  }
  if (rpat_data->train_duration <= 0)
  {
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
 mdl_valid_complex_name:
    Check that no molecule or named reaction pathway exists with the suplpied
    name.

 In: parse_state: parser state
     name: name for the new species
 Out: 0 if the name would be valid, 1 if not
**************************************************************************/
int mdl_valid_complex_name(struct mdlparse_vars *parse_state, char *name)
{
  if (retrieve_sym(name, parse_state->vol->rxpn_sym_table) != NULL)
  {
    mdlerror_fmt(parse_state, "There is already a named reaction pathway called '%s'.  Please change one of the names.", name);
    free(name);
    return 1;
  }
  else if (retrieve_sym(name, parse_state->vol->mol_sym_table) != NULL)
  {
    mdlerror_fmt(parse_state, "There is already a molecule or complex called '%s'.  Please change one of the names.", name);
    free(name);
    return 1;
  }

  return 0;
}

/**************************************************************************
 mdl_new_mol_species:
    Create a new species. There must not yet be a molecule or named reaction
    pathway with the supplied name.

 In: parse_state: parser state
     name: name for the new species
 Out: symbol for the species, or NULL if an error occurred
**************************************************************************/
struct sym_table *mdl_new_mol_species(struct mdlparse_vars *parse_state, char *name)
{
  struct sym_table *sym = NULL;
  if (retrieve_sym(name, parse_state->vol->mol_sym_table) != NULL)
  {
    mdlerror_fmt(parse_state, "Molecule already defined: %s", name);
    free(name);
    return NULL;
  }
  else if (retrieve_sym(name, parse_state->vol->rxpn_sym_table) != NULL)
  {
    mdlerror_fmt(parse_state, "Molecule already defined as a named reaction pathway: %s", name);
    free(name);
    return NULL;
  }
  else if ((sym = store_sym(name, MOL, parse_state->vol->mol_sym_table, NULL)) == NULL)
  {
    mdlerror_fmt(parse_state, "Out of memory while creating molecule: %s", name);
    free(name);
    return NULL;
  }

  free(name);
  return sym;
}

/**************************************************************************
 mdl_ensure_rdstep_tables_built:
    Build the r_step/d_step tables if they haven't been built yet.

 In: parse_state: parser state
 Out: 0 on success, 1 on failure
**************************************************************************/
static int mdl_ensure_rdstep_tables_built(struct mdlparse_vars *parse_state)
{
  if (parse_state->vol->r_step != NULL  &&
      parse_state->vol->r_step_surface != NULL  &&
      parse_state->vol->d_step != NULL)
    return 0;

  if (parse_state->vol->r_step==NULL)
  {
    if ((parse_state->vol->r_step = init_r_step(parse_state->vol->radial_subdivisions)) == NULL)
    {
      mdlerror(parse_state, "Out of memory while creating r_step data for molecule");
      return 1;
    }
  }

  if (parse_state->vol->r_step_surface==NULL)
  {
    parse_state->vol->r_step_surface = init_r_step_surface(parse_state->vol->radial_subdivisions);
    if (parse_state->vol->r_step_surface==NULL)
    {
      mdlerror(parse_state, "Cannot store r_step_surface data.");
      return 1;
    }
  }

  if (parse_state->vol->d_step==NULL)
  {
    if ((parse_state->vol->d_step = init_d_step(parse_state->vol->radial_directions,&parse_state->vol->num_directions))==NULL)
    {
      mdlerror(parse_state, "Out of memory while creating d_step data for molecule");
      return 1;
    }

    /* Num directions, rounded up to the nearest 2^n - 1 */
    parse_state->vol->directions_mask = parse_state->vol->num_directions;
    parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 1);
    parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 2);
    parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 4);
    parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 8);
    parse_state->vol->directions_mask |= (parse_state->vol->directions_mask >> 16);
    if (parse_state->vol->directions_mask > (1<<18))
    {
      mdlerror(parse_state, "Internal error: bad number of default RADIAL_DIRECTIONS (max 131072).");
      return 1;
    }
  }

  return 0;
}

/**************************************************************************
 mdl_assemble_mol_species:
    Assemble a molecule species from its component pieces.

 In: parse_state: parser state
     sym:  symbol for the species
     D_ref: reference diffusion constant
     D:    diffusion constant
     is_2d: 1 if the species is a 2D molecule, 0 if 3D
     time_step: time_step for the molecule (< 0.0 for a custom space step, >0.0
                for custom timestep, 0.0 for default timestep)
     target_only: 1 if the molecule cannot initiate reactions
 Out: the species, or NULL if an error occurred
 NOTE: This is just a thin wrapper around assemble_mol_species
**************************************************************************/
struct species *mdl_assemble_mol_species(struct mdlparse_vars *parse_state,
                                         struct sym_table *sym,
                                         double D_ref,
                                         double D,
                                         int is_2d,
                                         double time_step,
                                         int target_only,
                                         double max_step_length)
{
  /* Can't define molecule before we have a time step */
  double global_time_unit = parse_state->vol->time_unit;
  if (global_time_unit == 0)
  {
    mdlerror_fmt(parse_state, "TIME_STEP not yet specified.  Cannot define molecule: %s", sym->name);
    return NULL;
  }

  struct species *specp = assemble_mol_species(parse_state->vol,
                                               sym,
                                               D_ref,
                                               D,
                                               is_2d,
                                               time_step,
                                               target_only,
                                               max_step_length);

  if (mdl_ensure_rdstep_tables_built(parse_state))
    return NULL;
  else
    return specp;
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
                   struct reaction_rate *rate)
{
  if (rate->rate_type == RATE_UNSET)
  {
    mdlerror_fmt(parse_state, "File %s, Line %d: Internal error: Rate is not set", __FILE__, __LINE__);
    return 1;
  }

  else if (rate->rate_type == RATE_CONSTANT)
  {
    if (rate->v.rate_constant < 0.0)
    {
      if (parse_state->vol->notify->neg_reaction==WARN_ERROR)
      {
        mdlerror(parse_state, "Error: reaction rates should be zero or positive.");
        return 1;
      }
      else if (parse_state->vol->notify->neg_reaction == WARN_WARN)
      {
        mcell_warn("Warning: negative reaction rate %f; setting to zero and continuing.", rate->v.rate_constant);
        rate->v.rate_constant = 0.0;
      }
    }
  }

  else if (rate->rate_type == RATE_FILE)
  {
    if (rate->v.rate_file == NULL)
    {
      mdlerror_fmt(parse_state, "File %s, Line %d: Internal error: Rate filename is not set", __FILE__, __LINE__);
      return 1;
    }
  }

  else if (rate->rate_type == RATE_COMPLEX)
  {
    if (rate->v.rate_complex == NULL)
    {
      mdlerror_fmt(parse_state, "File %s, Line %d: Internal error: Complex rate structure is not set", __FILE__, __LINE__);
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
static struct species_opt_orient *mdl_new_reaction_player(struct mdlparse_vars *parse_state,
                                                          struct species_opt_orient *spec)
{
  struct species_opt_orient *new_spec;
  if ((new_spec = (struct species_opt_orient *) CHECKED_MEM_GET(parse_state->mol_data_list_mem, "molecule type")) == NULL)
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
                                  struct species_opt_orient_list *list,
                                  struct species_opt_orient *spec)
{
  struct species_opt_orient *player = mdl_new_reaction_player(parse_state, spec);
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
                            struct species_opt_orient_list *list,
                            struct species_opt_orient *spec)
{
  struct species_opt_orient *player = mdl_new_reaction_player(parse_state, spec);
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
                               struct sym_table *symp)
{
  switch (symp->sym_type)
  {
    case DBL:
      rate->rate_type = RATE_CONSTANT;
      rate->v.rate_constant = *(double *) symp->value;
      break;

    case STR:
      rate->rate_type = RATE_FILE;
      if ((rate->v.rate_file = mdl_strdup((char *) symp->value)) == NULL)
        return 1;
      break;

    default:
      mdlerror(parse_state, "Invalid variable used for reaction rate: must be a number or a filename");
      return 1;
  }
  return 0;
}

/**************************************************************************
 mdl_reaction_rate_complex:
    Set a complex reaction rate.

 In: parse_state: parser state
     rate: reaction rate struct
     symp: symbol for complex
     tbl:  rate rule table name
 Out: 0 on success, 1 on failure
**************************************************************************/
int mdl_reaction_rate_complex(struct mdlparse_vars *parse_state,
                              struct reaction_rate *rate,
                              struct sym_table *symp,
                              char *tbl)
{
  struct species *complex_species = (struct species *) symp->value;
  if (! (complex_species->flags & IS_COMPLEX))
  {
    mdlerror_fmt(parse_state,
                 "The molecule '%s' specified in the complex reaction rate is not a complex.",
                 symp->name);
    free(tbl);
    return 1;
  }

  rate->rate_type = RATE_COMPLEX;
  if ((rate->v.rate_complex = macro_lookup_ruleset((struct complex_species *) complex_species, tbl)) == NULL)
  {
    mdlerror_fmt(parse_state,
                 "The complex '%s' has no rate rule named '%s'.",
                 symp->name, tbl);
    free(tbl);
    return 1;
  }
  free(tbl);
  return 0;
}

/*************************************************************************
 create_rx_name:
    Assemble reactants alphabetically into a reaction name string.

 In:  p: reaction pathway whose reaction name we are to create
 Out: a string to be used as a symbol name for the reaction
*************************************************************************/
static char *create_rx_name( struct pathway *p)
{

  struct species *reagents[3];
  int n_reagents = 0;
  int is_complex = 0;

  /* Store reagents in an array. */
  reagents[0] = p->reactant1;
  reagents[1] = p->reactant2;
  reagents[2] = p->reactant3;

  /* Count non-null reagents. */
  for (n_reagents = 0; n_reagents < 3; ++ n_reagents)
    if (reagents[n_reagents] == NULL)
      break;
    else if (p->is_complex[n_reagents])
      is_complex = 1;

  /* Sort reagents. */
  for (int i = 0; i<n_reagents; ++i)
  {
    for (int j = i+1; j<n_reagents; ++ j)
    {
      /* If 'i' is a subunit, 'i' wins. */
      if (p->is_complex[i])
        break;

      /* If 'j' is a subunit, 'j' wins. */
      else if (p->is_complex[j])
      {
        struct species *tmp = reagents[j];
        reagents[j] = reagents[i];
        reagents[i] = tmp;
      }

      /* If 'j' precedes 'i', 'j' wins. */
      else if (strcmp(reagents[j]->sym->name, reagents[i]->sym->name) < 0)
      {
        struct species *tmp = reagents[j];
        reagents[j] = reagents[i];
        reagents[i] = tmp;
      }
    }
  }

  /* Now, produce a name! */
  if (is_complex)
  {
    switch (n_reagents)
    {
      case 1: return alloc_sprintf("(%s)", reagents[0]->sym->name);
      case 2: return alloc_sprintf("(%s)+%s", reagents[0]->sym->name, reagents[1]->sym->name);
      case 3: return alloc_sprintf("(%s)+%s+%s", reagents[0]->sym->name, reagents[1]->sym->name, reagents[2]->sym->name);
      default:
        mcell_internal_error("Invalid number of reagents in reaction pathway (%d).", n_reagents);
        return NULL;
    }
  }
  else
  {
    switch (n_reagents)
    {
      case 1: return alloc_sprintf("%s", reagents[0]->sym->name);
      case 2: return alloc_sprintf("%s+%s", reagents[0]->sym->name, reagents[1]->sym->name);
      case 3: return alloc_sprintf("%s+%s+%s", reagents[0]->sym->name, reagents[1]->sym->name, reagents[2]->sym->name);
      default:
        mcell_internal_error("Invalid number of reagents in reaction pathway (%d).", n_reagents);
        return NULL;
    }
  }
}

/*************************************************************************
 sort_product_list_compare:
    Comparison function for products to be sorted when generating the product
    signature.

 In:  list_item: first item to compare
      new_item:  second item to compare
 Out: -1 if list_item < new_item, 1 if list_item > new_item, 0 if they are
      equal
*************************************************************************/
static int sort_product_list_compare(struct product *list_item,
                                     struct product *new_item)
{
  int cmp = list_item->is_complex - new_item->is_complex;
  if (cmp != 0)
    return cmp;

  cmp = strcmp(list_item->prod->sym->name, new_item->prod->sym->name);
  if (cmp == 0)
  {
    if (list_item->orientation > new_item->orientation)
      cmp = -1;
    else if (list_item->orientation < new_item->orientation)
      cmp = 1;
    else
      cmp = 0;
  }
  return cmp;
}

/*************************************************************************
 sort_product_list:
    Sorts product_head in alphabetical order, and descending orientation order.
    Current algorithm uses insertion sort.

 In:  product_head: list to sort
 Out: the new list head
*************************************************************************/
static struct product *sort_product_list(struct product *product_head)
{
  struct product *next;             /* Saved next item (next field in product is overwritten) */
  struct product *iter;             /* List iterator */
  struct product *result = NULL;    /* Sorted list */
  int cmp;

  /* Use insertion sort to sort the list of products */
  for (struct product *current = product_head;
       current != NULL;
       current = next)
  {
    next = current->next;

    /* First item added always goes at the head */
    if (result == NULL)
    {
      current->next = result;
      result = current;
      continue;
    }

    /* Check if the item belongs at the head */
    cmp = sort_product_list_compare(result, current);
    if (cmp >= 0)
    {
      current->next = result;
      result = current;
      continue;
    }

    /* Otherwise, if it goes after the current entry, scan forward to find the insert point */
    else
    {
      /* locate the node before the point of insertion */
      iter = result;
      while (iter->next != NULL  &&  sort_product_list_compare(iter, current) < 0)
        iter = iter->next;

      current->next = iter->next;
      iter->next = current;
    }
  }

  return result;
}

/*************************************************************************
 create_prod_signature:
    Returns a string containing all products in the product_head list,
    separated by '+', and sorted in alphabetical order by name and descending
    orientation order.

 In:  product_head: list of products
 Out: product signature as a string.  *product_head list is sorted in
      alphabetical order by name, and descending order by orientation.  Returns
      NULL on failure.
*************************************************************************/
static char *create_prod_signature( struct product **product_head)
{

  /* points to the head of the sorted alphabetically list of products */
  char *prod_signature = NULL;

  *product_head = sort_product_list(*product_head);

  /* create prod_signature string */
  struct product *current = *product_head;
  prod_signature = CHECKED_STRDUP(current->prod->sym->name, "product name");

  /* Concatenate to create product signature */
  char *temp_str = NULL;
  while (current->next != NULL)
  {
    temp_str = prod_signature;
    prod_signature = CHECKED_SPRINTF("%s+%s",
                                     prod_signature,
                                     current->next->prod->sym->name);

    if (prod_signature == NULL)
    {
      if (temp_str != NULL) free(temp_str);
      return NULL;
    }
    if (temp_str != NULL) free(temp_str);

    current = current->next;
  }

  return prod_signature;
}

/*************************************************************************
 concat_rx_name:
    Concatenates reactants onto a reaction name.  Reactants which are subunits
    in macromolecular complexes will have their names parenthesized.

 In:  parse_state: parser state
      name1: name of first reactant (or first part of reaction name)
      is_complex1: 0 unless the first reactant is a subunit in a complex
      name2: name of second reactant (or second part of reaction name)
      is_complex2: 0 unless the second reactant is a subunit in a complex
 Out: reaction name as a string, or NULL if an error occurred
*************************************************************************/
static char *concat_rx_name(struct mdlparse_vars *parse_state,
                            char *name1,
                            int is_complex1,
                            char *name2,
                            int is_complex2)
{
  char *rx_name;

  /* Make sure they aren't both subunits  */
  if (is_complex1  &&  is_complex2)
  {
    mdlerror_fmt(parse_state, "File '%s', Line %ld: Internal error -- a reaction cannot have two reactants which are subunits of a macromolecule.", __FILE__, (long)__LINE__);
    return NULL;
  }

  /* Sort them */
  if (is_complex2  ||  strcmp(name2, name1) <= 0)
  {
    char *nametmp = name1;
    int is_complextmp = is_complex1;
    name1 = name2;
    is_complex1 = is_complex2;
    name2 = nametmp;
    is_complex2 = is_complextmp;
    assert(is_complex2 == 0);
  }

  /* Build the name */
  if (is_complex1)
    rx_name = CHECKED_SPRINTF("(%s)+%s", name1, name2);
  else
    rx_name = CHECKED_SPRINTF("%s+%s", name1, name2);

  /* Die if we failed to allocate memory */
  if (rx_name == NULL)
    return NULL;

  return rx_name;
}

/***********************************************************************
 invert_current_reaction_pathway:
    Creates a new reversed pathway, where the reactants of new pathway are the
    products of the current pathway, and the products of new pathway are the
    reactants of the current pathway.

 In:  parse_state: parser state
      pathp: pathway to invert
      reverse_rate: the reverse reaction rate
 Out: Returns 1 on error and 0 - on success.  The new pathway is added to the
      linked list of the pathways for the current reaction.
***********************************************************************/
static int invert_current_reaction_pathway(struct mdlparse_vars *parse_state,
                                           struct pathway *pathp,
                                           struct reaction_rate *reverse_rate)
{
  struct rxn *rx;
  struct pathway *path;
  struct product *prodp;
  struct sym_table *sym;
  char *inverse_name;
  int nprods;  /* number of products */
  int all_3d;  /* flag that tells whether all products are volume_molecules */
  int num_surf_products = 0;
  int num_grid_mols = 0;
  int num_vol_mols = 0;


  /* flag that tells whether there is a surface_class
     among products in the direct reaction */
  int is_surf_class = 0;

  all_3d=1;
  for (nprods=0,prodp=pathp->product_head ; prodp!=NULL ; prodp=prodp->next)
  {
    nprods++;
    if ((prodp->prod->flags&NOT_FREE)!=0) all_3d=0;
    if ((prodp->prod->flags & IS_SURFACE) != 0)
    {
           is_surf_class = 1;
    }
  }

  if (nprods==0)
  {
    mdlerror(parse_state, "Can't create a reverse reaction with no products");
    return 1;
  }
  if (nprods==1 && (pathp->product_head->prod->flags&IS_SURFACE))
  {
    mdlerror(parse_state, "Can't create a reverse reaction starting from only a surface");
    return 1;
  }
  if (nprods>3)
  {
    mdlerror(parse_state, "Can't create a reverse reaction involving more than three products. Please note that surface_class from the reaction reactant side also counts as a product.");
    return 1;
  }

  if (pathp->pathname != NULL)
  {
    mdlerror(parse_state, "Can't name bidirectional reactions--write each reaction and name them separately");
    return 1;
  }
  if (all_3d)
  {
    if ((pathp->reactant1->flags&NOT_FREE)!=0) all_3d = 0;
    if (pathp->reactant2!=NULL && (pathp->reactant2->flags&NOT_FREE)!=0) all_3d = 0;
    if (pathp->reactant3!=NULL && (pathp->reactant3->flags&NOT_FREE)!=0) all_3d = 0;

    if (!all_3d)
    {
      mdlerror(parse_state, "Cannot reverse orientable reaction with only volume products");
      return 1;
    }
  }

  prodp = pathp->product_head;
  if (nprods==1)
  {
    if (prodp->is_complex)
    {
      inverse_name = CHECKED_SPRINTF("(%s)",
                                     prodp->prod->sym->name);
    }
    else
      inverse_name = mdl_strdup(prodp->prod->sym->name);
    if (inverse_name == NULL)
      return 1;
  }
  else if (nprods == 2)
  {
    inverse_name = concat_rx_name(parse_state, prodp->prod->sym->name, prodp->is_complex, prodp->next->prod->sym->name, prodp->next->is_complex);
  }
  else
  {
    if (prodp->is_complex || prodp->next->is_complex || prodp->next->next->is_complex)
    {
      mdlerror(parse_state, "MCell does not currently support trimolecular reactions for macromolecules");
      return 1;
    }
    inverse_name = concat_rx_name(parse_state, prodp->prod->sym->name, 0, prodp->next->prod->sym->name, 0);
    inverse_name = concat_rx_name(parse_state, inverse_name, 0, prodp->next->next->prod->sym->name, 0);
  }
  if (inverse_name==NULL)
  {
    mdlerror(parse_state, "Out of memory forming reaction name");
    return 1;
  }

  sym = retrieve_sym(inverse_name, parse_state->vol->rxn_sym_table);
  if (sym==NULL)
  {
    sym = store_sym(inverse_name,RX,parse_state->vol->rxn_sym_table, NULL);
    if (sym==NULL)
    {
      mdlerror_fmt(parse_state, "File '%s', Line %ld: Out of memory while storing reaction pathway.", __FILE__, (long)__LINE__);
      return 1;
    }
  }
  free(inverse_name);
  rx = (struct rxn*)sym->value;
  rx->n_reactants = nprods;
  rx->n_pathways++;

  if ((path = (struct pathway*) CHECKED_MEM_GET(parse_state->path_mem, "reaction pathway")) == NULL)
      return 1;
  path->pathname=NULL;
  path->flags = 0;
  path->reactant1=prodp->prod;
  if ((path->reactant1->flags & NOT_FREE) == 0)
  {
         ++ num_vol_mols;
  }
  else 
  {
     if (path->reactant1->flags & ON_GRID)
     {
         ++ num_grid_mols;
     }
  }
  path->is_complex[0] = prodp->is_complex;
  path->is_complex[1] = 0;
  path->is_complex[2] = 0;
  path->orientation1 = prodp->orientation;
  path->reactant2=NULL;
  path->reactant3=NULL;
  path->prod_signature = NULL;
  if (nprods > 1)
  {
      path->reactant2 = prodp->next->prod;
      if ((path->reactant2->flags & NOT_FREE) == 0)
      {
         ++ num_vol_mols;
      }
      else 
      {
         if (path->reactant2->flags & ON_GRID)
         {
           ++ num_grid_mols;
         }
      }
      path->orientation2 = prodp->next->orientation;
      path->is_complex[1] = prodp->next->is_complex;
  }
  if (nprods > 2)
  {
      path->reactant3 = prodp->next->next->prod;
      if ((path->reactant3->flags & NOT_FREE) == 0)
      {
         ++ num_vol_mols;
      }
      else 
      {
         if (path->reactant3->flags & ON_GRID)
         {
           ++ num_grid_mols;
         }
      }
      path->orientation3 = prodp->next->next->orientation;
  }

  switch (reverse_rate->rate_type)
  {
    case RATE_UNSET:
      mdlerror_fmt(parse_state, "File %s, Line %d: Internal error: Reverse rate is not set", __FILE__, __LINE__);
      return 1;

    case RATE_CONSTANT:
      path->km = reverse_rate->v.rate_constant;
      path->km_filename = NULL;
      path->km_complex = NULL;
      break;

    case RATE_FILE:
      path->km = 0.0;
      path->km_filename = mdl_find_include_file(parse_state, reverse_rate->v.rate_file, parse_state->vol->curr_file);
      free(reverse_rate->v.rate_file);
      path->km_complex = NULL;
      break;

    case RATE_COMPLEX:
      path->km = 0.0;
      path->km_filename = NULL;
      path->km_complex = reverse_rate->v.rate_complex;
      break;

    default: UNHANDLED_CASE(reverse_rate->rate_type);
  }

  path->product_head = (struct product*) CHECKED_MEM_GET(parse_state->prod_mem, "reaction product");
  if (path->product_head == NULL)
    return 1;

  path->product_head->orientation = pathp->orientation1;
  path->product_head->prod = pathp->reactant1;
  path->product_head->is_complex = pathp->is_complex[0];
  path->product_head->next = NULL;
  if (path->product_head->prod->flags & ON_GRID)
    ++ num_surf_products;

  if ((pathp->reactant2!=NULL) && ((pathp->reactant2->flags & IS_SURFACE) == 0))
  {
    path->product_head->next = (struct product*) CHECKED_MEM_GET(parse_state->prod_mem, "reaction product");
    if (path->product_head->next == NULL)
      return 1;
    path->product_head->next->orientation = pathp->orientation2;
    path->product_head->next->prod = pathp->reactant2;
    path->product_head->next->is_complex = pathp->is_complex[1];
    path->product_head->next->next = NULL;
    if (path->product_head->next->prod->flags & ON_GRID)
      ++ num_surf_products;
  }

  if ((pathp->reactant3!=NULL) && ((pathp->reactant3->flags & IS_SURFACE) == 0))
  {
    path->product_head->next->next = (struct product*) CHECKED_MEM_GET(parse_state->prod_mem, "reaction product");
    if (path->product_head->next->next == NULL)
      return 1;
    path->product_head->next->next->orientation = pathp->orientation3;
    path->product_head->next->next->prod = pathp->reactant3;
    path->product_head->next->next->is_complex = pathp->is_complex[2];
    path->product_head->next->next->next = NULL;
    if (path->product_head->next->next->prod->flags & ON_GRID)
      ++ num_surf_products;
  }

  path->prod_signature = create_prod_signature(&path->product_head);
  if (path->prod_signature == NULL)
  {
      mdlerror(parse_state, "Error creating 'prod_signature' field for reaction pathway.");
      return 1;
  }


  if ((parse_state->vol->vacancy_search_dist2 == 0) &&
     (num_surf_products > num_grid_mols))
  {
      /* the case with one volume molecule reacting with the surface
         and producing one grid molecule is excluded */
      if (!((num_grid_mols == 0) && (num_vol_mols == 1)))
      {
          mdlerror(parse_state, "Error: number of surface products exceeds number of surface reactants, but VACANCY_SEARCH_DISTANCE is not specified or set to zero.");
           return 1;
      }
  }

  /* Now go back to the original reaction and if there is a "surface_class"
     among products - remove it.  We do not need it now on the product side
     of the reaction */
   if (is_surf_class)
   {
     prodp = pathp->product_head;
     if (prodp->prod->flags & IS_SURFACE)
     {
       pathp->product_head = prodp->next;
       prodp->next = NULL;
       mem_put(parse_state->prod_mem, (void *)prodp);
     }
     else if (prodp->next->prod->flags & IS_SURFACE)
     {
       struct product *temp = prodp->next;
       prodp->next = prodp->next->next;
       mem_put(parse_state->prod_mem, temp);
     }
     else
     {
       struct product *temp = prodp->next->next;
       prodp->next->next = prodp->next->next->next;
       mem_put(parse_state->prod_mem, temp);
     }
   }

  path->next = rx->pathway_head;
  rx->pathway_head = path;
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
struct rxn *mdl_assemble_reaction(struct mdlparse_vars *parse_state,
                                  struct species_opt_orient *reactants,
                                  struct species_opt_orient *surface_class,
                                  struct reaction_arrow *react_arrow,
                                  struct species_opt_orient *products,
                                  struct reaction_rates *rate,
                                  struct sym_table *pathname)
{
  char *rx_name;
  struct sym_table *symp;
  int bidirectional = 0;
  int catalytic = -1;
  int surface = -1;
  unsigned int oriented_count = 0;
  unsigned int num_surfaces = 0;
  int num_surf_products = 0;
  int num_grid_mols = 0;
  int num_vol_mols = 0;
  int num_complex_reactants = 0;
  int all_3d = 1;
  struct rxn *rxnp;
  struct pathway *pathp;

  /* Create pathway */
  pathp = (struct pathway*) CHECKED_MEM_GET(parse_state->path_mem, "reaction pathway");
  if (pathp == NULL)
    return NULL;
  memset(pathp, 0, sizeof(struct pathway));

  // test reaction addition
  mcell_add_reaction(parse_state->vol, reactants, react_arrow, surface_class);


  /* Scan reactants, copying into the new pathway */
  struct species_opt_orient *current_reactant;
  int reactant_idx = 0;
  parse_state->complex_type = 0;
  for (current_reactant = reactants;
       reactant_idx < 3 && current_reactant != NULL;
       ++ reactant_idx, current_reactant = current_reactant->next)
  {
    /* Extract orientation and species */
    short orient = current_reactant->orient_set ? current_reactant->orient : 0;
    struct species *reactant_species = (struct species *) current_reactant->mol_type->value;

    /* Count the type of this reactant */
    if (current_reactant->orient_set)
      ++ oriented_count;
    if (reactant_species->flags & NOT_FREE)
    {
      all_3d = 0;
      if (reactant_species->flags & ON_GRID)
        ++ num_grid_mols;
    }
    else
      ++ num_vol_mols;

    /* Sanity check this reactant */
    if (current_reactant->is_subunit)
    {
      if ((reactant_species->flags & NOT_FREE) == 0)
        parse_state->complex_type = TYPE_3D;
      else if (reactant_species->flags & ON_GRID)
        parse_state->complex_type = TYPE_GRID;
      else
      {
        mdlerror(parse_state, "Only a molecule may be used as a macromolecule subunit in a reaction.");
        return NULL;
      }
      ++ num_complex_reactants;
    }

    else if (reactant_species->flags & IS_SURFACE)
    {
      mdlerror(parse_state, "Surface class can be listed only as the last reactant on the left-hand side of the reaction with the preceding '@' sign.");
      return NULL;
    }

    /* Copy in reactant info */
    pathp->is_complex[reactant_idx] = current_reactant->is_subunit;
    switch (reactant_idx)
    {
      case 0:
        pathp->reactant1 = reactant_species;
        pathp->orientation1 = orient;
        break;

      case 1:
        pathp->reactant2 = reactant_species;
        pathp->orientation2 = orient;
        break;

      case 2:
        pathp->reactant3 = reactant_species;
        pathp->orientation3 = orient;
        break;

      default: UNHANDLED_CASE(reactant_idx);
    }
  }
  mem_put_list(parse_state->mol_data_list_mem, reactants);

  /* Only one complex reactant allowed */
  if (num_complex_reactants > 1)
  {
    mdlerror(parse_state, "Reaction may not include more than one reactant which is a subunit in a complex.");
    return NULL;
  }

  /* If we have reactants left over here, we had more than 3 */
  if (current_reactant != NULL)
  {
    mdlerror(parse_state, "Too many reactants--maximum number is three.");
    return NULL;
  }

  /* Grab info from the arrow */
  if (react_arrow->flags & ARROW_BIDIRECTIONAL)
    bidirectional = 1;
  if (react_arrow->flags & ARROW_CATALYTIC)
  {
    struct species *catalyst_species = (struct species *) react_arrow->catalyst.mol_type->value;
    short orient = react_arrow->catalyst.orient_set ? react_arrow->catalyst.orient : 0;

    /* XXX: Should surface class be allowed inside a catalytic arrow? */
    if (catalyst_species->flags & IS_SURFACE)
    {
      mdlerror(parse_state, "A surface classes may not appear inside a catalytic arrow");
      return NULL;
    }

    /* Count the type of this reactant */
    if (react_arrow->catalyst.orient_set)
      ++ oriented_count;
    if (catalyst_species->flags & NOT_FREE)
    {
      all_3d = 0;
      if (catalyst_species->flags & ON_GRID)
        ++ num_grid_mols;
    }
    else
      ++ num_vol_mols;

    if (reactant_idx >= 3)
    {
      mdlerror(parse_state, "Too many reactants--maximum number is three.");
      return NULL;
    }

    /* Copy in catalytic reactant */
    switch (reactant_idx)
    {
      case 1:
        pathp->reactant2 = (struct species*) react_arrow->catalyst.mol_type->value;
        pathp->orientation2 = orient;
        break;

      case 2:
        pathp->reactant3 = (struct species*) react_arrow->catalyst.mol_type->value;
        pathp->orientation3 = orient;
        break;

      case 0:
      default:
        mcell_internal_error("Catalytic reagent ended up in an invalid slot (%d).", reactant_idx);
        return NULL;
    }
    catalytic = reactant_idx;
    ++ reactant_idx;
  }

  /* If a surface was specified, include it */
  if (surface_class->mol_type != NULL)
  {
    short orient = surface_class->orient_set ? surface_class->orient : 0;
    if (surface_class->orient_set)
      ++ oriented_count;

    /* Copy reactant into next available slot */
    switch (reactant_idx)
    {
      case 0:
        mdlerror(parse_state, "Before defining reaction surface class at least one reactant should be defined.");
        return NULL;

      case 1:
        pathp->reactant2 = (struct species*) surface_class->mol_type->value;
        pathp->orientation2 = orient;
        break;

      case 2:
        pathp->reactant3 = (struct species*) surface_class->mol_type->value;
        pathp->orientation3 = orient;
        break;

      default:
        mdlerror(parse_state, "Too many reactants--maximum number is two plus reaction surface class.");
        return NULL;
    }

    surface = reactant_idx;
    ++ reactant_idx;
    ++ num_surfaces;
    all_3d = 0;
  }

  /* Create a reaction name for the pathway we're creating */
  rx_name = create_rx_name(pathp);
  if (rx_name==NULL)
  {
    mdlerror(parse_state, "Out of memory while creating reaction.");
    return NULL;
  }

  /* If this reaction doesn't exist, create it */
  if ((symp = retrieve_sym(rx_name, parse_state->vol->rxn_sym_table)) != NULL)
  {
    /* do nothing */
  }
  else if ((symp = store_sym(rx_name,RX,parse_state->vol->rxn_sym_table, NULL)) == NULL)
  {
    mdlerror(parse_state, "Out of memory while creating reaction.");
    free(rx_name);
    return NULL;
  }
  free(rx_name);

  rxnp = (struct rxn*) symp->value;
  rxnp->n_reactants = reactant_idx;
  ++ rxnp->n_pathways;

  /* Check for invalid reaction specifications */
  if (num_surfaces > 1)
  {
    /* Shouldn't happen */
    mcell_internal_error("Too many surfaces--reactions can take place on at most one surface.");
    return NULL;
  }
  if (num_surfaces == rxnp->n_reactants)
  {
    mdlerror(parse_state, "Reactants cannot consist entirely of surfaces.  Use a surface release site instead!");
    return NULL;
  }
  if ((num_vol_mols == 2) && (num_surfaces == 1))
  {
    mdlerror(parse_state, "Reaction between two volume molecules and a surface is not defined.");
    return NULL;
  }
  if (all_3d)
  {
    if (oriented_count != 0)
    {
      if (parse_state->vol->notify->useless_vol_orient==WARN_ERROR)
      {
        mdlerror(parse_state, "Error: orientation specified for molecule in reaction in volume");
        return NULL;
      }
      else if (parse_state->vol->notify->useless_vol_orient==WARN_WARN)
      {
        mdlerror(parse_state, "Warning: orientation specified for molecule in reaction in volume");
      }
    }
  }
  else
  {
    if (rxnp->n_reactants != oriented_count)
    {
      if (parse_state->vol->notify->missed_surf_orient==WARN_ERROR)
      {
        mdlerror_fmt(parse_state, "Error: orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
        return NULL;
      }
      else if (parse_state->vol->notify->missed_surf_orient==WARN_WARN)
      {
        mdlerror_fmt(parse_state, "Warning: orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
      }
    }
  }

  /* Add catalytic reagents to the product list.
   *    - For unidirectional catalytic reactions - copy catalyst to products
   *      only if catalyst is not a surface_clas.
   *    - For bidirectional catalytic reactions always copy catalyst to
   *      products and take care that surface_class will not appear in the
   *      products later after inverting the reaction
   */
  if (catalytic >= 0)
  {
    struct species *catalyst;
    short catalyst_orient;
    switch (catalytic)
    {
      case 0: catalyst = pathp->reactant1; catalyst_orient = pathp->orientation1; break;
      case 1: catalyst = pathp->reactant2; catalyst_orient = pathp->orientation2; break;
      case 2: catalyst = pathp->reactant3; catalyst_orient = pathp->orientation3; break;
      default:
        mcell_internal_error("Catalytic reagent index is invalid.");
        return NULL;
    }

    if (bidirectional || (! (catalyst->flags & IS_SURFACE)))
    {
      struct product *prodp;
      prodp = (struct product*) CHECKED_MEM_GET(parse_state->prod_mem,
                                                 "reaction product");
      if (prodp == NULL)
        return NULL;

      prodp->is_complex = 0;
      prodp->prod = catalyst;
      if (all_3d) prodp->orientation = 0;
      else prodp->orientation = catalyst_orient;
      prodp->next = pathp->product_head;
      pathp->product_head = prodp;
    }
  }

  /* Add in all products */
  struct species_opt_orient *current_product;
  int num_complex_products = 0;
  for (current_product = products;
       current_product != NULL;
       current_product = current_product->next)
  {
    /* Nothing to do for NO_SPECIES */
    if (current_product->mol_type == NULL)
      continue;

    /* Create new product */
    struct product *prodp = (struct product*) mem_get(parse_state->prod_mem);
    if (prodp == NULL)
    {
      mdlerror_fmt(parse_state,
                   "Out of memory while creating reaction: %s -> ... ",
                   rxnp->sym->name);
      return NULL;
    }

    /* Set product species and orientation */
    prodp->prod = (struct species *) current_product->mol_type->value;
    if (all_3d) prodp->orientation = 0;
    else prodp->orientation = current_product->orient;

    /* Disallow surface as product unless reaction is bidirectional */
    if (! bidirectional)
    {
      if (prodp->prod->flags & IS_SURFACE)
      {
        mdlerror_fmt(parse_state,
                     "Surface_class '%s' is not allowed to be on the product side of the reaction.",
                     prodp->prod->sym->name);
        return NULL;
      }
    }

    /* Copy over complex-related state for product */
    prodp->is_complex = current_product->is_subunit;
    if (current_product->is_subunit)
    {
      ++ num_complex_products;
      if ((prodp->prod->flags & NOT_FREE) != 0)
      {
        if (parse_state->complex_type == TYPE_3D)
        {
          mdlerror_fmt(parse_state,
                       "Volume subunit cannot become a surface subunit '%s' in a macromolecular reaction.",
                       prodp->prod->sym->name);
          return NULL;
        }
      }
      else if ((prodp->prod->flags & ON_GRID) == 0)
      {
        if (parse_state->complex_type == TYPE_GRID)
        {
          mdlerror_fmt(parse_state,
                       "Surface subunit cannot become a volume subunit '%s' in a macromolecular reaction.",
                       prodp->prod->sym->name);
          return NULL;
        }
      }
      else
      {
        mdlerror(parse_state, "Only a molecule may be used as a macromolecule subunit in a reaction.");
        return NULL;
      }
    }

    /* Append product to list */
    prodp->next = pathp->product_head;
    pathp->product_head = prodp;

    if (prodp->prod->flags & ON_GRID)
      ++ num_surf_products;

    /* Add product if it isn't a surface */
    if (! (prodp->prod->flags&IS_SURFACE))
    {
      if (all_3d == 0)
      {
        if (! current_product->orient_set)
        {
          if (parse_state->vol->notify->missed_surf_orient==WARN_ERROR)
          {
            mdlerror_fmt(parse_state, "Error: product orientation not specified in reaction with orientation\n  (use ; or ', or ,' for random orientation)");
            return NULL;
          }
          else if (parse_state->vol->notify->missed_surf_orient==WARN_WARN)
          {
            mdlerror_fmt(parse_state, "Warning: product orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
          }
        }
      }
      else
      {
        if ((prodp->prod->flags&NOT_FREE)!=0)
        {
          mdlerror(parse_state, "Reaction has only volume reactants but is trying to create a surface product");
          return NULL;
        }
        if (current_product->orient_set)
        {
          if (parse_state->vol->notify->useless_vol_orient==WARN_ERROR)
          {
            mdlerror(parse_state, "Error: orientation specified for molecule in reaction in volume");
            return NULL;
          }
          else if (parse_state->vol->notify->useless_vol_orient==WARN_WARN)
          {
            mdlerror(parse_state, "Warning: orientation specified for molecule in reaction at surface");
          }
        }
      }
    }
  }
  mem_put_list(parse_state->mol_data_list_mem, products);

  /* Subunits can neither be created nor destroyed */
  if (num_complex_reactants != num_complex_products)
  {
    mdlerror_fmt(parse_state, "Reaction must include the same number of complex-subunits on each side of the reaction (have %d reactants vs. %d products)", num_complex_reactants, num_complex_products);
    return NULL;
  }

  /* Attach reaction pathway name, if we have one */
  if (pathname != NULL)
  {
    struct rxn_pathname *rxpnp = (struct rxn_pathname *) pathname->value;
    rxpnp->rx = rxnp;
    pathp->pathname = rxpnp;
  }

  if (pathp->product_head != NULL)
  {
    pathp->prod_signature = create_prod_signature(&pathp->product_head);
    if (pathp->prod_signature == NULL)
    {
      mdlerror(parse_state, "Error creating 'prod_signature' field for the reaction pathway.");
      return NULL;
    }
  }
  else
    pathp->prod_signature = NULL;

  /* Copy in forward rate */
  switch (rate->forward_rate.rate_type)
  {
    case RATE_UNSET:
      mdlerror_fmt(parse_state, "File %s, Line %d: Internal error: Rate is not set", __FILE__, __LINE__);
      return NULL;

    case RATE_CONSTANT:
      pathp->km = rate->forward_rate.v.rate_constant;
      pathp->km_filename = NULL;
      pathp->km_complex = NULL;
      break;

    case RATE_FILE:
      pathp->km = 0.0;
      pathp->km_filename = mdl_find_include_file(parse_state, rate->forward_rate.v.rate_file, parse_state->vol->curr_file);
      free(rate->forward_rate.v.rate_file);
      pathp->km_complex = NULL;
      break;

    case RATE_COMPLEX:
      pathp->km = 0.0;
      pathp->km_filename = NULL;
      pathp->km_complex = rate->forward_rate.v.rate_complex;
      break;

    default: UNHANDLED_CASE(rate->forward_rate.rate_type);
  }

  /* Add the pathway to the list for this reaction */
  if (rate->forward_rate.rate_type == RATE_FILE)
  {
    struct pathway *tpp;
    if (rxnp->pathway_head == NULL)
    {
      rxnp->pathway_head = pathp;
      pathp->next = NULL;
    }
    else  /* Move varying reactions to the end of the list */
    {
      for (tpp = rxnp->pathway_head;
            tpp->next != NULL && tpp->next->km_filename==NULL;
            tpp = tpp->next) {}
      pathp->next = tpp->next;
      tpp->next = pathp;
    }
  }
  else
  {
    pathp->next = rxnp->pathway_head;
    rxnp->pathway_head = pathp;
  }

  /* If we're doing 3D releases, set up array so we can release reversibly */
  if (parse_state->vol->r_step_release == NULL  &&  all_3d  &&  pathp->product_head != NULL)
  {
    parse_state->vol->r_step_release = init_r_step_3d_release(parse_state->vol->radial_subdivisions);
    if (parse_state->vol->r_step_release == NULL)
    {
      mdlerror(parse_state,"Out of memory building r_step array.");
      return NULL;
    }
  }

  /* If the vacancy search distance is zero and this reaction produces more
   * grid molecules than it comsumes, it can never succeed, except if it is a
   * volume molecule hitting the surface and producing a single grid molecule.
   * Fail with an error message.
   */
  if ((parse_state->vol->vacancy_search_dist2 == 0)  &&
      (num_surf_products > num_grid_mols))
  {
    /* The case with one volume molecule reacting with the surface and
     * producing one grid molecule is okay.
     */
    if (num_grid_mols == 0 && num_vol_mols == 1 && num_surf_products == 1)
    {
      /* do nothing */
    }
    else
    {
      mdlerror(parse_state, "Error: number of surface products exceeds number of surface reactants, but VACANCY_SEARCH_DISTANCE is not specified or set to zero.");
      return NULL;
    }
  }

  /* A non-reversible reaction may not specify a reverse reaction rate */
  if (rate->backward_rate.rate_type != RATE_UNSET && ! bidirectional)
  {
    mdlerror(parse_state, "Reverse rate specified but the reaction isn't reversible.");
    return NULL;
  }

  /* Create reverse reaction if we need to */
  if (bidirectional)
  {
    /* A bidirectional reaction must specify a reverse rate */
    if (rate->backward_rate.rate_type == RATE_UNSET)
    {
      mdlerror(parse_state, "Reversible reaction indicated but no reverse rate supplied.");
      return NULL;
    }

    /* if "surface_class" is present on the reactant side of the reaction copy
     * it to the product side of the reaction.
     *
     * Reversible reaction of the type:
     *    A' @ surf' <---> C''[>r1,<r2]
     *
     * is equivalent now to the two reactions:
     *    A' @ surf' ---> C'' [r1]
     *    C'' @ surf' ----> A' [r2]
     *
     * Reversible reaction of the type:
     *    A' + B' @ surf' <---> C'' + D'' [>r1,<r2]
     *
     * is equivalent now to the two reactions:
     *    A' + B @ surf' ---> C'' + D'' [r1]
     *    C'' + D'' @ surf' ----> A' + B' [r2]
     */
    if (surface != -1  &&  surface != catalytic)
    {
      struct product *prodp;
      prodp = (struct product *) CHECKED_MEM_GET(parse_state->prod_mem,
                                                  "reaction product");
      if (prodp == NULL)
      {
        mem_put(parse_state->prod_mem, prodp);
        return NULL;
      }

      switch (surface)
      {
        case 1:
          prodp->prod = pathp->reactant2;
          prodp->orientation = pathp->orientation2;
          break;

        case 2:
          prodp->prod = pathp->reactant3;
          prodp->orientation = pathp->orientation3;
          break;

        case 0:
        default:
          mcell_internal_error("Surface appears in invalid reactant slot in reaction (%d).", surface);
          break;
      }
      prodp->next = pathp->product_head;
      pathp->product_head = prodp;
    }

    /* Invert the current reaction pathway */
    if (invert_current_reaction_pathway(parse_state, pathp, &rate->backward_rate))
      return NULL;
  }

  return rxnp;
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
struct rxn *mdl_assemble_surface_reaction(struct mdlparse_vars *parse_state,
                                          int reaction_type,
                                          struct species *surface_class,
                                          struct sym_table *reactant_sym,
                                          short orient)
{
  struct species *reactant = (struct species *) reactant_sym->value;
  struct product *prodp;
  struct rxn *rxnp;
  struct pathway *pathp;
  struct name_orient *no;

  /* Make sure the other reactant isn't a surface */
  if (reactant->flags == IS_SURFACE)
  {
    mdlerror_fmt(parse_state,
                 "Illegal reaction between two surfaces in surface reaction: %s -%s-> ...",
                 reactant_sym->name,
                 surface_class->sym->name);
    return NULL;
  }

  /* Build reaction name */
  char *rx_name = concat_rx_name(parse_state, surface_class->sym->name, 0, reactant_sym->name, 0);
  if (rx_name == NULL)
  {
    mdlerror_fmt(parse_state,
                 "Out of memory while parsing surface reaction: %s -%s-> ...",
                 surface_class->sym->name,
                 reactant_sym->name);
    return NULL;
  }

  /* Find or create reaction */
  struct sym_table *reaction_sym;
  if ((reaction_sym = retrieve_sym(rx_name, parse_state->vol->rxn_sym_table)) != NULL)
  {
    /* do nothing */
  }
  else if ((reaction_sym = store_sym(rx_name, RX, parse_state->vol->rxn_sym_table, NULL)) == NULL)
  {
    free(rx_name);
    mdlerror_fmt(parse_state,
                 "Out of memory while creating surface reaction: %s -%s-> ...",
                 reactant_sym->name,
                 surface_class->sym->name);
    return NULL;
  }
  free(rx_name);

  /* Create pathway */
  if ((pathp = (struct pathway *) mem_get(parse_state->path_mem)) == NULL)
  {
    mdlerror_fmt(parse_state,
                 "Out of memory while creating surface reaction: %s -%s-> ...",
                 reactant_sym->name,
                 surface_class->sym->name);
    return NULL;
  }

  rxnp = (struct rxn *)reaction_sym->value;
  rxnp->n_reactants = 2;
  ++ rxnp->n_pathways;
  pathp->pathname = NULL;
  pathp->reactant1 = surface_class;
  pathp->reactant2 = (struct species *) reactant_sym->value;
  pathp->reactant3 = NULL;
  pathp->is_complex[0] = pathp->is_complex[1] = pathp->is_complex[2] = 0;
  pathp->km = GIGANTIC;
  pathp->km_filename = NULL;
  pathp->km_complex = NULL;
  pathp->prod_signature = NULL;
  pathp->flags=0;

  pathp->orientation1 = 1;
  pathp->orientation3 = 0;
  if (orient == 0)
  {
    pathp->orientation2 = 0;
  }
  else
  {
    pathp->orientation2 = (orient < 0) ? -1 : 1;
  }

  no = CHECKED_MALLOC_STRUCT(struct name_orient, "struct name_orient");
  no->name = CHECKED_STRDUP(reactant->sym->name, "reactant name");
  if (orient == 0)
  {
    no->orient = 0;
  }
  else 
  {
    no->orient = (orient < 0) ? -1 : 1;
  }

  switch (reaction_type)
  {
    case RFLCT:
      if ((prodp=(struct product *)mem_get(parse_state->prod_mem))==NULL)
      {
        mdlerror_fmt(parse_state,
                     "Out of memory while creating surface reaction: %s -%s-> ...",
                     reactant_sym->name, surface_class->sym->name);
        return NULL;
      }
      pathp->flags |= PATHW_REFLEC;
      prodp->prod = pathp->reactant2;
      prodp->orientation = 1;
      prodp->next = NULL;
      pathp->product_head = prodp;
      if (pathp->product_head != NULL)
      {
        pathp->prod_signature = create_prod_signature(&pathp->product_head);
        if (pathp->prod_signature == NULL)
        {
          mdlerror(parse_state, "Error creating 'prod_signature' field for the reaction pathway.");
          return NULL;
        }
      }
      if (surface_class->refl_mols == NULL)
      {
         no->next = NULL;
         surface_class->refl_mols = no;
      }
      else 
      {
         no->next = surface_class->refl_mols;
         surface_class->refl_mols = no;
      }

      break;
    case TRANSP:
      if ((prodp = (struct product *)mem_get(parse_state->prod_mem))==NULL)
      {
        mdlerror_fmt(parse_state,
                     "Out of memory while creating surface reaction: %s -%s-> ...",
                     reactant_sym->name, surface_class->sym->name);
        return NULL;
      }
      pathp->flags |= PATHW_TRANSP;
      prodp->prod = pathp->reactant2;
      prodp->orientation = -1;
      prodp->next = NULL;
      pathp->product_head = prodp;
      if (pathp->product_head != NULL)
      {
        pathp->prod_signature = create_prod_signature(&pathp->product_head);
        if (pathp->prod_signature == NULL)
        {
          mdlerror(parse_state, "Error creating 'prod_signature' field for the reaction pathway.");
          return NULL;
        }
      }
      if (surface_class->transp_mols == NULL)
      {
         no->next = NULL;
         surface_class->transp_mols = no;
      }
      else 
      {
         no->next = surface_class->transp_mols;
         surface_class->transp_mols = no;
      }
      break;
    case SINK:
      pathp->flags |= PATHW_ABSORP;
      pathp->product_head = NULL;
      if (surface_class->absorb_mols == NULL)
      {
         no->next = NULL;
         surface_class->absorb_mols = no;
      }
      else 
      {
         no->next = surface_class->absorb_mols;
         surface_class->absorb_mols = no;
      }
      break;
    default:
      mdlerror(parse_state, "Unknown special surface type.");
      return NULL;
      break;
  }

  pathp->next = rxnp->pathway_head;
  rxnp->pathway_head = pathp;

#ifdef DEBUG
  no_printf("Surface reaction defined:\n");
  no_printf("  %s[%d] -%s[%d]->",
            rxnp->pathway_head->reactant2->sym->name,
            rxnp->pathway_head->orientation2,
            rxnp->pathway_head->reactant1->sym->name,
            rxnp->pathway_head->orientation1);
  for (prodp = rxnp->pathway_head->product_head;
       prodp != NULL;
       prodp = prodp->next)
  {
    if (prodp != rxnp->pathway_head->product_head)
      no_printf(" +");
    no_printf(" %s[%d]", prodp->prod->sym->name, prodp->orientation);
  }
  no_printf(" [%.9g]\n",rxnp->pathway_head->km);
#endif
  return rxnp;
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
struct rxn *mdl_assemble_concentration_clamp_reaction(struct mdlparse_vars *parse_state,
                                                      struct species *surface_class,
                                                      struct sym_table *mol_sym,
                                                      short orient,
                                                      double conc)
{
  struct rxn *rxnp;
  struct pathway *pathp;
  struct sym_table *stp3;
  struct species *specp = (struct species *) mol_sym->value;
  struct name_orient *no;

  if (specp->flags == IS_SURFACE)
  {
    mdlerror_fmt(parse_state,
                 "Illegal reaction between two surfaces in surface reaction: %s -%s-> ...",
                 mol_sym->name, surface_class->sym->name);
    return NULL;
  }
  if (specp->flags & ON_GRID)
  {
    mdlerror(parse_state, "Concentration clamp does not work on surface molecules.");
    return NULL;
  }
  if (specp->flags&NOT_FREE || specp->D <= 0.0)
  {
    mdlerror(parse_state, "Concentration clamp must be applied to molecule diffusing in 3D");
    return NULL;
  }
  if (conc < 0)
  {
    mdlerror(parse_state, "Concentration can only be clamped to positive values.");
    return NULL;
  }

  char *rx_name = concat_rx_name(parse_state, surface_class->sym->name, 0, mol_sym->name, 0);
  if (rx_name == NULL)
  {
    mdlerror_fmt(parse_state,
                 "Memory allocation error: %s -%s-> ...",
                 surface_class->sym->name, mol_sym->name);
    return NULL;
  }
  if ((stp3=retrieve_sym(rx_name, parse_state->vol->rxn_sym_table)) !=NULL)
  {
    /* do nothing */
  }
  else if ((stp3=store_sym(rx_name,RX,parse_state->vol->rxn_sym_table, NULL)) ==NULL)
 {
    free(rx_name);
    mdlerror_fmt(parse_state,
                 "Cannot store surface reaction: %s -%s-> ...",
                 mol_sym->name, surface_class->sym->name);
    return NULL;
  }
  free(rx_name);
  if ((pathp=(struct pathway *)mem_get(parse_state->path_mem))==NULL)
  {
    mdlerror_fmt(parse_state,
                 "Cannot store surface reaction: %s -%s-> ...",
                 mol_sym->name, surface_class->sym->name);
    return NULL;
  }
  rxnp = (struct rxn *)stp3->value;
  rxnp->n_reactants = 2;
  ++ rxnp->n_pathways;
  pathp->pathname = NULL;
  pathp->reactant1 = surface_class;
  pathp->reactant2 = (struct species *) mol_sym->value;
  pathp->reactant3 = NULL;
  pathp->is_complex[0] = pathp->is_complex[1] = pathp->is_complex[2] = 0;
  pathp->flags = 0;

  pathp->flags |= PATHW_CLAMP_CONC;

  pathp->km = conc;
  pathp->km_filename = NULL;
  pathp->km_complex = NULL;

  pathp->orientation1 = 1;
  pathp->orientation3 = 0;
  if (orient == 0)
  {
    pathp->orientation2 = 0;
  }
  else
  {
    pathp->orientation2 = (orient < 0) ? -1 : 1;
  }

  pathp->product_head = NULL;
  pathp->prod_signature = NULL;

  pathp->next = rxnp->pathway_head;
  rxnp->pathway_head = pathp;

  no = CHECKED_MALLOC_STRUCT(struct name_orient, "struct name_orient");
  no->name = CHECKED_STRDUP(mol_sym->name, "molecule name");
  no->orient = pathp->orientation2;

  if (surface_class->clamp_conc_mols == NULL)
  {
    no->next = NULL;
    surface_class->clamp_conc_mols = no;
  }
  else 
  {
    no->next = surface_class->clamp_conc_mols;
    surface_class->clamp_conc_mols = no;
  }


  return rxnp;
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
                             struct sym_table *symp)
{
  struct species *specp = (struct species *) symp->value;
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
     symp: symbol for the surface class
 Out: 0 on success, 1 on failure
**************************************************************************/
void mdl_finish_surface_class(struct mdlparse_vars *parse_state,
                              struct sym_table *symp)
{
  parse_state->current_surface_class = NULL;
}

/**************************************************************************
 mdl_new_effector_data:
    Create a new effector data for surface molecule initialization.

 In: parse_state: parser state
     eff_info:
     quant: the amount of surface molecules to release
 Out: 0 on success, 1 on failure
**************************************************************************/
struct eff_dat *mdl_new_effector_data(struct mdlparse_vars *parse_state,
                                      struct species_opt_orient *eff_info,
                                      double quant)
{
  struct eff_dat *effdp;
  struct species *specp = (struct species *) eff_info->mol_type->value;
  if (! (specp->flags & ON_GRID))
  {
    mdlerror_fmt(parse_state, "Cannot initialize surface with non-surface molecule '%s'", specp->sym->name);
    return NULL;
  }
  else if (mdl_check_valid_molecule_release(parse_state, eff_info))
  {
    return NULL;
  }

  if ((effdp = CHECKED_MALLOC_STRUCT(struct eff_dat, "surface molecule data")) == NULL)
    return NULL;

  effdp->next = NULL;
  effdp->eff = specp;
  effdp->quantity_type = 0;
  effdp->quantity = quant;
  effdp->orientation = eff_info->orient;
  return effdp;
}

/*************************************************************************
 * Macromolecules
 *************************************************************************/

/*************************************************************************
 macro_relation_index:
    Find a relation by name in the relations table.

    In:  struct subunit_relation const *relations - array of subunit relations
         int num_Relations - length of subunit relations array
         char const *name - the relation name to find
    Out: 0-based index into relations, or -1 if relation not found
*************************************************************************/
static int macro_relation_index(struct subunit_relation const *relations,
                                int num_relations,
                                char const *name)
{
  int i;
  for (i=0; i<num_relations; ++i)
    if (! strcmp(relations[i].name, name))
      return i;
  return -1;
}

/*************************************************************************
 macro_convert_clause:
    Fill in a single row in the rate rule table.

    In: struct macro_rate_clause *clauses - linked list of conditions to match
                                            for this row in the rate rule table
        struct subunit_relation const *relations - array of subunit relations
        int num_relations - length of relations array
        struct species **nptr - pointer to "neighbors" table containing species
                                            for this clause in the rate rule
        int *iptr - pointer to "invert" table which is 1 if this condition
                                            should be inverted (i.e. != instead
                                            of ==)
        signed char *optr - NULL, or pointer to "orient" table, which contains
                            -1, 0, or 1 for each condition in this clause
    Out: Nothing.  nptr, iptr, and optr table row are filled in
*************************************************************************/
static void macro_convert_clause(struct macro_rate_clause *clauses,
                                 struct subunit_relation const *relations,
                                 int num_relations,
                                 struct species **nptr,
                                 int *iptr,
                                 signed char *optr)
{
  for (; clauses != NULL; clauses = clauses->next)
  {
    int relation_index = macro_relation_index(relations, num_relations, clauses->name);
    if (relation_index == -1)
    {
      /* XXX: Sanity check - signal an error? */
      /* shouldn't occur because we already checked the clauses earlier */
      nptr[relation_index] = NULL;
      iptr[relation_index] = 0;
      if (optr) optr[relation_index] = 0;
    }
    else
    {
      nptr[relation_index] = clauses->species;
      iptr[relation_index] = clauses->invert;
      if (optr)
      {
        if (clauses->orient > 0)
          optr[relation_index] = 1;
        else if (clauses->orient < 0)
          optr[relation_index] = -1;
        else
          optr[relation_index] = 0;
      }
    }
  }
}

/*************************************************************************
 macro_build_rate_table:
    Convert the parse-time data structures to a run-time rate table.

 In:  cr: the complex rate structure to fill in
      relations: array of subunit relations
      num_relations: length of relations array
      rules: reverse-ordered list of rules to put into table
      is_surface: flag indicating whether to include orientations in the table
 Out: 0 on success, 1 if allocation fails.  'cr' structure is filled in.
*************************************************************************/
static int macro_build_rate_table(struct complex_rate *cr,
                                  struct subunit_relation const *relations,
                                  int num_relations,
                                  struct macro_rate_rule *rules,
                                  int is_surface)
{

  /* Count the rules */
  struct macro_rate_rule *rules_temp;
  int rule_count = 0, rule_index;
  for (rules_temp = rules; rules_temp != NULL; rules_temp = rules_temp->next)
    ++ rule_count;

  /* Allocate and zero the tables */
  cr->num_rules = rule_count;
  cr->num_neighbors = num_relations;
  cr->neighbors = NULL;
  cr->invert = NULL;
  cr->rates = NULL;
  cr->orientations = NULL;
  cr->neighbors = CHECKED_MALLOC_ARRAY(struct species *,
                                        rule_count * num_relations,
                                        "macromolecule rate rule neighbors table");
  if (cr->neighbors == NULL) return 1;
  cr->invert    = CHECKED_MALLOC_ARRAY(int,
                                        rule_count * num_relations,
                                        "macromolecule rate rule invert table");
  if (cr->invert == NULL) return 1;
  cr->rates     = CHECKED_MALLOC_ARRAY(double,
                                        rule_count * num_relations,
                                        "macromolecule rate rule rates");
  if (cr->rates == NULL) return 1;
  if (is_surface)
  {
    cr->orientations = CHECKED_MALLOC_ARRAY(signed char,
                                             rule_count * num_relations,
                                             "macromolecule rate rule orientations");
    if (cr->orientations == NULL) return 1;
  }
  memset(cr->neighbors, 0, rule_count * num_relations * sizeof(struct species *));
  memset(cr->invert, 0, rule_count * num_relations * sizeof(int));
  if (is_surface)
    memset(cr->orientations, 0, rule_count * num_relations * sizeof(signed char));

  /* Set up table pointers. */
  struct species **nptr = cr->neighbors + num_relations * (rule_count - 1);
  int *iptr = cr->invert + num_relations * (rule_count - 1);
  signed char *optr = NULL;
  if (is_surface) optr = cr->orientations + num_relations * (rule_count - 1);

  /* Now, dump all of the rules into the table */
  for (rule_index = rule_count - 1; rules != NULL; --rule_index, rules = rules->next)
  {
    /* Convert this row in the table */
    cr->rates[rule_index] = rules->rate;
    macro_convert_clause(rules->clauses, relations, num_relations, nptr, iptr, optr);

    /* Advance to next row in table */
    iptr -= num_relations;
    nptr -= num_relations;
    if (optr) optr -= num_relations;
  }

  return 0;
}

/*************************************************************************
 macro_free_runtime_rate_tables:
    Free the rate tables stored in a species.

 In:  cs: species whose rates we mean to clea
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_runtime_rate_tables(struct complex_species *cs)
{
  while (cs->rates != NULL)
  {
    struct complex_rate *cr = cs->rates;
    cs->rates = cs->rates->next;
    free(cr->neighbors);
    free(cr->invert);
    free(cr->rates);
    if (cr->orientations) free(cr->orientations);
    free(cr);
  }
}

/*************************************************************************
 macro_build_rate_tables:
    Convert the parse-time data structures to run-time rate tables.

 In:  parse_state: parser state
      cs: the species to receive the tables
      rates: a linked list of rate rulesets to be converted to run-time format
 Out: 0 on success, 1 if allocation fails.  'cs' structure is filled in.
*************************************************************************/
static int macro_build_rate_tables(struct mdlparse_vars *parse_state,
                                   struct complex_species *cs,
                                   struct macro_rate_ruleset *rates)
{
  /* Check whether we need to include orientation information */
  int is_surface = 1;
  if ((cs->base.flags & NOT_FREE) == 0) is_surface = 0;

  /* Now, convert each rate table to run-time usable format. */
  for (; rates != NULL; rates = rates->next)
  {
    struct complex_rate *cr = CHECKED_MALLOC_STRUCT(struct complex_rate,
                                                     "macromolecular rate table");
    if (cr == NULL)
      return 1;

    cr->next = cs->rates;
    cr->name = rates->name;
    cs->rates = cr;
    if (macro_build_rate_table(cr,
                               cs->relations,
                               cs->num_relations,
                               rates->rules,
                               is_surface))
      return 1;
  }

  return 0;
}

/*************************************************************************
 macro_new_subunit_assignment:
    Allocate a new subunit assignment data structure for a macromolecule.  A
    subunit assignment is the combination of:
      - a location specification, which is a cartesian product of contiguous
        ranges, one item for each dimension in the topology
      - a subunit species for this portion of the complex
      - an orientation for this portion of the complex, if the complex is a
        surface complex
    This data structure is used only at parse-time and is transferred to an
    array and freed before the simulation starts.

 In:  where: location specification
      mol: subunit type specification
      orient: subunit orientation specification
 Out: Freshly allocated structure, or NULL if we failed to allocate
*************************************************************************/
static struct macro_subunit_assignment *macro_new_subunit_assignment(struct macro_subunit_spec *where,
                                                                     struct species *mol,
                                                                     short orient)
{

  struct macro_subunit_assignment *a;
  a = CHECKED_MALLOC_STRUCT(struct macro_subunit_assignment, "macromolecule subunit assignments");
  if (a == NULL)
    return NULL;
  a->next = NULL;
  a->head = where;
  a->what = mol;
  a->orient = orient;
  return a;
}

/*************************************************************************
 macro_linear_array_index_to_string:
    This converts a linear subunit array index to an n-dimensional subunit
    index, in the form of a string.  This is useful for providing useful error
    messages.

 In:  linear_index: the index to convert
      topo: the macromol topology
 Out: a char * allocated to hold the array index
*************************************************************************/
static char *macro_linear_array_index_to_string(int linear_index,
                                                struct macro_topology *topo)
{

  int subunit_index_rem = linear_index;
  int dim_index;
  int coords[topo->head->dimensionality];

  /* Convert linear-coord to an ordered array of coordinates */
  for (dim_index = topo->head->dimensionality - 1; dim_index >= 0; -- dim_index)
    coords[dim_index] = 1 + (subunit_index_rem % topo->head->dimensions[dim_index]);

  /* Turn the array into a string */
  char *message = NULL;
  for (dim_index = 0; dim_index < topo->head->dimensionality; ++ dim_index)
  {
    char *newmessage;
    if (dim_index < topo->head->dimensionality - 1)
    {
      subunit_index_rem /= topo->head->dimensions[dim_index];
      if (dim_index == 0)
        newmessage = CHECKED_SPRINTF("[%d, ", coords[dim_index]);
      else
        newmessage = CHECKED_SPRINTF("%s%d, ", message, coords[dim_index]);
    }
    else
      newmessage = CHECKED_SPRINTF("%s%d]", message, coords[dim_index]);
    if (message != NULL) free(message);
    message = newmessage;
  }

  return message;
}


/*************************************************************************
 mdl_assemble_complex_relation_state:
    Allocate a new relation state structure for a rule table (presently, either
    a rate table, or a subunit counting table).

 In:  parse_state: the parser state for error logging
      rel_idx: the index of the relation in the relation table
      invert: the invert flag for the relation
      state: the referenced species with optional orientation
 Out: A freshly allocated relation state struct, or NULL if allocation fails.
*************************************************************************/
struct macro_relation_state *mdl_assemble_complex_relation_state(struct mdlparse_vars *parse_state,
                                                                 int rel_idx,
                                                                 int invert,
                                                                 struct species_opt_orient *state)
{
  struct species *mol = (struct species *) state->mol_type->value;
  struct macro_relation_state *relstate;

  if ((parse_state->current_complex->base.flags & NOT_FREE) == 0)
  {
    if (state->orient_set)
    {
      if (parse_state->vol->notify->useless_vol_orient==WARN_ERROR)
      {
        mdlerror_fmt(parse_state, "Error: orientation specified for subunit relation of volume complex '%s' in count statement", parse_state->current_complex->base.sym->name);
        return NULL;
      }
      else if (parse_state->vol->notify->useless_vol_orient==WARN_WARN)
        mdlerror_fmt(parse_state, "Warning: orientation specified for subunit relation of volume complex '%s' in count statement", parse_state->current_complex->base.sym->name);
    }
  }

  relstate = CHECKED_MALLOC_STRUCT(struct macro_relation_state, "macromolecule subunit relation state");
  if (relstate == NULL)
    return NULL;

  relstate->next = NULL;
  relstate->relation = rel_idx;
  relstate->invert = invert;
  relstate->mol = mol;
  relstate->orient = state->orient;
  return relstate;
}

/*************************************************************************
 macro_check_subunit_spec:
    Check the subunit spec for correctness.  If it is incorrect, a message will
    be logged, and 1 will be returned.

 In:  parse_state: the parser state for error logging
      topo: the macromol topology
      spec: the spec to check
 Out: 0 if the spec is well-formed, 1 if not.  On error, the details are written
      to the error log.
*************************************************************************/
static int macro_check_subunit_spec(struct mdlparse_vars *parse_state,
                                    struct macro_topology *topo,
                                    struct macro_subunit_spec *spec)
{
  /* N.B. The elements of 'spec' are in reverse order */
  struct macro_topology_element *topo_el = topo->head;
  int dim_index;

  /* Check each dimension of the array indices against the array dimensionality */
  for (dim_index = topo_el->dimensionality; dim_index > 0; -- dim_index)
  {
    /* Did we run out of array indices early? */
    if (spec == NULL)
    {
      mdlerror_fmt(parse_state,
                   "In subunit assignment, %d array indices must be specified (%d %s specified)",
                   topo_el->dimensionality,
                   topo_el->dimensionality - dim_index,
                   (topo_el->dimensionality - dim_index) == 1 ? "was" : "were");
      return 1;
    }

    /* Did the user specify the range coordinates in the wrong order? */
    if (spec->from > spec->to)
    {
      mdlerror_fmt(parse_state,
                   "In SUBUNIT[..., %d:%d, ...], the range must be increasing (i.e. %d:%d)",
                   spec->from,
                   spec->to,
                   spec->to,
                   spec->from);
      return 1;
    }

    /* Did the user specify indices <= 0? */
    if (spec->from <= 0  ||  spec->to <= 0)
    {
      mdlerror_fmt(parse_state,
                   "In subunit assignment, all array indices must be > 0");
      return 1;
    }

    /* Did the user specify indices > the corresponding array dimension?  */
    if (spec->to > topo_el->dimensions[dim_index - 1]  ||  spec->from > topo_el->dimensions[dim_index - 1])
    {
      if (spec->to == spec->from)
        mdlerror_fmt(parse_state, "In subunit assignment, array index %d must be in the range 1..%d (%d was specified)",
                     dim_index,
                     topo_el->dimensions[dim_index - 1],
                     spec->to);
      else
        mdlerror_fmt(parse_state, "In subunit assignment, array index %d must be entirely in 1..%d (%d:%d was specified)",
                     dim_index,
                     topo_el->dimensions[dim_index - 1],
                     spec->from,
                     spec->to);
      return 1;
    }

    spec = spec->next;
  }

  /* Did we have array indices left over afterwards? */
  if (spec != NULL)
  {
    while (spec != NULL)
    {
      ++ dim_index;
      spec = spec->next;
    }

    mdlerror_fmt(parse_state, "In subunit assignment, only %d array indices may be specified (%d were specified)",
                 topo_el->dimensionality,
                 dim_index + topo_el->dimensionality);
    return 1;
  }

  return 0;
}

/*************************************************************************
 macro_free_subunit_specs:
    Free a linked list of subunit coordinate specifications.  This will free
    the memory associated with the structures.

 In:  specs: the list to free
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_subunit_specs(struct macro_subunit_spec *specs)
{
  while (specs != NULL)
  {
    struct macro_subunit_spec *s = specs;
    specs = specs->next;
    free(s);
  }
}

/*************************************************************************
 mdl_assemble_complex_subunit_assignment:
    Assemble a subunit assignment for one or more subunits within a complex.
    These will be the species with which these subunits will be initialized
    when a complex of the containing type is created.

 In:  parse_state: parser state
      su:   specification of the subunit/subunits
      spec: species and (opt.) orientation
 Out: the subunit assignment, or NULL if an error occurred
*************************************************************************/
struct macro_subunit_assignment *mdl_assemble_complex_subunit_assignment(struct mdlparse_vars *parse_state,
                                                                         struct macro_subunit_spec *su,
                                                                         struct species_opt_orient *spec)
{
  struct macro_subunit_assignment *the_assignment;
  struct species *sp = (struct species *) spec->mol_type->value;
  if (sp->flags & ON_GRID)
  {
    if (parse_state->complex_type == 0)
      parse_state->complex_type = TYPE_GRID;
    else if (parse_state->complex_type != TYPE_GRID)
    {
      mdlerror_fmt(parse_state, "Subunit type '%s' is not a surface molecule, but the complex has other subunits which are surface molecules.", sp->sym->name);
      macro_free_subunit_specs(su);
      return NULL;
    }
  }
  else if (! (sp->flags & NOT_FREE))
  {
    if (parse_state->complex_type == 0)
      parse_state->complex_type = TYPE_3D;
    else if (parse_state->complex_type != TYPE_3D)
    {
      mdlerror_fmt(parse_state, "Subunit type '%s' is not a volume molecule, but the complex has other subunits which are volume molecules.", sp->sym->name);
      macro_free_subunit_specs(su);
      return NULL;
    }
  }
  else
  {
    mdlerror_fmt(parse_state, "Subunit type '%s' is not a molecule.", sp->sym->name);
    macro_free_subunit_specs(su);
    return NULL;
  }

  if (macro_check_subunit_spec(parse_state, parse_state->complex_topo, su))
  {
    macro_free_subunit_specs(su);
    return NULL;
  }

  if ((the_assignment = macro_new_subunit_assignment(su, sp, spec->orient)) == NULL)
    macro_free_subunit_specs(su);

  return the_assignment;
}

/*************************************************************************
 macro_check_subunit_geometry:
    Check a geometry assignment for a subunit.  In reality, this only needs to
    check that the subunit specification refers to an actual element within the
    topology of the complex.  If it does not, an error message will be logged,
    and 1 returned.

 In:  parse_state: the parser state for error logging
      topo: the macromol topology
      indices: the subunit indices to check
 Out: 0 if the subunit has well-formed geometry, 1 if not.  On error, the
      details are written to the error log.
*************************************************************************/
static int macro_check_subunit_geometry(struct mdlparse_vars *parse_state,
                                        struct macro_topology *topo,
                                        struct num_expr_list *indices)
{
  struct macro_topology_element *topo_el = topo->head;
  int dim_index;

  /* Check each dimension of the array indices against the array dimensionality */
  for (dim_index = 0; dim_index < topo_el->dimensionality; ++ dim_index)
  {
    /* Did we run out of array indices early? */
    if (indices == NULL)
    {
      mdlerror_fmt(parse_state,
                   "In subunit location assignment, %d array indices must be specified (%d %s specified)",
                   topo_el->dimensionality,
                   dim_index,
                   (dim_index == 1) ? "was" : "were");
      return 1;
    }

    /* Is the value out of range? */
    int value = (int) indices->value;
    if (value <= 0  ||  value > topo_el->dimensions[dim_index])
    {
      mdlerror_fmt(parse_state, "In subunit location assignment, array index %d must be in the range 1..%d (%d was specified)",
                   dim_index + 1,
                   topo_el->dimensions[dim_index],
                   value);
      return 1;
    }

    indices = indices->next;
  }

  /* Did we have array indices left over afterwards? */
  if (indices != NULL)
  {
    while (indices != NULL)
    {
      ++ dim_index;
      indices = indices->next;
    }

    mdlerror_fmt(parse_state, "In subunit location assignment, only %d array indices may be specified (%d were specified)",
                 topo_el->dimensionality,
                 dim_index + 1);
    return 1;
  }

  return 0;
}

/*************************************************************************
 macro_new_geometry:
    Allocate a new geometry data structure for a macromolecule.  Geometry is in
    the form of 3d offsets relative to a fictional "center" point of the
    complex.  This data structure is used only at parse-time and is transferred
    to an array and freed before the simulation starts.

 In:  parse_state: parser state for error reporting
      coord_indices: coordinates for which to specify location
      where: the location for this subunit
 Out: Freshly allocated structure, or NULL if we failed to allocate
*************************************************************************/
static struct macro_geometry *macro_new_geometry(struct mdlparse_vars *parse_state,
                                                 struct num_expr_list *coord_indices,
                                                 struct vector3 *where)
{
  struct macro_geometry *mg = CHECKED_MALLOC_STRUCT(struct macro_geometry,
                                                     "macromolecule geometry");
  if (mg == NULL)
    return NULL;

  mg->next = NULL;
  mg->index = coord_indices;
  mg->location = *where;
  mg->location.x *= parse_state->vol->r_length_unit;
  mg->location.y *= parse_state->vol->r_length_unit;
  mg->location.z *= parse_state->vol->r_length_unit;
  return mg;
}

/*************************************************************************
 mdl_assemble_complex_geometry:
    Assemble a complex geometry structure.

 In:  parse_state: parser state
      topo: topology of this complex
      coords: coordinates of the subunit within the complex
      pos: spatial position of the subunit
 Out: a valid macro_geometry structure, or NULL if there was an error.
*************************************************************************/
struct macro_geometry *mdl_assemble_complex_geometry(struct mdlparse_vars *parse_state,
                                                     struct macro_topology *topo,
                                                     struct num_expr_list_head *coords,
                                                     struct vector3 *pos)
{
  if (! macro_check_subunit_geometry(parse_state, topo, coords->value_head))
  {
    struct macro_geometry *geom = macro_new_geometry(parse_state, coords->value_head, pos);
    if (geom)
    {
      /* Do not free coords -- it's been added to the geom data structure */
      free(pos);
      return geom;
    }
  }

  mdl_free_numeric_list(coords->value_head);
  free(pos);
  return NULL;
}

/*************************************************************************
 macro_check_complex_geometry:
    Check geometry assignments for an entire complex.  This checks that each
    subunit's position is specified exactly once.

 In:  parse_state: the parser state for error logging
      topo: the macromol topology
      geom: the geometry to check
 Out: 0 if the complex has well-formed geometry, 1 if not.  On error, the
      details are written to the error log.
*************************************************************************/
static int macro_check_complex_geometry(struct mdlparse_vars *parse_state,
                                        struct macro_topology *topo,
                                        struct macro_geometry *geom)
{
  /* The "subunits" array contains flags indicating whether we've seen a
   * particular subunit
   */
  int subunits[topo->total_subunits];
  int i;

  /* Initialize all subunit flags to "unseen" */
  for (i = 0; i<topo->total_subunits; ++i)
    subunits[i] = 0;

  /* Mark each subunit as "seen" as we find its geometry */
  for (; geom != NULL; geom = geom->next)
  {
    int subunit_index_linear = 0;
    struct num_expr_list *nel;
    int dim_index = 0;
    for (nel = geom->index; nel != NULL; ++ dim_index, nel = nel->next)
    {
      if (dim_index > 0)
        subunit_index_linear *= topo->head->dimensions[dim_index];
      subunit_index_linear += ((int) nel->value) - 1;
    }

    ++ subunits[subunit_index_linear];
  }

  /* Check how many times we specified the location of each subunit */
  for (i = 0; i<topo->total_subunits; ++i)
  {
    /* Check if we've already seen this one */
    if (subunits[i] != 1)
    {
      char *ai = macro_linear_array_index_to_string(i, topo);
      if (subunits[i] == 0)
        mdlerror_fmt(parse_state,
                     "In complex species '%s', no position was specified for subunit %s",
                     parse_state->complex_name,
                     ai);
      else
        mdlerror_fmt(parse_state,
                     "In complex species '%s', %d separate positions were specified for subunit %s",
                     parse_state->complex_name,
                     subunits[i],
                     ai);
      free(ai);
      return 1;
    }
  }

  return 0;
}

/*************************************************************************
 macro_free_geometry:
    Free a linked list of geometry.  This will free the memory associated with
    the geometry structures except for the array of indices, which cannot be
    safely freed, as they may be shared.  (Currently, the parser does not
    duplicate arrays if they are assigned to a variable and multiply
    referenced, and no reference counting is done on the lists.)

 In:  rels: the relationships to free
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_geometry(struct macro_geometry *geom)
{
  while (geom != NULL)
  {
    struct macro_geometry *g = geom;

    /* JW: Memory cannot safely be freed. */
#if 0
    while (g->index != NULL)
    {
      struct num_expr_list *nel = g->index;
      g->index = g->index->next;
      free(nel);
    }
#endif

    geom = geom->next;
    free(g);
  }
}

/*************************************************************************
 mdl_validate_complex_geometry:
    Validate the geometry of a macromolecular complex, freeing it if it is
    invalid.

 In:  parse_state: parser state
      topo: macromolecule topology
      geom: macromolecule geometry
 Out: 0 on success, 1 on failure
*************************************************************************/
int mdl_validate_complex_geometry(struct mdlparse_vars *parse_state,
                                  struct macro_topology *topo,
                                  struct macro_geometry *geom)
{
  if (macro_check_complex_geometry(parse_state, topo, geom))
  {
    macro_free_geometry(geom);
    return 1;
  }
  return 0;
}

/*************************************************************************
 macro_check_complex_relationships:
    Check the relationship definitions for a complex.  At present, this just
    checks that no two relationships have the same name.

 In:  parse_state: the parser state for error logging
      topo: the macromol topology
      rels: linked list of relationships
 Out: 0 if the relationships are uniquely named, 1 if not.  On error, the
      details are written to the error log.
*************************************************************************/
static int macro_check_complex_relationships(struct mdlparse_vars *parse_state,
                                             struct macro_topology *topo,
                                             struct macro_relationship *rels)
{
  if (! rels  ||  ! rels->next)
    return 0;

  /* Check for duplicate names among the relationships */
  /* N.B. This is O(n^2), but unlikely to ever be a bottleneck. */
  /*      If it becomes a bottleneck, hash-based techniques can */
  /*      improve performance here. */
  struct macro_relationship *rel_other;
  for (rel_other = rels->next; rel_other != NULL; rel_other = rel_other->next)
  {
    if (! strcmp(rels->name, rel_other->name))
    {
      mdlerror_fmt(parse_state,
                   "In complex species '%s', relationship '%s' is specified more than once",
                   parse_state->complex_name,
                   rels->name);
      return 1;
    }
  }

  return macro_check_complex_relationships(parse_state, topo, rels->next);
}

/*************************************************************************
 macro_free_relationships:
    Free a linked list of relationships.  This will free the memory associated
    with the relationship structures except for the array of indices, which
    cannot be safely freed, as they may be shared.  (Currently, the parser does
    not duplicate arrays if they are assigned to a variable and multiply
    referenced, and no reference counting is done on the lists.)

 In:  rels: the relationships to free
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_relationships(struct macro_relationship *rels)
{
  while (rels != NULL)
  {
    struct macro_relationship *r = rels;

    /* JW: Memory cannot safely be freed. */
#if 0
    while (r->indices != NULL)
    {
      struct num_expr_list *nel = r->indices;
      r->indices = r->indices->next;
      free(nel);
    }
#endif

    rels = rels->next;
    if (r->name) free(r->name);
    free(r);
  }
}


/*************************************************************************
 mdl_validate_complex_relationships:
    Checks the parsed relationships for errors, freeing them if errors are
    found.  rels are freed if they are invalid.

 In:  parse_state: parser state
      topo: macromol topology
      rels: macromol relationships
 Out: 0 on success, 1 if there is a problem
*************************************************************************/
int mdl_validate_complex_relationships(struct mdlparse_vars *parse_state,
                                       struct macro_topology *topo,
                                       struct macro_relationship *rels)
{
  if (macro_check_complex_relationships(parse_state, topo, rels))
  {
    macro_free_relationships(rels);
    return 1;
  }
  return 0;
}

/*************************************************************************
 macro_check_relationship:
    Check a relationship definition for a complex.  It checks that all of the
    required parameters were set, and that each parameter has the right
    dimensionality.

 In:  parse_state: the parser state for error logging
      topo: the macromol topology
      rel_name: name of the relation
      indices: the coordinate offsets within the topology for this relationship
 Out: 0 if the relationship is well-formed, 1 if not.  On error, the details
      are written to the error log.
*************************************************************************/
static int macro_check_relationship(struct mdlparse_vars *parse_state,
                                    struct macro_topology *topo,
                                    char const *rel_name,
                                    struct num_expr_list *indices)
{
  int dim_index = 0;
  struct num_expr_list *this_index;

  /* Check all offsets in this relationship against the limits provided by the topology */
  for (dim_index = 0, this_index = indices;
       dim_index < topo->head->dimensionality;
       ++ dim_index, this_index = this_index->next)
  {
    /* Check for too few coordinates */
    if (this_index == NULL)
    {
      mdlerror_fmt(parse_state,
                   "In complex species '%s', relationship '%s', exactly %d array indices must be provided (%d %s provided)",
                   parse_state->complex_name,
                   rel_name,
                   topo->head->dimensionality,
                   dim_index,
                   (dim_index == 1) ? "was" : "were");
      return 1;
    }

    /* Check for offsets out of range */
    int su_offset = (int) this_index->value;
    if (su_offset <= - topo->head->dimensions[dim_index]  ||  su_offset >= topo->head->dimensions[dim_index])
    {
      mdlerror_fmt(parse_state,
                   "In complex species '%s', relationship '%s', offsets in array index %d must be between %d and %d, inclusive (value was %d)",
                   parse_state->complex_name,
                   rel_name,
                   dim_index + 1,
                   - (topo->head->dimensions[dim_index] - 1),
                   topo->head->dimensions[dim_index] - 1,
                   su_offset);
      return 1;
    }
  }

  /* Check for too many coordinates */
  if (this_index != NULL)
  {
    for (; this_index != NULL; this_index = this_index->next)
      ++ dim_index;

    mdlerror_fmt(parse_state,
                 "In complex species '%s', relationship '%s', exactly %d array indices must be provided (%d were provided)",
                 parse_state->complex_name,
                 rel_name,
                 topo->head->dimensionality,
                 dim_index);
    return 1;
  }

  return 0;
}

/*************************************************************************
 macro_new_relationship:
    Allocate a new relationship data structure for a macromolecule.  A
    relationship is currently the combination of a name and an offset within
    each dimension of the topology.  The offset is applied modulo the size of
    the dimension under consideration.  Put another way, the topology is
    currently always a discrete n-toroid, and the relationship specifies an
    offset within each of the n dimensions.  This data structure is transferred
    to a table and freed before the simulation starts;

 In:  name: the name of this relationship
      indices: offset within each dimension for this relationship
 Out: Freshly allocated structure, or NULL if we failed to allocate
*************************************************************************/
static struct macro_relationship *macro_new_relationship(char *name,
                                                         struct num_expr_list *indices)
{

  struct macro_relationship *mr;
  mr = CHECKED_MALLOC_STRUCT(struct macro_relationship, "macromolecule subunit relationship");
  if (mr == NULL)
    return NULL;
  mr->next = NULL;
  mr->name = name;
  mr->indices = indices;
  return mr;
}

/*************************************************************************
 mdl_assemble_complex_relationship:
    Assemble a complex relationship.

 In:  parse_state: parser state
      topo: macromolecule topology
      name: name of the relationship
      rel:  relative offset of related subunit from reference subunit in each
            dimension
 Out: the relationship, or NULL if an error occurred
*************************************************************************/
struct macro_relationship *mdl_assemble_complex_relationship(struct mdlparse_vars *parse_state,
                                                             struct macro_topology *topo,
                                                             char *name,
                                                             struct num_expr_list_head *rel)
{
  struct macro_relationship *relate;
  if (macro_check_relationship(parse_state, topo, name, rel->value_head))
  {
    if (! rel->shared)
      mdl_free_numeric_list(rel->value_head);
    free(name);
    return NULL;
  }

  if ((relate = macro_new_relationship(name, rel->value_head)) == NULL)
  {
    if (! rel->shared)
      mdl_free_numeric_list(rel->value_head);
    free(name);
    return NULL;
  }

  return relate;
}

/*************************************************************************
 macro_check_rates:
    Check the rate rule definitions for a complex.  At present, this just
    checks that no two complex rates have the same name.

 In:  parse_state: the parser state for error logging
      topo: the macromol topology
      rulesets: linked list of rulesets
 Out: 0 if the rate rules are uniquely named, 1 if not.  On error, the details
      are written to the error log.
*************************************************************************/
static int macro_check_rates(struct mdlparse_vars *parse_state,
                             struct macro_rate_ruleset *rulesets)
{
  if (! rulesets  ||  ! rulesets->next)
    return 0;

  /* Check for duplicate names among the relationships */
  /* N.B. See comment on O(n^2) algo in macro_check_complex_relationships */
  struct macro_rate_ruleset *rules_other;
  for (rules_other = rulesets->next; rules_other != NULL; rules_other = rules_other->next)
  {
    if (! strcmp(rulesets->name, rules_other->name))
    {
      mdlerror_fmt(parse_state,
                   "In complex species '%s', ruleset '%s' is specified more than once",
                   parse_state->complex_name,
                   rulesets->name);
      return 1;
    }
  }

  return macro_check_rates(parse_state, rulesets->next);
}

/*************************************************************************
 macro_free_rate_clauses:
    Free a linked list of rate clauses.  This will free the memory associated
    with the rate clause structures.

 In:  clauses: the clauses to free
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_rate_clauses(struct macro_rate_clause *clauses)
{
  while (clauses != NULL)
  {
    struct macro_rate_clause *cl = clauses;
    clauses = clauses->next;
    free(cl->name);
    free(cl);
  }
}

/*************************************************************************
 macro_free_rate_rules:
    Free a linked list of rate rules.  This will free the memory associated
    with the rate rule structures.

 In:  rules: the rules to free
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_rate_rules(struct macro_rate_rule *rules)
{
  while (rules != NULL)
  {
    struct macro_rate_rule *ru = rules;
    macro_free_rate_clauses(ru->clauses);
    rules = rules->next;
    free(ru);
  }
}


/*************************************************************************
 macro_free_rates:
    Free a linked list of rate rulesets.  This will free the memory associated
    with the rate ruleset structures.

 In:  rates: the rulesets to free
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_rates(struct macro_rate_ruleset *rates)
{
  while (rates != NULL)
  {
    struct macro_rate_ruleset *r = rates;
    macro_free_rate_rules(r->rules);
    rates = rates->next;
    free(r);
  }
}


/*************************************************************************
 mdl_validate_complex_rates:
    Validate the rate rules, freeing them if invalid.

 In:  parse_state: parser state
      rates: rate rules
 Out: 0 on success, 1 if there is a problem
*************************************************************************/
int mdl_validate_complex_rates(struct mdlparse_vars *parse_state,
                               struct macro_rate_ruleset *rates)
{
  if (macro_check_rates(parse_state, rates))
  {
    macro_free_rates(rates);
    return 1;
  }
  return 0;
}

/*************************************************************************
 macro_new_ruleset:
    Allocate a new rate ruleset for a macromolecule.  A rate ruleset is a named
    set of rules for associating a particular state of a set of subunits to a
    reaction rate.  This data structure is transferred to a table and freed
    before the simulation starts.

 In:  name: the name for this ruleset
      rules: a linked list of rules for this ruleset
 Out: Freshly allocated structure, or NULL if we failed to allocate
*************************************************************************/
static struct macro_rate_ruleset *macro_new_ruleset(char *name,
                                                    struct macro_rate_rule *rules)
{

  struct macro_rate_ruleset *mrr;
  mrr = CHECKED_MALLOC_STRUCT(struct macro_rate_ruleset, "macromolecular rate ruleset");
  if (mrr == NULL)
    return NULL;
  mrr->next = NULL;
  mrr->name = name;
  mrr->rules = rules;
  return mrr;
}

/*************************************************************************
 mdl_assemble_complex_ruleset:
    Assemble a macromolecular complex rate rule set.

 In:  parse_state: parser state
      name: name of the ruleset
      rules: rules for this ruleset
 Out: ruleset, or NULL if an error occurred
*************************************************************************/
struct macro_rate_ruleset *mdl_assemble_complex_ruleset(struct mdlparse_vars *parse_state,
                                                        char *name,
                                                        struct macro_rate_rule *rules)
{
  struct macro_rate_ruleset *ruleset;
  if ((ruleset = macro_new_ruleset(name, rules)) == NULL)
  {
    free(name);
    macro_free_rate_rules(rules);
    return NULL;
  }

  return ruleset;
}

/*************************************************************************
 macro_new_rule:
    Allocate a new rate rule for a macromolecule.  A rate rule associates a
    particular state of related subunits to a particular reaction rate.  This
    is used only at parse time and is transferred to a table before the
    simulation starts.

 In:  clauses: linked list of conditions which must be met
      rate: the reaction rate
 Out: Freshly allocated structure, or NULL if we failed to allocate
*************************************************************************/
static struct macro_rate_rule *macro_new_rule(struct macro_rate_clause *clauses,
                                              double rate)
{

  struct macro_rate_rule *mrr;
  mrr = CHECKED_MALLOC_STRUCT(struct macro_rate_rule, "macromolecular rate rule");
  if (mrr == NULL)
    return NULL;
  mrr->next = NULL;
  mrr->clauses = clauses;
  mrr->rate = rate;
  return mrr;
}

/*************************************************************************
 mdl_assemble_complex_rate_rule:
    Assemble a macromolecular complex rate rule.

 In:  parse_state: parser state
      clauses: the clauses for the rate rule
      rate: resultant rate if this clause is matched
 Out: rate rule, or NULL if an error occurred
*************************************************************************/
struct macro_rate_rule *mdl_assemble_complex_rate_rule(struct mdlparse_vars *parse_state,
                                                       struct macro_rate_clause *clauses,
                                                       double rate)
{
  struct macro_rate_rule *rule;
  if ((rule = macro_new_rule(clauses, rate)) == NULL)
  {
    macro_free_rate_clauses(clauses);
    return NULL;
  }
  return rule;
}

/*************************************************************************
 macro_check_rate_clause:
    Check the rate rule clause for validity.  This checks the types of species
    referenced in the clause, and also ensures that the relationship specified
    actually exists.

 In:  parse_state: the parser state for error logging
      topo: the macromol topology
      rel_name: the relationship name
      invert: the invert flag for the clause
      mol_type: the referenced species
 Out: 0 if the rate clause is valid, 1 if not.  On error, the details are
      written to the error log.
*************************************************************************/
static int macro_check_rate_clause(struct mdlparse_vars *parse_state,
                                   struct macro_relationship *rels,
                                   char const *rel_name,
                                   int invert,
                                   struct species *mol_type)
{
  UNUSED(invert);

  /* Disallow other complex molecules */
  if (mol_type->flags & IS_COMPLEX)
  {
    mdlerror_fmt(parse_state,
                 "In complex species '%s', ruleset references a complex molecule '%s'",
                 parse_state->complex_name,
                 mol_type->sym->name);
    return 1;
  }

  /* Check that the relationship exists */
  struct macro_relationship *the_rel;
  for (the_rel = rels; the_rel != NULL; the_rel = the_rel->next)
  {
    if (! strcmp(the_rel->name, rel_name))
      break;
  }

  /* If the_rel is null, we didn't find a matching relationship */
  if (the_rel == NULL)
  {
    mdlerror_fmt(parse_state,
                 "In complex species '%s', ruleset references a non-existent relationship '%s'",
                 parse_state->complex_name,
                 rel_name);
    return 1;
  }

  return 0;
}

/*************************************************************************
 macro_new_rate_clause:
    Allocate a new rate rule clause for a macromolecule.  A rate rule clause
    currently represents a condition of the form:

        RELATED_SUBUNIT == species [OPTIONAL_ORIENT]
      or:
        RELATED_SUBUNIT != species [OPTIONAL_ORIENT]

    This data structure is used only at parse time and is transferred to a
    table and freed before the simulation starts.

 In:  name: name of relationship for subunit
      invert: if 0, ==, else !=
      mol: species for RHS of rule
      orient: orient for RHS of rule
 Out: Freshly allocated structure, or NULL if we failed to allocate
*************************************************************************/
static struct macro_rate_clause *macro_new_rate_clause(char *name,
                                                       int invert,
                                                       struct species *mol,
                                                       short orient)
{

  struct macro_rate_clause *mrc;
  mrc = CHECKED_MALLOC_STRUCT(struct macro_rate_clause, "macromolecular rate rule clause");;
  if (mrc == NULL)
    return NULL;
  mrc->next = NULL;
  mrc->name = name;
  mrc->invert = invert;
  mrc->species = mol;
  mrc->orient = orient;
  return mrc;
}

/*************************************************************************
 mdl_assemble_complex_rate_rule_clause:
    Assemble a macromolecular rate rule clause.

 In:  parse_state: parser state
      rels: defined relations for this macromolecule
      relation_name: name of the relation for this clause
      invert: is this an == or a != clause?
      target: target type of molecule for this clause
 Out: newly allocated rate clause, or NULL if an error occurred
*************************************************************************/
struct macro_rate_clause *mdl_assemble_complex_rate_rule_clause(struct mdlparse_vars *parse_state,
                                                                struct macro_relationship *rels,
                                                                char *relation_name,
                                                                int invert,
                                                                struct species_opt_orient *target)
{
  struct macro_rate_clause *clause;
  struct species *subunit_type = (struct species *) target->mol_type->value;
  short orient = target->orient;
  if (macro_check_rate_clause(parse_state, rels, relation_name, invert, subunit_type))
  {
    free(relation_name);
    return NULL;
  }

  if (orient > 0) orient = 1;
  else if (orient < 0) orient = -1;

  if ((clause = macro_new_rate_clause(relation_name, invert, subunit_type, orient)) == NULL)
  {
    free(relation_name);
    return NULL;
  }

  return clause;
}

/*************************************************************************
 mdl_assemble_topology:
    Allocate a new topology data structure for a macromolecule.  Currently, a
    topology is always the cartesian product of 1 or more finite sets.  This
    data structure is used only at parse time, and then it is freed.

 In:  parse_state: parser state for error reporting
      dims: size of each dimension in the topology
 Out: Freshly allocated structure, or NULL if we failed to allocate
*************************************************************************/
struct macro_topology *mdl_assemble_topology(struct mdlparse_vars *parse_state,
                                             struct num_expr_list_head *dims)
{
  int cur_dim;
  struct num_expr_list *nel;
  struct macro_topology *mtopo;

  /* Allocate and initialize the topology structure */
  if ((mtopo = CHECKED_MALLOC_STRUCT(struct macro_topology, "macromolecule topology")) == NULL)
  {
    if (! dims->shared)
      mdl_free_numeric_list(dims->value_head);
    return NULL;
  }
  mtopo->total_subunits = 0;

  /* Allocate a single topology element (i.e. a single MxNx...xP array.  We may
   * later support other topologies, but this is the only one for now.
   */
  if ((mtopo->head = CHECKED_MALLOC_STRUCT(struct macro_topology_element, "macromolecule topology element")) == NULL)
  {
    free(mtopo);
    if (! dims->shared)
      mdl_free_numeric_list(dims->value_head);
    return NULL;
  }
  mtopo->head->next = NULL;
  mtopo->head->name = NULL;

  /* Count dimensionality, allocate space, and copy dimensions into array */
  mtopo->head->dimensionality = dims->value_count;
  if ((mtopo->head->dimensions =  CHECKED_MALLOC_ARRAY(int, mtopo->head->dimensionality,
                                                       "macromolecule topology")) == NULL)
  {
    free(mtopo->head);
    free(mtopo);
    if (! dims->shared)
      mdl_free_numeric_list(dims->value_head);
    return NULL;
  }
  for (cur_dim = 0, nel = dims->value_head; nel != NULL; nel = nel->next)
  {
    int dim = (int) nel->value;
    if (dim <= 0)
    {
      free(mtopo->head->dimensions);
      free(mtopo->head);
      free(mtopo);
      mdlerror_fmt(parse_state, "Dimensions for subunit topology must be greater than 0");
      if (! dims->shared)
        mdl_free_numeric_list(dims->value_head);
      return NULL;
    }
    mtopo->head->dimensions[cur_dim ++] = dim;
  }

  /* Compute the total number of subunits in this topology */
  struct macro_topology_element *mtopo_el;
  for (mtopo_el = mtopo->head; mtopo_el != NULL; mtopo_el = mtopo_el->next)
  {
    int nitems = 1;
    for (cur_dim = 0; cur_dim < mtopo_el->dimensionality; ++ cur_dim)
      nitems *= mtopo_el->dimensions[cur_dim];
    mtopo->total_subunits += nitems;
  }
  if (! dims->shared)
    mdl_free_numeric_list(dims->value_head);
  return mtopo;
}

/*************************************************************************
 mdl_assemble_subunit_spec_component:
    Allocate and populate a new component in a subunit coordinate
    specification.  A linked list of these, with one item corresponding each
    coordinate dimension, will specify a specific subset of all subunits.

 In: from: lowest value to include for this coordinate
     to: highest value to include for this coordinate
 Out: Freshly allocated structure, or NULL if no structure could be allocated.
*************************************************************************/
struct macro_subunit_spec *mdl_assemble_subunit_spec_component(int from, int to)
{

  struct macro_subunit_spec *cmp;
  cmp = CHECKED_MALLOC_STRUCT(struct macro_subunit_spec, "macromolecular subunit specification");
  if (cmp == NULL)
    return NULL;

  cmp->next = NULL;
  cmp->from = from;
  cmp->to = to;
  return cmp;
}

/*************************************************************************
 macro_assign_initial_subunits_helper:
    Helper function to assign the initial subunits.  This recursive helper
    function is used to iterate over the (arbitrarily many) dimensions of the
    molecule topology.

    In: struct complex_species *cs - species to receive species assignments
        struct macro_topology *topo - the topology of this species
        struct macro_subunit_spec *cur_coord - the specification of the
                                   coordinates for the dimension currently under
                                   consideration
        struct species *subunit_species - the species to assign to the selected
                                   coordinates.
        short subunit_orient - the orientation to assign to the selected
                                   coordinates
        int cur_dim_idx - the 0-based index of the coordinate under
                                   consideration
        int cur_stride - how many linear elements do we skip for every
                                   increment or decrement of a coordinate in
                                   the current dimension
        int cur_accum - what is the base linear index for the coordinates in
                                   all of the "earlier" dimensions.
    Out: Nothing.  Species data structure is filled in with subunit types.

    N.B. subunit coordinates are in reverse order in 'cur_coord'
*************************************************************************/
static void macro_assign_initial_subunits_helper(struct complex_species *cs,
                                                 struct macro_topology *topo,
                                                 struct macro_subunit_spec *cur_coord,
                                                 struct species *subunit_species,
                                                 short subunit_orient,
                                                 int cur_dim_idx,
                                                 int cur_stride,
                                                 int cur_accum)
{
  int this_coord;
  cur_accum += (cur_coord->from-1) * cur_stride;
  if (cur_coord->next == NULL)
  {
    for (this_coord = cur_coord->from-1; this_coord < cur_coord->to; ++ this_coord)
    {
      cs->subunits[cur_accum] = subunit_species;
      if (cs->orientations)
        cs->orientations[cur_accum] = subunit_orient;
      cur_accum += cur_stride;
    }
  }
  else
  {
    int next_stride = cur_stride * topo->head->dimensions[cur_dim_idx];
    for (this_coord = cur_coord->from-1; this_coord < cur_coord->to; ++ this_coord)
    {
      macro_assign_initial_subunits_helper(cs, topo, cur_coord->next, subunit_species, subunit_orient, cur_dim_idx - 1, next_stride, cur_accum);
      cur_accum += cur_stride;
    }
  }
}

/*************************************************************************
 macro_assign_initial_subunits:
    Store the initial subunit assignments in the species.

 In: cs: species to receive species assignments
     topo: the topology of this species
     assignments: list of assignments of subunit species to subsets of the
                  topology
 Out: Nothing.  Species data structure is filled in with subunit types.

 N.B. subunit coordinates are in reverse order in 'assignments'
*************************************************************************/
static void macro_assign_initial_subunits(struct complex_species *cs,
                                          struct macro_topology *topo,
                                          struct macro_subunit_assignment *assignments)
{
  while (assignments)
  {
    macro_assign_initial_subunits_helper(cs,
                                         topo,
                                         assignments->head,
                                         assignments->what,
                                         assignments->orient,
                                         topo->head->dimensionality - 1,
                                         1,
                                         0);
    assignments = assignments->next;
  }
}

/*************************************************************************
 macro_assign_geometry:
    Convert parse-time geometry structures to run-time geometry structures.

 In: cs: species to receive geometry
     topo: the topology of this species
     geom: list of geometry to put in the table
 Out: Nothing.  Subunit locations are filled in cs data structure.
*************************************************************************/
static void macro_assign_geometry(struct complex_species *cs,
                                  struct macro_topology *topo,
                                  struct macro_geometry *geom)
{
  for (; geom != NULL; geom = geom->next)
  {
    int dim_index, linear_index;
    struct num_expr_list *nel = geom->index;
    for (dim_index = 0, linear_index = 0; dim_index < topo->head->dimensionality; ++ dim_index)
    {
      if (dim_index != 0)
        linear_index *= topo->head->dimensions[dim_index];
      linear_index += (int) nel->value - 1;
      nel = nel->next;
    }

    cs->rel_locations[linear_index] = geom->location;
  }
}

/*************************************************************************
 macro_free_runtime_relationships:
    Free the relationships table from the given complex species.

 In:  cs: the species
 Out: table is freed and its pointer cleared
*************************************************************************/
static void macro_free_runtime_relationships(struct complex_species *cs)
{
  if (cs->relations)
  {
    int i;
    for (i=0; i<cs->num_relations; ++ i)
    {
      if (cs->relations[i].target)
        free((void *) cs->relations[i].target);
      if (cs->relations[i].inverse)
        free((void *) cs->relations[i].inverse);
    }
    free((void *) cs->relations);
  }
  cs->relations = NULL;
  cs->num_relations = 0;
}

/*************************************************************************
 macro_assign_relationships:
    Convert parse-time relationship structures to run-time relationship
    structures.

 In: cs:   species to receive relationship table
     topo: the topology of this species
     rels: list of relationships to put in the
                                       table
 Out: 0 on success, 1 if allocation fails.  Relationship table is filled in cs
      data structure.
*************************************************************************/
static int macro_assign_relationships(struct complex_species *cs,
                                      struct macro_topology *topo,
                                      struct macro_relationship *rels)
{

  cs->relations = NULL;

  /* Count relations */
  struct macro_relationship *rels_temp;
  int count = 0;
  for (rels_temp = rels; rels_temp != NULL; rels_temp = rels_temp->next)
    ++ count;

  /* Allocate space for relations */
  struct subunit_relation *final_relations;
  final_relations = CHECKED_MALLOC_ARRAY(struct subunit_relation,
                                          count,
                                          "macromolecule relation table");
  if (final_relations == NULL)
    return 1;
  memset(final_relations, 0, count*sizeof(struct subunit_relation));
  cs->num_relations = count;
  cs->relations = final_relations;

  /* Compute scaling factors for each array index */
  int scales[ topo->head->dimensionality + 1 ];
  scales[ topo->head->dimensionality ] = 1;
  for (int i = topo->head->dimensionality; i > 0; -- i)
    scales[i-1] = scales[i] * topo->head->dimensions[i - 1];

  /* Copy out all relationships */
  for (count = 0; rels != NULL; ++ count, rels = rels->next)
  {
    int coords[ topo->head->dimensionality ];
    int *targets = CHECKED_MALLOC_ARRAY(int, topo->total_subunits, "macromolecule subunit relations");
    if (targets == NULL)
      return 1;
    for (int i = 0; i < topo->head->dimensionality; ++i)
      coords[i] = 0;

    /* Fill in the current entry in the relation table. */
    final_relations[count].name = rels->name;
    final_relations[count].target = targets;
    final_relations[count].inverse = NULL;
    rels->name = NULL;

    /* Step through the coordinate space in order... */
    struct num_expr_list *nel;
    int linear_idx = 0, linear_target = 0;
    while (coords[0] < topo->head->dimensions[0])
    {
      linear_target = linear_idx;
      int target_coord_index;
      /* Step through each coordinate */
      for (target_coord_index = 0, nel = rels->indices; nel != NULL; ++ target_coord_index, nel = nel->next)
      {
        int offset = (int) nel->value;

        /* If the coordinate is offset in the negative direction... */
        if (offset < 0)
        {
          /* If the offset causes wrap... */
          if (coords[target_coord_index] + offset < 0)
            linear_target += scales[target_coord_index] + offset*scales[target_coord_index + 1];
          else
            linear_target += offset*scales[target_coord_index + 1];
        }

        /* If the coordinate is offset in the positive direction... */
        else if (offset > 0)
        {
          /* If the offset causes wrap... */
          if (coords[target_coord_index] + offset >= topo->head->dimensions[target_coord_index])
            linear_target -= scales[target_coord_index] - offset*scales[target_coord_index + 1];
          else
            linear_target += offset*scales[target_coord_index + 1];
        }
      }

      /* Store the computed target */
      targets[linear_idx] = linear_target;

      /* Advance the coordinates */
      ++ linear_idx;
      for (int i = topo->head->dimensionality - 1; i >= 0; --i)
      {
        if (++ coords[i] == topo->head->dimensions[i]  &&  i > 0)
          coords[i] = 0;
        else
          break;
      }
    }
  }

  return 0;
}

/*************************************************************************
 macro_assign_inverse_relationships:
    Assign the inverse relationship tables for a species.  This is mainly used
    for special macromolecule counting at present.

 In: cs:   species to receive inverse relation table
     topo: topology for complex
 Out: 0 on success, 1 if allocation fails.  cs data structure has inverse
      relation table filled in
*************************************************************************/
static int macro_assign_inverse_relationships(struct complex_species *cs,
                                              struct macro_topology *topo)
{

  struct subunit_relation *final_relations = (struct subunit_relation *) cs->relations;
  for (int rel_index = 0; rel_index < cs->num_relations; ++ rel_index)
  {
    int su_index;

    /* Allocate inverse table and initialize entries to -1 */
    int *inverse = CHECKED_MALLOC_ARRAY(int, topo->total_subunits, "macromolecule subunit inverse relations");
    if (inverse == NULL)
      return 1;

    final_relations[rel_index].inverse = inverse;
    for (su_index = 0; su_index < topo->total_subunits; ++ su_index)
      inverse[su_index] = -1;

    /* Now, invert such that if target[i] == j, inverse[j] == i */
    /* Obviously, this relies on a one-to-one mapping to be useful */
    for (su_index = 0; su_index < topo->total_subunits; ++ su_index)
    {
      /* -1 indicates no relation.  this can't currently happen */
      int target = final_relations[rel_index].target[su_index];
      if (target != -1)
        inverse[target] = su_index;
    }
  }
  return 0;
}

/*************************************************************************
 macro_check_empty_subunits:
    Check that all subunits in the complex have initial species assignments.
    If not, an error message will be logged, and 1 returned.

 In:  parse_state: the parser state for error logging
      topo: the macromol topology
      cs: the complex to check
 Out: 0 if the complex has fully-specified initial state, 1 if not.  On error,
      the details are written to the error log.
*************************************************************************/
static int macro_check_empty_subunits(struct mdlparse_vars *parse_state,
                                      struct macro_topology *topo,
                                      struct complex_species *cs)
{
  int subunit_index;
  for (subunit_index=0; subunit_index<cs->num_subunits; ++subunit_index)
  {
    if (cs->subunits[subunit_index] == NULL)
    {
      char *ai = macro_linear_array_index_to_string(subunit_index, topo);
      mdlerror_fmt(parse_state,
                   "In complex species '%s', no subunit type is specified for subunit %s",
                   parse_state->complex_name,
                   ai);
      free(ai);
      return 1;
    }
  }

  return 0;
}

/*************************************************************************
 macro_free_assignments:
    Free a linked list of macromolecule initial subunit assignments.  This will
    free the memory associated with the structures.

 In:  assignments: the list to free
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_assignments(struct macro_subunit_assignment *assignments)
{
  while (assignments != NULL)
  {
    struct macro_subunit_assignment *a = assignments;
    macro_free_subunit_specs(assignments->head);
    assignments = assignments->next;
    free(a);
  }
}

/*************************************************************************
 macro_free_topology:
    Free a linked list of macromolecule topology.  This will free the memory
    associated with the structures.

 In:  topo: the list to free
 Out: Nothing.  The memory is freed.
*************************************************************************/
static void macro_free_topology(struct macro_topology *topo)
{
  while (topo->head != NULL)
  {
    struct macro_topology_element *el = topo->head;
    topo->head = topo->head->next;
    free(el->dimensions);
    free(el);
  }
  free(topo);
}

/*************************************************************************
 mdl_assemble_complex_species:
    Assemble a complex species, adding it to the symbol table.

 In:  parse_state: parser state
      name: species name
      topo: species topology
      assignments: species initial subunit assignments
      geometry: species geometry
      rels: species relationships
      rates: species rate rules
 Out: 0 on success, 1 on failure; species is added to symbol table
*************************************************************************/
int mdl_assemble_complex_species(struct mdlparse_vars *parse_state,
                                 char *name,
                                 struct macro_topology *topo,
                                 struct macro_subunit_assignment *assignments,
                                 struct macro_geometry *geom,
                                 struct macro_relationship *rels,
                                 struct macro_rate_ruleset *rates)
{
  struct complex_species *cs = new_complex_species(topo->total_subunits, parse_state->complex_type);

  /* Copy parser data into our species */
  macro_assign_initial_subunits(cs, topo, assignments);
  macro_assign_geometry(cs, topo, geom);
  if (macro_assign_relationships(cs, topo, rels))
    goto failure;
  if (macro_assign_inverse_relationships(cs, topo))
    goto failure;
  if (macro_build_rate_tables(parse_state, cs, rates))
    goto failure;

  /* Check that no subunits were left empty */
  if (macro_check_empty_subunits(parse_state, topo, cs))
    goto failure;

  /* Free parser data */
  macro_free_rates(rates);
  macro_free_relationships(rels);
  macro_free_geometry(geom);
  macro_free_assignments(assignments);
  macro_free_topology(topo);
  parse_state->complex_topo = NULL;
  parse_state->complex_relations = NULL;

  cs->base.sym = store_sym(name, MOL, parse_state->vol->mol_sym_table, cs);
  if (cs->base.sym == NULL)
  {
    mdlerror_fmt(parse_state, "Failed to store the complex species '%s' in the symbol table", name);
    free(name);
    goto failure;
  }

  free(name);
  parse_state->complex_name = NULL;
  return 0;

failure:
  macro_free_runtime_rate_tables(cs);
  macro_free_runtime_relationships(cs);
  macro_free_rates(rates);
  macro_free_relationships(rels);
  macro_free_geometry(geom);
  macro_free_assignments(assignments);
  macro_free_topology(topo);
  parse_state->complex_topo = NULL;
  parse_state->complex_relations = NULL;
  return 1;
}

/*************************************************************************
 * prepare_reactions and related machinery
 *************************************************************************/

/*************************************************************************
 equivalent_geometry_for_two_reactants:

 In: o1a: orientation of the first reactant from first reaction
     o1b: orientation of the second reactant from first reaction
     o2a: orientation of the first reactant from second reaction
     o2b: orientation of the second reactant from second reaction
 Out: Returns 1 if the two pathways (defined by pairs o1a-o1b and o2a-o2b)
      have equivalent geometry, 0 otherwise.
*************************************************************************/
static int equivalent_geometry_for_two_reactants(int o1a, int o1b, int o2a, int o2b)
{

    /* both reactants for each pathway are in the same
       orientation class and parallel one another */
    if ((o1a == o1b) && (o2a == o2b))
    {
       return 1;
    /* both reactants for each pathway are in the same
       orientation class and opposite one another */
    }
    else if ((o1a == -o1b) && (o2a == -o2b))
    {
       return 1;
    }
    /* reactants are not in the same orientation class */
    if (abs(o1a) != abs(o1b))
    {
       if ((abs(o2a) != abs(o2b)) || ((o2a == 0) && (o2b == 0)))
       {
          return 1;
       }
    }
    if (abs(o2a) != abs(o2b))
    {
       if ((abs(o1a) != abs(o1b)) || ((o1a == 0) && (o1b == 0)))
       {
          return 1;
       }
    }

    return 0;
}

/*************************************************************************
 equivalent_geometry:

 In: p1, p2: pathways to compare
     n: The number of reactants for the pathways
 Out: Returns 1 if the two pathways are the same (i.e. have equivalent
      geometry), 0 otherwise.
*************************************************************************/
static int equivalent_geometry(struct pathway *p1, struct pathway *p2, int n)
{

  short o11,o12,o13,o21,o22,o23; /* orientations of individual reactants */
  /* flags for 3-reactant reactions signaling whether molecules orientations
   * are parallel one another and molecule and surface orientaions are parallel
   * one another
   */
  int mols_parallel_1 = SHRT_MIN + 1; /* for first pathway */
  int mols_parallel_2 = SHRT_MIN + 2; /* for second pathway */
  int mol_surf_parallel_1 = SHRT_MIN + 3; /* for first pathway */
  int mol_surf_parallel_2 = SHRT_MIN + 4; /* for second pathway */

  if (memcmp(p1->is_complex, p2->is_complex, 3))
    return 0;

  if (n < 2)
  {
     /* one reactant case */
     /* RULE: all one_reactant pathway geometries are equivalent */

      return 1;

  }
  else if (n < 3)
  {
    /* two reactants case */

    /* RULE - Two pathways have equivalent geometry when:
       1) Both pathways have exactly the same number of reactants;
       2) There exists an identity mapping between reactants from Pathway 1 and
          Pathway 2 such that for each pair of reactants, r1a and r1b from Pathway
          1, and r2a, and r2b from Pathway 2:
         - r1a is the same species as r2a (likewise for r1b and r2b);
         - r1a and r1b have the same orientation in the same orientation class
           if and only if r2a and r2b do;
         - r1a and r1b have the opposite orientation in the same orientation
           class if and only if r2a and r2b do;
         - r1a and r1b are not in the same orientation class, either because
           they have different orientation classes or both are in the zero
           orientation class, if and only if r2a and r2b are not in the same
           orientation class or both are in the zero orientation class
     */

    o11 = p1->orientation1;
    o12 = p1->orientation2;
    o21 = p2->orientation1;
    o22 = p2->orientation2;

    o13 = o23 = 0;


    return equivalent_geometry_for_two_reactants(o11, o12, o21, o22);

  }
  else if (n < 4)
  {
     /* three reactants case */

    o11 = p1->orientation1;
    o12 = p1->orientation2;
    o13 = p1->orientation3;
    o21 = p2->orientation1;
    o22 = p2->orientation2;
    o23 = p2->orientation3;

    /* special case: two identical reactants */
      if ((p1->reactant1 == p1->reactant2)
          && (p2->reactant1 == p2->reactant2))
      {

       /* Case 1: two molecules and surface are in the same orientation class */
        if ((abs(o11) == abs(o12)) && (abs(o11) == abs(o13)))
        {
          if (o11 == o12) mols_parallel_1 = 1;
          else mols_parallel_1 = 0;

          if (mols_parallel_1)
          {
            if ((o11 == -o13) || (o12 == -o13))
            {
               mol_surf_parallel_1 = 0;
            }
            else 
            {
               mol_surf_parallel_1 = 1;
            }
          }
          else 
          {
               mol_surf_parallel_1 = 0;
          }

          if ((abs(o21) == abs(o22)) && (abs(o21) == abs(o23)))
          {
             if (o21 == o22) mols_parallel_2 = 1;
             else mols_parallel_2 = 0;

             if (mols_parallel_2)
             {
               if ((o21 == -o23) || (o22 == -o23))
               {
                  mol_surf_parallel_2 = 0;
               }
               else 
               {
                  mol_surf_parallel_2 = 1;
               }
             }
             else 
             {
                  mol_surf_parallel_2 = 0;
             }

          }

          if ((mols_parallel_1 == mols_parallel_2) &&
              (mol_surf_parallel_1 == mol_surf_parallel_2))
              {
                 return 1;
          }

       } /* end case 1 */

       /* Case 2: one molecule and surface are in the same orientation class */
       else if ((abs(o11) == abs(o13)) || (abs(o12) == abs(o13)))
       {
          if ((o11 == o13) || (o12 == o13)) mol_surf_parallel_1 = 1;
          else mol_surf_parallel_1 = 0;

          /* check that pathway2 is also in the case2 */

          if ((abs(o21) != abs(o23)) || (abs(o22) != abs(o23)))
          {
             if ((abs(o21) == abs(o23)) || (abs(o22) == abs(o23)))
             {
                if ((o21 == o23) || (o22 == o23)) mol_surf_parallel_2 = 1;
                else mol_surf_parallel_2 = 0;

             }
          }
          if (mol_surf_parallel_1 == mol_surf_parallel_2)
          {
             return 1;
          }

       } /* end case 2 */

       /* Case 3: two molecules but not surface are in the same
                  orientation class */
       else if ((abs(o11) == abs(o12)) && (abs(o11) != abs(o13)))
       {
          if (o11 == o12) mols_parallel_1 = 1;
          else mols_parallel_1 = 0;

          if ((abs(o21) == abs(o22)) && (abs(o21) != abs(o23)))
          {
             if (o21 == o22) mols_parallel_2 = 1;
             else mols_parallel_2 = 0;
          }
          if (mols_parallel_1 == mols_parallel_2)
          {
                 return 1;
          }

       }
       /* Case 4: all molecules and surface are in different orientation classes */
       else if ((abs(o11) != abs(o13)) && (abs(o12) != abs(o13)) &&
                 (abs(o11) != abs(o12)))
                 {

          if ((abs(o21) != abs(o23)) && (abs(o22) != abs(o23)) &&
                 (abs(o21) != abs(o22)))
                 {
               return 1;
          }
       } /* end all cases */

    }
    else { /* no identical reactants */

       if ((equivalent_geometry_for_two_reactants(o11, o12, o21, o22))
           && (equivalent_geometry_for_two_reactants(o12, o13, o22, o23))
           && (equivalent_geometry_for_two_reactants(o11, o13, o21, o23)))
           {
                return 1;
       }

    }

  } // end if (n < 4)


  return 0;
}

/*************************************************************************
 create_sibling_reaction:
    Create a sibling reaction to the given reaction -- a reaction into which
    some of the pathways may be split by split_reaction.

 In:  rx:   reaction for whom to create sibling
 Out: sibling reaction, or NULL on error
*************************************************************************/
static struct rxn *create_sibling_reaction(struct rxn *rx)
{

  struct rxn *reaction = CHECKED_MALLOC_STRUCT(struct rxn, "reaction");
  if (reaction == NULL)
    return NULL;
  reaction->next = NULL;
  reaction->sym = rx->sym;
  reaction->n_reactants = rx->n_reactants;
  reaction->n_pathways = 0;
  reaction->cum_probs = NULL;
  reaction->product_idx = NULL;
  reaction->rates = NULL;
  reaction->max_fixed_p = 0.0;
  reaction->min_noreaction_p = 0.0;
  reaction->pb_factor = 0.0;
  reaction->players = NULL;
  reaction->geometries = NULL;
  reaction->is_complex = NULL;
  reaction->n_occurred = 0;
  reaction->n_skipped = 0.0;
  reaction->prob_t = NULL;
  reaction->pathway_head = NULL;
  reaction->info = NULL;
  return reaction;
}

/*************************************************************************
 split_reaction:
 In:  parse_state: parser state
      rx: reaction to split
 Out: Returns head of the linked list of reactions where each reaction
      contains only geometrically equivalent pathways
*************************************************************************/
static struct rxn *split_reaction(struct mdlparse_vars *parse_state, struct rxn *rx)
{
  struct rxn  *curr_rxn_ptr = NULL,  *head = NULL, *end = NULL;
  struct rxn *reaction;
  struct pathway *to_place, *temp;

  /* keep reference to the head of the future linked_list */
  head = end = rx;
  to_place = head->pathway_head->next;
  head->pathway_head->next = NULL;
  head->n_pathways = 1;
  while (to_place != NULL)
  {
    if (to_place->flags & (PATHW_TRANSP | PATHW_REFLEC | PATHW_ABSORP | PATHW_CLAMP_CONC))
    {
      reaction = create_sibling_reaction(rx);
      if (reaction == NULL)
        return NULL;

      reaction->pathway_head = to_place;
      to_place = to_place->next;
      reaction->pathway_head->next = NULL;
      ++ reaction->n_pathways;

      end->next = reaction;
      end = reaction;
    }
    else
    {
      for (curr_rxn_ptr = head; curr_rxn_ptr != NULL; curr_rxn_ptr = curr_rxn_ptr->next)
      {
        if (curr_rxn_ptr->pathway_head->flags & (PATHW_TRANSP | PATHW_REFLEC | PATHW_ABSORP))
          continue;
        if (equivalent_geometry(to_place, curr_rxn_ptr->pathway_head, curr_rxn_ptr->n_reactants))
          break;
      }

      if (! curr_rxn_ptr)
      {
        reaction = create_sibling_reaction(rx);
        if (reaction == NULL)
          return NULL;

        end->next = reaction;
        end = reaction;

        curr_rxn_ptr = end;
      }

      temp = to_place;
      to_place = to_place->next;

      temp->next = curr_rxn_ptr->pathway_head;
      curr_rxn_ptr->pathway_head = temp;
      ++ curr_rxn_ptr->n_pathways;
    }
  }

  return head;
}

/*************************************************************************
 check_reaction_for_duplicate_pathways:
 In:  head: head of linked list of pathways
 Out: Sorts linked list of pathways in alphabetical order according to the
      "prod_signature" field.  Checks for the duplicate pathways.  Prints error
      message and exits simulation if duplicates found.
 Note: This function is called after 'split_reaction()' function so all
       pathways have equivalent geometry from the reactant side.  Here we check
       whether relative orientation of all players (both reactants and
       products) is the same for the two seemingly identical pathways.
 RULE: Two reactions pathways are duplicates if and only if
        (a) they both have the same number and species of reactants;
        (b) they both have the same number and species of products;
        (c) there exists a bijective mapping between the reactants and products
            of the two pathways such that reactants map to reactants, products
            map to products, and the two pathways have equivalent geometry
            under mapping.
            Two pathways R1 and R2 have an equivalent geometry under a mapping
            M if and only if for every pair of players "i" and "j" in R1, the
            corresponding players M(i) and M(j) in R2 have the same orientation
            relation as do "i" and "j" in R1.
            Two players "i" and "j" in a reaction pathway have the following
            orientation:
              parallel - if both "i" and "j" are in the same nonzero orientation
              class with the same sign;
              antiparallel (opposite) - if they are both in the same nonzero
              orientation class but have opposite sign;
              independent - if they are in different orientation classes or both
              in the zero orientation class.

 PostNote: In this function we check only the validity of Rule (c) since
           conditions of Rule (a) and (b) are already satisfied when the
           function is called.
*************************************************************************/
static void check_reaction_for_duplicate_pathways(struct pathway **head)
{

   struct pathway *result = NULL; /* build the sorted list here */
   struct pathway *null_result = NULL; /* put pathways with NULL
                                          prod_signature field here */
   struct pathway *current, *next, **pprev;
   struct product *iter1, *iter2;
   int pathways_equivalent;  /* flag */
   int i, j;
   int num_reactants; /* number of reactants in the pathway */
   int num_products; /* number of products in the pathway */
   int num_players; /* total number of reactants and products in the pathway */
   int *orient_players_1, *orient_players_2; /* array of orientations of players */
   int o1a, o1b, o2a,o2b;

   /* extract  pathways with "prod_signature" field equal to NULL
     into "null_result" list */
   current = *head;
   pprev = head;
   while (current != NULL)
   {
     if (current->prod_signature == NULL)
     {
       *pprev = current->next;
       current->next = null_result;
       null_result = current;
       current = *pprev;
     }
     else
     {
       pprev = &current->next;
       current = current->next;
     }
   }

   /* check for duplicate pathways in null_result */
     current = null_result;
     if ((current != NULL) && (current->next != NULL))
     {
       /* From the previously called function "split_reaction()"
          we know that reactant-reactant pairs in two pathways
          are equivalent. Because there are no products the pathways
          are duplicates.
          RULE: There may be no more than one pathway with zero (--->NULL)
                products in the reaction->pathway_head
                after calling the function "split_reaction()"
       */
       if (current->reactant2 == NULL)
         mcell_error("Exact duplicates of reaction %s  ----> NULL are not allowed.  Please verify that orientations of reactants are not equivalent.",
                     current->reactant1->sym->name);
       else if (current->reactant3 == NULL)
         mcell_error("Exact duplicates of reaction %s + %s  ----> NULL are not allowed.  Please verify that orientations of reactants are not equivalent.",
                     current->reactant1->sym->name,
                     current->reactant2->sym->name);
       else
         mcell_error("Exact duplicates of reaction %s + %s + %s  ----> NULL are not allowed.  Please verify that orientations of reactants are not equivalent.",
                     current->reactant1->sym->name,
                     current->reactant2->sym->name,
                     current->reactant3->sym->name);
     }

    /* now sort the remaining pathway list by "prod_signature" field
       and check for the duplicates */
     current = *head;

  while(current != NULL)
  {
     next = current->next;

     /* insert in sorted order into the "result" */
     if (result == NULL || (strcmp(result->prod_signature, current->prod_signature) >= 0))
     {
        current->next = result;
        result = current;
     }
     else 
     {
        struct pathway *iter = result;
        while(iter->next != NULL && (strcmp(iter->next->prod_signature, current->prod_signature) < 0))
        {
             iter = iter->next;
        }
        current->next = iter->next;
        iter->next = current;
     }

     /* move along the original list */
     current = next;
  }

   /* Now check for the duplicate pathways */
   /* Since the list is sorted we can proceed down the list
      and compare the adjacent nodes */

   current = result;

   if (current != NULL)
   {
     while(current->next != NULL) 
     {
       if (strcmp(current->prod_signature, current->next->prod_signature) == 0)
       {

         pathways_equivalent  = 1;
         /* find total number of players in the pathways */
         num_reactants = 0;
         num_products = 0;
         if (current->reactant1 != NULL) num_reactants++;
         if (current->reactant2 != NULL) num_reactants++;
         if (current->reactant3 != NULL) num_reactants++;

         iter1 = current->product_head;
         while(iter1 != NULL)
         {
           num_products++;
           iter1 = iter1->next;
         }

         num_players = num_reactants + num_products;

         /* create arrays of players orientations */
         orient_players_1 = CHECKED_MALLOC_ARRAY(int, num_players, "reaction player orientations");
         if (orient_players_1 == NULL)
           mcell_die();
         orient_players_2 = CHECKED_MALLOC_ARRAY(int, num_players, "reaction player orientations");
         if (orient_players_2 == NULL)
           mcell_die();

         if (current->reactant1!=NULL) orient_players_1[0]=current->orientation1;
         if (current->reactant2!=NULL) orient_players_1[1]=current->orientation2;
         if (current->reactant3!=NULL) orient_players_1[2]=current->orientation3;
         if (current->next->reactant1!=NULL) orient_players_2[0]=current->next->orientation1;
         if (current->next->reactant2!=NULL) orient_players_2[1]=current->next->orientation2;
         if (current->next->reactant3!=NULL) orient_players_2[2]=current->next->orientation3;


         iter1 = current->product_head;
         iter2 = current->next->product_head;

         for (i = num_reactants; i < num_players; i++)
         {
           orient_players_1[i] = iter1->orientation;
           orient_players_2[i] = iter2->orientation;
           iter1 = iter1->next;
           iter2 = iter2->next;
         }


         /* below we will compare only reactant-product
            and product-product combinations
            because reactant-reactant combinations
            were compared previously in the function
            "equivalent_geometry()"
            */

         /* Initial assumption - pathways are equivalent.
            We check whether this assumption is
            valid by  comparing pairs as described
            above */

         i = 0;
         while((i < num_players) && (pathways_equivalent))
         {
           if (i < num_reactants)
           {
             j = num_reactants;
           }
           else 
           {
             j = i + 1;
           }
           for (; j < num_players; j++)
           {
             o1a = orient_players_1[i];
             o1b = orient_players_1[j];
             o2a = orient_players_2[i];
             o2b = orient_players_2[j];
             if (!equivalent_geometry_for_two_reactants(o1a, o1b, o2a, o2b))
             {
               pathways_equivalent = 0;
               break;
             }
           }
           i++;
         }

         if (pathways_equivalent)
         {
           if (current->reactant2 == NULL)
             mcell_error("Exact duplicates of reaction %s  ----> %s are not allowed.  Please verify that orientations of reactants are not equivalent.",
                         current->reactant1->sym->name,
                         current->prod_signature);
           else if (current->reactant3 == NULL)
             mcell_error("Exact duplicates of reaction %s + %s  ----> %s are not allowed.  Please verify that orientations of reactants are not equivalent.",
                         current->reactant1->sym->name,
                         current->reactant2->sym->name,
                         current->prod_signature);
           else
             mcell_error("Exact duplicates of reaction %s + %s + %s  ----> %s are not allowed.  Please verify that orientations of reactants are not equivalent.",
                         current->reactant1->sym->name,
                         current->reactant2->sym->name,
                         current->reactant3->sym->name,
                         current->prod_signature);
         }
       }

       current = current->next;
     }
   }

    if (null_result == NULL)
    {
       *head = result;
    }
    else if (result == NULL)
    {
       *head = null_result;
    }
    else
    {
       current = result;
       while(current->next != NULL)
       {
          current = current->next;
       }
       current->next = null_result;
       null_result->next = NULL;

       *head = result;
    }

}

/*************************************************************************
 reaction_has_complex_rates:
    Check if a reaction has any complex rates.


    In:  struct rxn *rx - the reaction to check
    Out: 1 if the reaction has complex pathways, 0 otherwise
*************************************************************************/
static int reaction_has_complex_rates(struct rxn *rx)
{
  struct pathway *path;
  for (path = rx->pathway_head;
       path != NULL;
       path = path->next)
  {
    if (path->km_complex)
      return 1;
  }

  return 0;
}

/*************************************************************************
 reorder_varying_pathways:
    Sort pathways so that all complex rates come at the end.  This allows us to
    quickly determine whether a reaction definitely occurs, definitely does not
    occur, or may occur depending on the states of the subunits in the complex.

    In:  struct rxn *rx - the reaction whose pathways to sort
    Out: 1 if the reaction has complex pathways, 0 otherwise

    XXX: Worthwhile sorting pathways by probability?
*************************************************************************/
static int reorder_varying_pathways(struct rxn *rx)
{

  int num_fixed = 0, num_varying = 0;
  int num_fixed_players = 0, num_varying_players = 0;
  int pathway_idx;
  int already_sorted = 1;

  /* If we have no rates, we're done */
  if (! rx->rates)
    return 0;

  /* Count fixed and varying pathways and players */
  for (pathway_idx = 0; pathway_idx < rx->n_pathways; ++ pathway_idx)
  {
    int player_count = rx->product_idx[pathway_idx+1] - rx->product_idx[pathway_idx];
    if (! rx->rates[pathway_idx])
    {
      ++ num_fixed;
      num_fixed_players += player_count;
      if (num_varying) already_sorted = 0;
    }
    else
    {
      ++ num_varying;
      num_varying_players += player_count;
    }
  }

  /* If all are fixed or all are varying, we're done */
  if (! num_fixed  || ! num_varying)
    return 0;

  /* If all fixed pathways already precede all varying pathways, we're done
   */
  if (already_sorted)
    return 0;

  /* Allocate space for sorted info */
  int pathway_mapping[rx->n_pathways];
  struct species **newplayers = NULL;
  short *newgeometries = NULL;
  unsigned char *new_is_complex = NULL;
  u_int *new_product_index = NULL;
  double *new_cum_probs = NULL;
  struct complex_rate **new_complex_rates = NULL;
  struct pathway_info *new_pathway_info = NULL;

  if ((newplayers = CHECKED_MALLOC_ARRAY(struct species*, rx->product_idx[rx->n_pathways], "reaction players array")) == NULL)
    goto failure;
  if ((newgeometries = CHECKED_MALLOC_ARRAY(short, rx->product_idx[rx->n_pathways], "reaction geometries array")) == NULL)
    goto failure;
  if (rx->is_complex)
    if ((new_is_complex = CHECKED_MALLOC_ARRAY(unsigned char, rx->product_idx[rx->n_pathways], "reaction 'is complex' flag array")) == NULL)
      goto failure;
  if ((new_product_index = CHECKED_MALLOC_ARRAY(u_int, rx->product_idx[rx->n_pathways] + 1, "reaction product index array")) == NULL)
    goto failure;
  if ((new_cum_probs = CHECKED_MALLOC_ARRAY(double, rx->n_pathways, "reaction cumulative probabilities array")) == NULL)
    goto failure;
  if ((new_complex_rates = CHECKED_MALLOC_ARRAY(struct complex_rate *, rx->n_pathways, "reaction complex rates array")) == NULL)
    goto failure;
  if ((new_pathway_info = CHECKED_MALLOC_ARRAY(struct pathway_info, rx->n_pathways, "reaction pathway info")) == NULL)
    goto failure;

  memcpy(newplayers, rx->players, rx->n_reactants * sizeof(struct species *));

  /* Now, step through the array until all fixed rates are at the beginning
   */
  int placed_fixed = 0, placed_varying = 0;
  int idx=0;
  int next_player_fixed = rx->n_reactants;
  int next_player_varying = rx->n_reactants + num_fixed_players;
  for (idx = 0; idx < rx->n_pathways; ++ idx)
  {
    int dest_player_idx;
    int dest_pathway;
    int num_players_to_copy = rx->product_idx[idx+1] - rx->product_idx[idx];

    /* Figure out where to put this pathway */
    if (! rx->rates[idx])
    {
      dest_player_idx = next_player_fixed;
      dest_pathway = placed_fixed;

      ++ placed_fixed;
      next_player_fixed += rx->product_idx[idx+1] - rx->product_idx[idx];
    }
    else
    {
      dest_player_idx = next_player_varying;
      dest_pathway = num_fixed + placed_varying;

      ++ placed_varying;
      next_player_varying += rx->product_idx[idx+1] - rx->product_idx[idx];
    }
    pathway_mapping[idx] = next_player_fixed;

    /* Copy everything in */
    memcpy(newplayers + dest_player_idx,
           rx->players + rx->product_idx[idx],
           sizeof(struct species *) * num_players_to_copy);
    memcpy(newgeometries + dest_player_idx,
           rx->geometries + rx->product_idx[idx],
           sizeof(short) * num_players_to_copy);
    if (new_is_complex)
      memcpy(new_is_complex + dest_player_idx,
             rx->is_complex + rx->product_idx[idx],
             sizeof(unsigned char) * num_players_to_copy);
    new_product_index[dest_pathway] = dest_player_idx;
    new_cum_probs[dest_pathway] = rx->cum_probs[idx];
    new_complex_rates[dest_pathway] = rx->rates[idx];
    new_pathway_info[dest_pathway].count = 0.0;
    new_pathway_info[dest_pathway].pathname = rx->info[idx].pathname;
    if (rx->info[idx].pathname) rx->info[idx].pathname->path_num = dest_pathway;
  }
  new_product_index[rx->n_pathways] = rx->product_idx[rx->n_pathways];

  /* Now, fix up varying rates */
  struct t_func *tf;
  for (tf = rx->prob_t; tf != NULL; tf = tf->next)
    tf->path = pathway_mapping[tf->path];

  /* Swap in newly ordered items */
  free(rx->players);
  free(rx->geometries);
  if (rx->is_complex)
    free(rx->is_complex);
  free(rx->product_idx);
  free(rx->cum_probs);
  free(rx->rates);
  free(rx->info);

  rx->players = newplayers;
  rx->geometries = newgeometries;
  rx->is_complex = new_is_complex;
  rx->product_idx = new_product_index;
  rx->cum_probs = new_cum_probs;
  rx->rates = new_complex_rates;
  rx->info = new_pathway_info;

  return 0;

failure:
  if (newplayers) free(newplayers);
  if (newgeometries) free(newgeometries);
  if (new_is_complex) free(new_is_complex);
  if (new_product_index) free(new_product_index);
  if (new_cum_probs) free(new_cum_probs);
  if (new_complex_rates) free(new_complex_rates);
  if (new_pathway_info) free(new_pathway_info);
  return 1;
}

/*************************************************************************
 set_reaction_player_flags:
    Set the reaction player flags for all participating species in this
    reaction.

 In:  rx: the reaction
 Out: Nothing.  Species flags may be updated.
*************************************************************************/
static void set_reaction_player_flags(struct rxn *rx)
{
  switch (rx->n_reactants)
  {
    case 1:
      /* do nothing */
      return;

    case 2:
      if (strcmp(rx->players[0]->sym->name, "ALL_MOLECULES")==0)
      {
          rx->players[0]->flags |= (CAN_MOLWALL|CAN_GRIDWALL);
      }
      else if (strcmp(rx->players[0]->sym->name, "ALL_VOLUME_MOLECULES")==0)
      {
          rx->players[0]->flags |= CAN_MOLWALL;
      }
      else if (strcmp(rx->players[0]->sym->name, "ALL_SURFACE_MOLECULES")==0)
      {
          rx->players[0]->flags |= CAN_GRIDWALL;
      }
      else if ((rx->players[0]->flags & NOT_FREE)==0)
      {
        /* two volume molecules */
        if ((rx->players[1]->flags & NOT_FREE)==0)
        {
          rx->players[0]->flags |= CAN_MOLMOL;
          rx->players[1]->flags |= CAN_MOLMOL;
        }
        /* one volume molecules and one wall */
        else if ((rx->players[1]->flags & IS_SURFACE)!=0)
        {
          rx->players[0]->flags |= CAN_MOLWALL;
        }
        /* one volume molecule and one grid molecule */
        else if ((rx->players[1]->flags & ON_GRID)!= 0)
        {
          rx->players[0]->flags |= CAN_MOLGRID;
        }
      }
      else if ((rx->players[0]->flags & IS_SURFACE)!=0)
      {
        /* one volume molecule and one wall */
        if ((rx->players[1]->flags & NOT_FREE)==0)
        {
          rx->players[1]->flags |= CAN_MOLWALL;
        }
        /* one grid molecule and one wall */
        else if ((rx->players[1]->flags & ON_GRID) != 0)
        {
          rx->players[1]->flags |= CAN_GRIDWALL;
        }
      }
      else if ((rx->players[0]->flags & ON_GRID)!= 0)
      {
        /* one volume molecule and one grid molecule */
        if ((rx->players[1]->flags & NOT_FREE)==0)
        {
          rx->players[1]->flags |= CAN_MOLGRID;
        }
        /* two grid molecules */
        else if ((rx->players[1]->flags & ON_GRID) != 0)
        {
          rx->players[0]->flags |= CAN_GRIDGRID;
          rx->players[1]->flags |= CAN_GRIDGRID;
        }
        /* one grid molecule and one wall */
        else if ((rx->players[1]->flags & IS_SURFACE) != 0)
        {
          rx->players[0]->flags |= CAN_GRIDWALL;
        }
      }
      break;

    case 3:
      if ((rx->players[2]->flags & IS_SURFACE) != 0)
      {
        /* two molecules and surface  */
        if ((rx->players[0]->flags & NOT_FREE)==0)
        {
          /* one volume molecule, one grid molecule, one surface */
          if ((rx->players[1]->flags & ON_GRID)!= 0)
          {
            rx->players[0]->flags |= CAN_MOLGRID;
          }
        }
        else if ((rx->players[0]->flags & ON_GRID)!= 0)
        {
          /* one volume molecule, one grid molecule, one surface */
          if ((rx->players[1]->flags & NOT_FREE)==0)
          {
            rx->players[1]->flags |= CAN_MOLGRID;
          }
          /* two grid molecules, one surface */
          else if ((rx->players[1]->flags & ON_GRID) != 0)
          {
            rx->players[0]->flags |= CAN_GRIDGRID;
            rx->players[1]->flags |= CAN_GRIDGRID;
          }
        }
      }
      else
      {
        if ((rx->players[0]->flags & NOT_FREE)==0)
        {
          if ((rx->players[1]->flags & NOT_FREE)==0)
          {
            /* three volume molecules */
            if ((rx->players[2]->flags & NOT_FREE)==0)
            {
              rx->players[0]->flags |= CAN_MOLMOLMOL;
              rx->players[1]->flags |= CAN_MOLMOLMOL;
              rx->players[2]->flags |= CAN_MOLMOLMOL;
            }
            /* two volume molecules and one grid molecule */
            else if ((rx->players[2]->flags & ON_GRID) !=0)
            {
              rx->players[0]->flags |= CAN_MOLMOLGRID;
              rx->players[1]->flags |= CAN_MOLMOLGRID;
            }
          }
          else if ((rx->players[1]->flags & ON_GRID) !=0)
          {
            /* two volume molecules and one grid molecule */
            if ((rx->players[2]->flags & NOT_FREE)==0)
            {
              rx->players[0]->flags |= CAN_MOLMOLGRID;
              rx->players[2]->flags |= CAN_MOLMOLGRID;
            }
            /* one volume molecules and two grid molecules */
            else if ((rx->players[2]->flags & ON_GRID) !=0)
            {
              rx->players[0]->flags |= CAN_MOLGRIDGRID;
            }
          }
        }
        else if ((rx->players[0]->flags & ON_GRID) != 0)
        {
          if ((rx->players[1]->flags & NOT_FREE)==0)
          {
            /* two volume molecules and one grid molecule */
            if ((rx->players[2]->flags & NOT_FREE)==0)
            {
              rx->players[1]->flags |= CAN_MOLMOLGRID;
              rx->players[2]->flags |= CAN_MOLMOLGRID;
            }
            /* one volume molecule and two grid molecules */
            else if ((rx->players[2]->flags & ON_GRID) !=0)
            {
              rx->players[1]->flags |= CAN_MOLGRIDGRID;
            }
          }
          else if ((rx->players[1]->flags & ON_GRID) !=0)
          {
            /* one volume molecule and two grid molecules */
            if ((rx->players[2]->flags & NOT_FREE)==0)
            {
              rx->players[2]->flags |= CAN_MOLGRIDGRID;
            }
            /* three grid molecules */
            else if ((rx->players[2]->flags & ON_GRID) !=0)
            {
              rx->players[0]->flags |= CAN_GRIDGRIDGRID;
              rx->players[1]->flags |= CAN_GRIDGRIDGRID;
              rx->players[2]->flags |= CAN_GRIDGRIDGRID;
            }
          }
        }
      }
      break;

    default:
      assert(0);
      break;
  }
}

/*************************************************************************
 build_reaction_hash_table:
    Scan the symbol table, copying all reactions found into the reaction hash.

 In:  parse_state: parser state
      num_rx: num reactions expected
 Out: 0 on success, 1 if we fail to allocate the table
*************************************************************************/
static int build_reaction_hash_table(struct mdlparse_vars *parse_state, int num_rx)
{
  struct rxn **rx_tbl = NULL;
  int rx_hash;
  for (rx_hash=2; rx_hash<=num_rx && rx_hash != 0; rx_hash <<= 1)
    ;
  rx_hash <<= 1;

  if (rx_hash == 0) rx_hash = MAX_RX_HASH;
  if (rx_hash > MAX_RX_HASH) rx_hash = MAX_RX_HASH;
#ifdef REPORT_RXN_HASH_STATS
  mcell_log("Num rxns: %d", num_rx);
  mcell_log("Size of hash: %d", rx_hash);
#endif

  /* Create the reaction hash table */
  parse_state->vol->rx_hashsize = rx_hash;
  rx_hash -= 1;
  rx_tbl = CHECKED_MALLOC_ARRAY(struct rxn*, parse_state->vol->rx_hashsize, "reaction hash table");
  if (rx_tbl==NULL)
     return 1;
  parse_state->vol->reaction_hash = rx_tbl;
  for (int i=0;i<=rx_hash;i++) rx_tbl[i] = NULL;

#ifdef REPORT_RXN_HASH_STATS
  int numcoll = 0;
#endif
  for (int i=0;i<parse_state->vol->rxn_sym_table->n_bins;i++)
  {
    for (struct sym_table *sym = parse_state->vol->rxn_sym_table->entries[i]; sym != NULL; sym = sym->next)
    {
      if (sym == NULL) continue;

      struct rxn *rx = (struct rxn*) sym->value;
      int table_slot;
      if (rx->n_reactants == 1)
      {
        table_slot = rx->players[0]->hashval & rx_hash;
      }
      else
      {
        table_slot = (rx->players[0]->hashval + rx->players[1]->hashval) & rx_hash;
      }


#ifdef REPORT_RXN_HASH_STATS
      if (rx_tbl[table_slot] != NULL)
      {
        mcell_log("Collision: %s and %s", rx_tbl[table_slot]->sym->name, sym->name);
        ++ numcoll;
      }
#endif
      parse_state->vol->n_reactions++;
      while (rx->next != NULL) rx = rx->next;
      rx->next = rx_tbl[table_slot];
      rx_tbl[table_slot] = (struct rxn*)sym->value;
    }
  }
#ifdef REPORT_RXN_HASH_STATS
  mcell_log("Num collisions: %d", numcoll);
#endif

  return 0;
}

/*************************************************************************
 prepare_reactions:
    Postprocess the parsed reactions, moving them to the reaction hash table,
    and transferring information from the pathway structures to a more compact,
    runtime-optimized form.

 In: parse_state: parser state
 Out: Returns 1 on error, 0 on success.
      Reaction hash table is built and geometries are set properly.  Unlike in
      the parser, reactions with different reactant geometries are _different
      reactions_, and are stored as separate struct rxns.

 Note: The user inputs _geometric equivalence classes_, but here we convert
       from that to _output construction geometry_.  A geometry of 0 means to
       choose a random orientation.  A geometry of k means to adopt the
       geometry of the k'th species in the list (reactants start at #1,
       products are in order after reactants).  A geometry of -k means to adopt
       the opposite of the geometry of the k'th species.  The first n_reactants
       products determine the fate of the reactants (NULL = destroyed), and the
       rest are real products.
 PostNote: The reactants are used for triggering, and those have equivalence
       class geometry even in here.
*************************************************************************/
int prepare_reactions(struct mdlparse_vars *parse_state)
{
  struct pathway *path;
  struct product *prod;
  struct rxn *rx;
  struct t_func *tp;
  //double D_tot, t_step;
  short geom;
  int k, kk;
  /* flags that tell whether reactant_1 is also on the product list,
     same for reactant_2 and reactant_3 */
  int recycled1, recycled2, recycled3;
  int num_rx, num_players;
  struct species *temp_sp;
  int n_prob_t_rxns; /* # of pathways with time-varying rates */
  struct rxn *reaction;

  num_rx = 0;

  parse_state->vol->vacancy_search_dist2 *= parse_state->vol->r_length_unit;         /* Convert units */
  parse_state->vol->vacancy_search_dist2 *= parse_state->vol->vacancy_search_dist2;  /* Take square */

  if (parse_state->vol->rx_radius_3d <= 0.0)
  {
    parse_state->vol->rx_radius_3d = 1.0/sqrt(MY_PI*parse_state->vol->grid_density);
  }
  parse_state->vol->tv_rxn_mem = create_mem(sizeof(struct t_func) , 100);
  if (parse_state->vol->tv_rxn_mem == NULL) return 1;

  for (int n_rxn_bin=0; n_rxn_bin<parse_state->vol->rxn_sym_table->n_bins; n_rxn_bin++)
  {
    for (struct sym_table *sym = parse_state->vol->rxn_sym_table->entries[n_rxn_bin];
         sym != NULL;
         sym = sym->next)
    {
      reaction = (struct rxn*)sym->value;
      reaction->next = NULL;

      for (path=reaction->pathway_head ; path != NULL ; path = path->next)
      {
        check_duplicate_special_reactions(path);

        /* if one of the reactants is a surface, move it to the last reactant.
         * Also arrange reactant1 and reactant2 in alphabetical order */
        if (reaction->n_reactants>1)
        {
          /* Put surface last */
          if ((path->reactant1->flags & IS_SURFACE) != 0)
          {
            temp_sp = path->reactant1;
            path->reactant1 = path->reactant2;
            path->reactant2 = temp_sp;
            geom = path->orientation1;
            path->orientation1 = path->orientation2;
            path->orientation2 = geom;
          }
          if (reaction->n_reactants>2)
          {
            if ((path->reactant2->flags & IS_SURFACE) != 0)
            {
              temp_sp = path->reactant3;
              path->reactant3 = path->reactant2;
              path->reactant2 = temp_sp;
              geom = path->orientation3;
              path->orientation3 = path->orientation2;
              path->orientation2 = geom;
            }
          }
          alphabetize_pathway(path, reaction);
        } /* end if (n_reactants > 1) */

      }  /* end for (path = reaction->pathway_head; ...) */


      /* if reaction contains equivalent pathways, split this reaction into a
       * linked list of reactions each containing only equivalent pathways.
       */

      rx = split_reaction(parse_state, reaction);

      /* set the symbol value to the head of the linked list of reactions */
      sym->value = (void *)rx;

      while (rx != NULL)
      {
        double pb_factor = 0.0;
        /* Check whether reaction contains pathways with equivalent product
         * lists.  Also sort pathways in alphabetical order according to the
         * "prod_signature" field.
         */
        check_reaction_for_duplicate_pathways(&rx->pathway_head);

        num_rx++;

        /* At this point we have reactions of the same geometry and can collapse them
         * and count how many non-reactant products are in each pathway. */

        /* Search for reactants that appear as products */
        /* Any reactants that don't appear are set to be destroyed. */
        rx->product_idx = CHECKED_MALLOC_ARRAY(u_int, rx->n_pathways+1, "reaction product index array");
        rx->cum_probs = CHECKED_MALLOC_ARRAY(double, rx->n_pathways, "reaction cumulative probabilities array");

        /* Note, that the last member of the array "rx->product_idx"
         * contains size of the array "rx->players" */

        if (rx->product_idx == NULL  || rx->cum_probs == NULL)
          return 1;

        if (reaction_has_complex_rates(rx))
        {
          int pathway_idx;
          rx->rates = CHECKED_MALLOC_ARRAY(struct complex_rate *, rx->n_pathways, "reaction complex rates array");
          if (rx->rates == NULL)
            return 1;
          for (pathway_idx = 0; pathway_idx < rx->n_pathways; ++ pathway_idx)
            rx->rates[pathway_idx] = NULL;
        }

        n_prob_t_rxns = 0;
        path = rx->pathway_head;

        for (int n_pathway=0; path!=NULL ; n_pathway++ , path = path->next)
        {

          rx->product_idx[n_pathway] = 0;
          if (rx->rates)
            rx->rates[n_pathway] = path->km_complex;

          /* Look for concentration clamp */
          if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE)!=0 &&
              path->km >= 0.0 && path->product_head==NULL && ((path->flags & PATHW_CLAMP_CONC) != 0))
          {
            struct ccn_clamp_data *ccd;

            if (n_pathway!=0 || path->next!=NULL)
              mcell_warn("Mixing surface modes with other surface reactions.  Please don't.");

            if (path->km>0)
            {
              ccd = CHECKED_MALLOC_STRUCT(struct ccn_clamp_data, "concentration clamp data");
              if (ccd==NULL)
                return 1;

              ccd->surf_class = path->reactant2;
              ccd->mol = path->reactant1;
              ccd->concentration = path->km;
              if (path->orientation1*path->orientation2==0)
              {
                ccd->orient = 0;
              }
              else
              {
                ccd->orient = (path->orientation1==path->orientation2) ? 1 : -1;
              }
              ccd->sides = NULL;
              ccd->next_mol = NULL;
              ccd->next_obj = NULL;
              ccd->objp = NULL;
              ccd->n_sides = 0;
              ccd->side_idx = NULL;
              ccd->cum_area = NULL;
              ccd->scaling_factor = 0.0;
              ccd->next = parse_state->vol->clamp_list;
              parse_state->vol->clamp_list = ccd;
            }
            path->km = GIGANTIC;
          }
          else if ((path->flags & PATHW_TRANSP) != 0)
          {
            rx->n_pathways = RX_TRANSP;
            if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE) &&
               (path->reactant1->flags & ON_GRID))
            {
                    path->reactant1->flags |= CAN_REGION_BORDER;
            }
          }
          else if ((path->flags & PATHW_REFLEC) != 0)
          {
            rx->n_pathways = RX_REFLEC;
            if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE) &&
               (path->reactant1->flags & ON_GRID))
            {
                    path->reactant1->flags |= CAN_REGION_BORDER;
            }
          }
          else if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE) && (path->reactant1->flags & ON_GRID) && (path->product_head==NULL) && (path->flags & PATHW_ABSORP))
          {
             rx->n_pathways = RX_ABSORB_REGION_BORDER;
             path->reactant1->flags |= CAN_REGION_BORDER;
          }
          else if ((strcmp(path->reactant1->sym->name, "ALL_SURFACE_MOLECULES") == 0))
          {
            if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE)  && (path->product_head==NULL) && (path->flags & PATHW_ABSORP))
            {
              rx->n_pathways = RX_ABSORB_REGION_BORDER;
              path->reactant1->flags |= CAN_REGION_BORDER;
            }
          }
          if (path->km_filename == NULL) rx->cum_probs[n_pathway] = path->km;
          else
          {
            rx->cum_probs[n_pathway]=0;
            n_prob_t_rxns++;
          }

          recycled1 = 0;
          recycled2 = 0;
          recycled3 = 0;

          for (prod=path->product_head ; prod != NULL ; prod = prod->next)
          {
            if (recycled1 == 0 && prod->prod == path->reactant1) recycled1 = 1;
            else if (recycled2 == 0 && prod->prod == path->reactant2) recycled2 = 1;
            else if (recycled3 == 0 && prod->prod == path->reactant3) recycled3 = 1;
            else rx->product_idx[n_pathway]++;
          }


        } /* end for (n_pathway=0,path=rx->pathway_head; ...) */

        /* Now that we know how many products there really are, set the index array */
        /* and malloc space for the products and geometries. */
        num_players = rx->n_reactants;
        kk = rx->n_pathways;
        if (kk<=RX_SPECIAL) kk = 1;
        for (int n_pathway=0;n_pathway<kk;n_pathway++)
        {
          k = rx->product_idx[n_pathway] + rx->n_reactants;
          rx->product_idx[n_pathway] = num_players;
          num_players += k;
        }
        rx->product_idx[kk] = num_players;

        rx->players = CHECKED_MALLOC_ARRAY(struct species*, num_players, "reaction players array");
        rx->geometries = CHECKED_MALLOC_ARRAY(short, num_players, "reaction geometries array");
        if (rx->pathway_head->is_complex[0] ||
            rx->pathway_head->is_complex[1] ||
            rx->pathway_head->is_complex[2])
        {
          rx->is_complex = CHECKED_MALLOC_ARRAY(unsigned char, num_players, "reaction 'is complex' flag");
          if (rx->is_complex == NULL)
            return 1;
          memset(rx->is_complex, 0, sizeof(unsigned char) * num_players);
        }
        else
          rx->is_complex = NULL;

        if (rx->players==NULL || rx->geometries==NULL)
          return 1;

        /* Load all the time-varying rates from disk (if any), merge them into */
        /* a single sorted list, and pull off any updates for time zero. */
        if (n_prob_t_rxns > 0)
        {
          path = rx->pathway_head;
          for (int n_pathway = 0; path!=NULL ; n_pathway++, path=path->next)
          {
            if (path->km_filename != NULL)
            {
              if (load_rate_file(parse_state, rx, path->km_filename, n_pathway))
                mcell_error("Failed to load rates from file '%s'.", path->km_filename);
            }
          }
          rx->prob_t = (struct t_func*) ae_list_sort((struct abstract_element*)rx->prob_t);

          while (rx->prob_t != NULL && rx->prob_t->time <= 0.0)
          {
            rx->cum_probs[ rx->prob_t->path ] = rx->prob_t->value;
            rx->prob_t = rx->prob_t->next;
          }
        } /* end if (n_prob_t_rxns > 0) */


        /* Set the geometry of the reactants.  These are used for triggering.                 */
        /* Since we use flags to control orientation changes, just tell everyone to stay put. */
        path = rx->pathway_head;
        rx->players[0] = path->reactant1;
        rx->geometries[0] = path->orientation1;
        if (rx->is_complex) rx->is_complex[0] = path->is_complex[0];
        if (rx->n_reactants > 1)
        {
          rx->players[1] = path->reactant2;
          rx->geometries[1] = path->orientation2;
          if (rx->is_complex) rx->is_complex[1] = path->is_complex[1];
          if (rx->n_reactants > 2)
          {
            rx->players[2] = path->reactant3;
            rx->geometries[2] = path->orientation3;
            if (rx->is_complex) rx->is_complex[2] = path->is_complex[2];
          }
        }

        /* maximum number of surface products */
        path = rx->pathway_head;
        int max_num_surf_products = set_product_geometries(path, rx, prod);

        pb_factor = compute_pb_factor(parse_state->vol, rx, max_num_surf_products);
        rx->pb_factor = pb_factor;
        path = rx->pathway_head;

        if (scale_probabilities(path, rx, parse_state, pb_factor))
          return 1;

        if (n_prob_t_rxns > 0)
        {
          for (tp = rx->prob_t ; tp != NULL ; tp = tp->next)
            tp->value *= pb_factor;
        }

        /* Move counts from list into array */
        if (rx->n_pathways > 0)
        {
          rx->info = CHECKED_MALLOC_ARRAY(struct pathway_info, rx->n_pathways, "reaction pathway info");
          if (rx->info == NULL)
            return 1;

          path = rx->pathway_head;
          for (int n_pathway=0; path!=NULL ; n_pathway++,path=path->next)
          {
            rx->info[n_pathway].count = 0;
            rx->info[n_pathway].pathname = path->pathname;    /* Keep track of named rxns */
            if (path->pathname!=NULL)
            {
              rx->info[n_pathway].pathname->path_num = n_pathway;
              rx->info[n_pathway].pathname->rx = rx;
            }
          }
        }
        else /* Special reaction, only one exit pathway */
        {
          rx->info = CHECKED_MALLOC_STRUCT(struct pathway_info,
                                               "reaction pathway info");
          if (rx->info == NULL)
            return 1;
          rx->info[0].count = 0;
          rx->info[0].pathname = rx->pathway_head->pathname;
          if (rx->pathway_head->pathname!=NULL)
          {
            rx->info[0].pathname->path_num = 0;
            rx->info[0].pathname->rx = rx;
          }
        }

        /* Sort pathways so all fixed pathways precede all varying pathways */
        if (rx->rates  &&  rx->n_pathways > 0)
          reorder_varying_pathways(rx);

        /* Compute cumulative properties */
        for (int n_pathway=1; n_pathway<rx->n_pathways; ++n_pathway)
          rx->cum_probs[n_pathway] += rx->cum_probs[n_pathway-1];
        if (rx->n_pathways > 0)
          rx->min_noreaction_p = rx->max_fixed_p = rx->cum_probs[rx->n_pathways - 1];
        else
          rx->min_noreaction_p = rx->max_fixed_p = 1.0;
        if (rx->rates)
          for (int n_pathway=0; n_pathway<rx->n_pathways; ++n_pathway)
            if (rx->rates[n_pathway])
              rx->min_noreaction_p += macro_max_rate(rx->rates[n_pathway], pb_factor);

        rx = rx->next;
      }
    }
  }

  if (parse_state->vol->grid_grid_reaction_flag  || parse_state->vol->grid_grid_grid_reaction_flag)
  {
    if (parse_state->vol->notify->reaction_probabilities==NOTIFY_FULL)
      mcell_log("For reaction between two (or three) surface molecules the upper probability limit is given. The effective reaction probability will be recalculated dynamically during simulation.");
  }

  if (build_reaction_hash_table(parse_state, num_rx))
    return 1;

  parse_state->vol->rx_radius_3d *= parse_state->vol->r_length_unit; /* Convert into length units */

  for (int n_rxn_bin=0;n_rxn_bin<parse_state->vol->rx_hashsize;n_rxn_bin++)
  {
    for (struct rxn *this_rx = parse_state->vol->reaction_hash[n_rxn_bin];
         this_rx != NULL;
         this_rx = this_rx->next)
    {
      /* Here we deallocate some memory used for creating pathways.
         Other pathways related memory will be freed in
         'mdlparse.y'.  */
      for (path = this_rx->pathway_head; path != NULL; path = path->next)
      {
        if (path->prod_signature != NULL) free(path->prod_signature);
      }

      set_reaction_player_flags(this_rx);
      this_rx->pathway_head = NULL;
    }
  }

  add_surface_reaction_flags(parse_state);

  if (parse_state->vol->notify->reaction_probabilities==NOTIFY_FULL)
    mcell_log_raw("\n");

  return 0;
}

/*************************************************************************
 warn_about_high_rates:
    If HIGH_REACTION_PROBABILITY is set to WARNING or ERROR, and the reaction
    probability is high, give the user a warning or error respectively.

 In: parse_state: parser state
     warn_file: The log/error file. Can be stdout/stderr 
     rate_warn: If 1, warn the user about high reaction rates (or give error)
     print_once: If the warning has been printed once, don't repeat it
 Out: print_once. Also print out reaction probabilities (with warning/error)
*************************************************************************/
int warn_about_high_rates(struct mdlparse_vars *parse_state, FILE *warn_file, int rate_warn, int print_once)
{
  if (rate_warn)
  {
    if (parse_state->vol->notify->high_reaction_prob==WARN_ERROR)
    {
      warn_file = mcell_get_error_file();
      if (!print_once)
      {
        fprintf(warn_file, "\n");
        fprintf(warn_file, "Reaction probabilities generated for the following reactions:\n");
        print_once = 1;
      }
      fprintf(warn_file,"\tError: High ");
    }
    else
    {
      if (!print_once)
      {
        fprintf(warn_file, "\n");
        fprintf(warn_file, "Reaction probabilities generated for the following reactions:\n");
        print_once = 1;
      }
      if (parse_state->vol->notify->high_reaction_prob==WARN_WARN) fprintf(warn_file,"\tWarning: High ");
      else fprintf(warn_file,"\t");
    }
  }
  else 
  {
      if (!print_once)
      {
        fprintf(warn_file, "\n");
        fprintf(warn_file, "Reaction probabilities generated for the following reactions:\n");
        print_once = 1;
      }
      fprintf(warn_file,"\t");
  }
  return print_once;
}

/*************************************************************************
 alphabetize_pathway:
    The reaction pathway (path) is alphabetized.

 In: path: Parse-time structure for reaction pathways
     reaction: Reaction pathways leading away from a given intermediate
 Out: Nothing. 
*************************************************************************/
void alphabetize_pathway(struct pathway *path, struct rxn *reaction)
{
  unsigned char temp_is_complex;
  short geom, geom2;
  struct species *temp_sp, *temp_sp2;

  /* Alphabetize if we have two molecules */
  if ((path->reactant2->flags&IS_SURFACE)==0)
  {
    if (strcmp(path->reactant1->sym->name, path->reactant2->sym->name) > 0)
    {
      temp_sp = path->reactant1;
      path->reactant1 = path->reactant2;
      path->reactant2 = temp_sp;
      geom = path->orientation1;
      path->orientation1 = path->orientation2;
      path->orientation2 = geom;
      temp_is_complex = path->is_complex[0];
      path->is_complex[0] = path->is_complex[1];
      path->is_complex[1] = temp_is_complex;
    }
    else if (strcmp(path->reactant1->sym->name, path->reactant2->sym->name) == 0)
    {
      if (path->orientation1 < path->orientation2)
      {
        geom = path->orientation1;
        path->orientation1 = path->orientation2;
        path->orientation2 = geom;
        temp_is_complex = path->is_complex[0];
        path->is_complex[0] = path->is_complex[1];
        path->is_complex[1] = temp_is_complex;
      }
    }
  }

  /* Alphabetize if we have three molecules */
  if (reaction->n_reactants == 3)
  {
    if ((path->reactant3->flags&IS_SURFACE)==0)
    {
      if (strcmp(path->reactant1->sym->name, path->reactant3->sym->name) > 0)
      {
         /* Put reactant3 at the beginning */
         temp_sp = path->reactant1;
         geom = path->orientation1;
         path->reactant1 = path->reactant3;
         path->orientation1 = path->orientation3;

         /* Put former reactant1 in place of reactant2 */
         temp_sp2 = path->reactant2;
         geom2 = path->orientation2;
         path->reactant2 = temp_sp;
         path->orientation2 = geom;

         /* Put former reactant2 in place of reactant3 */
         path->reactant3 = temp_sp2;
         path->orientation3 = geom2;
         /* XXX: Update to deal with macromolecules? */

      }
      else if (strcmp(path->reactant2->sym->name, path->reactant3->sym->name) > 0)
      {

         /* Put reactant3 after reactant1 */
         temp_sp = path->reactant2;
         path->reactant2 = path->reactant3;
         path->reactant3 = temp_sp;
         geom = path->orientation2;
         path->orientation2 = path->orientation3;
         path->orientation3 = geom;

      }
    } /*end */
  }
}

/*************************************************************************
 check_duplicate_special_reactions:
   Check for duplicate special reaction pathways (e.g. TRANSPARENT = molecule).

 In: path: Parse-time structure for reaction pathways
 Out: Nothing. 
 Note: I'm not sure if this code is ever actually called.
*************************************************************************/
void check_duplicate_special_reactions(struct pathway *path)
{
  /* if it is a special reaction - check for the duplicates pathways */
  if (path->next != NULL)
  {
    if ((path->flags & PATHW_TRANSP) && (path->next->flags & PATHW_TRANSP))
    {
      if ((path->orientation2 == path->next->orientation2) ||
         (path->orientation2 == 0) || (path->next->orientation2 == 0))
      {
         mcell_error("Exact duplicates of special reaction TRANSPARENT = %s are not allowed.  Please verify the contents of DEFINE_SURFACE_CLASS statement.", path->reactant2->sym->name);
      }
    }

    if ((path->flags & PATHW_REFLEC) && (path->next->flags & PATHW_REFLEC))
    {
      if ((path->orientation2 == path->next->orientation2) ||
         (path->orientation2 == 0) || (path->next->orientation2 == 0))
      {
         mcell_error("Exact duplicates of special reaction REFLECTIVE = %s are not allowed.  Please verify the contents of DEFINE_SURFACE_CLASS statement.", path->reactant2->sym->name);
      }
    }
    if ((path->flags & PATHW_ABSORP) && (path->next->flags & PATHW_ABSORP))
    {
      if ((path->orientation2 == path->next->orientation2) ||
         (path->orientation2 == 0) || (path->next->orientation2 == 0))
      {
        mcell_error("Exact duplicates of special reaction ABSORPTIVE = %s are not allowed.  Please verify the contents of DEFINE_SURFACE_CLASS statement.", path->reactant2->sym->name);
      }
    }
  }
}


/*************************************************************************
 set_product_geometries:
 
  Walk through the list, setting the geometries of each of the products. We do
  this by looking for an earlier geometric match and pointing there or we just
  point to 0 if there is no match.

 In: path: Parse-time structure for reaction pathways
     rx: Pathways leading away from a given intermediate
     prod: Parse-time structure for products of reaction pathways
 Out: max_num_surf_products: Maximum number of surface products 
*************************************************************************/
int set_product_geometries(struct pathway *path, struct rxn *rx, struct product *prod)
{
  int recycled1, recycled2, recycled3;
  int k, kk, k2;
  short geom;
  struct product *prod2;
  int max_num_surf_products;         /* maximum number of surface products */
  int num_surf_products_per_pathway;

  max_num_surf_products = 0;
  for (int n_pathway=0; path!=NULL ; n_pathway++ , path = path->next)
  {
    recycled1 = 0;
    recycled2 = 0;
    recycled3 = 0;
    k = rx->product_idx[n_pathway] + rx->n_reactants;
    num_surf_products_per_pathway = 0;
    for (prod=path->product_head ; prod != NULL ; prod = prod->next)
    {
      if (recycled1==0 && prod->prod == path->reactant1)
      {
        recycled1 = 1;
        kk = rx->product_idx[n_pathway] + 0;
      }
      else if (recycled2==0 && prod->prod == path->reactant2)
      {
        recycled2 = 1;
        kk = rx->product_idx[n_pathway] + 1;
      }
      else if (recycled3==0 && prod->prod == path->reactant3)
      {
        recycled3 = 1;
        kk = rx->product_idx[n_pathway] + 2;
      }
      else
      {
        kk = k;
        k++;
      }

      if (prod->prod->flags & ON_GRID) num_surf_products_per_pathway++;

      rx->players[kk] = prod->prod;
      if (rx->is_complex) rx->is_complex[kk] = prod->is_complex;

      if ((prod->orientation+path->orientation1)*(prod->orientation-path->orientation1)==0 && prod->orientation*path->orientation1!=0)
      {
        if (prod->orientation == path->orientation1) rx->geometries[kk] = 1;
        else rx->geometries[kk] = -1;
      }
      else if (rx->n_reactants > 1 &&
                (prod->orientation+path->orientation2)*(prod->orientation-path->orientation2)==0 && prod->orientation*path->orientation2!=0
             )
      {
        if (prod->orientation == path->orientation2) rx->geometries[kk] = 2;
        else rx->geometries[kk] = -2;
      }
      else if (rx->n_reactants > 2 &&
                (prod->orientation+path->orientation3)*(prod->orientation-path->orientation3)==0 && prod->orientation*path->orientation3!=0
             )
      {
        if (prod->orientation == path->orientation3) rx->geometries[kk] = 3;
        else rx->geometries[kk] = -3;
      }
      else
      {
        k2 = 2*rx->n_reactants + 1;  /* Geometry index of first non-reactant product, counting from 1. */
        geom = 0;
        for (prod2=path->product_head ; prod2!=prod && prod2!=NULL && geom==0 ; prod2 = prod2->next)
        {
          if ((prod2->orientation+prod->orientation)*(prod2->orientation-prod->orientation)==0 && prod->orientation*prod2->orientation!=0)
          {
            if (prod2->orientation == prod->orientation) geom = 1;
            else geom = -1;
          }
          else geom = 0;

          if (recycled1 == 1)
          {
            if (prod2->prod == path->reactant1)
            {
              recycled1 = 2;
              geom *= rx->n_reactants+1;
            }
          }
          else if (recycled2==1)
          {
            if (prod2->prod == path->reactant2)
            {
              recycled2 = 2;
              geom *= rx->n_reactants+2;
            }
          }
          else if (recycled3==1)
          {
            if (prod2->prod == path->reactant3)
            {
              recycled3 = 2;
              geom *= rx->n_reactants+3;
            }
          }
          else
          {
            geom *= k2;
            k2++;
          }
        }
        rx->geometries[kk] = geom;
      }
      if (num_surf_products_per_pathway > max_num_surf_products) max_num_surf_products = num_surf_products_per_pathway;
    }

    k = rx->product_idx[n_pathway];
    if (recycled1==0) rx->players[k] = NULL;
    if (recycled2==0 && rx->n_reactants>1) rx->players[k+1] = NULL;
    if (recycled3==0 && rx->n_reactants>2) rx->players[k+2] = NULL;
  } /* end for (n_pathway = 0, ...) */
  return max_num_surf_products;
}


/*************************************************************************
 scale_probabilities:
 
  Scale probabilities, notifying and warning as appropriate.

 In: path: Parse-time structure for reaction pathways
     rx: Pathways leading away from a given intermediate
     parse_state: parser state
     pb_factor:
 Out: Return 1 if rates are high and HIGH_REACTION_PROBABILITY is set to ERROR
 Note: This does not work properly right now. Even if rates are high and
       HIGH_REACTION_PROBABILITY is set to ERROR, the error is ignored
*************************************************************************/
int scale_probabilities(struct pathway *path, struct rxn *rx, struct mdlparse_vars *parse_state, double pb_factor)
{
  int print_once = 0;  /* flag */
  FILE *warn_file;
  int is_gigantic;
  double rate;

  for (int n_pathway=0;path != NULL;n_pathway++, path = path->next)
  {
    int rate_notify=0, rate_warn=0;
    if (rx->cum_probs[n_pathway]==GIGANTIC) is_gigantic=1;
    else is_gigantic=0;

    /* automatic surface reactions will be printed out from 'init_sim()'. */
    if (is_gigantic) continue;

    if (! rx->rates  ||  ! rx->rates[n_pathway])
      rate = pb_factor*rx->cum_probs[n_pathway];
    else
      rate = 0.0;
    rx->cum_probs[n_pathway] = rate;

    if ((parse_state->vol->notify->reaction_probabilities==NOTIFY_FULL && ((rate>=parse_state->vol->notify->reaction_prob_notify) || (parse_state->vol->notify->reaction_prob_notify==0.0))))
      rate_notify = 1;
    if ((parse_state->vol->notify->high_reaction_prob != WARN_COPE && ((rate>=parse_state->vol->notify->reaction_prob_warn) || ((parse_state->vol->notify->reaction_prob_warn==0.0)))))
      rate_warn = 1;

    if ((rate > 1.0) && (!parse_state->vol->reaction_prob_limit_flag))
    {
      parse_state->vol->reaction_prob_limit_flag = 1;
    }


    if (rate_warn || rate_notify)
    {

      warn_file = mcell_get_log_file();

      print_once = warn_about_high_rates(parse_state, warn_file, rate_warn, print_once);

      if (rx->rates  &&  rx->rates[n_pathway])
        fprintf(warn_file,"Varying probability \"%s\" set for ", rx->rates[n_pathway]->name);
      else
        fprintf(warn_file,"Probability %.4e set for ",rate);
      if (rx->n_reactants==1) fprintf(warn_file,"%s{%d} -> ",rx->players[0]->sym->name,rx->geometries[0]);
      else if (rx->n_reactants == 2)
      {
        if (rx->players[1]->flags & IS_SURFACE)
        {
          fprintf(warn_file,"%s{%d} @ %s{%d} -> ",
                  rx->players[0]->sym->name,rx->geometries[0],
                  rx->players[1]->sym->name,rx->geometries[1]);
         }
         else
         {
           fprintf(warn_file,"%s{%d} + %s{%d} -> ",
                   rx->players[0]->sym->name,rx->geometries[0],
                   rx->players[1]->sym->name,rx->geometries[1]);
         }
      }
      else
      {
        if (rx->players[2]->flags & IS_SURFACE)
        {
          fprintf(warn_file,"%s{%d} + %s{%d}  @ %s{%d} -> ",
                  rx->players[0]->sym->name,rx->geometries[0],
                  rx->players[1]->sym->name,rx->geometries[1],
                  rx->players[2]->sym->name,rx->geometries[2]);
        }
        else
        {
          fprintf(warn_file,"%s{%d} + %s{%d}  + %s{%d} -> ",
                  rx->players[0]->sym->name,rx->geometries[0],
                  rx->players[1]->sym->name,rx->geometries[1],
                  rx->players[2]->sym->name,rx->geometries[2]);
        }
      }
      if (path->product_head == NULL)
      {
        fprintf(warn_file,"NULL ");
      }
      else
      {
        for (struct product *prod = path->product_head ; prod != NULL ; prod = prod->next)
        {
         fprintf(warn_file,"%s{%d} ",prod->prod->sym->name, prod->orientation);
        }
      }

      fprintf(warn_file,"\n");

      if (rate_warn && parse_state->vol->notify->high_reaction_prob==WARN_ERROR)
        return 1;
    }
  }
  return 0;
}


/*************************************************************************
 add_surface_reaction_flags:
 
 In: parse_state: parser state
 Out: Nothing
*************************************************************************/
void add_surface_reaction_flags(struct mdlparse_vars *parse_state)
{
  struct species *temp_sp;

  /* Add flags for surface reactions with ALL_MOLECULES */
  if (parse_state->vol->all_mols->flags & (CAN_MOLWALL|CAN_GRIDWALL))
  {
    for (int n_mol_bin=0; n_mol_bin<parse_state->vol->mol_sym_table->n_bins; n_mol_bin++)
    {
      for (struct sym_table *symp = parse_state->vol->mol_sym_table->entries[n_mol_bin];
           symp != NULL;
           symp = symp->next)
      {
        temp_sp = (struct species*) symp->value;
        if (temp_sp == parse_state->vol->all_mols) continue;
        if (temp_sp == parse_state->vol->all_volume_mols) continue;
        if (temp_sp == parse_state->vol->all_surface_mols) continue;

        if (((temp_sp->flags & NOT_FREE) == 0) && ((temp_sp->flags & CAN_MOLWALL) == 0))
        {
          temp_sp->flags |= CAN_MOLWALL;
        }
        else if ((temp_sp->flags & ON_GRID) && ((temp_sp->flags & CAN_REGION_BORDER) == 0))
        {
          temp_sp->flags |= CAN_REGION_BORDER;
        }
      }
    }
  }

  /* Add flags for surface reactions with ALL_VOLUME_MOLECULES */
  if (parse_state->vol->all_volume_mols->flags & CAN_MOLWALL)
  {
    for (int n_mol_bin=0; n_mol_bin<parse_state->vol->mol_sym_table->n_bins; n_mol_bin++)
    {
      for (struct sym_table *symp = parse_state->vol->mol_sym_table->entries[n_mol_bin];
           symp != NULL;
           symp = symp->next)
      {
        temp_sp = (struct species*) symp->value;
        if (temp_sp == parse_state->vol->all_mols) continue;
        if (temp_sp == parse_state->vol->all_volume_mols) continue;
        if (temp_sp == parse_state->vol->all_surface_mols) continue;
        if (((temp_sp->flags & NOT_FREE) == 0) && ((temp_sp->flags & CAN_MOLWALL) == 0))
        {
          temp_sp->flags |= CAN_MOLWALL;
        }
      }
    }
  }

  /* Add flags for surface reactions with ALL_SURFACE_MOLECULES */
  if (parse_state->vol->all_surface_mols->flags & CAN_GRIDWALL)
  {
    for (int n_mol_bin=0; n_mol_bin<parse_state->vol->mol_sym_table->n_bins; n_mol_bin++)
    {
      for (struct sym_table *symp = parse_state->vol->mol_sym_table->entries[n_mol_bin];
           symp != NULL;
           symp = symp->next)
      {
        temp_sp = (struct species*) symp->value;
        if (temp_sp == parse_state->vol->all_mols) continue;
        if (temp_sp == parse_state->vol->all_volume_mols) continue;
        if (temp_sp == parse_state->vol->all_surface_mols) continue;
        if (((temp_sp->flags & ON_GRID) && ((temp_sp->flags & CAN_REGION_BORDER) == 0)))
        {
          temp_sp->flags |= CAN_REGION_BORDER;
        }
      }
    }
  }
}
