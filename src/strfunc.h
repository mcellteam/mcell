#ifndef STRFUNC_H
#define STRFUNC_H

#include <stdarg.h>

/*  Header file for character string handling functions */

char *alloc_vsprintf(char const *fmt, va_list args)
  PRINTF_FORMAT_V(1);
char *alloc_sprintf(char const *fmt, ...)
  PRINTF_FORMAT(1);
char *my_strcat(char const *s1, char const *s2);
char *my_strclump(char **slist);
char *strip_quotes(char const *s);

#endif
