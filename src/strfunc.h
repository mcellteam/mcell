#ifndef STRFUNC_H
#define STRFUNC_H

#include <stdarg.h>

/*  Header file for character string handling functions */

char *alloc_vsprintf(char const *fmt, va_list args);
char *alloc_sprintf(char const *fmt, ...)
  __attribute__((format (printf, 1, 2)));
char *my_strcat(char *s1, char *s2);
char *my_strclump(char **slist);
char *strip_quotes(char *s);

#endif
