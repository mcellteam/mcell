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

#pragma once

#include <stdarg.h>

/*  Header file for character string handling functions */

char *alloc_vsprintf(char const *fmt, va_list args) PRINTF_FORMAT_V(1);
char *alloc_sprintf(char const *fmt, ...) PRINTF_FORMAT(1);
char *my_strcat(char const *s1, char const *s2);
char *strip_quotes(char const *s);
