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

/* *nix-specific includes and defines */

#ifndef MCELL_CONFIG_NIX_H
#define MCELL_CONFIG_NIX_H

/* Macro for eliminating "unused variable" or "unused parameter" warnings. */
#define UNUSED(p) ((void)(p))

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

#define PRINTF_FORMAT(arg)                                                     \
  __attribute__((__format__(printf, arg, arg + 1))) /* specifies that a        \
                                                       function is uses a      \
                                                       printf-like format -    \
                                                       only used for           \
                                                       compile-time warnings   \
                                                       */
#define PRINTF_FORMAT_V(arg)                                                   \
  __attribute__((__format__(printf, arg, 0))) /* same but with va_list         \
                                                 argument */

#endif
