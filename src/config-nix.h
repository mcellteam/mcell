/* *nix-specific includes and defines */

#ifndef MCELL_CONFIG_NIX_H
#define MCELL_CONFIG_NIX_H

/* Macro for eliminating "unused variable" or "unused parameter" warnings. */
#define UNUSED(p) ((void) (p))

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

#define PRINTF_FORMAT(arg) __attribute__((__format__(printf, arg, arg+1))) /* specifies that a function is uses a printf-like format - only used for compile-time warnings */
#define PRINTF_FORMAT_V(arg) __attribute__((__format__(printf, arg, 0))) /* same but with va_list argument */

#endif

