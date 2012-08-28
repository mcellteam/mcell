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

#endif

