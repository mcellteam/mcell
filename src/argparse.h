#ifndef ARGPARSE_H
#define ARGPARSE_H

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

struct volume;

void print_usage(FILE *f, char const *argv0);

void argerror(struct volume *vol, char const *s, ...)
  __attribute__((format (printf, 2, 3)));

int argparse_init(int argc, char * const argv[], struct volume *vol);

#endif
