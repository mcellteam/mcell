#ifndef ARGPARSE_H
#define ARGPARSE_H

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

void argerror(struct volume *vol, char const *s, ...)
  __attribute__((format (printf, 2, 3)));

int argparse_init(int argc, char * const argv[], struct volume *vol);

#endif
