#ifndef ARGPARSE_H
#define ARGPARSE_H

#include "mcell_structs.h"

struct argparse_vars {
  int ival;
  double rval;
  char *cval;
  char *arg_err_msg;
  struct volume *vol;
};

int argparse(void *p);
void argerror(char *s,...);
int argparse_init(int argc, char *argv[], struct volume *vol);

#endif
