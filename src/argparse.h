#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <stdio.h>

struct volume;

void print_usage(FILE *f, char const *argv0);

int argparse_init(int argc, char * const argv[], struct volume *vol);

#endif
