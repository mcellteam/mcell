#ifndef INCLUDED_VERSION_INFO_H
#define INCLUDED_VERSION_INFO_H

#include <stdio.h>

/* MCell version as a string */
extern char const mcell_version[];

/* Write the credits to a file handle */
void print_credits(FILE *f);

/* Write the version info to a file handle */
void print_version(FILE *f);

/* Write the version info to a file handle */
void print_full_version(FILE *f);

#endif /* INCLUDED_VERSION_INFO_H */
