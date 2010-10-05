#include "version_info.h"
#include "version.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>

const char mcell_version[] = MCELL_VERSION;

/*
 * Prints out authors' institutions.  Ordering between the two institutions is
 * chosen randomly.
 *
 *    f: file handle to which to write
 */
void print_credits(FILE *f)
{
  char const *institute[2];

  srand(time(NULL) / 2);
  if (rand() & 1)
  {
    institute[0]="The Salk Institute for Biological Studies";
    institute[1]="Pittsburgh Supercomputing Center, Carnegie Mellon University";
  }
  else
  {
    institute[0]="Pittsburgh Supercomputing Center, Carnegie Mellon University";
    institute[1]="The Salk Institute for Biological Studies";
  }

  fprintf(f,"  Copyright (C) 2006 - 2010 by\n    %s and\n    %s\n\n", institute[0], institute[1]);
}

/*
 * Prints out brief version info.
 *
 *    f: file handle to which to write
 */
void print_version(FILE *f)
{
  char hostname[256];
  char curtime[128];
  time_t now = time(NULL);

  /* Print the version line */
  fprintf(f, "MCell %s", MCELL_VERSION);
  if (MCELL_REVISION != -1)
    fprintf(f, " (revision %d/%s)", MCELL_REVISION, MCELL_REVISION_DATE);
  if (! MCELL_REVISION_COMMITTED)
    fprintf(f, " [unofficial revision]");
  fprintf(f, "\n");

  /* Print the current machine details */
  gethostname(hostname, 256);
  fprintf(f,"  Running on %s at %s\n", hostname, ctime_r(&now, curtime));
  print_credits(f);
}

/*
 * Prints out full version info.
 *
 *    f: file handle to which to write
 */
void print_full_version(FILE *f)
{
  /* Print version line */
  print_version(f);

  /* Print build info */
  fprintf(f, "  Built at %s on %s by %s\n",
             MCELL_BUILDDATE,
             MCELL_BUILDHOST,
             MCELL_BUILDUSER);
  fprintf(f, "    Src directory: %s\n", MCELL_SRCDIR);
  fprintf(f, "    Build directory: %s\n", MCELL_BUILDDIR);
  fprintf(f, "    Machine info: %s\n", MCELL_BUILDUNAME);
  fprintf(f, "\n");

  /* Print detailed build options */
  fprintf(f, "  Build configuration:\n");
  fprintf(f, "    Compiler:  %s [%s]\n", MCELL_CC, MCELL_CC_PATH);
  fprintf(f, "      version: %s\n", MCELL_CC_VERSION);
  fprintf(f, "    Linker:    %s [%s]\n", MCELL_LD, MCELL_LD_PATH);
  fprintf(f, "      version: %s\n", MCELL_LD_VERSION);
  fprintf(f, "    lex:       %s [%s]\n", MCELL_FLEX, MCELL_FLEX_PATH);
  fprintf(f, "      version: %s\n", MCELL_FLEX_VERSION);
  fprintf(f, "    yacc:      %s [%s]\n", MCELL_BISON, MCELL_BISON_PATH);
  fprintf(f, "      version: %s\n", MCELL_BISON_VERSION);
  fprintf(f, "\n");
  if (MCELL_CFLAGS[0] != '\0'  || MCELL_LDFLAGS[0] != '\0'  ||
      MCELL_LFLAGS[0] != '\0'  || MCELL_YFLAGS[0] != '\0')
  {
    fprintf(f, "  Build flags:\n");
    if (MCELL_CFLAGS[0] != '\0')
      fprintf(f, "    Compiler flags:  %s\n", MCELL_CFLAGS);
    if (MCELL_LDFLAGS[0] != '\0')
      fprintf(f, "    Linker flags:    %s\n", MCELL_LDFLAGS);
    if (MCELL_LFLAGS[0] != '\0')
      fprintf(f, "    lex flags:       %s\n", MCELL_LFLAGS);
    if (MCELL_YFLAGS[0] != '\0')
      fprintf(f, "    yacc flags:      %s\n", MCELL_YFLAGS);
    fprintf(f, "\n");
  }
}
