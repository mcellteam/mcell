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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _MSC_VER
#include <unistd.h>
#include <sys/time.h>
#endif
#include <time.h>

#include "config.h"
#include "version_info.h"
#include "version.h"

const char mcell_version[] = MCELL_VERSION;

/*
 * Prints out authors' institutions.
 *
 *    f: file handle to which to write
 */
void print_credits(FILE *f) {
  fprintf(
      f,
      "  Copyright (C) 2006-2018 by\n"
      "    The National Center for Multiscale Modeling of Biological Systems,\n"
      "    The Salk Institute for Biological Studies, and\n"
      "    Pittsburgh Supercomputing Center, Carnegie Mellon University,\n\n\n"
      "**********************************************************************\n"
      "MCell development is supported by the NIGMS-funded (P41GM103712)\n"
      "National Center for Multiscale Modeling of Biological Systems (MMBioS).\n"
      "Please acknowledge MCell in your publications.\n"
      "**********************************************************************"
      "\n\n");
}

/*
 * Prints out brief version info.
 *
 *    f: file handle to which to write
 */
void print_version(FILE *f) {
  char hostname[256];
  char curtime[128];
  time_t now = time(NULL);

  /* Print the version line */
  fprintf(f, "MCell %s", MCELL_VERSION);
  if (strncmp(MCELL_REVISION, "", strlen(MCELL_REVISION)) != 0)
    fprintf(f, " (commit: %s  date: %s)", MCELL_REVISION, MCELL_REVISION_DATE);
  if (!MCELL_REVISION_COMMITTED)
    fprintf(f, " [unofficial revision]");
  fprintf(f, "\n");

#ifndef _WIN64
  /* Print the current machine details */
  gethostname(hostname, sizeof(hostname));
  fprintf(f, "  Running on %s at %s\n", hostname, ctime_r(&now, curtime));
#endif
  print_credits(f);
}

/*
 * Prints out full version info.
 *
 *    f: file handle to which to write
 */
void print_full_version(FILE *f) {
  /* Print version line */
  print_version(f);

  /* Print build info */
  fprintf(f, "  Built at %s on %s by %s\n", MCELL_BUILDDATE, MCELL_BUILDHOST,
          MCELL_BUILDUSER);
  fprintf(f, "    Src directory: %s\n", MCELL_SRCDIR);
  fprintf(f, "    Build directory: %s\n", MCELL_BUILDDIR);
  fprintf(f, "    Branch: %s\n", MCELL_REVISION_BRANCH);
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
  if (MCELL_CFLAGS[0] != '\0' || MCELL_LDFLAGS[0] != '\0' ||
      MCELL_LFLAGS[0] != '\0' || MCELL_YFLAGS[0] != '\0') {
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
