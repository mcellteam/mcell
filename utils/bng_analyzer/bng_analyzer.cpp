/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#include <getopt.h>
#include <iostream>
#include <cassert>
#include <string>

#include "nfsim_species_unifier.h"


using namespace std;

/* Command-line arguments structure:
 *     long arg name
 *     has argument
 *     pointer to flag (should always be 0)
 *     short argument letter
 */
static const option long_options[] = {
    { "help", 0, 0, 'h' },
    { "version", 0, 0, 'v' },
    { nullptr, 0, 0, 0 }
};


void print_usage(const char* argv0) {
  cout << "TODO\n";
}


void print_version(const char* argv0) {
  cout << "TODO\n";
}


const int ARG_PARSE_ERROR = 1;
const int ARG_PARSE_QUIT = 0;
const int ARG_PARSE_OK = -1;
// returns ARG_PARSE_OK if execution should continue
// ARG_PARSE_END to end with exit code 0,
// ARG_PARSE_ERROR to end with exit code 1,
int process_args(
    const int argc, char* argv[],
    string& input_file
) {
  input_file = "";

  assert(argc > 0);
  while (1) {

    // get the next argument
    int c = getopt_long_only(argc, argv, "hvgo:", long_options, nullptr);
    if (c == -1)
      break;

    switch (c) {
      case 'h':
        print_usage(argv[0]);
        return ARG_PARSE_QUIT;
      case 'v':
        print_version(argv[0]);
        return ARG_PARSE_QUIT;
    }
  }

  if (optind < argc) {
    if (argc - optind > 1) {
      cerr << "Only one input file can be specified.\n";
      return ARG_PARSE_ERROR;
    }
    input_file = argv[optind];
  }
  else {
    cerr << "Input file was not specified.\n";
    return ARG_PARSE_ERROR;
  }

  return ARG_PARSE_OK;
}


int main(const int argc, char* argv[]) {

  string input_file;

  int arg_process_res = process_args(argc, argv, input_file);
  if (arg_process_res != ARG_PARSE_OK) {
    return arg_process_res;
  }

  MCell::NFSimSpeciesUnifier nfsim_unifier;
  bool ok = nfsim_unifier.read_species_file(input_file);
  if (!ok) {
    cerr << "There was an error while reading " << input_file << ".\n";
    return 1;
  }
  nfsim_unifier.print_unified_species();
  return 0;
}
