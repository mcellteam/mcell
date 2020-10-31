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
    { "dat", 0, 0, 'd' },
    { "output", 1, 0, 'o' },
    { nullptr, 0, 0, 0 }
};


void print_usage(const char* argv0) {
  cout <<
      " [-d] INPUT_FILE [-o OUTPUT_FILE]\n"
      "  -d  Optional argument that specifies that the input is a .dat ascii file \n"
      "      for MCell visualization, NFSim .species file otherwise\n";
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
    string& input_file,
    string& output_file,
    bool& is_dat
) {
  input_file = "";
  output_file = "";
  is_dat = false;

  assert(argc > 0);
  while (1) {

    // get the next argument
    int c = getopt_long_only(argc, argv, "hvo:", long_options, nullptr);
    if (c == -1)
      break;

    switch (c) {
      case 'h':
        print_usage(argv[0]);
        return ARG_PARSE_QUIT;
      case 'v':
        print_version(argv[0]);
        return ARG_PARSE_QUIT;
      case 'o':
        output_file = optarg;
        break;
      case 'd':
        is_dat = true;
        break;
      default:
        cerr << "Invalid arguments.\n";
        return ARG_PARSE_ERROR;
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
  string output_file;
  bool is_dat;

  int arg_process_res = process_args(argc, argv, input_file, output_file, is_dat);
  if (arg_process_res != ARG_PARSE_OK) {
    return arg_process_res;
  }

  MCell::NFSimSpeciesUnifier nfsim_unifier;
  bool ok;
  ok = nfsim_unifier.read_input_file(input_file, is_dat);
  if (!ok) {
    cerr << "There was an error while reading " << input_file << ".\n";
    return 1;
  }
  ok = nfsim_unifier.print_unified_species(output_file);
  return ok ? 0 : 1;
}
