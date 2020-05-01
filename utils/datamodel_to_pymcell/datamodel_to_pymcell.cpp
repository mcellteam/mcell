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

#include "pymcell_converter.h"

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
    { "output_file_prefix", 1, 0, 'o'},
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
    string& input_file,
    string& output_files_prefix
) {
  assert(argc > 0);
  while (1) {

    // get the next argument
    int c = getopt_long_only(argc, argv, "?h", long_options, nullptr);
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
        output_files_prefix = optarg;
        break;
    }
  }

  if (optind < argc) {
    if (argc - optind > 1) {
      cerr << "Only one input datamodel file can be specified.\n";
      return ARG_PARSE_ERROR;
    }
    input_file = argv[optind];
  }
  else {
    cerr << "Input datamodel file file was not specified.\n";
    return ARG_PARSE_ERROR;
  }


  if (output_files_prefix == "") {
    cerr << "Output files prefix must specified.\n";
    return ARG_PARSE_ERROR;
  }

  return ARG_PARSE_OK;
}


int main(const int argc, char* argv[]) {

  string input_file;
  string output_files_prefix;

  int arg_process_res = process_args(argc, argv, input_file, output_files_prefix);
  if (arg_process_res != ARG_PARSE_OK) {
    return arg_process_res;
  }

  MCell::PymcellConverter converter;
  bool ok = converter.convert(input_file, output_files_prefix);

  if (!ok) {
    cerr << "There was an error while converting " << input_file << " to pymcell code.\n";
    return 1;
  }
  return 0;
}
