/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef _MSC_VER
#include <getopt.h>
#else
#include "win_getopt/win_getopt.h"
#endif

#include <iostream>
#include <cassert>
#include <string>
#include <ctype.h>

#include "mcell4_generator.h"
#include "generator_structs.h"

#include "version.h"

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
    { "debug", 0, 0, 'g' },
    { "testing", 0, 0, 't' },
    { "checkpoint_iters", 1, 0, 'k' },
    { "bng", 0, 0, 'b' },
    { "not_overridable_python_params", 0, 0, 'p'},
    { "output_file_prefix", 1, 0, 'o'},
    { nullptr, 0, 0, 0 }
};


void print_usage(const char* argv0) {
  cout << argv0 << " options:\n";

  size_t i = 0;
  while (long_options[i].name != nullptr) {
    cout << "  -" << (char)long_options[i].val << "| -" << long_options[i].name << "\n";
    i++;
  }
}


void print_version(const char* argv0) {
  cout << argv0 << " version:" << MCELL_VERSION << "\n";
}


// iters_arg is a list of integers separated by comma
std::vector<int> parse_checkpoint_iters(char* iters_arg) {
  assert(iters_arg != nullptr);

  std::vector<int> res;
  string current;
  char* c = iters_arg;

  while (*c != '\0') {
    if (*c == ' ') {
      // continue
    }
    if (isdigit(*c)) {
      current += *c;
    }
    else if (*c == ',') {
      if (current == "") {
        cerr <<
            "Could not parse comma-separated list of integers for argument -k: '" <<
            iters_arg << "'.\n";
        exit(1);
      }
      // conversion is safe because we are in
      res.push_back(stoi(current));
      current = "";
    }
    else {
      cerr <<
          "Could not parse comma-separated list of integers for argument -k: '" <<
          iters_arg << "'.\n";
      exit(1);
    }

    c++;
  }

  if (current != "") {
    res.push_back(stoi(current));
  }

  return res;
}


const int ARG_PARSE_ERROR = 1;
const int ARG_PARSE_QUIT = 0;
const int ARG_PARSE_OK = -1;
// returns ARG_PARSE_OK if execution should continue
// ARG_PARSE_END to end with exit code 0,
// ARG_PARSE_ERROR to end with exit code 1,
int process_args(
    const int argc, char* argv[],
    MCell::SharedGenData& opts
) {
  opts.reset();

  assert(argc > 0);
  while (1) {

    // get the next argument
    int c = getopt_long_only(argc, argv, "hvgtk:cbo:", long_options, nullptr);
    if (c == -1)
      break;

    switch (c) {
      case 'h':
        print_usage(argv[0]);
        return ARG_PARSE_QUIT;
      case 'v':
        print_version(argv[0]);
        return ARG_PARSE_QUIT;
      case 'g':
        opts.debug_mode = true;
        break;
      case 't':
        opts.testing_mode = true;
        break;
      case 'p':
        opts.not_overridable_python_params = true;
        break;
      case 'k':
        opts.checkpoint_iterations = parse_checkpoint_iters(optarg);
        opts.testing_mode = true;
        break;
      case 'b':
        opts.bng_mode = true;
        break;
      case 'o':
        opts.output_files_prefix = optarg;
        break;
      default:
        cerr << "Invalid arguments.\n";
        print_usage(argv[0]);
        return ARG_PARSE_ERROR;
    }
  }

  if (optind < argc) {
    if (argc - optind > 1) {
      cerr << "Only one input data model file can be specified.\n";
      return ARG_PARSE_ERROR;
    }
    opts.input_file = argv[optind];
  }
  else {
    cerr << "Input data model file file was not specified.\n";
    return ARG_PARSE_ERROR;
  }

  if (opts.output_files_prefix == "Untitled") {
    cout << "Ignoring files prefix 'Untitled'.\n";
    opts.output_files_prefix = "";
  }

  return ARG_PARSE_OK;
}


int main(const int argc, char* argv[]) {
  MCell::SharedGenData opts;

  int arg_process_res = process_args(argc, argv, opts);
  if (arg_process_res != ARG_PARSE_OK) {
    return arg_process_res;
  }

  MCell::MCell4Generator converter;
  bool ok = converter.generate(opts);

  if (!ok) {
    cerr << "There was an error while converting " << opts.input_file << " to pymcell code.\n";
    return 1;
  }
  return 0;
}
