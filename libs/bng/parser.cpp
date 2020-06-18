
#include <iostream>
#include <cstdlib>
#include <stdio.h>

#include "parser.h"

#include "bngl_parser.hpp"
#include "bng/species.h"
#include "bng/semantic_analyzer.h"
#include "bng/bng_data.h"

using namespace std;

extern FILE *bnglin;

namespace BNG {

int parse_bngl_file(const std::string& file_name, BNGData& bng_data) {
  bng_data.clear();

  create_parser_context();
  ParserContext* ctx = get_parser_context();

  FILE *infile = fopen(file_name.c_str(), "r");
  if (infile == nullptr) {
    cerr << "Could not open input file.\n";
    return 1;
  }

  ctx->set_current_file_name(file_name.c_str());
  // set input file for flex
  bnglin = infile;
  // run bison parser
  int res = bnglparse();
  if (res != 0) {
    // parse error, do not continue
    return 1;
  }

  BNG::SemanticAnalyzer sema;
  sema.check_and_convert(ctx, &bng_data);

  ctx->print_error_report();

  int errors = ctx->get_error_count();

  BNG::delete_parser_context();

  return errors;
}

} /* namespace BNG */
