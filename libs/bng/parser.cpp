
#include <iostream>
#include <cstdlib>
#include <stdio.h>

#include "parser.h"

#include "bngl_scanner.hpp"
#include "bngl_parser.hpp"
#include "bng/species.h"
#include "bng/semantic_analyzer.h"
#include "bng/bng_data.h"

using namespace std;

extern FILE *bnglin;

namespace BNG {

int parse_bngl_file(
    const std::string& file_name,
    BNGData& bng_data,
    const std::map<std::string, float_t>& parameter_overrides) {

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
    BNG::delete_parser_context();
    return 1;
  }

  BNG::SemanticAnalyzer sema;
  sema.check_and_convert_parsed_file(ctx, &bng_data, parameter_overrides);

  ctx->print_error_report();

  int errors = ctx->get_error_count();

  BNG::delete_parser_context();
  fflush(bnglin); // valgrind reports memory leak without this flush call
  fclose(bnglin);
  bngllex_destroy();

  return errors;
}

// bng_data are not cleared and one can continue with adding
// complexes gradually
int parse_single_cplx_string(
    const std::string& cplx_string, BNGData& bng_data,
    Cplx& res_cplx
) {

  create_parser_context();
  ParserContext* ctx = get_parser_context();


  // form input for parser, the !CPLX switches it to a mode where it parses a single string
  char* input = new char[cplx_string.size() + 32];
  strcpy(input, "!CPLX ");
  strcat(input, cplx_string.c_str());

  bngl_scan_string(input);

  ctx->set_current_file_name(cplx_string.c_str());

  // run bison parser
  int res = bnglparse();

  delete [] input;

  if (res != 0) {
    // parse error, do not continue
    BNG::delete_parser_context();
    return 1;
  }

  BNG::SemanticAnalyzer sema;

  CplxVector cplx_vec;

  sema.check_and_convert_single_cplx(ctx, &bng_data, res_cplx);

  ctx->print_error_report();

  int errors = ctx->get_error_count();

  BNG::delete_parser_context();
  bngllex_destroy();

  return errors;
}

} /* namespace BNG */
