#include <iostream>
#include <cstdlib>

#include "bngl_parser.hpp"
#include "bng/species.h"
#include "bng/semantic_analyzer.h"
#include "bng/bng_engine.h"

using namespace std;

extern FILE *bnglin;

// returns 0 if everything was ok
int parse_bngl(char const *name, const bool dump_ast, const bool dump_bng_data) {

  BNG::create_parser_context();
  BNG::ParserContext* ctx = BNG::get_parser_context();

  FILE *infile = fopen(name, "r");
  if (infile == nullptr) {
    cerr << "Could not open input file.\n";
    return 1;
  }

  ctx->set_current_file_name(name);
  bnglin = infile;
  int res = bnglparse();
  if (res != 0) {
    cerr << "\nParse error code is " << res << ".\n";
  }

  if (dump_ast) {
    ctx->dump();
  }

  BNG::SemanticAnalyzer sema;

  BNG::BNGConfig bng_config; // we do not care about the config values here
  BNG::BNGEngine bng_engine(bng_config);
  sema.check_and_convert_parsed_file(ctx, &bng_engine.get_data());

  if (dump_bng_data) {
    bng_engine.get_data().dump();
  }

  ctx->print_error_report();

  int errors = ctx->get_error_count();

  BNG::delete_parser_context();

  return errors != 0;
}

extern int bngldebug;

int main(int argc, const char* argv[]) {

  if ((argc == 2 && strcmp(argv[1], "-h") == 0) || (argc != 2 && argc != 3)) {
    cerr << "Expected input file as argument, second optional arg (-a or -b) enables AST or BNG dump, or -g to enable parsing debug\n";
    return 1;
  }

  bool dump_ast = false;
  bool dump_bng = false;
  if (argc == 3) {
    if (string(argv[2]) == "-a") {
      dump_ast = true;
    }
    else if (string(argv[2]) == "-b") {
      dump_bng = true;
    }
    else if (argc == 3 && string(argv[2]) == "-g") {
      bngldebug = 1;
    }
    else {
      cerr << "Invalid argument.\n";
      return 1;
    }
  }

  return parse_bngl(argv[1], dump_ast, dump_bng);
}
