#include <iostream>
#include <cstdlib>

#include "bngl_parser.hpp"
#include "semantic_analyzer.h"
#include "bng_engine.h"

using namespace std;

extern FILE *bnglin;

int parse_bngl(char const *name, const bool dump) {

  BNG::create_ast_context();
  BNG::ASTContext* ctx = BNG::get_ast_context();

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

  if (dump) {
    ctx->dump();
  }

  BNG::SemanticAnalyzer sema;

  BNG::BNGData bng_data;
  sema.check_and_convert(ctx, &bng_data);

  ctx->print_error_report();

  int errors = ctx->get_error_count();

  BNG::delete_ast_context();

  return errors == 0;
}

int main(int argc, const char* argv[]) {

  // TODO - use file as input
  if (argc != 2 && argc != 3) {
    cerr << "Expected input file as argument, second optional arg enables AST dump\n";
  }

  bool dump = false;
  if (argc == 3 && string(argv[2]) == "-d") {
    dump = true;
  }

  return parse_bngl(argv[1], dump);
}
