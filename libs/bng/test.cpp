#include <iostream>
#include <cstdlib>

#include "bngl_parser.hpp"
#include "semantic_analyzer.h"
#include "bng_engine.h"
#include "cplx_species.h"

using namespace std;

extern FILE *bnglin;

// returns 0 if everything was ok
int parse_bngl(char const *name, const bool dump_ast, const bool dump_bng_data) {

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

  if (dump_ast) {
    ctx->dump();
  }

  BNG::SemanticAnalyzer sema;

  BNG::BNGEngine<BNG::CplxSpecies> bng_engine;
  sema.check_and_convert(ctx, &bng_engine.get_data());

  if (dump_bng_data) {
    bng_engine.get_data().dump();
  }

  ctx->print_error_report();

  int errors = ctx->get_error_count();

  BNG::delete_ast_context();

  return errors != 0;
}

int main(int argc, const char* argv[]) {

  // TODO - use file as input
  if (argc != 2 && argc != 3) {
    cerr << "Expected input file as argument, second optional arg enables AST dump\n";
  }

  bool dump_ast = false;
  bool dump_bng = false;
  if (argc == 3 && string(argv[2]) == "-a") {
    dump_ast = true;
  }
  if (argc == 3 && string(argv[2]) == "-b") {
    dump_bng = true;
  }

  return parse_bngl(argv[1], dump_ast, dump_bng);
}
