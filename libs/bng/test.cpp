#include <iostream>
#include <cstdlib>

#include "bngl_parser.hpp"

using namespace std;


void parse_bngl() {

  BNG::create_ast_context();

  int res = bnglparse();
  if (res != 0) {
    cerr << "\nParse error code is " << res << ".\n";
  }

  BNG::ASTContext* ctx = BNG::get_ast_context();
  ctx->dump();

  ctx->print_error_report();

  BNG::delete_ast_context();
}

int main(int argc, const char* argv[]) {

  // TODO - use file as input
  /*if (argc != 2) {
    cerr << "Expected input file as argument\n";
  }*/

  parse_bngl();

  return 0;
}
