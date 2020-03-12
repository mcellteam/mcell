#include <iostream>
#include <cstdlib>

#include "bngl_parser.hpp"

using namespace std;

int main(int argc, const char* argv[]) {

  // TODO - use file as input
  /*if (argc != 2) {
    cerr << "Expected input file as argument\n";
  }*/

  int res = bnglparse();



  if (res != 0) {
    cerr << "\nParse error code is " << res << ".\n";
  }
  return res;
}
