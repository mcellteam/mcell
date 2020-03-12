
#include <cstdlib>
#include <cerrno>
#include <string>


#include "parser_utils.h"
#include "bngl_parser.hpp"

using namespace std;

extern int bngllineno;

namespace BNG {

ostream& errs() {
  cerr << "Line " << bngllineno << ": ";
  return cerr;
}

double convert_to_dbl(const char* str) {
  char* end;
  double res;

  errno = 0; // note: errno is thread-local
  res = strtod(str, &end);
  if (errno != 0 || *end != '\0') {
    errs() << "Could not convert floating point value '" << str << "'.\n";
  }

  return res;
}

long long convert_dec_to_llong(const char* str) {
  char* end;
  long long int res;

  errno = 0; // note: errno is thread-local
  res = strtoll(str, &end, 10);
  if (errno != 0 || *end != '\0') {
    errs() << "Could not convert integer value '" << str << "'.\n";
  }

  return res;
}

} // namespace BNG
