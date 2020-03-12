
#include <cstdlib>
#include <cerrno>
#include <string>

#include "parser_utils.h"

using namespace std;

namespace BNGL {

double convert_to_dbl(const char* str) {
  char* end;
  double res;

  errno = 0; // note: errno is thread-local
  res = strtod(str, &end);
  if (errno != 0 || *end != '\0') {
    bnglerror_fmt("Could not convert floating point value '%s'.", str);
  }

  return res;
}

long long convert_dec_to_llong(const char* str) {
  char* end;
  long long int res;

  errno = 0; // note: errno is thread-local
  res = strtoll(str, &end, 10);
  if (errno != 0 || *end != '\0') {
    bnglerror_fmt("Could not convert integer value '%s'.", str);
  }

  return res;
}

} // namespace BNGL
