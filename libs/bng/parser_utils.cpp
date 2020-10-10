
#include <cstdlib>
#include <cerrno>
#include <string>


#include "bng/parser_utils.h"
#include "bngl_parser.hpp"

using namespace std;

extern int bngllineno;
extern BNG::ParserContext* g_ctx;

namespace BNG {

ostream& errs_loc() {
  assert(g_ctx != nullptr);
  cerr <<
      g_ctx->get_current_file_name() << ":" << bngllineno <<
      ": error: ";
  return cerr;
}


ostream& errs_loc(const ASTBaseNode* loc) {
  assert(loc != nullptr);
  assert(loc->has_loc);
  assert(loc->file != nullptr);
  cerr <<
      loc->file << ":" << loc->line <<
      ": error: ";
  return cerr;
}


double convert_to_dbl(const char* str) {
  char* end;
  double res;

  errno = 0; // note: errno is thread-local
  res = strtod(str, &end);
  if (errno != 0 || *end != '\0') {
    errs_loc() << "Could not convert floating point value '" << str << "'.\n";
  }

  return res;
}


long long convert_dec_to_llong(const char* str) {
  char* end;
  long long int res;

  errno = 0; // note: errno is thread-local
  res = strtoll(str, &end, 10);
  if (errno != 0 || *end != '\0') {
    errs_loc() << "Could not convert integer value '" << str << "'.\n";
  }

  return res;
}


char *strdup_new(const char *src)
{
  char *str;
  size_t size = strlen(src) + 1;

  str = new char[size];
  if (str != nullptr) {
    memcpy(str, src, size);
  }
  return str;
}

} // namespace BNG
