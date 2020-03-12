#ifndef LIBS_BNG_PARSER_UTILS_H_
#define LIBS_BNG_PARSER_UTILS_H_

#include <ostream>

namespace BNG {

std::ostream& errs();
double convert_to_dbl(const char* str);
long long convert_dec_to_llong(const char* str);

}

#endif // LIBS_BNG_PARSER_UTILS_H_
