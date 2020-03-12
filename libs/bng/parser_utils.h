

// defined in bngl_parser.y
void bnglerror(const char *s);
void bnglerror_fmt(const char *s, ...);

namespace BNGL {

double convert_to_dbl(const char* str);
long long convert_dec_to_llong(const char* str);

}
