#include "defines.h"

namespace mcell {

std::ostream & operator<<(std::ostream &out, const vec3_t &a) {
    out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
    return out;
}

} // namespace mcell
