
#ifndef API_GEOMETRY_OBJECT_H
#define API_GEOMETRY_OBJECT_H

#include <string>

#include "../generated/gen_geometry_object.h"
#include "common.h"

namespace MCell {
namespace API {

class GeometryObject: public GenGeometryObject {
public:
  GEOMETRY_OBJECT_CTOR()
};


} // namespace API
} // namespace MCell

#endif // API_GEOMETRY_OBJECT_H
