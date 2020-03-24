
#ifndef API_INSTANTIATION_DATA_H
#define API_INSTANTIATION_DATA_H

#include <string>

#include "../generated/gen_instantiation_data.h"
#include "common.h"

namespace MCell {
namespace API {

class InstantiationData: public GenInstantiationData {
public:

  void add_release_site(const Species* s) override {}
  ReleaseSite* find_release_site(const std::string& name) override {return nullptr;}
  void add_geometry_object(const GeometryObject* o) override {}
  void find_geometry_object(const std::string& name) override {}
};


} // namespace API
} // namespace MCell

#endif // API_INSTANTIATION_DATA_H
