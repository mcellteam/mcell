
#ifndef API_RELEASE_SITE_H
#define API_RELEASE_SITE_H

#include <string>

#include "../examples/gen_release_site.h"
#include "common.h"

namespace MCell {
namespace API {

class ReleaseSite: public GenReleaseSite {
public:
  RELEASE_SITE_CTOR()

  // actual manual implementation of a semantic check
  SemRes check_semantics(std::ostream& out) const override {
    SemRes base_res = GenReleaseSite::check_semantics(out);
    if (base_res != SemRes::OK) {
      return base_res;
    }

    if (is_set(site_diameter) && is_set(site_radius)) {
      out << "Only either 'site_diameter' or 'site_radius' can be set.\n";
      return SemRes::ERROR;
    }

    return SemRes::OK;
  }
};



} // namespace API
} // namespace MCell

#endif // API_RELEASE_SITE_H
