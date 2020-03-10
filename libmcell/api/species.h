
#ifndef API_SPECIES_H
#define API_SPECIES_H

#include <string>

#include "../generated/gen_species.h"
#include "common.h"

namespace MCell {
namespace API {

class Species: public GenSpecies {
public:
  SPECIES_CTOR()

  // actual manual implementation of a semantic check
  SemRes check_semantics(std::ostream& out) const override {
    SemRes base_res = GenSpecies::check_semantics(out);
    if (base_res != SemRes::OK) {
      return base_res;
    }

    if (is_set(diffusion_constant_2d) && is_set(diffusion_constant_3d)) {
      out << get_object_name() << "Only either 'diffusion_constant_2d' or 'diffusion_constant_3d' can be set.";
      return SemRes::ERROR;
    }

    return SemRes::OK;
  }
};

} // namespace API
} // namespace MCell

#endif // API_SPECIES_H
