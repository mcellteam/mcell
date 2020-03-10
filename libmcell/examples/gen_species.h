/**
 * This is an attempt to define a C++ class that holds information on model.
 * Such classes are then initialized from Python code.
 *
 * Semantic check code is contained
 *
 * Separate file is generated for each class
 */
#ifndef API_GEN_SPECIES_H
#define API_GEN_SPECIES_H

#include <string>
#include "../api/mcell.h"

namespace MCell {
namespace API {

// this macro must be instantiated in the derived class
#define SPECIES_CTOR() \
    Species( \
        const std::string& name_, \
        const float_t diffusion_constant_2d_ = FLT_UNSET, \
        const float_t diffusion_constant_3d_ = FLT_UNSET \
    ) { \
      class_name = "Species"; \
      name = name_; \
      diffusion_constant_2d = diffusion_constant_2d_; \
      diffusion_constant_3d = diffusion_constant_3d_; \
    }

// do not instantiate this class directly, use its superclass Species
class GenSpecies: public BaseDataClass {
public:
  SemRes check_semantics(std::ostream& out) const override;
  std::string to_str() const override;

  // --- attributes ---
  float_t diffusion_constant_2d;
  // get/set

  float_t diffusion_constant_3d;
  // get/set
};

void define_binding_Species(py::module& m);

} // namespace API
} // namespace MCell

#endif // API_GEN_SPECIES_H
