/**
 * This is an attempt to define a C++ class that holds information on model.
 * Such classes are then initialized from Python code.
 *
 * Semantic check code is contained
 *
 * Separate file is generated for each class
 */
#ifndef API_GEN_RELEASE_SITE_H
#define API_GEN_RELEASE_SITE_H

#include "../api/mcell.h"

namespace MCell {
namespace API {

#define RELEASE_SITE_CTOR() \
    ReleaseSite( \
        const std::string& name_, \
        const std::string& shape_, \
        const Species* molecule_, \
        const Vec3& location_ = VEC3_UNSET, \
        const float_t site_diameter_ = FLT_UNSET, \
        const float_t site_radius_ = FLT_UNSET, \
        const float_t release_probability_ = FLT_UNSET \
    ) { \
      class_name = "ReleaseSite"; \
      name = name_; \
      shape = shape_; \
      molecule = molecule_; \
      location = location_; \
      site_diameter = site_diameter_; \
      site_radius = site_radius_; \
      release_probability = release_probability_; \
    }

// do not instantiate this class directly, use its superclass ReleaseSite
class GenReleaseSite: public BaseDataClass {
public:
  SemRes check_semantics(std::ostream& out) const override;
  std::string to_str() const override;

  std::string shape;
  // get/set

  Vec3 location;
  const Species* molecule;
  // get/set

  float_t site_diameter;
  virtual void set_site_diameter(float_t v) {
    site_diameter = v;
  }
  virtual float_t get_site_diameter() const {
    return site_diameter;
  }

  float_t site_radius;
  // get/set

  float_t release_probability;
  // get/set

  /*
  uint number_to_release;
  float_t concentration;
  float_t density;

  // Gaussian release number:
  float_t mean_number;
  float_t standard_deviation;

  // release pattern
  std::string release_pattern_name;
  float_t release_pattern_delay;
  float_t release_pattern_release_interval;
  float_t release_pattern_train_duration;
  float_t release_pattern_train_interval;
  */
};

void define_binding_ReleaseSite(py::module& m);

} // namespace API
} // namespace MCell

#endif // API_GEN_RELEASE_SITE_H
