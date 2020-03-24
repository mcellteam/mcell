
#ifndef API_MODEL_H
#define API_MODEL_H

#include <string>

#include "../generated/gen_model.h"
#include "common.h"
#include "subsystem.h"
#include "instantiation_data.h"

namespace MCell {
namespace API {

class Model: public GenModel, public Subsystem, public InstantiationData {
public:

  void run_iterations(const long iterations) override {}
  void add_subsystem(const Subsystem* subsystem) override {}
  void add_instantiation_data(const InstantiationData* instantiation_data) override {}


};


} // namespace API
} // namespace MCell

#endif // API_MODEL_H
