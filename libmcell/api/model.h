
#ifndef API_MODEL_H
#define API_MODEL_H

#include <string>

#include "../generated/gen_model.h"
#include "common.h"
#include "subsystem.h"
#include "instantiation_data.h"

namespace MCell {
namespace API {

// multiple inheritance...
class Model: public GenModel, public Subsystem, public InstantiationData {
public:

};


} // namespace API
} // namespace MCell

#endif // API_MODEL_H
