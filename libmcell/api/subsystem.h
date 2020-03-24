
#ifndef API_SUBSYSTEM_H
#define API_SUBSYSTEM_H

#include <string>

#include "../generated/gen_subsystem.h"
#include "common.h"

namespace MCell {
namespace API {

class Subsystem: public GenSubsystem {
public:
  void add_species(const Species* s) override {}
  Species* find_species(const std::string& name) override {return nullptr;}
};


} // namespace API
} // namespace MCell

#endif // API_SUBSYSTEM_H
