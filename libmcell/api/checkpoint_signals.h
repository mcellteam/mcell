/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef LIBMCELL_API_CHECKPOINT_SIGNALS_H_
#define LIBMCELL_API_CHECKPOINT_SIGNALS_H_

#include "defines.h"
#include <string>

namespace MCell {
namespace API {

class Model;

const int SIGNO_NOT_SIGNALED = -1;

// used as a template argument for CustomFunctionCallEvent
struct CheckpointSaveEventContext {
  Model* model;
  std::string dir_prefix;
  bool append_it_to_dir;
};

// may be called only once
void set_checkpoint_signals(Model* model);

// may be safely called multiple times
void unset_checkpoint_signals(Model* model);

// called from CustomFunctionCallEvent
void save_checkpoint_func(const double time, CheckpointSaveEventContext ctx);

} // namespace API
} // namespace MCell

#endif // LIBMCELL_API_CHECKPOINT_SIGNALS_H_
