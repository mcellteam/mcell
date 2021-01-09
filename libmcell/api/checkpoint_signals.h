/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#ifndef LIBMCELL_API_CHECKPOINT_SIGNALS_H_
#define LIBMCELL_API_CHECKPOINT_SIGNALS_H_

#include <string>

namespace MCell {
namespace API {

class Model;

// may be called only once
void set_checkpoint_signals(Model* model);

// may be safely called multiple times
void unset_checkpoint_signals(Model* model);

// used as a template argument for CustomFunctionCallEvent
struct CheckpointSaveEventContext {
  Model* model;
  std::string dir_prefix;
  bool append_it_to_dir;
};

} // namespace API
} // namespace MCell

#endif // LIBMCELL_API_CHECKPOINT_SIGNALS_H_
