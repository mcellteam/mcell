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

#ifndef LIBMCELL_API_BINDINGS_H
#define LIBMCELL_API_BINDINGS_H

namespace MCell {
namespace API {

// used to detect Ctrl-C signal from the user,
// throws an exception if so
void check_ctrl_c(const double current_time, World* arg);

}
}


#endif // LIBMCELL_API_BINDINGS_H
