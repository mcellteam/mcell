/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef UTILS_DATA_MODEL_GEOMETRY_H_
#define UTILS_DATA_MODEL_GEOMETRY_H_

#include "json/json.h"

namespace MCell {

// throws ConversionError if computation was not successful (may possibly throw other exceptions)
void compute_volume_and_area(Json::Value& model_object, double& volume, double& area);

} // namespace MCell

#endif /* UTILS_DATA_MODEL_GEOMETRY_H_ */
