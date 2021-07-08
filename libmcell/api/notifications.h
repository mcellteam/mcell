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

#ifndef API_NOTIFICATIONS_H
#define API_NOTIFICATIONS_H

#include "generated/gen_notifications.h"
#include "api/api_common.h"

namespace MCell {
namespace API {

class Notifications: public GenNotifications {
public:
  NOTIFICATIONS_CTOR()
};

} // namespace API
} // namespace MCell

#endif // API_NOTIFICATIONS_H
