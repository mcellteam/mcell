/*
 * callback_data.h
 *
 *  Created on: Oct 14, 2019
 *      Author: ahusar
 */

#ifndef SRC4_CALLBACK_INFO_H_
#define SRC4_CALLBACK_INFO_H_

#ifndef SWIG
#include "callback_structs.h"
#endif

typedef void (*wall_hit_callback_func)(const MCell::WallHitInfo& info, void*);

#ifndef NOSWIG
// to be used only from SWIG .i file
void py_callback_wall_hit(const MCell::WallHitInfo& info, void*);
#endif

#endif /* SRC4_CALLBACK_INFO_H_ */
