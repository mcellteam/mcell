#ifndef __DUMP_STATE_H__
#define __DUMP_STATE_H__

#include "mcell_structs.h"


#define DUMP_EVERYTHING 0xFFFFFFFF

#ifdef __cplusplus
extern "C"
#endif
void dump_volume(struct volume* s, const char* comment, unsigned int selected_details);


#endif
