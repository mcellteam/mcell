#ifndef MCELL_COUNT_UTIL
#define MCELL_COUNT_UTIL

#include "mcell_structs.h"

void update_collision_count(struct species *sp,struct region_list *rl,int direction,int crossed);
void count_me_by_region(struct abstract_molecule *me,int n);

#endif
