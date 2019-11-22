// diverse debug macros

// for mcell3 - crossing memory partitions causes reordering in diffusion and it is not possible to
// compare results anymore

#include <iostream>

#include "dump_state.h"

#define MCELL3_IDENTICAL
//MCell4 check

#ifndef MCELL3_IDENTICAL
// enable several things that make comparison with mcell4 easier
//#define MCELL3_ONLY_ONE_MEMPART
//#define MCELL3_SORTED_VIZ_OUTPUT
//#define MCELL3_SORTED_WALLS_FOR_COLLISION
#define ASSERT_FOR_MCELL4(...) assert(__VA_ARGS__)
//#define ASSERT_FOR_MCELL4(...) do { } while(0)

#else

#define ASSERT_FOR_MCELL4(...) do { } while(0)

#endif

//#define DUMP_ALWAYS
#define DUMP_NEVER


#if (!defined(NDEBUG) || defined(DUMP_ALWAYS)) && !defined(DUMP_NEVER)

//#define DEBUG_WALL_COLLISIONS

#define DEBUG_DYNAMIC_GEOMETRY
#define DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
#define DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS

#define DEBUG_CLOSEST_INTERIOR_POINT

//#define DEBUG_EDGE_INITIALIZATION

//#define DEBUG_SCHEDULER

//#define DEBUG_DEFRAGMENTATION

// cannot be conditioned by iterations
#define DEBUG_RNG_CALLS

// does not generate the same dump as mcell3
//#define DEBUG_SUBPARTITIONS

#define DEBUG_DIFFUSION
#define DEBUG_COLLISIONS
//#define DEBUG_COLLISIONS_WALL_EXTRA
#define DEBUG_REACTIONS

//#define DEBUG_GRIDS

#define FROM_ITERATION 6//250

#define TO_ITERATION 9

#define DUMP_CONDITION3(code) do { if ((int)world->current_iterations >= (int)FROM_ITERATION && (int)world->current_iterations <= (int)TO_ITERATION) { code; } } while (0)
#define DUMP_CONDITION4(code) do { if ((int)world->get_current_iteration() >= (int)FROM_ITERATION && (int)world->get_current_iteration() <= (int)TO_ITERATION) { code; } } while (0)

#ifdef DEBUG_SCHEDULER
//#define DUMP_LOCAL_SCHEDULE_HELPER
#endif

#endif
