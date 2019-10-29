// diverse debug macros

// for mcell3 - crossing memory partitions causes reordering in diffusion and it is not possible to
// compare results anymore

#define MCELL3_IDENTICAL
//MCell4 check

#ifndef MCELL3_IDENTICAL
// enable several things that make comparison with mcell4 easier
#define MCELL3_ONLY_ONE_MEMPART
#define MCELL3_SORTED_VIZ_OUTPUT
#define ASSERT_FOR_MCELL4(...) assert(__VA_ARGS__)

#else

#define ASSERT_FOR_MCELL4(...) do { } while(0)

#endif

//#define DUMP_ALWAYS
//#define DUMP_NEVER

//#define DEBUG_WALL_COLLISIONS

#define DEBUG_DYNAMIC_GEOMETRY
//#define DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS

#if (!defined(NDEBUG) || defined(DUMP_ALWAYS)) && !defined(DUMP_NEVER)

//#define DEBUG_SCHEDULER

#define DEBUG_DEFRAGMENTATION

// cannot be conditioned by iterations
#define DEBUG_RNG_CALLS

// does not generate the same dump as mcell3
//#define DEBUG_SUBPARTITIONS

#define DEBUG_DIFFUSION
#define DEBUG_COLLISIONS
#define DEBUG_REACTIONS

//#define DEBUG_GRIDS

#define FROM_ITERATION 0//250

#define DUMP_CONDITION3(code) do { if ((int)world->current_iterations >= (int)FROM_ITERATION) { code; } } while (0)
#define DUMP_CONDITION4(code) do { if ((int)world->current_iteration >= (int)FROM_ITERATION) { code; } } while (0)

#ifdef DEBUG_SCHEDULER
//#define DUMP_LOCAL_SCHEDULE_HELPER
#endif

#endif
