// diverse debug macros

// for mcell3 - crossing memory partitions causes reordering in diffusion and it is not possible to
// copare results anymore
#define MCELL3_ONLY_ONE_MEMPART

//#define DUMP_ALWAYS
//#define DUMP_NEVER

#if (!defined(NDEBUG) || defined(DUMP_ALWAYS)) && !defined(DUMP_NEVER)

//#define DEBUG_SCHEDULER

#ifdef DEBUG_SCHEDULER
//#define DUMP_LOCAL_SCHEDULE_HELPER
#endif

#define DEBUG_DEFRAGMENTATION

// cannot be conditioned by iterations
#define DEBUG_RNG_CALLS

// does not generate the same dump as mcell3
//#define DEBUG_SUBPARTITIONS

#define DEBUG_DIFFUSION
#define DEBUG_COLLISIONS
#define DEBUG_REACTIONS


#define FROM_ITERATION 0

#define DUMP_CONDITION3(code) do { if ((int)world->current_iterations >= (int)FROM_ITERATION) { code; } } while (0)
#define DUMP_CONDITION4(code) do { if ((int)world->current_iteration >= (int)FROM_ITERATION) { code; } } while (0)

#endif
