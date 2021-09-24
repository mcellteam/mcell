/******************************************************************************
 *
 * Copyright (C) 2019,2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/


// diverse debug macros

#ifndef __DEBUG_CONFIG_H__
#define __DEBUG_CONFIG_H__

// TODO: make the dumping system integrated and build it all for debug builds,
//       enabled by cmdline arguments

// for mcell3 - crossing memory partitions causes reordering in diffusion and it is not possible to
// compare results anymore

#include <iostream>

// use gperftools to print a snapshot of heap, can be displayed with
// pprof mcell.so mem<N>.dump -gv
// MCell must be built with -DENABLE_GPERFTOOLS=ON otherwise
// linking fails with undefined symbol MallocExtension::instance()
#ifdef WITHGPERFTOOLS
//#define PROFILE_MEMORY
#endif

// when enabled, mcell3 produces identical result to the mcell master branch
#define MCELL3_IDENTICAL
#define MCELL4_IDENTICAL

//#define MCELL3_SORTED_MOLS_ON_RUN_TIMESTEP
//#define MCELL3_SORTED_VIZ_OUTPUT

// define when comparin mcell4 and pymcell4 outputs
//#define PYMCELL4_TESTING

// bug in mcellr - exact disk
// most probably mcell3r ever had support for reactive surfaces,
// but keeping this as a macro to maintain compatibility
// if needed
#define FIX_EXTERNAL_SPECIES_WO_RXS_IN_EXACT_DISK

// number placements are in reverse order, but density placements are in correct order,
// enabling this macro unifies it
#define MCELL3_REVERSE_INITIAL_SURF_MOL_PLACEMENT_BY_NUM


// in MCell3 the newly created particles that have long time steps gradually increase
// their timestep to the full value
// we cannot have individual timesteps for each molecule in MCell4 because there would be no way to
// parallelize it
// enabling this macro disables logic in safe_diffusion_step and set_inertness_and_maxtime
//#define MCELL3_MOLECULE_MOVES_WITH_MAXIMUM_TIMESTEP

//#define MCELL3_NO_COMPARMTENT_IN_PRINTOUTS

#ifndef MCELL3_IDENTICAL

// ---- MCell4 macros to match MCell3R ----
#ifndef MCELL4_IDENTICAL
// MCell3R always creates new products when there is the same species
// on the reactants and products side (unlike in normal MCell)
// this switch is for validation of MCell4 BNG against MCell3R
// One can detect this automatically or with an option, but both seemed as ugly solution
// testsuite for mcell4 won't pass
#define MCELL4_DO_NOT_REUSE_REACTANT

// MCell3R seems to sort reaction products by name
//#define MCELL4_SORT_RXN_PRODUCTS_BY_NAME
//#define MCELL4_SORT_RXN_PRODUCTS_BY_NAME_REV
#define MCELL4_SORT_RXN_PRODUCTS_BY_LENGTH_DESC

// tentative sorting of reactions in a rxn class to match MCell3R
#define MCELL4_REVERSED_RXNS_IN_RXN_CLASS

// do not call random generator when probability of an unimol rxn is 0
// required for compatibility with MCell3R that ignores such reactions, 
// however MCell3 call rng even if the prob is 0
// testsuite for mcell4 won't pass
#define MCELL4_NO_RNG_FOR_UNIMOL_RXN_P_0

// when a volume product is created in a surface reaction, it is bumped from the surface,
// when is this macro enabled, it is bumped only by EPS same as in MCell3, not 16*EPS
#define MCELL4_VOL_PROD_BUMP_FROM_SURFACE_ONLY_ONE_EPS

// sort molecules in schedule helper according to ID before a new timestep begins
// testsuite for mcell4 won't pass
#define MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID

#define MCELL3_UNIMOL_RX_ABSORB_NO_RNG

// do not reuse molecule IDs
#define MCELL3_DO_NOT_REUSE_MOL_ID_UNIMOL_RXN
#endif

// ^^^^ MCell4 macros to match MCell3R ^^^^

// enable several things that make comparison with mcell4 easier
#define MCELL3_ONLY_ONE_MEMPART
#define MCELL3_SORTED_VIZ_OUTPUT
#define MCELL3_SORTED_WALLS_FOR_COLLISION

#ifndef MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID
// sort molecules when run_timestep is started, replaced by better MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID that
// however changes ordering for mcell4
#define MCELL3_SORTED_MOLS_ON_RUN_TIMESTEP
#endif

//#define MCELL3_4_SAFE_DIFF_STEP_RETURNS_CONSTANT

//#define MCELL3_NEXT_BARRIER_IS_THE_NEXT_TIMESTEP // do not diffuse more than until the end of the timestep

// messes up original ordering when all releases are planned for 0, not sure why yet
// use only really when needed
//#define MCELL3_RELEASE_ACCORDING_TO_EVENT_TIME

//#define MCELL3_ROUND_TSTEPS

//#define MCELL3_ALWAYS_DIFFUSE // non-diffusable molecules are scheduled differently when there is a unimol reaction and ??
#define ASSERT_FOR_MCELL4(...) assert(__VA_ARGS__)

#else

#define ASSERT_FOR_MCELL4(...) do { } while(0)

#endif


#define DUMP4_PRECISION_DEFAULT 5

#ifdef PYMCELL4_TESTING
// needed for easier diff between mcell4 and pymcell4
#define SORT_MCELL4_SPECIES_BY_NAME

// testsuite for mcell4 won't pass
// #define ORDER_RXNS_IN_RXN_CLASS_BY_NAME

#define DUMP4_PRECISION 17
#else
#define DUMP4_PRECISION DUMP4_PRECISION_DEFAULT
#endif

//#define DEBUG_EXTRA_CHECKS



//#define DUMP_ALWAYS
#define DUMP_NEVER

#if (!defined(NDEBUG) || defined(DUMP_ALWAYS)) && !defined(DUMP_NEVER)

/*
#define MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID
#define MCELL3_SORTED_VIZ_OUTPUT
#define MCELL3_ONLY_ONE_MEMPART
#define MCELL3_4_SAFE_DIFF_STEP_RETURNS_CONSTANT

#define MCELL3_SORTED_WALLS_FOR_COLLISION
#define MCELL3_RELEASE_ACCORDING_TO_EVENT_TIME
#define MCELL3_UNIMOL_RX_ABSORB_NO_RNG
*/
//#define DUMP_LOCAL_SCHEDULE_HELPER

#define TRACK_MOL 1 // true
#define TRACKED_MOL_ID 61254

#define FROM_ITERATION 100000
#define TO_ITERATION 1000000

#define DUMP_NONDIFFUSING_VMS

#if 1
#define DEBUG_DIFFUSION
#define DEBUG_COLLISIONS
//#define NODEBUG_WALL_COLLISIONS
#endif

#define DEBUG_RXNS

//#define DEBUG_COMPARTMENTS
//#define DEBUG_RNG_CALLS // cannot be conditioned by iterations

//#define DEBUG_WALL_COLLISIONS

//#define DEBUG_DYNAMIC_GEOMETRY
//#define DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
//#define DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS

//#define DEBUG_CLOSEST_INTERIOR_POINT


//#define DEBUG_POLY_EDGE_INITIALIZATION
//#define DEBUG_EDGE_INITIALIZATION

//#define DEBUG_SCHEDULER_ACTION
//#define DEBUG_SCHEDULER

//#define DEBUG_DEFRAGMENTATION

//#define DEBUG_EXACT_DISK


//#define DEBUG_RELEASES // cannot be conditioned by iterations


// does not generate the same dump as mcell3
//#define DEBUG_SUBPARTITIONS

//#define DEBUG_COUNTED_VOLUMES
//#define DEBUG_TRANSPARENT_SURFACES

//#define DEBUG_TIMING

//#define DEBUG_DIFFUSION_EXTRA
//#define DEBUG_COLLISIONS_WALL_EXTRA

//#define DEBUG_COUNTED_VOLUMES

//#define DEBUG_REACTION_PROBABILITIES  // cannot be conditioned by iterations

//#define DEBUG_GRIDS


#define DUMP_CONDITION3(code) do { if ((int)world->current_iterations >= (int)FROM_ITERATION && (int)world->current_iterations <= (int)TO_ITERATION) { code; } } while (0)
#define DUMP_CONDITION4(id, code) do { \
  if (\
      ((int)world->get_current_iteration() >= (int)FROM_ITERATION && (int)world->get_current_iteration() <= (int)TO_ITERATION) || \
      (TRACK_MOL && id == TRACKED_MOL_ID)) { \
    code; } } while (0)

#define DUMP_CONDITION4P(code) do { if ((int)p.stats.get_current_iteration() >= (int)FROM_ITERATION && (int)p.stats.get_current_iteration() <= (int)TO_ITERATION) { code; } } while (0)

#ifdef DEBUG_SCHEDULER
#define DUMP_LOCAL_SCHEDULE_HELPER
#endif

#endif

#endif // __DEBUG_CONFIG_H__
