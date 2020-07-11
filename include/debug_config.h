/******************************************************************************
 *
 * Copyright (C) 2019,2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/


// diverse debug macros

#ifndef DEBUG_CONFIG_H
#define DEBUG_CONFIG_H

// TODO: do not include in defines.h for faster compilation

// TODO: make the dumping system integrated and build it all for debug builds,
//       enabled by cmdline arguments

// for mcell3 - crossing memory partitions causes reordering in diffusion and it is not possible to
// compare results anymore

#include <iostream>

#include "dump_state.h"

// when enabled, mcell3 produces identical result to the mcell master branch
//#define MCELL3_IDENTICAL

// define when comparin mcell4 and pymcell4 outputs
//#define PYMCELL4_TESTING

// bug in mcellr - exact disk
// most probably mcell3r ever had support for reactive surfaces,
// but keeping this as a macro to maintain compatibility
// if needed
#define FIX_EXTERNAL_SPECIES_WO_RXS_IN_EXACT_DISK

// number placement are in reverse order, but densitty placement are in correct order, 
// enabling this macro unifies it
#define MCELL3_REVERSE_INITIAL_SURF_MOL_PLACEMENT_BY_NUM

#ifndef MCELL3_IDENTICAL

// sort molecules in schedule helper according to ID before a new timestep begins
// testsuite for mcell4 won't pass
//#define MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID

// enable several things that make comparison with mcell4 easier
#define MCELL3_ONLY_ONE_MEMPART
#define MCELL3_SORTED_VIZ_OUTPUT
#define MCELL3_SORTED_WALLS_FOR_COLLISION

#ifndef MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID
// sort molecules when run_timestep is started, replaced by better MCELL3_4_ALWAYS_SORT_MOLS_BY_TIME_AND_ID that
// however changes ordering for mcell4
#define MCELL3_SORTED_MOLS_ON_RUN_TIMESTEP
#endif

#define MCELL3_NEXT_BARRIER_IS_THE_NEXT_TIMESTEP // do not diffuse more than until the end of the timestep

#define MCELL3_RELEASE_ACCORDING_TO_EVENT_TIME

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

//#define COLLECT_SUBPARTS_LEGACY


//#define DUMP_ALWAYS
//#define DUMP_NEVER

#if (!defined(NDEBUG) || defined(DUMP_ALWAYS)) && !defined(DUMP_NEVER)

#define FROM_ITERATION 0
#define TO_ITERATION 7

#define DUMP_NONDIFFUSING_VMS

#if 1
#define DEBUG_DIFFUSION
#define DEBUG_COLLISIONS
#define NODEBUG_WALL_COLLISIONS
#endif

#define DEBUG_RXNS

//#define DEBUG_RNG_CALLS // cannot be conditioned by iterations

//#define DEBUG_CPLX_RXNS
//#define DEBUG_CPLX_MATCHING
//#define DEBUG_CPLX_MATCHING_EXTRA_COMPARE

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
#define DUMP_CONDITION4(code) do { if ((int)world->get_current_iteration() >= (int)FROM_ITERATION && (int)world->get_current_iteration() <= (int)TO_ITERATION) { code; } } while (0)
#define DUMP_CONDITION4P(code) do { if ((int)p.stats.get_current_iteration() >= (int)FROM_ITERATION && (int)p.stats.get_current_iteration() <= (int)TO_ITERATION) { code; } } while (0)

#ifdef DEBUG_SCHEDULER
#define DUMP_LOCAL_SCHEDULE_HELPER
#endif

#endif

#endif // DEBUG_CONFIG_H
