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

// TODO: make the dumping system integrated and build it all for debug builds,
//       enabled by cmdline arguments

// for mcell3 - crossing memory partitions causes reordering in diffusion and it is not possible to
// compare results anymore

#include <iostream>

#include "dump_state.h"

//#define MCELL3_IDENTICAL

#ifndef MCELL3_IDENTICAL
// enable several things that make comparison with mcell4 easier
#define MCELL3_ONLY_ONE_MEMPART
#define MCELL3_SORTED_VIZ_OUTPUT
#define MCELL3_SORTED_WALLS_FOR_COLLISION
#define MCELL3_SORTED_MOLS_ON_RUN_TIMESTEP  // sort molecules in schedule helper according to ID before a new timestep begins
#define ASSERT_FOR_MCELL4(...) assert(__VA_ARGS__)

#else

#define ASSERT_FOR_MCELL4(...) do { } while(0)

#endif

//#define DUMP_ALWAYS
//#define DUMP_NEVER


#if (!defined(NDEBUG) || defined(DUMP_ALWAYS)) && !defined(DUMP_NEVER)

//#define DEBUG_WALL_COLLISIONS

#define DEBUG_DYNAMIC_GEOMETRY
#define DEBUG_DYNAMIC_GEOMETRY_MCELL4_ONLY
//#define DEBUG_DYNAMIC_GEOMETRY_COLLISION_DETECTIONS

//#define DEBUG_CLOSEST_INTERIOR_POINT

//#define DEBUG_EDGE_INITIALIZATION

//#define DEBUG_SCHEDULER

//#define DEBUG_DEFRAGMENTATION


//#define DEBUG_RNG_CALLS // cannot be conditioned by iterations

// does not generate the same dump as mcell3
//#define DEBUG_SUBPARTITIONS

#if 1
#define DEBUG_DIFFUSION
#define DEBUG_COLLISIONS
#define DEBUG_REACTIONS
#endif


//#define DEBUG_TIMING

//#define DEBUG_DIFFUSION_EXTRA
//#define DEBUG_COLLISIONS_WALL_EXTRA

//
//#define DEBUG_REACTION_PROBABILITIES  // cannot be conditioned by iterations

//#define DEBUG_GRIDS

#define FROM_ITERATION 0

#define TO_ITERATION 20

#define DUMP_CONDITION3(code) do { if ((int)world->current_iterations >= (int)FROM_ITERATION && (int)world->current_iterations <= (int)TO_ITERATION) { code; } } while (0)
#define DUMP_CONDITION4(code) do { if ((int)world->get_current_iteration() >= (int)FROM_ITERATION && (int)world->get_current_iteration() <= (int)TO_ITERATION) { code; } } while (0)
#define DUMP_CONDITION4P(code) do { if ((int)p.stats.get_current_iteration() >= (int)FROM_ITERATION && (int)p.stats.get_current_iteration() <= (int)TO_ITERATION) { code; } } while (0)

#ifdef DEBUG_SCHEDULER
#define DUMP_LOCAL_SCHEDULE_HELPER
#endif

#endif

#endif // DEBUG_CONFIG_H
