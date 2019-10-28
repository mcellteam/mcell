/******************************************************************************
 *
 * Copyright (C) 2019 by
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

#ifndef __DUMP_STATE_H__
#define __DUMP_STATE_H__

#include "mcell_structs.h"
#include "grid_util.h"

#include "rng.h"

#define DUMP_EVERYTHING 0xFFFFFFFF

void dump_volume(struct volume* s, const char* comment, unsigned int selected_details);

// the functions below are used to generate log that can be diffed with mcell4's log
void dump_collisions(struct collision* shead);

void dump_surface_molecule(
    struct surface_molecule* amp,
    const char* ind,
    bool for_diff,
    const char* extra_comment,
    unsigned long long iteration,
    double time,
    bool print_position
);

void dump_volume_molecule(
    struct volume_molecule* amp,
    const char* ind,
    bool for_diff,
    const char* extra_comment,
    unsigned long long iteration,
    double time,
    bool print_position
);

void dump_object_list(geom_object* obj, const char* name, const char* comment, const char* ind);

void dump_vector2(struct vector2 vec, const char* extra_comment);

void dump_vector3(struct vector3 vec, const char* extra_comment);

void dump_rng_call_info(struct isaac64_state* rng, const char* extra_comment);

void dump_tile_neighbors_list(struct tile_neighbor *tile_nbr_head, const char* extra_comment, const char* ind);

void dump_wall(wall* w, const char* ind);

void dump_schedule_helper(
    struct schedule_helper* shp,
    const char* name,
    const char* comment,
    const char* ind,
    bool simplified_for_vm);

#endif
