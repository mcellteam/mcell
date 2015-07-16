/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
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

#ifndef OUTBOUND_MOLECULES_H
#define OUTBOUND_MOLECULES_H

#include "config.h"

#include "util.h"
#include "vector.h"

enum {
  MOLECULE_QUEUE_LENGTH = 256
};

/* Data structures to describe release events */
struct release_event_queue {
  struct release_event_queue *next;
  double event_time;                     /* Time of the release */
  struct release_site_obj *release_site; /* What to release, where to release it, etc */
  double t_matrix[4][4];                 /* transformation matrix for location of release site */
  int train_counter;                     /* counts executed trains */
  double train_high_time;                /* time of the train's start */
};

typedef struct transmitted_molecule
{
  struct volume_molecule         *molecule;
  struct subvolume               *target;
  struct vector3                  disp_remainder;
  double                          time_remainder;
} transmitted_molecule_t;

typedef struct transmitted_molecules
{
  struct transmitted_molecules   *next;
  int                             fill;
  transmitted_molecule_t          molecules[MOLECULE_QUEUE_LENGTH];
} transmitted_molecules_t;

typedef struct delayed_release
{
  struct delayed_release         *next;
  struct release_event_queue      event;
  struct magic_list              *incantation;
} delayed_release_t;

typedef struct outbound_molecules
{
  transmitted_molecules_t       *molecule_queue;
  delayed_release_t             *release_queue;
} outbound_molecules_t;

typedef int transmitted_molecule_iter_t;

/* Initialize an outbound molecules queue. */
#define outbound_molecules_init(om) do {                    \
  (om)->molecule_queue = NULL;                              \
  (om)->release_queue  = NULL;                              \
} while (0)

/* Add a volume molecule to the outbound queue. */
void outbound_molecules_add_molecule(outbound_molecules_t *queue,
                                     struct volume_molecule *mol,
                                     struct subvolume      *target,
                                     struct vector3        *disp,
                                     double                 t_remain);

/* Add a reaction-triggered release to the outbound queue. */
void outbound_molecules_add_release(outbound_molecules_t *queue,
                                    struct release_event_queue *event,
                                    struct magic_list *incantation);

/* Play all delayed molecule placements. */
struct volume;
void outbound_molecules_play(struct volume *world,
                             outbound_molecules_t *queue);

/* Begin iterating over all outbound molecules (destructively). */
#define outbound_molecules_begin(q, iter) do {              \
    if ((q)->molecule_queue == NULL)                        \
      *(iter) = -1;                                         \
    else                                                    \
      *(iter) = 0;                                          \
} while (0)

/* Check if we've finished iterating. */
#define outbound_molecules_finished(q, iter) ((*(iter)) < 0)

/* Get the next outbound molecule (destructively). */
transmitted_molecule_t *outbound_molecules_next(outbound_molecules_t *queue,
                                                transmitted_molecule_iter_t *iter);

#endif
