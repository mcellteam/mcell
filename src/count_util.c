/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
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

/**************************************************************************\
** File: count_util.c                                                     **
**                                                                        **
** Purpose: Handles counting of interesting events                        **
**                                                                        **
** Testing status: untested.                                              **
\**************************************************************************/

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"
#include "logging.h"
#include "grid_util.h"
#include "wall_util.h"
#include "vol_util.h"
#include "count_util.h"
#include "react_output.h"
//#include "util.h"
#include "sym_table.h"
#include "dyngeom_parse_extras.h"

/* Instantiate a request to track a particular quantity */
static int instantiate_count_request(
  int dyn_geom_flag, struct output_request *request, int count_hashmask,
  struct counter **count_hash, struct mem_helper *trig_request_mem,
  double *elapsed_time, struct mem_helper *counter_mem);

/* Create a new counter data structure */
static struct counter *create_new_counter(struct region *where, void *who,
  byte what, struct periodic_image *img, struct mem_helper *counter_mem);

/* Pare down the region lists, annihilating any regions which appear in both
 * lists. */
static void clean_region_lists(struct subvolume *my_sv,
                               struct region_list **p_all_regs,
                               struct region_list **p_all_antiregs);

/* Check if the object corresponding to a particular symbol has been referenced
 * directly or indirectly from one of the INSTANTIATE blocks in the model.
 */
/*static int is_object_instantiated(struct sym_entry *entry,*/
/*                                  struct object *root_instance);*/

/* Find the list of regions enclosing a particular point. given a particular
 * starting point and starting region list. */
static int find_enclosing_regions(struct volume *world, struct vector3 *loc,
                                  struct vector3 *start,
                                  struct region_list **rlp,
                                  struct region_list **arlp,
                                  struct mem_helper *rmem);

void count_region_list(
    struct volume *world,
    struct region_list *regions,
    struct surface_molecule *sm,
    struct vector3 *where,
    int count_hashmask,
    int inc,
    struct periodic_image *previous_box);

/*************************************************************************
eps_equals:
   In: two doubles
   Out: 1 if they are equal to within some small tolerance, 0 otherwise
*************************************************************************/
static int eps_equals(double x, double y) {
  double mag = fabs(x);
  if (mag < fabs(y))
    mag = fabs(y);

  double diff = fabs(x - y);
  return diff < EPS_C * (mag + 1.0);
}

/*************************************************************************
dup_region_list:
   In: a list of regions
       memory handler to use for duplicated regions
   Out: The duplicated list of regions, or NULL on a memory allocation
        error.
*************************************************************************/
static struct region_list *dup_region_list(struct region_list *r,
                                           struct mem_helper *mh) {
  struct region_list *nr, *rp, *r0;

  if (r == NULL)
    return NULL;

  r0 = rp = NULL;
  while (r != NULL) {
    nr = (struct region_list *)CHECKED_MEM_GET(mh, "region list entry");
    nr->next = NULL;
    nr->reg = r->reg;
    if (rp == NULL)
      r0 = rp = nr;
    else {
      rp->next = nr;
      rp = nr;
    }

    r = r->next;
  }

  return r0;
}

/*************************************************************************
region_listed:
   In: list of regions
       one specific region we're interested in
   Out: 1 if the region is in the list.  0 if not.
*************************************************************************/
int region_listed(struct region_list *rl, struct region *r) {
  while (rl != NULL) {
    if (rl->reg == r)
      return 1;
    rl = rl->next;
  }
  return 0;
}

/*************************************************************************
count_region_update:
   In: world: simulation state 
       sp: species of thing that hit
       periodic_box: periodic box of molecule being counted
       rl: region list for the wall we hit
       dir: direction of impact relative to surface normal for volume molecule,
         or relative to the region border for surface molecule (inside out = 1,
         ouside in = 0)
       crossed: whether we crossed or not
       loc: location of the hit (for triggers)
       t: time of the hit (for triggers)
   Out: Returns none.
        Appropriate counters are updated, that is, hit counters are updated
        according to which side was hit, and crossings counters and counts
        within enclosed regions are updated if the surface was crossed.
*************************************************************************/
void count_region_update(
    struct volume *world,
    struct species *sp,
    struct periodic_image *periodic_box,
    struct region_list *rl,
    int direction,
    int crossed,
    struct vector3 *loc,
    double t) {

  int count_hits = 0;
  double hits_to_ccn = 0;
  if ((sp->flags & COUNT_HITS) && ((sp->flags & NOT_FREE) == 0)) {
    count_hits = 1;
    hits_to_ccn = sp->time_step *
                  2.9432976599069717358e-3 / /* 1e6*sqrt(MY_PI)/(1e-15*N_AV) */
                  (sp->space_step * world->length_unit * world->length_unit *
                   world->length_unit);
  }

  struct counter *hit_count = NULL;
  for (; rl != NULL; rl = rl->next) {
    if (!(rl->reg->flags & COUNT_SOME_MASK)) {
      continue;
    }

    if (!(rl->reg->flags & sp->flags & (COUNT_HITS | COUNT_CONTENTS | COUNT_ENCLOSED))) {
      continue;
    }

    int hash_bin = (rl->reg->hashval + sp->hashval) & world->count_hashmask;
    for (hit_count = world->count_hash[hash_bin]; hit_count != NULL;
         hit_count = hit_count->next) {

      if (!(hit_count->reg_type == rl->reg && hit_count->target == sp)) {
        continue;
      }

      // count only in the relevant periodic box
      if (world->periodic_box_obj && !world->periodic_traditional) {
        if (!periodic_boxes_are_identical(periodic_box, hit_count->periodic_box)) {
          continue;
        }
        else {
          struct vector3 pos_output = {0.0, 0.0, 0.0};
          convert_relative_to_abs_PBC_coords(
              world->periodic_box_obj,
              periodic_box,
              world->periodic_traditional,
              loc,
              &pos_output);
          loc->x = pos_output.x;   
          loc->y = pos_output.y;   
          loc->z = pos_output.z;   
        }
      }

      if (crossed) {
        if (direction == 1) {
          if (hit_count->counter_type & TRIG_COUNTER) {
            hit_count->data.trig.t_event = (double)world->current_iterations + t;
            hit_count->data.trig.orient = 0;
            if (rl->reg->flags & sp->flags & COUNT_HITS) {
              fire_count_event(world, hit_count, 1, loc,
                               REPORT_FRONT_HITS | REPORT_TRIGGER);

              fire_count_event(world, hit_count, 1, loc,
                               REPORT_FRONT_CROSSINGS | REPORT_TRIGGER);
            }
            if (rl->reg->flags & sp->flags & COUNT_CONTENTS) {
              fire_count_event(world, hit_count, 1, loc,
                               REPORT_ENCLOSED | REPORT_CONTENTS |
                                   REPORT_TRIGGER);
            }
          } else {
            if (rl->reg->flags & sp->flags & COUNT_HITS) {
              hit_count->data.move.front_hits++;
              hit_count->data.move.front_to_back++;
            }
            if (rl->reg->flags & sp->flags & COUNT_CONTENTS) {
              hit_count->data.move.n_enclosed++;
            }
          }
        } else {
          if (hit_count->counter_type & TRIG_COUNTER) {
            hit_count->data.trig.t_event = (double)world->current_iterations + t;
            hit_count->data.trig.orient = 0;
            if (rl->reg->flags & sp->flags & COUNT_HITS) {
              fire_count_event(world, hit_count, 1, loc,
                               REPORT_BACK_HITS | REPORT_TRIGGER);
              fire_count_event(world, hit_count, 1, loc,
                               REPORT_BACK_CROSSINGS | REPORT_TRIGGER);
            }
            if (rl->reg->flags & sp->flags & COUNT_CONTENTS) {
              fire_count_event(world, hit_count, -1, loc,
                               REPORT_ENCLOSED | REPORT_CONTENTS |
                                   REPORT_TRIGGER);
            }
          } else {
            if (rl->reg->flags & sp->flags & COUNT_HITS) {
              hit_count->data.move.back_hits++;
              hit_count->data.move.back_to_front++;
            }
            if (rl->reg->flags & sp->flags & COUNT_CONTENTS) {
              hit_count->data.move.n_enclosed--;
            }
          }
        }
      } else if (rl->reg->flags & sp->flags & COUNT_HITS) {
      /* Didn't cross, only hits might update */
        if (direction == 1) {
          if (hit_count->counter_type & TRIG_COUNTER) {
            hit_count->data.trig.t_event = (double)world->current_iterations + t;
            hit_count->data.trig.orient = 0;
            fire_count_event(world, hit_count, 1, loc,
                             REPORT_FRONT_HITS | REPORT_TRIGGER);
          } else {
            hit_count->data.move.front_hits++;
          }
        } else {
          if (hit_count->counter_type & TRIG_COUNTER) {
            hit_count->data.trig.t_event = (double)world->current_iterations + t;
            hit_count->data.trig.orient = 0;
            fire_count_event(world, hit_count, 1, loc,
                             REPORT_BACK_HITS | REPORT_TRIGGER);
          } else
            hit_count->data.move.back_hits++;
        }
      }
      if ((count_hits && rl->reg->area != 0.0) &&
          ((sp->flags & NOT_FREE) == 0)) {
        if ((hit_count->counter_type & TRIG_COUNTER) == 0) {
          hit_count->data.move.scaled_hits += hits_to_ccn / rl->reg->area;
        }
      }
    }
  }
}

/**************************************************************************
count_region_border_update:
  In: world: simulation state
      sp: species of thing that hit
      hd_info: information about the hit (linked list of "hit_data")
  Out: Returns none.
       Appropriate counters are updated, that is,
       hit counters are updated according to which side was hit,
       and crossings counters  are updated if the region border was crossed.
       We consider FRONT_HITS/FRONT_CROSSINGS when the molecule
       crosses region border "inside out" and BACK_HITS/BACK_CROSSINGS
       when the molecule hits region border "outside in".
**************************************************************************/
void count_region_border_update(struct volume *world, struct species *sp,
                                struct hit_data *hd_info) {

  assert((sp->flags & NOT_FREE) != 0);

  for (struct hit_data *hd = hd_info; hd != NULL; hd = hd->next) {
    for (struct region_list *rl = hd->count_regions; rl != NULL; rl = rl->next) {
      if ((rl->reg->flags & COUNT_SOME_MASK) &&
          (rl->reg->flags & sp->flags & COUNT_HITS)) {

        int hash_bin = (rl->reg->hashval + sp->hashval) & world->count_hashmask;

        for (struct counter *hit_count = world->count_hash[hash_bin];
          hit_count != NULL; hit_count = hit_count->next) {

          if (((hit_count->reg_type != rl->reg) || (hit_count->target != sp))) {
            continue;
          }

          if ((hit_count->orientation != ORIENT_NOT_SET) &&
              (hit_count->orientation != hd->orientation) &&
              (hit_count->orientation != 0)) {
            continue;
          }

          if (hit_count->counter_type & TRIG_COUNTER) {
            hit_count->data.trig.t_event = hd->t;
            hit_count->data.trig.orient = 0;
            if (hd->direction == 1) {
              fire_count_event(world, hit_count, 1, &(hd->loc),
                REPORT_FRONT_HITS | REPORT_TRIGGER);
              if (hd->crossed) {
                fire_count_event(world, hit_count, 1, &(hd->loc),
                  REPORT_FRONT_CROSSINGS | REPORT_TRIGGER);
              }
            } else {
              fire_count_event(world, hit_count, 1, &(hd->loc),
                REPORT_BACK_HITS | REPORT_TRIGGER);
              if (hd->crossed) {
                fire_count_event(world, hit_count, 1, &(hd->loc),
                  REPORT_BACK_CROSSINGS | REPORT_TRIGGER);
              }
            }
          } else {
            if (hd->direction == 1) {
              hit_count->data.move.front_hits++;
              if (hd->crossed) {
                hit_count->data.move.front_to_back++;
              }
            } else {
              hit_count->data.move.back_hits++;
              if (hd->crossed) {
                hit_count->data.move.back_to_front++;
              }
            }
          }
        }
      }
    }
  } /* end for (hd...) */
}

/*************************************************************************
count_region_from_scratch:
   In: world: simulation state 
       am: molecule to count, or NULL
       rxpn: reaction pathname to count, or NULL
       n: number of these to count
       loc: location at which to count them (may be NULL)
       my_wall: wall at which this happened (may be NULL)
       t: time of the hit (for triggers)
       periodic_box:
   Out: Returns zero on success and 1 on failure.
        Appropriate counters are updated and triggers are fired.
   Note: At least one of molecule or rxn pathname must be non-NULL; if other
         inputs are NULL, sensible values will be guessed (which may themselves
         be NULL). This routine is not super-fast for volume counts (enclosed
         counts) since it has to dynamically create and test lists of enclosing
         regions.
*************************************************************************/
void count_region_from_scratch(struct volume *world,
                               struct abstract_molecule *am,
                               struct rxn_pathname *rxpn,
                               int n,
                               struct vector3 *loc,
                               struct wall *my_wall,
                               double t,
                               struct periodic_image *periodic_box) {
  struct region_list *rl, *arl, *nrl, *narl; /*a=anti n=new*/
  struct counter *c;
  void *target; /* what we're counting: am->properties or rxpn */
  int hashval;  /* Hash value of what we're counting */
  double t_hit, t_sv_hit;
  struct vector3 delta, hit; /* For raytracing */
  struct vector3 xyz_loc;          /* Computed location of mol if loc==NULL */
  byte count_flags;
  int pos_or_neg;        /* Sign of count (neg for antiregions) */
  int orient = SHRT_MIN; /* orientation of the molecule also serves as a flag
                            for triggering reactions  */

  /* Set up values and fill in things the calling function left out */
  if (rxpn != NULL) {
    hashval = rxpn->hashval;
    target = rxpn;
    count_flags = REPORT_RXNS;
  } else {
    hashval = am->properties->hashval;
    target = am->properties;
    count_flags = REPORT_CONTENTS;
    if (loc == NULL) {
      if (am->properties->flags & ON_GRID) {
        uv2xyz(&(((struct surface_molecule *)am)->s_pos),
               ((struct surface_molecule *)am)->grid->surface, &xyz_loc);
        loc = &xyz_loc;
      } else
        loc = &(((struct volume_molecule *)am)->pos);
    }

    if (my_wall == NULL && (am->properties->flags & ON_GRID) != 0) {
      my_wall = ((struct surface_molecule *)am)->grid->surface;
    }

    if (am->properties->flags & ON_GRID) {
      orient = ((struct surface_molecule *)am)->orient;
    } else {
      orient = 0;
    }
  }

  /* Count surface molecules and reactions on surfaces--easy */
  if (my_wall != NULL && (my_wall->flags & COUNT_CONTENTS) != 0) {
    for (rl = my_wall->counting_regions; rl != NULL; rl = rl->next) {
      int hash_bin = (hashval + rl->reg->hashval) & world->count_hashmask;
      for (c = world->count_hash[hash_bin]; c != NULL; c = c->next) {
        if (c->target == target && c->reg_type == rl->reg &&
            (c->counter_type & ENCLOSING_COUNTER) == 0) {
          if (c->counter_type & TRIG_COUNTER) {
            c->data.trig.t_event = t;
            c->data.trig.orient = orient;
            // XXX: may need to convert loc for PBCs
            fire_count_event(world, c, n, loc, count_flags | REPORT_TRIGGER);
          } else if (rxpn == NULL) {
            if (am->properties->flags & ON_GRID) {
              if ((c->orientation == ORIENT_NOT_SET) ||
                  (c->orientation == orient) || (c->orientation == 0)) {
                // count only in the relevant periodic box
                if (periodic_boxes_are_identical(am->periodic_box, c->periodic_box)) {
                  c->data.move.n_at += n;
                }
              }
            } else {
              c->data.move.n_at += n;
            }
          } else if ((rxpn != NULL) && (periodic_boxes_are_identical(periodic_box, c->periodic_box))) {
            c->data.rx.n_rxn_at += n;
          }
        }
      }
    }
  }

  /* Waypoints must have been placed in order for the following code to work. */
  if (!world->place_waypoints_flag) {
    assert (am != NULL && (am->properties->flags & COUNT_ENCLOSED) == 0 &&
      (am->properties->flags & NOT_FREE) != 0);
    return;
  }

  /* Count volume molecules, vol reactions, and surface stuff that is
   * enclosed--hard!!*/
  if (am == NULL || (am->properties->flags & COUNT_ENCLOSED) != 0 ||
      (am->properties->flags & NOT_FREE) == 0) {
    const int px = bisect(world->x_partitions, world->nx_parts, loc->x);
    const int py = bisect(world->y_partitions, world->ny_parts, loc->y);
    const int pz = bisect(world->z_partitions, world->nz_parts, loc->z);
    const int this_sv =
        pz + (world->nz_parts - 1) * (py + (world->ny_parts - 1) * px);
    struct waypoint *wp = &(world->waypoints[this_sv]);
    struct subvolume *my_sv = &(world->subvol[this_sv]);

    struct vector3 here = {.x = wp->loc.x, .y = wp->loc.y, .z = wp->loc.z};

    struct region_list *all_regs = NULL;
    struct region_list *all_antiregs = NULL;

    /* Copy all the potentially relevant regions from the nearest waypoint */
    for (rl = wp->regions; rl != NULL; rl = rl->next) {
      if (rl->reg == NULL)
        continue;
      int hash_bin = (hashval + rl->reg->hashval) & world->count_hashmask;
      if (world->count_hash[hash_bin] == NULL)
        continue; /* Won't count on this region so ignore it */

      nrl = (struct region_list *)CHECKED_MEM_GET(
          my_sv->local_storage->regl, "list of enclosing regions for count");
      nrl->reg = rl->reg;
      nrl->next = all_regs;
      all_regs = nrl;
    }

    /* And all the antiregions (regions crossed from inside to outside only) */
    for (arl = wp->antiregions; arl != NULL; arl = arl->next) {
      int hash_bin = (hashval + arl->reg->hashval) & world->count_hashmask;
      if (world->count_hash[hash_bin] == NULL)
        continue; /* Won't count on this region so ignore it */

      narl = (struct region_list *)CHECKED_MEM_GET(
          my_sv->local_storage->regl, "list of enclosing regions for count");
      narl->reg = arl->reg;
      narl->next = all_antiregs;
      all_antiregs = narl;
    }

    /* Raytrace across any walls from waypoint to us and add to region lists */
    for (struct subvolume *sv = &(world->subvol[this_sv]); sv != NULL;
         sv = next_subvol(&here, &delta, sv, world->x_fineparts,
                          world->y_fineparts, world->z_fineparts,
                          world->ny_parts, world->nz_parts)) {
      delta.x = loc->x - here.x;
      delta.y = loc->y - here.y;
      delta.z = loc->z - here.z;

      t_sv_hit = collide_sv_time(&here, &delta, sv, world->x_fineparts,
                                 world->y_fineparts, world->z_fineparts);
      if (t_sv_hit > 1.0)
        t_sv_hit = 1.0;

      for (struct wall_list *wl = sv->wall_head; wl != NULL; wl = wl->next) {
        /* Skip wall that we are on unless we're a volume molecule */
        if (my_wall == wl->this_wall &&
            (am == NULL || (am->properties->flags & NOT_FREE))) {
          continue;
        }

        if (wl->this_wall->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) {
          int hit_code = collide_wall(&here, &delta, wl->this_wall, &t_hit,
                                      &hit, 0, world->rng, world->notify,
                                      &(world->ray_polygon_tests));
          if (hit_code == COLLIDE_MISS) {
            continue;
          }

          world->ray_polygon_colls++;
          if (t_hit <= t_sv_hit && (hit.x - loc->x) * delta.x +
            (hit.y - loc->y) * delta.y + (hit.z - loc->z) * delta.z < 0) {
            for (rl = wl->this_wall->counting_regions; rl != NULL;
                 rl = rl->next) {
              if ((rl->reg->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) != 0) {
                int hash_bin =
                    (hashval + rl->reg->hashval) & world->count_hashmask;
                if (world->count_hash[hash_bin] == NULL) {
                  continue; /* Won't count on this region so ignore it */
                }
                nrl = (struct region_list *)CHECKED_MEM_GET(
                    my_sv->local_storage->regl,
                    "list of enclosing regions for count");
                nrl->reg = rl->reg;
                if (hit_code == COLLIDE_FRONT) {
                  nrl->next = all_regs;
                  all_regs = nrl;
                } else if (hit_code == COLLIDE_BACK) {
                  nrl->next = all_antiregs;
                  all_antiregs = nrl;
                }
              }
            }
          }
        }
      }
    }

    /* Clean up region lists */
    if (all_regs != NULL && all_antiregs != NULL)
      clean_region_lists(my_sv, &all_regs, &all_antiregs);

    /* Actually check the regions here */
    count_flags |= REPORT_ENCLOSED;

    for (nrl = all_regs; nrl != NULL; nrl = (nrl == all_regs) ? all_antiregs : NULL)
    /* Trick so we don't need to duplicate this code */ {
      if (nrl == all_regs) {
        pos_or_neg = 1;
      } else {
        pos_or_neg = -1;
      }
      for (rl = nrl; rl != NULL; rl = rl->next) {
        int hash_bin = (hashval + rl->reg->hashval) & world->count_hashmask;
        for (c = world->count_hash[hash_bin]; c != NULL; c = c->next) {
          if (am != NULL
              && !periodic_boxes_are_identical(c->periodic_box, periodic_box)) {
            continue;
          }
          if (rxpn != NULL
              && !periodic_boxes_are_identical(c->periodic_box, periodic_box)) {
            continue;
          }
          if (c->target == target && c->reg_type == rl->reg &&
              ((c->counter_type & ENCLOSING_COUNTER) != 0 ||
               (am != NULL && (am->properties->flags & ON_GRID) == 0)) &&
              (my_wall == NULL ||
               (am != NULL && (am->properties->flags & NOT_FREE) == 0) ||
               !region_listed(my_wall->counting_regions, rl->reg))) {
            if (c->counter_type & TRIG_COUNTER) {
              c->data.trig.t_event = t;
              c->data.trig.orient = orient;
              // Don't count triggers after a dynamic geometry event
              if (!world->dynamic_geometry_flag) {
                // XXX: may need to convert loc for PBCs
                fire_count_event(world, c, n * pos_or_neg, loc,
                                 count_flags | REPORT_TRIGGER);
              }
            } else if (rxpn == NULL) {
              if (am->properties->flags & ON_GRID) {
                if ((c->orientation == ORIENT_NOT_SET) ||
                    (c->orientation == orient) || (c->orientation == 0)) {
                  c->data.move.n_enclosed += n * pos_or_neg;
                }
              } else {
                c->data.move.n_enclosed += n * pos_or_neg;
              }
            } else {
              c->data.rx.n_rxn_enclosed += n * pos_or_neg;
            }
          }
        }
      }
    }

    /* Free region memory */
    if (all_regs != NULL)
      mem_put_list(my_sv->local_storage->regl, all_regs);
    if (all_antiregs != NULL)
      mem_put_list(my_sv->local_storage->regl, all_antiregs);
  }
}

/*************************************************************************
count_moved_surface_mol:
   In: world: simulation state 
       sm: molecule to count
       sg: new grid for molecule
       loc: new location on that grid
       count_hashmask:
       count_hash:
       ray_polygon_colls:
   Out: Returns zero on success and 1 on failure.
        Appropriate counters are updated and triggers are fired.
   Note: This routine is not super-fast for enclosed counts for
         surface molecules since it raytraces without using waypoints.
*************************************************************************/
void count_moved_surface_mol(
    struct volume *world,
    struct surface_molecule *sm,
    struct surface_grid *sg,
    struct vector2 *loc,
    int count_hashmask,
    struct counter **count_hash,
    long long *ray_polygon_colls,
    struct periodic_image *previous_box) {

  struct vector3 origin;
  struct vector3 target;
  struct region_list *pos_regs = NULL;
  struct region_list *neg_regs = NULL;
  struct region_list *nrl = NULL;
  struct region_list *prl = NULL;
  struct region_list *rl = NULL;
  struct region_list *rl2 = NULL;
  struct storage *stor = sm->grid->surface->birthplace;
  // Different grids implies different walls, so we might have changed regions 
  /*if (sm->grid != sg) {*/
  if ((sm->grid != sg) ||
      (previous_box != NULL && !world->periodic_traditional)) {
    int delete_me = 0;
    if ((sm->grid->surface->flags & COUNT_CONTENTS) != 0 &&
      (sg->surface->flags & COUNT_CONTENTS) != 0) {

      delete_me = 1;
      nrl = sm->grid->surface->counting_regions;
      prl = sg->surface->counting_regions;
      while (prl != NULL && nrl != NULL) {
        if (prl->reg == nrl->reg) {

          rl = (struct region_list *)CHECKED_MEM_GET(stor->regl, "region list entry");
          rl->next = pos_regs;
          rl->reg = prl->reg;
          pos_regs = rl;
          prl = prl->next;

          rl2 = (struct region_list *)CHECKED_MEM_GET(stor->regl, "region list entry");
          rl2->next = neg_regs;
          rl2->reg = nrl->reg;
          neg_regs = rl2;
          nrl = nrl->next;

          continue;
        }
        while (prl != NULL && prl->reg < nrl->reg) { /* Entering these regions */
          rl = (struct region_list *)CHECKED_MEM_GET(stor->regl, "region list entry");
          rl->next = pos_regs;
          rl->reg = prl->reg;
          pos_regs = rl;
          prl = prl->next;
        }
        while (nrl != NULL && (prl == NULL || nrl->reg < prl->reg)) { /* Leaving these regions */
          rl = (struct region_list *)CHECKED_MEM_GET(stor->regl, "region list entry");
          rl->next = neg_regs;
          rl->reg = nrl->reg;
          neg_regs = rl;
          nrl = nrl->next;
        }
      }
      /* If we exhaust all negative regions before we've exhausted all
       * positive regions, the above loop will terminate, leaving some
       * regions uncounted. */
      while (prl != NULL) {
        rl = (struct region_list *)CHECKED_MEM_GET(stor->regl,
                                                   "region list entry");
        rl->next = pos_regs;
        rl->reg = prl->reg;
        pos_regs = rl;
        prl = prl->next;
      }

      /* I don't think this can happen, but it could potentially happen
       * if, say, prl started off NULL (i.e. one of the grids belonged
       * to no counting regions at all). */
      while (nrl != NULL) {
        rl = (struct region_list *)CHECKED_MEM_GET(stor->regl,
                                                   "region list entry");
        rl->next = neg_regs;
        rl->reg = nrl->reg;
        neg_regs = rl;
        nrl = nrl->next;
      }
    } else if (sm->grid->surface->flags & COUNT_CONTENTS) {
      neg_regs = sm->grid->surface->counting_regions;
    } else if (sg->surface->flags & COUNT_CONTENTS) {
      pos_regs = sg->surface->counting_regions;
    }

    // Different walls and possibly different periodic box
    if (pos_regs != NULL) {
      uv2xyz(loc, sg->surface, &target);
      count_region_list(
          world, pos_regs, sm, &target, count_hashmask, 1, previous_box);
    }
    if (neg_regs != NULL) {
      uv2xyz(&(sm->s_pos), sm->grid->surface, &origin);
      count_region_list(
          world, neg_regs, sm, &origin, count_hashmask, -1, previous_box);
    }

    if (delete_me) {
      if (pos_regs != NULL)
        mem_put_list(stor->regl, pos_regs);
      if (neg_regs != NULL)
        mem_put_list(stor->regl, neg_regs);
    }
  }

  /* Have to raytrace */
  uv2xyz(&(sm->s_pos), sm->grid->surface, &origin);
  uv2xyz(loc, sg->surface, &target);
  if ((sm->properties->flags & COUNT_ENCLOSED) &&
      (periodic_boxes_are_identical(previous_box, sm->periodic_box))) {

    pos_regs = neg_regs = NULL;
    struct vector3 delta = {target.x - origin.x, target.y - origin.y, target.z - origin.z};
    struct vector3 here = origin;

    /* Collect all the relevant regions we pass through */
    for (struct subvolume *sv = find_subvolume(world, &origin, NULL); sv != NULL;
         sv = next_subvol(&here, &delta, sv, world->x_fineparts, world->y_fineparts,
          world->z_fineparts, world->ny_parts, world->nz_parts)) {

      int j = 0;
      for (struct wall_list *wl = sv->wall_head; wl != NULL; wl = wl->next) {
        if (wl->this_wall == sm->grid->surface || wl->this_wall == sg->surface) {
          continue; // ignore the origin and target walls
        }

        struct vector3 hit = {0.0, 0.0, 0.0};
        double t = 0.0;
        j = collide_wall(&here, &delta, wl->this_wall, &t, &hit, 0, world->rng,
          world->notify, &(world->ray_polygon_tests));

        /* we only consider the collision if it happens in the current subvolume.
           Otherwise we may double count collision for walls that span multiple
           subvolumes */
        if (!inside_subvolume(&hit, sv, world->x_fineparts, world->y_fineparts,
          world->z_fineparts)) {
          continue;
        }

        if (j != COLLIDE_MISS) {
          (*ray_polygon_colls)++;
        }

        /* check that hit is encountered before we reach the target */
        if (j != COLLIDE_MISS && (hit.x - target.x) * delta.x +
          (hit.y - target.y) * delta.y + (hit.z - target.z) * delta.z < 0) {
          for (rl = wl->this_wall->counting_regions; rl != NULL; rl = rl->next) {
            if ((rl->reg->flags & COUNT_ENCLOSED) == 0) {
              continue; /* Only ENCLOSED counted here */
            }

            if (j == COLLIDE_FRONT) {
              prl = (struct region_list *)CHECKED_MEM_GET(stor->regl, "region list entry");
              prl->reg = rl->reg;
              prl->next = pos_regs;
              pos_regs = prl;
            } else if (j == COLLIDE_BACK) {
              nrl = (struct region_list *)CHECKED_MEM_GET(stor->regl, "region list entry");
              nrl->reg = rl->reg;
              nrl->next = neg_regs;
              neg_regs = nrl;
            }
          }
        }
      }
    }

    if (pos_regs != NULL) {
      pos_regs = (struct region_list *)void_list_sort((struct void_list *)pos_regs);
    }
    if (neg_regs != NULL) {
      neg_regs = (struct region_list *)void_list_sort((struct void_list *)neg_regs);
    }

    prl = pos_regs;
    nrl = neg_regs;
    struct vector3 *where = NULL;
    int n;
    while (prl != NULL || nrl != NULL) {
      if (prl == NULL) {
        rl = nrl;
        nrl = nrl->next;
        n = -1;
        where = &origin;
      } else if (nrl == NULL) {
        rl = prl;
        prl = prl->next;
        n = 1;
        where = &target;
      } else if (prl->reg < nrl->reg) {
        rl = prl;
        prl = prl->next;
        n = 1;
        where = &origin;
      } else if (nrl->reg < prl->reg) {
        rl = nrl;
        nrl = nrl->next;
        n = -1;
        where = &target;
      } else {
        rl = NULL;
        prl = prl->next;
        nrl = nrl->next;
      }

      if (rl == NULL) {
        continue;
      }

      int hash_bin = (sm->properties->hashval + rl->reg->hashval) & count_hashmask;
      for (struct counter *c = count_hash[hash_bin]; c != NULL; c = c->next) {
        if (c->target == sm->properties && c->reg_type == rl->reg &&
            (c->counter_type & ENCLOSING_COUNTER) != 0) {

          assert(!region_listed(sm->grid->surface->counting_regions, rl->reg));
          assert(!region_listed(sg->surface->counting_regions, rl->reg));

          if (c->counter_type & TRIG_COUNTER) {
            c->data.trig.t_event = sm->t;
            c->data.trig.orient = sm->orient;
            fire_count_event(world, c, n, where, REPORT_CONTENTS | REPORT_ENCLOSED |
              REPORT_TRIGGER);
          } else if ((c->orientation == ORIENT_NOT_SET) ||
                     (c->orientation == sm->orient) ||
                     (c->orientation == 0)) {
            /*c->data.move.n_enclosed += n;*/
            if (periodic_boxes_are_identical(c->periodic_box, sm->periodic_box)) {
              c->data.move.n_enclosed += n;
            }
          }
        }
      }
    }

    if (pos_regs != NULL)
      mem_put_list(stor->regl, pos_regs);
    if (neg_regs != NULL)
      mem_put_list(stor->regl, neg_regs);
  }
  else if ((sm->properties->flags & COUNT_ENCLOSED) &&
      (!periodic_boxes_are_identical(previous_box, sm->periodic_box))) {
    // Increment count of where we are going now (target)
    count_region_from_scratch(world, (struct abstract_molecule *)sm, NULL, 1, &target, NULL, 1.0, sm->periodic_box);
    // Decrement count of where we were before (origin)
    count_region_from_scratch(world, (struct abstract_molecule *)sm, NULL, -1, &origin, NULL, 1.0, previous_box);
  }
}

/*************************************************************************
fire_count_event:
   In: world: simulation state
       event: counter of thing that just happened (trigger of some sort)
       n: number of times that thing happened (or hit direction for triggers)
       where: location where it happened
       what: what happened (Report Type Flags)
   Out: None
*************************************************************************/
void fire_count_event(struct volume *world, struct counter *event, int n,
                      struct vector3 *where, byte what) {

  short flags;
  if ((what & REPORT_TYPE_MASK) == REPORT_RXNS)
    flags = TRIG_IS_RXN;
  else if ((what & REPORT_TYPE_MASK) == REPORT_CONTENTS)
    flags = 0;
  else
    flags = TRIG_IS_HIT;

  byte whatelse = what;
  if ((what & REPORT_TYPE_MASK) == REPORT_FRONT_HITS)
    whatelse = (what - REPORT_FRONT_HITS) | REPORT_ALL_HITS;
  else if ((what & REPORT_TYPE_MASK) == REPORT_BACK_HITS)
    whatelse = (what - REPORT_BACK_HITS) | REPORT_ALL_HITS;
  else if ((what & REPORT_TYPE_MASK) == REPORT_FRONT_CROSSINGS)
    whatelse = (what - REPORT_FRONT_CROSSINGS) | REPORT_ALL_CROSSINGS;
  else if ((what & REPORT_TYPE_MASK) == REPORT_BACK_CROSSINGS)
    whatelse = (what - REPORT_BACK_CROSSINGS) | REPORT_ALL_CROSSINGS;

  struct trigger_request *tr;
  for (tr = event->data.trig.listeners; tr != NULL; tr = tr->next) {
    if (tr->ear->report_type == what) {
      memcpy(&(event->data.trig.loc), where, sizeof(struct vector3));
      if ((what & REPORT_TYPE_MASK) == REPORT_FRONT_HITS ||
          (what & REPORT_TYPE_MASK) == REPORT_FRONT_CROSSINGS) {
        add_trigger_output(world, event, tr->ear, n, flags);
      } else if ((what & REPORT_TYPE_MASK) == REPORT_BACK_HITS ||
                 (what & REPORT_TYPE_MASK) == REPORT_BACK_CROSSINGS) {
        add_trigger_output(world, event, tr->ear, -n, flags);
      } else {
        add_trigger_output(world, event, tr->ear, n, flags);
      }

    } else if (tr->ear->report_type == whatelse) {
      memcpy(&(event->data.trig.loc), where, sizeof(struct vector3));
      if ((what & REPORT_TYPE_MASK) == REPORT_FRONT_HITS ||
          (what & REPORT_TYPE_MASK) == REPORT_FRONT_CROSSINGS) {
        add_trigger_output(world, event, tr->ear, n, flags);
      } else {
        add_trigger_output(world, event, tr->ear, -n, flags);
      }
    }
  }
}

/*************************************************************************
find_enclosing_regions:
   In: world: simulation state 
       loc: location we want to end up
       start: starting position
       rlp: list of regions we're inside at the starting position
       arlp: list of inside-out regions we're "outside" at the starting position
       rmem: memory handler to store lists of regions
   Out: 0 on success, 1 on memory allocation error.  The region and
        inside-out region lists are updated to be correct at the ending
        position.
*************************************************************************/
static int find_enclosing_regions(struct volume *world,
                                  struct vector3 *loc,
                                  struct vector3 *start,
                                  struct region_list **rlp,
                                  struct region_list **arlp,
                                  struct mem_helper *rmem) {
  struct vector3 outside, delta, hit;
  struct wall_list *wl;
  struct region_list *rl, *arl;
  struct region_list *trl, *tarl, *xrl, *yrl, *nrl;
  struct wall_list dummy;

  rl = *rlp;
  arl = *arlp;

  if (start == NULL || (distinguishable(loc->x, start->x, EPS_C)) ||
      (distinguishable(loc->y, start->y, EPS_C)) || loc->z < start->z) {
    outside.x = loc->x;
    outside.y = loc->y;
    outside.z = (world->z_partitions[0] + world->z_partitions[1]) / 2;
  } else {
    outside.x = start->x;
    outside.y = start->y;
    outside.z = start->z;
  }

  delta.x = 0.0;
  delta.y = 0.0;
  delta.z = loc->z - outside.z;

  struct subvolume *sv = find_subvolume(world, &outside, NULL);
  struct subvolume *svt = find_subvolume(world, loc, NULL);
  int traveling = 1;

  double t;
  while (traveling) {
    tarl = trl = NULL;
    double t_hit_sv = collide_sv_time(&outside, &delta, sv, world->x_fineparts,
                                      world->y_fineparts, world->z_fineparts);

    for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
      int hit_code =
          collide_wall(&outside, &delta, wl->this_wall, &t, &hit, 0, world->rng,
                       world->notify, &(world->ray_polygon_tests));

      if ((hit_code != COLLIDE_MISS) &&
          (world->notify->final_summary == NOTIFY_FULL)) {
        world->ray_polygon_colls++;
      }

      if (hit_code == COLLIDE_REDO) {
        while (trl != NULL) {
          xrl = trl->next;
          mem_put(rmem, trl);
          trl = xrl;
        }
        while (tarl != NULL) {
          xrl = tarl->next;
          mem_put(rmem, tarl);
          tarl = xrl;
        }
        dummy.next = sv->wall_head;
        wl = &dummy;
        continue; /* Trick to restart for loop */
      } else if (hit_code == COLLIDE_MISS || !(t >= 0 && t < 1.0) ||
                 t > t_hit_sv ||
                 (wl->this_wall->flags &
                  (COUNT_CONTENTS | COUNT_RXNS | COUNT_ENCLOSED)) == 0 ||
                 (hit.x - outside.x) * delta.x + (hit.y - outside.y) * delta.y +
                         (hit.z - outside.z) * delta.z <
                     0)
        continue;
      else {
        for (xrl = wl->this_wall->counting_regions; xrl != NULL;
             xrl = xrl->next) {
          if ((xrl->reg->flags &
               (COUNT_CONTENTS | COUNT_RXNS | COUNT_ENCLOSED)) != 0) {
            nrl = (struct region_list *)CHECKED_MEM_GET(rmem,
                                                        "region list entry");
            nrl->reg = xrl->reg;

            if (hit_code == COLLIDE_BACK) {
              nrl->next = tarl;
              tarl = nrl;
            } else {
              nrl->next = trl;
              trl = nrl;
            }
          }
        }
      }
    }

    xrl = trl;
    while (trl != NULL) {
      nrl = NULL;
      yrl = arl;
      while (yrl != NULL) {
        if (xrl->reg == yrl->reg) {
          if (nrl == NULL) {
            arl = yrl->next;
            mem_put(rmem, yrl);
          } else {
            nrl->next = yrl->next;
            mem_put(rmem, yrl);
          }
          trl = trl->next;
          mem_put(rmem, xrl);
          xrl = NULL;
          break;
        } else {
          nrl = yrl;
          yrl = yrl->next;
        }
      }
      if (xrl != NULL) {
        trl = trl->next;
        xrl->next = rl;
        rl = xrl;
        xrl = trl;
      } else
        xrl = trl;
    }

    xrl = tarl;
    while (tarl != NULL) {
      nrl = NULL;
      yrl = rl;
      while (yrl != NULL) {
        if (xrl->reg == yrl->reg) {
          if (nrl == NULL) {
            rl = yrl->next;
            mem_put(rmem, yrl);
          } else {
            nrl->next = yrl->next;
            mem_put(rmem, yrl);
          }
          tarl = tarl->next;
          mem_put(rmem, xrl);
          xrl = NULL;
          break;
        } else {
          nrl = yrl;
          yrl = yrl->next;
        }
      }
      if (xrl != NULL) {
        tarl = tarl->next;
        xrl->next = arl;
        arl = xrl;
        xrl = tarl;
      } else
        xrl = tarl;
    }

    if (sv == svt)
      traveling = 0;
    else {
      sv = next_subvol(&outside, &delta, sv, world->x_fineparts,
                       world->y_fineparts, world->z_fineparts, world->ny_parts,
                       world->nz_parts);
      delta.x = loc->x - outside.x;
      delta.y = loc->y - outside.y;
      delta.z = loc->z - outside.z;

      if (sv == NULL) {
        if ((delta.x * delta.x + delta.y * delta.y + delta.z * delta.z) <
            EPS_C * EPS_C) {
          mcell_log("Didn't quite reach waypoint target, fudging.");
          traveling = 0;
        } else {
          mcell_log("Couldn't reach waypoint target.");
          sv = find_subvolume(world, &outside, NULL);
        }
      }
    }
  }

  *rlp = rl;
  *arlp = arl;

  return 0;
}

/*************************************************************************
place_waypoints:
   In: world: simulation state
   Out: Returns 1 if malloc fails, 0 otherwise.
        Allocates waypoints to SSVs, if any are needed.
   Note: you must have initialized SSVs before calling this routine!
*************************************************************************/
int place_waypoints(struct volume *world) {
  int waypoint_in_wall = 0;

/* Being exactly in the center of a subdivision can be bad. */
/* Define "almost center" positions for X, Y, Z */
#define W_Xa (0.5 + 0.0005 * MY_PI)
#define W_Ya (0.5 + 0.0002 * MY_PI *MY_PI)
#define W_Za (0.5 - 0.00007 * MY_PI *MY_PI *MY_PI)
#define W_Xb (1.0 - W_Xa)
#define W_Yb (1.0 - W_Ya)
#define W_Zb (1.0 - W_Za)

  /* Probably ought to check for whether you really need waypoints */

  if (world->waypoints != NULL)
    free(world->waypoints);
  world->n_waypoints = world->n_subvols;
  world->waypoints =
      CHECKED_MALLOC_ARRAY(struct waypoint, world->n_waypoints, "waypoints");
  memset(world->waypoints, 0, world->n_waypoints * sizeof(struct waypoint));

  for (int px = 0; px < world->nx_parts - 1; px++) {
    for (int py = 0; py < world->ny_parts - 1; py++) {
      for (int pz = 0; pz < world->nz_parts - 1; pz++) {
        const int this_sv =
            pz + (world->nz_parts - 1) * (py + (world->ny_parts - 1) * px);
        struct waypoint *wp = &(world->waypoints[this_sv]);

        struct subvolume *sv = &(world->subvol[this_sv]);

        /* Place waypoint near center of subvolume (W_#a=W_#b=0.5 gives center)
         */
        wp->loc.x = W_Xa * world->x_fineparts[sv->llf.x] +
                    W_Xb * world->x_fineparts[sv->urb.x];
        wp->loc.y = W_Ya * world->y_fineparts[sv->llf.y] +
                    W_Yb * world->y_fineparts[sv->urb.y];
        wp->loc.z = W_Za * world->z_fineparts[sv->llf.z] +
                    W_Zb * world->z_fineparts[sv->urb.z];

        do {
          waypoint_in_wall = 0;
          struct wall_list *wl;
          for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
            double d = dot_prod(&(wp->loc), &(wl->this_wall->normal));
            if (eps_equals(d, wl->this_wall->d)) {
              waypoint_in_wall++;
              d = EPS_C * (double)((rng_uint(world->rng) & 0xF) - 8);
              if (!distinguishable(d, 0, EPS_C))
                d = 8 * EPS_C;
              wp->loc.x += d * wl->this_wall->normal.x;
              wp->loc.y += d * wl->this_wall->normal.y;
              wp->loc.z += d * wl->this_wall->normal.z;
              break;
            }
          }
        } while (waypoint_in_wall);

        if (pz > 0) {
          if (world->waypoints[this_sv - 1].regions != NULL) {
            wp->regions = dup_region_list(world->waypoints[this_sv - 1].regions,
                                          sv->local_storage->regl);
            if (wp->regions == NULL)
              return 1;
          } else
            wp->regions = NULL;

          if (world->waypoints[this_sv - 1].antiregions != NULL) {
            wp->antiregions =
                dup_region_list(world->waypoints[this_sv - 1].antiregions,
                                sv->local_storage->regl);
            if (wp->antiregions == NULL)
              return 1;
          } else
            wp->antiregions = NULL;

          if (find_enclosing_regions(
                  world, &(wp->loc), &(world->waypoints[this_sv - 1].loc),
                  &(wp->regions), &(wp->antiregions), sv->local_storage->regl))
            return 1;
        } else {
          wp->regions = NULL;
          wp->antiregions = NULL;
          if (find_enclosing_regions(world, &(wp->loc), NULL, &(wp->regions),
                                     &(wp->antiregions),
                                     sv->local_storage->regl))
            return 1;
        }
      }
    }
  }

  return 0;
#undef W_Zb
#undef W_Yb
#undef W_Xb
#undef W_Za
#undef W_Ya
#undef W_Xa
}

/******************************************************************
prepare_counters:
  In: world: simulation state
  Out: 0 if counter statements are correct, 1 otherwise.
  Note: A statement is incorrect if a non-closed manifold region
        tries to count a freely diffusing molecule.  Fixes up all
        count requests to point at the data we care about.
********************************************************************/
int prepare_counters(struct volume *world) {
  /* First give everything a sensible name, if needed */
  for (struct output_block *block = world->output_block_head; block != NULL;
       block = block->next) {
    for (struct output_set *set = block->data_set_head; set != NULL;
         set = set->next) {
      if (set->header_comment == NULL)
        continue;
      for (struct output_column *column = set->column_head; column != NULL;
           column = column->next) {
        if (column->expr->title == NULL)
          column->expr->title = oexpr_title(column->expr);
        if (column->expr->title == NULL)
          mcell_allocfailed("Unable to create title for reaction data output.");
      }
    }
  }

  /* Then convert all requests to real counts */
  for (struct output_request *request = world->output_request_head;
       request != NULL; request = request->next) {
    /* check whether the "count_location" refers to the instantiated
       object or region */
    if (request->count_location != NULL) {
      char *name = request->count_location->name;
      if (!((is_object_instantiated(request->count_location, world->root_instance)) ||
          ((world->dg_parse != NULL ) &&
          ((retrieve_sym(name, world->dg_parse->reg_sym_table) != NULL) ||
          (retrieve_sym(name, world->dg_parse->obj_sym_table) != NULL)))))

        mcell_error("The object/region name '%s' in the COUNT/TRIGGER "
                    "statement is not fully referenced.\n"
                    "  This occurs when a count is requested on an object "
                    "which has not been referenced\n"
                    "  (directly or indirectly) from an INSTANTIATE block in "
                    "the MDL file.", name);
    }

    if (request->count_target->sym_type == MOL) {
      struct species *sp = (struct species *)(request->count_target->value);

      /* For volume molecules: */
      if ((sp->flags & ON_GRID) == 0) {
        /* Make sure orientation is not set */
        if (request->count_orientation != ORIENT_NOT_SET) {
          switch (world->notify->useless_vol_orient) {
          case WARN_WARN:
            mcell_warn("An orientation has been given for the molecule '%s', "
                       "which is a volume molecule.\n"
                       "  Orientation is valid only for surface molecules, and "
                       "will be ignored.",
                       request->count_target->name);
          /* Fall through */
          case WARN_COPE:
            request->count_orientation = ORIENT_NOT_SET;
            break;

          case WARN_ERROR:
            mcell_error("An orientation has been given for the molecule '%s', "
                        "which is a volume molecule.\n"
                        "  Orientation is valid only for surface molecules.",
                        request->count_target->name);
            /*break;*/
          }
        }
      }
    }

    if (request->count_location != NULL &&
        request->count_location->sym_type == OBJ) {
      if (expand_object_output(request,
                               (struct object *)request->count_location->value,
                               world->reg_sym_table))
        mcell_error("Failed to expand request to count on object.");
    }

    if (instantiate_count_request(
        world->dynamic_geometry_flag, request, world->count_hashmask,
        world->count_hash, world->trig_request_mem, &world->elapsed_time,
        world->counter_mem)) {
      mcell_error("Failed to instantiate count request.");
    }
  }

  return 0;
}

/******************************************************************
is_object_instantiated:
  In: entry: symbol table entry to check
      root_instance: symbol_table entry against which the object is tested
  Out: 1 if the name of the object or one of its descendants matches the name
       of the symbol passed, 0 otherwise.
  Note: Checking is performed for all instantiated objects
********************************************************************/
int is_object_instantiated(struct sym_entry *entry,
                           struct object *root_instance) {
  struct object *obj = NULL;
  if (entry->sym_type == REG) {
    struct region *reg = entry->value;
    obj = ((struct region *)(entry->value))->parent;
    if (region_listed(obj->regions, reg)) {
      return 1; 
    }
    else {
      return 0; 
    }
  }
  else if (entry->sym_type == OBJ && entry->count != 0)
    obj = ((struct object *)(entry->value));
  else
    return 0;

  for (; obj != NULL; obj = obj->parent) {
    if (obj == root_instance)
      return 1;
  }

  return 0;
}

/*************************************************************************
check_counter_geometry:
   In: count_hashmask:
       count_hash:
       place_waypoints_flag:
   Out: 0 on success, 1 on failure.
        Checks all counters to make sure that if they are ENCLOSING,
        they count on closed regions.  If not, the function prints out
        the offending region name and returns 1.
*************************************************************************/
int check_counter_geometry(int count_hashmask, struct counter **count_hash,
                           byte *place_waypoints_flag) {
  /* Check to make sure what we've created is geometrically sensible */
  for (int i = 0; i < count_hashmask + 1; i++) {
    for (struct counter *cp = count_hash[i]; cp != NULL; cp = cp->next) {
      if ((cp->counter_type & ENCLOSING_COUNTER) != 0) {
        struct region *rp = cp->reg_type;

        if (rp->manifold_flag == MANIFOLD_UNCHECKED) {
          int count_regions_flag = 1;
          if (is_manifold(rp, count_regions_flag))
            rp->manifold_flag = IS_MANIFOLD;
          else
            rp->manifold_flag = NOT_MANIFOLD;
        }

        if (rp->manifold_flag == NOT_MANIFOLD)
          mcell_error("Cannot count molecules or events inside non-manifold "
                      "object region '%s'.  Please make sure that all "
                      "objects/regions used to count 3D molecules are "
                      "closed/watertight.",
                      rp->sym->name);

        (*place_waypoints_flag) = 1;
      }
    }
  }

  return 0;
}

/*************************************************************************
expand_object_output:
   In: request: request for a count
       obj: object upon which the request is made.
   Out: 0 on success, 1 on failure (memory allocation only?).
        Request is split into a separate request for each BOX and POLY
        object's ALL region that is a child of this object.  The result
        is then added up here.
   Note: This is probably broken for concentration.  It may also not be
         the most intuitive interpretation when used inside a large
         object with multiple layers of nesting--if one molecule is
         inside three sub-objects, it will be counted three times!
   PostNote: Checks that COUNT/TRIGGER statements are not allowed for
             metaobjects and release objects.
*************************************************************************/
int expand_object_output(
    struct output_request *request,
    struct object *obj,
    struct sym_table_head *reg_sym_table) {
#ifdef ALLOW_COUNTS_ON_METAOBJECT
  int n_expanded;
#endif

  switch (obj->object_type) {
  case REL_SITE_OBJ:
    mcell_error("COUNT and TRIGGER statements on release object '%s' are not "
                "allowed.\n",
                obj->sym->name);
    /*break;*/

  case META_OBJ:
/* XXX: Should this really be disabled?  Some comments by Tom lead me to
 * believe that, despite the potential confusion for users, this should
 * not be disabled.
 */
#ifndef ALLOW_COUNTS_ON_METAOBJECT
    mcell_error(
        "COUNT and TRIGGER statements on metaobject '%s' are not allowed.\n",
        obj->sym->name);
#else
#error "Support for counting in/on a metaobject doesn't work right now."
    n_expanded = 0;
    for (struct object *child = obj->first_child; child != NULL;
         child = child->next) {
      if (!object_has_geometry(child))
        continue; /* NOTE -- for objects nested N deep, we check this
                     N+(N-1)+...+2+1 times (slow) */
      if (n_expanded > 0) {
        struct output_request *new_request =
            (struct output_request *)mem_get(world->outp_request_mem);
        struct output_expression *oe = request->requester;
        struct output_expression *oel = new_output_expr(world->oexpr_mem);
        struct output_expression *oer = new_output_expr(world->oexpr_mem);
        if (new_request == NULL || oel == NULL || oer == NULL)
          mcell_allocfailed("Failed to expand count expression on object %s.",
                            obj->sym->name);
        oel->column = oer->column = oe->column;
        oel->expr_flags = oer->expr_flags = oe->expr_flags;
        oel->up = oer->up = oe;
        oel->left = request;
        oer->left = new_request;
        oel->oper = oer->oper = '#';
        oe->expr_flags = (oe->expr_flags & OEXPR_TYPE_MASK) | OEXPR_LEFT_OEXPR |
                         OEXPR_RIGHT_OEXPR;
        oe->left = oel;
        oe->right = oer;
        oe->oper = '+';

        new_request->report_type = request->report_type;
        new_request->count_target = request->count_target;
        new_request->requester = oer;
        request->requester = oel;
        new_request->next = request->next;
        request->next = new_request;
        request = new_request;
      }
      if (expand_object_output(request, child, reg_sym_table))
        return 1;
      ++n_expanded;
    }
    if (n_expanded == 0)
      mcell_error("Trying to count on object %s but it has no geometry.",
                  obj->sym->name);
#endif
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    request->count_location = NULL;
    for (struct region_list *rl = obj->regions; rl != NULL; rl = rl->next) {
      if (is_reverse_abbrev(",ALL", rl->reg->sym->name)) {
        request->count_location = rl->reg->sym;
        break;
      }
    }
    char *region_name = CHECKED_SPRINTF("%s,ALL", obj->sym->name);
    if (request->count_location == NULL)
      request->count_location = retrieve_sym(region_name, reg_sym_table);
    free(region_name);
    if (request->count_location == NULL)
      mcell_internal_error("ALL region missing on object %s", obj->sym->name);
    break;

  default:
    UNHANDLED_CASE(obj->object_type);
    /*return 1;*/
  }
  return 0;
}

/*************************************************************************
object_has_geometry:
   In: obj: object (instantiated in world)
   Out: 0 if there are no geometrical objects within that object (and it
        is not a geometrical object itself).  1 if there are such object.
*************************************************************************/
int object_has_geometry(struct object *obj) {
  struct object *child;
  switch (obj->object_type) {
  case BOX_OBJ:
  case POLY_OBJ:
    return 1;
    /*break;*/

  case META_OBJ:
    for (child = obj->first_child; child != NULL; child = child->next) {
      if (object_has_geometry(child))
        return 1;
    }
    break;

  case REL_SITE_OBJ:
  case VOXEL_OBJ:
  default:
    return 0;
    /*break;*/
  }
  return 0;
}

/*************************************************************************
instantiate_count_request:
   In: request: request for a count
       count_hashmask:
       count_hash:
       trig_request_mem:
       elapsed_time:
       counter_mem:
   Out: 0 on success, 1 on failure (memory allocation only?).
        Requesting output tree gets appropriate node pointed to the
        memory location where we will be collecting data.
*************************************************************************/
static int instantiate_count_request(
  int dyn_geom_flag, struct output_request *request, int count_hashmask,
  struct counter **count_hash, struct mem_helper *trig_request_mem,
  double *elapsed_time, struct mem_helper *counter_mem) {

  int request_hash = 0;
  struct rxn_pathname *rxpn_to_count;
  struct rxn *rx_to_count = NULL;
  struct species *mol_to_count = NULL;
  struct region *reg_of_count;
  struct counter *count = NULL;
  struct trigger_request *trig_req;
  u_int report_type_only;
  byte count_type;
  int is_enclosed;

  /* Set up and figure out hash value */
  void *to_count = request->count_target->value;
  switch (request->count_target->sym_type) {
  case MOL:
    rxpn_to_count = NULL;
    rx_to_count = NULL;
    mol_to_count = (struct species *)to_count;
    if ((mol_to_count->flags & NOT_FREE) == 0 &&
        (request->report_type & REPORT_TYPE_MASK) == REPORT_CONTENTS) {
      request->report_type |= REPORT_ENCLOSED;
    }
    request_hash = mol_to_count->hashval;
    break;

  case RXPN:
    rxpn_to_count = (struct rxn_pathname *)to_count;
    rx_to_count = rxpn_to_count->rx;
    mol_to_count = NULL;
    if ((rx_to_count->players[0]->flags & NOT_FREE) == 0 &&
        (rx_to_count->n_reactants == 1 ||
         (rx_to_count->players[1]->flags & NOT_FREE) == 0)) {
      request->report_type |= REPORT_ENCLOSED;
    }
    request_hash = rxpn_to_count->hashval;
    break;

  default:
    UNHANDLED_CASE(request->count_target->sym_type);
    /*return 1;*/
  }

  if (request->count_location != NULL) {
    if (request->count_location->sym_type != REG)
      mcell_internal_error(
          "Non-region location symbol (type=%d) in count request.",
          request->count_location->sym_type);
    reg_of_count = (struct region *)request->count_location->value;

    request_hash += reg_of_count->hashval;

  } else
    reg_of_count = NULL;
  request_hash &= count_hashmask;

  /* Now create count structs and set output expression to point to data */
  report_type_only = request->report_type & REPORT_TYPE_MASK;
  if (!dyn_geom_flag) {
    request->requester->expr_flags &= ~OEXPR_LEFT_REQUEST;
  }
  if ((request->report_type & REPORT_TRIGGER) == 0 &&
      request->count_location == NULL) /* World count is easy! */
  {
    request->report_type &= ~REPORT_ENCLOSED;
    switch (report_type_only) {
    case REPORT_CONTENTS:
      request->requester->expr_flags |= OEXPR_LEFT_INT;
      request->requester->left = (void *)&(mol_to_count->population);
      break;
    case REPORT_RXNS:
      assert(rx_to_count != NULL);
      request->requester->expr_flags |= OEXPR_LEFT_DBL;
      request->requester->left =
          (void *)&(rx_to_count->info[rxpn_to_count->path_num].count);
      break;
    default:
      mcell_internal_error("Invalid report type 0x%x in count request.",
                           report_type_only);
      /*return 1;*/
    }
  } else /* Triggered count or count on region */
  {
    /* Set count type flags */
    if (report_type_only == REPORT_RXNS)
      count_type = RXN_COUNTER;
    else
      count_type = MOL_COUNTER;
    if (request->report_type & REPORT_ENCLOSED) {
      assert(reg_of_count != NULL);
      reg_of_count->flags |= COUNT_ENCLOSED;
      count_type |= ENCLOSING_COUNTER;
      if (mol_to_count != NULL)
        mol_to_count->flags |= COUNT_ENCLOSED;
    }
    if (request->report_type & REPORT_TRIGGER) {
      assert(reg_of_count != NULL);
      count_type |= TRIG_COUNTER;
      reg_of_count->flags |= COUNT_TRIGGER;
    }

    for (count = count_hash[request_hash]; count != NULL; count = count->next) {
      if (count->reg_type == reg_of_count && count->target == to_count &&
          count_type == count->counter_type &&
          count->orientation == request->count_orientation &&
          periodic_boxes_are_identical(count->periodic_box, request->periodic_box))
        break;
    }
    if (count == NULL) {
      count = create_new_counter(reg_of_count, request->count_target->value,
                                 count_type, request->periodic_box, counter_mem);
      if (request->count_orientation != ORIENT_NOT_SET) {
        count->orientation = request->count_orientation;
      }

      count->next = count_hash[request_hash];
      count_hash[request_hash] = count;
    }

    is_enclosed = ((request->report_type & REPORT_ENCLOSED) != 0);

    /* set periodic box */
    count->periodic_box = request->periodic_box;

    /* Point appropriately */
    if (request->report_type & REPORT_TRIGGER) {
      trig_req = (struct trigger_request *)CHECKED_MEM_GET(
          trig_request_mem, "trigger notification request");
      trig_req->next = count->data.trig.listeners;
      count->data.trig.listeners = trig_req;
      trig_req->ear = request;

      request->requester->expr_flags |= OEXPR_TYPE_TRIG;

      if (mol_to_count != NULL)
        mol_to_count->flags |= COUNT_TRIGGER;
      assert(reg_of_count != NULL);
      switch (report_type_only) {
      case REPORT_CONTENTS:
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_CONTENTS;
        reg_of_count->flags |= COUNT_CONTENTS;
        break;
      case REPORT_RXNS:
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_RXNS;
        reg_of_count->flags |= COUNT_RXNS;
        break;
      case REPORT_FRONT_HITS:
      case REPORT_BACK_HITS:
      case REPORT_FRONT_CROSSINGS:
      case REPORT_BACK_CROSSINGS:
      case REPORT_ALL_HITS:
      case REPORT_ALL_CROSSINGS:
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_HITS;
        reg_of_count->flags |= COUNT_HITS;
        break;
      case REPORT_CONCENTRATION:
        if (mol_to_count != NULL) {
          if (mol_to_count->flags & ON_GRID) {
            mcell_error("ESTIMATE_CONC counting on regions is implemented only "
                        "for volume molecules, while %s is a surface molecule.",
                        mol_to_count->sym->name);
          } else {
            mol_to_count->flags |= COUNT_HITS;
            reg_of_count->flags |= COUNT_HITS;
          }
        }
        break;
      default:
        UNHANDLED_CASE(report_type_only);
        /*return 1;*/
      }
    } else /* Not trigger--set up for regular count */
    {
      assert(reg_of_count != NULL);
      request->requester->expr_flags |= OEXPR_LEFT_DBL; /* Assume double */
      switch (report_type_only) {
      case REPORT_CONTENTS:
        request->requester->expr_flags -= OEXPR_LEFT_DBL;
        request->requester->expr_flags |= OEXPR_LEFT_INT;

        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_CONTENTS;
        reg_of_count->flags |= COUNT_CONTENTS;
        if (!is_enclosed)
          request->requester->left = (void *)&(count->data.move.n_at);
        else
          request->requester->left = (void *)&(count->data.move.n_enclosed);
        break;
      case REPORT_RXNS:
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_RXNS;
        reg_of_count->flags |= COUNT_RXNS;
        if (!is_enclosed)
          request->requester->left = (void *)&(count->data.rx.n_rxn_at);
        else
          request->requester->left = (void *)&(count->data.rx.n_rxn_enclosed);
        break;
      case REPORT_FRONT_HITS:
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_HITS;
        reg_of_count->flags |= COUNT_HITS;
        request->requester->left = (void *)&(count->data.move.front_hits);
        break;
      case REPORT_BACK_HITS:
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_HITS;
        reg_of_count->flags |= COUNT_HITS;
        request->requester->left = (void *)&(count->data.move.back_hits);
        break;
      case REPORT_FRONT_CROSSINGS:
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_HITS;
        reg_of_count->flags |= COUNT_HITS;
        request->requester->left = (void *)&(count->data.move.front_to_back);
        break;
      case REPORT_BACK_CROSSINGS:
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_HITS;
        reg_of_count->flags |= COUNT_HITS;
        request->requester->left = (void *)&(count->data.move.back_to_front);
        break;
      case REPORT_ALL_HITS:
        request->requester->expr_flags |= OEXPR_RIGHT_DBL;
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_HITS;
        reg_of_count->flags |= COUNT_HITS;
        request->requester->left = (void *)&(count->data.move.front_hits);
        request->requester->right = (void *)&(count->data.move.back_hits);
        break;
      case REPORT_ALL_CROSSINGS:
        request->requester->expr_flags |= OEXPR_RIGHT_DBL;
        reg_of_count->flags |= COUNT_HITS;
        if (mol_to_count != NULL)
          mol_to_count->flags |= COUNT_HITS;
        request->requester->left = (void *)&(count->data.move.front_to_back);
        request->requester->right = (void *)&(count->data.move.back_to_front);
        break;
      case REPORT_CONCENTRATION:
        request->requester->expr_flags |= OEXPR_RIGHT_DBL;
        if (mol_to_count != NULL) {
          if (mol_to_count->flags & ON_GRID) {
            mcell_error("ESTIMATE_CONC counting on regions is implemented only "
                        "for volume molecules, while %s is a surface molecule.",
                        mol_to_count->sym->name);
          } else {
            mol_to_count->flags |= COUNT_HITS;
          }
        }
        reg_of_count->flags |= COUNT_HITS;
        request->requester->left = (void *)&(count->data.move.scaled_hits);
        request->requester->right = (void *)(elapsed_time);
        request->requester->oper = '/';
        break;

      default:
        UNHANDLED_CASE(report_type_only);
        /*return 1;*/
      }
    }
  }

  return 0;
}

/*************************************************************************
create_new_counter:
   In: where: region upon which to count
       who: target we're going to count (species or rxn pathname)
       what: what to count (*_COUNTER flags)
       img: in what periodic image to count
       counter_mem:
   Out: Newly allocated counter initialized with the given region and
        target, or NULL if there is a memory allocation error.
   Note: memory is allocated from world->counter_mem using mem_get,
         not from the global heap using malloc.
*************************************************************************/
static struct counter *create_new_counter(struct region *where, void *who,
  byte what, struct periodic_image *img, struct mem_helper *counter_mem) {

  struct counter *c;

  c = (struct counter *)CHECKED_MEM_GET(counter_mem, "counter");
  c->next = NULL;
  c->reg_type = where;
  c->target = who;
  c->orientation = ORIENT_NOT_SET;
  c->counter_type = what;
  c->periodic_box = img;
  if (what & TRIG_COUNTER) {
    c->data.trig.t_event = 0.0;
    c->data.trig.loc.x = c->data.trig.loc.y = c->data.trig.loc.z = 0.0;
    c->data.trig.orient = SHRT_MIN;
    c->data.trig.listeners = NULL;
  } else if (what & RXN_COUNTER) {
    c->data.rx.n_rxn_at = c->data.rx.n_rxn_enclosed = 0.0;
  } else if (what & MOL_COUNTER) {
    c->data.move.n_at = c->data.move.n_enclosed = 0;
    c->data.move.front_hits = c->data.move.back_hits = 0.0;
    c->data.move.front_to_back = c->data.move.back_to_front = 0.0;
    c->data.move.scaled_hits = 0.0;
  }
  return c;
}

/*************************************************************************
clean_region_lists:
   Cleans the region and antiregion lists, annihilating any items which appear
   on both lists.

   In: my_sv: subvolume containing waypoint
       p_all_regs: pointer to receive list of regions
       p_all_antiregs: pointer to receive list of antiregions
   Out: None
*************************************************************************/
static void clean_region_lists(struct subvolume *my_sv,
                               struct region_list **p_all_regs,
                               struct region_list **p_all_antiregs) {
  if ((*p_all_regs)->next != NULL || (*p_all_antiregs)->next != NULL) {
    struct region_list pre_sentry, pre_antisentry;
    struct region_list *prl, *parl, *rl, *arl;

    /* Sort by memory address to make mutual annihilation faster */
    if ((*p_all_regs)->next != NULL)
      *p_all_regs =
          (struct region_list *)void_list_sort((struct void_list *)*p_all_regs);
    if ((*p_all_antiregs)->next != NULL)
      *p_all_antiregs = (struct region_list *)void_list_sort(
          (struct void_list *)*p_all_antiregs);

    /* Need previous entry to fix up list, so we'll make an imaginary one for
     * 1st list element */
    pre_sentry.next = *p_all_regs;
    pre_antisentry.next = *p_all_antiregs;
    prl = &pre_sentry;
    parl = &pre_antisentry;

    /* If we cross a region both ways, throw both out (once) */
    for (rl = *p_all_regs, arl = *p_all_antiregs; rl != NULL && arl != NULL;
         prl = rl, rl = rl->next, parl = arl, arl = arl->next) {
      if (rl->reg == arl->reg) /* Mutual annihilation */
      {
        prl->next = rl->next;
        parl->next = arl->next;
        mem_put(my_sv->local_storage->regl, rl);
        mem_put(my_sv->local_storage->regl, arl);
        rl = prl;
        arl = parl;
      }
    }
    *p_all_regs = pre_sentry.next;
    *p_all_antiregs = pre_antisentry.next;
  } else if ((*p_all_regs)->reg == (*p_all_antiregs)->reg) {
    /* Crossed one region both ways, toss them */
    mem_put(my_sv->local_storage->regl, *p_all_regs);
    mem_put(my_sv->local_storage->regl, *p_all_antiregs);
    *p_all_regs = NULL;
    *p_all_antiregs = NULL;
  }
}

/*
 * function updating the hit counts during diffusion of a 2d
 * molecule if the latter hits a counted on region border on
 * the target wall.
 *
 * in:
 * ----
 *
 * hd_head  : head to linked list of hit_data for target region
 * current  : wall we are currently on
 * target   : wall we are hitting and which is counted on
 * sm        : surface molecule which is diffusing
 * direction: direction in which we are hitting the region border
 *            (0: outside in, 1: inside out)
 * crossed  : indicates if we crossed the region border or not
 *            (0: did not cross, 1: crossed)
 *
 * out:
 * ----
 *
 * nothing
 *
 *
 * side effects:
 * -------------
 *
 * a new hit_data structure is created and appended to the linked
 * list hd_head
 *
 * */
void update_hit_data(struct hit_data **hd_head, struct wall *current,
                     struct wall *target, struct surface_molecule *sm,
                     struct vector2 boundary_pos, int direction, int crossed) {

  struct hit_data *hd;

  hd = CHECKED_MALLOC_STRUCT(struct hit_data, "hit_data");
  hd->count_regions = target->counting_regions;
  hd->direction = direction;
  hd->crossed = crossed;
  hd->orientation = sm->orient;
  uv2xyz(&boundary_pos, current, &(hd->loc));
  hd->t = sm->t;
  if (*hd_head == NULL) {
    hd->next = NULL;
    *hd_head = hd;
  } else {
    hd->next = *hd_head;
    *hd_head = hd;
  }
}


/* count_regions_list updates COUNTS and TRIGGERS for surface_molecule sm
 * for all regions in the provided region_list */
void count_region_list(
    struct volume *world,
    struct region_list *regions,
    struct surface_molecule *sm,
    struct vector3 *where,
    int count_hashmask,
    int inc,
    struct periodic_image *previous_box) {

  struct counter **count_hash = world->count_hash;
  for (struct region_list *rl = regions; rl != NULL; rl = rl->next) {
    int hash_bin = (sm->properties->hashval + rl->reg->hashval) & count_hashmask;
    for (struct counter *c = count_hash[hash_bin]; c != NULL; c = c->next) {
      if (c->target == sm->properties && c->reg_type == rl->reg &&
          (c->counter_type & ENCLOSING_COUNTER) == 0) {
        if (c->counter_type & TRIG_COUNTER) {
          c->data.trig.t_event = sm->t;
          c->data.trig.orient = sm->orient;
          fire_count_event(world, c, inc, where, REPORT_CONTENTS | REPORT_TRIGGER);
        } else if ((c->orientation == ORIENT_NOT_SET) ||
                   (c->orientation == sm->orient) || (c->orientation == 0)) {
          if ((inc == 1) && (periodic_boxes_are_identical(
              sm->periodic_box, c->periodic_box))) {
            c->data.move.n_at++;
          }
          else if ((inc == -1) && (previous_box != NULL) &&
                   (periodic_boxes_are_identical(previous_box, c->periodic_box))) {
            c->data.move.n_at--;
          }
        }
      }
    }
  }
}
