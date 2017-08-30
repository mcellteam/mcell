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

#include "config.h"

#include "volume_output.h"
#include "logging.h"
#include "mcell_structs.h"
#include "sched_util.h"
#include "vol_util.h"
#include "strfunc.h"
#include "util.h"

#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

static int produce_item_header(FILE *out_file, struct volume_output_item *vo);

static int produce_mol_counts(struct volume *wrld, FILE *out_file,
                              struct volume_output_item *vo);

static int find_species_in_array(struct species **mols, int num_mols,
                                 struct species *ptr);

static int reschedule_volume_output_item(struct volume *wrld,
                                         struct volume_output_item *vo);

/*
 * Output a block of volume data as requested by the 'vo' object.
 */
int update_volume_output(struct volume *wrld, struct volume_output_item *vo) {
  int failure = 0;
  char *filename;

  switch (wrld->notify->volume_output_report) {
  case NOTIFY_NONE:
    break;

  case NOTIFY_BRIEF:
  case NOTIFY_FULL:
    mcell_log("Updating volume output '%s' scheduled at time %.15g on "
              "iteration %lld.",
              vo->filename_prefix, vo->t, wrld->current_iterations);
    break;

  default:
    UNHANDLED_CASE(wrld->notify->volume_output_report);
  }

  /* build the filename */
  filename = CHECKED_SPRINTF("%s.%lld.dat", vo->filename_prefix, wrld->current_iterations);

  /* Try to make the directory if it doesn't exist */
  if (make_parent_dir(filename)) {
    free(filename);
    return 1;
  }

  /* Output the volume item */
  failure = output_volume_output_item(wrld, filename, vo);
  free(filename);

  /* Reschedule this volume item, if appropriate */
  if (!failure)
    failure = reschedule_volume_output_item(wrld, vo);

  /* Should we return failure if we can't create the file?  Doing so will bring
   * down the entire sim...
   */
  return failure;
}

/*
 * Produce the output for a volume item.
 */
int output_volume_output_item(struct volume *wrld, char const *filename,
                              struct volume_output_item *vo) {
  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    mcell_perror_nodie(errno, "Couldn't open volume output file '%s'.",
                       filename);
    return 1;
  }

  if (produce_item_header(f, vo))
    goto failure;

  if (produce_mol_counts(wrld, f, vo))
    goto failure;

  fclose(f);
  return 0;

failure:
  fclose(f);
  return 1;
}

/*
 * Write the molecule counts to the file.
 *
 * XXX: Update this code to be smarter about the tradeoff between large slab
 * size (= large mem usage) and small slab size (= worse cache usage).
 */
static int produce_mol_counts(struct volume *wrld, FILE *out_file,
                              struct volume_output_item *vo) {
  struct volume_molecule *curmol;
  int *counters, *countersptr;
  int k, u, v;
  double z = vo->location.z, y = vo->location.y, x = vo->location.x;
  double x_lim = x + vo->voxel_size.x * (double)vo->nvoxels_x;
  double y_lim = y + vo->voxel_size.y * (double)vo->nvoxels_y;
  struct subvolume *cur_partition_z;
  double z_lim_part;

  /* Allocate memory for counters. */
  counters =
      CHECKED_MALLOC_ARRAY(int, vo->nvoxels_x * vo->nvoxels_y, "voxel slab");

  cur_partition_z = find_subvolume(wrld, &vo->location, NULL);
  if (cur_partition_z == NULL) {
    free(counters);
    mcell_internal_error(
        "While counting at [%g, %g, %g]: point isn't within a partition.", x, y,
        z);
    /*return 1;*/
  }

  z_lim_part = wrld->z_fineparts[cur_partition_z->urb.z];

  /* For each slab: */
  double r_voxsz_x = 1.0 / vo->voxel_size.x;
  double r_voxsz_y = 1.0 / vo->voxel_size.y;
  for (k = 0; k < vo->nvoxels_z; ++k) {
    double z_lim_slab = z + vo->voxel_size.z;
    struct subvolume *cur_partition_y;

    /* reset counters for this slab */
    memset(counters, 0, sizeof(int) * vo->nvoxels_x * vo->nvoxels_y);

  /* Loop over relevant partitions */
  keep_counting:
    cur_partition_y = cur_partition_z;
    while (cur_partition_y != NULL &&
           wrld->y_fineparts[cur_partition_y->llf.y] < y_lim) {
      struct subvolume *cur_partition = cur_partition_y;
      while (cur_partition != NULL &&
             wrld->x_fineparts[cur_partition->llf.x] < x_lim) {
        /* Count molecules */
        struct per_species_list *psl;
        int check_nonreacting = 0;
        int i = 0;
        for (i = 0; i < vo->num_molecules; ++i) {
          if (!(vo->molecules[i]->flags & CAN_VOLVOL)) {
            check_nonreacting = 1;
            break;
          }
        }
        for (psl = cur_partition->species_head; psl != NULL; psl = psl->next) {
          if (psl->properties == NULL) {
            if (!check_nonreacting)
              continue;
            else {
              for (curmol = psl->head; curmol != NULL;
                   curmol = curmol->next_v) {
                /* See if we're interested in this molecule */
                if (vo->num_molecules == 1) {
                  if (*vo->molecules != curmol->properties)
                    continue;
                } else {
                  if ((find_species_in_array(vo->molecules, vo->num_molecules,
                                             curmol->properties)) == -1)
                    continue;
                }

                /* Skip molecules not in our slab */
                if (curmol->pos.z < z || curmol->pos.z >= z_lim_slab)
                  continue;

                /* Skip molecules outside our domain */
                if (curmol->pos.x < x || curmol->pos.x >= x_lim ||
                    curmol->pos.y < y || curmol->pos.y >= y_lim)
                  continue;

                /* We've got a winner!  Add one to the appropriate voxel. */
                ++counters[((int)floor((curmol->pos.y - y) * r_voxsz_y)) *
                               vo->nvoxels_x +
                           (int)floor((curmol->pos.x - x) * r_voxsz_x)];
              }
            }
          } else {
            /* See if we're interested in this molecule */
            if (vo->num_molecules == 1) {
              if (*vo->molecules != psl->properties)
                continue;
            } else {
              if ((find_species_in_array(vo->molecules, vo->num_molecules,
                                         psl->properties)) == -1)
                continue;
            }

            for (curmol = psl->head; curmol != NULL; curmol = curmol->next_v) {
              /* Skip molecules not in our slab */
              if (curmol->pos.z < z || curmol->pos.z >= z_lim_slab)
                continue;

              /* Skip molecules outside our domain */
              if (curmol->pos.x < x || curmol->pos.x >= x_lim ||
                  curmol->pos.y < y || curmol->pos.y >= y_lim)
                continue;

              /* We've got a winner!  Add one to the appropriate voxel. */
              ++counters[((int)floor((curmol->pos.y - y) * r_voxsz_y)) *
                             vo->nvoxels_x +
                         (int)floor((curmol->pos.x - x) * r_voxsz_x)];
            }
          }
        }

        /* Advance to next x-partition */
        cur_partition =
            traverse_subvol(cur_partition, X_POS, wrld->ny_parts,
                            wrld->nz_parts);
      }

      /* Advance to next y-partition */
      cur_partition_y =
          traverse_subvol(cur_partition_y, Y_POS, wrld->ny_parts,
                          wrld->nz_parts);
    }

    /* If the slab crosses a Z boundary, keep on truckin' */
    if (z_lim_slab >= z_lim_part) {
      /* If we can get to the next partition, don't update slab and don't
       * spill!
       */
      cur_partition_z =
          traverse_subvol(cur_partition_z, Z_POS, wrld->ny_parts,
                          wrld->nz_parts);

      if (cur_partition_z != NULL) {
        z_lim_part = wrld->z_fineparts[cur_partition_z->urb.z];
        goto keep_counting;
      }
    } else {
      z = z_lim_slab;
      /* z_lim_slab will be updated on the next pass through the loop. */
    }

    /* Spill our counts */
    countersptr = counters;
    for (u = 0; u < vo->nvoxels_y; ++u) {
      for (v = 0; v < vo->nvoxels_x; ++v)
        fprintf(out_file, "%d ", *countersptr++);
      fprintf(out_file, "\n");
    }

    /* Extra newline to put visual separation between slabs */
    fprintf(out_file, "\n");
  }

  free(counters);
  return 0;
}

/*
 * Binary search for a pointer in an array of pointers.
 */
static int find_species_in_array(struct species **mols, int num_mols,
                                 struct species *ptr) {
  int lo = 0, hi = num_mols;
  while (hi - lo > 1) {
    int mid = (hi + lo) / 2;
    if (mols[mid] > ptr)
      hi = mid;
    else if (mols[mid] < ptr)
      lo = mid;
    else
      return mid;
  }

  if (mols[lo] == ptr)
    return lo;
  else
    return -1;
}

/*
 * Write the item header to the file.
 */
static int produce_item_header(FILE *out_file, struct volume_output_item *vo) {
  if (fprintf(out_file, "# nx=%d ny=%d nz=%d time=%g\n", vo->nvoxels_x,
              vo->nvoxels_y, vo->nvoxels_z, vo->t) < 0) {
    mcell_perror_nodie(errno, "Couldn't write header of volume output file.");
    return 1;
  }

  return 0;
}

/*
 * Reschedule a volume output item, if necessary.
 */
static int reschedule_volume_output_item(struct volume *wrld,
                                         struct volume_output_item *vo) {
  /* Find the next time */
  if (vo->timer_type == OUTPUT_BY_STEP)
    vo->t += (vo->step_time / wrld->time_unit);
  else {
    double time_scale = 0.0;

    /* Check if we're done */
    if (vo->next_time == vo->times + vo->num_times) {
      free(vo->filename_prefix);
      free(vo->molecules);
      free(vo->times);
      free(vo);
      return 0;
    }

    /* Compute the next time and advance the next_time ptr */
    if (vo->timer_type == OUTPUT_BY_ITERATION_LIST)
      time_scale = 1.0;
    else
      time_scale = 1.0 / wrld->time_unit;
    vo->t = (*vo->next_time++) * time_scale;
  }

  switch (wrld->notify->volume_output_report) {
  case NOTIFY_NONE:
  case NOTIFY_BRIEF:
    break;

  case NOTIFY_FULL:
    mcell_log("  Next output scheduled for time %.15g.",
              vo->t * wrld->time_unit);
    break;

  default:
    UNHANDLED_CASE(wrld->notify->volume_output_report);
  }

  /* Add to the schedule */
  if (schedule_add(wrld->volume_output_scheduler, vo))
    mcell_allocfailed("Failed to add volume output request to scheduler.");

  return 0;
}
