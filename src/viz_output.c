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

#include <stdarg.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

/*
#include "isaac64.h"
#include "rng.h"
*/

#include "logging.h"
#include "mcell_structs.h"
#include "grid_util.h"
#include "sched_util.h"
#include "viz_output.h"
#include "strfunc.h"
#include "util.h"
#include "vol_util.h"
#include "sym_table.h"

/* Output frame types. */
static int output_ascii_molecules(struct volume *world,
                                  struct viz_output_block *,
                                  struct frame_data_list *fdlp);

static int output_cellblender_molecules(struct volume *world,
                                        struct viz_output_block *,
                                        struct frame_data_list *fdlp);

/* == viz-specific Utilities == */

/*************************************************************************
frame_iteration:
    Gets the iteration number for a given time/iteration value and "type".

        In:  double iterval - the time/iteration value
             int type - the type of value
        Out: the frame time as an iteration number
**************************************************************************/
static long long frame_iteration(struct volume *world, double iterval,
                                 int type) {
  switch (type) {
  case OUTPUT_BY_ITERATION_LIST:
    return (long long)iterval;

  case OUTPUT_BY_TIME_LIST:
    if (world->chkpt_seq_num == 1) {
      return (long long)(iterval / world->time_unit + ROUND_UP);
    } else {
      if (iterval >= world->simulation_start_seconds) {
        return (long long) convert_seconds_to_iterations(
            world->start_iterations, world->time_unit,
            world->simulation_start_seconds, iterval) + ROUND_UP;
      } else {
        /* This iteration_time was in the past - just return flag.
           We do this because TIME_STEP may have been changed between
           checkpoints */
        return INT_MIN;
      }
    }

  default:
    mcell_internal_error("Invalid frame_data_list list_type (%d).", type);
    /*return -1;*/
  }
}

/*************************************************************************
sort_molecules_by_species:
    Scans over all molecules, sorting them into arrays by species.

        In:  struct abstract_molecule ****viz_molpp
             u_int  **viz_mol_countp
             int include_volume - should the lists include vol mols?
             int include_grid - should the lists include surface mols?
        Out: 0 on success, 1 on error; viz_molpp and viz_mol_countp arrays are
             allocated and filled with sorted data.
**************************************************************************/
static int sort_molecules_by_species(struct volume *world,
                                     struct viz_output_block *vizblk,
                                     struct abstract_molecule ****viz_molpp,
                                     u_int **viz_mol_countp, int include_volume,
                                     int include_grid) {
  struct storage_list *slp;
  u_int *counts;
  int species_index;

  /* XXX: May leave memory allocated on failure */
  if ((*viz_molpp = (struct abstract_molecule ***)allocate_ptr_array(
           world->n_species)) == NULL)
    return 1;
  if ((counts = *viz_mol_countp = allocate_uint_array(world->n_species, 0)) ==
      NULL)
    return 1;

  /* Walk through the species allocating arrays for all molecules of that
   * species */
  for (species_index = 0; species_index < world->n_species; ++species_index) {
    int mol_count;
    u_int spec_id = world->species_list[species_index]->species_id;

    if (vizblk->species_viz_states[species_index] == EXCLUDE_OBJ)
      continue;

    if (world->species_list[species_index]->flags & IS_SURFACE)
      continue;

    if (!include_grid && (world->species_list[species_index]->flags & ON_GRID))
      continue;

    if (!include_volume &&
        !(world->species_list[species_index]->flags & ON_GRID))
      continue;

    mol_count = world->species_list[species_index]->population;
    if (mol_count <= 0)
      continue;

    if (((*viz_molpp)[spec_id] =
             (struct abstract_molecule **)allocate_ptr_array(mol_count)) ==
        NULL)
      return 1;
  }

  /* Sort molecules by species id */
  for (slp = world->storage_head; slp != NULL; slp = slp->next) {
    struct storage *sp = slp->store;
    struct schedule_helper *shp;
    struct abstract_molecule *amp;
    int sched_slot_index;
    for (shp = sp->timer; shp != NULL; shp = shp->next_scale) {
      for (sched_slot_index = -1; sched_slot_index < shp->buf_len;
           ++sched_slot_index) {
        for (amp = (struct abstract_molecule *)((sched_slot_index < 0)
                                                    ? shp->current
                                                    : shp->circ_buf_head
                                                          [sched_slot_index]);
             amp != NULL; amp = amp->next) {
          u_int spec_id;
          if (amp->properties == NULL)
            continue;

          spec_id = amp->properties->species_id;
          if (vizblk->species_viz_states[spec_id] == EXCLUDE_OBJ)
            continue;

          if (!include_grid && (amp->flags & TYPE_MASK) != TYPE_VOL)
            continue;

          if (!include_volume && (amp->flags & TYPE_MASK) == TYPE_VOL)
            continue;

          if (counts[spec_id] < amp->properties->population)
            (*viz_molpp)[spec_id][counts[spec_id]++] = amp;
          else {
            mcell_warn("Molecule count disagreement!\n"
                       "  Species %s  population = %d  count = %d",
                       amp->properties->sym->name, amp->properties->population,
                       counts[spec_id]);
          }
        }
      }
    }
  }

  return 0;
}

/*************************************************************************
reset_time_values:
    Scan over all frame data elements, resetting the "next" iteration state to
    the soonest iteration following or equal to the specified iteration.

        In:  struct frame_data_list *fdlp - the head of the frame data list
             long long curiter - the minimum iteration number
        Out: 0 on success, 1 on failure
**************************************************************************/
static int reset_time_values(struct volume *world, struct frame_data_list *fdlp,
                             long long curiter) {
  /* If we've loaded a checkpoint, don't output on the first iter */
  if (curiter != 0)
    ++curiter;

  for (; fdlp != NULL; fdlp = fdlp->next) {
    fdlp->curr_viz_iteration = fdlp->iteration_list;
    fdlp->viz_iteration = -1;

    /* Scan for first iteration >= curiter */
    while (fdlp->curr_viz_iteration != NULL) {
      if (frame_iteration(world, fdlp->curr_viz_iteration->value,
                          fdlp->list_type) >= curiter)
        break;
      fdlp->curr_viz_iteration = fdlp->curr_viz_iteration->next;
    }

    /* If we had an iteration, use it to set viz_iteration */
    if (fdlp->curr_viz_iteration != NULL)
      fdlp->viz_iteration = frame_iteration(
          world, fdlp->curr_viz_iteration->value, fdlp->list_type);
  }

  return 0;
}

/*************************************************************************
count_time_values:
    Scan over all frame data elements, counting iterations.  The number of
    distinct iterations on which viz data will be output is returned.  In the
    process, each frame data elements n_viz_iterations is set, as is the "final
    viz iteration".  The total number of distinct iterations is returned.

        In:  struct frame_data_list *fdlp - the head of the frame data list
        Out: num distinct iterations
**************************************************************************/
static int count_time_values(struct volume *world,
                             struct frame_data_list *const fdlp) {
  int time_values = 0;
  long long curiter = -1;
  struct frame_data_list *fdlpcur = NULL;
  for (fdlpcur = fdlp; fdlpcur != NULL; fdlpcur = fdlpcur->next) {
    fdlpcur->curr_viz_iteration = fdlpcur->iteration_list;
    fdlpcur->n_viz_iterations = 0;
    fdlpcur->viz_iteration = -1;
  }

  while (1) {
    curiter = -1;

    /* Find the next iteration */
    for (fdlpcur = fdlp; fdlpcur != NULL; fdlpcur = fdlpcur->next) {
      long long thisiter;
      if (fdlpcur->curr_viz_iteration == NULL)
        continue;

      thisiter = frame_iteration(world, fdlpcur->curr_viz_iteration->value,
                                 fdlpcur->list_type);

      if (curiter == -1)
        curiter = thisiter;
      else if (thisiter < curiter)
        curiter = thisiter;
    }

    /* If we found no new iteration, quit */
    if (curiter == -1)
      break;

    /* We won't create any more output frames after the apocalypse. */
    if (curiter > world->iterations)
      break;

    /* We won't create any more output frames after we checkpoint. */
    if (world->chkpt_iterations != 0 &&
        curiter > world->start_iterations + world->chkpt_iterations)
      break;

    /* We found at least one more.  Note that the only time we will output at
     * iteration == start_iterations is when start_iterations is zero.  This is because we
     * do not output on the first iteration after we resume.
     */
    if (curiter > world->start_iterations)
      ++time_values;
    else if ((world->start_iterations | curiter) == 0)
      ++time_values;

    /* Advance any frame data items which are set to this iteration */
    for (fdlpcur = fdlp; fdlpcur != NULL; fdlpcur = fdlpcur->next) {
      if (fdlpcur->curr_viz_iteration == NULL)
        continue;

      if (curiter > world->start_iterations || (world->start_iterations | curiter) == 0) {
        if (frame_iteration(world, fdlpcur->curr_viz_iteration->value,
                            fdlpcur->list_type) == curiter)
          ++fdlpcur->n_viz_iterations;
      }

      while (fdlpcur->curr_viz_iteration &&
             frame_iteration(world, fdlpcur->curr_viz_iteration->value,
                             fdlpcur->list_type) == curiter)
        fdlpcur->curr_viz_iteration = fdlpcur->curr_viz_iteration->next;
    }
  }

  return time_values;
}

/************************************************************************
output_ascii_molecules:
In: vizblk: VIZ_OUTPUT block for this frame list
    a frame data list (internal viz output data structure)
Out: 0 on success, 1 on failure.  The positions of molecules are output
     in exponential floating point notation (with 8 decimal places)
*************************************************************************/
static int output_ascii_molecules(struct volume *world,
                                  struct viz_output_block *vizblk,
                                  struct frame_data_list *fdlp) {
  FILE *custom_file;
  char *cf_name;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  struct surface_molecule *gmp;
  short orient = 0;

  int ndigits, i;
  long long lli;

  struct vector3 where, norm;

  no_printf("Output in ASCII mode (molecules only)...\n");

  if ((fdlp->type == ALL_MOL_DATA) || (fdlp->type == MOL_POS)) {
    lli = 10;
    for (ndigits = 1; lli <= world->iterations && ndigits < 20;
         lli *= 10, ndigits++) {
    }
    cf_name =
        CHECKED_SPRINTF("%s.ascii.%.*lld.dat", vizblk->file_prefix_name,
                        ndigits, fdlp->viz_iteration);
    if (cf_name == NULL)
      return 1;
    if (make_parent_dir(cf_name)) {
      free(cf_name);
      mcell_error(
          "Failed to create parent directory for ASCII-mode VIZ output.");
      /*return 1;*/
    }
    custom_file = open_file(cf_name, "w");
    if (!custom_file)
      mcell_die();
    else {
      no_printf("Writing to file %s\n", cf_name);
    }
    free(cf_name);
    cf_name = NULL;

    for (slp = world->storage_head; slp != NULL; slp = slp->next) {
      for (shp = slp->store->timer; shp != NULL; shp = shp->next_scale) {
        for (i = -1; i < shp->buf_len; i++) {
          for (aep = (i < 0) ? shp->current : shp->circ_buf_head[i];
               aep != NULL; aep = aep->next) {
            amp = (struct abstract_molecule *)aep;
            if (amp->properties == NULL)
              continue;

            int id = vizblk->species_viz_states[amp->properties->species_id];
            if (id == EXCLUDE_OBJ)
              continue;

            if ((amp->properties->flags & NOT_FREE) == 0) {
              mp = (struct volume_molecule *)amp;
              where.x = mp->pos.x;
              where.y = mp->pos.y;
              where.z = mp->pos.z;
              norm.x = 0;
              norm.y = 0;
              norm.z = 0;
            } else if ((amp->properties->flags & ON_GRID) != 0) {
              gmp = (struct surface_molecule *)amp;
              uv2xyz(&(gmp->s_pos), gmp->grid->surface, &where);
              orient = gmp->orient;
              norm.x = orient * gmp->grid->surface->normal.x;
              norm.y = orient * gmp->grid->surface->normal.y;
              norm.z = orient * gmp->grid->surface->normal.z;
            } else
              continue;

            where.x *= world->length_unit;
            where.y *= world->length_unit;
            where.z *= world->length_unit;
            /*
                        fprintf(custom_file,"%d %15.8e %15.8e %15.8e
               %2d\n",id,where.x,where.y,where.z,orient);
            */
            if (id == INCLUDE_OBJ) {
              /* write name of molecule */
              fprintf(custom_file, "%s %lu %.9g %.9g %.9g %.9g %.9g %.9g\n",
                      amp->properties->sym->name, amp->id, where.x, where.y,
                      where.z, norm.x, norm.y, norm.z);
            } else {
              /* write state value of molecule */
              fprintf(custom_file, "%d %lu %.9g %.9g %.9g %.9g %.9g %.9g\n", id,
                      amp->id, where.x, where.y, where.z, norm.x, norm.y,
                      norm.z);
            }
          }
        }
      }
    }
    fclose(custom_file);
  }

  return 0;
}


typedef struct external_mol_viz_struct {
  char mol_type;  /* s = surface, v = volume, n = none (don't display?) */
  float pos_x, pos_y, pos_z;
  float norm_x, norm_y, norm_z;
  struct external_mol_viz_struct *next_mol;
} external_mol_viz;

typedef struct external_mol_viz_by_name_struct {
  char *mol_name;
  external_mol_viz *mol_list;
  struct external_mol_viz_by_name_struct *next_name;
} external_mol_viz_by_name;



typedef struct external_mol_viz_entry_struct {
  char *mol_id_string;  /* Typically the NAUTY string (although that will be the key) */
} external_mol_viz_entry;


static struct sym_table_head *graph_pattern_table = NULL;
static long next_molcomp_id = 0L;

typedef struct external_molcomp_loc_struct {
  bool is_mol;
  bool has_coords;
  bool is_final;
  double x, y, z;
  double kx, ky, kz;
  char *name;
  char *graph_string;
  int num_peers;
  int *peers;
} external_molcomp_loc;

typedef struct molcomp_list_struct {
  external_molcomp_loc *molcomp_array;
  int num_molcomp_items;
  long molcomp_id;
} molcomp_list;

static void dump_molcomp_array ( external_molcomp_loc *molcomp_array, int num_parts ) {
  int i, j;
  fprintf ( stdout, "%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%\n" );
  for (i=0; i<num_parts; i++) {
    fprintf ( stdout, "[%d] = %s (", i, molcomp_array[i].name );
    if (molcomp_array[i].is_mol) {
      fprintf ( stdout, "m" );
    } else {
      fprintf ( stdout, "c" );
    }
    fprintf ( stdout, ") at  (%g, %g, %g) with peers [", molcomp_array[i].x, molcomp_array[i].y, molcomp_array[i].z );
    for (j=0; j<molcomp_array[i].num_peers; j++) {
      fprintf ( stdout, "%d", molcomp_array[i].peers[j] );
      if (j < molcomp_array[i].num_peers - 1) {
        fprintf ( stdout, "," );
      }
    }
    fprintf ( stdout, "]\n" );
  }
  fprintf ( stdout, "%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%=%%\n" );
}

/*
static void dump_molcomp_list ( molcomp_list *mcl ) {
  dump_molcomp_array ( mcl->molcomp_array, mcl->num_molcomp_items );
}
*/

double clipped_cos ( double angle ) {
  double v = cos(angle);
  if ((v > -1e-6) && (v < 1e-6)) v = 0;
  return ( v );
}

double clipped_sin ( double angle ) {
  double v = sin(angle);
  if ((v > -1e-6) && (v < 1e-6)) v = 0;
  return ( v );
}

/*
static void set_molcomp_positions_2D ( external_molcomp_loc *molcomp_array, int num_parts ) {
  // Compute positions for all molecules/components in a molcomp_array
  //#### fprintf ( stdout, "Begin Building molcomp_positions for:\n" );
  //#### dump_molcomp_array ( molcomp_array, num_parts );

  double scale = 0.02;

  int i, j;

  // Start by initializing all molecules and components to defaults
  for (i=0; i<num_parts; i++) {
    molcomp_array[i].x = 0;
    molcomp_array[i].y = 0;
    molcomp_array[i].z = 0;
    molcomp_array[i].has_coords = 1;
    molcomp_array[i].is_final = 0;
  }

  // Next assign default locations for binding sites around each molecule
  for (i=0; i<num_parts; i++) {
    if (molcomp_array[i].is_mol) {
      // This is a molecule
      int num_peers = molcomp_array[i].num_peers;
      for (j=0; j<num_peers; j++) {
        double angle = j * 2 * MY_PI / num_peers;
        molcomp_array[molcomp_array[i].peers[j]].x = scale * clipped_cos(angle);
        molcomp_array[molcomp_array[i].peers[j]].y = scale * clipped_sin(angle);
        molcomp_array[molcomp_array[i].peers[j]].z = 0;
      }
    }
  }

  // Start by setting all molecules at the origin and assigning default binding site locations
  for (i=0; i<num_parts; i++) {
    if (molcomp_array[i].has_coords == 0) {
      molcomp_array[i].x = 0;
      molcomp_array[i].y = 0;
      molcomp_array[i].z = 0;
      molcomp_array[i].has_coords = 1;
    }
    molcomp_array[i].is_final = 0;
  }

  //#### fprintf ( stdout, "Finished assigning molecule-centric coordinates:\n" );
  //#### dump_molcomp_array ( molcomp_array, num_parts );


  //#### fprintf ( stdout, "Done Building molcomp_positions for:\n" );
  //#### dump_molcomp_array ( molcomp_array, num_parts );
}
*/

/*
static void set_component_positions_2D ( struct volume *world, external_molcomp_loc *mc, int num_parts ) {
  double scale = 0.02;
  int mi;
  for (mi=0; mi<num_parts; mi++) {
    if (mc[mi].is_mol) {
      //#### fprintf ( stdout, "Setting component positions for %s\n", mc[mi].name );
      for (int ci=0; ci<mc[mi].num_peers; ci++) {
        double angle = 2 * MY_PI * ci / mc[mi].num_peers;
        mc[mc[mi].peers[ci]].x = scale * cos(angle);
        mc[mc[mi].peers[ci]].y = scale * sin(angle);
        //#### fprintf ( stdout, "  Component %s is at (%g,%g)\n", mc[mc[mi].peers[ci]].name, mc[mc[mi].peers[ci]].x, mc[mc[mi].peers[ci]].y );
      }
    }
  }
}
*/

static void set_component_positions_by_table ( struct volume *world, external_molcomp_loc *mc, int num_parts ) {
  double scale = 0.02;
  int mi;
  //fprintf ( stdout, "Setting positions by table.\n" );
  // Dump the table just to verify how to read it
  //fprintf ( stdout, "==============================================\n" );
  //dump_symtab(world->mol_ss_sym_table);
  //fprintf ( stdout, "==============================================\n" );

  for (mi=0; mi<num_parts; mi++) {
    if (mc[mi].is_mol) {
      // Look up this molecule in the table
      //fprintf ( stdout, "  Looking up component positions for %s\n", mc[mi].name );
      struct sym_entry *sp;
      sp = retrieve_sym(mc[mi].name, world->mol_ss_sym_table);
      if (sp != NULL) {
        //fprintf ( stdout, "    Got an entry with symbol type %d and name %s\n", sp->sym_type, sp->name );
        // Fill in all of the component positions in this molecule from the list
        if (sp->sym_type == MOL_SS) {
          //fprintf ( stdout, "       It's a Spatially Structured Molecule!!\n" );

          // Set the has_coords flag to false on each component of this molecule
          for (int ci=0; ci<mc[mi].num_peers; ci++) {
            mc[mc[mi].peers[ci]].has_coords = false;
          }

          // Walk through the list of component positions and for each one, find a matching mol_comp entry without coordinates

          struct mol_ss *mol_ss_ptr = (struct mol_ss *)(sp->value);
          struct mol_comp_ss *mc_ptr = mol_ss_ptr->mol_comp_ss_head;
          int comp_count = 0;
          // char *translations[5] = { "COINCIDENT", "XYZ", "XYZA", "XYZRef", "XYZVA" };
          while (mc_ptr != NULL) {
            //fprintf ( stdout, "         Component %d is \"%s\" of type %s at (%g,%g,%g).\n", comp_count, mc_ptr->name, translations[mc_ptr->spatial_type], mc_ptr->loc_x, mc_ptr->loc_y, mc_ptr->loc_z );
            for (int ci=0; ci<mc[mi].num_peers; ci++) {
              if ( (!mc[mc[mi].peers[ci]].has_coords) && (strcmp(mc[mc[mi].peers[ci]].name, mc_ptr->name) == 0) ) {
                mc[mc[mi].peers[ci]].x = mc_ptr->loc_x;
                mc[mc[mi].peers[ci]].y = mc_ptr->loc_y;
                mc[mc[mi].peers[ci]].z = mc_ptr->loc_z;
                mc[mc[mi].peers[ci]].kx = mc_ptr->rot_axis_x;  // These are currently key locations rather than rotation axis locations
                mc[mc[mi].peers[ci]].ky = mc_ptr->rot_axis_y;  // These are currently key locations rather than rotation axis locations
                mc[mc[mi].peers[ci]].kz = mc_ptr->rot_axis_z;  // These are currently key locations rather than rotation axis locations
                mc[mc[mi].peers[ci]].has_coords = true;
                if (world->dump_level >= 20) {
                  fprintf ( stdout, "    Component %s is at (%g,%g,%g)\n", mc[mc[mi].peers[ci]].name, mc[mc[mi].peers[ci]].x, mc[mc[mi].peers[ci]].y, mc[mc[mi].peers[ci]].z );
                  fprintf ( stdout, "       Ref key for %s is at (%g,%g,%g)\n", mc[mc[mi].peers[ci]].name, mc[mc[mi].peers[ci]].kx, mc[mc[mi].peers[ci]].ky, mc[mc[mi].peers[ci]].kz );
                }
                break;
              }
            }
            comp_count += 1;
            mc_ptr = mc_ptr->next;
          }

          //fprintf ( stdout, "       This molecule has %d components.\n", comp_count );

          // Just to be sure, set any unset coordinates to 0

          for (int ci=0; ci<mc[mi].num_peers; ci++) {
            if (!mc[mc[mi].peers[ci]].has_coords) {
              mc[mc[mi].peers[ci]].x = 0.0;
              mc[mc[mi].peers[ci]].y = 0.0;
              mc[mc[mi].peers[ci]].z = 0.0;
              mc[mc[mi].peers[ci]].kx = 0.0;
              mc[mc[mi].peers[ci]].ky = 0.0;
              mc[mc[mi].peers[ci]].kz = 0.0;
              mc[mc[mi].peers[ci]].has_coords = true;
            }
          }

        }
      } else {
        fprintf ( stdout, "    No entry found for %s, using default. This is unexpected!!\n", mc[mi].name );
        for (int ci=0; ci<mc[mi].num_peers; ci++) {
          double angle = 2 * MY_PI * ci / mc[mi].num_peers;
          mc[mc[mi].peers[ci]].x = scale * cos(angle);
          mc[mc[mi].peers[ci]].y = scale * sin(angle);
          mc[mc[mi].peers[ci]].z = 0.0;
          mc[mc[mi].peers[ci]].kx = 0.0;
          mc[mc[mi].peers[ci]].ky = 0.0;
          mc[mc[mi].peers[ci]].kz = scale;
          //#### fprintf ( stdout, "  Component %s is at (%g,%g)\n", mc[mc[mi].peers[ci]].name, mc[mc[mi].peers[ci]].x, mc[mc[mi].peers[ci]].y );
        }
      }
    }
  }
}


static void bind_molecules_at_components ( struct volume *world, external_molcomp_loc *mc, int num_parts, int fixed_comp_index, int var_comp_index, bool as3D, bool with_rot ) {
  // Bind these two molecules by aligning their axes and shifting to align their components
  //#### fprintf ( stdout, "########## Binding %s to %s\n", mc[fixed_comp_index].name, mc[var_comp_index].name );
  //#### dump_molcomp_array(mc,num_parts);

  int fixed_mol_index = mc[fixed_comp_index].peers[0];
  int var_mol_index = mc[var_comp_index].peers[0];

  double fixed_vec[3];  // This will either hold 2 or 3 values, but since it's temporary, the allocation difference doesn't accumulate.
  double   var_vec[3];
  fixed_vec[0] = mc[fixed_comp_index].x - mc[fixed_mol_index].x;
  fixed_vec[1] = mc[fixed_comp_index].y - mc[fixed_mol_index].y;
  var_vec[0]   = mc[  var_comp_index].x - mc[  var_mol_index].x;
  var_vec[1]   = mc[  var_comp_index].y - mc[  var_mol_index].y;
  if (as3D) {
    fixed_vec[2] = mc[fixed_comp_index].z - mc[fixed_mol_index].z;
    var_vec[2]   = mc[  var_comp_index].z - mc[  var_mol_index].z;
  }

  double fixed_mag;
  double   var_mag;
  if (as3D) {
    fixed_mag = sqrt ( (fixed_vec[0]*fixed_vec[0]) + (fixed_vec[1]*fixed_vec[1]) + (fixed_vec[2]*fixed_vec[2]) );
      var_mag = sqrt ( (  var_vec[0]*  var_vec[0]) + (  var_vec[1]*  var_vec[1]) + (  var_vec[2]*  var_vec[2]) );
  } else {
    fixed_mag = sqrt ( (fixed_vec[0]*fixed_vec[0]) + (fixed_vec[1]*fixed_vec[1]) );
      var_mag = sqrt ( (  var_vec[0]*  var_vec[0]) + (  var_vec[1]*  var_vec[1]) );
  }

  double dot_prod;
  if (as3D) {
    dot_prod = (fixed_vec[0] * var_vec[0]) + (fixed_vec[1] * var_vec[1]) + (fixed_vec[2] * var_vec[2]);
  } else {
    dot_prod = (fixed_vec[0] * var_vec[0]) + (fixed_vec[1] * var_vec[1]);
  }

  // In general, the magnitudes should be checked for 0. However, in this case, they were generated as non-zero.
  double norm_dot_prod = dot_prod / ( fixed_mag * var_mag );

  // Ensure that the dot product is a legal argument for the "acos" function:
  if (norm_dot_prod >  1) { norm_dot_prod =  1; }
  if (norm_dot_prod < -1) { norm_dot_prod = -1; }

  //#### fprintf ( stdout, "norm_dot_prod = %g\n", norm_dot_prod );
  double angle = acos ( norm_dot_prod );
  //#### fprintf ( stdout, "Angle (from acos) = %g\n", angle );

  if (as3D) {
    // This seems to be required to get everything right:
    angle = -angle;
  } else {
    // Try using the cross product to fix the direction issue
    //#### fprintf ( stdout, "Cross of (%g,%g) X (%g,%g) = %g\n", fixed_vec[0], fixed_vec[1], var_vec[0], var_vec[1], (fixed_vec[0] * var_vec[1]) - (fixed_vec[1] * var_vec[0]) );
    if ( ( (fixed_vec[0] * var_vec[1]) - (fixed_vec[1] * var_vec[0]) ) > 0 ) {
      angle = -angle;
    }
  }

  // Reverse the direction since we want the components attached to each other
  angle = MY_PI + angle;

  // Normalize between -PI and PI
  while (angle > MY_PI) {
    angle = angle - (2 * MY_PI);
  }
  while (angle <= -MY_PI) {
    angle = angle + (2 * MY_PI);
  }

  //angle = -angle;

  //#### fprintf ( stdout, "Final corrected angle = %g\n", angle );

  //#### fprintf ( stdout, "Binding between f(%s: %g,%g) and v(%s: %g,%g) is at angle %g deg\n", mc[fixed_comp_index].name, fixed_vec[0], fixed_vec[1], mc[var_comp_index].name, var_vec[0], var_vec[1], 180*angle/MY_PI );

  // Rotate all of the components of the var_mol_index by the angle
  double cos_angle = cos(angle);
  double sin_angle = sin(angle);

  if (as3D) {
    double cross_prod[3];
    cross_prod[0] = (fixed_vec[1] * var_vec[2]) - (fixed_vec[2] * var_vec[1]);
    cross_prod[1] = (fixed_vec[2] * var_vec[0]) - (fixed_vec[0] * var_vec[2]);
    cross_prod[2] = (fixed_vec[0] * var_vec[1]) - (fixed_vec[1] * var_vec[0]);

    double xpx, xpy, xpz;
    xpx = cross_prod[0] / (fixed_mag * var_mag);
    xpy = cross_prod[1] / (fixed_mag * var_mag);
    xpz = cross_prod[2] / (fixed_mag * var_mag);

    double axis_length = sqrt ( (xpx*xpx) + (xpy*xpy) + (xpz*xpz) );

    double R[3][3] = { { 1, 0, 0 },
                       { 0, 1, 0 },
                       { 0, 0, 1 } };
    if (axis_length < 1e-30) {
      // Can't compute a meaningful unit vector for the rotation matrix ... make it identity
      if (norm_dot_prod < 0) {
        // R is fine as defined
      } else {
        // Change the sign on the diagonal of R
        for (int i=0; i<3; i++) {
          R[i][i] = -1;
        }
      }
    } else {
      // Build the rotation matrix R
      double ux = xpx / axis_length;
      double uy = xpy / axis_length;
      double uz = xpz / axis_length;
      double omca = 1 - cos_angle;

      R[0][0] = cos_angle + (ux*ux*omca);
      R[0][1] = (ux*uy*omca) - (uz*sin_angle);
      R[0][2] = (ux*uz*omca) + (uy*sin_angle);

      R[1][0] = (uy*ux*omca) + (uz*sin_angle);
      R[1][1] = cos_angle + (uy*uy*omca);
      R[1][2] = (uy*uz*omca) - (ux*sin_angle);
                
      R[2][0] = (uz*ux*omca) - (uy*sin_angle);
      R[2][1] = (uz*uy*omca) + (ux*sin_angle);
      R[2][2] = cos_angle + (uz*uz*omca);
    }

    for (int ci=0; ci<mc[var_mol_index].num_peers; ci++) {
      // Rotate the component locations
      double x = mc[mc[var_mol_index].peers[ci]].x;
      double y = mc[mc[var_mol_index].peers[ci]].y;
      double z = mc[mc[var_mol_index].peers[ci]].z;
      mc[mc[var_mol_index].peers[ci]].x = (R[0][0]*x) + (R[0][1]*y) + (R[0][2]*z);
      mc[mc[var_mol_index].peers[ci]].y = (R[1][0]*x) + (R[1][1]*y) + (R[1][2]*z);
      mc[mc[var_mol_index].peers[ci]].z = (R[2][0]*x) + (R[2][1]*y) + (R[2][2]*z);

      // Rotate the rotation key locations
      x = mc[mc[var_mol_index].peers[ci]].kx;
      y = mc[mc[var_mol_index].peers[ci]].ky;
      z = mc[mc[var_mol_index].peers[ci]].kz;
      mc[mc[var_mol_index].peers[ci]].kx = (R[0][0]*x) + (R[0][1]*y) + (R[0][2]*z);
      mc[mc[var_mol_index].peers[ci]].ky = (R[1][0]*x) + (R[1][1]*y) + (R[1][2]*z);
      mc[mc[var_mol_index].peers[ci]].kz = (R[2][0]*x) + (R[2][1]*y) + (R[2][2]*z);
    }

  } else {
    // Note that the 2D branch is not expected to be used in MCell, so it may be out of date.

    //#### fprintf ( stdout, "Rotating component positions for %s by %g\n", mc[var_mol_index].name, 180*angle/MY_PI );
    for (int ci=0; ci<mc[var_mol_index].num_peers; ci++) {
      //#### fprintf ( stdout, "  Component %s before is at (%g,%g)\n", mc[mc[var_mol_index].peers[ci]].name, mc[mc[var_mol_index].peers[ci]].x, mc[mc[var_mol_index].peers[ci]].y );
      double x = mc[mc[var_mol_index].peers[ci]].x;
      double y = mc[mc[var_mol_index].peers[ci]].y;
      mc[mc[var_mol_index].peers[ci]].x = (x * cos_angle) - (y * sin_angle);
      mc[mc[var_mol_index].peers[ci]].y = (x * sin_angle) + (y * cos_angle);
      //#### fprintf ( stdout, "  Component %s after  is at (%g,%g)\n", mc[mc[var_mol_index].peers[ci]].name, mc[mc[var_mol_index].peers[ci]].x, mc[mc[var_mol_index].peers[ci]].y );
    }
  }


  // Now the molecules are aligned as they should be except for rotation along their bonding axis

  if ( as3D && with_rot) {
    // Rotate the variable molecule along its bonding axis to align based on the rotation key angle



    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    double fixed_req_bond_angle = world->bond_angle / 2; // Radians ... should be a function of the bond!!
    double var_req_bond_angle = world->bond_angle / 2;   // Radians ... should be a function of the bond!!
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////



    // fixed_vcomp (fvc) will be the vector from the fixed molecule to the fixed component
    double fvc[3];
    fvc[0] = mc[fixed_comp_index].x - mc[fixed_mol_index].x;
    fvc[1] = mc[fixed_comp_index].y - mc[fixed_mol_index].y;
    fvc[2] = mc[fixed_comp_index].z - mc[fixed_mol_index].z;

    // var_vcomp (vvc) will be the vector from the var molecule to the var component
    double vvc[3];
    vvc[0] = mc[var_comp_index].x - mc[var_mol_index].x;
    vvc[1] = mc[var_comp_index].y - mc[var_mol_index].y;
    vvc[2] = mc[var_comp_index].z - mc[var_mol_index].z;


    // fixed_vkey (fvk) will be the vector from the fixed molecule to the fixed key
    double fvk[3];
    fvk[0] = mc[fixed_comp_index].kx - mc[fixed_mol_index].x;
    fvk[1] = mc[fixed_comp_index].ky - mc[fixed_mol_index].y;
    fvk[2] = mc[fixed_comp_index].kz - mc[fixed_mol_index].z;

    // var_vkey (vvk) will be the vector from the var molecule to the var key
    double vvk[3];
    vvk[0] = mc[var_comp_index].kx - mc[var_mol_index].x;
    vvk[1] = mc[var_comp_index].ky - mc[var_mol_index].y;
    vvk[2] = mc[var_comp_index].kz - mc[var_mol_index].z;

    if (world->dump_level >= 20) {
      fprintf ( stdout, "  Fixed vcomp = [ %g %g %g ]\n", fvc[0], fvc[1], fvc[2] );
      fprintf ( stdout, "  Var   vcomp = [ %g %g %g ]\n", vvc[0], vvc[1], vvc[2] );
      fprintf ( stdout, "  Fixed vkey  = [ %g %g %g ]\n", fvk[0], fvk[1], fvk[2] );
      fprintf ( stdout, "  Var vkey    = [ %g %g %g ]\n", vvk[0], vvk[1], vvk[2] );
    }

    // Use the cross product to get the normal to the fixed molecule-component-key plane
    double fixed_normal[3];
    fixed_normal[0] = (fvc[1] * fvk[2]) - (fvc[2] * fvk[1]);
    fixed_normal[1] = (fvc[2] * fvk[0]) - (fvc[0] * fvk[2]);
    fixed_normal[2] = (fvc[0] * fvk[1]) - (fvc[1] * fvk[0]);

    // Use the cross product to get the normal to the variable molecule-component-key plane
    double var_normal[3];
    var_normal[0] = (vvc[1] * vvk[2]) - (vvc[2] * vvk[1]);
    var_normal[1] = (vvc[2] * vvk[0]) - (vvc[0] * vvk[2]);
    var_normal[2] = (vvc[0] * vvk[1]) - (vvc[1] * vvk[0]);

    // Get the magnitudes of the two vectors for normalization
    double fixed_norm_mag = sqrt ( (fixed_normal[0]*fixed_normal[0]) + (fixed_normal[1]*fixed_normal[1]) + (fixed_normal[2]*fixed_normal[2]) );
    double var_norm_mag   = sqrt ( (  var_normal[0]*  var_normal[0]) + (  var_normal[1]*  var_normal[1]) + (  var_normal[2]*  var_normal[2]) );

    // Calculate unit vectors
    double fixed_unit[3];
    fixed_unit[0] = fixed_normal[0] / fixed_norm_mag;
    fixed_unit[1] = fixed_normal[1] / fixed_norm_mag;
    fixed_unit[2] = fixed_normal[2] / fixed_norm_mag;
    double var_unit[3];
    var_unit[0] = var_normal[0] / var_norm_mag;
    var_unit[1] = var_normal[1] / var_norm_mag;
    var_unit[2] = var_normal[2] / var_norm_mag;

    if (world->dump_level >= 20) {
      fprintf ( stdout, "  Fixed unit = [ %g %g %g ]\n", fixed_unit[0], fixed_unit[1], fixed_unit[2] );
      fprintf ( stdout, "  Var unit = [ %g %g %g ]\n", var_unit[0], var_unit[1], var_unit[2] );
    }

    double norm_dot_prod_again;
    norm_dot_prod_again = (fixed_unit[0] * var_unit[0]) + (fixed_unit[1] * var_unit[1]) + (fixed_unit[2] * var_unit[2]);

    // Ensure that the dot product is a legal argument for the "acos" function:
    if (norm_dot_prod_again >  1) {
      if (world->dump_level >= 20) {
        fprintf ( stdout, "Numerical Warning: normalized dot product %g was greater than 1\n", norm_dot_prod_again );
      }
      norm_dot_prod_again =  1;
    }
    if (norm_dot_prod_again < -1) {
      if (world->dump_level >= 20) {
        fprintf ( stdout, "Numerical Warning: normalized dot product %g was less than -1\n", norm_dot_prod_again );
      }
      norm_dot_prod_again = -1;
    }
    if (world->dump_level >= 20) {
      fprintf ( stdout, "  Normalized Dot Product between fixed and var is %g\n", norm_dot_prod_again );
    }

    // Compute the amount of rotation to bring the planes into alignment offset by the requested bond angles
    double cur_key_plane_angle = acos ( norm_dot_prod_again );

    if (world->dump_level >= 20) {
      fprintf ( stdout, "Current key plane angle = %g\n", (180*cur_key_plane_angle/MY_PI) );
    }

    double cross_prod[3];

    cross_prod[0] = (fixed_unit[1] * var_unit[2]) - (fixed_unit[2] * var_unit[1]);
    cross_prod[1] = (fixed_unit[2] * var_unit[0]) - (fixed_unit[0] * var_unit[2]);
    cross_prod[2] = (fixed_unit[0] * var_unit[1]) - (fixed_unit[1] * var_unit[0]);

    double dot_cross_rot = (cross_prod[0] * vvc[0]) + (cross_prod[1] * vvc[1]) + (cross_prod[2] * vvc[2]);
    if (dot_cross_rot > 0) {
      cur_key_plane_angle = (2*MY_PI) - cur_key_plane_angle;
    }

    if (world->dump_level >= 20) {
      fprintf ( stdout, "Current key plane angle = %g,  dot_cross_rot = %g\n", (180*cur_key_plane_angle/MY_PI), dot_cross_rot );
    }

    double composite_rot_angle = MY_PI + (var_req_bond_angle+fixed_req_bond_angle) + cur_key_plane_angle;  // The "MY_PI" adds 180 degrees to make the components "line up"

    if (world->dump_level >= 20) {
      fprintf ( stdout, "  Fixed angle                is = %g degrees\n", 180 * fixed_req_bond_angle / MY_PI );
      fprintf ( stdout, "  Var angle                  is = %g degrees\n", 180 * var_req_bond_angle / MY_PI );
      fprintf ( stdout, "  Current angle between keys is = %g degrees\n", 180 * cur_key_plane_angle / MY_PI );
      fprintf ( stdout, "  Composite rotation angle   is = %g degrees\n", 180 * composite_rot_angle / MY_PI );
    }

    // Build a 3D rotation matrix along the axis of the molecule to the component
    double var_vcomp_mag = sqrt ( (vvc[0]*vvc[0]) + (vvc[1]*vvc[1]) + (vvc[2]*vvc[2]) );

    double var_rot_unit[3];
    var_rot_unit[0] = vvc[0] / var_vcomp_mag;
    var_rot_unit[1] = vvc[1] / var_vcomp_mag;
    var_rot_unit[2] = vvc[2] / var_vcomp_mag;

    double ux = var_rot_unit[0];
    double uy = var_rot_unit[1];
    double uz = var_rot_unit[2];

    // Build the rotation matrix directly

    double cca = cos(composite_rot_angle);
    double sca = sin(composite_rot_angle);
    double omcca = 1 - cca;

    double R[3][3] = { { 1, 0, 0 },
                       { 0, 1, 0 },
                       { 0, 0, 1 } };

    R[0][0] = cca + (ux*ux*omcca);
    R[0][1] = (ux*uy*omcca) - (uz*sca);
    R[0][2] = (ux*uz*omcca) + (uy*sca);

    R[1][0] = (uy*ux*omcca) + (uz*sca);
    R[1][1] = cca + (uy*uy*omcca);
    R[1][2] = (uy*uz*omcca) - (ux*sca);

    R[2][0] = (uz*ux*omcca) - (uy*sca);
    R[2][1] = (uz*uy*omcca) + (ux*sca);
    R[2][2] = cca + (uz*uz*omcca);

    // Apply the rotation matrix after subtracting the molecule center location from all components and keys

    for (int ci=0; ci<mc[var_mol_index].num_peers; ci++) {
      // Rotate the component locations
      double x = mc[mc[var_mol_index].peers[ci]].x - mc[var_mol_index].x;
      double y = mc[mc[var_mol_index].peers[ci]].y - mc[var_mol_index].y;
      double z = mc[mc[var_mol_index].peers[ci]].z - mc[var_mol_index].z;
      mc[mc[var_mol_index].peers[ci]].x = (R[0][0]*x) + (R[0][1]*y) + (R[0][2]*z) + mc[var_mol_index].x;
      mc[mc[var_mol_index].peers[ci]].y = (R[1][0]*x) + (R[1][1]*y) + (R[1][2]*z) + mc[var_mol_index].y;
      mc[mc[var_mol_index].peers[ci]].z = (R[2][0]*x) + (R[2][1]*y) + (R[2][2]*z) + mc[var_mol_index].z;

      // Rotate the rotation key locations
      x = mc[mc[var_mol_index].peers[ci]].kx - mc[var_mol_index].x;
      y = mc[mc[var_mol_index].peers[ci]].ky - mc[var_mol_index].y;
      z = mc[mc[var_mol_index].peers[ci]].kz - mc[var_mol_index].z;
      mc[mc[var_mol_index].peers[ci]].kx = (R[0][0]*x) + (R[0][1]*y) + (R[0][2]*z) + mc[var_mol_index].x;
      mc[mc[var_mol_index].peers[ci]].ky = (R[1][0]*x) + (R[1][1]*y) + (R[1][2]*z) + mc[var_mol_index].y;
      mc[mc[var_mol_index].peers[ci]].kz = (R[2][0]*x) + (R[2][1]*y) + (R[2][2]*z) + mc[var_mol_index].z;
    }

  }

  //#### dump_molcomp_array(mc,num_parts);

  // Shift the var molecule location and the locations of all of its components by the difference of the binding components

  double dx = mc[fixed_comp_index].x - mc[var_comp_index].x;
  double dy = mc[fixed_comp_index].y - mc[var_comp_index].y;
  double dz = mc[fixed_comp_index].z - mc[var_comp_index].z;

  //#### fprintf ( stdout, "Shifting molecule and component positions for %s\n", mc[var_mol_index].name );
  mc[var_mol_index].x += dx;
  mc[var_mol_index].y += dy;
  mc[var_mol_index].z += dz;
  for (int ci=0; ci<mc[var_mol_index].num_peers; ci++) {
    // Shift the component locations
    mc[mc[var_mol_index].peers[ci]].x += dx;
    mc[mc[var_mol_index].peers[ci]].y += dy;
    mc[mc[var_mol_index].peers[ci]].z += dz;
    // Shift the rotation key locations
    mc[mc[var_mol_index].peers[ci]].kx += dx;
    mc[mc[var_mol_index].peers[ci]].ky += dy;
    mc[mc[var_mol_index].peers[ci]].kz += dz;
    //#### fprintf ( stdout, "  Component %s is at (%g,%g)\n", mc[mc[var_mol_index].peers[ci]].name, mc[mc[var_mol_index].peers[ci]].x, mc[mc[var_mol_index].peers[ci]].y );
  }

  //#### dump_molcomp_array(mc,num_parts);

  //#### fprintf ( stdout, "########## Done Binding %s to %s\n", mc[fixed_comp_index].name, mc[var_comp_index].name );
}


static void bind_all_molecules ( struct volume *world, external_molcomp_loc *molcomp_array, int num_parts, bool as3D, bool with_rot ) {
  // Compute positions for all molecules/components in a molcomp_array
  int mi=0;
  int pi=0;

  // Find first molecule
  for (mi=0; mi<num_parts; mi++) {
    if (molcomp_array[mi].is_mol) break;
  }
  if (molcomp_array[mi].is_mol) {
    // Set this first molecule and all of its components to final
    molcomp_array[mi].is_final = true;
    for (int ci=0; ci<molcomp_array[mi].num_peers; ci++) {
      molcomp_array[molcomp_array[mi].peers[ci]].is_final = true;
    }
    int done = 0;
    while (done == 0) {
      // Look for a bond between a non-final and a final component
      done = 1;
      for (mi=0; mi<num_parts; mi++) {
        if (!molcomp_array[mi].is_mol) {
          // Only search components for bonds
          if (molcomp_array[mi].num_peers > 1) {
            // This component has bonds, so search them
            for (int ci=1; ci<molcomp_array[mi].num_peers; ci++) {
              pi = molcomp_array[mi].peers[ci];  // Peer index
              if (molcomp_array[mi].is_final != molcomp_array[pi].is_final) {
                done = 0;
                // One of these is final and the other is not so join them and make them all final
                int fci, vci;  // Fixed comp index and Variable comp index
                int vmi;  // Fixed mol index and Variable mol index
                // Figure out which is fixed and which is not
                if (molcomp_array[mi].is_final) {
                  fci = mi;
                  vci = pi;
                } else {
                  fci = pi;
                  vci = mi;
                }
                // Set the molecule index values for the bond
                vmi = molcomp_array[vci].peers[0];

                // Set the initial (relative) positions of the components with each molecule at (0,0)
                // set_component_positions_2D ( world, molcomp_array, num_parts );

                // Perform the bond (changes the locations)
                bind_molecules_at_components ( world, molcomp_array, num_parts, fci, vci, as3D, with_rot );

                // Set the variable molecule and its components to final
                molcomp_array[vmi].is_final = true;
                for (int vmici=0; vmici<molcomp_array[vmi].num_peers; vmici++) {
                  molcomp_array[molcomp_array[vmi].peers[vmici]].is_final = true;
                }

              }
            }
          }
        }
      }
    }
  }

}

/*
static void set_molcomp_positions_2D_first_try ( struct volume *world, external_molcomp_loc *molcomp_array, int num_parts ) {
  // Compute positions for all molecules/components in a molcomp_array
  // This might be done recursively, but it's being done iteratively here first.
  // Note that this procedure does not work properly
  //#### fprintf ( stdout, "Begin Building molcomp_positions for:\n" );
  //#### dump_molcomp_array ( molcomp_array, num_parts );

  int all_done = 0;
  if (num_parts > 0) {
    // Note that a simple molecule will fall in this case and be at the origin
    molcomp_array[0].x = 0;
    molcomp_array[0].y = 0;
    molcomp_array[0].z = 0;
    molcomp_array[0].has_coords = 1;
  }
  while (all_done == 0) {
    // Continually go through the list and try to assign coordinates to each part
    double scale = 0.02;
    int i, j;
    all_done = 1;
    //#### fprintf ( stdout, "Top of loop through all molcomps\n" );
    for (i=0; i<num_parts; i++) {
      //#### fprintf ( stdout, "  Working on index %d named %s\n", i, molcomp_array[i].name );
      //#### if (molcomp_array[i].has_coords != 0) fprintf ( stdout, "    Index %d named %s has coordinates (%g,%g,%g).\n", i, molcomp_array[i].name, molcomp_array[i].x, molcomp_array[i].y, molcomp_array[i].z );
      if (molcomp_array[i].has_coords == 0) {
        // This molecule or component does NOT have any coordinates.
        // If it's a component, check to see if its molecule has coordinates
        // If it's a molecule, check to see if any of its components have coordinates
        if (molcomp_array[i].is_mol == 0) {
          //#### fprintf ( stdout, "    Index %d named %s has no coordinates.\n", i, molcomp_array[i].name );
          // This is a component. Components can only have 1 (mol) or 2 (mol,bond) peer entries
          if (molcomp_array[molcomp_array[i].peers[0]].has_coords != 0) {
            //#### fprintf ( stdout, "      Index %d named %s has a molecule with coordinates.\n", i, molcomp_array[i].name );
            // The molecule associated with this component has coordinates so use them first
            // Find out which "spoke" of the molecule contains this component
            int parent_mol_index = molcomp_array[i].peers[0];
            int comp_num_in_mol = 0;
            while (molcomp_array[parent_mol_index].peers[comp_num_in_mol] != i) {
              comp_num_in_mol += 1;
            }
            double angle = comp_num_in_mol * 2 * MY_PI / molcomp_array[parent_mol_index].num_peers;
            molcomp_array[i].x = molcomp_array[parent_mol_index].x + (scale * clipped_cos(angle));
            molcomp_array[i].y = molcomp_array[parent_mol_index].y + (scale * clipped_sin(angle));
            molcomp_array[i].z = molcomp_array[parent_mol_index].z;
            molcomp_array[i].has_coords = 1;
            //#### fprintf ( stdout, "        Index %d comp named %s is now at coordinates (%g,%g,%g)\n", i, molcomp_array[i].name, molcomp_array[i].x, molcomp_array[i].y, molcomp_array[i].z );
          } else if (molcomp_array[i].num_peers > 1) {
            // This component is bound to another, so check if those coordinates are defined
            //#### fprintf ( stdout, "      Index %d named %s is bound to a peer.\n", i, molcomp_array[i].name );
            if (molcomp_array[molcomp_array[i].peers[1]].has_coords != 0) {
              //#### fprintf ( stdout, "        Index %d named %s is bound to a peer with coordinates (%g,%g,%g).\n", i, molcomp_array[i].name, molcomp_array[molcomp_array[i].peers[1]].x, molcomp_array[molcomp_array[i].peers[1]].y, molcomp_array[molcomp_array[i].peers[1]].z );
              // The bound component has coordinates, so use those same coords for this component
              molcomp_array[i].x = molcomp_array[molcomp_array[i].peers[1]].x;
              molcomp_array[i].y = molcomp_array[molcomp_array[i].peers[1]].y;
              molcomp_array[i].z = molcomp_array[molcomp_array[i].peers[1]].z;
              molcomp_array[i].has_coords = 1;
            }
          }
        } else {
          // This is a molecule. Molecules can have any number of peer entries which are all components
          //#### fprintf ( stdout, "    Index %d named %s is a molecule with %d components.\n", i, molcomp_array[i].name, molcomp_array[i].num_peers );
          for (j=0; j<molcomp_array[i].num_peers; j++) {
            if (molcomp_array[molcomp_array[i].peers[j]].has_coords != 0) {
              // The jth component of this molecule has coordinates so use them
              int child_comp_index = molcomp_array[i].peers[j];
              //#### fprintf ( stdout, "      Index %d named %s has a %dth component with coordinates (%g,%g,%g).\n", i, molcomp_array[i].name, j, molcomp_array[child_comp_index].x, molcomp_array[child_comp_index].y, molcomp_array[child_comp_index].z );
              double angle = j * 2 * MY_PI / molcomp_array[i].num_peers;
              // Use the negative of the offset to locate this molecule relative to its component
              molcomp_array[i].x = molcomp_array[child_comp_index].x - (scale * clipped_cos(angle));
              molcomp_array[i].y = molcomp_array[child_comp_index].y - (scale * clipped_sin(angle));
              molcomp_array[i].z = molcomp_array[child_comp_index].z;
              molcomp_array[i].has_coords = 1;
              //#### fprintf ( stdout, "      Index %d mol named %s is now at coordinates (%g,%g,%g).\n", i, molcomp_array[i].name, molcomp_array[i].x, molcomp_array[i].y, molcomp_array[i].z );
              break;
            }
          }
        }
        all_done = false;

        //double scale = 0.048;
        // Temporarily assign random values
        //molcomp_array[i].x = scale * (drand48()-0.5) * .70710678118654752440;
        //molcomp_array[i].y = scale * (drand48()-0.5) * .70710678118654752440;
        //molcomp_array[i].z = scale * (drand48()-0.5) * .70710678118654752440;
        //molcomp_array[i].has_coords = 1;
      }
    }
  }

  //#### fprintf ( stdout, "Done Building molcomp_positions for:\n" );
  //#### dump_molcomp_array ( molcomp_array, num_parts );
}
*/

static external_molcomp_loc *build_molcomp_array ( struct volume *world, char **graph_strings ) {
  int part_num;
  char *next_part;

  /* Had trouble using the rng_state mechanism, skip for now
  if (rng == NULL) {
    rng = (rng_state *) malloc ( sizeof (struct rng_state) );
    rng_init ( rng, 12345 );
  }
  */

  part_num = 0;
  next_part = graph_strings[part_num];
  while (next_part != NULL) {
    part_num++;
    next_part = graph_strings[part_num];
  }

  // Allocate the entire block at once ...
  external_molcomp_loc *molcomp_loc_array = (external_molcomp_loc *) malloc ( part_num * sizeof(external_molcomp_loc) );

  // Copy the data into the block
  part_num = 0;
  next_part = graph_strings[part_num];
  while (next_part != NULL) {
    molcomp_loc_array[part_num].has_coords = 0;
    molcomp_loc_array[part_num].graph_string = (char *) malloc ( 1 + strlen(next_part) );
    strcpy ( molcomp_loc_array[part_num].graph_string, next_part );
    molcomp_loc_array[part_num].x = 0;
    molcomp_loc_array[part_num].y = 0;
    molcomp_loc_array[part_num].z = 0;
    if (strstr(next_part,"m:") == next_part) {
      // This is a molecule
      molcomp_loc_array[part_num].is_mol = 1;
      // For molecules, the name ends with ! or possibly the end of the string
      if (index(next_part,'!') == NULL) {
        molcomp_loc_array[part_num].name = (char *) malloc ( 1 + strlen(next_part) - 2 );
        strcpy ( molcomp_loc_array[part_num].name, &next_part[2] );
      } else {
        char *end_point = index(next_part,'!');
        *end_point = '\0';
        molcomp_loc_array[part_num].name = (char *) malloc ( 1 + strlen(next_part) - 2 );
        strcpy ( molcomp_loc_array[part_num].name, &next_part[2] );
        *end_point = '!';
      }
      // Remove any @ portions if they exist
      char *at_sign = index(molcomp_loc_array[part_num].name, '@');
      if (at_sign != NULL) {
        // Make a copy up to that point
        *at_sign = '\0';
        char *shorter_name = (char *) malloc ( 1 + strlen(molcomp_loc_array[part_num].name) );
        strcpy ( shorter_name, molcomp_loc_array[part_num].name );
        *at_sign = '@';
        free ( molcomp_loc_array[part_num].name );
        molcomp_loc_array[part_num].name = shorter_name;
      }
      // Get the molecule's neighbors which should all be components
      molcomp_loc_array[part_num].num_peers = 0;
      molcomp_loc_array[part_num].peers = NULL;
      char *next_excl = index(next_part,'!');
      while (next_excl != NULL) {
        molcomp_loc_array[part_num].num_peers++;
        next_excl++;
        next_excl = index(next_excl,'!');
      }
      if (molcomp_loc_array[part_num].num_peers > 0) {
        molcomp_loc_array[part_num].peers = (int *) malloc ( molcomp_loc_array[part_num].num_peers * sizeof(int) );
        next_excl = index(next_part,'!');
        int peer_num = 0;
        int comp_index;
        while (next_excl != NULL) {
          next_excl++;
          comp_index = atoi(next_excl);
          molcomp_loc_array[part_num].peers[peer_num] = comp_index;
          peer_num++;
          next_excl = index(next_excl,'!');
        }
      }
    } else {
      // This is a component
      molcomp_loc_array[part_num].is_mol = 0;
      // For components, the name ends with ~ or ! or possibly the end of the string
      char *first_exc = index(next_part,'!');
      char *first_til = index(next_part,'~');
      char *end_point;
      if ( (first_exc != NULL) && (first_til != NULL) ) {
        // Use whichever comes first
        if (first_exc < first_til) {
          end_point = first_exc;
        } else {
          end_point = first_til;
        }
      } else if (first_exc != NULL) {
        // Name ends at exc
        end_point = first_exc;
      } else if (first_til != NULL) {
        // Name ends at til
        end_point = first_til;
      } else {
        // Name ends at end of string
        end_point = index(next_part,'\0');
      }
      char previous_end = *end_point;
      *end_point = '\0';
      molcomp_loc_array[part_num].name = (char *) malloc ( 1 + strlen(next_part) - 2 );
      strcpy ( molcomp_loc_array[part_num].name, &next_part[2] );
      *end_point = previous_end;

      // Get the component's neighbors (the first will be the molecule)
      molcomp_loc_array[part_num].num_peers = 0;
      molcomp_loc_array[part_num].peers = NULL;
      char *next_excl = index(next_part,'!');
      while (next_excl != NULL) {
        molcomp_loc_array[part_num].num_peers++;
        next_excl++;
        next_excl = index(next_excl,'!');
      }
      if (molcomp_loc_array[part_num].num_peers > 0) {
        molcomp_loc_array[part_num].peers = (int *) malloc ( molcomp_loc_array[part_num].num_peers * sizeof(int) );
        next_excl = index(next_part,'!');
        int peer_num = 0;
        int comp_index;
        while (next_excl != NULL) {
          next_excl++;
          comp_index = atoi(next_excl);
          molcomp_loc_array[part_num].peers[peer_num] = comp_index;
          peer_num++;
          next_excl = index(next_excl,'!');
        }
      }
    }
    part_num++;
    next_part = graph_strings[part_num];
  }
  // set_molcomp_positions_2D_first_try ( world, molcomp_loc_array, part_num );
  // Set the initial (relative) positions of the components with each molecule at (0,0)

  // set_component_positions_2D ( world, molcomp_loc_array, part_num );
  set_component_positions_by_table ( world, molcomp_loc_array, part_num );

  bind_all_molecules ( world, molcomp_loc_array, part_num, true, true );

  if (world->dump_level >= 20) {
    fprintf ( stdout, ">>>>>>>>>>>>>>>>>>>>>>> Final molcomp_loc_array <<<<<<<<<<<<<<<<<<<\n" );
    dump_molcomp_array ( molcomp_loc_array, part_num );
    fprintf ( stdout, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n" );
  }

  return molcomp_loc_array;
}


static char **get_graph_strings ( char *nauty_string ) {
  // Parse the graph pattern
  // This code assumes that the graph pattern ends with a comma!!
  int num_parts;
  int part_num;
  char **graph_parts;
  char *first, *last;

  // Start by just counting the parts
  part_num = 0;
  first = nauty_string;
  last = strchr ( first, ',' );
  while (last != NULL) {
    first = last+1;
    last = strchr ( first, ',' );
    part_num++;
  }
  num_parts = part_num;

  // Allocate the array and mark the end
  graph_parts = (char **) malloc ( (num_parts+1) * sizeof(char *) );
  graph_parts[num_parts] = NULL; // Marks the end of the "list"

  // Copy the parts into the array of strings
  part_num = 0;
  first = nauty_string;
  last = strchr ( first, ',' );
  while (last != NULL) {
    *last = '\0';
    char *s = (char *) malloc ( strlen(first) + 1 );
    strcpy ( s, first );
    graph_parts[part_num] = s;
    *last = ',';
    first = last+1;
    last = strchr ( first, ',' );
    part_num++;
  }

  return graph_parts;
}



static void free_graph_parts ( char **graph_parts ) {
  if (graph_parts != NULL) {
    int part_num = 0;
    char *next_part = graph_parts[part_num];
    while (next_part != NULL) {
      free ( next_part );
      part_num++;
      next_part = graph_parts[part_num];
    }
    free ( graph_parts );
  }
}

/************************************************************************
output_cellblender_molecules:
In: vizblk: VIZ_OUTPUT block for this frame list
    a frame data list (internal viz output data structure)
Out: 0 on success, 1 on failure.  The names and positions of molecules are
     output in binary format designed for fast visualization in CellBlender

     Format of binary file is:
       Header:
         A single four-byte u_int containing version number of binary file
         format. This is value less than or equal to 16777215, which is
         0x00ffffff.  This allows automagic detection of ASCII and binary
         format molecule viz files, as well as the endianness of the binary
         format files during molecule viz in CellBlender.
       Molecule Viz Data:
         Version 0x00000001 files contain a block of binary data for each
         molecule species structured as follows:

           A single byte containing the length of the ASCII string of the
           molecule species name or state value, not including the terminating
           NULL. 32 chars max.

           The ASCII string containing the molecule species name or state value,
           not including the terminating NULL. 32 chars max.

           A single byte containing the molecule species type.
           Type 0 means volume molecule.  Type 1 means surface surface molecule.

           A four-byte u_int, N, whose value is 3 times the number of molecules
           of this species contained in the block, i.e. the number of floats
           in the block for the x,y,z coordinates of the molecule positions.

           A block of N four-byte floats containing the x,y,z coordinates of
           the molecule positions.

           If the molecules are surface molecules then block of x,y,z positions
           is followed by another block of N four-byte floats containing the
           i,j,k components of the orientation vector of the surface molecules.

         Note that the end of the file is indicated by the usual EOF only.

*************************************************************************/
static int output_cellblender_molecules(struct volume *world,
                                        struct viz_output_block *vizblk,
                                        struct frame_data_list *fdlp) {

  no_printf("Output in CELLBLENDER mode (molecules only)...\n");

  if ( (world->dump_level >= 0) && (world->viz_options != 0) ) {
    fprintf ( stdout, "vizblk->file_prefix_name = \"%s\"\n", vizblk->file_prefix_name );
  }
  if ( (world->dump_level >= 50) && (world->viz_options != 0) ) {
    fprintf ( stdout, "Visualization Options = 0x%lx\n", world->viz_options );
  }
  if (world->dump_level >= 50) {
    fprintf ( stdout, ">>>>>>>>>>>>>>>>>>>>>>> Top of MolViz Output <<<<<<<<<<<<<<<<<<<\n" );
  }

  // Unfortunately, the file name in vizblk->file_prefix_name has the "Scene" attached to the end.
  // This makes it difficult to use to build other paths that don't have a "Scene" in them.
  // For example: vizblk->file_prefix_name = "./viz_data/seed_00001/Scene"
  // So this will take some string "monkey business" to fix
  // We want "./viz_data/seed_#..#/ without the "Scene" attached
  char *file_prefix_no_Scene = NULL;
  char *file_prefix_usually_Scene = NULL;

  // Get the location of the last separator
  char *last_sep = strrchr ( vizblk->file_prefix_name, '/' );
  // Use the last_sep to copy the file part (usually Scene)
  file_prefix_usually_Scene = my_strcat ( last_sep+1, NULL );
  // Also use the last sep to copy the path part
  // Set it to \0 to mark the end of the string
  *last_sep = '\0';
  // Copy it with the desired prefix
  file_prefix_no_Scene = my_strcat ( vizblk->file_prefix_name, NULL );
  // Restore the vizblk->file_prefix_name
  *last_sep = '/';

  /*
  vizblk->file_prefix_name = "./viz_data/seed_00001/Scene"
  file_prefix_no_Scene = ./viz_data/seed_00001
  file_prefix_usually_Scene = Scene
  */

fprintf ( stdout, "path without file = %s\n", file_prefix_no_Scene );
fprintf ( stdout, "file without path = %s\n", file_prefix_usually_Scene );

  if ((fdlp->type == ALL_MOL_DATA) || (fdlp->type == MOL_POS)) {
    long long lli = 10;
    int ndigits = 1;
    for (; lli <= world->iterations && ndigits < 20;
         lli *= 10, ndigits++) {
    }
    char *cf_name =
        CHECKED_SPRINTF("%s.cellbin.%.*lld.dat", vizblk->file_prefix_name,
                        ndigits, fdlp->viz_iteration);
    if (cf_name == NULL)
      return 1;
    if (make_parent_dir(cf_name)) {
      free(cf_name);
      mcell_error(
          "Failed to create parent directory for CELLBLENDER-mode VIZ output.");
      /*return 1;*/
    }
    FILE *custom_file = open_file(cf_name, "wb");
    if (!custom_file)
      mcell_die();
    else {
      no_printf("Writing to file %s\n", cf_name);
    }
    free(cf_name);
    cf_name = NULL;

    FILE *space_struct_file = NULL;
    if (world->viz_options & 0x10L) {
      if (world->dump_level >= 20) {
        fprintf ( stdout, "Spatially Structured Option = 0x%lx\n", world->viz_options & 0x10L );
      }
      // Create the spatially structured mol file to hold the instances
      cf_name =
          CHECKED_SPRINTF("%s/viz_bngl/%s.bnglviz.%.*lld.dat", file_prefix_no_Scene, file_prefix_usually_Scene,
                          ndigits, fdlp->viz_iteration);
      if (cf_name == NULL)
        return 1;
      if (make_parent_dir(cf_name)) {
        free(cf_name);
        mcell_error(
            "Failed to create parent directory for SPATIAL-mode VIZ output.");
        /*return 1;*/
      }
      space_struct_file = open_file(cf_name, "wb");
      if (!space_struct_file) {
        mcell_die();
      } else {
        no_printf("Writing to file %s\n", cf_name);
      }
      free(cf_name);
      cf_name = NULL;
    }

    /* Get a list of molecules sorted by species. */
    u_int *viz_mol_count = NULL;
    struct abstract_molecule ***viz_molp = NULL;
    if (sort_molecules_by_species(
        world, vizblk, &viz_molp, &viz_mol_count, 1, 1)) {
      fclose(custom_file);
      custom_file = NULL;
      return 1;
    }

    /* Write file header */
    u_int cellbin_version = 1;
    fwrite(&cellbin_version, sizeof(cellbin_version), 1, custom_file);

    /* Write all the molecules whether EXTERNAL_SPECIES or not (for now) */

    for (int species_idx = 0; species_idx < world->n_species; species_idx++) {
      const unsigned int this_mol_count = viz_mol_count[species_idx];
      if (this_mol_count == 0)
        continue;

      const int id = vizblk->species_viz_states[species_idx];
      if (id == EXCLUDE_OBJ)
        continue;

      struct abstract_molecule **const mols = viz_molp[species_idx];
      if (mols == NULL)
        continue;

      /* Write species name: */
      struct abstract_molecule *amp = mols[0];
      char mol_name[33];
      if (id == INCLUDE_OBJ) {
        /* encode name of species as ASCII string, 32 chars max */
        snprintf(mol_name, 33, "%s", amp->properties->sym->name);
      } else {
        /* encode state value of species as ASCII string, 32 chars max */
        snprintf(mol_name, 33, "%d", id);
      }
      byte name_len = strlen(mol_name);
      fwrite(&name_len, sizeof(name_len), 1, custom_file);
      fwrite(mol_name, sizeof(char), name_len, custom_file);

      /* Write species type: */
      byte species_type = 0;
      if ((amp->properties->flags & ON_GRID) != 0) {
        species_type = 1;
      }
      fwrite(&species_type, sizeof(species_type), 1, custom_file);

      /* write number of x,y,z floats for mol positions to follow: */
      u_int n_floats = 3 * this_mol_count;
      fwrite(&n_floats, sizeof(n_floats), 1, custom_file);

      /* Write positions of volume and surface molecules: */
      float pos_x = 0.0;
      float pos_y = 0.0;
      float pos_z = 0.0;
      for (unsigned int n_mol = 0; n_mol < this_mol_count; ++n_mol) {
        amp = mols[n_mol];
        if ((amp->properties->flags & NOT_FREE) == 0) {
          struct volume_molecule *mp = (struct volume_molecule *)amp;
          struct vector3 pos_output = {0.0, 0.0, 0.0};
          if (!convert_relative_to_abs_PBC_coords(
              world->periodic_box_obj,
              mp->periodic_box,
              world->periodic_traditional,
              &mp->pos,
              &pos_output)) {
            pos_x = pos_output.x;   
            pos_y = pos_output.y;   
            pos_z = pos_output.z;   
          }
          else {
            pos_x = mp->pos.x; 
            pos_y = mp->pos.y; 
            pos_z = mp->pos.z; 
          }

        } else if ((amp->properties->flags & ON_GRID) != 0) {
          struct surface_molecule *gmp = (struct surface_molecule *)amp;
          struct vector3 where;
          uv2xyz(&(gmp->s_pos), gmp->grid->surface, &where);
          struct vector3 pos_output = {0.0, 0.0, 0.0};
          if (!convert_relative_to_abs_PBC_coords(
              world->periodic_box_obj,
              gmp->periodic_box,
              world->periodic_traditional,
              &where,
              &pos_output)) {
            pos_x = pos_output.x;   
            pos_y = pos_output.y;   
            pos_z = pos_output.z;   
          }
          else {
            pos_x = where.x; 
            pos_y = where.y; 
            pos_z = where.z; 
          }
        }

        pos_x *= world->length_unit;
        pos_y *= world->length_unit;
        pos_z *= world->length_unit;

        fwrite(&pos_x, sizeof(pos_x), 1, custom_file);
        fwrite(&pos_y, sizeof(pos_y), 1, custom_file);
        fwrite(&pos_z, sizeof(pos_z), 1, custom_file);
      }

      /* Write orientations of surface surface molecules: */
      amp = mols[0];
      if ((amp->properties->flags & ON_GRID) != 0) {
        for (unsigned int n_mol = 0; n_mol < this_mol_count; ++n_mol) {
          struct surface_molecule *gmp = (struct surface_molecule *)mols[n_mol];
          short orient = gmp->orient;
          float norm_x = orient * gmp->grid->surface->normal.x;
          float norm_y = orient * gmp->grid->surface->normal.y;
          float norm_z = orient * gmp->grid->surface->normal.z;

          if (world->periodic_box_obj && !(world->periodic_traditional)) {
            if (gmp->periodic_box->x % 2 != 0) {
              norm_x *= -1;
            }
            if (gmp->periodic_box->y % 2 != 0) {
              norm_y *= -1;
            }
            if (gmp->periodic_box->z % 2 != 0) {
              norm_z *= -1;
            }
          }

          fwrite(&norm_x, sizeof(norm_x), 1, custom_file);
          fwrite(&norm_y, sizeof(norm_y), 1, custom_file);
          fwrite(&norm_z, sizeof(norm_z), 1, custom_file);
        }
      }
    }

    /* Add additional Viz blocks for all EXTERNAL_SPECIES molecules */
    /* Note that this could be done while processing normal molecules, but separating makes code clearer. */

    external_mol_viz_by_name *mol_name_list = NULL;

    for (int species_idx = 0; species_idx < world->n_species; species_idx++) {
      const unsigned int this_mol_count = viz_mol_count[species_idx];

      if (this_mol_count == 0)
        continue;

      const int id = vizblk->species_viz_states[species_idx];
      if (id == EXCLUDE_OBJ)
        continue;

      struct abstract_molecule **const mols = viz_molp[species_idx];
      if (mols == NULL)
        continue;

      /* Get species name: */
      struct abstract_molecule *amp;
      amp = mols[0];
      char mol_name[33];
      if (id == INCLUDE_OBJ) {
        /* encode name of species as ASCII string, 32 chars max */
        snprintf(mol_name, 33, "%s", amp->properties->sym->name);
      } else {
        /* encode state value of species as ASCII string, 32 chars max */
        snprintf(mol_name, 33, "%d", id);
      }

      /* Get species type: */
      /*byte species_type = 0;*/
      /*if ((amp->properties->flags & ON_GRID) != 0) {*/
      /*  species_type = 1;*/
      /*}*/

      /* Get and save positions of EXTERNAL_SPECIES volume and surface molecules: */
      for (unsigned int n_mol = 0; n_mol < this_mol_count; ++n_mol) {
        amp = mols[n_mol];

        float pos_x = 0.0;
        float pos_y = 0.0;
        float pos_z = 0.0;
        float norm_x = 0.0;
        float norm_y = 0.0;
        float norm_z = 0.0;
        if ((amp->properties->flags & NOT_FREE) == 0) {
          struct volume_molecule *mp = (struct volume_molecule *)amp;
          pos_x = mp->pos.x;
          pos_y = mp->pos.y;
          pos_z = mp->pos.z;
        } else if ((amp->properties->flags & ON_GRID) != 0) {
          struct surface_molecule *gmp = (struct surface_molecule *)amp;
          struct vector3 where;
          uv2xyz(&(gmp->s_pos), gmp->grid->surface, &where);
          pos_x = where.x;
          pos_y = where.y;
          pos_z = where.z;
        }

        pos_x *= world->length_unit;
        pos_y *= world->length_unit;
        pos_z *= world->length_unit;

        char mol_type = 'v';
        if ((amp->properties->flags & ON_GRID) != 0) {
          mol_type = 's';
          struct surface_molecule *gmp = (struct surface_molecule *)mols[n_mol];
          short orient = gmp->orient;
          norm_x = orient * gmp->grid->surface->normal.x;
          norm_y = orient * gmp->grid->surface->normal.y;
          norm_z = orient * gmp->grid->surface->normal.z;
        }

        float x_offset = 0.0;

        if ((amp->properties->flags & EXTERNAL_SPECIES) != 0) {
          /* This is complex molecule, so add a new viz molecule for each molecule in the complex */
          /* The graph pattern will be something like: */
          /*    c:SH2~NO_STATE!5,c:U~NO_STATE!5!3,c:a~NO_STATE!6,c:b~Y!6!1,c:g~Y!6,m:Lyn@PM!0!1,m:Rec@PM!2!3!4, */
          char *next_mol = amp->graph_data->graph_pattern;

/* BEGIN NEW PROCESSING */
// fprintf ( stdout, "=============== BEGIN NEW PROCESSING ===============\n" );

          if (graph_pattern_table == NULL) {
            graph_pattern_table = init_symtab ( 10 );
          }

          struct sym_entry *sp;
          sp = retrieve_sym(next_mol, graph_pattern_table);

          if (sp == NULL) {

            // This pattern has not been saved yet, so parse it, save it, and possibly print it.

            char **graph_parts = get_graph_strings ( next_mol );

            // Print for use with external tools like the SpatialMols2D.java
            if (world->dump_level >= 10) {
              fprintf ( stdout, "=#= New Graph Pattern: %s\n", next_mol );
            }

            // Count the number of parts in graph_parts linked list
            int num_parts = 0;
            char *next_part = graph_parts[num_parts];
            while (next_part != NULL) {
              //#### fprintf ( stdout, "  Graph Part %d: %s\n", num_parts, next_part );
              num_parts++;
              next_part = graph_parts[num_parts];
            }

            external_molcomp_loc *molcomp_array = build_molcomp_array ( world, graph_parts );

            if (world->dump_level >= 10) {
              fprintf ( stdout, "=============== molcomp_array ===============\n" );
              dump_molcomp_array ( molcomp_array, num_parts );
              fprintf ( stdout, "=============================================\n" );
            }

            next_molcomp_id += 1;

            molcomp_list *mcl = (molcomp_list *) malloc ( sizeof(molcomp_list) );
            mcl->molcomp_array = molcomp_array;
            mcl->num_molcomp_items = num_parts;
            mcl->molcomp_id = next_molcomp_id;


            //#### fprintf ( stdout, "=============== molcomp_list ===============\n" );
            //#### dump_molcomp_list ( mcl );
            //#### fprintf ( stdout, "=============================================\n" );

            sp = store_sym ( next_mol, VOID_PTR, graph_pattern_table, mcl );

            //#### fprintf ( stdout, "=============== graph_pattern_table ===============\n" );
            //#### dump_symtab ( graph_pattern_table );
            //#### fprintf ( stdout, "===================================================\n" );

            free_graph_parts ( graph_parts );
          }

          molcomp_list *mcl = NULL;
          if (sp != NULL) {
            mcl = (molcomp_list *) sp->value;
            // fprintf ( stdout, "     sp: %s\n", mcl->name );
          }

          if (mcl != NULL) {

            // This should be the normal path

            // Note that this logic assumes that each part is either a molecule or a component
            // That's fine at the time of this design, but be aware that adding anything else
            //   that's other than a molecule or component may break this logic!!

            int part_num;

            for (part_num = 0; part_num<mcl->num_molcomp_items; part_num++) {
              // fprintf ( stdout, "    Mol Viz part_num %d\n", part_num );
              if ( mcl->molcomp_array[part_num].is_mol || (world->viz_options>0) ) {
                // fprintf ( stdout, "    mcl %s\n", mcl->molcomp_array[part_num].name );

                // Choose the name based on the type of item and the viz flags
                char *name_to_find_or_add = NULL;

                if (mcl->molcomp_array[part_num].is_mol) {

                  // Handle Molecule Viz

                  name_to_find_or_add = (char *) malloc (1+strlen(mcl->molcomp_array[part_num].name));
                  strcpy ( name_to_find_or_add, mcl->molcomp_array[part_num].name );

                } else {

                  // Handle Component Viz

                  int name_bits = world->viz_options & 0xf;
                  if (name_bits == 1) {
                    // Build a global component name alone (same glyph for all components and all molecules)
                    name_to_find_or_add = (char *) malloc (1+strlen("component"));
                    strcpy ( name_to_find_or_add, "component" );
                  } else if (name_bits == 2) {
                    // Build the name from the component name alone (same glyph for all components with this name across all molecules)
                    name_to_find_or_add = (char *) malloc (1+strlen("comp_")+strlen(mcl->molcomp_array[part_num].name));
                    strcpy ( name_to_find_or_add, "comp_" );
                    strcpy ( &name_to_find_or_add[strlen("comp_")], mcl->molcomp_array[part_num].name );
                  } else if (name_bits == 3) {
                    // Build the name from the molecule name and the component name (glyph only applies to this mol/comp combination)
                    char *last_mol_name = NULL;
                    if (mcl->molcomp_array[part_num].num_peers < 1) {
                      // This shouldn't happen ...
                      last_mol_name = (char *) malloc (1+strlen("unknown_"));
                      strcpy ( last_mol_name, "unknown_" );
                    } else {
                      // Copy the molecule name
                      char *name_ptr = mcl->molcomp_array[mcl->molcomp_array[part_num].peers[0]].name;
                      last_mol_name = (char *) malloc (1+strlen(name_ptr));
                      strcpy ( last_mol_name, name_ptr );
                    }
                    name_to_find_or_add = (char *) malloc (1+strlen(last_mol_name)+strlen("_comp_")+strlen(mcl->molcomp_array[part_num].name));
                    strcpy ( name_to_find_or_add,                               last_mol_name );
                    strcpy ( &name_to_find_or_add[strlen(name_to_find_or_add)], "_comp_" );
                    strcpy ( &name_to_find_or_add[strlen(name_to_find_or_add)], mcl->molcomp_array[part_num].name );

                    if (last_mol_name != NULL) {
                      free ( last_mol_name );
                    }
                  }

                }

                /* Check to see if this name is already in the mol_name_list */
                external_mol_viz_by_name *next_mol_name = mol_name_list;
                int found = 0;
                // Check for the actual mol name being in the list
                do {
                  if (next_mol_name == NULL) {
                    break;
                  }
                  if (strcmp(name_to_find_or_add, next_mol_name->mol_name) == 0) {
                    found = 1;
                    break;
                  }
                  next_mol_name = next_mol_name->next_name;
                } while ( found == 0 );

                if (found == 0) {
                  /* This molecule or component name is not in the list, so add a new name to the front */
                  next_mol_name = (external_mol_viz_by_name *) malloc ( sizeof(external_mol_viz_by_name) );
                  next_mol_name->mol_name = name_to_find_or_add;  /* This takes "ownership" of the allocated "name_to_find_or_add" memory */
                  next_mol_name->mol_list = NULL;
                  next_mol_name->next_name = mol_name_list;
                  mol_name_list = next_mol_name;
                }

                /* next_mol_name now points to the list of molecules by this name */

                /* Make a new molecule viz item to store this location */
						    external_mol_viz *new_mol_viz_item = (external_mol_viz *) malloc ( sizeof(external_mol_viz) );

                /* Set its values */
                new_mol_viz_item->mol_type = mol_type;

                new_mol_viz_item->pos_x = pos_x + mcl->molcomp_array[part_num].x;
                new_mol_viz_item->pos_y = pos_y + mcl->molcomp_array[part_num].y;
                new_mol_viz_item->pos_z = pos_z + mcl->molcomp_array[part_num].z;

                // fprintf ( stdout, "=MVM= Mol Viz Mol at: (%g,%g,%g)\n", new_mol_viz_item->pos_x, new_mol_viz_item->pos_y, new_mol_viz_item->pos_z );

                new_mol_viz_item->norm_x = norm_x;
                new_mol_viz_item->norm_y = norm_y;
                new_mol_viz_item->norm_z = norm_z;

                /* Add it to the top of the molecule list */
                new_mol_viz_item->next_mol = next_mol_name->mol_list;
                next_mol_name->mol_list = new_mol_viz_item;

                next_mol += 1;

              }
            }

// fprintf ( stdout, "=============== END NEW PROCESSING ===============\n" );

/* END NEW PROCESSING */

          } else {

            // This was the old way of making mols directly from the NAUTY strings.
            // This code either placed all molecule parts at the same location or
            //   spaced them evenly along the "x" axis (see x_offset below).
            // This should not generally get executed since the previous processing was added.
            // It's being kept here because it is computationally simpler and may be useful.

            while ((next_mol = strstr(next_mol,"m:")) != NULL ) {
              /* Pull the next actual molecule name out of the graph pattern */
              char *end_mol = strpbrk ( next_mol, "@!,(~" );
              if (end_mol == NULL) {
                end_mol = next_mol + strlen(next_mol);
              }
              int ext_name_len = end_mol - next_mol;
              char *ext_name = (char *) malloc ( ext_name_len + 1 );
              strncpy ( ext_name, next_mol+2, ext_name_len-2 );
              ext_name[ext_name_len-2] = '\0';

              /* Check to see if this name is already in the list */
              external_mol_viz_by_name *next_mol_name = mol_name_list;
              int found = 0;
              do {
                if (next_mol_name == NULL) {
                  break;
                }
                if (strcmp(ext_name, next_mol_name->mol_name) == 0) {
                  found = 1;
                  break;
                }
                next_mol_name = next_mol_name->next_name;
              } while ( found == 0 );

              if (found == 0) {
                /* This molecule name is not in the list, so add a new name to the front */
                next_mol_name = (external_mol_viz_by_name *) malloc ( sizeof(external_mol_viz_by_name) );
                next_mol_name->mol_name = ext_name;  /* This takes "ownership" of the allocated name memory */
                next_mol_name->mol_list = NULL;
                next_mol_name->next_name = mol_name_list;
                mol_name_list = next_mol_name;
              } else {
                /* This molecule name is already in the list and next_mol_name points to it, so just free the name. */
                free ( ext_name );
              }

              /* next_mol_name now points to the list of molecules by this name */

              /* Make a new molecule viz item to store this location */
						  external_mol_viz *new_mol_viz_item = (external_mol_viz *) malloc ( sizeof(external_mol_viz) );

              /* Set its values */
              new_mol_viz_item->mol_type = mol_type;

              new_mol_viz_item->pos_x = pos_x + x_offset;  x_offset += 0.008;
              new_mol_viz_item->pos_y = pos_y;
              new_mol_viz_item->pos_z = pos_z;

              new_mol_viz_item->norm_x = norm_x;
              new_mol_viz_item->norm_y = norm_y;
              new_mol_viz_item->norm_z = norm_z;

              /* Add it to the top of the molecule list */
              new_mol_viz_item->next_mol = next_mol_name->mol_list;
              next_mol_name->mol_list = new_mol_viz_item;

              next_mol += 1;
            }

          }

        }
      }
    }

    /* Write out the molecules with their proper names */
    external_mol_viz_by_name *nl = mol_name_list;
    external_mol_viz *mv;

    while (nl != NULL) {

      /* Write the name length and name */
      byte name_len = strlen(nl->mol_name);
      fwrite(&name_len, sizeof(name_len), 1, custom_file);
      fwrite(nl->mol_name, sizeof(char), name_len, custom_file);

      /* Write species type: */
      byte species_type = 0;
      if (nl->mol_list != NULL) {
        if (nl->mol_list->mol_type == 's') {
          species_type = 1;
        }
      }
      fwrite(&species_type, sizeof(species_type), 1, custom_file);

      /* write number of x,y,z floats for mol positions to follow: */
      u_int n_floats = 0;
      mv = nl->mol_list;
      while (mv != NULL) {
        n_floats += 3;
        mv = mv->next_mol;
      }
      fwrite(&n_floats, sizeof(n_floats), 1, custom_file);

      /* Write positions of volume and surface molecules: */
      mv = nl->mol_list;
      while (mv != NULL) {
        float pos_x = mv->pos_x;
        float pos_y = mv->pos_y;
        float pos_z = mv->pos_z;
        fwrite(&pos_x, sizeof(pos_x), 1, custom_file);
        fwrite(&pos_y, sizeof(pos_y), 1, custom_file);
        fwrite(&pos_z, sizeof(pos_z), 1, custom_file);
        mv = mv->next_mol;
      }
      /* Write orientations of surface surface molecules: */
      mv = nl->mol_list;
      if (mv->mol_type == 's') {
        while (mv != NULL) {
          float norm_x = mv->norm_x;
          float norm_y = mv->norm_y;
          float norm_z = mv->norm_z;
          fwrite(&norm_x, sizeof(norm_x), 1, custom_file);
          fwrite(&norm_y, sizeof(norm_y), 1, custom_file);
          fwrite(&norm_z, sizeof(norm_z), 1, custom_file);
          mv = mv->next_mol;
        }
      }
      nl = nl->next_name;
    }

    /* Free the structures used to build the molecule lists from the complexes */
    while (mol_name_list != NULL) {
      nl = mol_name_list;
      /* Free the name block */
      free ( nl->mol_name );
      /* Free the list of molecule instances */
      while (nl->mol_list != NULL) {
        mv = nl->mol_list;
        nl->mol_list = mv->next_mol;
        free ( mv );
      }
      mol_name_list = nl->next_name;
      free ( nl );
    }

    fclose(custom_file);
    custom_file = NULL;

    free_ptr_array((void **)viz_molp, world->n_species);
    viz_molp = NULL;
    free(viz_mol_count);
    viz_mol_count = NULL;
  }

  free ( file_prefix_no_Scene );
  free ( file_prefix_usually_Scene );

  if (world->dump_level >= 50) {
    fprintf ( stdout, ">>>>>>>>>>>>>>>>>>>>>>> Bottom of MolViz Output <<<<<<<<<<<<<<<<<<<\n" );
  }

  return 0;
}

/*********************************************************************
init_frame_data_list:

   In: vizblk: the VIZ_OUTPUT block to initialize
   Out: 0 on success, 1 on error.
        Initializes frame_data_list structure.
        Sets the value of the current iteration step to the start value.
        Sets the number of iterations.
***********************************************************************/
int init_frame_data_list(struct volume *world,
                         struct viz_output_block *vizblk) {
  int mol_orient_frame_present = 0;
  int mol_pos_frame_present = 0;
  struct frame_data_list *fdlp;

  if (vizblk->frame_data_head == NULL)
    return 0;

  switch (vizblk->viz_mode) {
  case NO_VIZ_MODE:
    count_time_values(world, vizblk->frame_data_head);
    if (reset_time_values(world, vizblk->frame_data_head, world->start_iterations))
      return 1;
    break;
  //      return 0;

  case ASCII_MODE:
    count_time_values(world, vizblk->frame_data_head);
    if (reset_time_values(world, vizblk->frame_data_head, world->start_iterations))
      return 1;
    break;

  case CELLBLENDER_MODE:
    count_time_values(world, vizblk->frame_data_head);
    if (reset_time_values(world, vizblk->frame_data_head, world->start_iterations))
      return 1;
    break;

  default:
    count_time_values(world, vizblk->frame_data_head);
    if (reset_time_values(world, vizblk->frame_data_head, world->start_iterations))
      return 1;
    break;
  }

  for (fdlp = vizblk->frame_data_head; fdlp != NULL; fdlp = fdlp->next) {
    if (fdlp->curr_viz_iteration == NULL)
      continue;

    switch (fdlp->type) {
    case MOL_ORIENT:
      mol_orient_frame_present = 1;
      break;

    case MOL_POS:
      mol_pos_frame_present = 1;
      break;

    case ALL_MOL_DATA:
      mol_pos_frame_present = 1;
      mol_orient_frame_present = 1;
      break;

    default:
      /* Do nothing */
      ;
    }
  } /* end while */

  /* Check that the user hasn't selected a useless set of output info */
  if ((mol_orient_frame_present) & (!mol_pos_frame_present))
    mcell_warn("The input file contains ORIENTATIONS but not POSITIONS "
               "statement in the MOLECULES block. The molecules cannot be "
               "visualized.");

  return 0;
}

/**************************************************************************
update_frame_data_list:
        In: vizblk: VIZ_OUTPUT block
        Out: 0 on success, 1 on failure.
             Calls output visualization functions if necessary.
             Updates value of the current iteration step and pointer
             to the current iteration in the linked list.
**************************************************************************/
int update_frame_data_list(struct volume *world,
                           struct viz_output_block *vizblk) {
  static char const *const FRAME_TYPES[NUM_FRAME_TYPES] = {
    "MOL_POS",  "MOL_ORIENT", "ALL_MOL_DATA", 
  };

  if (vizblk == NULL)
    return 0;
  if (vizblk->frame_data_head == NULL)
    return 0;

  switch (world->notify->viz_output_report) {
  case NOTIFY_NONE:
    break;

  case NOTIFY_BRIEF:
  case NOTIFY_FULL:
    mcell_log("Updating viz output on iteration %lld.", world->current_iterations);
    break;

  default:
    UNHANDLED_CASE(world->notify->viz_output_report);
  }

  /* Scan over all frames, producing appropriate output. */
  for (struct frame_data_list *fdlp = vizblk->frame_data_head; fdlp != NULL;
       fdlp = fdlp->next) {
    if (world->current_iterations != fdlp->viz_iteration)
      continue;

    if (world->notify->viz_output_report == NOTIFY_FULL) {
      if (fdlp->type >= NUM_FRAME_TYPES)
        mcell_warn("  Updating data frame of unknown type %d.", fdlp->type);
      else
        mcell_log("  Updating data frame of type %s.", FRAME_TYPES[fdlp->type]);
    }

    switch (vizblk->viz_mode) {
    case ASCII_MODE:
      if (output_ascii_molecules(world, vizblk, fdlp))
        return 1;
      break;

    case CELLBLENDER_MODE:
      if (output_cellblender_molecules(world, vizblk, fdlp))
        return 1;
      break;

    case NO_VIZ_MODE:
    default:
      /* Do nothing for vizualization */
      break;
    }

    while (fdlp->curr_viz_iteration != NULL &&
           fdlp->viz_iteration == world->current_iterations) {
      fdlp->curr_viz_iteration = fdlp->curr_viz_iteration->next;
      if (fdlp->curr_viz_iteration)
        fdlp->viz_iteration = frame_iteration(
            world, fdlp->curr_viz_iteration->value, fdlp->list_type);
    }
    if (world->notify->viz_output_report == NOTIFY_FULL)
      mcell_log("  Next update on iteration %lld.", fdlp->viz_iteration);
  }
  return 0;
}

/**************************************************************************
finalize_viz_output:
        In: vizblk: VIZ_OUTPUT block
        Out: Returns 1 on error and zero otherwise. Writes final information
             into visualization output files.
**************************************************************************/
int finalize_viz_output(struct volume *world, struct viz_output_block *vizblk) {
  if (vizblk == NULL)
    return 0;

  switch (vizblk->viz_mode) {
  case NO_VIZ_MODE:
  case ASCII_MODE:
  default:
    /* Do nothing for vizualization */
    break;
  }

  return 0;
}
