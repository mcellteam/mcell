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


typedef struct external_molcomp_loc_struct {
  bool is_mol;
  double x, y, z;
  char *name;
  char *graph_string;
  int num_peers;
  int *peers;
} external_molcomp_loc;

static void dump_molcomp_array ( external_molcomp_loc *molcomp_array, int num_parts ) {
  int i, j;
  for (i=0; i<num_parts; i++) {
    if (molcomp_array[i].is_mol) {
      fprintf ( stdout, "m" );
    } else {
      fprintf ( stdout, "c" );
    }
    fprintf ( stdout, "[%d]: %s", i, molcomp_array[i].name );
    for (j=0; j<molcomp_array[i].num_peers; j++) {
      fprintf ( stdout, "   [%d]", molcomp_array[i].peers[j] );
    }
    fprintf ( stdout, "\n" );
  }
}

static external_molcomp_loc *build_molcomp_array ( char **graph_strings ) {
  int part_num;
  char *next_part;

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
    molcomp_loc_array[part_num].x = 0;
    molcomp_loc_array[part_num].y = 0;
    molcomp_loc_array[part_num].z = 0;
    molcomp_loc_array[part_num].graph_string = (char *) malloc ( 1 + strlen(next_part) );
    strcpy ( molcomp_loc_array[part_num].graph_string, next_part );
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

          if (graph_pattern_table == NULL) {
            graph_pattern_table = init_symtab ( 10 );
          }

          struct sym_entry *sp;
          sp = retrieve_sym(next_mol, graph_pattern_table);

          if (sp == NULL) {

            // This pattern has not been saved yet, so parse it, print it, and save it

            char **graph_parts = get_graph_strings ( next_mol );

            fprintf ( stdout, "#=# Graph Pattern: %s\n", next_mol );

            int part_num = 0;
            char *next_part = graph_parts[part_num];
            while (next_part != NULL) {
              fprintf ( stdout, "  Graph Part %d: %s\n", part_num, next_part );
              part_num++;
              next_part = graph_parts[part_num];
            }

            struct sym_entry *stored_mol = store_sym ( next_mol, VOID_PTR, graph_pattern_table, NULL );

            fprintf ( stdout, "=============== graph_pattern_table ===============\n" );
            dump_symtab ( graph_pattern_table );
            fprintf ( stdout, "===================================================\n" );

            external_molcomp_loc *molcomp_array = build_molcomp_array ( graph_parts );

            fprintf ( stdout, "=============== molcomp_array ===============\n" );
            dump_molcomp_array ( molcomp_array, part_num );
            fprintf ( stdout, "=============================================\n" );

            free_graph_parts ( graph_parts );
          }

/* END NEW PROCESSING */

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
