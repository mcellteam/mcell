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

#include "mcell_structs.h"
#include "grid_util.h"
#include "sched_util.h"
#include "viz_output.h"
#include "strfunc.h"
#include "util.h"

static const int endian_test_word = 0x04030201;

extern struct volume *world;

/* == viz-specific Utilities == */

/*************************************************************************
frame_iteration:
    Gets the iteration number for a given time/iteration value and "type".

        In:  double iterval - the time/iteration value
             int type - the type of value
        Out: the frame time as an iteration number
**************************************************************************/
static long long frame_iteration(double iterval, int type)
{
  switch (type)
  {
    case OUTPUT_BY_ITERATION_LIST:
      return (long long) iterval;

    case OUTPUT_BY_TIME_LIST:
      return (long long) (iterval / world->time_unit + ROUND_UP);

    default:
      fprintf(world->err_file, "File '%s', Line %ld: error - wrong frame_data_list list_type %d\n", __FILE__, (long)__LINE__, type);
      return -1;
  }
}

/*************************************************************************
sort_molecules_by_species:
    Scans over all molecules, sorting them into arrays by species.

        In:  struct abstract_molecule ****viz_molpp
             u_int  **viz_mol_countp
             int include_volume - should the lists include vol mols?
             int include_grid - should the lists include grid mols?
        Out: 0 on success, 1 on error; viz_molpp and viz_mol_countp arrays are
             allocated and filled with sorted data.
**************************************************************************/
static int sort_molecules_by_species(struct abstract_molecule ****viz_molpp,
                                     u_int **viz_mol_countp,
                                     int include_volume,
                                     int include_grid)
{
  struct storage_list *slp;
  u_int *counts;
  int species_index;

  /* XXX: May leave memory allocated on failure */
  if ((*viz_molpp = (struct abstract_molecule ***) allocate_ptr_array(world->n_species)) == NULL)
    return 1;
  if ((counts = *viz_mol_countp = allocate_uint_array(world->n_species, 0)) == NULL)
    return 1;

  /* Walk through the species allocating arrays for all molecules of that species */
  for (species_index = 0; species_index < world->n_species; ++ species_index)
  {
    int mol_count;
    u_int spec_id = world->species_list[species_index]->species_id;

    if (world->species_list[species_index]->viz_state == EXCLUDE_OBJ)
      continue;

    if (world->species_list[species_index]->flags & IS_SURFACE)
      continue;

    if (! include_grid  &&  (world->species_list[species_index]->flags & ON_GRID))
      continue;

    if (! include_volume  &&  ! (world->species_list[species_index]->flags & ON_GRID))
      continue;

    mol_count = world->species_list[species_index]->population;
    if (mol_count <= 0)
      continue;

    if (((*viz_molpp)[spec_id] = (struct abstract_molecule **) allocate_ptr_array(mol_count)) == NULL)
      return 1;
  }

  /* Sort molecules by species id */
  for (slp = world->storage_head;
       slp != NULL;
       slp = slp->next)
  {
    struct storage *sp = slp->store;
    struct schedule_helper *shp;
    struct abstract_molecule *amp;
    int sched_slot_index;
    for (shp = sp->timer;
         shp != NULL;
         shp = shp->next_scale)
    {
      for (sched_slot_index = -1; sched_slot_index < shp->buf_len; ++ sched_slot_index)
      {
        for (amp = (struct abstract_molecule *) ((sched_slot_index < 0) ? shp->current : shp->circ_buf_head[sched_slot_index]);
             amp != NULL;
             amp = amp->next)
        {
          u_int spec_id;
          if (amp->properties == NULL)
            continue;

          if (amp->properties->viz_state == EXCLUDE_OBJ)
            continue;

          if (! include_grid  &&  (amp->flags & TYPE_MASK) != TYPE_3D)
            continue;

          if (! include_volume  &&  (amp->flags & TYPE_MASK) == TYPE_3D)
            continue;

          spec_id = amp->properties->species_id;
          if (counts[spec_id] < amp->properties->population)
            (*viz_molpp)[spec_id][counts[spec_id]++] = amp;
          else {
            fprintf(world->log_file,"MCell: molecule count disagreement!!\n");
            fprintf(world->log_file,"  Species %s  population = %d  count = %d\n",
                    amp->properties->sym->name,
                    amp->properties->population,
                    counts[spec_id]);
          }
        }
      }
    }
  }

  return 0;
}

/*************************************************************************
active_this_iteration:
    Check if any visualization is to be done this iteration.

        In: struct frame_data_list *fdlp - the list of frame data
        Out: 0 on success, 1 on failure
**************************************************************************/
static int active_this_iteration(struct frame_data_list *fdlp)
{
  for (; fdlp != NULL; fdlp = fdlp->next)
  {
    if (fdlp->viz_iteration == world->it_time)
      return 1;
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
static int reset_time_values(struct frame_data_list *fdlp,
                             long long curiter)
{
  /* If we've loaded a checkpoint, don't output on the first iter */
  if (curiter != 0)
    ++ curiter;

  for (; fdlp != NULL; fdlp = fdlp->next)
  {
    fdlp->curr_viz_iteration = fdlp->iteration_list;
    fdlp->viz_iteration = -1;

    /* Scan for first iteration >= curiter */
    while (fdlp->curr_viz_iteration != NULL)
    {
      if (frame_iteration(fdlp->curr_viz_iteration->value, fdlp->list_type) >= curiter)
        break;
      fdlp->curr_viz_iteration = fdlp->curr_viz_iteration->next;
    }

    /* If we had an iteration, use it to set viz_iteration */
    if (fdlp->curr_viz_iteration != NULL)
      fdlp->viz_iteration = frame_iteration(fdlp->curr_viz_iteration->value, fdlp->list_type);
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
static int count_time_values(struct frame_data_list * const fdlp)
{
  int time_values = 0;
  long long curiter = -1;
  struct frame_data_list *fdlpcur = NULL;
  for (fdlpcur = fdlp; fdlpcur != NULL; fdlpcur = fdlpcur->next)
  {
    fdlpcur->curr_viz_iteration = fdlpcur->iteration_list;
    fdlpcur->n_viz_iterations = 0;
    fdlpcur->viz_iteration = -1;
  }

  while (1)
  {
    curiter = -1;

    /* Find the next iteration */
    for (fdlpcur = fdlp; fdlpcur != NULL; fdlpcur = fdlpcur->next)
    {
      long long thisiter;
      if (fdlpcur->curr_viz_iteration == NULL)
        continue;

      thisiter = frame_iteration(fdlpcur->curr_viz_iteration->value, fdlpcur->list_type);

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
    if (world->chkpt_iterations != 0  &&  curiter > world->start_time + world->chkpt_iterations)
      break;

    /* Optimistically, store this as the "final" iteration */
    if (curiter > world->viz_state_info.final_iteration)
      world->viz_state_info.final_iteration = curiter;

    /* We found at least one more.  Note that the only time we will output at
     * iteration == start_time is when start_time is zero.  This is because we
     * do not output on the first iteration after we resume.
     */
    if (curiter > world->start_time)
      ++ time_values;
    else if ((world->start_time | curiter) == 0)
      ++ time_values;

    /* Advance any frame data items which are set to this iteration */
    for (fdlpcur = fdlp; fdlpcur != NULL; fdlpcur = fdlpcur->next)
    {
      if (fdlpcur->curr_viz_iteration == NULL)
        continue;

      if (curiter > world->start_time  ||  (world->start_time | curiter) == 0)
      {
        if (frame_iteration(fdlpcur->curr_viz_iteration->value, fdlpcur->list_type) == curiter)
          ++ fdlpcur->n_viz_iterations;
      }

      while (fdlpcur->curr_viz_iteration  &&
             frame_iteration(fdlpcur->curr_viz_iteration->value, fdlpcur->list_type) == curiter)
        fdlpcur->curr_viz_iteration = fdlpcur->curr_viz_iteration->next;
    }
  }

  return time_values;
}

/*************************************************************************
initialize_iteration_counters:
    Sets up iteration counter buffers, allocating space for iteration and time
    data, which is needed for DREAMM output.

        In:  int time_values_total - the maximum number of time/iteration
                                     elements
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int initialize_iteration_counters(int time_values_total)
{
  memset(&world->viz_state_info.output_times, 0, sizeof(world->viz_state_info.output_times));
  memset(&world->viz_state_info.mesh_output_iterations, 0, sizeof(world->viz_state_info.mesh_output_iterations));
  memset(&world->viz_state_info.vol_mol_output_iterations, 0, sizeof(world->viz_state_info.vol_mol_output_iterations));
  memset(&world->viz_state_info.grid_mol_output_iterations, 0, sizeof(world->viz_state_info.grid_mol_output_iterations));

  if (initialize_iteration_counter(&world->viz_state_info.output_times, time_values_total))
    goto failure;
  if (initialize_iteration_counter(&world->viz_state_info.mesh_output_iterations, time_values_total))
    goto failure;
  if (initialize_iteration_counter(&world->viz_state_info.vol_mol_output_iterations, time_values_total))
    goto failure;
  if (initialize_iteration_counter(&world->viz_state_info.grid_mol_output_iterations, time_values_total))
    goto failure;

  return 0;

failure:
  destroy_iteration_counter(&world->viz_state_info.output_times);
  destroy_iteration_counter(&world->viz_state_info.mesh_output_iterations);
  destroy_iteration_counter(&world->viz_state_info.vol_mol_output_iterations);
  destroy_iteration_counter(&world->viz_state_info.grid_mol_output_iterations);
  return 1;
}

/*************************************************************************
collect_objects:
    Collect all objects being visualized.  Pointers to objects will be sorted
    by address so that a binary search may be used on them.

        In:  struct viz_obj *root - the head of the viz_obj linked list
             struct object ***objs - pointer to receive array of object *s
             int *num_objects - int to receive object count
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int collect_objects(struct viz_obj *root,
                           struct object ***objs,
                           int *num_objects)
{
  struct viz_obj *vizp;

  /* Count up the objects */
  int count = 0;
  for (vizp = root; vizp != NULL; vizp = vizp->next)
  {
    struct viz_child *vcp;
    for (vcp = vizp->viz_child_head;
         vcp != NULL;
         vcp = vcp->next)
      ++ count;
  }

  /* Allocate a flattened array for the objects */
  if ((*objs = (struct object **) allocate_ptr_array(count)) == NULL)
    return 1;

  *num_objects = count;

  /* Traverse the tree again and fill the objects array */
  struct object **curobjptr = *objs;
  for (vizp = root;
       vizp != NULL  &&  count > 0;
       vizp = vizp->next)
  {
    struct viz_child *vcp;
    for (vcp = vizp->viz_child_head;
         vcp != NULL;
         vcp = vcp->next)
    {
      *(curobjptr++) = vcp->obj;
      if (-- count == 0)
        break;
    }
  }

  qsort(*objs, *num_objects, sizeof(struct object *), &void_ptr_compare);

  return 0;
}

/*************************************************************************
collect_species:
    Collect all objects being visualized.

        In:  struct species ***vol_species - pointer to receive array of
                                             species *s for vol mols
             int *n_vol_species - int to receive vol mol species count
             struct species ***grid_species - pointer to receive array of
                                              species *s for grid mols
             int *n_grid_species - int to receive grid mol species count
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int collect_species(struct species ***vol_species,
                           int *n_vol_species,
                           struct species ***grid_species,
                           int *n_grid_species)
{
  int vcount=0, gcount=0;

  /* Count species */
  int spec_id;
  for (spec_id = 0; spec_id < world->n_species; ++ spec_id)
  {
    struct species *specp = world->species_list[spec_id];
    if ((specp->flags & IS_SURFACE) != 0) continue;
    if (strcmp(specp->sym->name, "GENERIC_MOLECULE") == 0) continue;
    if (specp->viz_state <= 0) continue;

    if (((specp->flags & NOT_FREE) != 0))
      ++ gcount;
    else
      ++ vcount;
  }

  /* Allocate arrays */
  if ((*vol_species = (struct species **) allocate_ptr_array(vcount)) == NULL)
    goto failure;
  if ((*grid_species = (struct species **) allocate_ptr_array(gcount)) == NULL)
    goto failure;
  *n_vol_species = vcount;
  *n_grid_species = gcount;

  /* Sort species into arrays */
  struct species **vol_cur = *vol_species, **grid_cur = *grid_species;
  for (spec_id = 0;
       spec_id < world->n_species  &&  (vcount || gcount);
       ++ spec_id)
  {
    struct species *specp = world->species_list[spec_id];
    if ((specp->flags & IS_SURFACE) != 0) continue;
    if (strcmp(specp->sym->name, "GENERIC_MOLECULE") == 0) continue;
    if (specp->viz_state <= 0) continue;

    if (((specp->flags & NOT_FREE) != 0))
    {
      *(grid_cur++) = specp;
      -- gcount;
    }
    else
    {
      *(vol_cur++) = specp;
      -- vcount;
    }
  }

  return 0;

failure:
  if (*vol_species) free(*vol_species);
  if (*grid_species) free(*grid_species);
  *vol_species = NULL;
  *grid_species = NULL;
  *n_vol_species = *n_grid_species = -1;
  return 1;
}

/*************************************************************************
convert_frame_data_to_iterations:
    Converts time-indexed frame to an iteration-indexed frame.

        In:  struct frame_data_list *fdlp - the frame to convert
        Out: 0 if successful
**************************************************************************/
static int convert_frame_data_to_iterations(struct frame_data_list *fdlp)
{
  struct num_expr_list *nel;
  for (nel = fdlp->iteration_list; nel != NULL; nel = nel->next)
    nel->value = (double) (long long)(nel->value / world->time_unit + ROUND_UP);
  fdlp->list_type = OUTPUT_BY_ITERATION_LIST;
  return 0;
}

/* == DX/DREAMM Utilities == */

/*************************************************************************
dx_output_worldfloat
    Writes a floating point world coordinate to the specified file.

        In:  FILE *f - file handle to receive output
             float fval - floating point value to write
        Out: num bytes written
**************************************************************************/
static int dx_output_worldfloat(FILE *f, double fval)
{
  float v = (float) (fval * world->length_unit);
  fwrite(&v, sizeof(v), 1, f);
  return sizeof(float);
}

/*************************************************************************
dx_output_vector3
    Writes a 3d vector to a DX output file in binary format.

        In:  FILE *f - file handle to receive output
             struct vector3 *v3 - vector to write
        Out: num bytes written
**************************************************************************/
static int dx_output_vector3(FILE *f, struct vector3 const *v3)
{
  int size = 0;
  size += dx_output_worldfloat(f, v3->x);
  size += dx_output_worldfloat(f, v3->y);
  size += dx_output_worldfloat(f, v3->z);
  return size;
}

/*************************************************************************
dx_output_oriented_normal
    Writes a 3d oriented normal vector to a DX output file in binary format.

        In:  FILE *f - file handle to receive output
             struct vector3 *v3 - vector to write
             short orient - the orientation of the normal (+1/-1)
        Out: num bytes written
**************************************************************************/
static int dx_output_oriented_normal(FILE *f, struct vector3 const *v3, short orient)
{
  float f1 = (float) (v3->x * (float) orient);
  float f2 = (float) (v3->y * (float) orient);
  float f3 = (float) (v3->z * (float) orient);
  fwrite(&f1, sizeof(f1), 1, f);
  fwrite(&f2, sizeof(f2), 1, f);
  fwrite(&f3, sizeof(f3), 1, f);
  return 3*sizeof(float);
}

/*************************************************************************
dx_output_vertices
    Writes the vertices for an object out to a DX output file in binary format.

        In:  FILE *f - file handle to receive output
             struct object *objp - the object whose vertices to write
        Out: num bytes written
**************************************************************************/
static int dx_output_vertices(FILE *f, struct object const *objp)
{
  int size = 0;
  int i;
  for (i=0; i<objp->n_verts; ++i)
    size += dx_output_vector3(f, &objp->verts[i]);
  return size;
}

/*************************************************************************
dx_output_wall_vertices
    Writes the vertex indices for a wall to a DX output file in binary format.

    XXX: Is this right?  sizeof(int) can vary between platforms...

        In:  FILE *f - file handle to receive output
             struct element_data *edp - the wall vertices to write
        Out: num bytes written
**************************************************************************/
static int dx_output_wall_vertices(FILE *f, struct element_data const *edp)
{
  int i;
  for (i=0; i<3; ++i)
  {
    int vi = edp->vertex_index[i];
    fwrite(&vi, sizeof(int), 1, f);
  }
  return 3*sizeof(int);
}

/*************************************************************************
dx_output_gridpt_and_normal
    Writes the point and oriented normal for a location on a grid.

        In:  FILE *f - file handle to receive output
             struct wall *w - the wall of the surface grid
             int grid_index - the grid triangle index
             short orient - the orientation (+1/-1)
        Out: num bytes written
**************************************************************************/
static int dx_output_gridpt_and_normal(FILE *f,
                                       struct wall *w,
                                       int grid_index,
                                       short orient)
{
  int size = 0;
  struct vector3 p0;

  grid2xyz(w->grid, grid_index, &p0);
  size += dx_output_vector3(f, &p0);
  size += dx_output_oriented_normal(f, &w->normal, orient);

  return size;
}

/* == DX output == */

/*************************************************************************
dx_open_file:
    Opens a DX output file.

        In:  char const *cls  - class of file (second part of filename)
             char const *prefix - prefix of file (first part of filename)
             long long iteration - iteration number
        Out: file handle on success, NULL on failure
**************************************************************************/
static FILE *dx_open_file(char const *cls,
                          char const *prefix,
                          long long iteration)
{
  char *filename_buffer = alloc_sprintf("%s.%s.%lld.dx",
                                        prefix,
                                        cls,
                                        iteration);
  if (filename_buffer == NULL)
  {
    fprintf(world->err_file,"File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return NULL;
  }

  FILE *f = open_file(filename_buffer, "wb");
  free(filename_buffer);
  return f;
}

/*************************************************************************
dx_output_walls
    Outputs all wall data for a given object to the specified files.

        In:  FILE *wall_verts_header
             FILE *wall_states_header
             struct object *objp
        Out: 0 on success, 1 on error; on success, the objects are written to
             the files
**************************************************************************/
static int dx_output_walls(FILE *wall_verts_header,
                           FILE *wall_states_header,
                           struct object const *objp)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";

  struct polygon_object *pop = (struct polygon_object *) objp->contents;
  struct ordered_poly *opp = (struct ordered_poly *) pop->polygon_data;
  struct element_data *edp = opp->element;
  int element_data_count = objp->n_walls_actual;
  int wall_index;

  if (wall_verts_header)
  {
    fprintf(wall_verts_header,
            "object \"%s.positions\" class array type float rank 1 shape 3 items %d %s binary data follows\n",
            objp->sym->name,
            objp->n_verts,
            my_byte_order);

    /* output polyhedron vertices */
    (void) dx_output_vertices(wall_verts_header, objp);
    fprintf(wall_verts_header,
            "\nattribute \"dep\" string \"positions\"\n#\n");

    fprintf(wall_verts_header,
            "object \"%s.connections\" class array type int rank 1 shape 3 items %d %s binary data follows\n",
            objp->sym->name, element_data_count, my_byte_order);
  }

  if (wall_states_header)
    fprintf(wall_states_header,
            "object \"%s.states\" class array type int rank 0 items %d ascii data follows\n",
            objp->sym->name, element_data_count);

  /* output polygon element connections */
  for (wall_index = 0; wall_index < objp->n_walls; ++ wall_index)
  {
    if (get_bit(pop->side_removed,wall_index))
      continue;

    if (wall_verts_header)
      dx_output_wall_vertices(wall_verts_header, &edp[wall_index]);

    if (wall_states_header)
      fprintf(wall_states_header, "%d\n", objp->viz_state[wall_index]);
  }

  if (wall_verts_header)
  {
    fprintf(wall_verts_header, "\nattribute \"ref\" string \"positions\"\n");
    fprintf(wall_verts_header, "attribute \"element type\" string \"triangles\"\n#\n");

    /* XXX: The following used to be wall_states_header, but I don't think that
     * was right...
     */
    fprintf(wall_verts_header,"\nattribute \"dep\" string \"connections\"\n#\n");
  }

  if (wall_states_header)
    fprintf(wall_states_header,"\nattribute \"dep\" string \"states\"\n#\n");

  return 0;
}

/*************************************************************************
dx_output_count_effectors
    Counts the total number of (non-excluded) effectors in a given object.

        In:  struct object *objp
        Out: 0 on success, 1 on error; on success, the objects are written to
             the files
**************************************************************************/
static int dx_output_count_effectors(struct object *objp)
{
  int wall_index;
  int n_eff = 0;
  for (wall_index = 0; wall_index < objp->n_walls; ++ wall_index)
  {
    int nremain;
    unsigned int tile_index;
    struct wall *w = objp->wall_p[wall_index];
    if (w == NULL  ||  w->grid == NULL)
      continue;

    nremain = w->grid->n_occupied;
    for (tile_index = 0; tile_index < w->grid->n_tiles  &&  nremain != 0; ++ tile_index)
    {
      if (w->grid->mol[tile_index] == NULL)
        continue;

      -- nremain;
      if (w->grid->mol[tile_index]->properties->viz_state != EXCLUDE_OBJ)
        ++ n_eff;
    }
  }

  return n_eff;
}

/*************************************************************************
dx_output_effectors_on_wall
    Writes effector data for a specific wall to the specified output files.

        In:  FILE *eff_pos_header
             FILE *eff_states_header
             struct wall *w
        Out: 0 on success, 1 on error; on success, the objects are written to
             the files
**************************************************************************/
static int dx_output_effectors_on_wall(FILE *eff_pos_header,
                                       FILE *eff_states_header,
                                       struct wall *w)
{
  int nremain;
  unsigned int tile_index;

  nremain = w->grid->n_occupied;
  for (tile_index = 0; tile_index<w->grid->n_tiles && nremain != 0; ++ tile_index)
  {
    struct grid_molecule *gmol=w->grid->mol[tile_index];
    if (gmol == NULL)
      continue;

    if (gmol->properties->viz_state == EXCLUDE_OBJ)
      continue;

    -- nremain;

    if (eff_states_header)
      fprintf(eff_states_header, "%d\n", gmol->properties->viz_state);

    if (eff_pos_header)
      (void) dx_output_gridpt_and_normal(eff_pos_header,
                                         w,
                                         tile_index,
                                         gmol->orient);
  }

  return 0;
}

/*************************************************************************
dx_output_effectors
    writes effector data to effector position/state files

        In:  FILE *eff_pos_header
             FILE *eff_states_header
             struct object *objp
        out: 0 on success, 1 on error; on success, the objects are written to
             the files
**************************************************************************/
static int dx_output_effectors(FILE *eff_pos_header,
                               FILE *eff_states_header,
                               struct object *objp)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";
  int n_eff=0;
  int wall_index;

  no_printf("Traversing walls in object %s\n", objp->sym->name);

  n_eff = dx_output_count_effectors(objp);

  no_printf("Dumping %d effectors...\n", n_eff);
  fflush(world->log_file);

  if (eff_pos_header) {
    if (n_eff) {
      fprintf(eff_pos_header,
              "object \"%s.pos_and_norm\" class array type float rank 2 shape 2 3 items %d %s binary data follows\n",
              objp->sym->name,
              n_eff,
              my_byte_order);
    }
    else {
      fprintf(eff_pos_header,
              "object \"%s.pos_and_norm\" array",
              objp->sym->name);
    }
  }

  if (eff_states_header) {
    if (n_eff) {
      fprintf(eff_states_header,
              "object \"%s.states\" class array type int rank 0 items %d ascii data follows\n",
              objp->sym->name,
              n_eff);
    }
    else {
      fprintf(eff_states_header,
              "object \"%s.states\" array\n",
              objp->sym->name);
    }
  }

  /* dump the effectors */
  for (wall_index = 0; wall_index < objp->n_walls; ++ wall_index)
  {
    struct wall *w = objp->wall_p[wall_index];
    if (w == NULL)
      continue;

    if (w->grid == NULL)
      continue;

    if (dx_output_effectors_on_wall(eff_pos_header, eff_states_header, w))
      return 1;
  }

  if (eff_pos_header)
    fprintf(eff_pos_header, "\n#\n");

  if (eff_states_header)
    fprintf(eff_states_header,
            "attribute \"dep\" string \"positions\"\n#\n");
  return 0;
}

/*************************************************************************
dx_output_walls_groups
    Writes group objects for walls

        In:  FILE *wall_verts_header
             FILE *wall_states_header
             struct viz_obj *vizp
        Out: 0 on success, 1 on error; on success, the objects are written to
             the files
**************************************************************************/
static int dx_output_walls_groups(FILE *wall_verts_header,
                                  FILE *wall_states_header,
                                  struct viz_obj *vizp)
{
  struct viz_child *vcp;

  /* output surface positions null objects and group objects */
  if (wall_verts_header) {
    fprintf(wall_verts_header, "object \"%s\" group\n", vizp->full_name);
    fprintf(wall_verts_header, "  member \"null_object (default)\" \"null_object\"\n");
  }

  /* surface states */
  if (wall_states_header) {
    fprintf(wall_states_header, "object \"%s\" group\n", vizp->full_name);
    fprintf(wall_states_header, "  member \"null_object (default)\" \"null_object\"\n");
  }

  /* output group object members */
  /* member name is full MCell child object name */
  for (vcp = vizp->viz_child_head;
       vcp != NULL;
       vcp = vcp->next)
  {
    struct object *objp = vcp->obj;

    /* surface positions */
    if (wall_verts_header) {
      fprintf(wall_verts_header,
              "  member \"%s\" \"%s.field\"\n",
              objp->sym->name,
              objp->sym->name);
    }

    /* surface states */
    if (wall_states_header) {
      fprintf(wall_states_header,
              "  member \"%s.states\" \"%s.states\"\n",
              objp->sym->name,
              objp->sym->name);
    }
  }

  return 0;
}

/*************************************************************************
dx_output_effectors_groups
    Writes group objects for effectors

        In:  FILE *eff_pos_header
             FILE *eff_states_header
             struct viz_obj *vizp
        Out: 0 on success, 1 on error; on success, the objects are written to
             the files
**************************************************************************/
static int dx_output_effectors_groups(FILE *eff_pos_header,
                                      FILE *eff_states_header,
                                      struct viz_obj *vizp)
{
  struct viz_child *vcp;
  /* output effector positions null objects and group objects */
  if (eff_pos_header) {
    fprintf(eff_pos_header, "object \"%s\" group\n", vizp->full_name);
    fprintf(eff_pos_header, "  member \"null_object (default)\" \"null_object\"\n");
  }

  /* effector states */
  if (eff_states_header) {
    fprintf(eff_states_header, "object \"%s\" group\n", vizp->full_name);
    fprintf(eff_states_header, "  member \"null_object (default)\" \"null_object\"\n");
  }

  /* output group object members */
  /* member name is full MCell child object name */
  for (vcp = vizp->viz_child_head;
       vcp != NULL;
       vcp = vcp->next)
  {
    struct object *objp = vcp->obj;

    /* effector positions */
    if (eff_pos_header) {
      fprintf(eff_pos_header,
              "  member \"%s\" \"%s.field\"\n",
              objp->sym->name,
              objp->sym->name);
    }

    /* effector states */
    if (eff_states_header) {
      fprintf(eff_states_header,
              "  member \"%s.states\" \"%s.states\"\n",
              objp->sym->name,
              objp->sym->name);
    }
  }

  return 0;
}

/*************************************************************************
dx_output_walls_null
    Writes null objects for walls

        In:  FILE *wall_verts_header
             FILE *wall_states_header
        Out: 0 on success, 1 on error; on success, the objects are written to
             the files
**************************************************************************/
static void dx_output_walls_null(FILE *wall_verts_header, FILE *wall_states_header)
{
  if (wall_verts_header)
  {
    fprintf(wall_verts_header, "object \"null_positions\" array\n\n");
    fprintf(wall_verts_header, "object \"null_connections\" array\n\n");
    fprintf(wall_verts_header, "object \"null_object\" field\n");
    fprintf(wall_verts_header, "  component \"positions\" \"null_positions\"\n");
    fprintf(wall_verts_header, "  component \"connections\" \"null_connections\"\n\n");
  }

  if (wall_states_header)
  {
    fprintf(wall_states_header, "object \"null_object\" array\n");
    fprintf(wall_states_header, "  attribute \"dep\" string \"connections\"\n\n");
  }
}

/*************************************************************************
dx_output_effectors_null
    Writes null objects for effectors

        In:  FILE *eff_pos_header
             FILE *eff_states_header
        Out: 0 on success, 1 on error; on success, objects are written to the
             files
**************************************************************************/
static void dx_output_effectors_null(FILE *eff_pos_header, FILE *eff_states_header)
{
  if (eff_pos_header)
  {
    fprintf(eff_pos_header, "object \"null_pos_and_norm\" array\n\n");
    fprintf(eff_pos_header, "object \"null_object\" field\n");
    fprintf(eff_pos_header, "  component \"data\" \"null_pos_and_norm\"\n\n");
  }

  if (eff_states_header)
  {
    fprintf(eff_states_header, "object \"null_object\" array\n");
    fprintf(eff_states_header, "  attribute \"dep\" string \"positions\"\n\n");
  }
}

/*************************************************************************
dx_output_walls_and_effectors_fields
    Writes field descriptors for wall/effector DX output files

        In:  FILE *wall_verts_header
             FILE *eff_pos_header
             struct viz_child *vcp
        Out: 0 on success, 1 on error; on success, the objects are written to
             the file
**************************************************************************/
static void dx_output_walls_and_effectors_fields(FILE *wall_verts_header,
                                                 FILE *eff_pos_header,
                                                 struct viz_child *vcp)
{
  while (vcp != NULL)
  {
    struct object *objp = vcp->obj;

    /* effector positions */
    if (eff_pos_header) {
      fprintf(eff_pos_header, "object \"%s.field\" field\n", objp->sym->name);
      fprintf(eff_pos_header, "  component \"data\" \"%s.pos_and_norm\"\n\n", objp->sym->name);
    }

    /* surface positions */
    if (wall_verts_header) {
      fprintf(wall_verts_header, "object \"%s.field\" field\n", objp->sym->name);
      fprintf(wall_verts_header, "  component \"positions\" \"%s.positions\"\n", objp->sym->name);
      fprintf(wall_verts_header, "  component \"connections\" \"%s.connections\"\n\n", objp->sym->name);
    }

    vcp = vcp->next;
  }
}

/*************************************************************************
dx_output_walls_and_effectors_single
    Writes a single set of wall/effector DX output files

        In:  struct frame_data_list *fdlp
             struct viz_obj *vizp
        Out: 0 on success, 1 on error; on success, all output files are created
**************************************************************************/
static int dx_output_walls_and_effectors_single(struct frame_data_list *fdlp,
                                                struct viz_obj *vizp)
{
  int viz_surf_pos     = ((fdlp->type == ALL_FRAME_DATA) || (fdlp->type == SURF_POS));
  int viz_surf_states  = ((fdlp->type == ALL_FRAME_DATA) || (fdlp->type == SURF_STATES));
  int viz_eff_pos      = ((fdlp->type == ALL_FRAME_DATA) || (fdlp->type == EFF_POS));
  int viz_eff_states   = ((fdlp->type == ALL_FRAME_DATA) || (fdlp->type == EFF_STATES));

  FILE *wall_verts_header = NULL;
  FILE *wall_states_header = NULL;
  FILE *eff_pos_header = NULL;
  FILE *eff_states_header = NULL;

  struct viz_child *vcp = NULL;

  /* XXX: May leave file handles open on failure */

  if (viz_surf_pos  &&
      (wall_verts_header = dx_open_file("mesh_elements", vizp->name, fdlp->viz_iteration)) == NULL)
    goto failure;

  if (viz_surf_states  &&
      (wall_states_header = dx_open_file("mesh_element_states", vizp->name, fdlp->viz_iteration)) == NULL)
    goto failure;

  if (viz_eff_pos  &&
      (eff_pos_header = dx_open_file("effector_site_positions", vizp->name, fdlp->viz_iteration)) == NULL)
    goto failure;

  if (viz_eff_states  &&
      (eff_states_header = dx_open_file("effector_site_states", vizp->name, fdlp->viz_iteration)) == NULL)
    goto failure;

  for (vcp = vizp->viz_child_head;
       vcp != NULL;
       vcp = vcp->next)
  {
    struct object *objp = vcp->obj;

    if ((viz_surf_pos || viz_surf_states)  &&
        (objp->object_type==POLY_OBJ || objp->object_type==BOX_OBJ))
    {
      if (dx_output_walls(wall_verts_header,
                          wall_states_header,
                          objp))
        goto failure;
    }

    if (viz_eff_pos  ||  viz_eff_states) {
      if (dx_output_effectors(eff_pos_header,
                              eff_states_header,
                              objp))
        goto failure;
    }
  }

  /* Output null objects */
  if (viz_surf_states || viz_surf_pos)
    dx_output_walls_null(wall_verts_header, wall_states_header);
  if (viz_eff_states || viz_eff_pos)
    dx_output_effectors_null(eff_pos_header, eff_states_header);

  /* output effector and surface positions field objects */
  if (viz_eff_pos  ||  viz_surf_pos)
  {
    dx_output_walls_and_effectors_fields(wall_verts_header,
                                         eff_pos_header,
                                         vizp->viz_child_head);
  }

  if (viz_eff_pos  ||  viz_eff_states)
    if (dx_output_effectors_groups(eff_pos_header, eff_states_header, vizp))
      goto failure;

  if (viz_surf_pos  ||  viz_surf_states)
    if (dx_output_walls_groups(wall_verts_header, wall_states_header, vizp))
      goto failure;

  if (wall_verts_header != NULL)  fclose(wall_verts_header);
  if (wall_states_header != NULL) fclose(wall_states_header);
  if (eff_pos_header != NULL)     fclose(eff_pos_header);
  if (eff_states_header != NULL)  fclose(eff_states_header);
  return 0;

failure:
  if (wall_verts_header != NULL)  fclose(wall_verts_header);
  if (wall_states_header != NULL) fclose(wall_states_header);
  if (eff_pos_header != NULL)     fclose(eff_pos_header);
  if (eff_states_header != NULL)  fclose(eff_states_header);
  return 1;
}

/*************************************************************************
dx_output_walls_and_effectors
    Writes wall and/or effector output files for the old DX output format.

        In:  struct frame_data_list *fdlp
        Out: 0 on success, 1 on error; on success, all output files are created
**************************************************************************/
static int dx_output_walls_and_effectors(struct frame_data_list *fdlp)
{
  struct viz_obj *vizp;
  for (vizp = world->viz_obj_head;
       vizp != NULL;
       vizp = vizp->next)
  {
    if (dx_output_walls_and_effectors_single(fdlp, vizp))
      return 1;
  }

  return 0;
}

/*************************************************************************
dx_output_molecules_position_fields_and_groups
    Writes field and group information for molecule positions output file in
    the old DX output format.

        In:  FILE *mol_pos_header
             int   mol_pos_index
        Out: 0 on success, 1 on error; data is appended to the provided file
             handle.
**************************************************************************/
static int dx_output_molecules_position_fields_and_groups(FILE *mol_pos_header, int mol_pos_index)
{
  int mol_pos_field_index = 0;
  int mol_pos_group_index = 0;
  int species_index;

  for (species_index = 0; species_index < world->n_species; ++ species_index)
  {
    if ((world->species_list[species_index]->flags & NOT_FREE)!=0)
      continue;

    fprintf(mol_pos_header,
            "object \"%d\" field\n",mol_pos_index+mol_pos_field_index);
    fprintf(mol_pos_header,
            "  component \"positions\" \"%d\"\n\n",1+mol_pos_field_index);
    ++ mol_pos_field_index;
  }
  fprintf(mol_pos_header, "object \"null_positions\" array\n\n");
  fprintf(mol_pos_header, "object \"null_object\" field\n");
  fprintf(mol_pos_header, "  component \"positions\" \"null_positions\"\n\n");
  fprintf(mol_pos_header, "object \"%d\" group\n", 2*mol_pos_index-1);
  fprintf(mol_pos_header, "  member \"null_object (default)\" \"null_object\"\n");

  for (species_index = 0; species_index < world->n_species; ++ species_index)
  {
    if ((world->species_list[species_index]->flags & NOT_FREE)!=0)
      continue;

    fprintf(mol_pos_header, "  member \"%s\" \"%d\"\n",
                            world->species_list[species_index]->sym->name,
                            mol_pos_index+mol_pos_group_index);
    ++ mol_pos_group_index;
  }

  return 0;
}

/*************************************************************************
dx_output_molecules_states_groups:
    Writes group information for molecule state output file in the old DX
    output format.

        In:  FILE *mol_states_header
             int   mol_states_index
        Out: 0 on success, 1 on error; data is appended to the provided file
             handle.
**************************************************************************/
static int dx_output_molecules_states_groups(FILE *mol_states_header, int mol_states_index)
{
  int mol_states_group_index = 0;
  int species_index;

  fprintf(mol_states_header, "object \"null_object\" array\n");
  fprintf(mol_states_header, "  attribute \"dep\" string \"positions\"\n\n");
  fprintf(mol_states_header, "object \"%d\" group\n", mol_states_index);
  fprintf(mol_states_header, "  member \"null_object (default)\" \"null_object\"\n");

  for (species_index = 0; species_index < world->n_species; ++ species_index)
  {
    if ((world->species_list[species_index]->flags & NOT_FREE) != 0)
      continue;

    fprintf(mol_states_header, "  member \"%s\" \"%d\"\n",
                               world->species_list[species_index]->sym->name,
                               mol_states_group_index + 1);
    ++ mol_states_group_index;
  }

  return 0;
}

/*************************************************************************
dx_output_molecules_position:
    Writes molecule positions to the molecule position DX file.

        In:  FILE *mol_pos_header
             struct volume_molecule **viz_molp
             u_int mol_count
             int   *mol_pos_index
        Out: 0 on success, 1 on error; data is appended to the provided file
             handle, and mol_states_index is incremented.
**************************************************************************/
static int dx_output_molecules_position(FILE *mol_pos_header, struct volume_molecule **viz_molp, u_int mol_count, int *mol_pos_index)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";

  if (mol_count > 0)
  {
    unsigned int mol_index;
    fprintf(mol_pos_header, "object \"%d\" class array type float rank 1 shape 3 items %d %s binary data follows\n",
            *mol_pos_index,
            mol_count,
            my_byte_order);

    for (mol_index = 0; mol_index<mol_count; ++ mol_index)
    {
      struct volume_molecule *molp = viz_molp[mol_index];
      dx_output_vector3(mol_pos_header, &molp->pos);
    }
    fprintf(mol_pos_header, "\n#\n");
  }
  else /* mol_count == 0 */
    fprintf(mol_pos_header, "object \"%d\" array\n\n", *mol_pos_index);
  ++ *mol_pos_index;
  return 0;
}

/*************************************************************************
dx_output_molecules_state:
    Writes state info to the molecule state DX file.

        In:  FILE *mol_states_header
             u_int mol_count
             int   state
             int   *mol_states_index
        Out: 0 on success, 1 on error; data is appended to the provided file
             handle, and mol_states_index is incremented.
**************************************************************************/
static int dx_output_molecules_state(FILE *mol_states_header, u_int mol_count, int state, int *mol_states_index)
{
  if (mol_count > 0)
  {
    fprintf(mol_states_header, "object \"%d\"\n  constantarray type int items %d\n",
                               *mol_states_index,
                               mol_count);
    fprintf(mol_states_header, "  data follows %d\n",
                               state);
  }
  else /* mol_count == 0 */
    fprintf(mol_states_header, "object \"%d\" array\n",
                               *mol_states_index);
  fprintf(mol_states_header, "  attribute \"dep\" string \"positions\"\n\n");
  ++ *mol_states_index;
  return 0;
}

/*************************************************************************
dx_output_molecules:
    Write out molecules for the old DX output format.

        In: struct frame_data_list *fdlp
        Out: 0 on success, 1 on error; output visualization files (*.dx)
             are written.
**************************************************************************/
static int dx_output_molecules(struct frame_data_list *fdlp)
{
  int retcode = 0;
  FILE *mol_pos_header = NULL;
  FILE *mol_states_header = NULL;
  byte viz_mol_pos=((fdlp->type==ALL_FRAME_DATA) || (fdlp->type==MOL_POS));
  byte viz_mol_states=((fdlp->type==ALL_FRAME_DATA) || (fdlp->type==MOL_STATES));
  int mol_pos_index=1;
  int mol_states_index=1;
  struct volume_molecule ***viz_molp = NULL;
  u_int *viz_mol_count = NULL;
  int species_index;

  /* XXX: May leave file handles open on failure */
  if (viz_mol_pos  &&
      (mol_pos_header = dx_open_file("molecule_positions",
                                     world->molecule_prefix_name,
                                     fdlp->viz_iteration)) == NULL)
    goto failure;

  if (viz_mol_states  &&
      (mol_states_header = dx_open_file("molecule_states",
                                        world->molecule_prefix_name,
                                        fdlp->viz_iteration)) == NULL)
      goto failure;

  if (sort_molecules_by_species((struct abstract_molecule ****) (void *) &viz_molp,
                                &viz_mol_count,
                                1, 0))
    goto failure;

  /* Iterate over species */
  for (species_index = 0; species_index<world->n_species; ++ species_index)
  {
    u_int spec_id = world->species_list[species_index]->species_id;
    int mol_count = viz_mol_count[spec_id];
    int state = world->species_list[species_index]->viz_state;
    if ((world->species_list[species_index]->flags & NOT_FREE) != 0)
      continue;

    if (viz_mol_pos  &&  dx_output_molecules_position(mol_pos_header,
                                                      viz_molp[spec_id],
                                                      mol_count,
                                                      &mol_pos_index))
      goto failure;

    if (viz_mol_states  &&  dx_output_molecules_state(mol_states_header,
                                                      mol_count,
                                                      state,
                                                      &mol_states_index))
      goto failure;
  }

  /* build fields and groups here */
  if (viz_mol_pos  &&  dx_output_molecules_position_fields_and_groups(mol_pos_header, mol_pos_index))
    goto failure;

  if (viz_mol_states  &&  dx_output_molecules_states_groups(mol_states_header, mol_states_index))
    goto failure;

  goto success;

failure:
  retcode = 1;

success:
  /* clean up */
  if (mol_states_header)
    fclose(mol_states_header);
  mol_states_header = NULL;

  if (mol_pos_header)
    fclose(mol_pos_header);
  mol_pos_header = NULL;

  if (viz_molp != NULL)
    free_ptr_array((void **) viz_molp, world->n_species);
  viz_molp = NULL;

  if (viz_mol_count != NULL)
    free (viz_mol_count);
  viz_mol_count = NULL;

  return retcode;
}

/*************************************************************************
output_dx_objects:
    Write out a frame of data in the old DX output format.

        In: struct frame_data_list *fdlp
        Out: 0 on success, 1 on error; output visualization files (*.dx)
             are written.
**************************************************************************/
int output_dx_objects(struct frame_data_list *fdlp)
{
  int viz_type;
  byte viz_eff, viz_mol, viz_surf;

  no_printf("Viz output in DX mode...\n");

  viz_type = fdlp->type;
  viz_eff  = ((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)  || (viz_type==EFF_STATES));
  viz_mol  = ((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)  || (viz_type==MOL_STATES));
  viz_surf = ((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS) || (viz_type==SURF_STATES));

  /* dump walls and effectors: */
  if (viz_surf  ||  viz_eff)
    if (dx_output_walls_and_effectors(fdlp))
      return 1;

  /* dump diffusible molecules: */
  if (viz_mol)
    if (dx_output_molecules(fdlp))
      return 1;

  return 0;
}

static const struct vector3 v3_unit_z = { 0.0, 0.0, 1.0 };

/* == DREAMM Generic Utilities == */

/*************************************************************************
dreamm_v3_generic_open_file:
    Open a file in a given directory, being careful to remove symlinks before
    doing so.

        In: char const *dir - directory name for file
            char const *fname - filename for file
            char const *mode - file access mode, as to fopen
        Out: file handle for file, NULL on error
**************************************************************************/
static FILE *dreamm_v3_generic_open_file(char const *dir, char const *fname, char const *mode)
{
  struct stat f_stat;
  FILE *f;

  /* concatenate dir and fname to get the mesh states file path */
  char *path = NULL;
  if (dir != NULL)
    path = alloc_sprintf("%s/%s", dir, fname);
  else
    path = strdup(fname);
  if (path == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    goto failure;
  }

  /* If the file exists and is a symlink, remove it */
  if (stat(path, &f_stat) == 0  &&
      (f_stat.st_mode & S_IFLNK) == S_IFLNK)
  {
    /* remove the symbolic link */
    if (unlink(path) != 0)
    {
      fprintf(world->err_file, "File %s, Line %ld: error %d removing file %s.\n", __FILE__, (long)__LINE__, errno, fname);
      goto failure;
    }
  }

  f = open_file(path, mode);
  free(path);
  return f;

failure:
  if (path) free(path);
  return NULL;
}

/*************************************************************************
dreamm_v3_generic_init:
    Initialize state which is common to both forms of DREAMM output.

        In:  Nothing
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_generic_init()
{
  /* Break file prefix name into basename and dirname */
  if (get_basename(world->file_prefix_name,
                   &world->viz_state_info.filename_prefix_basename))
    return 1;
  if (get_dirname(world->file_prefix_name,
                  &world->viz_state_info.filename_prefix_dirname))
    return 1;

  /* test whether a directory structure created by the user exists */
  if (world->viz_state_info.filename_prefix_dirname)
  {
    switch (dir_exists(world->viz_state_info.filename_prefix_dirname))
    {
      case 0:
        if (mkdirs(world->viz_state_info.filename_prefix_dirname, world->err_file))
        {
          fprintf(world->err_file, "File %s, Line %ld: viz_output data directory '%s' could not be created.\n", __FILE__, (long)__LINE__, world->viz_state_info.filename_prefix_dirname);
          return 1;
        }
        break;

      case -1:
        return 1;

      case 1:
        /* All's well */
        break;

      default:
        fprintf(world->err_file, "File %s, Line %ld: unexpected return value from dir_exists (internal error).\n", __FILE__, (long)__LINE__);
        return 1;
    }
  }

  /* Collect all mesh objects to be visualized */
  if (collect_objects(world->viz_obj_head,
                      &world->viz_state_info.viz_objects,
                      &world->viz_state_info.n_viz_objects))
    return 1;

  /* Collect all species to be visualized */
  if (collect_species(&world->viz_state_info.vol_species,
                      &world->viz_state_info.n_vol_species,
                      &world->viz_state_info.grid_species,
                      &world->viz_state_info.n_grid_species))
    return 1;

  /*
   * Check here if MESHES or MOLECULES blocks are not supplied.  The
   * corresponding files should be empty in order to prevent unintentional
   * mixing of pre-existing and new files.
   */
  if (world->viz_state_info.n_viz_objects == 0)
    fprintf(world->log_file, "\nMESHES keyword is absent or commented.\nEmpty 'meshes' output files are created.\n\n");
  if ((world->viz_state_info.n_vol_species == 0) && (world->viz_state_info.n_grid_species == 0))
    fprintf(world->log_file, "\nMOLECULES keyword is absent or commented.\nEmpty 'molecules' output files are created.\n\n");

  return 0;
}

/*************************************************************************
dreamm_v3_generic_merge_frame_data:
    Merges the data frame 'src' into 'dest'.  No checking is done, so it's
    important to make sure that src and dest are the same type of frame before
    calling this.  It's used before DREAMM output to merge frames of the same
    type so that we can convert paired MOL_POS and MOL_ORIENT frames into
    ALL_MOL_DATA frames.  The source frame is not freed, but its iteration list
    will be NULL after this method is called.

        In:  struct frame_data_list *dest - destination frame
             struct frame_data_list *src - source frame
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_generic_merge_frame_data(struct frame_data_list *dest,
                                              struct frame_data_list *src)
{
  struct num_expr_list *nelSrc = src->iteration_list,
                       *nelDest = dest->iteration_list,
                       *outTail = NULL;
  src->iteration_list = NULL;
  dest->iteration_list = NULL;

  while (nelSrc != NULL  &&  nelDest != NULL)
  {
    if (nelSrc->value < nelDest->value)
    {
      if (outTail)
        outTail->next = nelSrc;
      else
        dest->iteration_list = nelSrc;
      outTail = nelSrc;
      nelSrc = nelSrc->next;
      outTail->next = NULL;
    }
    else
    {
      if (outTail)
        outTail->next = nelDest;
      else
        dest->iteration_list = nelDest;
      outTail = nelDest;
      nelDest = nelDest->next;
      outTail->next = NULL;
    }
  }

  if (nelSrc)
  {
    if (outTail)
      outTail->next = nelSrc;
    else
      dest->iteration_list = nelSrc;
  }
  else if (nelDest)
  {
    if (outTail)
      outTail->next = nelSrc;
    else
      dest->iteration_list = nelSrc;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_generic_merge_coincident_frames:
    Merges coincident frames from two lists into a third list.  This is used to
    turn MOL_POS+MOL_ORIENT into ALL_MOL_DATA, or MESH_GEOMETRY+REG_DATA into
    ALL_MESH_DATA.  This greatly simplifies DREAMM processing.

        In:  struct frame_data_list **both - target list for merged frames
             struct frame_data_list **first - first list to merge
             struct frame_data_list **second - second list to merge
             struct frame_data_list **discard - any empty frames are put here
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_generic_merge_coincident_frames(struct frame_data_list **both,
                                                     struct frame_data_list **first,
                                                     struct frame_data_list **second,
                                                     struct frame_data_list **discard)
{
  struct num_expr_list **nelBoth, **nelFirst, **nelSecond;
  nelBoth = &(*both)->iteration_list;
  nelFirst = &(*first)->iteration_list;
  nelSecond = &(*second)->iteration_list;

  while (*nelFirst != NULL  &&  *nelSecond != NULL)
  {
    if (fabs((*nelFirst)->value - (*nelSecond)->value) < EPS_C)
    {
      struct num_expr_list *temp;

      /* Find insertion point */
      while ((*nelBoth)  &&  (*nelBoth)->value < (*nelFirst)->value)
        nelBoth = &(*nelBoth)->next;

      /* Splice time into "both" list */
      temp = *nelFirst;
      *nelFirst = (*nelFirst)->next;
      temp->next = *nelBoth;
      *nelBoth = temp;

      /* Elide time from "second" list" */
      temp = (*nelSecond);
      *nelSecond = (*nelSecond)->next;
    }

    /* Advance list 1 */
    else if ((*nelFirst)->value < (*nelSecond)->value)
      nelFirst = &(*nelFirst)->next;

    /* Advance list 2 */
    else
      nelSecond = &(*nelSecond)->next;
  }

  if ((*first)->iteration_list == NULL)
  {
    (*first)->next = *discard;
    *discard = *first;
    *first = NULL;
  }

  if ((*second)->iteration_list == NULL)
  {
    (*second)->next = *discard;
    *discard = *second;
    *second = NULL;
  }

  if ((*both)->iteration_list == NULL)
  {
    (*both)->next = *discard;
    *discard = *both;
    *both = NULL;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_generic_discard_frames:
    Discards all frames in the list, freeing any memory associated with them.

        In:  struct frame_data_list *discard - frame discard pile
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_generic_discard_frames(struct frame_data_list *discard)
{
  struct frame_data_list *dnext;
  while (discard != NULL)
  {
    dnext = discard->next;
    free_num_expr_list(discard->iteration_list);
    free(discard);
    discard = dnext;
  }
  return 0;
}

/*************************************************************************
dreamm_v3_generic_new_frame:
    Allocates a new frame of a given type.  This is used when merging molecule
    frames into ALL_MOL_DATA frames or mesh frames into ALL_MESH_DATA frames.

        In:  int type - the frame type (ALL_MOL_DATA, etc.)
        Out: the frame, or NULL if allocation failed
**************************************************************************/
static struct frame_data_list *dreamm_v3_generic_new_frame(int type)
{
  struct frame_data_list *fdlp = (struct frame_data_list *) malloc(sizeof(struct frame_data_list));
  if (fdlp == NULL)
  {
    fprintf(world->err_file,"File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return NULL;
  }

  fdlp->next = NULL;
  fdlp->list_type = OUTPUT_BY_ITERATION_LIST;
  fdlp->type = type;
  fdlp->viz_iteration = 0;
  fdlp->n_viz_iterations = 0;
  fdlp->iteration_list = NULL;
  fdlp->curr_viz_iteration = NULL;
  return fdlp;
}

/*************************************************************************
dreamm_v3_generic_preprocess_frame_data:
    Preprocesses the frame data, normalizing it before the frames begin to be
    rendered.  This greatly simplifies the output logic.  The normalizations
    that are done are:

        1. Convert all frames to use iteration number rather than a mixture of
           iteration numbers and times

        2. Coalesce different frames of the same type into individual frames.

        3. Merge coincident MOL_POS and MOL_ORIENT frames into ALL_MOL_DATA
           frames, and coincident MESH_GEOMETRY and REG_DATA frames into
           ALL_MESH_DATA frames.

        4. Remove duplicated iteration numbers from all lists.

        5. Discard any frames which have no iterations.

        In:  struct frame_data_list *fdlpp - pointer to head of frame_data_list
        Out: 0 on success, 1 on failure.  *fdlpp will be updated if successful
**************************************************************************/
static int dreamm_v3_generic_preprocess_frame_data(struct frame_data_list **fdlpp)
{
  struct frame_data_list *mol_pos = NULL,
                         *mol_orient = NULL,
                         *all_mol_data = NULL,
                         *surf_pos = NULL,
                         *region_data = NULL,
                         *all_mesh_data = NULL,
                         *discard = NULL;
  struct frame_data_list *fdlp, *fdlpnext;

  /* Parse out all frames */
  for (fdlp = *fdlpp;
       fdlp != NULL;
       fdlp = fdlpnext)
  {
    fdlpnext = fdlp->next;
    fdlp->next = NULL;

    /* Convert frames to output by iteration list */
    if (fdlp->list_type == OUTPUT_BY_TIME_LIST)
      convert_frame_data_to_iterations(fdlp);

    /* Sort frames into one of our 6 lists */
    switch (fdlp->type)
    {
      case MOL_POS:
        if (mol_pos == NULL)
          mol_pos = fdlp;
        else
        {
          dreamm_v3_generic_merge_frame_data(mol_pos, fdlp);
          fdlp->next = discard;
          discard = fdlp;
        }
        break;

      case MOL_ORIENT:
        if (mol_orient == NULL)
          mol_orient = fdlp;
        else
        {
          dreamm_v3_generic_merge_frame_data(mol_orient, fdlp);
          fdlp->next = discard;
          discard = fdlp;
        }
        break;

      case ALL_MOL_DATA:
        if (all_mol_data == NULL)
          all_mol_data = fdlp;
        else
        {
          dreamm_v3_generic_merge_frame_data(all_mol_data, fdlp);
          fdlp->next = discard;
          discard = fdlp;
        }
        break;

      case MESH_GEOMETRY:
        if (surf_pos == NULL)
          surf_pos = fdlp;
        else
        {
          dreamm_v3_generic_merge_frame_data(surf_pos, fdlp);
          fdlp->next = discard;
          discard = fdlp;
        }
        break;

      case REG_DATA:
        if (region_data == NULL)
          region_data = fdlp;
        else
        {
          dreamm_v3_generic_merge_frame_data(region_data, fdlp);
          fdlp->next = discard;
          discard = fdlp;
        }
        break;

      case ALL_MESH_DATA:
        if (all_mesh_data == NULL)
          all_mesh_data = fdlp;
        else
        {
          dreamm_v3_generic_merge_frame_data(all_mesh_data, fdlp);
          fdlp->next = discard;
          discard = fdlp;
        }
        break;

      default:
        /* Unrecognized frame type.  Leave it in the original list */
        *fdlpp = fdlp;
        fdlpp = &fdlp->next;
    }
  }
  *fdlpp = NULL;

  /* Merge coincident frames  */
  if (mol_pos  &&  mol_orient)
  {
    if (all_mol_data == NULL)
      if ((all_mol_data = dreamm_v3_generic_new_frame(ALL_MOL_DATA)) == NULL)
        return 1;
    dreamm_v3_generic_merge_coincident_frames(&all_mol_data, &mol_pos, &mol_orient, &discard);
  }
  if (surf_pos  &&  region_data)
  {
    if (all_mesh_data == NULL)
      if ((all_mesh_data = dreamm_v3_generic_new_frame(ALL_MESH_DATA)) == NULL)
        return 1;
    dreamm_v3_generic_merge_coincident_frames(&all_mesh_data, &surf_pos, &region_data, &discard);
  }

  /* Remove duplicated iteration numbers */
  if (all_mol_data) uniq_num_expr_list(all_mol_data->iteration_list);
  if (mol_pos) uniq_num_expr_list(mol_pos->iteration_list);
  if (mol_orient) uniq_num_expr_list(mol_orient->iteration_list);
  if (all_mesh_data) uniq_num_expr_list(all_mesh_data->iteration_list);
  if (surf_pos) uniq_num_expr_list(surf_pos->iteration_list);
  if (region_data) uniq_num_expr_list(region_data->iteration_list);

  /* Discard now-empty frames */
  if (discard)
    dreamm_v3_generic_discard_frames(discard);

  /* Chain remaining frames together */
  if (all_mesh_data)
  {
    *fdlpp = all_mesh_data;
    fdlpp = &all_mesh_data->next;
  }
  if (surf_pos)
  {
    *fdlpp = surf_pos;
    fdlpp = &surf_pos->next;
  }
  if (region_data)
  {
    *fdlpp = region_data;
    fdlpp = &region_data->next;
  }
  if (all_mol_data)
  {
    *fdlpp = all_mol_data;
    fdlpp = &all_mol_data->next;
  }
  if (mol_pos)
  {
    *fdlpp = mol_pos;
    fdlpp = &mol_pos->next;
  }
  if (mol_orient)
  {
    *fdlpp = mol_orient;
    fdlpp = &mol_orient->next;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_generic_scan_for_frame
    Scan a list of frame data objects to see if there is an upcoming frame for
    the specified iteration number.

        In: struct frame_data_list *fdlp - the list to scan
            long long iterno - the iteration number
        Out: 1 if a matching frame is found, 0 otherwise.
**************************************************************************/
static int dreamm_v3_generic_scan_for_frame(struct frame_data_list *fdlp, long long iterno)
{
  for (; fdlp != NULL; fdlp = fdlp->next)
    if (fdlp->viz_iteration == iterno)
      return 1;
  return 0;
}

/*************************************************************************
dreamm_v3_generic_dump_time_values:
    Writes the time values to the time values data file.

        In:  char const *viz_data_dir - directory to receive time data
             char const *time_values_name - name of time data file
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_dump_time_values(char const *viz_data_dir,
                                              char const *time_values_name)
{
  FILE *time_values_data = NULL;
  int time_value_index;

  /* Open time values data file */
  if ((time_values_data = dreamm_v3_generic_open_file(viz_data_dir, time_values_name, "wb")) == NULL)
    return 1;

  /* Write out time values */
  for (time_value_index = 0;
       time_value_index < world->viz_state_info.output_times.n_iterations;
       ++ time_value_index)
  {
    double t_value = world->viz_state_info.output_times.iterations[time_value_index] * world->time_unit;
    fwrite(&t_value, sizeof(t_value), 1, time_values_data);
  }

  fclose(time_values_data);
  return 0;
}

/*************************************************************************
dreamm_v3_generic_dump_iteration_numbers:
    Writes the iteration numbers to the iteration numbers data file.

        In:  char const *viz_data_dir - directory to receive time data
             char const *iteration_numbers_name - name of iteration data
             u_int iteration_numbers_count - maximum number of iterations in
                                   any of the data

        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_dump_iteration_numbers(char const *viz_data_dir,
                                                    char const *iteration_numbers_name,
                                                    u_int iteration_numbers_count)
{
  FILE *iteration_numbers_data = NULL;
  long long iteration_index;

  /* Open iteration numbers data file */
  if ((iteration_numbers_data = dreamm_v3_generic_open_file(viz_data_dir, iteration_numbers_name, "wb")) == NULL)
    return 1;

  /* Write out iteration data */
  int last_mesh = -1, last_vol_mol = -1, last_surf_mol = -1;
  for (iteration_index = 0; iteration_index < iteration_numbers_count; iteration_index++)
  {
    /* Write meshes iteration */
    if (iteration_index < world->viz_state_info.mesh_output_iterations.n_iterations)
      last_mesh = world->viz_state_info.mesh_output_iterations.iterations[iteration_index];
    fwrite(&last_mesh, sizeof(last_mesh), 1, iteration_numbers_data);

    /* Write vol mols iteration */
    if (iteration_index < world->viz_state_info.vol_mol_output_iterations.n_iterations)
      last_vol_mol = world->viz_state_info.vol_mol_output_iterations.iterations[iteration_index];
    fwrite(&last_vol_mol, sizeof(last_vol_mol), 1, iteration_numbers_data);

    /* Write surface mols iteration */
    if (iteration_index < world->viz_state_info.grid_mol_output_iterations.n_iterations)
      last_surf_mol = world->viz_state_info.grid_mol_output_iterations.iterations[iteration_index];
    fwrite(&last_surf_mol, sizeof(last_surf_mol), 1, iteration_numbers_data);
  }

  fclose(iteration_numbers_data);
  return 0;
}

/*************************************************************************
dreamm_v3_generic_write_time_info:
    Writes the timing info to the master header file.

        In:  FILE *master_header - the master header to receive index info
             char const *iteration_numbers_name - name of iteration data
             char const *time_values_name - name of time data
             char const *dreamm3mode - the dreamm3mode as a string
             int dreamm3mode_number - the dreamm3mode as an integer
             u_int iteration_numbers_count - number of iteration data
             u_int time_values_count - number of time data
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_write_time_info(FILE *master_header,
                                             char const *iteration_numbers_name,
                                             char const *time_values_name,
                                             char const *dreamm3mode,
                                             int dreamm3mode_number,
                                             u_int iteration_numbers_count,
                                             u_int time_values_count)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";

  /* Write iteration object to header */
  fprintf(master_header,
          "object \"iteration_numbers\" class array "
          "type unsigned int "
          "rank 1 shape 3 items %u "
          "%s binary data "
          "file %s,0\n",
          iteration_numbers_count,
          my_byte_order,
          iteration_numbers_name);
  fprintf(master_header,
          "\tattribute \"dreamm3mode\" number %d\t#%s#\n",
          dreamm3mode_number,
          dreamm3mode);
  fprintf(master_header, "\n\n");

  /* If we have time data, write time object to header file */
  if (time_values_count > 0)
  {
    fprintf(master_header,
            "object \"time_values\" class array "
            "type double "
            "rank 0 items %u "
            "%s binary data "
            "file %s,0\n",
            time_values_count,
            my_byte_order,
            time_values_name);
    fprintf(master_header,
            "\tattribute \"dreamm3mode\" number %d\t#%s#\n",
            dreamm3mode_number,
            dreamm3mode);
    fprintf(master_header, "\n\n");
  }

  return 0;
}

/*************************************************************************
dreamm_v3_generic_write_mesh_fields:
    Writes the mesh fields to the header file.

        In:  FILE *meshes_header - the header to receive index info
             int field_index_base - base index for field objects
             int surf_index - base index for mesh data object numbers
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_write_mesh_fields(struct frame_data_list const * const fdlp,
                                               FILE *meshes_header,
                                               int field_index_base,
                                               int surf_index)
{
  byte viz_surf_pos_flag    = (fdlp->type == ALL_MESH_DATA || fdlp->type == MESH_GEOMETRY);
  byte viz_surf_states_flag = (viz_surf_pos_flag  &&  (world->viz_output_flag & VIZ_SURFACE_STATES) != 0);
  byte viz_region_data_flag = (fdlp->type == ALL_MESH_DATA || fdlp->type == REG_DATA);

  int obj_index;
  for (obj_index = 0; obj_index < world->viz_state_info.n_viz_objects; obj_index++)
  {
    fprintf(meshes_header,
            "object %d field # %s #\n",
            field_index_base ++,
            world->viz_state_info.viz_objects[obj_index]->sym->name);
    if (viz_surf_pos_flag)
      fprintf(meshes_header,
              "\tcomponent \"positions\" value %d\n",
              surf_index ++);
    if (viz_surf_pos_flag)
      fprintf(meshes_header,
              "\tcomponent \"connections\" value %d\n",
              surf_index ++);
    if (viz_surf_states_flag)
      fprintf(meshes_header,
              "\tcomponent \"state_values\" value %d\n",
              surf_index ++);
    if (viz_region_data_flag)
    {
      struct region_list *rlp;
      for (rlp = world->viz_state_info.viz_objects[obj_index]->regions;
           rlp != NULL;
           rlp = rlp->next)
      {
        if (! strcmp(rlp->reg->region_last_name, "ALL"))
          continue;

        fprintf(meshes_header,
                "\tcomponent \"%s\" value %d\n",
                rlp->reg->region_last_name,
                surf_index ++);
      }
    }
    fprintf(meshes_header, "\n");
  }

  return 0;
}

/*************************************************************************
dreamm_v3_generic_write_rank0_int_array_index:
    Writes index info for a rank 0 integer array to a header file.

        In:  FILE *header - the header to receive index info
             int obj_index - the index number for the object
             int array_length - number of items in array
             char const *filename - filename with data
             long file_offset - offset within file for data
             char const *symname - name for symbol
             char const *objtype - type of object (for comment)
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_write_rank0_int_array_index(FILE *header,
                                                         int obj_index,
                                                         int array_length,
                                                         char const *filename,
                                                         long file_offset,
                                                         char const *symname,
                                                         char const *objtype)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";

  fprintf(header,
          "object %d class array type int "
          "rank 0 items %d "
          "%s binary data "
          "file %s,%ld # %s.%s #\n",
          obj_index,
          array_length,
          my_byte_order,
          filename,
          file_offset,
          symname,
          objtype);
  return 0;
}

/*************************************************************************
dreamm_v3_generic_write_rank1_int_array_index:
    Writes index info for a rank 1 integer array to a header file.

        In:  FILE *header - the header to receive index info
             int obj_index - the index number for the object
             int array_length - number of items in array
             char const *filename - filename with data
             long file_offset - offset within file for data
             char const *symname - name for symbol
             char const *objtype - type of object (for comment)
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_write_rank1_int_array_index(FILE *header,
                                                         int obj_index,
                                                         int array_length,
                                                         char const *filename,
                                                         long file_offset,
                                                         char const *symname,
                                                         char const *objtype)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";

  fprintf(header,
          "object %d class array type int "
          "rank 1 shape 3 items %d "
          "%s binary data "
          "file %s,%ld # %s.%s #\n",
          obj_index,
          array_length,
          my_byte_order,
          filename,
          file_offset,
          symname,
          objtype);
  return 0;
}

/*************************************************************************
dreamm_v3_generic_write_float_array_index:
    Writes index info for a float array to a header file.

        In:  FILE *header - the header to receive index info
             int obj_index - the index number for the object
             int array_length - number of items in array
             char const *filename - filename with data
             long file_offset - offset within file for data
             char const *symname - name for symbol
             char const *objtype - type of object (for comment)
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_write_float_array_index(FILE *header,
                                                     int obj_index,
                                                     int array_length,
                                                     char const *filename,
                                                     long file_offset,
                                                     char const *symname,
                                                     char const *objtype)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";

  if (array_length <= 0)
    fprintf(header, "object %d array # %s.%s #\n",
            obj_index,
            symname,
            objtype);
  else
    fprintf(header,
            "object %d class array type float "
            "rank 1 shape 3 items %d "
            "%s binary data "
            "file %s,%ld # %s.%s #\n",
            obj_index,
            array_length,
            my_byte_order,
            filename,
            file_offset,
            symname,
            objtype);
  return 0;
}

/*************************************************************************
dreamm_v3_generic_dump_mesh_data:
    Writes the mesh data to mesh data files, and appropriate index info to the
    header file.  Mesh object indices are assigned, in order, to position,
    connections, states, and region info, omitting whichever indices are not
    needed.

        In:  struct frame_data_list const * const fdlp - the frame to write
             FILE *meshes_header - the header to receive index info
             char const *dirname - the directory to receive data
             char const *mesh_pos_filename - filename for mesh pos data
             char const *mesh_states_filename - filename for mesh state data
             char const *region_data_filename - filename for region data
             int *main_index_base - ptr to index for allocating obj numbers
             int *surf_index_base - ptr to rcv index for mesh data objects
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_dump_mesh_data(struct frame_data_list const * const fdlp,
                                            FILE *meshes_header,
                                            char const *dirname,
                                            char const *mesh_pos_filename,
                                            char const *mesh_states_filename,
                                            char const *region_data_filename,
                                            int *meshes_main_index)
{
  /* File handles */
  FILE *mesh_pos_data = NULL;
  FILE *mesh_states_data = NULL;
  FILE *region_data = NULL;

  /* Index for iteration over all mesh objects */
  int obj_index;

  /* Control flags */
  byte viz_surf_pos_flag    = (fdlp->type == ALL_MESH_DATA || fdlp->type == MESH_GEOMETRY);
  byte viz_surf_states_flag = (viz_surf_pos_flag  &&  (world->viz_output_flag & VIZ_SURFACE_STATES) != 0);
  byte viz_region_data_flag = (fdlp->type == ALL_MESH_DATA || fdlp->type == REG_DATA);

  /* Open Surface Positions file, if necessary */
  if (viz_surf_pos_flag  &&  (mesh_pos_data = dreamm_v3_generic_open_file(dirname, mesh_pos_filename, "ab")) == NULL)
    goto failure;

  /* Open Surface States file, if necessary */
  if (viz_surf_states_flag  &&  (mesh_states_data = dreamm_v3_generic_open_file(dirname, mesh_states_filename, "ab")) == NULL)
    goto failure;

  /* Open Region Data file, if necessary */
  if (viz_region_data_flag  &&  (region_data = dreamm_v3_generic_open_file(dirname, region_data_filename, "ab")) == NULL)
    goto failure;

  /* Traverse all visualized objects and output mesh/region data */
  for (obj_index = 0; obj_index < world->viz_state_info.n_viz_objects; ++ obj_index)
  {
    struct object *objp = world->viz_state_info.viz_objects[obj_index];
    if (objp->viz_state == NULL) continue;

    if (objp->object_type != POLY_OBJ  &&  objp->object_type != BOX_OBJ)
      continue;

    struct polygon_object *pop = (struct polygon_object *) objp->contents;
    struct ordered_poly *opp = (struct ordered_poly *) pop->polygon_data;
    struct element_data *edp = opp->element;
    int element_data_count = objp->n_walls_actual;

    if (viz_surf_pos_flag)
    {
      int wall_index;
      dreamm_v3_generic_write_float_array_index(meshes_header,
                                                (*meshes_main_index) ++,
                                                objp->n_verts,
                                                mesh_pos_filename,
                                                ftell(mesh_pos_data),
                                                objp->sym->name,
                                                "positions");
      fprintf(meshes_header,
              "\tattribute \"dep\" string \"positions\"\n\n");

      /* output polyhedron vertices */
      dx_output_vertices(mesh_pos_data, objp);

      /* output polygon element connections */
      dreamm_v3_generic_write_rank1_int_array_index(meshes_header,
                                                    (*meshes_main_index) ++,
                                                    element_data_count,
                                                    mesh_pos_filename,
                                                    ftell(mesh_pos_data),
                                                    objp->sym->name,
                                                    "connections");
      fprintf(meshes_header,
              "\tattribute \"ref\" string \"positions\"\n");
      fprintf(meshes_header,
              "\tattribute \"element type\" string \"triangles\"\n\n");

      for (wall_index = 0; wall_index < objp->n_walls; ++ wall_index)
        if (! get_bit(pop->side_removed, wall_index))
          dx_output_wall_vertices(mesh_pos_data, &edp[ wall_index]);
    }

    if (viz_surf_states_flag)
    {
      int wall_index;
      dreamm_v3_generic_write_rank0_int_array_index(meshes_header,
                                                    (*meshes_main_index) ++,
                                                    element_data_count,
                                                    mesh_states_filename,
                                                    ftell(mesh_states_data),
                                                    objp->sym->name,
                                                    "states");
      fprintf(meshes_header,
              "\tattribute \"dep\" string \"connections\"\n\n");

      /* XXX: Why write this as binary? */
      for (wall_index = 0; wall_index < objp->n_walls; ++ wall_index)
      {
        if (! get_bit(pop->side_removed, wall_index))
        {
          int state=objp->viz_state[wall_index];
          fwrite(&state, sizeof (state), 1, mesh_states_data);
        }
      }
    }

    if (viz_region_data_flag && (objp->num_regions > 1))
    {
      struct region_list *rlp;

      for(rlp = objp->regions; rlp != NULL; rlp = rlp->next)
      {
        int wall_index;
        struct region *rp = rlp->reg;
        if (strcmp(rp->region_last_name, "ALL") == 0) continue; 

        /* number of walls in the region */
        int region_walls_number = 0; 
        long pos = ftell(region_data);

        for(wall_index = 0; wall_index < objp->n_walls; ++ wall_index)
        {
          int n = objp->wall_p[wall_index]->side;
          if (get_bit(rp->membership,n))
          {
            fwrite(&n, sizeof (n), 1, region_data);
            region_walls_number++;
          }
        }

        dreamm_v3_generic_write_rank0_int_array_index(meshes_header,
                                                      (*meshes_main_index) ++,
                                                      region_walls_number,
                                                      region_data_filename,
                                                      pos,
                                                      objp->sym->name,
                                                      "region_data");
        fprintf(meshes_header,
                "\tattribute \"ref\" string \"connections\"\n");
        fprintf(meshes_header,
                "\tattribute \"identity\" string \"region_indices\"\n");
        fprintf(meshes_header,
                "\tattribute \"name\" string \"%s\"\n",
                rp->region_last_name);
        if (rp->region_viz_value > 0)
          fprintf(meshes_header,
                  "\tattribute \"viz_value\" number %d\n",
                  rp->region_viz_value);

        fprintf(meshes_header, "\n\n");
      } /* end for */
    } /* end if (region_data_flag) */
  }

  if (mesh_pos_data) fclose(mesh_pos_data);
  if (mesh_states_data) fclose(mesh_states_data);
  if (region_data) fclose(region_data);
  return 0;

failure:
  if (mesh_pos_data) fclose(mesh_pos_data);
  if (mesh_states_data) fclose(mesh_states_data);
  if (region_data) fclose(region_data);
  return 1;
}

/*************************************************************************
dreamm_v3_generic_write_molecule_fields:
    Writes the molecule fields to the header file.

        In:  FILE *mol_header - the header to receive index info
             struct species **specs - all relevant species
             int num_molecules - num relevant species
             int field_index - base index for field objects
             int mol_data_index - base index for mol data objects
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_write_molecule_fields(struct frame_data_list const * const fdlp,
                                                   FILE *mol_header,
                                                   struct species **specs,
                                                   int num_molecules,
                                                   int field_idx,
                                                   int mol_data_index)
{
  /* Control flags */
  byte viz_mol_pos_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_POS);
  byte viz_mol_orient_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_ORIENT);
  byte viz_mol_states_flag = (viz_mol_pos_flag  &&  (world->viz_output_flag & VIZ_MOLECULES_STATES));

  int mol_index;
  /* Build fields for molecules here */
  for (mol_index=0; mol_index < num_molecules; ++ mol_index)
  {
    fprintf(mol_header,
            "object %d field # %s #\n",
            field_idx ++,
            specs[mol_index]->sym->name);
    if (viz_mol_pos_flag)
      fprintf(mol_header,
              "\tcomponent \"positions\" value %d\n",
              mol_data_index ++);
    if (viz_mol_orient_flag)
      fprintf(mol_header,
              "\tcomponent \"data\" value %d # orientations #\n",
              mol_data_index ++);
    if (viz_mol_states_flag)
      fprintf(mol_header,
              "\tcomponent \"state_values\" value %d\n",
              mol_data_index ++);
    fprintf(mol_header, "\n");
  }
  return 0;
}

/*************************************************************************
dreamm_v3_generic_write_vol_orientations_index:
    Write the orientations constant array for volume molecules to the index.

        In:  FILE *mol_header - the header to receive index info
             int obj_index - object index number
             int count - number of molecules
             char const *filename - filename for data
             long file_offset - offset within file for data
             char const *symname - symbol name for molecule
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_write_vol_orientations_index(FILE *mol_header,
                                                          int obj_index,
                                                          int count,
                                                          char const *filename,
                                                          long file_offset,
                                                          char const *symname)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";

  if (count > 0)
  {
    fprintf(mol_header,
            "object %d class constantarray type float "
            "rank 1 shape 3 items %d "
            "%s binary data "
            "file %s,%ld # %s.orientations #\n",
            obj_index,
            count,
            my_byte_order,
            filename,
            file_offset,
            symname);
    fprintf(mol_header,
            "\tattribute \"dep\" string \"positions\"\n\n");
  }
  else
    fprintf(mol_header,
            "object %d array # %s.orientations #\n",
            obj_index,
            symname);
  return 0;
}

/*************************************************************************
dreamm_v3_generic_write_state_array_index:
    Write the state index constant array to the header

        In:  FILE *mol_header - the header to receive index info
             int obj_index - object index number
             int count - number of molecules
             char const *filename - filename for data
             long file_offset - offset within file for data
             char const *symname - symbol name for molecule
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_write_state_array_index(FILE *mol_header,
                                                     int obj_index,
                                                     int count,
                                                     char const *filename,
                                                     long file_offset,
                                                     char const *symname)
{
  char const *my_byte_order = (*(unsigned char *)&endian_test_word == 1) ? "lsb" : "msb";

  if (count > 0)
  {
    fprintf(mol_header,
            "object %d class constantarray type int "
            "items %d "
            "%s binary data "
            "file %s,%ld # %s.states #\n",
            obj_index,
            count,
            my_byte_order,
            filename,
            file_offset,
            symname);
    fprintf(mol_header,
            "\tattribute \"dep\" string \"positions\"\n\n");
  }
  else
    fprintf(mol_header,
            "object %d array # %s.states #\n",
            obj_index,
            symname);
  return 0;
}

/*************************************************************************
dreamm_v3_generic_dump_grid_molecule_data:
    Writes the grid molecule data to appropriate data files and index info to
    the header file.  Object numbers are assigned first to position, then to
    orientation, and finally to state data objects.

        In:  struct frame_data_list const * const fdlp - frame to write
             FILE *surf_mol_header - the header to receive index info
             char const *dirname - directory for data files
             char const *mol_pos_name - name for position data file
             char const *mol_orient_name - name for orientation data file
             char const *mol_states_name - name for state data file
             int *main_index - ptr to index for allocating obj numbers
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_dump_grid_molecule_data(struct frame_data_list const * const fdlp,
                                                     FILE *surf_mol_header,
                                                     char const *dirname,
                                                     char const *mol_pos_name,
                                                     char const *mol_orient_name,
                                                     char const *mol_states_name,
                                                     int *main_index)
{
  /* File handles */
  FILE *surf_mol_pos_data = NULL;
  FILE *surf_mol_states_data = NULL;
  FILE *surf_mol_orient_data = NULL;

  /* Grid molecules, sorted into species */
  struct grid_molecule ***grid_mols_by_species = NULL;
  u_int *grid_mol_counts_by_species = NULL;

  /* Iteration variables */
  int species_index;

  /* Control flags */
  byte viz_mol_pos_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_POS);
  byte viz_mol_orient_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_ORIENT);
  byte viz_mol_states_flag = (viz_mol_pos_flag  &&  (world->viz_output_flag & VIZ_MOLECULES_STATES));

  /* Open surface molecules position data file */
  if (viz_mol_pos_flag  &&  (surf_mol_pos_data = dreamm_v3_generic_open_file(dirname, mol_pos_name, "ab")) == NULL)
    goto failure;

  /* Open surface molecules orientation data file */
  if (viz_mol_orient_flag  &&  (surf_mol_orient_data =  dreamm_v3_generic_open_file(dirname, mol_orient_name, "ab")) == NULL)
    goto failure;

  /* Open surface molecules states data file */
  if (viz_mol_states_flag  &&  (surf_mol_orient_data =  dreamm_v3_generic_open_file(dirname, mol_states_name, "ab")) == NULL)
    goto failure;

  /* Get a list of molecules sorted by species. */
  if (sort_molecules_by_species((struct abstract_molecule ****) (void *) &grid_mols_by_species,
                                &grid_mol_counts_by_species,
                                0, 1))
    goto failure;

  /* Emit all molecules for each species */
  for (species_index = 0; species_index < world->n_species; ++ species_index)
  {
    long fpos_pos = 0;
    long fpos_orient = 0;
    long fpos_states = 0;
    struct species *specp = world->species_list[species_index];

    /* Skip objects marked for exclusion and non-grid molecules */
    if (specp->viz_state == EXCLUDE_OBJ  ||  ! (specp->flags & ON_GRID))
      continue;

    /* Save offsets of data for current species */
    if (grid_mol_counts_by_species[species_index] > 0)
    {
      if (viz_mol_pos_flag) fpos_pos = ftell(surf_mol_pos_data);
      if (viz_mol_orient_flag) fpos_orient = ftell(surf_mol_orient_data);
      if (viz_mol_states_flag) fpos_states = ftell(surf_mol_states_data);
    }

    /* Iterate over specific molecules in this species */
    int count = 0;
    unsigned int mol_index;
    for (mol_index = 0; mol_index < grid_mol_counts_by_species[species_index]; ++ mol_index)
    {
      struct grid_molecule *gmol = (grid_mols_by_species[species_index])[mol_index];
      struct wall *w = gmol->grid->surface;
      int objidx = void_array_search((void **) world->viz_state_info.viz_objects,
                                     world->viz_state_info.n_viz_objects,
                                     w->parent_object);

      /* If the molecule is on an object we don't care about, skip it */
      if (objidx == -1)
        continue;

      /* Keep count of the items we write */
      ++ count;

      /* Write positions information */
      if (viz_mol_pos_flag)
      {
        struct vector3 p0;
        grid2xyz(gmol->grid, gmol->grid_index, &p0);
        dx_output_vector3(surf_mol_pos_data, &p0);
      }

      /* Write orientations information */
      if (viz_mol_orient_flag)
        dx_output_oriented_normal(surf_mol_orient_data, &w->normal, gmol->orient);
    }

    /* Write data for surface molecule states */
    if (viz_mol_states_flag  &&  count > 0)
    {
      int state = specp->viz_state;
      fwrite(&state, sizeof state, 1, surf_mol_states_data);
    }

    /* Emit array of mol positions */
    if (viz_mol_pos_flag)
    {
      dreamm_v3_generic_write_float_array_index(surf_mol_header,
                                                (*main_index) ++,
                                                count,
                                                mol_pos_name,
                                                fpos_pos,
                                                specp->sym->name,
                                                "positions");
      if (count > 0) fprintf(surf_mol_header, "\tattribute \"dep\" string \"positions\"\n\n");
    }

    /* Emit array of mol orientations */
    if (viz_mol_orient_flag)
    {
      dreamm_v3_generic_write_float_array_index(surf_mol_header,
                                                (*main_index) ++,
                                                count,
                                                mol_orient_name,
                                                fpos_orient,
                                                specp->sym->name,
                                                "orientations");
      if (count > 0) fprintf(surf_mol_header, "\tattribute \"dep\" string \"positions\"\n\n");
    }

    /* Emit array of mol states */
    if (viz_mol_states_flag)
    {
      dreamm_v3_generic_write_state_array_index(surf_mol_header,
                                                (*main_index) ++,
                                                count,
                                                mol_states_name,
                                                fpos_states,
                                                specp->sym->name);
    }

    /* For readability, add an extra newline */
    if (count <= 0)
      fprintf(surf_mol_header, "\n");
  }

  free_ptr_array((void **) grid_mols_by_species, world->n_species);
  free(grid_mol_counts_by_species);
  if (surf_mol_pos_data) fclose(surf_mol_pos_data);
  if (surf_mol_orient_data) fclose(surf_mol_orient_data);
  if (surf_mol_states_data) fclose(surf_mol_states_data);
  return 0;

failure:
  if (grid_mols_by_species) free_ptr_array((void **) grid_mols_by_species, world->n_species);
  if (grid_mol_counts_by_species) free(grid_mol_counts_by_species);
  if (surf_mol_pos_data) fclose(surf_mol_pos_data);
  if (surf_mol_orient_data) fclose(surf_mol_orient_data);
  if (surf_mol_states_data) fclose(surf_mol_states_data);
  return 1;
}

/*************************************************************************
dreamm_v3_generic_dump_volume_molecule_data:
    Writes the volume molecule data to appropriate data files and index info to
    the header file.  Object id numbers are allocated, for each molecule, first
    to position, then to orientation, then to state, omitting any which are not
    desired in this output frame.

        In:  struct frame_data_list const * const fdlp - frame to write
             FILE *vol_mol_header - the header to receive index info
             char const *dirname - directory for data files
             char const *mol_pos_name - name for position data file
             char const *mol_orient_name - name for orientation data file
             char const *mol_states_name - name for state data file
             int *main_index - ptr to index for allocating obj numbers
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_generic_dump_volume_molecule_data(struct frame_data_list const * const fdlp,
                                                       FILE *vol_mol_header,
                                                       char const *dirname,
                                                       char const *mol_pos_name,
                                                       char const *mol_orient_name,
                                                       char const *mol_states_name,
                                                       int *main_index)
{
  /* File handles */
  FILE *vol_mol_pos_data = NULL;
  FILE *vol_mol_states_data = NULL;
  FILE *vol_mol_orient_data = NULL;

  /* All volume molecules, sorted into species */
  struct volume_molecule ***viz_molp = NULL;
  u_int *viz_mol_count = NULL;

  /* Iteration variables */
  int species_index;

  /* Control flags */
  byte viz_mol_pos_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_POS);
  byte viz_mol_orient_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_ORIENT);
  byte viz_mol_states_flag = (viz_mol_pos_flag  &&  (world->viz_output_flag & VIZ_MOLECULES_STATES));

  /* Prepare for position output */
  if (viz_mol_pos_flag  &&  (vol_mol_pos_data = dreamm_v3_generic_open_file(dirname, mol_pos_name, "ab")) == NULL)
    goto failure;

  /* Prepare for orientation output */
  if (viz_mol_orient_flag  &&  (vol_mol_orient_data = dreamm_v3_generic_open_file(dirname, mol_orient_name, "ab")) == NULL)
    goto failure;

  /* Prepare for states output */
  if (viz_mol_states_flag  &&  (vol_mol_states_data = dreamm_v3_generic_open_file(dirname, mol_states_name, "ab")) == NULL)
    goto failure;

  /* Get a list of molecules sorted by species. */
  if (sort_molecules_by_species((struct abstract_molecule ****) (void *) &viz_molp,
                                &viz_mol_count,
                                1, 0))
    goto failure;

  /* Process all volume mols */
  for (species_index=0; species_index<world->n_species; ++ species_index)
  {
    struct species *specp = world->species_list[species_index];

    /* Skip this species if it is marked as excluded, or not a vol mol */
    if (specp->viz_state == EXCLUDE_OBJ  ||  (specp->flags & NOT_FREE))
      continue;

    /* Check that our molecule count agrees with the species population count */
    if (viz_mol_count[species_index] != specp->population)
    {
      fprintf(world->log_file, "MCell: molecule count disagreement!!\n");
      fprintf(world->log_file, "  Species %s  population = %d  count = %d\n",
              specp->sym->name,
              specp->population,
              viz_mol_count[species_index]);
    }

    /* Emit an array of molecule positions */
    if (viz_mol_pos_flag)
    {
      unsigned int mol_index;
      dreamm_v3_generic_write_float_array_index(vol_mol_header,
                                                (*main_index) ++,
                                                viz_mol_count[species_index],
                                                mol_pos_name,
                                                ftell(vol_mol_pos_data),
                                                specp->sym->name,
                                                "positions");
      if (viz_mol_count[species_index] > 0)
      {
        fprintf(vol_mol_header,
                "\tattribute \"dep\" string \"positions\"\n\n");
        for (mol_index = 0; mol_index < viz_mol_count[species_index]; ++ mol_index)
          dx_output_vector3(vol_mol_pos_data, &(viz_molp[species_index][mol_index])->pos);
      }
    }

    /* Emit a constantarray of molecule orientations (constant [0, 0, 1]) */
    if (viz_mol_orient_flag)
    {
      dreamm_v3_generic_write_vol_orientations_index(vol_mol_header,
                                                     (*main_index) ++,
                                                     viz_mol_count[species_index],
                                                     mol_orient_name,
                                                     ftell(vol_mol_orient_data),
                                                     specp->sym->name);
      if (viz_mol_count[species_index] > 0)
        dx_output_oriented_normal(vol_mol_orient_data, &v3_unit_z, 1);
    }

    /* Emit molecule state */
    if (viz_mol_states_flag)
    {
        /* write molecule states information. */ 
        dreamm_v3_generic_write_state_array_index(vol_mol_header,
                                                  (*main_index) ++,
                                                  viz_mol_count[species_index],
                                                  mol_states_name,
                                                  ftell(vol_mol_states_data),
                                                  specp->sym->name);
      if (viz_mol_count[species_index] > 0)
      {
        int state = specp->viz_state;
        fwrite(&state,sizeof state,1,vol_mol_states_data);
      }
    }

    /* For readability, add an extra newline */
    if (viz_mol_count[species_index] <= 0)
      fprintf(vol_mol_header, "\n");
  }

  free_ptr_array((void **) viz_molp, world->n_species);
  free(viz_mol_count);
  if (vol_mol_pos_data) fclose(vol_mol_pos_data);
  if (vol_mol_orient_data) fclose(vol_mol_orient_data);
  if (vol_mol_states_data) fclose(vol_mol_states_data);
  return 0;

failure:
  if (viz_molp != NULL) free_ptr_array((void **) viz_molp, world->n_species);
  if (viz_mol_count != NULL) free(viz_mol_count);
  if (vol_mol_pos_data) fclose(vol_mol_pos_data);
  if (vol_mol_orient_data) fclose(vol_mol_orient_data);
  if (vol_mol_states_data) fclose(vol_mol_states_data);
  return 1;
}

/* == DREAMM V3 Output (non-grouped) == */

/* Filenames for DREAMM V3 (non-grouped) output */
static char const * const DREAMM_MESH_POS_NAME = "mesh_positions.bin";
static char const * const DREAMM_MESH_STATES_NAME = "mesh_states.bin";
static char const * const DREAMM_REGION_VIZ_DATA_NAME = "region_indices.bin";
static char const * const DREAMM_MESHES_HEADER_NAME = "meshes.dx";
static char const * const DREAMM_VOL_MOL_POS_NAME = "volume_molecules_positions.bin";
static char const * const DREAMM_SURF_MOL_POS_NAME = "surface_molecules_positions.bin";
static char const * const DREAMM_VOL_MOL_ORIENT_NAME = "volume_molecules_orientations.bin";
static char const * const DREAMM_SURF_MOL_ORIENT_NAME = "surface_molecules_orientations.bin";
static char const * const DREAMM_VOL_MOL_STATES_NAME = "volume_molecules_states.bin";
static char const * const DREAMM_SURF_MOL_STATES_NAME = "surface_molecules_states.bin";
static char const * const DREAMM_VOL_MOL_HEADER_NAME = "volume_molecules.dx";
static char const * const DREAMM_SURF_MOL_HEADER_NAME = "surface_molecules.dx";


/*************************************************************************
dreamm_v3_init:
    Initialize state for DREAMM_V3 output.

        In:  struct frame_data_list *fdlp - the head of the frame data list
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_init(struct frame_data_list *fdlp)
{
  if (dreamm_v3_generic_init())
    return 1;

  /* create viz_data dir filename */  
  char *viz_data_dir = my_strcat(world->file_prefix_name, "_viz_data");
  if (viz_data_dir == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  /* make viz_data dir */  
  if (mkdirs(viz_data_dir, world->err_file))
  {
    free(viz_data_dir);
    return 1;
  }

  /* create frame data dir filename */  
  world->viz_state_info.frame_data_dir = my_strcat(viz_data_dir, "/frame_data");
  free(viz_data_dir);
  if (world->viz_state_info.frame_data_dir == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  /* make directory for "frame data" */
  if (mkdirs(world->viz_state_info.frame_data_dir, world->err_file))
    return 1;

  /* Prepare iteration counters and timing info for frame_data_list */
  int time_values_total = count_time_values(fdlp);
  if (reset_time_values(fdlp, world->start_time))
    return 1;
  return initialize_iteration_counters(time_values_total);
}

/*************************************************************************
dreamm_v3_remove_file:
    Remove a file in the DREAMM V3 output directory.  Used by DREAMM V3 output
    to get rid of files at the beginning of each iteration.  We eliminate files
    which we are going to write to because some of these files may be symlinks.

        In:  char const *fname - filename for file
        Out: 0 on success, 1 on failure
**************************************************************************/
static int dreamm_v3_remove_file(char const *fname)
{
  struct stat f_stat;

  /* concatenate dir and fname to get the mesh states file path */
  char *path = alloc_sprintf("%s/%s", world->viz_state_info.iteration_number_dir, fname);
  if (path == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    goto failure;
  }

  /* If the file exists, remove it */
  if (stat(path, &f_stat) == 0)
  {
    if (unlink(path) != 0)
    {
      fprintf(world->err_file, "File %s, Line %ld: error %d removing file %s.\n", __FILE__, (long)__LINE__, errno, fname);
      goto failure;
    }
  }

  free(path);
  return 0;

failure:
  if (path) free(path);
  return 1;
}

/*************************************************************************
dreamm_v3_clean_files:
    Clean up any files which are going to be created during this iteration.  We
    need to remove these files because an old run may have created symlinks
    here, which will cause the files to overwrite data in older iterations.

        In: struct frame_data_list *fdlp - the frame data list
        Out: 0 on success, 1 on failure
**************************************************************************/
static int dreamm_v3_clean_files(struct frame_data_list *fdlp)
{
  if (! active_this_iteration(fdlp))
    return 0;

  /* Free old directory */
  if (world->viz_state_info.iteration_number_dir)
    free(world->viz_state_info.iteration_number_dir);
  world->viz_state_info.iteration_number_dir = NULL;

  /* Make new directory name */
  world->viz_state_info.iteration_number_dir = alloc_sprintf("%s/iteration_%lld",
                                                             world->viz_state_info.frame_data_dir,
                                                             world->it_time);
  if (world->viz_state_info.iteration_number_dir == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  /* If new directory doesn't exist, create it and return */
  if (! is_dir(world->viz_state_info.iteration_number_dir))
    return mkdirs(world->viz_state_info.iteration_number_dir, world->err_file);

  /* Directory already existed.  Clear out any files we're going to write */
  for (; fdlp != NULL; fdlp = fdlp->next)
  {
    if (fdlp->viz_iteration != world->it_time)
      continue;

    switch (fdlp->type)
    {
      case ALL_MOL_DATA:
        dreamm_v3_remove_file(DREAMM_VOL_MOL_POS_NAME);
        dreamm_v3_remove_file(DREAMM_VOL_MOL_ORIENT_NAME);
        dreamm_v3_remove_file(DREAMM_VOL_MOL_STATES_NAME);
        dreamm_v3_remove_file(DREAMM_SURF_MOL_POS_NAME);
        dreamm_v3_remove_file(DREAMM_SURF_MOL_ORIENT_NAME);
        dreamm_v3_remove_file(DREAMM_SURF_MOL_STATES_NAME);
        break;

      case MOL_POS:
        dreamm_v3_remove_file(DREAMM_VOL_MOL_POS_NAME);
        dreamm_v3_remove_file(DREAMM_VOL_MOL_STATES_NAME);
        dreamm_v3_remove_file(DREAMM_SURF_MOL_POS_NAME);
        dreamm_v3_remove_file(DREAMM_SURF_MOL_STATES_NAME);
        break;

      case MOL_ORIENT:
        dreamm_v3_remove_file(DREAMM_VOL_MOL_ORIENT_NAME);
        dreamm_v3_remove_file(DREAMM_SURF_MOL_ORIENT_NAME);
        break;

      case ALL_MESH_DATA:
        dreamm_v3_remove_file(DREAMM_MESH_POS_NAME);
        dreamm_v3_remove_file(DREAMM_MESH_STATES_NAME);
        dreamm_v3_remove_file(DREAMM_REGION_VIZ_DATA_NAME);
        break;

      case MESH_GEOMETRY:
        dreamm_v3_remove_file(DREAMM_MESH_POS_NAME);
        dreamm_v3_remove_file(DREAMM_MESH_STATES_NAME);
        break;

      case REG_DATA:
        dreamm_v3_remove_file(DREAMM_REGION_VIZ_DATA_NAME);
        break;

      default:
        fprintf(world->err_file, "File %s, Line %ld: Unexpected frame type in DREAMM V3 output mode (type=%d).\n", __FILE__, (long)__LINE__, fdlp->type);
        return 1;
    }
  }

  return 0;
}

/*************************************************************************
dreamm_v3_update_last_iteration_info:
    Update the "last iteration" information for meshes and molecules.  This
    information is used to create the symlinks needed in the DREAMM V3 output
    format.

        In: struct frame_data_list *fdlp - the frame data list
        Out: 0 on success, 1 on failure
**************************************************************************/
static int dreamm_v3_update_last_iteration_info(struct frame_data_list *fdlp)
{
  for (;
       fdlp != NULL;
       fdlp = fdlp->next)
  {
    if (world->it_time !=fdlp->viz_iteration)
      continue;

    if ((fdlp->type == ALL_MESH_DATA) ||
        (fdlp->type == REG_DATA) ||
        (fdlp->type == MESH_GEOMETRY))
    {
      if (fdlp->viz_iteration > world->viz_state_info.last_meshes_iteration)
        world->viz_state_info.last_meshes_iteration = fdlp->viz_iteration;
    }
    else if ((fdlp->type == ALL_MOL_DATA) ||
             (fdlp->type == MOL_POS) ||
             (fdlp->type == MOL_ORIENT))
    {
      if (fdlp->viz_iteration > world->viz_state_info.last_mols_iteration)
        world->viz_state_info.last_mols_iteration = fdlp->viz_iteration;
    }
  }

  return 0;
}

/*************************************************************************
dreamm_v3_create_empty_file:
    Create an empty file in the DREAMM V3 output directory.  This may be
    redundant now, since I'm explicitly removing files which, if left over from
    a previous run, would cause trouble.

        In: char const *fname - filename for new file
        Out: 0 upon success, 1 upon failure
**************************************************************************/
static int dreamm_v3_create_empty_file(char const *fname)
{
  FILE *f = NULL;
  char *path = NULL;
  if (dreamm_v3_remove_file(fname))
    return 1;

  if ((path = alloc_sprintf("%s/%s", world->viz_state_info.iteration_number_dir, fname)) == NULL)
    return 1;

  if ((f = fopen(path, "w")) == NULL)
  {
    free(path);
    return 1;
  }

  fclose(f);
  free(path);
  return 0;
}

/*************************************************************************
dreamm_v3_create_symlink:
    Create a symlink to a previous iteration appropriate for the DREAMM data
    format.

        In:  long long newiter - iteration needing symlink
             long long lastiter - iteration to which to make symlink
             char const *filename - filename to link from old dir
        Out: 0 upon succes, 1 upon failure.
**************************************************************************/
static int dreamm_v3_create_symlink(long long newiter,
                                    long long lastiter,
                                    char const *filename)
{
  char *newpath = NULL;
  char *oldpath = NULL;
  char *effoldpath = NULL;
  struct stat f_stat;

  /* Create old path for 'stat' */
  effoldpath = alloc_sprintf("%s/iteration_%lld/%s", world->viz_state_info.frame_data_dir, lastiter, filename);
  if (effoldpath == NULL)
    goto failure;

  /* If the old file doesn't exist, don't make the symlink. */
  if (stat(effoldpath, &f_stat) != 0)
  {
    free(effoldpath);
    return 0;
  }

  /* Old path name is a relative path */
  oldpath = alloc_sprintf("../iteration_%lld/%s", lastiter, filename);
  if (oldpath == NULL)
    goto failure;

  /* New path name is in new iteration directory */
  newpath = alloc_sprintf("%s/iteration_%lld/%s", world->viz_state_info.frame_data_dir, newiter, filename);
  if (newpath == NULL)
    goto failure;

  /* Try to create symlink */
  if (symlink(oldpath, newpath) == -1)
  {
    /* If we failed because the file exists */
    if (errno == EEXIST)
    {
      /* Remove the symlink and try again */
      if (unlink(newpath) == -1)
      {
        fprintf(world->err_file, "File %s, Line %ld: error %d removing symlink %s.\n", __FILE__, (long)__LINE__, errno, newpath);
        goto failure;
      }

      /* Try again to create symlink */
      if (symlink(oldpath, newpath) == -1)
      {
        fprintf(world->err_file, "File %s, Line %ld: error %d creating symlink %s.\n", __FILE__, (long)__LINE__, errno, newpath);
        goto failure;
      }
    }

    /* If we failed to make the symlink for some other reason */
    else
    {
      fprintf(world->err_file, "File %s, Line %ld: error %d creating symlink %s.\n", __FILE__, (long)__LINE__, errno, newpath);
      goto failure;
    }
  }

  /* Clean up and return success */
  free(effoldpath);
  free(oldpath);
  free(newpath);
  return 0;

failure:
  if (effoldpath != NULL) free(effoldpath);
  if (oldpath != NULL) free(oldpath);
  if (newpath != NULL) free(newpath);
  return 1;
}

/*************************************************************************
dreamm_v3_scan_for_mol_frames
    Scan a list of frame data objects to see if there is an upcoming mol frame
    for the specified iteration number.

        In: struct frame_data_list *fdlp - the list to scan
            long long iterno - the iteration number
        Out: 1 if a matching frame is found, 0 otherwise.
**************************************************************************/
static int dreamm_v3_scan_for_mol_frames(struct frame_data_list *fdlp,
                                         long long iterno)
{
  for (; fdlp != NULL; fdlp = fdlp->next)
  {
    /* XXX: This used to return as soon as it found an fdlp->viz_iteration
     * unequal to iterno...  I think that was incorrect.
     */
    if (fdlp->viz_iteration != iterno)
      continue;

    if ((fdlp->type == ALL_MOL_DATA) ||
        (fdlp->type == MOL_POS)      ||
        (fdlp->type == MOL_ORIENT))
      return 1;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_create_molecule_symlinks:
    If appropriate, creates molecule symlinks for a given iteration of DREAMM
    V3 output.

        In:  struct frame_data_list const *fdlp - the frame for which to create
                                                 links
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_create_molecule_symlinks(struct frame_data_list const *fdlp)
{
  long long lastiter = world->viz_state_info.last_mols_iteration;
  int mol_frame_found  = dreamm_v3_scan_for_mol_frames(fdlp->next, fdlp->viz_iteration);

  if (!mol_frame_found         &&                   /* No upcoming mol frames this iteration */
      lastiter >= 0            &&                   /* There is an old mol frame */
      fdlp->viz_iteration > lastiter)               /* The old mol frame was not created during this iteration */
  {
    /* Create a whole mess of symlinks */
    if (dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_SURF_MOL_HEADER_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_SURF_MOL_ORIENT_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_SURF_MOL_POS_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_SURF_MOL_STATES_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_VOL_MOL_HEADER_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_VOL_MOL_ORIENT_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_VOL_MOL_POS_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_VOL_MOL_STATES_NAME))
      return 1;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_scan_for_mesh_frames
    Scan a list of frame data objects to see if there is an upcoming mesh frame
    for the specified iteration number.

        In: struct frame_data_list *fdlp - the list to scan
            long long iterno - the iteration number
        Out: 1 if a matching frame is found, 0 otherwise.
**************************************************************************/
static int dreamm_v3_scan_for_mesh_frames(struct frame_data_list *fdlp,
                                          long long iterno)
{
  for (; fdlp != NULL; fdlp = fdlp->next)
  {
    if (fdlp->viz_iteration != iterno)
      continue;

    if ((fdlp->type == ALL_MESH_DATA) ||
        (fdlp->type == MESH_GEOMETRY) ||
        (fdlp->type == REG_DATA))
      return 1;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_create_mesh_symlinks:
    If appropriate, creates mesh symlinks for a given iteration of DREAMM V3
    output.

        In:  struct frame_data_list const *fdlp - the frame for which to create
                                                 links
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_create_mesh_symlinks(struct frame_data_list const *fdlp)
{
  long long lastiter = world->viz_state_info.last_meshes_iteration;
  int mesh_frame_found = dreamm_v3_scan_for_mesh_frames(fdlp->next, fdlp->viz_iteration);

  if (! mesh_frame_found        &&                  /* No more mesh frames this iteration */
      lastiter >= 0 &&                              /* There is an old meshes frame */
      fdlp->viz_iteration > lastiter)               /* The old meshes frame was not created during this iteration */
  {
    /* Create a whole mess of symlinks */
    if (dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_MESHES_HEADER_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_MESH_POS_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_MESH_STATES_NAME)
        || dreamm_v3_create_symlink(fdlp->viz_iteration, lastiter, DREAMM_REGION_VIZ_DATA_NAME))
      return 1;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_write_time_info:
    Writes the master header on the very last frame of the very last iteration
    of a DREAMM V3 viz output run.

        In:  char const *viz_data_dir - dir for master header
             char const *master_header_name - name of master header
             char const *iteration_numbers_name - name of iteration data
             char const *time_values_name - name of time data
             u_int iteration_numbers_count - number of iteration data
             u_int time_values_count - number of time data
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_write_time_info(char const *viz_data_dir,
                                     char const *master_header_name,
                                     char const *iteration_numbers_name,
                                     char const *time_values_name,
                                     u_int iteration_numbers_count,
                                     u_int time_values_count)
{
  FILE *master_header = NULL;

  /* Open master header file. */
  if ((master_header = dreamm_v3_generic_open_file(viz_data_dir, master_header_name, "w")) == NULL)
    goto failure;

  if (dreamm_v3_generic_write_time_info(master_header,
                                        iteration_numbers_name,
                                        time_values_name,
                                        "DREAMM_V3_MODE",
                                        1,
                                        iteration_numbers_count,
                                        time_values_count))
    goto failure;

  if (master_header) fclose(master_header);
  return 0;

failure:
  if (master_header) fclose(master_header);
  return 1;
}

/*************************************************************************
dreamm_v3_make_time_info_filename:
    Make output filename for DREAMM V3 (ungrouped) time/iteration data.

        In:  char const *typename - "iteration_numbers" or "time_values"
        Out: the filename, or NULL if allocation fails
**************************************************************************/
static char *dreamm_v3_make_time_info_filename(char const *typename)
{
  char *filename = NULL;
  if (world->chkpt_flag)
    filename = alloc_sprintf("%s.%s.%d.bin",
                             world->viz_state_info.filename_prefix_basename,
                             typename,
                             world->chkpt_seq_num);
  else
    filename = alloc_sprintf("%s.%s.bin",
                             typename,
                             world->viz_state_info.filename_prefix_basename);
  if (filename == NULL)
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);

  return filename;
}

/*************************************************************************
dreamm_v3_dump_time_info:
    Writes the master header, as well as the iteration and time data files on
    the very last frame of the very last iteration of a DREAMM V3 viz output
    run.

        In:  None
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_dump_time_info(void)
{
  char *time_values_name = NULL;
  char *iteration_numbers_name = NULL;
  char *master_header_name = NULL;
  long long iteration_numbers_count = 0;

  /* Build viz data dir name */
  char *viz_data_dir = my_strcat(world->file_prefix_name, "_viz_data");
  if (viz_data_dir == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    goto failure;
  }

  /* Build iteration numbers filename */
  if ((iteration_numbers_name = dreamm_v3_make_time_info_filename("iteration_numbers")) == NULL)
    goto failure;

  /* Build time values filename */
  if ((time_values_name = dreamm_v3_make_time_info_filename("time_values")) == NULL)
    goto failure;

  /* Build master header filename */
  if (world->chkpt_flag)
    master_header_name = alloc_sprintf("%s.%u.dx",
                                       world->viz_state_info.filename_prefix_basename,
                                       world->chkpt_seq_num);
  else
    master_header_name = alloc_sprintf("%s.dx",
                                       world->viz_state_info.filename_prefix_basename);
  if (master_header_name == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    goto failure;
  }

  /* Find maximum iteration numbers count */
  if (world->viz_state_info.mesh_output_iterations.n_iterations >
      world->viz_state_info.grid_mol_output_iterations.n_iterations)
    iteration_numbers_count = world->viz_state_info.mesh_output_iterations.n_iterations;
  else
    iteration_numbers_count = world->viz_state_info.grid_mol_output_iterations.n_iterations;
  if (world->viz_state_info.vol_mol_output_iterations.n_iterations > iteration_numbers_count)
    iteration_numbers_count = world->viz_state_info.vol_mol_output_iterations.n_iterations;

  /* Write iteration numbers file */
  if (dreamm_v3_generic_dump_iteration_numbers(viz_data_dir,
                                               iteration_numbers_name,
                                               iteration_numbers_count))
    goto failure;

  /* write "time_values" object. */
  if (world->viz_state_info.output_times.n_iterations > 0  &&
      dreamm_v3_generic_dump_time_values(viz_data_dir, time_values_name))
    goto failure;
  if (dreamm_v3_write_time_info(viz_data_dir,
                                master_header_name,
                                iteration_numbers_name,
                                time_values_name,
                                iteration_numbers_count,
                                world->viz_state_info.output_times.n_iterations))
    goto failure;

  free(viz_data_dir);
  free(iteration_numbers_name);
  free(time_values_name);
  free(master_header_name);
  return 0;

failure:
  if (viz_data_dir) free(viz_data_dir);
  if (iteration_numbers_name) free(iteration_numbers_name);
  if (time_values_name) free(time_values_name);
  if (master_header_name) free(master_header_name);
  return 1;
}

/*************************************************************************
dreamm_v3_write_mesh_group:
    Write a group containing all meshes to the header file.

        In:  FILE *master_header - file to which to write group
             int groupidx - the group index
             int field_idx_base - object number for first field
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_write_mesh_group(FILE *header,
                                      int groupidx,
                                      int fieldsidx)
{
  int obj_index;

  /* Create group object */
  fprintf(header,
          "object %d group # meshes #\n",
          groupidx);

  /* Add member declarations to group */
  for (obj_index = 0; obj_index < world->viz_state_info.n_viz_objects; ++ obj_index)
  {
    fprintf(header,
            "\tmember \"%s\" value %d\n",
            world->viz_state_info.viz_objects[obj_index]->sym->name,
            fieldsidx ++);
  }
  fprintf(header, "\n");

  return 0;
}

/*************************************************************************
dreamm_v3_dump_mesh_data:
    Dump all mesh data to the mesh output files for this iteration, writing
    index information to the header.

        In:  struct frame_data_list *fdlp - frame for which to write meshes
             FILE *meshes_header - header to receive index info
             char const *iteration_dir - directory for this iteration
             int *meshes_main_index - ptr to index for allocating obj numbers
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_dump_mesh_data(struct frame_data_list const * const fdlp,
                                    FILE *meshes_header,
                                    char const *iteration_dir,
                                    int *meshes_main_index)
{
  return dreamm_v3_generic_dump_mesh_data(fdlp,
                                          meshes_header,
                                          iteration_dir,
                                          DREAMM_MESH_POS_NAME,
                                          DREAMM_MESH_STATES_NAME,
                                          DREAMM_REGION_VIZ_DATA_NAME,
                                          meshes_main_index);
}

/*************************************************************************
dreamm_v3_dump_meshes:
    Dump all mesh data to the mesh output files for this iteration.

        In:  struct frame_data_list *fdlp - frame for which to write meshes
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_dump_meshes(struct frame_data_list const * const fdlp)
{
  int meshes_main_index = 1;
  FILE *meshes_header = NULL;

  /* If desired, open meshes header */
  if ((meshes_header = dreamm_v3_generic_open_file(world->viz_state_info.iteration_number_dir,
                                                   DREAMM_MESHES_HEADER_NAME, "w")) == NULL)
    return 1;

  /* Dump mesh info */
  if (dreamm_v3_dump_mesh_data(fdlp,
                               meshes_header,
                               world->viz_state_info.iteration_number_dir,
                               &meshes_main_index))
    goto failure;

  /* Create field objects for meshes */
  if (world->viz_state_info.n_viz_objects > 0)
  {
    int field_idx_base = meshes_main_index;
    meshes_main_index += world->viz_state_info.n_viz_objects;

    int group_idx = meshes_main_index ++;

    if (dreamm_v3_generic_write_mesh_fields(fdlp,
                                            meshes_header,
                                            field_idx_base,
                                            1))
      goto failure;

    /* Create a group object for all meshes */
    if (dreamm_v3_write_mesh_group(meshes_header, group_idx, field_idx_base))
      goto failure;

    /* Store iteration_number for meshes */
    if (add_to_iteration_counter(&world->viz_state_info.mesh_output_iterations, fdlp->viz_iteration))
      goto failure;

    /* Put value of viz_iteration into the time values */
    if (add_to_iteration_counter_monotonic(&world->viz_state_info.output_times, fdlp->viz_iteration))
      goto failure;
  }

  if (meshes_header) fclose(meshes_header);
  return 0;

failure:
  if (meshes_header) fclose(meshes_header);
  return 0;
}

/*************************************************************************
dreamm_v3_write_molecule_group:
    Write a group containing all molecules of the specified "type" (either
    "surface" or "volume") to the header file.

        In:  FILE *master_header - file to which to write group
             char const *desc - type of mol ("surface molecules" or "volume
                                molecules")
             struct species **all_species - species in group
             int n_species - number of species in group
             int field_idx_base - object number for first field
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_write_molecule_group(FILE *header,
                                          char const *desc,
                                          struct species **specs,
                                          int num_specs,
                                          int groupidx,
                                          int fieldsidx)
{
  int species_index;

  /* Create group object */
  fprintf(header,
          "object %d group # %s #\n",
          groupidx,
          desc);

  /* Add member declarations to group */
  for (species_index = 0; species_index < num_specs; ++ species_index)
    fprintf(header,
            "\tmember \"%s\" value %d\n",
            specs[species_index]->sym->name,
            fieldsidx ++);
  fprintf(header, "\n");

  return 0;
}

/*************************************************************************
dreamm_v3_dump_grid_molecule_data:
    Dump desired data for all grid molecules to grid molecule files, writing
    index information to the header file for this iteration.

        In:  struct frame_data_list *fdlp - frame for which to write mols
             FILE *surf_mol_header - the header file for grid molecules
             char const *iteration_number_dir - directory for this iteration
             int *surf_mol_main_index - pointer to counter of object numbers
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_dump_grid_molecule_data(struct frame_data_list const * const fdlp,
                                             FILE *surf_mol_header,
                                             char const *iteration_number_dir,
                                             int *surf_mol_main_index)
{
  return dreamm_v3_generic_dump_grid_molecule_data(fdlp,
                                                   surf_mol_header,
                                                   iteration_number_dir,
                                                   DREAMM_SURF_MOL_POS_NAME,
                                                   DREAMM_SURF_MOL_ORIENT_NAME,
                                                   DREAMM_SURF_MOL_STATES_NAME,
                                                   surf_mol_main_index);
}

/*************************************************************************
dreamm_v3_dump_grid_molecules:
    Write the desired grid molecule data to the grid molecule data files for
    this iteration.

        In:  struct frame_data_list *fdlp - frame for which to write mols
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_dump_grid_molecules(struct frame_data_list const * const fdlp)
{
  int surf_mol_main_index = 1;

  FILE  *surf_mol_header = NULL;

  /* Open surface molecules header file */
  if ((surf_mol_header = dreamm_v3_generic_open_file(world->viz_state_info.iteration_number_dir,
                                                     DREAMM_SURF_MOL_HEADER_NAME, "w")) == NULL)
    return 1;

  if (dreamm_v3_dump_grid_molecule_data(fdlp,
                                        surf_mol_header,
                                        world->viz_state_info.iteration_number_dir,
                                        &surf_mol_main_index))
    goto failure;

  if (world->viz_state_info.n_grid_species > 0)
  {
    int field_idx_base = surf_mol_main_index;
    surf_mol_main_index += world->viz_state_info.n_grid_species;

    int group_idx = surf_mol_main_index ++;

    /* Build fields for grid molecules here */
    if (dreamm_v3_generic_write_molecule_fields(fdlp,
                                                surf_mol_header,
                                                world->viz_state_info.grid_species,
                                                world->viz_state_info.n_grid_species,
                                                field_idx_base,
                                                1))
      goto failure;
    if (dreamm_v3_write_molecule_group(surf_mol_header,
                                       "surface molecules",
                                       world->viz_state_info.grid_species,
                                       world->viz_state_info.n_grid_species,
                                       group_idx,
                                       field_idx_base))
      goto failure;

    /* Store iteration_number for surface molecules */
    if (add_to_iteration_counter(&world->viz_state_info.grid_mol_output_iterations, fdlp->viz_iteration))
      goto failure;

    /* Put value of viz_iteration into the time values */
    if (add_to_iteration_counter_monotonic(&world->viz_state_info.output_times, fdlp->viz_iteration))
      goto failure;
  }

  if (surf_mol_header) fclose(surf_mol_header);
  return 0;

failure:
  if (surf_mol_header) fclose(surf_mol_header);
  return 1;
}

/*************************************************************************
dreamm_v3_dump_volume_molecule_data:
    Dump desired data for all volume molecules to volume molecule files,
    writing index information to the header file for this iteration.

        In:  struct frame_data_list *fdlp - frame for which to write mols
             FILE *vol_mol_header - the header file for volume molecules
             char const *iteration_number_dir - directory for this iteration
             int *vol_mol_main_index - pointer to counter of object numbers
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_dump_volume_molecule_data(struct frame_data_list const * const fdlp,
                                               FILE *vol_mol_header,
                                               char const *iteration_number_dir,
                                               int *vol_mol_main_index)
{
  return dreamm_v3_generic_dump_volume_molecule_data(fdlp,
                                                     vol_mol_header,
                                                     iteration_number_dir,
                                                     DREAMM_VOL_MOL_POS_NAME,
                                                     DREAMM_VOL_MOL_ORIENT_NAME,
                                                     DREAMM_VOL_MOL_STATES_NAME,
                                                     vol_mol_main_index);
}

/*************************************************************************
dreamm_v3_dump_volume_molecules:
    Write the desired volume molecule data to the volume molecule data files
    for this iteration.

        In:  struct frame_data_list *fdlp - frame for which to write mols
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_dump_volume_molecules(struct frame_data_list const * const fdlp)
{
  int vol_mol_main_index = 1;
  FILE *vol_mol_header = NULL;

  /* Open volume molecules header */
  if ((vol_mol_header = dreamm_v3_generic_open_file(world->viz_state_info.iteration_number_dir,
                                                    DREAMM_VOL_MOL_HEADER_NAME, "w")) == NULL)
    return 1;

  if (dreamm_v3_dump_volume_molecule_data(fdlp,
                                          vol_mol_header,
                                          world->viz_state_info.iteration_number_dir,
                                          &vol_mol_main_index))
    goto failure;

  if (world->viz_state_info.n_vol_species > 0)
  {
    int field_idx_base = vol_mol_main_index;
    vol_mol_main_index += world->viz_state_info.n_vol_species;

    int group_idx = vol_mol_main_index ++;

    /* Build fields for grid molecules here */
    if (dreamm_v3_generic_write_molecule_fields(fdlp,
                                                vol_mol_header,
                                                world->viz_state_info.vol_species,
                                                world->viz_state_info.n_vol_species,
                                                field_idx_base,
                                                1))
      goto failure;

    /* Create group objects for molecules */
    if (dreamm_v3_write_molecule_group(vol_mol_header,
                                       "volume molecules",
                                       world->viz_state_info.vol_species,
                                       world->viz_state_info.n_vol_species,
                                       group_idx,
                                       field_idx_base))
      goto failure;

    /* Store iteration_number for volume_molecules */
    if (add_to_iteration_counter(&world->viz_state_info.vol_mol_output_iterations, fdlp->viz_iteration))
      goto failure;

    /* Put value of viz_iteration into the time values */
    if (add_to_iteration_counter_monotonic(&world->viz_state_info.output_times, fdlp->viz_iteration))
      goto failure;
  }

  if (vol_mol_header) fclose(vol_mol_header);
  return 0;

failure:
  if (vol_mol_header) fclose(vol_mol_header);
  return 1;
}

/*************************************************************************
output_dreamm_objects:
        In: struct frame_data_list *fdlp
        Out: 0 on success, 1 on error; output visualization files (*.dx)
             in dreamm group format are written.
**************************************************************************/
int output_dreamm_objects(struct frame_data_list const * const fdlp)
{
  /* What do we output? */
  byte viz_meshes = 0;
  byte viz_mols = 0;

  no_printf("Viz output in DREAMM_V3 mode...\n");

  /* Set flags based on fdlp->type */
  if (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_POS  ||  fdlp->type == MOL_ORIENT)
    viz_mols = 1;
  if (fdlp->type == ALL_MESH_DATA  ||  fdlp->type == MESH_GEOMETRY  ||  fdlp->type == REG_DATA)
    viz_meshes = 1;

  /* Create empty files, if appropriate */
  /* XXX: Do we still need to do this?  We're not creating links to these files, so... */
  if (world->viz_state_info.n_viz_objects == 0)
  {
    dreamm_v3_create_empty_file(DREAMM_MESH_POS_NAME);
    dreamm_v3_create_empty_file(DREAMM_MESH_STATES_NAME);
    dreamm_v3_create_empty_file(DREAMM_REGION_VIZ_DATA_NAME);
    dreamm_v3_create_empty_file(DREAMM_MESHES_HEADER_NAME);
  }
  if (world->viz_state_info.n_vol_species == 0  &&  world->viz_state_info.n_grid_species == 0)
  {
    dreamm_v3_create_empty_file(DREAMM_VOL_MOL_POS_NAME);
    dreamm_v3_create_empty_file(DREAMM_VOL_MOL_ORIENT_NAME);
    dreamm_v3_create_empty_file(DREAMM_VOL_MOL_STATES_NAME);
    dreamm_v3_create_empty_file(DREAMM_VOL_MOL_HEADER_NAME);

    dreamm_v3_create_empty_file(DREAMM_SURF_MOL_POS_NAME);
    dreamm_v3_create_empty_file(DREAMM_SURF_MOL_ORIENT_NAME);
    dreamm_v3_create_empty_file(DREAMM_SURF_MOL_STATES_NAME);
    dreamm_v3_create_empty_file(DREAMM_SURF_MOL_HEADER_NAME);
  }

  /* Dump meshes */
  if (viz_meshes  &&  dreamm_v3_dump_meshes(fdlp))
      return 1;

  /* Dump molecules */
  if (viz_mols)
  {
    /* Dump grid molecules. */
    if (dreamm_v3_dump_grid_molecules(fdlp))
      return 1;

    /* dump 3D molecules: */
    if (dreamm_v3_dump_volume_molecules(fdlp))
      return 1;
  }

  /* If we're on the final iteration... */
  if (world->it_time == world->viz_state_info.final_iteration)
  {
    /* Check for other frames this iteration */
    if (! dreamm_v3_generic_scan_for_frame(fdlp->next, fdlp->viz_iteration) &&
        dreamm_v3_dump_time_info())
      return 1;
  }

  /*
   * Now check whether there is a need to create symlinks to the
   * meshes/molecules files saved previously in the previous frame directories.
   * Create link in this "frame_#" folder to the last existing "meshes" files.
   */
  if (viz_meshes  &&  dreamm_v3_create_molecule_symlinks(fdlp))
    return 1;

  if (viz_mols  &&  dreamm_v3_create_mesh_symlinks(fdlp))
    return 1;

  return 0;
}

/* == DREAMM V3 Output (grouped) == */

/* Filename base for DREAMM V3 (grouped) output */
static char const * const DREAMM_GROUPED_MOL_POS_NAME = "molecule_positions";
static char const * const DREAMM_GROUPED_MOL_ORIENT_NAME = "molecule_orientations";
static char const * const DREAMM_GROUPED_MOL_STATES_NAME = "molecule_states";
static char const * const DREAMM_GROUPED_MESH_POS_NAME = "mesh_positions";
static char const * const DREAMM_GROUPED_MESH_STATES_NAME = "mesh_states";
static char const * const DREAMM_GROUPED_REGION_VIZ_DATA_NAME = "region_indices";

/*************************************************************************
dreamm_v3_grouped_get_master_header_name:
    Get the name of the master header file for DREAMM V3 grouped output.

        In:  None
        Out: The path, or NULL if an error occurs
**************************************************************************/
static char *dreamm_v3_grouped_get_master_header_name()
{
  char *master_header_file_path = NULL;
  if (world->chkpt_flag)
    master_header_file_path = alloc_sprintf("%s.%d.dx", world->file_prefix_name, world->chkpt_seq_num);
  else
    master_header_file_path = alloc_sprintf("%s.dx", world->file_prefix_name);

  if (master_header_file_path == NULL)
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);

  return master_header_file_path;
}

/*************************************************************************
dreamm_v3_grouped_create_filepath:
    Creates a filepath for a given type of output file.

        In:  char const *kind - the type of output file
             char **path - pointer to receive allocated filepath
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_grouped_create_filepath(char const *kind,
                                             char **path)
{
  if (world->chkpt_flag)
    *path = alloc_sprintf("%s.%s.%d.bin",
                          world->file_prefix_name,
                          kind,
                          world->chkpt_seq_num);
  else
    *path = alloc_sprintf("%s.%s.bin",
                          world->file_prefix_name,
                          kind);

  if (*path == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_grouped_create_filename:
    Creates a filename for a given type of output file.

        In:  char const *kind - the type of output file
             char **path - pointer to receive allocated filename
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_grouped_create_filename(char const *kind,
                                             char **name)
{
  if (world->chkpt_flag)
    *name = alloc_sprintf("%s.%s.%d.bin",
                          world->viz_state_info.filename_prefix_basename,
                          kind,
                          world->chkpt_seq_num);
  else
    *name = alloc_sprintf("%s.%s.bin",
                          world->viz_state_info.filename_prefix_basename,
                          kind);

  if (*name == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_grouped_remove_file:
    Remove a file for DREAMM V3 Grouped output.

        In:  char const *kind - the "type" of file to remove (used to build
                                filename)
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_grouped_remove_file(char const *kind)
{
  char *filename = NULL;
  if (dreamm_v3_grouped_create_filepath(kind, &filename))
    return 1;

  if (unlink(filename)  &&  errno != ENOENT)
  {
    fprintf(world->err_file, "File %s, Line %ld: failed to unlink %s file: %s\n", __FILE__, (long)__LINE__, kind, filename);
    free(filename);
    return 1;
  }

  free(filename);
  return 0;
}

/*************************************************************************
dreamm_v3_grouped_clean_files:
    Clean files for DREAMM V3 Grouped output.  This will remove files which are
    going to be created by the DREAMM V3 Grouped output routines.  This should
    be called before the first frame of output is processed.

        In:  struct frame_data_list *fdlp - the head of the frame data list
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_grouped_clean_files(struct frame_data_list *fdlp)
{
  /* Delete master header */
  char *filename = dreamm_v3_grouped_get_master_header_name();
  if (filename == NULL)
    return 1;
  if (unlink(filename)  &&  errno != ENOENT)
  {
    fprintf(world->err_file, "File %s, Line %ld: failed to unlink master header file: %s\n", __FILE__, (long)__LINE__, filename);
    free(filename);
    return 1;
  }
  free(filename);

  /* Delete any subfiles which we'll be creating */
  for (; fdlp != NULL; fdlp = fdlp->next)
  {
    switch (fdlp->type)
    {
      case ALL_MOL_DATA:
        if (dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MOL_POS_NAME)     ||
            dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MOL_ORIENT_NAME)  ||
            dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MOL_STATES_NAME))
          return 1;
        break;

      case MOL_POS:
        if (dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MOL_POS_NAME)  ||
            dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MOL_STATES_NAME))
          return 1;
        break;

      case MOL_ORIENT:
        if (dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MOL_ORIENT_NAME))
          return 1;
        break;

      case ALL_MESH_DATA:
        if (dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MESH_POS_NAME)  ||
            dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MESH_STATES_NAME)     ||
            dreamm_v3_grouped_remove_file(DREAMM_GROUPED_REGION_VIZ_DATA_NAME))
          return 1;
        break;

      case MESH_GEOMETRY:
        if (dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MESH_POS_NAME)  ||
            dreamm_v3_grouped_remove_file(DREAMM_GROUPED_MESH_STATES_NAME))
          return 1;
        break;

      case REG_DATA:
        if (dreamm_v3_grouped_remove_file(DREAMM_GROUPED_REGION_VIZ_DATA_NAME))
          return 1;
        break;

      default:
        fprintf(world->err_file, "File %s, Line %ld: Unexpected frame type in DREAMM V3 grouped output mode (type=%d).\n", __FILE__, (long)__LINE__, fdlp->type);
        return 1;
    }
  }

  return 0;
}

/*************************************************************************
dreamm_v3_grouped_init:
    Initialize state for DREAMM_V3_GROUPED output.

        In:  struct frame_data_list *fdlp - the head of the frame data list
        Out: 0 if successful, 1 if failed
**************************************************************************/
static int dreamm_v3_grouped_init(struct frame_data_list *fdlp)
{
  if (dreamm_v3_generic_init())
    return 1;

  if (dreamm_v3_grouped_clean_files(fdlp))
    return 1;

  world->viz_state_info.dx_main_object_index = 1;
  world->viz_state_info.dreamm_last_iteration_meshes = -1;
  world->viz_state_info.dreamm_last_iteration_vol_mols = -1;
  world->viz_state_info.dreamm_last_iteration_surf_mols = -1;

  int time_values_total = count_time_values(fdlp);
  if (reset_time_values(fdlp, world->start_time))
    return 1;
  if (initialize_iteration_counters(time_values_total))
    return 1;

  if (initialize_string_buffer(&world->viz_state_info.combined_group_members,
                               time_values_total))
    return 1;

  return 0;
}

/*************************************************************************
dreamm_v3_grouped_write_time_info:
    Writes time index information to the master header on the very last frame
    of the very last iteration of a DREAMM V3 Grouped viz output run.

        In:  FILE *master_header - header file for time index info
             char const *iteration_numbers_name - name of iteration data
             char const *time_values_name - name of time data
             u_int iteration_numbers_count - number of iteration data
             u_int time_values_count - number of time data
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_write_time_info(FILE *master_header,
                                             char const *iteration_numbers_name,
                                             char const *time_values_name,
                                             u_int iteration_numbers_count,
                                             u_int time_values_count)
{
  return dreamm_v3_generic_write_time_info(master_header,
                                           iteration_numbers_name,
                                           time_values_name,
                                           "DREAMM_V3_GROUPED_MODE",
                                           2,
                                           iteration_numbers_count,
                                           time_values_count);
}

/*************************************************************************
dreamm_v3_grouped_dump_time_info:
    Writes time info to the master header, and creates the iteration and time
    data files on the very last frame of the very last iteration of a DREAMM V3
    Grouped viz output run.

        In:  FILE *master_header - the master header file
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_dump_time_info(FILE *master_header)
{
  char *time_values_name = NULL;
  char *iteration_numbers_name = NULL;
  long long iteration_numbers_count = 0;

  /* Create filenames */
  if (dreamm_v3_grouped_create_filename("iteration_numbers",
                                        &iteration_numbers_name))
    goto failure;
  if (dreamm_v3_grouped_create_filename("time_values",
                                        &time_values_name))
    goto failure;

  /* Find maximum iteration numbers count */
  if (world->viz_state_info.mesh_output_iterations.n_iterations > world->viz_state_info.grid_mol_output_iterations.n_iterations)
    iteration_numbers_count = world->viz_state_info.mesh_output_iterations.n_iterations;
  else
    iteration_numbers_count = world->viz_state_info.grid_mol_output_iterations.n_iterations;
  if (world->viz_state_info.vol_mol_output_iterations.n_iterations > iteration_numbers_count)
    iteration_numbers_count = world->viz_state_info.vol_mol_output_iterations.n_iterations;

  /* Write iteration numbers file */
  if (dreamm_v3_generic_dump_iteration_numbers(world->viz_state_info.filename_prefix_dirname,
                                               iteration_numbers_name,
                                               iteration_numbers_count))
    goto failure;

  /* write "time_values" object. */
  if (world->viz_state_info.output_times.n_iterations > 0  &&
      dreamm_v3_generic_dump_time_values(world->viz_state_info.filename_prefix_dirname, time_values_name))
    goto failure;

  /* Write header details */
  if (dreamm_v3_grouped_write_time_info(master_header,
                                        iteration_numbers_name,
                                        time_values_name,
                                        iteration_numbers_count,
                                        world->viz_state_info.output_times.n_iterations))
    goto failure;

  if (iteration_numbers_name) free(iteration_numbers_name);
  if (time_values_name) free(time_values_name);
  return 0;

failure:
  if (iteration_numbers_name) free(iteration_numbers_name);
  if (time_values_name) free(time_values_name);
  return 0;
}

/*************************************************************************
dreamm_v3_grouped_write_mesh_group:
    Write a group containing all meshes to the header file.

        In:  FILE *master_header - file to which to write group
             struct frame_data_list *fdlp - frame for which to write group
             int field_idx_base - object number for first field
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_write_mesh_group(FILE *master_header,
                                              struct frame_data_list const * const fdlp,
                                              int field_idx_base)
{
  int obj_index;
  fprintf(master_header, "object \"meshes_%lld\" group # meshes #\n", fdlp->viz_iteration);
  for (obj_index = 0; obj_index < world->viz_state_info.n_viz_objects; ++obj_index)
    fprintf(master_header,
            "\tmember \"%s\" value %d\n",
            world->viz_state_info.viz_objects[obj_index]->sym->name,
            field_idx_base ++);
  fprintf(master_header, "\n");
  return 0;
}

/*************************************************************************
dreamm_v3_grouped_dump_mesh_data:
    Dump all mesh data to appropriate output files, storing index information
    in the header file.

        In:  struct frame_data_list *fdlp - frame for which to write meshes
             FILE *master_header - file to which to write group
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_dump_mesh_data(struct frame_data_list const * const fdlp,
                                            FILE *master_header)
{
  /* Control flags */
  byte viz_surf_pos_flag    = (fdlp->type == ALL_MESH_DATA || fdlp->type == MESH_GEOMETRY);
  byte viz_region_data_flag = (fdlp->type == ALL_MESH_DATA || fdlp->type == REG_DATA);
  byte viz_surf_states_flag = (viz_surf_pos_flag  &&  (world->viz_output_flag & VIZ_SURFACE_STATES) != 0);

  /* Filenames */
  char *mesh_pos_name = NULL;
  char *mesh_states_name = NULL;
  char *region_viz_data_name = NULL;

  /* Build filenames */
  if (viz_surf_pos_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_MESH_POS_NAME,
                                                               &mesh_pos_name))
    goto failure;
  if (viz_surf_states_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_MESH_STATES_NAME,
                                                                  &mesh_states_name))
    goto failure;
  if (viz_region_data_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_REGION_VIZ_DATA_NAME,
                                                                  &region_viz_data_name))
    goto failure;

  /* Output mesh info */
  if (dreamm_v3_generic_dump_mesh_data(fdlp,
                                       master_header,
                                       world->viz_state_info.filename_prefix_dirname,
                                       mesh_pos_name,
                                       mesh_states_name,
                                       region_viz_data_name,
                                       &world->viz_state_info.dx_main_object_index))
    goto failure;

  if (mesh_pos_name) free(mesh_pos_name);
  if (mesh_states_name) free(mesh_states_name);
  if (region_viz_data_name) free(region_viz_data_name);
  return 0;

failure:
  if (mesh_pos_name) free(mesh_pos_name);
  if (mesh_states_name) free(mesh_states_name);
  if (region_viz_data_name) free(region_viz_data_name);
  return 1;
}

/*************************************************************************
dreamm_v3_grouped_dump_meshes:
    Dump all mesh data to appropriate output files and write index information
    into the header file, creating fields and groups as appropriate.

        In:  struct frame_data_list *fdlp - frame for which to write meshes
             FILE *master_header - file to which to write group
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_dump_meshes(struct frame_data_list const * const fdlp, FILE *master_header)
{
  int surf_index = world->viz_state_info.dx_main_object_index;

  if (dreamm_v3_grouped_dump_mesh_data(fdlp,
                                       master_header))
    return 1;

  if (world->viz_state_info.n_viz_objects > 0)
  {
    /* Allocate ids for fields */
    int field_idx_base = world->viz_state_info.dx_main_object_index;
    world->viz_state_info.dx_main_object_index += world->viz_state_info.n_viz_objects;

    /* Create field objects */
    if (dreamm_v3_generic_write_mesh_fields(fdlp,
                                            master_header,
                                            field_idx_base,
                                            surf_index))
      return 1;

    /* Create a group object for all meshes */
    world->viz_state_info.dreamm_last_iteration_meshes = fdlp->viz_iteration;
    if (dreamm_v3_grouped_write_mesh_group(master_header, fdlp, field_idx_base))
      return 1;

    /* Store iteration_number for meshes */
    if (add_to_iteration_counter(&world->viz_state_info.mesh_output_iterations,
                                 fdlp->viz_iteration))
      return 1;
  }

  return 0;
}

/*************************************************************************
dreamm_v3_grouped_write_molecule_group:
    Write a group containing all molecules of the specified "type" (either
    "surface" or "volume") to the header file.

        In:  FILE *master_header - file to which to write group
             struct frame_data_list *fdlp - frame for which to write group
             char const *moltype - type of mol ("surface" or "volume")
             struct species **all_species - species in group
             int n_species - number of species in group
             int field_idx_base - object number for first field
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_write_molecule_group(FILE *master_header,
                                                  struct frame_data_list const * const fdlp,
                                                  char const *moltype,
                                                  struct species **all_species,
                                                  int n_species,
                                                  int field_idx_base)
{
  int mol_index;
  fprintf(master_header,
          "object \"%s_molecules_%lld\" group # %s molecules #\n",
          moltype,
          fdlp->viz_iteration,
          moltype);
  for (mol_index = 0; mol_index < n_species; ++ mol_index)
    fprintf(master_header,
            "\tmember \"%s\" value %d\n",
            all_species[mol_index]->sym->name,
            field_idx_base++);
  return 0;
}

/*************************************************************************
dreamm_v3_grouped_dump_grid_molecule_data:
    Dump desired data for all grid molecules to grid molecule files, writing
    index information to the header file for this iteration.

        In:  struct frame_data_list *fdlp - frame for which to write mols
             FILE *master_header - the header file for index info
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_dump_grid_molecule_data(struct frame_data_list const * const fdlp,
                                                     FILE *master_header)
{

  /* Control flags */
  byte viz_mol_pos_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_POS);
  byte viz_mol_orient_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_ORIENT);
  byte viz_mol_states_flag = (viz_mol_pos_flag  &&  (world->viz_output_flag & VIZ_MOLECULES_STATES));

  /* Filenames */
  char *mol_pos_name = NULL;
  char *mol_orient_name = NULL;
  char *mol_states_name = NULL;

  /* Build filenames */
  if (viz_mol_pos_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_MOL_POS_NAME,
                                                              &mol_pos_name))
    goto failure;
  if (viz_mol_orient_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_MOL_ORIENT_NAME,
                                                                 &mol_orient_name))
    goto failure;
  if (viz_mol_states_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_MOL_STATES_NAME,
                                                                 &mol_states_name))
    goto failure;

  /* Output molecule info */
  if (dreamm_v3_generic_dump_grid_molecule_data(fdlp,
                                                master_header,
                                                world->viz_state_info.filename_prefix_dirname,
                                                mol_pos_name,
                                                mol_orient_name,
                                                mol_states_name,
                                                &world->viz_state_info.dx_main_object_index))
    goto failure;

  if (mol_pos_name) free(mol_pos_name);
  if (mol_orient_name) free(mol_orient_name);
  if (mol_states_name) free(mol_states_name);
  return 0;

failure:
  if (mol_pos_name) free(mol_pos_name);
  if (mol_orient_name) free(mol_orient_name);
  if (mol_states_name) free(mol_states_name);
  return 1;
}

/*************************************************************************
dreamm_v3_grouped_dump_grid_molecules:
    Write the desired grid molecule data to the volume molecule data files and
    write appropriate index information into the master header.

        In:  struct frame_data_list *fdlp - frame for which to write mols
             FILE *master_header - file to which to write group
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_dump_grid_molecules(struct frame_data_list const * const fdlp,
                                                 FILE *master_header)
{
  int eff_index_base = world->viz_state_info.dx_main_object_index;
  if (dreamm_v3_grouped_dump_grid_molecule_data(fdlp,
                                                master_header))
    return 1;

  if (world->viz_state_info.n_grid_species > 0)
  {
    /* Allocate field indices for grid molecule data */
    int field_idx_base = world->viz_state_info.dx_main_object_index;
    world->viz_state_info.dx_main_object_index += world->viz_state_info.n_grid_species;

    /* Build fields for grid molecules here */
    if (dreamm_v3_generic_write_molecule_fields(fdlp,
                                                master_header,
                                                world->viz_state_info.grid_species,
                                                world->viz_state_info.n_grid_species,
                                                field_idx_base,
                                                eff_index_base))
      return 1;

    /* Create groups for effectors */
    world->viz_state_info.dreamm_last_iteration_surf_mols = fdlp->viz_iteration;
    if (dreamm_v3_grouped_write_molecule_group(master_header,
                                               fdlp,
                                               "surface",
                                               world->viz_state_info.grid_species,
                                               world->viz_state_info.n_grid_species,
                                               field_idx_base))
      return 1;

    /* Store iteration_number for surface molecules */
    if (add_to_iteration_counter(&world->viz_state_info.grid_mol_output_iterations,
                                 fdlp->viz_iteration))
      return 1;
  }
  fprintf(master_header, "\n"); 
  return 0;
}

/*************************************************************************
dreamm_v3_grouped_dump_volume_molecule_data:
    Dump desired data for all volume molecules to volume molecule files,
    writing index information to the header file for this iteration.  Object id
    numbers are allocated first to pos, then to orientation, then to state for
    each object.

        In:  struct frame_data_list *fdlp - frame for which to write mols
             FILE *master_header - the header file for index info
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_dump_volume_molecule_data(struct frame_data_list const * const fdlp,
                                                       FILE *master_header)
{
  /* Control flags */
  byte viz_mol_pos_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_POS);
  byte viz_mol_orient_flag = (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_ORIENT);
  byte viz_mol_states_flag = (viz_mol_pos_flag  &&  (world->viz_output_flag & VIZ_MOLECULES_STATES));

  /* Filenames */
  char *mol_pos_name = NULL;
  char *mol_orient_name = NULL;
  char *mol_states_name = NULL;

  /* Build filenames */
  if (viz_mol_pos_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_MOL_POS_NAME,
                                                              &mol_pos_name))
    goto failure;
  if (viz_mol_orient_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_MOL_ORIENT_NAME,
                                                                 &mol_orient_name))
    goto failure;
  if (viz_mol_states_flag  &&  dreamm_v3_grouped_create_filename(DREAMM_GROUPED_MOL_STATES_NAME,
                                                                 &mol_states_name))
    goto failure;

  /* Output molecule info */
  if (dreamm_v3_generic_dump_volume_molecule_data(fdlp,
                                                  master_header,
                                                  world->viz_state_info.filename_prefix_dirname,
                                                  mol_pos_name,
                                                  mol_orient_name,
                                                  mol_states_name,
                                                  &world->viz_state_info.dx_main_object_index))
    goto failure;

  if (mol_pos_name) free(mol_pos_name);
  if (mol_orient_name) free(mol_orient_name);
  if (mol_states_name) free(mol_states_name);
  return 0;

failure:
  if (mol_pos_name) free(mol_pos_name);
  if (mol_orient_name) free(mol_orient_name);
  if (mol_states_name) free(mol_states_name);
  return 1;
}

/*************************************************************************
dreamm_v3_grouped_dump_volume_molecules:
    Write the desired volume molecule data to the volume molecule data files
    and write appropriate index information into the master header.

        In:  struct frame_data_list *fdlp - frame for which to write mols
             FILE *master_header - file to which to write group
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_dump_volume_molecules(struct frame_data_list const * const fdlp,
                                                   FILE *master_header)
{
  int mol_index_base = world->viz_state_info.dx_main_object_index;
  if (dreamm_v3_grouped_dump_volume_molecule_data(fdlp,
                                                  master_header))
    return 1;

  if (world->viz_state_info.n_vol_species > 0)
  {
    int field_idx_base = world->viz_state_info.dx_main_object_index;
    world->viz_state_info.dx_main_object_index += world->viz_state_info.n_vol_species;

    /* Build fields for volume mols here */
    if (dreamm_v3_generic_write_molecule_fields(fdlp,
                                                master_header,
                                                world->viz_state_info.vol_species,
                                                world->viz_state_info.n_vol_species,
                                                field_idx_base,
                                                mol_index_base))
      return 1;

    /* Create group objects for volume molecules */
    world->viz_state_info.dreamm_last_iteration_vol_mols = fdlp->viz_iteration;
    if (dreamm_v3_grouped_write_molecule_group(master_header,
                                               fdlp,
                                               "volume",
                                               world->viz_state_info.vol_species,
                                               world->viz_state_info.n_vol_species,
                                               field_idx_base))
      return 1;

    /* Store iteration_number for volume_molecules */
    if (add_to_iteration_counter(&world->viz_state_info.vol_mol_output_iterations,
                                 fdlp->viz_iteration))
      return 1;
  }
  fprintf(master_header, "\n"); 
  return 0;
}

/*************************************************************************
dreamm_v3_grouped_write_combined_group:
    Write a "combined group" to the header file.  The combined group will
    contain the most recent version of both the mesh and the molecule data, as
    of the current iteration.  Also updates the combined_group_members data
    structure so that the combined group will be properly written into the
    frame_data series object.

        In: FILE *master_header - file to which to write group
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_write_combined_group(FILE *master_header,
                                                  long long viz_iteration)
{
  int combined_group_index = world->viz_state_info.dx_main_object_index ++;

  fprintf(master_header, "object %d group\n", combined_group_index);

  if (world->viz_state_info.dreamm_last_iteration_meshes != -1)
    fprintf(master_header,
            "\tmember \"meshes\" value \"meshes_%lld\"\n",
            world->viz_state_info.dreamm_last_iteration_meshes);

  if (world->viz_state_info.dreamm_last_iteration_vol_mols != -1)
    fprintf(master_header,
            "\tmember \"volume_molecules\" value \"volume_molecules_%lld\"\n",
            world->viz_state_info.dreamm_last_iteration_vol_mols); 

  if (world->viz_state_info.dreamm_last_iteration_surf_mols != -1)
    fprintf(master_header,
            "\tmember \"surface_molecules\" value \"surface_molecules_%lld\"\n",
            world->viz_state_info.dreamm_last_iteration_surf_mols);

  fprintf(master_header,"\n");

  /* create an entry into a 'frame_data' object. */
  char *str = alloc_sprintf("\tmember %d value %d position %lld\n",
                            world->viz_state_info.combined_group_members.n_strings,
                            combined_group_index,
                            viz_iteration);
  if (str == NULL)
  {
    fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  if (add_string_to_buffer(&world->viz_state_info.combined_group_members, str))
  {
    free(str);
    return 1;
  }

  fprintf(master_header, "\n\n");
  return 0;
}

/*************************************************************************
dreamm_v3_grouped_write_frame_series:
    Write the frame_data series object into the file.  The frame_data series
    object ties together all of the combined groups which have been created at
    the end of each iteration.

        In: FILE *master_header - file to which to write series
        Out: 0 on success, 1 on error
**************************************************************************/
static int dreamm_v3_grouped_write_frame_series(FILE *master_header)
{
  fprintf(master_header, "object \"frame_data\" class series\n");
  int frame_data_index;
  for (frame_data_index = 0;
       frame_data_index < world->viz_state_info.combined_group_members.n_strings;
       ++ frame_data_index)
    fprintf(master_header,
            "\t%s",
            world->viz_state_info.combined_group_members.strings[frame_data_index]);
  fprintf(master_header, "\n\n");
  return 0;
}

/*************************************************************************
output_dreamm_objects_grouped:
        In: struct frame_data_list *fdlp
        Out: 0 on success, 1 on error; output visualization files (*.dx)
             in dreamm  group format are written.
**************************************************************************/
int output_dreamm_objects_grouped(struct frame_data_list const * const fdlp)
{
  FILE *master_header = NULL;

  /* Flags */
  byte viz_mols = 0;
  byte viz_meshes = 0;

  no_printf("Viz output in DREAMM_V3_GROUPED mode...\n");

  /* Initialize flags */
  if (fdlp->type == ALL_MOL_DATA  ||  fdlp->type == MOL_POS  ||  fdlp->type == MOL_ORIENT)
    viz_mols = 1;
  if (fdlp->type == ALL_MESH_DATA  ||  fdlp->type == MESH_GEOMETRY  ||  fdlp->type == REG_DATA)
    viz_meshes = 1;

  /* Open master header file. */
  {
    char *master_header_file_path = dreamm_v3_grouped_get_master_header_name();
    if (master_header_file_path == NULL)
      return 1;

    if ((master_header = open_file(master_header_file_path, "a")) == NULL)
    {
      free(master_header_file_path);
      return 1;
    }
    free(master_header_file_path);
  }

  /* dump walls */
  if (viz_meshes)
    if (dreamm_v3_grouped_dump_meshes(fdlp, master_header))
      goto failure;

  if (viz_mols)
  {
    /* Dump grid molecules. */
    if (dreamm_v3_grouped_dump_grid_molecules(fdlp, master_header))
      goto failure;

    /* Dump 3D molecules: */
    if (dreamm_v3_grouped_dump_volume_molecules(fdlp, master_header))
      goto failure;
  }

  /* Create combined group if we're the last frame this iteration */
  if (! dreamm_v3_generic_scan_for_frame(fdlp->next, fdlp->viz_iteration))
  {
    if (dreamm_v3_grouped_write_combined_group(master_header, fdlp->viz_iteration))
      goto failure;

    /* Put value of viz_iteration into the time values */
    if (add_to_iteration_counter_monotonic(&world->viz_state_info.output_times,
                                           fdlp->viz_iteration))
      goto failure;

    /* If we're on the final iteration... */
    if (world->it_time == world->viz_state_info.final_iteration)
    {
      if (dreamm_v3_grouped_dump_time_info(master_header))
        goto failure;

      if (dreamm_v3_grouped_write_frame_series(master_header))
        goto failure;
    }
  }

  if (master_header != NULL) fclose(master_header);
  return 0;

failure:
  if (master_header != NULL) fclose(master_header);
  return 1;
}

/************************************************************************ 
output_rk_custom:
Rex Kerr's personal visualization mode output function 
*************************************************************************/

int output_rk_custom(struct frame_data_list *fdlp)
{
  FILE *custom_file;
  char cf_name[1024];
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  struct grid_molecule *gmp;
  short orient = 0;
  struct species *target;
  
  int i,j,k;
  double d;

  int id;
  struct vector3 where;
  
  no_printf("Output in CUSTOM_RK mode...\n");
  
  if (world->rk_mode_var==NULL) return output_ascii_molecules(fdlp);
  
  if ((fdlp->type==ALL_FRAME_DATA) || (fdlp->type==MOL_POS) || (fdlp->type==MOL_STATES))
  {
    sprintf(cf_name,"%s.rk.dat",world->molecule_prefix_name);
    custom_file = fopen(cf_name,(world->rk_mode_var->n_written)?"a+":"w");
    if (!custom_file)
    {
      fprintf(world->log_file,"Couldn't open file %s for viz output.\n",cf_name);
      return 1;
    }
    else { no_printf("Writing to file %s\n",cf_name); }
    
    world->rk_mode_var->n_written++;
    fprintf(custom_file,"%lld",fdlp->viz_iteration);
    for (j=0;j<world->n_species;j++)
    {
      target = world->species_list[j];
      if (target==NULL) continue;
      if (target->viz_state==EXCLUDE_OBJ) continue;
      for (k=0;k<world->rk_mode_var->n_bins;k++) world->rk_mode_var->bins[k] = 0;
      id = target->viz_state;
      
      fprintf(custom_file,"\t%d",id);
      
      for (slp = world->storage_head ; slp != NULL ; slp = slp->next)
      {
	for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale)
	{
	  for (i=-1;i<shp->buf_len;i++)
	  {
	    if (i<0) aep = shp->current;
	    else aep = shp->circ_buf_head[i];
	    
	    for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next)
	    {
	      amp = (struct abstract_molecule*)aep;
	      if (amp->properties != target) continue;
  
	      if ((amp->properties->flags & NOT_FREE)==0)
	      {
		mp = (struct volume_molecule*)amp;
		where.x = mp->pos.x;
		where.y = mp->pos.y;
		where.z = mp->pos.z;
	      }
	      else if ((amp->properties->flags & ON_GRID)!=0)
	      {
		gmp = (struct grid_molecule*)amp;
		uv2xyz(&(gmp->s_pos),gmp->grid->surface,&where);
		orient = gmp->orient;
	      }
	      else continue;
	      
	      d = dot_prod(&where , world->rk_mode_var->direction);
	      
	      k = bin(world->rk_mode_var->parts,world->rk_mode_var->n_bins-1,d);
	      world->rk_mode_var->bins[k]++;
	    }
	  }
	}
      }
      for (i=k=0 ; k<world->rk_mode_var->n_bins ; k++) i+=world->rk_mode_var->bins[k];
      if (i!=target->population) printf("Wanted to bin %d but found %d instead\n",target->population,i);
      for (k=0 ; k<world->rk_mode_var->n_bins ; k++) fprintf(custom_file," %d",world->rk_mode_var->bins[k]);
    }
    fprintf(custom_file,"\n");
    fclose(custom_file);
  }
  
  return 0;
}


/************************************************************************ 
output_ascii_molecules:
In: a frame data list (internal viz output data structure)
Out: 0 on success, 1 on failure.  The positions of molecules are output
     in exponential floating point notation (with 8 decimal places)
*************************************************************************/

int output_ascii_molecules(struct frame_data_list *fdlp)
{
  FILE *custom_file;
  char cf_name[1024];
  char cf_format[256];
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  struct grid_molecule *gmp;
  short orient = 0;
  
  int ndigits,i;
  long long lli;

  int id;
  struct vector3 where;
  
  no_printf("Output in ASCII mode (molecules only)...\n");
  
  if ((fdlp->type==ALL_FRAME_DATA) || (fdlp->type==MOL_POS) || (fdlp->type==MOL_STATES))
  {
    lli = 10;
    for (ndigits = 1 ; lli <= world->iterations && ndigits<20 ; lli*=10 , ndigits++) {}
    sprintf(cf_format,"%%s.ascii.%%0%dlld.dat",ndigits);
    sprintf(cf_name,cf_format,world->molecule_prefix_name,fdlp->viz_iteration);
    custom_file = fopen(cf_name,"w");
    if (!custom_file)
    {
      fprintf(world->log_file,"Couldn't open file %s for viz output.\n",cf_name);
      return 1;
    }
    else { no_printf("Writing to file %s\n",cf_name); }
    
    for (slp = world->storage_head ; slp != NULL ; slp = slp->next)
    {
      for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale)
      {
        for (i=-1;i<shp->buf_len;i++)
        {
          if (i<0) aep = shp->current;
          else aep = shp->circ_buf_head[i];
          
          for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next)
          {
            amp = (struct abstract_molecule*)aep;
            if (amp->properties == NULL) continue;
            if (amp->properties->viz_state == EXCLUDE_OBJ) continue;

            id = amp->properties->viz_state;
            
            if ((amp->properties->flags & NOT_FREE)==0)
            {
              mp = (struct volume_molecule*)amp;
              where.x = mp->pos.x;
              where.y = mp->pos.y;
              where.z = mp->pos.z;
            }
            else if ((amp->properties->flags & ON_GRID)!=0)
            {
              gmp = (struct grid_molecule*)amp;
              uv2xyz(&(gmp->s_pos),gmp->grid->surface,&where);
              orient = gmp->orient;
            }
            else continue;
            
            where.x *= world->length_unit;
            where.y *= world->length_unit;
            where.z *= world->length_unit;
            fprintf(custom_file,"%d %15.8e %15.8e %15.8e %2d\n",id,where.x,where.y,where.z,orient);
          }
        }
      }
    }
    fclose(custom_file);
  }
  
  return 0;
}

/********************************************************************* 
init_frame_data_list:
   In: struct frame_data_list* fdlp
   Out: 0 on success, 1 on error.
        Initializes frame_data_list structure. 
        Sets the value of the current iteration step to the start value.
        Sets the number of iterations. 
        Initializes state used in output_dreamm_objects and
                     output_dreamm_objects_grouped
***********************************************************************/
int init_frame_data_list(struct frame_data_list **fdlpp)
{
  int mol_orient_frame_present = 0;
  int mol_pos_frame_present = 0;
  int reg_data_frame_present = 0;
  int mesh_geometry_frame_present = 0;
  struct frame_data_list *fdlp;

  if (fdlpp == NULL) return 0;

  fdlp = *fdlpp;
  if (fdlp == NULL) return 0;

  switch (world->viz_mode)
  {
    case DREAMM_V3_MODE:
      if (dreamm_v3_generic_preprocess_frame_data(fdlpp))
        return 1;
      fdlp = *fdlpp;
      if (dreamm_v3_init(fdlp))
        return 1;
      break;

    case DREAMM_V3_GROUPED_MODE:
      if (dreamm_v3_generic_preprocess_frame_data(fdlpp))
        return 1;
      fdlp = *fdlpp;
      if (dreamm_v3_grouped_init(fdlp))
        return 1;
      break;

    default:
      count_time_values(fdlp);
      if (reset_time_values(fdlp, world->start_time))
        return 1;
      break;
  }

  for (; fdlp != NULL; fdlp = fdlp->next)
  {
    if (fdlp->curr_viz_iteration == NULL)
      continue;

    switch (fdlp->type)
    {
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

      case MESH_GEOMETRY:
        mesh_geometry_frame_present = 1;
        break;

      case REG_DATA:
        reg_data_frame_present = 1;
        break;

      case ALL_MESH_DATA:
        mesh_geometry_frame_present = 1;
        reg_data_frame_present = 1;
        break;

      default:
        /* Do nothing */
        ;
    }
  } /* end while */

  /* Check that the user hasn't selected a useless set of output info */
  if ((mol_orient_frame_present) & (!mol_pos_frame_present))
    fprintf(world->log_file, "The input file contains ORIENTATIONS but not POSITIONS statement in the MOLECULES block. The molecules cannot be visualized.\n");
  if ((reg_data_frame_present) & (!mesh_geometry_frame_present))
    fprintf(world->log_file, "The input file contains REGION_DATA but not GEOMETRY statement in the MESHES block. The meshes cannot be visualized.\n");

  return 0;
}

/**************************************************************************
update_frame_data_list:
        In: struct frame_data_list * fdlp
        Out: Nothing. Calls output visualization functions if necessary.
             Updates value of the current iteration step and pointer
             to the current iteration in the linked list
**************************************************************************/
void update_frame_data_list(struct frame_data_list *fdlp)
{
  if (fdlp == NULL) return;

  /* These statements need to precede handling of any frames to make sure
   * symlinks are created properly.
   */
  if (world->viz_mode == DREAMM_V3_MODE)
  {
    dreamm_v3_clean_files(fdlp);
    dreamm_v3_update_last_iteration_info(fdlp);
  }

  /* Scan over all frames, producing appropriate output. */
  for (; fdlp != NULL; fdlp = fdlp->next)
  {
    if (world->it_time!=fdlp->viz_iteration)
      continue;

    switch (world->viz_mode)
    {
      case DX_MODE:
        output_dx_objects(fdlp);
        break;

      case DREAMM_V3_MODE:
        output_dreamm_objects(fdlp);
        break;

      case DREAMM_V3_GROUPED_MODE:
        output_dreamm_objects_grouped(fdlp);
        break;

      case RK_MODE:
        output_rk_custom(fdlp);
        break;

      case ASCII_MODE:
        output_ascii_molecules(fdlp);
        break;

      case NO_VIZ_MODE:
      default:
        /* Do nothing for vizualization */
        break;
    }

    fdlp->curr_viz_iteration = fdlp->curr_viz_iteration->next;
    if (fdlp->curr_viz_iteration != NULL)
      fdlp->viz_iteration = frame_iteration(fdlp->curr_viz_iteration->value, fdlp->list_type);
  }
}

