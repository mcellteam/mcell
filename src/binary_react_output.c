
/* C headers */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* local includes */
#include "binary_react_output.h"
#include "logging.h"


/********************************************************************
 *
 * this function initializes the binary reaction data output,
 * opens the binary data files, and writes the initial file
 * headers.
 *
 ********************************************************************/
int init_binary_reaction_data(struct output_block *block_data,
                              struct volume *world)
{
  /* make sure the compression level is sane */
  if (block_data->binary_out->compression_level < 0
      || block_data->binary_out->compression_level > 9)
  {
    mcell_log("The compression level for binary output files has "
              "to be between 0 (off) and 9 (maximal)");
    return 1;
  }


  if (create_binary_output_file(block_data)) return 1;
  if (write_binary_header(block_data, world->time_unit, world->iterations,
    world->chkpt_iterations))
    return 1;
  if (write_binary_data_info(block_data)) return 1;

  return 0;
}



/********************************************************************
 *
 * function creating the binary reaction data output files
 *
 ********************************************************************/
int
create_binary_output_file(struct output_block *block_data)
{
  /* create directory if requested */
  char dir_name[BUFSIZ] = "";
  if (block_data->binary_out->directory != NULL)
  {
    if (mkdirs(block_data->binary_out->directory) != 0)
      return 1;

    strncpy(dir_name, block_data->binary_out->directory, BUFSIZ-2);
    strncat(dir_name, "/", 2);
  }

  /* create output file */
  if (block_data->binary_out->compression_level == 0)
  {
    FILE *output_file = fopen(block_data->binary_out->filename, "w");
    if (output_file == NULL) return 1;
    block_data->binary_out->output_file = output_file;
  }
  else
  {
    /* initialize compressed gzip format */
    if (block_data->binary_out->compression_type == COMPRESS_GZIP)
    {
      gzFile compressed_output_file =
        gzopen OF((block_data->binary_out->filename, "w"));
      if (compressed_output_file == NULL) return 1;

      block_data->binary_out->gz_compressed_output_file = compressed_output_file;
      gzsetparams OF((compressed_output_file,
                      block_data->binary_out->compression_level,
                      Z_DEFAULT_STRATEGY));
    }
    /* initialize compressed bzip format */
    else if (block_data->binary_out->compression_type == COMPRESS_BZIP2)
    {
      FILE *output_file = fopen(block_data->binary_out->filename, "wb");
      if (output_file == NULL) return 1;

      block_data->binary_out->output_file = output_file;

      int error = 0;
      BZFILE* compressed_output_file =
        BZ2_bzWriteOpen(&error, block_data->binary_out->output_file,
                        block_data->binary_out->compression_level, 0, 0);
      if (error != BZ_OK) return 1;

      block_data->binary_out->bz_compressed_output_file = compressed_output_file;
    }
    else
    {
      mcell_log("Unknown compression method for binary reaction output.");
      return 1;
    }
  }

  return 0;
}


/***********************************************************************
 *
 * create the header file for the binary data file
 *
 * the header has the following structure:
 *
 * const char* api_tag
 * uint16_t output_type        1 = STEP;  2 = TIME_LIST/ITERATION_LIST
 * uint16_t time_list_length   length of timelist/iteration list
 * double [] time_list         list of output times/iterations; for
 *                             output type STEP only contains a single
 *                             double equal to the output time step
 * uint32_t buffersize         size of the output buffer used to write
 *                             the output (needed for parsing the output)
 *
 **********************************************************************/
int write_binary_header(struct output_block *block_data,
                        double time_step, long long iterations,
                        long long chkpt_iterations)
{
  char api_tag[] = "MCELL_BINARY_API_2";
  BINARY_WRITE(api_tag, sizeof(api_tag), block_data);

  /* write type of output (STEP, TIMELIST, ITERATIONS) */
  uint64_t num_data_items = 0;
  uint64_t time_list_length = 0;
  uint16_t output_type = 0;
  if (block_data->timer_type == OUTPUT_BY_STEP)
  {
    output_type = 1;
    time_list_length = 1;

    double output_step = block_data->step_time;

    // the number of data items is determined by the smaller of iterations
    // or chkpt_iterations
    long long is = (chkpt_iterations < iterations) ? chkpt_iterations : iterations;
    num_data_items = (uint64_t)(is*time_step/output_step)+1;

    /* write info */
    BINARY_WRITE(&output_type, sizeof(output_type), block_data);
    BINARY_WRITE(&num_data_items, sizeof(num_data_items), block_data);
    BINARY_WRITE(&time_list_length, sizeof(time_list_length), block_data);
    BINARY_WRITE(&output_step, sizeof(output_step), block_data);
  }
  else if (block_data->timer_type == OUTPUT_BY_TIME_LIST
        || block_data->timer_type == OUTPUT_BY_ITERATION_LIST)
  {
    if (block_data->timer_type == OUTPUT_BY_TIME_LIST)
    {
      output_type = 2;
    }
    else
    {
      output_type = 3;
    }

    time_list_length = 0;
    struct num_expr_list *time_list = block_data->time_list_head;
    while (time_list != NULL)
    {
      ++time_list_length;
      time_list = time_list->next;
    }

    num_data_items = time_list_length;

    /* write info */
    BINARY_WRITE(&output_type, sizeof(output_type), block_data);
    BINARY_WRITE(&num_data_items, sizeof(num_data_items), block_data);
    BINARY_WRITE(&time_list_length, sizeof(time_list_length), block_data);
    time_list = block_data->time_list_head;
    while (time_list != NULL)
    {
      BINARY_WRITE(&(time_list->value), sizeof(double), block_data);
      time_list = time_list->next;
    }
  }

  /* write the output block size; the reader will need this info
   * to recover the data from the file */
  uint64_t buffersize = (uint64_t)(block_data->buffersize);
  BINARY_WRITE(&buffersize, sizeof(uint64_t), block_data);

  return 0;
}



/***********************************************************************
 *
 * create part of the header for the binary data file containing
 * information about the contained data sets
 *
 * this part of the header has the following structure:
 *
 * uint64_t num_data_files     number of data files contained in file
 *
 * --- repeat num_data_files times ----
 * char *  filename            name of the contained file
 * uint16_t num_columns        number of data columns in file
 *     ---- repeat num_columns times -----
 *     uint16_t column_type    data type contained in column
 *                             0 = int, 1 = double
 *
 **********************************************************************/
int write_binary_data_info(struct output_block *block_data)
{
  /* determine the total number of data files and write to file*/
  uint64_t num_data_files = 0;
  for (struct output_set *set = block_data->data_set_head; set != NULL;
       set=set->next)
    ++num_data_files;
  BINARY_WRITE(&num_data_files, sizeof(num_data_files), block_data);

  /* for each output file we extract the filename, the number of
   * columns and the type stored in each column (int/double) */
  for (struct output_set *set = block_data->data_set_head; set != NULL;
       set=set->next)
  {
    /* extract filename from path and write to file */
    char *name = set->outfile_name;
    char *filename_start = strrchr(name,'/');
    if (filename_start == NULL)
      filename_start = name;

    size_t filename_length = strlen(name)-(name-filename_start);
    char filename[filename_length];
    strncpy(filename, filename_start+1, filename_length);
    BINARY_WRITE(filename, strlen(filename)+1, block_data);

    /* extract the number of columns and their datatype */
    uint64_t num_columns = 0;
    struct output_column *col_ptr;
    for (col_ptr = set->column_head; col_ptr != NULL; col_ptr=col_ptr->next)
      ++num_columns;
    BINARY_WRITE(&num_columns, sizeof(num_columns), block_data);

    for (col_ptr = set->column_head; col_ptr != NULL; col_ptr=col_ptr->next)
    {
      uint16_t data_type = 0;
      if (col_ptr->data_type == COUNT_INT)
      {
        data_type = 0;
      }
      else if (col_ptr->data_type == COUNT_DBL)
      {
        data_type = 1;
      }
      else
      {
        mcell_log("Encountered unsupported output type (such as TRIGGER)\n"
                  "in binary reaction output. Please use regular reaction\n"
                  "data output instead.");
        return 1;
      }

      BINARY_WRITE(&data_type, sizeof(data_type), block_data);
    }
  }
  return 0;
}



/**************************************************************************
 * write_binary_reaction_output:
 *
 * this function writes reaction output data to the binary data file.
 *
 * In: the output_set we want to write to disk
 *     the flag that signals an end to the scheduled reaction outputs
 * Out: 0 on success, 1 on failure.
 *      The reaction output buffer is flushed and written to disk.
 *      Indices are not reset; that's the job of the calling function.
**************************************************************************/
int write_binary_reaction_output(struct output_block *block_data,
    struct output_set *set)
{
  struct output_column *column;
  u_int n_output;
  u_int i;

  n_output = set->block->buffersize;
  if (set->block->buf_index < set->block->buffersize)
    n_output = set->block->buf_index;

  /* Write data */
  double value;
  for (i=0;i<n_output;i++)
  {
    for (column=set->column_head; column!=NULL; column=column->next)
    {
      switch (column->data_type)
      {
        case COUNT_INT:
          value = (double)((int*)column->buffer)[i];
          BINARY_WRITE(&value, sizeof(double), block_data);
          break;

        case COUNT_DBL:
          BINARY_WRITE(&(((double*)column->buffer)[i]), sizeof(double),
                       block_data);
          break;

        default:
          mcell_error("Error in write_binary_reaction_output: "
                      "unknown reaction data type encountered.");
          break;
      }
    }
  }

  set->chunk_count++;
  return 0;
}



/*********************************************************************
 *
 * function for closing the binary reaction data ouput file and
 * doing any other necessary cleanup.
 *
 ********************************************************************/
void close_binary_reaction_data(struct output_block *block_data)
{
  free(block_data->binary_out->filename);
  if (block_data->binary_out->compression_level == 0)
  {
    fclose(block_data->binary_out->output_file);
  }
  else
  {
    if (block_data->binary_out->compression_type == COMPRESS_GZIP)
    {
      gzclose OF((block_data->binary_out->gz_compressed_output_file));
    }
    else
    {
      unsigned int bytes_in = 0;
      unsigned int bytes_out = 0;
      int error = 0;
      BZ2_bzWriteClose(&error, block_data->binary_out->bz_compressed_output_file,
                        0, &bytes_in, &bytes_out);
      if (error != BZ_OK)
        mcell_error("Error closing the bzlib compressed file.");

      /* also close the underlying file pointer */
      fclose(block_data->binary_out->output_file);
    }
  }
}
