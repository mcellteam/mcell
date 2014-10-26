#ifndef BINARY_REACT_OUTPUT_H
#define BINARY_REACT_OUTPUT_H

#include "mcell_structs.h"

/* macro used for unifying regular and compressed writing of
 * binary reaction data */
#define BINARY_WRITE(data, size, block_data) \
  if (block_data->binary_out->compression_level == 0) { \
    fwrite(data, size, 1, block_data->binary_out->output_file); \
  } else { \
    if (block_data->binary_out->compression_type == COMPRESS_GZIP) { \
      int status = gzwrite OF((block_data->binary_out->gz_compressed_output_file, data, size)); \
      if (status == 0) { \
        mcell_error("Failed to write compressed binary output file %s", \
                    block_data->binary_out->filename);                       \
      } \
    } \
    else if (block_data->binary_out->compression_type == COMPRESS_BZIP2) { \
      int error = 0; \
      BZ2_bzWrite(&error, block_data->binary_out->bz_compressed_output_file, data, size); \
      if (error != BZ_OK) \
        mcell_error("Failed to write compressed binary output file %s", \
                    block_data->binary_out->filename);                       \
    }\
  } \


/* function initializing the binary reaction data output file */
int init_binary_reaction_data(struct output_block *block_data,
                               struct volume* world);
int create_binary_output_file(struct output_block *block_data);
int write_binary_header(struct output_block *block_data,
                        double time_step, long long iterations,
                        long long chkpt_iterations,
                        u_int checkpt_seq_number,
                        long long start_time);
int write_binary_data_info(struct output_block *block_data);

// functions for writing binary reaction data output
int write_binary_reaction_output(struct output_block* block_data,
                                 struct output_set *set);


// function for closing binary output file and doing
// any remaining cleanup
void close_binary_reaction_data(struct output_block *block_data);



#endif
