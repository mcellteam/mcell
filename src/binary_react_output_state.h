#ifndef BINARY_REACT_OUTPUT_STATE_H
#define BINARY_REACT_OUTPUT_STATE_H

/* required for binary data reader */
#include <zlib.h>
#include <bzlib.h>


/* some helper macros */
#define COMPRESS_GZIP 1
#define COMPRESS_BZIP2 2


/* 
 * struct keeping track of mcell simulation state required for 
 * binary reaction data output 
 */
struct binary_output
{
  /* number of reaction data items present in output block */
  int num_react_data;

  /* array holding the type of stored data (int or double) in each 
   * data block */
  char *type_array;

  /* output filename and directory */
  char *filename;
  char *directory;

  /* file pointer to output file for uncompressed output */
  FILE *output_file;

  /* file pointer to gzip compressed file stream and compression
   * level */
  gzFile gz_compressed_output_file;
  BZFILE* bz_compressed_output_file;
  int compression_type;
  int compression_level;
};

#endif
