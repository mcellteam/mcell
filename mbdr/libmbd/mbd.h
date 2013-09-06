/***************************************************************
 *
 * the class IO provides all the functionality to access a 
 * binary mcell3 reaction data output file.
 *
 * (c) 2008 Markus Dittrich 
 *
 * This program is free software; you can redistribute it 
 * and/or modify it under the terms of the GNU General Public 
 * License Version 3 as published by the Free Software Foundation. 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License Version 3 for more details.
 *
 * You should have received a copy of the GNU General Public 
 * License along with this program; if not, write to the Free 
 * Software Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *
 ****************************************************************/

#ifndef IO_H
#define IO_H

/* C header */
#include <stdint.h>

/* output types */
const int STEP = 1;
const int TIME_LIST = 2;
const int ITERATION_LIST = 3;

/* data types */
const int INT_DATA = 0;
const int DOUBLE_DATA = 1;

/*
 * C++ interface
 */
#ifdef __cplusplus

/* STL headers */
#include <vector>
#include <string>
#include <map>


/** global macros */
const std::string NL("\n");
const std::string PRE("reader> ");
const std::string VERSION("0.2");
const std::string AUTHOR("Markus Dittrich");
const std::string DATE("2013");


/**********************************************************
 * struct entryType keeps track of the type (int/double)
 * in each data block as well as the start and end of
 * the data in the binary file
 **********************************************************/
struct EntryType
{
  std::string name;
  uint64_t num_cols;
  uint64_t offset;
  std::vector<uint16_t> data_type;
};


/**********************************************************
 * class IO  
 *********************************************************/
struct IO
{

public:

  explicit IO(std::string file_name);
  ~IO();

  /* initialize data (parse headers, etc.) */
  bool initialize();

  /* getter functions for data properties */
  uint64_t get_num_datablocks() const;
  double get_stepsize() const;
  uint64_t get_blocksize() const;
  uint64_t get_num_columns(int id) const;
  uint16_t get_output_type() const;
  std::vector<std::string> get_blocknames() const;
  std::vector<double> get_iteration_list() const;
  int get_id_from_name(std::string data_name) const;

  bool get_data(std::vector<std::vector<double> >& data_list,
                std::vector<uint16_t>& data_types,
                int id) const;

  bool get_data(std::vector<std::vector<double> >& data_list,
                std::vector<uint16_t>& data_types,
                std::string name) const;

  /** print basic file info extracted from header */
  void info() const;


private:

  /** member functions */
  bool open_gz_file_(FILE* inputFile);
  bool open_bz_file_(FILE* inputFile);
  bool open_binary_file_(FILE* inputFile);
  bool open_();
  bool parse_header_();

  uint64_t parse_api_tag_(uint64_t charCounter);
  uint64_t parse_block_information_(uint64_t charCounter);
  uint64_t parse_block_names_(uint64_t charCounter);
 
  /** member data */
  char* buffer_;     // main data buffer 
  std::string file_name_;
  uint64_t num_blocks_;
  uint16_t output_type_;
  uint64_t blocksize_;
  
  double stepsize_;

  uint64_t buffer_size_;
  uint64_t total_columns_;
  uint64_t fileoffset_to_data_;
  
  std::vector<double> time_it_list_;
  std::vector<std::string> block_names_;
  std::map<std::string,int> block_names_map_;
  std::vector<EntryType> block_info_;
	
  std::string needed_api_tag_;
};
#else
struct IO;
#endif


/*
 * C interface
 */

typedef struct IO* MBDIO;

#ifdef __cplusplus
extern "C" 
{
#endif

MBDIO mbd_new(const char* file_name);
void mbd_delete(MBDIO obj);
int mbd_initialize(MBDIO obj);

/* getter functions for data properties */
uint64_t mbd_get_num_datablocks(MBDIO obj);
double mbd_get_stepsize(MBDIO obj);
uint64_t mbd_get_blocksize(MBDIO obj);
uint16_t mbd_get_output_type(MBDIO obj);
int mbd_get_id_from_name(MBDIO obj, char* name);


/* getter functions for metadata */
void mbd_get_blocknames(MBDIO obj, char** block_names, size_t length);
void mbd_get_iteration_list(MBDIO obj, double* iteration_list); 
uint64_t mbd_get_num_columns_by_id(MBDIO obj, int id);


/* getter functions for data */
int mbd_get_data_by_id(MBDIO obj, double** data_list, uint16_t* data_types,
                       int id);

int mbd_get_data_by_name(MBDIO obj, double** data_list, uint16_t* data_types,
                         char* name);

/** print basic file info extracted from header */
void mbd_info(MBDIO obj);

#ifdef __cplusplus
}
#endif

#endif
