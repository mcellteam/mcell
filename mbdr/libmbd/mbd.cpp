/*****************************************************************
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

/* C headers */
#include <bzlib.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <zlib.h>

/* STL headers */
#include <iostream>


/* local headers */
#include "mbd.h"


/** pull in some things from namespace std */
using std::vector;
using std::cout;
using std::endl;
using std::string;


/*****************************************************************
 * constructor
 *****************************************************************/
IO::IO(string afile_name) 
  :
    file_name_(afile_name),
    num_blocks_(0),
    output_type_(0),
    blocksize_(0),
    stepsize_(0.0),
    needed_api_tag_("MCELL_BINARY_API_2")
{}


/******************************************************************
 * destructor 
 ******************************************************************/
IO::~IO()
{
  free(buffer_);
}


/*****************************************************************
 *
 * public member functions
 *
 *****************************************************************/

/******************************************************************
 * member function initizializing the object, vie opening
 * the file and parsing the header 
 ******************************************************************/
bool 
IO::initialize()
{
  /* open file, don't parse headers if that fails already */
  bool status = open_();

  if (!status)
  {
    return false;
  }

  /* parse headers */
  status = parse_header_();

  return status;
}


/*******************************************************************
 * getter function returning the number of datablocks 
 *******************************************************************/
uint64_t 
IO::get_num_datablocks() const
{
  return num_blocks_;
}


/*******************************************************************
 * getter function returning the STEP
 *******************************************************************/
double 
IO::get_stepsize() const
{
  return stepsize_;
}


/******************************************************************
 * getter function returning the block size
 ******************************************************************/
uint64_t 
IO::get_blocksize() const
{
  return blocksize_;
}



/******************************************************************
 * getter function for number of columns contained in data at id
 *****************************************************************/
uint64_t 
IO::get_num_columns(int id) const
{
  if (id < 0 || id > num_blocks_-1)
    return -1;

  return block_info_[id].num_cols;
}



/******************************************************************
 * getter function returning the block size
 ******************************************************************/
vector<string> 
IO::get_blocknames() const
{
  return block_names_;
}



/*******************************************************************
 * getter function returning the type of output
 ******************************************************************/
uint16_t 
IO::get_output_type() const
{
  return output_type_;
}


/*******************************************************************
 * getter function to return the iteration list
 ******************************************************************/
std::vector<double>
IO::get_iteration_list() const
{
  return time_it_list_;
}


/*******************************************************************
 * returns the id corresponding to the data set of the given name
 *******************************************************************/
int 
IO::get_id_from_name(std::string data_name) const
{
  std::map<std::string, int>::const_iterator id = 
    block_names_map_.find(data_name);

  if (id == block_names_map_.end())
    return -1;

  return id->second;
}



/*******************************************************************
 * getter function for the raw data either selected by name or
 * by id
 * 
 * two types of data are returned via references:
 * 1) data_list: a vector of data rows containing doubles
 * 2) data_types: a vector of ints describing the underlying datatype 
 *                contained in each data_list row.
 *
 * NOTE: this method will not cast the data according to its true
 *       data type as returned by data_type. It's up to the caller
 *       to decide if this is necessary.
 ******************************************************************/
bool 
IO::get_data(std::vector<std::vector<double> >& data_list, 
  std::vector<uint16_t>& data_types, string name) const
{
  int id = get_id_from_name(name);
  if (id == -1)
    return false;

  /* call via block id */
  return get_data(data_list, data_types, id);
}


bool 
IO::get_data(std::vector<std::vector<double> >& data_list, 
  std::vector<uint16_t>& data_types, int id) const
{
  if (id < 0 || id > num_blocks_-1)
    return false;

  EntryType entry = block_info_[id];
  data_list.reserve(entry.num_cols);
  for (uint64_t i=0; i < entry.num_cols; ++i) 
  {
    std::vector<double> column;
    column.reserve(blocksize_ + 1);
    data_list.push_back(column);
    data_types.push_back(entry.data_type[i]);
  }

  /* fast forward to starting location of data */
  char* chunkLocation = buffer_ +
    fileoffset_to_data_ + sizeof(double) * buffer_size_ * entry.offset;

  uint64_t numElements = 0;
  uint64_t stream_block = 1;
  while (numElements < blocksize_)
  {
    /* check if we need to fast forward to next block 
      * NOTE: We need a special case for the last and potentiall incomplete
      *       stream block with less than buffer_size elements 
      */
    if (numElements >= (stream_block * buffer_size_))
    {
      chunkLocation = buffer_ + fileoffset_to_data_ + sizeof(double) 
        * (buffer_size_ * stream_block) * total_columns_;
      if ((blocksize_ - numElements) >= buffer_size_)
      {
        chunkLocation += sizeof(double) * buffer_size_ * entry.offset;
      } 
      else 
      {
        chunkLocation += sizeof(double) * (blocksize_ - numElements) 
          * entry.offset;
      }
      ++stream_block;
    }

    for (uint64_t i=0; i < entry.num_cols; ++i)
    {
      double* temp = reinterpret_cast<double*>(chunkLocation);
      data_list[i].push_back(*temp);
      chunkLocation += sizeof(double);
    }
    ++numElements;
  }
  return true;
}



/*******************************************************************
 * member function printing basic file info as extracted from
 * the headers
 *******************************************************************/
void 
IO::info() const
{
  cout << "MCell3 binary reader v" << VERSION << " (C) " 
       << DATE << " " << AUTHOR << NL
       << "----------------------------------------------------" << NL
       << PRE << "found " << num_blocks_ << " datablocks with " << blocksize_
       << " items each." << NL;
  
  if (output_type_ == STEP) 
  {
    cout << PRE << "count STEP was " << stepsize_ << "[s]" << NL;
  }
  else if (output_type_ == TIME_LIST)
  {
    cout << PRE << "Output was generated via TIME_LIST" << NL;
  }
  else if (output_type_ == ITERATION_LIST)
  {
    cout << PRE << "Output was generated via ITERATION_LIST" << NL;
  }
  else
  {
          cout << PRE << "Unknown output type" << NL;
  }

  cout << endl;
}


/******************************************************************
 * 
 * private member functions
 *
 ******************************************************************/

/*******************************************************************
 * member function opening the binary file 
 *******************************************************************/
bool
IO::open_()
{
  FILE* input_file = fopen(file_name_.c_str(),"r");
  if (input_file == NULL)
  {
    return false;
  }

  /* call the proper parsing code depending on if the file is uncompressed 
   * data of a compressed bzip file */
  string::size_type location = file_name_.find_last_of(".");
  string extension = file_name_.substr(location);

  bool status = false;
  if (extension == ".gz") {
    status = open_gz_file_(input_file);
  }
  else if (extension == ".bz2")
  {
    status = open_bz_file_(input_file);
  }
  // try it as uncompressed binary file
  else
  {
    status = open_binary_file_(input_file);
  }

  /* done with file */
  fclose(input_file);

  return status;
}


/*****************************************************************
 * member function parsing the header 
 *****************************************************************/
bool 
IO::parse_header_()
{
  /* keep track of where we are */
  uint64_t charCounter = 0;

  charCounter = parse_api_tag_(charCounter);
  if (charCounter == 0)
          return false;

  charCounter = parse_block_information_(charCounter);
  charCounter = parse_block_names_(charCounter);

  // store the total number of columns and offset to begin of data
  fileoffset_to_data_ = charCounter;

  return true;
}



/*********************************************************************
 * parse the api tag
 ********************************************************************/
uint64_t
IO::parse_api_tag_(uint64_t charCounter)
{
  /** read API tag and make sure we have the proper output file type */
  string receivedAPItag;
  for (unsigned int i=0; i < needed_api_tag_.size(); ++i)
  {
    receivedAPItag.push_back(*(buffer_ + charCounter));
    ++charCounter;
  }
  ++charCounter;

  if (receivedAPItag != needed_api_tag_) 
  {
    return 0;
  }

  return charCounter;
}


/*********************************************************************
 * parse the block information 
 ********************************************************************/
uint64_t 
IO::parse_block_information_(uint64_t charCounter)
{
  output_type_ = *reinterpret_cast<uint16_t*>(buffer_ + charCounter);
  charCounter += sizeof(uint16_t);

  blocksize_ = *reinterpret_cast<uint64_t*>(buffer_ + charCounter);
  charCounter += sizeof(uint64_t);

  uint64_t length = 0;
  if (output_type_ == STEP)          // STEP
  {
    length = *reinterpret_cast<uint64_t*>(buffer_ + charCounter);
    charCounter += sizeof(uint64_t);

    stepsize_ = *reinterpret_cast<double*>(buffer_ + charCounter);
    charCounter += sizeof(double); 
  }
  else if (output_type_ == TIME_LIST || output_type_ == ITERATION_LIST) 
  {
    length = *reinterpret_cast<uint64_t*>(buffer_ + charCounter);
    charCounter += sizeof(uint64_t);

    for (uint64_t i = 0; i < length; ++i)
    {
      double temp;
      temp = *reinterpret_cast<double*>(buffer_ + charCounter);
      charCounter += sizeof(double); 

      time_it_list_.push_back(temp);
    }
  }
  else 
  {
    std::cout << "Unknown count output format" << std::endl;
    return false;
  }

  /* parse output buffer size */
  buffer_size_ = *reinterpret_cast<uint64_t*>(buffer_ + charCounter);
  charCounter += sizeof(uint64_t);

  /* parse number of blocks */
  num_blocks_ = *reinterpret_cast<uint64_t*>(buffer_ + charCounter);
  charCounter += sizeof(uint64_t);

  return charCounter;
}



/*********************************************************************
 * parse the block names
 ********************************************************************/
uint64_t 
IO::parse_block_names_(uint64_t charCounter)
{
  uint64_t block_counter = 0;
  uint64_t offset = 0;
  string buffer;
  while ((block_counter != num_blocks_))
  {
    char currentChar = *(buffer_ + charCounter);

    /** read until we hit the end of a substring; this is then copied to a
     ** new string and appended to available_blocks */
    if (currentChar == '\0')
    {
      EntryType entry;

      block_names_.push_back(buffer);
      block_names_map_[buffer] = block_counter;
      entry.name = buffer;
      buffer.clear();

      ++charCounter;

      entry.offset = offset;
      entry.num_cols = *reinterpret_cast<uint64_t*>(buffer_ + charCounter);
      offset += entry.num_cols;
      charCounter += sizeof(uint64_t);

      for(uint64_t i=0; i<entry.num_cols; ++i)
      {
        entry.data_type.push_back(*reinterpret_cast<uint16_t*> (buffer_ + charCounter));
        charCounter += sizeof(uint16_t);
      }

      block_info_.push_back(entry);

      if (++block_counter > num_blocks_) 
      {
        return false;
      }
    }
    else 
    {
      buffer.push_back(currentChar);
      ++charCounter;
    }
  }

  total_columns_ = offset;

  return charCounter;
}



/**********************************************************************
 * uncompress and read already opened binary file into a buffer_
 *********************************************************************/
bool 
IO::open_gz_file_(FILE* inputFile)
{
  // get the file descriptor of the FILE object
  int fdKeep = fileno(inputFile);
  int fd = dup(fdKeep);
  if (fd == -1) {
    return false;
  }

  gzFile myFile = gzdopen OF((fd, "r")); 
  if ( myFile == NULL ) {
    return false;
  }

  /* keep reading data in chunks of 1M */
  const long chunkSize = 1048576L;
  int chunkCount = 1;
  int status = 0;
  char* nextChunkLocation = NULL;
  buffer_ = NULL;
  
  do
  {
    /* request next chunk */
    buffer_ = (char*)realloc(buffer_, chunkSize*chunkCount/sizeof(char) );

    if ( buffer_ == NULL )
    {
      return false;
    }

    /* read it */
    nextChunkLocation = buffer_ + chunkSize * (chunkCount-1);

    status = gzread OF((myFile, nextChunkLocation, chunkSize));
    if ( status == -1 ) {
      return false;
    }

    ++chunkCount;
  } while ( status != 0 );

  gzclose(myFile);  

  return true;
}


/******************************************************************
 * uncompress and read already opened binary file into a Buffer_
 ******************************************************************/
bool
IO::open_bz_file_(FILE* inputFile)
{
  int bzError;
  BZFILE* myFile;
  myFile = BZ2_bzReadOpen(&bzError, inputFile, 0, 0, NULL, 0);

  if ( bzError != BZ_OK )
  {
    return false;
  }

  /* keep reading data in chunks of 1M */
  const long chunkSize = 1048576L;
  int chunkCount = 1;
  char* nextChunkLocation = NULL;
  buffer_ = NULL;

  do
  {
    /* request next chunk */
    buffer_ = (char*)realloc(buffer_, chunkSize*chunkCount/sizeof(char) );

    if ( buffer_ == NULL )
    {
      return false;
    }

    /* read it */
    nextChunkLocation = buffer_ + chunkSize * (chunkCount-1);
    BZ2_bzRead(&bzError, myFile, nextChunkLocation, chunkSize);
 
    if ( ( bzError != BZ_OK ) && ( bzError != BZ_STREAM_END ) )
    {
      return false;
    }

    ++chunkCount;
  } while ( bzError != BZ_STREAM_END );

  /* release libbzip2's memory */
  BZ2_bzReadClose(&bzError, myFile);

  if ( bzError != BZ_OK )
  {
    return false;
  }

  return true;
}


/******************************************************************
 * read already opened binary file into a buffer_
 ******************************************************************/
bool IO::open_binary_file_(FILE* inputFile)
{
  /* keep reading data in chunks of 1M */
  const unsigned int chunkSize = 1048576;
  int chunkCount = 1;
  size_t num_read = 0;
  buffer_ = NULL;
  
  do
  {
    /* request next chunk */
    buffer_ = (char*)realloc(buffer_, chunkSize*chunkCount*sizeof(char) );

    if ( buffer_ == NULL )
    {
      return false;
    }

    /* read it */
    char* nextChunkLocation = buffer_ + chunkSize * (chunkCount-1) * 
    sizeof(char);
    num_read = fread(nextChunkLocation, sizeof(char), chunkSize, inputFile);
    ++chunkCount;

  } while ( num_read == chunkSize );

  if ((feof(inputFile) != 0) && (ferror(inputFile) != 0 ) )
  {
    return false;
  }

  return true;
}



/*************************************************************************
 * 
 * implementation of C interface 
 *
 ************************************************************************/

/*
 * construct a new MBDIO data parsing object for the binary reaction
 * data file file_name
 */
MBDIO
mbd_new(const char* file_name)
{
  return reinterpret_cast<MBDIO>(new IO(file_name));
}


/*
 * delete an existing MBDIO data parsing object 
 */
void
mbd_delete(MBDIO obj)
{
  delete obj;
  return;
}


/*
 * initialize the file reader object
 */
int 
mbd_initialize(MBDIO obj) 
{
  return (obj->initialize() ? 1 : -1);
}


/*
 * returns the number of distinct data items (data_blocks)
 * in the data set handled by the data reader
 */
uint64_t 
mbd_get_num_datablocks(MBDIO obj)
{
  return obj->get_num_datablocks();
}


/*
 * returns the blocksize, i.e. the number of data items per
 * data set.
 */
uint64_t
mbd_get_blocksize(MBDIO obj)
{
  return obj->get_blocksize();
}


/*
 * returns the type of output (STEP, TIME_LIST, ITERATION_LIST)
 */
uint16_t 
mbd_get_output_type(MBDIO obj)
{
  return obj->get_output_type();
}
  

/*
 * returns the step size for output type STEP
 */
double
mbd_get_stepsize(MBDIO obj)
{
  return obj->get_stepsize();
}


/*
 * returns the names of each datablock and their count
 * NOTE: the caller is responsible for providing an array block_names
 *       of the proper size.
 */
void 
mbd_get_blocknames(MBDIO obj, char** block_names, size_t length)
{
  vector<string> names = obj->get_blocknames();
  for (int i = 0; i < names.size(); ++i)
  {
    strncpy(block_names[i], names[i].c_str(), length);
  }

  return;
}



/*
 * returns an array of timepoints/iterations and the count
 * NOTE: the caller is responsible for providing an array iteration_list
 *       of proper size (blocksize).
 */
void 
mbd_get_iteration_list(MBDIO obj, double* iteration_list) 
{
  vector<double> iter_list = obj->get_iteration_list();

  for (int i = 0; i < iter_list.size(); ++i)
  {
    iteration_list[i] = iter_list[i];
  }

  return;
}



/* 
 * return the number of available data columns in data set
 * with id.
 */
uint64_t
mbd_get_num_columns_by_id(MBDIO obj, int id)
{
  return obj->get_num_columns(id);
}



/*
 * return the data values corresponding to the data set with the
 * given id. The id is identical to the index of a given data set
 * in the list of blocknames.
 * NOTE: The called is responsible for providing data_list and
 *       data_types arrays of the correct size.
 */
int
mbd_get_data_by_id(MBDIO obj, double** data_list, uint16_t* data_types, 
                       int id)
{
  vector<vector<double> > raw_data;
  vector<uint16_t> raw_data_types;
  if (!obj->get_data(raw_data, raw_data_types, id))
    return -1;
 
  size_t col_count = raw_data_types.size();
  size_t row_count = obj->get_blocksize();

  /* extract data */
  for (int col = 0; col < col_count; ++col)
  {
    for (int row = 0; row < row_count; ++row)
    {
      data_list[col][row] = raw_data[col][row];
    }
  }

  /* extract data types */
  for (int col = 0; col < col_count; ++col) 
  {
    data_types[col] = raw_data_types[col];
  }

  return 0;
}
 


/*
 * return the id corresponding to the dataset with the given
 * name. If the name doesn't match any dataset -1 is returned.
 */
int 
mbd_get_id_from_name(MBDIO obj, char* name)
{
  return obj->get_id_from_name(string(name));
}



/*
 * return the data values corresponding to the data set with the
 * given name. The id is identical to the index of a given data set
 * in the list of blocknames.
 * NOTE: The called is responsible for providing data_list and
 *       data_types arrays of the correct size.
 */
int
mbd_get_data_by_name(MBDIO obj, double** data, uint16_t* data_types, 
                     char* data_name)
{
  int id = obj->get_id_from_name(string(data_name));
  if (id == -1)
    return 1;

  return mbd_get_data_by_id(obj, data, data_types, id);
}


/*
 * prints basic information regarding the data contained in the
 * reader object.
 */
void 
mbd_info(MBDIO obj)
{
  return obj->info();
}
