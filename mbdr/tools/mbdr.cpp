/**************************************************************
 *
 * mbdr is a wrapper to extract information from binary mcell data files
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
 **************************************************************/

/** standerd C headers */
#include <cstdio>
#include <cstdlib>

/** STL headers */
#include <iostream>

/** boost headers */
#include <boost/foreach.hpp>
#include <boost/regex.hpp>

/** application headers */
#include "mbd.h"
#include "mbdr.hpp"

/* some convenience stuff */
#define foreach BOOST_FOREACH

/* pull some symbols into namespace */
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;


//----------------------------------------------------------------
//  main
//----------------------------------------------------------------
int main(int argc, char** argv)
{
  /** tokens for command line switches */
  bool listNames = false;
  bool extractData = false;
  bool extractDataByID = false;
  bool extractDataByName = false;
  bool extractDataByRegex = false;
  bool listInfo = false;
  bool write_file = false;
  bool print_times = false;
  std::string selectionString;

  /** read command line */
  int c;

  if ( argc == 1 ) usage();

  while ( (c = getopt(argc, argv,"lhN:I:R:iwt" )) != -1 )
  {
    switch (c)
    {
      case 'l':
	listNames = true;
	break;
      case 'h':
	usage();
	break;
      case 'I':
	extractData = true;
        extractDataByID = true;
	selectionString = optarg;
	break;
      case 'N':
        extractData = true;
        extractDataByName = true;
        selectionString = optarg;
        break;
      case 'R':
        extractData = true;
        extractDataByRegex = true;
        selectionString = optarg;
        break;
      case 'i':
	listInfo = true;
	break;
      case 'w':
        write_file = true;
        break;
      case 't':
        print_times = true;
        break;
      default:
	break;
    }
  }

  /* selecting by name, ID, or regex are mutually exclusive */
  if (  (extractDataByID && extractDataByName) 
     || (extractDataByID && extractDataByRegex) 
     || (extractDataByName && extractDataByRegex) )
  {
    usage();
  }

  /* check if at least a single filename was provided on the command line */
  if ( argc <= optind ) 
  {
    usage();
  }

  /* loop over all files privded on the command line and parse */
  while ( optind < argc )
  {

    /* generate IO object for file access and parse headers */
    std::string filename(*(argv+optind));
    IO file_object(filename);
 
    if ( !file_object.initialize() )
    {
      std::cout << "ERROR: Could not initialize " << filename << std::endl;
      exit(EXIT_FAILURE);
    }
  
    /* list info */
    if ( listInfo )
    {
      display_file_info(file_object);

      /* next file */
      ++optind;
      continue;
    }
 
    /* list block names */
    if ( listNames ) 
    {
      display_block_names(file_object);

      /* next file */
      ++optind;
      continue;
    }

    
    /* extract data block(s) by (range of) id(s) */
    if ( extractData && extractDataByID )
    {
      display_data_by_id(file_object, selectionString, print_times, 
          write_file);

      /* next file */
      ++optind;
      continue;
    }

    /* extract data block by data name */
    if ( extractData && extractDataByName )
    {
      display_data_by_name(file_object, selectionString, print_times, 
          write_file);

      /* next file */
      ++optind;
      continue;
    }

    /* extract data block(s) by regular expression */
    if ( extractData && extractDataByRegex )
    {
      display_data_by_regex(file_object, selectionString, print_times,
          write_file);

      /* next file */
      ++optind;
      continue;
    }


    /* next file */
    ++optind;
  }

  /* done */
  return EXIT_SUCCESS;
}



/****************************************************************
 * return a vector with the time points for each iteration
 * based on the time step 
 *****************************************************************/
vector<double> setup_time_list(const IO& file_object)
{
  vector<double> timeList;
  uint16_t type = file_object.get_output_type();
  if (type == STEP)
  {
    double timeStep = file_object.get_stepsize();
    double dataLength = file_object.get_blocksize();

    for ( int counter = 0; counter < dataLength; ++counter)
      timeList.push_back(timeStep*counter);
  }
  else if (type == TIME_LIST)
  {
    timeList = file_object.get_iteration_list();
  }
  else 
  {
    cerr << "Error: no time values available for this data set" << endl;
    throw;
  }

  return timeList;
}


//----------------------------------------------------------------
// return a vector with the iteration values for each iteration
// based on the time step 
//----------------------------------------------------------------
vector<int> setup_iteration_list(const IO& file_object)
{
  vector<double> temp = file_object.get_iteration_list();
  vector<int> iteration_list;
  foreach(double item, temp) 
  {
    iteration_list.push_back(static_cast<int>(item));
  }
  
  return iteration_list;
}


//----------------------------------------------------------------
// utility function printing out usage + commandline flags 
//----------------------------------------------------------------
void usage()
{
  cout << "usage: reader [options] <binary count file> " << NL << NL
       << "options:" << NL
       << "\t -i                   provide general info " << NL
       << "\t -l                   list available data blocks"  << NL
       << "\t -I num[:num1]        extract data by their index" << NL 
       << "\t                      or range of indices" << NL
       << "\t -N <data set name>   extract data by data set name" << NL
       << "\t -R <regexp>          extract data by regular expression" << NL
       << "\t -w                   write extracted filedata to ASCII file" << NL
       << "\t                      with original name" << NL
       << "\t -t                   add iteration times" << NL
       << "\t -h                   help" << NL
       << endl;
  exit(1);
}


//----------------------------------------------------------------
// display file info 
//----------------------------------------------------------------
void display_file_info(const IO& file_object)
{
  file_object.info();
}


//----------------------------------------------------------------
// display block names
//----------------------------------------------------------------
void display_block_names(const IO& file_object)
{
  //int numDataBlocks = file_object.get_num_datablocks();
  vector<std::string> blockNames = file_object.get_blocknames(); 	
  
   /* list block names */
   int counter = 0;
   vector<string>::iterator myPos;
   for ( myPos = blockNames.begin(); myPos != blockNames.end(); ++myPos)
   {
     cout << "[" << counter++ << "]" << " " << (*myPos) << NL;
   }
   cout << endl;
}


//------------------------------------------------------------------
// display data set selected by name
//------------------------------------------------------------------
void display_data_by_name(const IO& file_object, string selectionString,
    bool print_times, bool write_file)
{
  vector<vector<double> > data;
  vector<uint16_t> data_types;
  if (!file_object.get_data(data, data_types, selectionString))
  {
    cerr << "ERROR: Could not retrieve data" << endl;
    exit(EXIT_FAILURE);
  }

  vector<double> times;
  vector<int> iterations;
  bool has_iterations = 
    (file_object.get_output_type() == ITERATION_LIST) ? true : false;
  if (print_times)
  {
    if (has_iterations) 
      iterations = setup_iteration_list(file_object);
    else
      times = setup_time_list(file_object);
  }


  if (write_file)
  {
    if (has_iterations) 
      write_data(data, data_types, iterations, selectionString);
    else
      write_data(data, data_types, times, selectionString);
  }
  else
  {
    if (has_iterations)
      write_data(data, data_types, iterations);
    else
      write_data(data, data_types, times);
  }
}


//------------------------------------------------------------------
// display data set(s) selected by regex
//------------------------------------------------------------------
void display_data_by_regex(const IO& file_object, string selectionString,
    bool print_times, bool write_file)
{
  /* set up regular expression */
  boost::regex searchRegex(selectionString);

  /* loop through all blocknames and display all matches */
  vector<string> blockNames = file_object.get_blocknames(); 	
  foreach(string blockName, blockNames)
  {
    if ( regex_match(blockName, searchRegex) )
    {
      display_data_by_name(file_object, blockName, print_times, write_file);
    }
  }
}


//------------------------------------------------------------------
// display data set(s) selected by ID
//------------------------------------------------------------------
void display_data_by_id(const IO& file_object, string selectionString,
    bool print_times, bool write_file)
{
  int numDataBlocks = file_object.get_num_datablocks();
  vector<std::string> blockNames = file_object.get_blocknames(); 	

  /* check if user provided a range */
  int start;
  int end;
  std::string::size_type position;
  if ((position = selectionString.find(":")) != std::string::npos )
  {
    start = atoi(selectionString.substr(0,position).c_str());
    end = atoi(selectionString.substr(position+1,selectionString.size()).c_str());

    if ( start > end )
    {
      std::cout << "ERROR: Invalid range" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else 
  {
    start = atoi(selectionString.c_str());
    end = start;
  }

  /* make sure we are in range with our data */
  if ((start < 0 || start >= numDataBlocks) 
   || (end < 0 || end >= numDataBlocks))
  {
    std::cout << "ERROR: Data out of range" << std::endl;
    exit(EXIT_FAILURE);
  }


  /* get time or iterations of output */
  vector<double> times;
  vector<int> iterations;
  bool has_iterations = 
    (file_object.get_output_type() == ITERATION_LIST) ? true : false;
  if (print_times)
  {
    if (has_iterations)
      iterations = setup_iteration_list(file_object);
    else
      times = setup_time_list(file_object);
  }

  /* loop over all data sets and write them to stdout or file */
  for ( int counter = start; counter <= end; ++counter)
  {
    vector<vector<double> > data;
    vector<uint16_t> data_types;
    if (!file_object.get_data(data, data_types, counter))
    {
      cerr << "ERROR: Could not retrieve data" << endl;
      exit(EXIT_FAILURE);
    }

    if (write_file)
    {
      if (has_iterations) 
        write_data(data, data_types, iterations, blockNames[counter]);
      else
        write_data(data, data_types, times, blockNames[counter]);
    }
    else
    {
      if (has_iterations)
        write_data(data, data_types, iterations);
      else
        write_data(data, data_types, times);
    }
  }
}



/*************************************************************
 * write a data array to stdout or disk 
 *************************************************************/
template <typename T> void 
write_data(const vector<vector<double> >& data, 
    const vector<uint16_t>& data_types, const vector<T>& times, 
    string file_name)
{
  if (data.size() == 0)
    return;

  std::streambuf* stream = cout.rdbuf();
  std::ofstream of;
  if (!file_name.empty())
  {
    of.open(file_name.c_str());
    stream = of.rdbuf();
  }
  std::ostream out(stream);
  out.precision(10);
 
  if (!times.empty())
  {
    for (size_t row = 0; row < data[0].size(); ++row) 
    {
      out << std::fixed << times[row] << "  ";
      for (size_t col = 0; col < data_types.size(); ++col)
      {
        if (data_types[col] == INT_DATA) 
          out << static_cast<int>(data[col][row]) << "  ";
        else
          out << data[col][row] << "  ";
      }
      out << NL;
    }
  }
  else
  {
    for (size_t row = 0; row < data[0].size(); ++row) 
    {
      for (size_t col = 0; col < data_types.size(); ++col)
      {
        if (data_types[col] == INT_DATA) 
          out << static_cast<int>(data[col][row]) << "  ";
        else
          out << data[col][row] << "  ";
      }
      out << NL;
    }
  }
}
