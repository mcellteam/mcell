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

/** STL headers */
#include <fstream>
#include <string>
#include <vector>

/* forward declarations */
class file_object;

/* main routine for writing a data set to stdout or disk */
template<typename T> void write_data(
    const std::vector<std::vector<double> >& data,
    const std::vector<uint16_t>& data_types,
    const std::vector<T>& times,
    std::string file_name = "");


/* return a vector with the time points or iterations numbers
 * for each iteration based on the time step */
std::vector<double> setup_time_list(const IO& file_object);
std::vector<int> setup_iteration_list(const IO& file_object);

/* punch usage information */
void usage();

/* punch file info */
void display_file_info(const IO& file_object);

/* display block names */
void display_block_names(const IO& file_object);

/* display data by ID */
void display_data_by_id(const IO& file_object, 
    std::string selectionString, bool print_times, bool write_file);

/* display data by name */
void display_data_by_name(const IO& file_object, 
    std::string selectionString, bool print_times, bool write_file);

/* display data by name */
void display_data_by_regex(const IO& file_object, 
    std::string selectionString, bool print_times, bool write_file);
