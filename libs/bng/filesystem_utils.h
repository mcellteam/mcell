/*
 * filesystem_utils.h
 *
 *  Created on: Jan 21, 2021
 *      Author: Adam
 */

#ifndef LIBS_BNG_FILESYSTEM_UTILS_H_
#define LIBS_BNG_FILESYSTEM_UTILS_H_

// no namespace?
#include <string>

// exists with error message if the directory could not be created
void make_dir_for_file_w_multiple_attempts(const std::string& file_path);

void make_dir_w_multiple_attempts(const std::string& dir_path);


#endif /* LIBS_BNG_FILESYSTEM_UTILS_H_ */
