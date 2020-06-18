/*
 * parser.h
 *
 *  Created on: Jun 17, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_PARSER_H_
#define LIBS_BNG_PARSER_H_

#include <string>

namespace BNG {

class BNGData;

// parses input BNGL file and performs semantic analysis
// returns number of errors encountered while parsing
// prints errors and warnings directly to the error output
// !! function is not reentrant
int parse_bngl_file(const std::string& file_name, BNGData& bng_data);

} /* namespace BNG */

#endif /* LIBS_BNG_PARSER_H_ */
