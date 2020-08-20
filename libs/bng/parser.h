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
class CplxInstance;

// parses input BNGL file and performs semantic analysis
// returns number of errors encountered while parsing
// prints errors and warnings directly to the error output
// !! function is not reentrant
int parse_bngl_file(const std::string& file_name, BNGData& bng_data);


// parses input BNGL complex instance string and performs
// semantic analysis
// bng_data are not cleared and one can continue with adding
// complexes gradually
// prints errors and warnings directly to the error output
// resulting complex is stored into the res_cplx
int parse_single_cplx_instance_file(
    const std::string& cplx_instance_string, BNGData& bng_data,
    CplxInstance& res_cplx
);


} /* namespace BNG */

#endif /* LIBS_BNG_PARSER_H_ */
