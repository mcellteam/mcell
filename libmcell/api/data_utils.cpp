/******************************************************************************
*
* Copyright (C) 2020 by
* The Salk Institute for Biological Studies
*
* Use of this source code is governed by an MIT-style
* license that can be found in the LICENSE file or at
* https://opensource.org/licenses/MIT.
*
******************************************************************************/

#include <fstream>
#include <cerrno>

#include "generated/gen_data_utils.h"

using namespace std;

namespace MCell {
namespace API {

namespace data_utils {

// trim from start (in place)
static void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}


// trim from end (in place)
static void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}


// trim from both ends (in place)
static void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}


std::vector<std::vector<double>> load_dat_file(const std::string& file_name) {
  std::vector<std::vector<double>> res;

  ifstream fin;
  fin.open(file_name, ios::in);

  int linenr = 1;
  string line;
  while(getline(fin, line)) {
    // is only whitespace? -> skip
    trim(line);
    if (line == "") {
      continue;
    }

    size_t pos = line.find_first_of(" \t");
    if (pos == string::npos) {
      throw RuntimeError("Invalid format of input file " + file_name +
          ", did not find a column separator on line " + to_string(linenr) + ".");
    }

    string col1 = line.substr(0, pos);
    string col2 = line.substr(pos + 1);

    char* end;
    errno = 0;
    double col1_value = strtod(col1.c_str(), &end);
    if (errno != 0 || *end != '\0') {
      throw RuntimeError("Invalid format of input file " + file_name +
          ", could not convert value in the first column " + col1 + " to a float, error on line " +
          to_string(linenr) + ".");
    }

    errno = 0;
    double col2_value = strtod(col2.c_str(), &end);
    if (errno != 0 || *end != '\0') {
      throw RuntimeError("Invalid format of input file " + file_name +
          ", could not convert value in the second column " + col2 + " to a float, error on line " +
          to_string(linenr) + ".");
    }

    res.push_back(std::vector<double>());
    res.back().push_back(col1_value);
    res.back().push_back(col2_value);

    linenr++;
  }

  return res;
}

} // namespace data_utils

} // namespace API
} // namespace MCell
