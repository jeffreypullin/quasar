/*
  Copyright (C) 2024-26 Jeffrey Pullin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 or 3 of the License
  (at your option).

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is available at
   http://www.r-project.org/Licenses/
*/

#include "Utils.hpp"

#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <unordered_set>
#include <iostream>

std::string format_with_commas(size_t number) {
    std::string numStr = std::to_string(number);
    
    int pos = numStr.length();
    while (pos > 3) {
        pos -= 3;
        numStr.insert(pos, ",");
    }

    return numStr;
}

// From regenie.
std::vector<std::string> string_split(std::string const& s, const char* delims) {

  std::vector<std::string> out;

  if (s.size() == 0) {
    return out;
  }

  const char* p = s.c_str();
  const char* q = strpbrk(p + 1, delims);

  for( ; q != NULL; q = strpbrk(p, delims)){
    out.push_back( std::string(p,q) );
    p = q + 1;
  }

  if (p && (p[0] != '\0')) {
    out.push_back(std::string(p));
  }

  return(out);
}

std::string trim_string(const std::string& s) {
  size_t start = 0;
  while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
    ++start;
  }
  size_t end = s.size();
  while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
    --end;
  }
  return s.substr(start, end - start);
}

std::vector<std::string> parse_comma_separated_names(const std::string& s) {
  std::vector<std::string> out;
  if (s.empty()) {
    return out;
  }

  std::vector<std::string> parts = string_split(s, ",");
  std::unordered_set<std::string> seen;
  for (const auto& part : parts) {
    std::string name = trim_string(part);
    if (name.empty()) {
      std::cerr << "Error: empty name in comma-separated --interaction list." << std::endl;
      std::exit(1);
    }
    if (!seen.insert(name).second) {
      std::cerr << "Error: duplicate interaction covariate '" << name << "'." << std::endl;
      std::exit(1);
    }
    out.push_back(name);
  }
  return out;
}

// From regenie.
void remove_carriage_return(std::string& str) {
  if (!str.empty() && str.back() == '\r') {
    str.pop_back();
  }
}

bool next_field(const char*& p, const char* end, const char* delims,
                const char*& field_begin, const char*& field_end) {
  if (p >= end) {
    return false;
  }
  field_begin = p;
  while (p < end && std::strchr(delims, *p) == nullptr) {
    ++p;
  }
  field_end = p;
  if (p < end) {
    ++p;  // consume one delimiter
  }
  return true;
}

double parse_double_field(const char* begin, const char* end) {
  if (begin == end) {
    std::cerr << "Error: Empty numeric field." << std::endl;
    exit(1);
  }
  char* endptr = nullptr;
  double value = std::strtod(begin, &endptr);
  if (endptr != end) {
    std::cerr << "Error: Failed to parse numeric field '"
              << std::string(begin, end) << "'." << std::endl;
    exit(1);
  }
  return value;
}

std::vector<std::string> intersection(std::vector<std::vector<std::string>> &vecs) {

  // Sort each vector.
  for (auto& vec : vecs) {
    std::sort(vec.begin(), vec.end());
  }
  
  auto last_intersection = vecs[0];
  std::vector<std::string> curr_intersection;
  for (std::size_t i = 1; i < vecs.size(); ++i) {
      std::set_intersection(last_intersection.begin(), last_intersection.end(),
          vecs[i].begin(), vecs[i].end(),
          std::back_inserter(curr_intersection));
      std::swap(last_intersection, curr_intersection);
      curr_intersection.clear();
  }
  return last_intersection;
}
