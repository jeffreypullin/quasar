/*
  Copyright (C) 2024-25 Jeffrey Pullin

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

// From regenie.
void remove_carriage_return(std::string& str) {
  if (!str.empty() && str.back() == '\r') {
    str.pop_back();
  }
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
