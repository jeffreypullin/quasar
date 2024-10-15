/*
 * MIT License
 * 
 * Copyright (c) 2024 Jeffrey Pullin
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "Utils.hpp"

#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <unordered_set>
#include <iostream>

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
