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

#ifndef QUASAR_H
#define QUASAR_H

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <string>

#include "Geno.hpp"
#include "Data.hpp"

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/StdVector"
#include "cxxopts.hpp"

struct Params {

    std::string plink_prefix;
    std::string grm_file;
    std::string cov_file;
    std::string bed_file;
    std::string resid_file;
    std::string fit_file;

    std::string mode; 
    std::string model;
    int window_size = 1000000;
    bool use_apl;

    std::string out;
    bool verbose = false;
};

#endif