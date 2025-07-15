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

#ifndef QTLMAPPING_UTILS_HPP
#define QTLMAPPING_UTILS_HPP

#include "Geno.hpp"
#include "Data.hpp"
#include "LMM.hpp"
#include "GLMM.hpp"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <vector>
#include <iostream>
#include <numeric>

std::vector<int> rank_vector(const std::vector<double>& v);
void rank_normalize(Eigen::MatrixXd& Y);

double ACAT(const std::vector<double>& pvals);

double pnorm(double x, bool lower);
double qnorm(double p, bool lower);
double qcauchy(double p, bool lower);
double pcauchy(double x, bool lower);

std::string make_variant_header_line(std::string& model);

#endif