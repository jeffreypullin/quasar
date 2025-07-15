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

#ifndef RESID_HPP
#define RESID_HPP

#include "ModelFit.hpp"
#include "Data.hpp"
#include "Geno.hpp"

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalXd;

void residualise(Params& params, ModelFit& model_fit, CovData& cov_data, PhenoData& pheno_data, GRM& grm);

#endif