/*
  Copyright (C) 2024-26 Jeffrey Pullin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 or 3 of the License
  (at your option).
*/

#ifndef SADDLEPOINT_H
#define SADDLEPOINT_H

#include <Eigen/Dense>

double poisson_saddlepoint_pvalue(
    const Eigen::VectorXd& mu,
    const Eigen::VectorXd& adjusted_genotype,
    double score,
    double score_variance,
    double normal_pvalue,
    double cutoff
);

#endif
