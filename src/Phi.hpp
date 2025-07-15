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

#ifndef PHI_H
#define PHI_H

#include <Eigen/Dense>

double estimate_phi_ml(Eigen::VectorXd y, Eigen::VectorXd mu, bool& phi_converged);
double estimate_phi_apl(Eigen::VectorXd, Eigen::VectorXd mu, const Eigen::MatrixXd& X, bool& phi_converged);
double theta_score(Eigen::VectorXd y, Eigen::VectorXd mu, double theta);
double theta_info(Eigen::VectorXd y, Eigen::VectorXd mu, double theta);

#endif
