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

#ifndef QTLMAPPING_UTILS_HPP
#define QTLMAPPING_UTILS_HPP

#include "Geno.hpp"
#include "Data.hpp"
#include "LMM.hpp"
#include "GLMM_GRM.hpp"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <vector>
#include <iostream>
#include <numeric>

std::vector<int> rank_vector(const std::vector<double>& v);
void rank_normalize(Eigen::MatrixXd& Y);
void rank_normalize_vec(Eigen::VectorXd& y);

Eigen::VectorXi assign_quantile_groups(const Eigen::VectorXd& x, int n_groups);

double ACAT(const std::vector<double>& pvals);

double pnorm(double x, bool lower);
double pt(double x, double df, bool lower);
double qnorm(double p, bool lower);
double qcauchy(double p, bool lower);
double pcauchy(double x, bool lower);

bool is_glmm(const std::string& model);
bool is_poisson_glmm(const std::string& model);
Eigen::MatrixXd donor_mean_covariates(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& ns
);
std::vector<Eigen::MatrixXd> group_donor_mean_covariates(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& ns,
    const std::vector<int>& cell_to_group,
    size_t n_groups
);
Eigen::MatrixXd between_donor_columns(const Eigen::MatrixXd& X);
Eigen::MatrixXd between_donor_columns(const Eigen::MatrixXd& X, const Eigen::VectorXi& donor_counts);

struct CochranQResult {
    double q;
    double pvalue;
    int df;
};
CochranQResult compute_cochran_q(const std::vector<double>& beta, const std::vector<double>& se);

struct WeightedTrendResult {
    double beta;
    double se;
    double pvalue;
};
WeightedTrendResult compute_weighted_trend(
    const std::vector<double>& beta,
    const std::vector<double>& se,
    const std::vector<double>& score
);

std::vector<double> make_group_linear_scores_values(const std::vector<double>& values);

std::string make_variant_header_line(const Params& params, const std::vector<std::string>& group_ids = {}, bool has_group_values = false);
std::string make_region_header_line(const Params& params, const std::vector<std::string>& group_ids = {}, bool has_group_values = false);

#endif