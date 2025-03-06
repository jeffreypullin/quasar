
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
Eigen::VectorXd standardise_vec(const Eigen::VectorXd& x);

double ACAT(const std::vector<double>& pvals);

double pnorm(double x, bool lower);
double qnorm(double p, bool lower);
double qcauchy(double p, bool lower);
double pcauchy(double x, bool lower);
double qchisq(double p, double df, bool lower);
double pchisq(double x, double df, bool lower);

#endif