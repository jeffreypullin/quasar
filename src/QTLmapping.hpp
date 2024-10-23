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

#ifndef QTLMAPPING_HPP
#define QTLMAPPING_HPP

#include "Geno.hpp"
#include "Data.hpp"
#include "Regions.hpp"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalXd;

void run_qtl_mapping_lmm(GenoData& geno_data, FeatData& feat_data, CovData& cov_data, PhenoData& pheno_data, GRM& grm, Regions& regions);
double estimate_r(const Eigen::MatrixXd& X, const Eigen::MatrixXd& GRM, const Eigen::MatrixXd& G, double sigma2_g, double delta);
std::vector<int> rank_vector(const std::vector<double>& v);
void rank_normalize(Eigen::MatrixXd& Y);

double ACAT(const std::vector<double>& pvals);

double pnorm(double x, bool lower);
double qnorm(double p, bool lower);
double qcauchy(double p, bool lower);
double pcauchy(double x, bool lower);

#endif