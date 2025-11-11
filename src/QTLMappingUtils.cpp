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

#include "QTLMappingUtils.hpp"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <vector>
#include <iostream>
#include <numeric>

double pnorm(double x, bool lower) {
    boost::math::normal N01(0.0, 1.0);
    if (lower) { 
        return boost::math::cdf(boost::math::complement(N01, x));
    } 
    return boost::math::cdf(N01, x);
}

double p_bd = 1e-300;
double q_bd = 3e+299;

double qnorm(double p, bool lower){
    boost::math::normal N01(0.0, 1.0);
    if (lower) { 
        return boost::math::quantile(boost::math::complement(N01, p));
    } 
    return boost::math::quantile(N01, p);
}

double qcauchy(double p, bool lower){
    p = p > p_bd ? p : p_bd;
    p = p < 1 - p_bd ? p : 1 - p_bd;

    boost::math::cauchy C01(0.0, 1.0);
    if (lower) {
        return boost::math::quantile(boost::math::complement(C01, p));
    }
    return boost::math::quantile(C01, p);
}

double pcauchy(double x, bool lower){
    x = x < q_bd ? x : q_bd;
    x = x > -q_bd ? x : -q_bd;

    boost::math::cauchy C01(0.0, 1.0);
    if (lower) {
        return boost::math::cdf(boost::math::complement(C01, x));
    }
    return boost::math::cdf(C01, x);
}

void rank_normalize(Eigen::MatrixXd& Y){
    double n = Y.rows();
    double p = Y.cols();
    
    std::vector<double> z((int) n);
    std::vector<double> rk((int) n);
    
    double mu = 0;
    double sd = 0;
    for (int i = 0; i < n; ++i){
        z[i] = qnorm(((double)i+1.0) / ((double)n+1.0), true);
        mu += z[i];
        sd += z[i] * z[i];
    }
    sd = std::sqrt(sd / (n - 1) - mu * mu / (n * (n - 1.0)));
    mu = mu / n;
    for (int i = 0; i < n; ++i){
        z[i] = (z[i] - mu) / sd;
    }

    for (int j = 0; j < p; ++j){

        std::vector<double> v(n);
        for (int i = 0; i < n; ++i){
            v[i] = Y(i, j);
        }

        std::vector<int> ranks = rank_vector(v);
        for (int i = 0; i < n; ++i){
            Y(i, j) = z[ranks[i] - 1];
        }
    }
}

std::vector<int> rank_vector(const std::vector<double>& v){
    
    std::vector<size_t> w(v.size());
    std::iota(w.begin(), w.end(), 0);
    std::sort(w.begin(), w.end(), [&v](size_t i, size_t j) { return v[i] < v[j]; });

    std::vector<int> r(w.size());
    for (size_t n, i = 0; i < w.size(); i += n) {
        n = 1;
        while (i + n < w.size() && v[w[i]] == v[w[i+n]]) ++n;
        for (size_t k = 0; k < n; ++k) {
            r[w[i+k]] = i + (n + 1) / 2.0;
        }
    }
    return r;
}

double ACAT(const std::vector<double>& pvals) {
    long double sum = 0.0;
    double n = pvals.size();
    for (const double& p: pvals) {
        if (p >= 1){
            sum += (qcauchy(1 - 1 / n, true) / n);
        } else if (p <= 0 || std::isnan(p)){
            continue;
        } else {
            sum += (qcauchy(p, true) / n);
        }
    }
    return pcauchy(sum, true);
}

std::string make_variant_header_line(std::string& model) {

    std::string line = "feature_id\tsnp_id\tchrom\tpos\tref\talt\tmaf\tbeta\tse\tpvalue";

    if (model == "p_glm") {
        line = line + "\tglm_converged";
    } else if (model == "nb_glm") {
        line = line + "\tglm_converged\tphi\tphi_converged";
    } else if (model == "p_glmm") {
        line = line + "\tglmm_converged";
    } else if (model == "nb_glmm") {
        line = line + "\tglmm_converged\tphi\tphi_converged";
    }

    line = line + "\n";
    return line;
}