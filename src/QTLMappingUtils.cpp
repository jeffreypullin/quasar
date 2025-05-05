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

double qchisq(double p, double df, bool lower){
    boost::math::chi_squared CHISQ(df);
    if (lower) { 
        return boost::math::quantile(boost::math::complement(CHISQ, p));
    }
    return boost::math::quantile(CHISQ, p);
}

double pchisq(double x, double df, bool lower){
    boost::math::chi_squared CHISQ(df);
    if (lower) { 
        return boost::math::cdf(boost::math::complement(CHISQ, x));
    }
    return boost::math::cdf(CHISQ, x);
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
            sum += qcauchy(1 - 1 / n, true);
        } else if (p <= 0){
            // FIXME: Handle this case better.
            continue;
        } else {
            sum += qcauchy(p, true);
        }
    }
    return pcauchy(sum, true);
}

void standardise_vec(Eigen::VectorXd& x) {
    double mean = x.mean();
    double sd = std::sqrt((x.array() - mean).square().sum() / (x.size() - 1));
    if (sd < 1e-10) {
        x.setZero();
    } else {
        x.array() -= mean;
        x.array() /= sd;
    }
}

std::string make_variant_header_line(std::string& model) {

    std::string line = "feature_id\tchrom\tpos\tref\talt\tbeta\tse\tpvalue";

    if (model == "nb_glm") {
        line = line + "\tphi\tglm_converged\tphi_converged";
    }

    line = line + "\n";
    return line;
}