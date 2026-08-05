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

#include "QTLMappingUtils.hpp"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <vector>
#include <iostream>
#include <numeric>

double pnorm(double x, bool lower) {
    boost::math::normal N01(0.0, 1.0);
    if (lower) { 
        return boost::math::cdf(N01, x);
    } 
    return boost::math::cdf(boost::math::complement(N01, x));
}

double p_bd = 1e-300;
double q_bd = 3e+299;

double qnorm(double p, bool lower){
    boost::math::normal N01(0.0, 1.0);
    if (lower) { 
        return boost::math::quantile(N01, p);
    } 
    return boost::math::quantile(boost::math::complement(N01, p));
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

void rank_normalize_vec(Eigen::VectorXd& y){
    double n = y.size();

    std::vector<double> z((int) n);

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

    std::vector<double> v(n);
    for (int i = 0; i < n; ++i){
        v[i] = y(i);
    }

    std::vector<int> ranks = rank_vector(v);
    for (int i = 0; i < n; ++i){
        y(i) = z[ranks[i] - 1];
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
    double nan_count = 0;
    for (const double& p: pvals) {
        if (p >= 1){
            sum += (qcauchy(1 - 1 / n, true) / n);
        } else if (p <= 0 || std::isnan(p)) {
            // We only want to throw NaN if all the p-values are NaN,
            // otherwise just calculate the ACAT on the non-NaN p-values.
            // This is because NaN's either arise for all variants
            // due to gene-level non-convergence or for single-variants
            // due to MAF/MAC issues.
            nan_count += 1;
            continue;
        } else {
            sum += (qcauchy(p, true) / n);
        }
    }
    if (std::abs(nan_count - n) < 1e-8) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return pcauchy(sum, true);
}

CochranQResult compute_cochran_q(const std::vector<double>& beta, const std::vector<double>& se) {

    CochranQResult res{std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN(),
                       0};

    if (beta.size() != se.size()) {
        return res;
    }

    double sum_w = 0.0;
    double sum_wb = 0.0;
    int k = 0;
    for (size_t i = 0; i < beta.size(); ++i) {
        if (std::isnan(beta[i]) || std::isnan(se[i]) || se[i] <= 0.0) {
            continue;
        }
        double w = 1.0 / (se[i] * se[i]);
        sum_w  += w;
        sum_wb += w * beta[i];
        k++;
    }

    if (k < 2 || sum_w <= 0.0) {
        return res;
    }

    double beta_fe = sum_wb / sum_w;
    double Q = 0.0;
    for (size_t i = 0; i < beta.size(); ++i) {
        if (std::isnan(beta[i]) || std::isnan(se[i]) || se[i] <= 0.0) {
            continue;
        }
        double w = 1.0 / (se[i] * se[i]);
        double d = beta[i] - beta_fe;
        Q += w * d * d;
    }

    res.q = Q;
    res.df = k - 1;
    boost::math::chi_squared chi(static_cast<double>(res.df));
    res.pvalue = boost::math::cdf(boost::math::complement(chi, Q));
    return res;
}

WeightedTrendResult compute_weighted_trend(
    const std::vector<double>& beta,
    const std::vector<double>& se,
    const std::vector<double>& score
) {
    WeightedTrendResult res{std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::quiet_NaN()};

    if (beta.size() != se.size() || beta.size() != score.size()) {
        return res;
    }

    double s0 = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    double sy = 0.0;
    double sxy = 0.0;
    int k = 0;

    for (size_t i = 0; i < beta.size(); ++i) {
        if (std::isnan(beta[i]) || std::isnan(se[i]) || std::isnan(score[i]) || se[i] <= 0.0) {
            continue;
        }
        double w = 1.0 / (se[i] * se[i]);
        s0 += w;
        s1 += w * score[i];
        s2 += w * score[i] * score[i];
        sy += w * beta[i];
        sxy += w * score[i] * beta[i];
        k++;
    }

    double det = s0 * s2 - s1 * s1;
    if (k < 2 || det <= 0.0 || std::isnan(det)) {
        return res;
    }

    res.beta = (s0 * sxy - s1 * sy) / det;
    res.se = std::sqrt(s0 / det);
    double z = res.beta / res.se;
    if (res.se <= 0.0 || std::isnan(z)) {
        return res;
    }
    res.pvalue = 2 * pnorm(std::abs(z), false);
    return res;
}

std::vector<double> make_group_linear_scores_values(const std::vector<double>& values) {
    std::vector<double> scores(values.size());
    double mean = 0.0;
    for (double v : values) mean += v;
    mean /= static_cast<double>(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        scores[i] = values[i] - mean;
    }
    return scores;
}

std::vector<double> make_group_quadratic_scores(const std::vector<double>& linear_scores) {
    std::vector<double> scores(linear_scores.size());
    double mean = 0.0;
    for (size_t i = 0; i < linear_scores.size(); ++i) {
        scores[i] = linear_scores[i] * linear_scores[i];
        mean += scores[i];
    }
    mean /= static_cast<double>(scores.size());
    for (double& score : scores) {
        score -= mean;
    }

    double lin_quad = 0.0;
    double lin_lin = 0.0;
    for (size_t i = 0; i < linear_scores.size(); ++i) {
        lin_quad += linear_scores[i] * scores[i];
        lin_lin += linear_scores[i] * linear_scores[i];
    }
    double slope = lin_lin > 0.0 ? lin_quad / lin_lin : 0.0;
    for (size_t i = 0; i < scores.size(); ++i) {
        scores[i] -= slope * linear_scores[i];
    }
    return scores;
}

std::string make_variant_header_line(const Params& params, const std::vector<std::string>& group_ids, bool has_group_values) {

    const std::string& model = params.model;

    std::string line = "feature_id\tsnp_id\tchrom\tpos\talt\tref\tmaf";
    
    if (params.do_interaction) {
        line = line + "\tsnp_beta\tsnp_se\tsnp_pvalue";
    } else {
        line = line + "\tbeta\tse\tpvalue";
    }

    if (model == "p_glm") {
        line = line + "\tglm_converged";
    } else if (model == "nb_glm") {
        line = line + "\tglm_converged\tphi\tphi_converged";
    } else if (model == "p_glmm" ||
               model == "p_glmm_grm" ||
               ((model == "p_glmm_sc") & !params.do_interaction)) {
        line = line + "\tglmm_converged\tsigma2";
    } else if ((model == "p_glmm_sc") & params.do_interaction) {
        line = line + "\tglmm_converged";
    } else if (model == "nb_glmm") {
        line = line + "\tglmm_converged\tsigma2\tphi\tphi_converged";
    }

    if (params.do_interaction) {
        std::string interaction_id = params.interaction_cov;

        std::string snake_case_id;
        snake_case_id.reserve(interaction_id.size());
        for (unsigned char c : interaction_id) {
            if (std::isalnum(c)) {
                snake_case_id.push_back(static_cast<char>(std::tolower(c)));
            } else if (!snake_case_id.empty() && snake_case_id.back() != '_') {
                snake_case_id.push_back('_');
            }
        }
        if (!snake_case_id.empty() && snake_case_id.back() == '_') {
            snake_case_id.pop_back();
        }

        line += "\tsnp_x_" + snake_case_id + "_beta";
        line += "\tsnp_x_" + snake_case_id + "_se";
        line += "\tsnp_x_" + snake_case_id + "_pvalue";
    }

    if (!group_ids.empty()) {
        for (const auto& gid : group_ids) {
            if (has_group_values) line += "\t" + gid + "_value";
            line += "\t" + gid + "_beta";
            line += "\t" + gid + "_se";
            line += "\t" + gid + "_pvalue";
        }
        line += "\tgroup_het_q\tgroup_het_pvalue";
        if (has_group_values) {
            line += "\tgroup_linear_beta\tgroup_linear_se\tgroup_linear_pvalue";
            line += "\tgroup_quadratic_beta\tgroup_quadratic_se\tgroup_quadratic_pvalue";
            line += "\tgroup_acat_pvalue";
        }
    }

    line = line + "\n";
    return line;
}