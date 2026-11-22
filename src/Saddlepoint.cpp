/*
  Copyright (C) 2024-26 Jeffrey Pullin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 or 3 of the License
  (at your option).
*/

#include "Saddlepoint.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

double normal_tail(double z) {
    return 0.5 * std::erfc(z / std::sqrt(2.0));
}

void cumulant_derivatives(
    double t,
    const Eigen::VectorXd& mu,
    const Eigen::VectorXd& g,
    double& first,
    double& second
) {
    first = 0.0;
    second = 0.0;
    for (Eigen::Index i = 0; i < mu.size(); ++i) {
        const double x = g(i) * t;
        if (x > 700.0) {
            first = std::copysign(std::numeric_limits<double>::infinity(), g(i));
            second = std::numeric_limits<double>::infinity();
            return;
        }
        const double exp_x = std::exp(x);
        first += mu(i) * g(i) * (exp_x - 1.0);
        second += mu(i) * g(i) * g(i) * exp_x;
    }
}

// Newton-Raphson from t = 0 with the step-halving safeguard of SPAtest/SAIGE
// (getroot_K1): if K'(t) - q changes sign with a jump no smaller than the
// previous one, halve the previous jump instead.
bool find_saddlepoint(
    double q,
    const Eigen::VectorXd& mu,
    const Eigen::VectorXd& g,
    double& root
) {
    const double tol = 1e-10;
    double t = 0.0;
    double k1 = 0.0;
    double k2 = 0.0;
    cumulant_derivatives(t, mu, g, k1, k2);
    k1 -= q;
    double prev_jump = std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < 100; ++iter) {
        double t_new = t - k1 / k2;
        if (!std::isfinite(t_new)) {
            return false;
        }
        if (std::abs(t_new - t) < tol) {
            root = t_new;
            return true;
        }

        double k1_new = 0.0;
        double k2_new = 0.0;
        cumulant_derivatives(t_new, mu, g, k1_new, k2_new);
        k1_new -= q;
        if (std::signbit(k1) != std::signbit(k1_new)) {
            if (std::abs(t_new - t) > prev_jump - tol) {
                t_new = t + std::copysign(prev_jump / 2.0, k1_new - k1);
                cumulant_derivatives(t_new, mu, g, k1_new, k2_new);
                k1_new -= q;
                prev_jump /= 2.0;
            } else {
                prev_jump = std::abs(t_new - t);
            }
        }
        t = t_new;
        k1 = k1_new;
        k2 = k2_new;
    }
    return false;
}

bool saddlepoint_tail(
    double q,
    const Eigen::VectorXd& mu,
    const Eigen::VectorXd& g,
    double& probability
) {
    double t = 0.0;
    if (!find_saddlepoint(q, mu, g, t)) {
        return false;
    }

    double cumulant = 0.0;
    double second_derivative = 0.0;
    for (Eigen::Index i = 0; i < mu.size(); ++i) {
        const double x = g(i) * t;
        const double exp_x = std::exp(x);
        cumulant += mu(i) * (std::expm1(x) - x);
        second_derivative += mu(i) * g(i) * g(i) * exp_x;
    }

    const double deviation = t * q - cumulant;
    if (!std::isfinite(deviation) || !std::isfinite(second_derivative) ||
        deviation <= 0.0 || second_derivative <= 0.0) {
        return false;
    }

    const double w = std::copysign(std::sqrt(2.0 * deviation), t);
    const double v = t * std::sqrt(second_derivative);
    const double ratio = v / w;
    if (!std::isfinite(w) || !std::isfinite(ratio) || w == 0.0 || ratio <= 0.0) {
        return false;
    }

    const double z = w + std::log(ratio) / w;
    if (!std::isfinite(z)) {
        return false;
    }
    probability = z > 0.0 ? normal_tail(z) : normal_tail(-z);
    return std::isfinite(probability);
}

} // namespace

double poisson_saddlepoint_pvalue(
    const Eigen::VectorXd& mu,
    const Eigen::VectorXd& adjusted_genotype,
    double score,
    double score_variance,
    double normal_pvalue,
    double cutoff
) {
    if (mu.size() == 0 || mu.size() != adjusted_genotype.size() ||
        !std::isfinite(score) || !std::isfinite(score_variance) ||
        score_variance <= 0.0) {
        return normal_pvalue;
    }

    const double standardised_score = std::abs(score) / std::sqrt(score_variance);
    if (standardised_score <= cutoff) {
        return normal_pvalue;
    }

    double poisson_variance = 0.0;
    for (Eigen::Index i = 0; i < mu.size(); ++i) {
        if (!std::isfinite(mu(i)) || mu(i) < 0.0 ||
            !std::isfinite(adjusted_genotype(i))) {
            return normal_pvalue;
        }
        poisson_variance += mu(i) * adjusted_genotype(i) * adjusted_genotype(i);
    }
    if (!std::isfinite(poisson_variance) || poisson_variance <= 0.0) {
        return normal_pvalue;
    }

    const double q = score * std::sqrt(poisson_variance / score_variance);
    double upper_tail = 0.0;
    double lower_tail = 0.0;
    if (!saddlepoint_tail(q, mu, adjusted_genotype, upper_tail) ||
        !saddlepoint_tail(-q, mu, adjusted_genotype, lower_tail)) {
        return normal_pvalue;
    }

    return std::min(1.0, upper_tail + lower_tail);
}
