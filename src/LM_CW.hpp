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

#ifndef LM_CW_H
#define LM_CW_H

#include <Eigen/Dense>
#include <brent_fmin.hpp>
#include <cmath>
#include <functional>

class LM_CW {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> ns;
        double df_resid;

    public:

        double delta;
        double sigma2;
        Eigen::VectorXd beta;
        Eigen::VectorXd w;

        LM_CW(const Eigen::Ref<Eigen::MatrixXd> X_,
              const Eigen::Ref<Eigen::VectorXd> y_,
              const Eigen::Ref<Eigen::VectorXd> ns_) :
            X(X_),
            y(y_),
            ns(ns_)
        {
            df_resid = X.rows() - X.cols();
        };

        double neg_ll_reml(double delta_arg) {
            Eigen::VectorXd w_vals = (ns.array().inverse() + delta_arg).inverse().matrix();
            Eigen::MatrixXd XtWX = X.transpose() * w_vals.asDiagonal() * X;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
            Eigen::VectorXd beta_hat = ldlt.solve(X.transpose() * w_vals.cwiseProduct(y));
            Eigen::VectorXd r = y - X * beta_hat;
            double sigma2_ = r.cwiseAbs2().dot(w_vals) / df_resid;
            double ll = 0.5 * (
                df_resid * std::log(sigma2_) -
                w_vals.array().log().sum() +
                ldlt.vectorD().array().log().sum()
            );
            return ll;
        };

        void fit() {

            std::function<double(double)> f = [this](double log_delta) {
                return neg_ll_reml(std::exp(log_delta));
            };
            delta = std::exp(Brent_fmin(std::log(1e-6), std::log(1e3), f, 1e-4));

            w = (ns.array().inverse() + delta).inverse().matrix();
            Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
            beta = ldlt.solve(X.transpose() * w.cwiseProduct(y));
            Eigen::VectorXd r = y - X * beta;
            sigma2 = r.cwiseAbs2().dot(w) / df_resid;
        }
};

#endif
