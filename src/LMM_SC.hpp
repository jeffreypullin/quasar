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

#ifndef LMM_SC_H
#define LMM_SC_H

#include <Eigen/Dense>
#include <brent_fmin.hpp>
#include <cmath>
#include <functional>

class LMM_SC {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> ns;
        double df_resid;

    public:

        Eigen::MatrixXd X_tilde;
        Eigen::VectorXd y_tilde;
        Eigen::MatrixXd XtX;
        Eigen::VectorXd Xty;
        double yty;

        double sigma2;
        double delta;
        Eigen::VectorXd beta;
        Eigen::VectorXd cum_ns;
        size_t N;
        size_t n;
        size_t p;

        Eigen::VectorXd y_out;
        Eigen::VectorXd mu_out;
        Eigen::MatrixXd XtWX_inv;
        Eigen::VectorXd Xty_res;
        Eigen::MatrixXd XtWZ;
        double r_approx;

        LMM_SC(const Eigen::Ref<Eigen::MatrixXd> X_, 
               const Eigen::Ref<Eigen::VectorXd> y_, 
               const Eigen::Ref<Eigen::VectorXd> ns_
        ) : 
            X(X_),
            y(y_),
            ns(ns_)
        {
            n = static_cast<size_t>(ns.size());
            N = X.rows();
            p = X.cols();
            df_resid = N - p;
            XtX = X.transpose() * X;
            Xty = X.transpose() * y;
            yty = y.dot(y);

            compute_cum_ns();

            X_tilde = Eigen::MatrixXd::Zero(n, p);
            for (int j = 0; j < X.cols(); ++j) {
                X_tilde.col(j) = orthogonal_collapse_vec(X.col(j));
            }
            y_tilde = orthogonal_collapse_vec(y);
       
        };

        // From GLMM_SC.hpp
        void compute_cum_ns() {
            cum_ns = Eigen::VectorXd::Zero(n);
            cum_ns(0) = 0;
            for (size_t i = 1; i < n; ++i) {
                cum_ns(i) = cum_ns(i - 1) + ns(i - 1);
            }
        }

        Eigen::VectorXd expand_orthogonal_vec(const Eigen::VectorXd& v) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                const Eigen::Index start = static_cast<Eigen::Index>(cum_ns(i));
                const Eigen::Index len = static_cast<Eigen::Index>(ns(i));
                out.segment(start, len).setConstant(v(i) / std::sqrt(ns(i)));
            }
            return out;
        }

        Eigen::VectorXd orthogonal_collapse_vec(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
            for (size_t i = 0; i < n; ++i) {
                out(i) = x.segment(cum_ns(i), ns(i)).sum() / std::sqrt(ns(i));
            }
            return out;
        }

        Eigen::VectorXd collapse_vec(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
            for (size_t i = 0; i < n; ++i) {
                out(i) = x.segment(cum_ns(i), ns(i)).sum();
            }
            return out;
        }

        double neg_ll_reml(double delta_arg) {
            Eigen::VectorXd psi_vals = delta_arg * ns.array() / (delta_arg * ns.array() + 1.0);
            Eigen::MatrixXd XtPsiX = X_tilde.transpose() * psi_vals.asDiagonal() * X_tilde;
            Eigen::VectorXd XtPsiy = X_tilde.transpose() * psi_vals.asDiagonal() * y_tilde;
            Eigen::MatrixXd XtVX = XtX - XtPsiX;
            Eigen::VectorXd XtVy = Xty - XtPsiy;
            
            Eigen::VectorXd beta_hat = XtVX.colPivHouseholderQr().solve(XtVy);

            double ytPsiy = y_tilde.dot(psi_vals.asDiagonal() * y_tilde);
            double sigma2_ = ((yty - ytPsiy) - XtVy.dot(beta_hat)) / df_resid;
            double ll = 0.5 * (
                df_resid * std::log(sigma2_) + 
                (1.0 + delta_arg * ns.array()).log().sum() +
                std::log(XtVX.determinant())
            );

            return ll;
        };

        void compute_r_approx() {

            double a, a1, a2, b;
            double tau = 1.0 / sigma2;
            Eigen::MatrixXd ZtX = Eigen::MatrixXd::Zero(n, p);
            for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(p); ++j) {
                ZtX.col(j) = collapse_vec(X.col(j));
            }

            Eigen::VectorXd psi_vals = delta * ns.array() / (delta * ns.array() + 1.0);
            Eigen::MatrixXd XtPsiX = X_tilde.transpose() * psi_vals.asDiagonal() * X_tilde;
            Eigen::MatrixXd XtSigma_invX_inv = (tau * (XtX - XtPsiX)).inverse();
            Eigen::VectorXd ns_sqrt_psi = ns.array().sqrt().matrix().cwiseProduct(psi_vals);
            Eigen::MatrixXd ZtSigma_invX = (tau * (ZtX - ns_sqrt_psi.asDiagonal() * X_tilde));

            a1 = tau * N - tau * (ns.array() * psi_vals.array()).sum();
            a2 = (ZtSigma_invX * XtSigma_invX_inv * ZtSigma_invX.transpose()).trace();
            a = a1 - a2;

            b = N - (ZtX * XtX.inverse() * ZtX.transpose()).trace();

            r_approx = a / b;
        }

        void fit() {

            std::function<double(double)> f = [this](double x) { 
                return neg_ll_reml(x);
            };
            delta = Brent_fmin(0.00, 10000, f, 2e-5);
            Eigen::VectorXd psi_vals = delta * ns.array() / (delta * ns.array() + 1.0);
            Eigen::MatrixXd XtPsiX = X_tilde.transpose() * psi_vals.asDiagonal() * X_tilde;
            Eigen::VectorXd XtPsiy = X_tilde.transpose() * psi_vals.asDiagonal() * y_tilde;
            Eigen::MatrixXd XtVX = XtX - XtPsiX;
            Eigen::VectorXd XtVy = Xty - XtPsiy;
            
            beta = XtVX.colPivHouseholderQr().solve(XtVy);
            double ytPsiy = y_tilde.dot(psi_vals.asDiagonal() * y_tilde);
            sigma2 = ((yty - ytPsiy) - XtVy.dot(beta)) / df_resid;

            Eigen::VectorXd y_res_1 = (y - X * beta) / sigma2;
            Eigen::VectorXd y_res_2 = expand_orthogonal_vec(
                psi_vals.asDiagonal() * y_tilde - psi_vals.asDiagonal() * X_tilde * beta
            ) / sigma2;
            Eigen::VectorXd y_res = y_res_1 - y_res_2;

            y_out = collapse_vec(y_res);
            mu_out = ns;
            XtWX_inv = XtX.inverse();
            Xty_res = X.transpose() * y_res;

            Eigen::MatrixXd X_collapsed = Eigen::MatrixXd::Zero(n, p);
            for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(p); ++j) {
                X_collapsed.col(j) = collapse_vec(X.col(j));
            }
            XtWZ = X_collapsed.transpose();
            compute_r_approx();
       
            return;
        }
};

#endif
