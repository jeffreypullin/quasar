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

#ifndef GLMM_SC_H
#define GLMM_SC_H

#include <Eigen/Dense>
#include <limits>
#include <iostream>
#include "GLM.hpp"
#include "Family.hpp"

class GLMM_SC {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> offset;
        std::unique_ptr<Family> family;
        const Eigen::Ref<Eigen::VectorXd> ns;
        size_t N;
        size_t c;
        size_t n;
        double tol = 1e-5;
        int max_iter = 50;

	public:

        Eigen::VectorXd cum_ns;
        Eigen::VectorXd w;
        Eigen::MatrixXd XtSigma_invX_inv;
        Eigen::MatrixXd XSigma_inv;
        Eigen::VectorXd y_tilde;
        std::vector<Eigen::MatrixXd> Sigma_invXs;
        std::vector<Eigen::MatrixXd> Xis;
        std::vector<Eigen::VectorXd> wis;
        Eigen::MatrixXd ZtSigma_invX;
        Eigen::VectorXd eta;
        Eigen::VectorXd mu;
        Eigen::VectorXd u;

        // Output.
        double r_approx;
        Eigen::VectorXd y_out;
        Eigen::VectorXd mu_out;
        Eigen::MatrixXd XtWX_inv;
        Eigen::VectorXd Xty_res;
        Eigen::MatrixXd XtWZ;

        // Parameters
        Eigen::VectorXd beta;
        Eigen::VectorXd beta_prev;
        double sigma2;
        double sigma2_prev;

        // Control parameters.
        bool glmm_converged;
        int iter;
        double step_size;

        void update_y_tilde() {
            Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);
            y_tilde = (eta - offset) + (y - mu).cwiseProduct(mu_eta_vec);
        }

        void update_eta() {
            eta = X * beta + offset + u;
            for (size_t i = 0; i < n; i++) {
                if (eta(i) < -30 || eta(i) > 30) {
                    eta(i) = std::numeric_limits<double>::epsilon();
                }
            }
        }

        void update_mu() {
            mu = family->invlink(eta);
        }

        void update_w() {
            Eigen::VectorXd v_vec = family->var(mu);
            Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);
            w = (v_vec.array() * mu_eta_vec.array().square()).inverse().matrix();
        }

        void compute_cum_ns() {
            cum_ns = Eigen::VectorXd::Zero(n);
            cum_ns(0) = 0;
            for (size_t i = 1; i < n; ++i) {
                cum_ns(i) = cum_ns(i - 1) + ns(i - 1);
            }
        }

        Eigen::VectorXd Sigma_inv_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            if (sigma2 <= 0.0) {
                return w.cwiseProduct(x);
            }
            const double tau = 1.0 / sigma2;
            for (size_t i = 0; i < n; ++i) {
                const Eigen::VectorXd w_i = wis[i];
                const Eigen::VectorXd x_i = x.segment(cum_ns(i), ns(i));
                const double w_sum = w_i.sum();
                const double w_dot_x = w_i.dot(x_i);
                out.segment(cum_ns(i), ns(i)) =
                    w_i.cwiseProduct(x_i) - (w_i * w_dot_x) / (tau + w_sum);
            }
            return out;
        }

        Eigen::MatrixXd Sigma_inv_X(Eigen::MatrixXd X) {
            Eigen::MatrixXd out = Eigen::MatrixXd::Zero(X.rows(), X.cols());
            for (size_t j = 0; j < c; ++j) {
                out.col(j) = Sigma_inv_x(X.col(j));
            }
            return out;
        }

        void compute_Xis() {
            Xis.clear();
            for (size_t i = 0; i < n; ++i) {
                Xis.push_back(X.block(cum_ns(i), 0, ns(i), X.cols()));
            }
        }
        
        void compute_wis() {
            wis.clear();
            for (size_t i = 0; i < n; ++i) {
                wis.push_back(w.segment(cum_ns(i), ns(i)));
            }
        }

        Eigen::MatrixXd Sigma_inv_X_i(size_t i) {
            Eigen::MatrixXd out = Eigen::MatrixXd::Zero(Xis[i].rows(), X.cols());
            for (size_t j = 0; j < c; ++j) {
                Eigen::VectorXd w_i = wis[i];
                Eigen::VectorXd x_ij = Xis[i].col(j);
                double w_sum = w_i.sum();
                if (sigma2 <= 0.0) {
                    out.col(j) = w_i.cwiseProduct(x_ij);
                    continue;
                }
                double tau = 1.0 / sigma2;
                double w_dot_x = w_i.dot(x_ij);
                out.col(j) = w_i.cwiseProduct(x_ij) - (w_i * w_dot_x) / (tau + w_sum);
            }
            return out;
        }

        void update_Sigma_invXs() {
            Sigma_invXs.clear();
            for (size_t i = 0; i < n; ++i) {
                Sigma_invXs.push_back(Sigma_inv_X_i(i));
            }
        }

        void update_XtSigma_invX_inv() {
            Eigen::MatrixXd XtSigma_invX = Eigen::MatrixXd::Zero(c, c);
            for (size_t i = 0; i < n; ++i) {
                XtSigma_invX += Xis[i].transpose() * Sigma_invXs[i];
            }
            XtSigma_invX_inv = XtSigma_invX.inverse();
        }

        Eigen::VectorXd P_x(Eigen::VectorXd x) {
            Eigen::VectorXd a = Sigma_inv_x(x);
            Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
            Eigen::VectorXd tmp = Eigen::VectorXd::Zero(c);
            for (size_t i = 0; i < n; ++i) {
                tmp += Sigma_invXs[i].transpose() * x.segment(cum_ns(i), ns(i));
            }
            tmp = XtSigma_invX_inv * tmp;
            for (size_t i = 0; i < n; ++i) {
                b.segment(cum_ns(i), ns(i)) = Sigma_invXs[i] * tmp;
            }
            Eigen::VectorXd result = a - b;
            return result;
        }

        Eigen::VectorXd Omega_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                double sum_i = x.segment(cum_ns(i), ns(i)).sum();
                Eigen::VectorXd tmp = Eigen::VectorXd::Constant(ns(i), sum_i);
                out.segment(cum_ns(i), ns(i)) = tmp;
            }
            return out;
        }

        Eigen::MatrixXd Omega_X(Eigen::MatrixXd X) {
            Eigen::MatrixXd out = Eigen::MatrixXd::Zero(X.rows(), X.cols());
            for (int j = 0; j < X.cols(); ++j) {
                out.col(j) = Omega_x(X.col(j));
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

        void update_ZtSigma_invX() {
            ZtSigma_invX = Eigen::MatrixXd::Zero(n, c);
            for (size_t i = 0; i < n; ++i) {
                ZtSigma_invX.row(i) = Sigma_invXs[i].colwise().sum();
            }
        }

        double compute_trPOmega() {
            double tr_Sigma_invOmega = 0;
            double tr_Sigma_invRest = 0;

            for (size_t i = 0; i < n; ++i) {
                double wi_sum = wis[i].sum();
                tr_Sigma_invOmega += wi_sum / (1 + sigma2 * wi_sum);
            }

            Eigen::MatrixXd tmp;
            tmp = (XtSigma_invX_inv) * ZtSigma_invX.transpose();
            tr_Sigma_invRest = (ZtSigma_invX * tmp).trace();

            return tr_Sigma_invOmega - tr_Sigma_invRest; 
        }

        void update_beta() {
            beta_prev = beta;
            beta = XtSigma_invX_inv * X.transpose() * (Sigma_inv_x(y_tilde));
        }

        void update_u() {
            u = sigma2 * Omega_x(Sigma_inv_x(y_tilde - (X * beta)));
        }

        void update_sigma2() {
            sigma2_prev = sigma2;

            double score, ai;
            Eigen::VectorXd Py_tilde = P_x(y_tilde);
            Eigen::VectorXd OmegaPy_tilde = Omega_x(Py_tilde);

            score = Py_tilde.dot(OmegaPy_tilde) - compute_trPOmega();
            ai = (OmegaPy_tilde.transpose() * P_x(OmegaPy_tilde));
            sigma2 += step_size * (score / ai);

            // Handle the case when sigma2 < 0 after the update.
            if (sigma2 < tol && sigma2_prev < tol) {
                sigma2 = 0.0;
            }

            double ss = step_size;
            while (sigma2 < 0.0) {
                ss = ss * 0.5;
                sigma2 = sigma2_prev + ss * (score / ai);

                if (sigma2 < tol && sigma2_prev < tol) {
                    sigma2 = 0.0;
                }
            }

            if (sigma2 < tol && sigma2_prev < tol) {
                sigma2 = 0.0;
            }  
        }

        void update_step_size() {
            if ((iter + 1) % 10 == 0) {
                step_size = 0.9 * step_size;
            }
        }

        void check_converge() {
            Eigen::VectorXd tol_vec_beta = Eigen::VectorXd::Constant(beta.size(), tol);
            double diff1 = ((beta - beta_prev).cwiseAbs().cwiseQuotient(
                (beta).cwiseAbs() + (beta_prev).cwiseAbs() + tol_vec_beta)).maxCoeff();
            double diff2 = std::abs(sigma2 - sigma2_prev) / (std::abs(sigma2) + std::abs(sigma2_prev) + tol);
            glmm_converged = (2 * std::max(diff1, diff2)) < tol;
        }

        double compute_trWPwOmega() {
            double a = w.sum();
            Eigen::MatrixXd WX = w.asDiagonal() * X;
            Eigen::MatrixXd XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
            Eigen::MatrixXd XtWOmegaWX = WX.transpose() * Omega_X(WX);
            double b = (XtWOmegaWX * XtWX_inv).trace();
            return a - b;
        }

        double compute_trPOmegaWPwOmega() {
            double a = 0;
            double b = 0;
            double c1 = 0;
            double c2 = 0;

            for (size_t i = 0; i < n; ++i) {
                double tau = 1 / sigma2;
                double wis_sum = wis[i].sum();
                a += pow(wis_sum, 2) - (pow(wis_sum, 3) / (tau + wis_sum));
            }

            Eigen::MatrixXd tmp = (XtSigma_invX_inv) * ZtSigma_invX.transpose();
            b = (ZtSigma_invX * tmp * collapse_vec(w).asDiagonal()).trace();

            Eigen::MatrixXd XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
            Eigen::MatrixXd OmegaWX = Omega_X(w.asDiagonal() * X);
            Eigen::MatrixXd Sigma_invOmegaWX = Sigma_inv_X(OmegaWX);
            Eigen::MatrixXd tmp1 = OmegaWX.transpose() * Sigma_invOmegaWX;
            c1 = (tmp1 * XtWX_inv).trace();

            Eigen::MatrixXd XtSigma_invOmegaWX = X.transpose() * Sigma_invOmegaWX;
            Eigen::MatrixXd tmp2 = XtSigma_invOmegaWX.transpose() * (XtSigma_invX_inv * XtSigma_invOmegaWX);
            c2 = (tmp2 * XtWX_inv).trace();

            return a - b - (c1 - c2);
        }

        double compute_trWPwOmegaWPwOmega() {
            double a = 0;

            for (size_t i = 0; i < n; ++i) {
                a += pow(wis[i].sum(), 2);
            } 

            Eigen::MatrixXd WX = w.asDiagonal() * X;
            Eigen::MatrixXd XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
            Eigen::MatrixXd tmp = (WX.transpose() * Omega_X(w.asDiagonal() * Omega_X(WX)));
            double b = (XtWX_inv * tmp).trace();

            Eigen::MatrixXd XtWOmegaWX = WX.transpose() * Omega_X(WX);
            double c = (XtWX_inv * XtWOmegaWX * XtWX_inv * XtWOmegaWX).trace();
            
            return a - 2 * b + c;
        }

        void compute_r_approx() {

            double tr_POmega = compute_trPOmega();
            double tr_WPwOmega = compute_trWPwOmega();
            
            double a = tr_POmega / tr_WPwOmega;
            double b = 2 * compute_trPOmegaWPwOmega()  / pow(tr_WPwOmega, 2);
            double c = ((2 * compute_trWPwOmegaWPwOmega()) * tr_POmega) / pow(tr_WPwOmega, 3);
            r_approx = a - b + c;
        }

        void compute_output() {
            mu_out = collapse_vec(mu);
            Eigen::VectorXd y_res = (y.array() - mu.array());
            y_out = collapse_vec(y_res);
            XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
            Xty_res = X.transpose()  * y_res;
            XtWZ = Eigen::MatrixXd::Zero(c, n);
            for (size_t i = 0; i < n; ++i) {
                XtWZ.col(i) = Xis[i].transpose() * wis[i];
            }
        }

        GLMM_SC(
            const Eigen::Ref<Eigen::MatrixXd> X_, 
            const Eigen::Ref<Eigen::VectorXd> y_, 
            const Eigen::Ref<Eigen::VectorXd> offset_,
            std::unique_ptr<Family> family_, 
            const Eigen::Ref<Eigen::VectorXd> ns_
        ) : 
            X(X_),
            y(y_),
            offset(offset_),
            family(std::move(family_)),
            ns(ns_)
        {   
            N = X.rows();
            c = X.cols();
            n = ns.size();
            init_params();
        };

        void init_params() {
            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLM glm(X, y, offset, std::move(poisson));
            glm.fit();
            beta = glm.beta;
            mu = glm.mu;
            compute_cum_ns();
            compute_Xis();
            u = Eigen::VectorXd::Zero(N);
            update_eta();
            update_w();
            compute_wis();
            update_y_tilde();
            sigma2 = 1;
            iter = 0;
            glmm_converged = false;
            step_size = 1;
        }

        void fit() {

            while (iter < max_iter) {
                
                update_Sigma_invXs();
                update_XtSigma_invX_inv();
                update_ZtSigma_invX();

                update_beta();
                update_u();
                update_eta();
                update_mu();

                update_sigma2();

                update_w();
                compute_wis();

                update_y_tilde();

                check_converge();
                if (glmm_converged) {
                    break;
                }
                update_step_size();
                iter += 1;
            }
            if (std::isnan(sigma2) || (beta.hasNaN())) {
                glmm_converged = false;
            }
            compute_r_approx();
            compute_output();
        }
};

#endif