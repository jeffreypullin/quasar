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

#ifndef SCGLMM_H
#define SCGLMM_H

#include <Eigen/Dense>
#include <limits>
#include <iostream>
#include "GLM.hpp"
#include "Family.hpp"

class scGLMM {
    
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
        Eigen::MatrixXd ZtSigma_invX;
        Eigen::VectorXd eta;
        Eigen::VectorXd mu;
        Eigen::VectorXd u;

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
            for (int i = 0; i < n; i++) {
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
            for (int i = 0; i < n ; i++) {
                cum_ns(i) = cum_ns.head(i).sum() + ns(i);
            }
        }

        // Multiply Sigma^-1 by a vector x.
        Eigen::VectorXd Sigma_inv_x(Eigen::VectorXd x) {
            double w_sum = w.sum();
            double tau = 1 / sigma2;
            w * x + (w * (w.dot(x))) / (tau + w_sum);
        }

        // Multiply Sigma^-1 by a matrix X.
        Eigen::MatrixXd Sigma_inv_X(Eigen::MatrixXd X) {
            Eigen::MatrixXd out = Eigen::MatrixXd::Zero(X.rows(), X.cols());
            for (int j = 0; j < X.cols(); ++j) {
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

        void update_Sigma_invXs() {
            Sigma_invXs.clear();
            for (size_t i = 0; i < n; ++i) {
                Sigma_invXs.push_back(Sigma_inv_X(Xis[i]));
            }
        } 

        // Update (X^T Sigma^-1 X)^-1
        void update_XtSigma_invX_inv() {
            Eigen::MatrixXd XtSigma_invX = Eigen::MatrixXd::Zero(c, c);
            for (size_t i = 0; i < n; ++i) {
                XtSigma_invX += Xis[i].transpose() * Sigma_invXs[i];
            }
            XtSigma_invX_inv = XtSigma_invX.inverse();
        }

        // Multiply P by a vector x.
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
            return a - b;
        }

        // Multiply Omega = ZZ^T by a vector x.
        Eigen::VectorXd Omega_x(Eigen::VectorXd x) {
            Eigen::VectorXd out = Eigen::VectorXd::Zero(N);
            for (size_t i = 0; i < n; ++i) {
                double sum_i = x.segment(cum_ns(i), ns(i)).sum();
                Eigen::VectorXd tmp = Eigen::VectorXd::Constant(ns(i), sum_i);
                out.segment(cum_ns(i), ns(i)) = tmp;
            }
            return out;
        }

        // Update Z^T Sigma^-1 X
        void update_ZtSigma_invX() {
            ZtSigma_invX = Eigen::MatrixXd::Zero(n, c);
            for (size_t i; i < c; i++) {
                ZtSigma_invX.col(i) = Sigma_invXs[i].rowwise().sum(); 
            }
        }

        // Compute tr(P * Omega).
        double compute_trPOmega() {
            double tr_Sigma_invOmega = 0;
            double tr_Sigma_invRest = 0;

            for (int i = 0; i < n; ++i) {
                double wi_sum = w.segment(cum_ns(i), ns(i)).sum();
                tr_Sigma_invOmega += wi_sum / (1 + sigma2 * wi_sum);
            }

            Eigen::MatrixXd tmp;
            tmp = (XtSigma_invX_inv) * ZtSigma_invX.transpose();
            tr_Sigma_invRest = (ZtSigma_invX * tmp).trace();

            return tr_Sigma_invOmega + tr_Sigma_invRest; 
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

        void init_params() {
            GLM glm(X, y, offset, std::move(family));
            glm.fit();
            beta = glm.beta;
            mu = glm.mu;

            compute_cum_ns();
            compute_Xis();

            u = Eigen::VectorXd::Zero(n);
            update_eta();
            update_w();
            update_y_tilde();
            sigma2 = 1;
            iter = 0;
            glmm_converged = false;
            step_size = 1;
        }

        scGLMM(
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

        void fit() {
            
            while (iter < max_iter) {

                update_beta();
                update_u();
                update_eta();
                update_mu();

                update_sigma2();
                
                update_w();
                
                update_y_tilde();

                check_converge();
                if (glmm_converged) {
                    break;
                }
                update_step_size();
                iter += 1;
            }
        }
};

#endif
