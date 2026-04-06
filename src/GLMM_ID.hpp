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

#ifndef GLMM_ID_H
#define GLMM_ID_H

#include <Eigen/Dense>
#include <limits>
#include <iostream>
#include "GLM.hpp"
#include "Family.hpp"

class GLMM_ID {
    
    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> offset;
        std::unique_ptr<Family> family;
        double tol = 1e-5;
        int max_iter = 50;
        int n;
        int p;

	public:

        // Diagonal elements of W.
        Eigen::VectorXd w;
        // Diagonal elements of Sigma.
        Eigen::VectorXd Sigma_diag;
        Eigen::VectorXd Sigma_diag_inv;
        Eigen::MatrixXd Sigma_invX;
        Eigen::MatrixXd XtSigma_invX;
        Eigen::MatrixXd XtSigma_invX_inv;
        Eigen::MatrixXd P;
        Eigen::VectorXd y_tilde;
        Eigen::VectorXd eta;
        Eigen::VectorXd mu;
        Eigen::VectorXd u;

        // Parameters
        Eigen::VectorXd beta;
        Eigen::VectorXd beta_prev;
        double sigma2;
        double sigma2_prev;

        Eigen::LDLT<Eigen::MatrixXd> XtSigma_invX_ldlt;

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

        void update_Sigma() {
            Sigma_diag = w.cwiseInverse().array() + sigma2;
            Sigma_diag_inv = Sigma_diag.cwiseInverse();
        }

        void update_P_components() {
            Sigma_invX = X.array().colwise() * Sigma_diag_inv.array();
            XtSigma_invX = X.transpose() * Sigma_invX;
            XtSigma_invX_ldlt.compute(XtSigma_invX);
            XtSigma_invX_inv = XtSigma_invX_ldlt.solve(Eigen::MatrixXd::Identity(p, p));
        }

        Eigen::VectorXd apply_P(Eigen::VectorXd x) {
            Eigen::VectorXd Sigma_inv_x = Sigma_diag_inv.array() * x.array();
            Eigen::VectorXd XtSigma_inv_x = X.transpose() * Sigma_inv_x;
            Eigen::VectorXd mid = XtSigma_invX_ldlt.solve(XtSigma_inv_x);
            Eigen::VectorXd Px = Sigma_inv_x - (Sigma_diag_inv.asDiagonal() * X) * mid;
            return Px;
        }

        void update_beta() {
            beta_prev = beta;
            Eigen::VectorXd Sigma_inv_y = Sigma_diag_inv.array() * y_tilde.array();
            Eigen::VectorXd XtSigma_inv_y = X.transpose() * Sigma_inv_y;
            beta = XtSigma_invX_ldlt.solve(XtSigma_inv_y);
        }
        
        void update_u() {
            Eigen::VectorXd resid = y_tilde - X * beta;
            u = sigma2 * (Sigma_diag_inv.array() * resid.array()).matrix();
        }

        void update_sigma2() {
            sigma2_prev = sigma2;
            
            double score, ai;
            Eigen::VectorXd Py_tilde = apply_P(y_tilde);
            Eigen::VectorXd KPy_tilde = Py_tilde;
            Eigen::VectorXd Sigma_inv2_diag = Sigma_diag_inv.array().square();
            Eigen::MatrixXd XtS_inv2X = X.transpose() * (Sigma_inv2_diag.asDiagonal() * X);
            double a = Sigma_diag_inv.sum();
            double b = (XtS_inv2X * XtSigma_invX_inv).trace();
            double tr = a - b;
            score = Py_tilde.dot(KPy_tilde) - tr;
            ai = (KPy_tilde.transpose() * apply_P(KPy_tilde));

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
            // Fit a GLM to initialise beta and mu.
            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLM glm(X, y, offset, std::move(poisson));
            glm.fit();
            beta = glm.beta;
            mu = glm.mu;

            u = Eigen::VectorXd::Zero(n);
            update_eta();
            update_w();
            update_y_tilde();
            sigma2 = 1;
            update_Sigma();
            update_P_components();
            iter = 0;
            glmm_converged = false;
            step_size = 1;
        }

        GLMM_ID(
            const Eigen::Ref<Eigen::MatrixXd> X_, 
            const Eigen::Ref<Eigen::VectorXd> y_, 
            const Eigen::Ref<Eigen::VectorXd> offset_,
            std::unique_ptr<Family> family_
        ) : 
            X(X_),
            y(y_),
            offset(offset_),
            family(std::move(family_))
        {
            n = X.rows();
            p = X.cols();
            init_params();
        };

        void fit() {
            
            while (iter < max_iter) {

                update_Sigma();
                update_P_components();
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
            if (std::isnan(sigma2) || (beta.hasNaN())) {
                glmm_converged = false;
            }
            P = Eigen::MatrixXd(Sigma_diag_inv.asDiagonal()) - Sigma_invX * XtSigma_invX_inv * Sigma_invX.transpose();
        }
};

#endif
