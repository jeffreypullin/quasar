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

#ifndef GLMM_H
#define GLMM_H

#include <Eigen/Dense>
#include <limits>
#include <iostream>
#include "GLM.hpp"
#include "Family.hpp"

class GLMM {
    
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
        Eigen::MatrixXd Sigma;
        Eigen::MatrixXd Sigma_inv;
        Eigen::MatrixXd Sigma_invX;
        Eigen::MatrixXd XtSigma_invX;
        Eigen::MatrixXd XtSigma_invX_inv;
        Eigen::MatrixXd P;
        Eigen::VectorXd y_tilde;
        Eigen::VectorXd eta;
        Eigen::VectorXd mu;
        Eigen::VectorXd u;

        // Relatedness matrix.
        Eigen::MatrixXd K;

        // Parameters
        Eigen::VectorXd beta;
        Eigen::VectorXd beta_prev;
        double sigma2;
        double sigma2_prev;

        Eigen::LDLT<Eigen::MatrixXd> Sigma_ldlt;
        Eigen::LDLT<Eigen::MatrixXd> XtSigmaX_ldlt;

        // Control parameters.
        bool is_converged;
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
            Sigma = w.asDiagonal().inverse();
            Sigma += sigma2 * K;
            Sigma_ldlt.compute(Sigma);
        }

        void update_P() {
            Sigma_invX = Sigma_ldlt.solve(X);
            XtSigma_invX = X.transpose() * Sigma_invX;
            XtSigmaX_ldlt.compute(XtSigma_invX);
            
            P = Sigma_ldlt.solve(Eigen::MatrixXd::Identity(n, n)) - 
                Sigma_invX * XtSigmaX_ldlt.solve(Sigma_invX.transpose());
        }

        void update_beta() {
            beta_prev = beta;
            beta = XtSigmaX_ldlt.solve(X.transpose() * Sigma_ldlt.solve(y_tilde));
        }
        
        void update_u() {
            u = sigma2 * K * Sigma_ldlt.solve(y_tilde - (X * beta));
        }

        void update_sigma2() {
            sigma2_prev = sigma2;
            
            double score, ai;
            Eigen::VectorXd Py_tilde = P * y_tilde;
            Eigen::VectorXd KPy_tilde = K * Py_tilde;
            score = Py_tilde.dot(KPy_tilde) - (P * K).trace();
            ai = (KPy_tilde.transpose() * P * KPy_tilde);

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
            is_converged = (2 * std::max(diff1, diff2)) < tol;
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
            update_P();
            iter = 0;
            is_converged = false;
            step_size = 1;
        }

        GLMM(
            const Eigen::Ref<Eigen::MatrixXd> X_, 
            const Eigen::Ref<Eigen::VectorXd> y_, 
            const Eigen::Ref<Eigen::VectorXd> offset_,
            std::unique_ptr<Family> family_, 
            const Eigen::MatrixXd K_
        ) : 
            X(X_),
            y(y_),
            offset(offset_),
            family(std::move(family_)),
            K(K_)
        {
            n = X.rows();
            p = X.cols();
            init_params();
        };

        void fit() {
            
            while (iter < max_iter) {

                std::cout << "GLMM iteration: " << iter << std::endl;
                update_Sigma();
                update_P();

                update_beta();
                update_u();
                update_eta();
                update_mu();

                update_sigma2();
                
                update_w();
                
                update_y_tilde();

                check_converge();
                if (is_converged) {
                    break;
                }
                update_step_size();
                iter += 1;
            }
        }
};

#endif
