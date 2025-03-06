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
#include "Family.hpp"

inline double variance(const Eigen::VectorXd& x) {
    double mean = x.mean();
    return (x.array() - mean).square().sum() / (x.size() - 1);
}

class GLMM {
    
    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        std::unique_ptr<Family> family;
        double tol = 1e-5;
        int max_iter = 10;
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
        // Relatedness matrix.
        Eigen::MatrixXd K;

        // Parameters
        Eigen::VectorXd beta;
        Eigen::VectorXd beta_prev;
        double sigma2;
        double sigma2_prev;

        // Control parameters.
        bool is_converged;
        int iter;

        void update_y_tilde() {
            Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);
            y_tilde = eta + (y - mu).cwiseProduct(mu_eta_vec);
        }

        void update_eta() {
            Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);
            eta = y_tilde + (Sigma_inv * (y_tilde - X * beta)).cwiseProduct(mu_eta_vec);
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
            Sigma += (sigma2 * K);
        }

        void update_P() {
            Sigma_invX = Sigma_inv * X;
            XtSigma_invX = X.transpose() * Sigma_invX;
            XtSigma_invX_inv = XtSigma_invX.inverse();
            P = Sigma_inv - Sigma_invX * XtSigma_invX_inv * Sigma_invX.transpose();
        }

        void update_beta() {
            beta_prev = beta;
            beta = XtSigma_invX_inv * Sigma_invX.transpose() * y_tilde;
        }

        void update_sigma2() {
            sigma2_prev = sigma2;
            double score;
            double ai;

            Eigen::VectorXd Py_tilde = P * y_tilde;
            Eigen::VectorXd KPy_tilde = K * Py_tilde;
            score = Py_tilde.dot(KPy_tilde) - (P.cwiseProduct(K)).sum();
            ai = (KPy_tilde.transpose() * KPy_tilde);

            sigma2 += score / ai;

            if (sigma2 < tol && sigma2_prev < tol) {
                sigma2 = 0.0;
            }

            while (sigma2 < 0.0) {
                sigma2 = sigma2_prev + score / ai;

                if (sigma2 < tol && sigma2_prev < tol) {
                    sigma2 = 0.0;
                } 
            }
            if (sigma2 < tol && sigma2_prev < tol) {
                sigma2 = 0.0;
            }  
        }

        void check_converge() {
            Eigen::VectorXd tol_vec_beta = Eigen::VectorXd::Constant(beta.size(), tol);
            double diff1 = ((beta - beta_prev).cwiseAbs().cwiseQuotient(
                (beta).cwiseAbs() + (beta_prev).cwiseAbs() + tol_vec_beta)).maxCoeff();
            double diff2 = std::abs(sigma2 - sigma2_prev) / (std::abs(sigma2) + std::abs(sigma2_prev) + tol);
            is_converged = (2 * std::max(diff1, diff2)) < tol;
        }

        void init_params(Eigen::Ref<Eigen::VectorXd> init_beta, const double init_sigma2) {
            beta = init_beta;
            sigma2 = init_sigma2;
            eta = X * beta;
            update_mu();
            update_w();
            update_y_tilde();
            sigma2 = std::min(0.9, variance(y_tilde));
            update_Sigma();
            update_P();
            Eigen::VectorXd Py_tilde = P * y_tilde;
            sigma2 += (Py_tilde.dot(K * Py_tilde) - P.cwiseProduct(K).sum()) * sigma2 * sigma2;
            sigma2 = std::max(0.0, sigma2 / n);
            iter = 0;
            is_converged = false;
        }

        GLMM(
            const Eigen::Ref<Eigen::MatrixXd> X_, 
            const Eigen::Ref<Eigen::VectorXd> y_, 
            std::unique_ptr<Family> family_, 
            const Eigen::MatrixXd K_, 
            const Eigen::Ref<Eigen::VectorXd> init_beta_,
            const double init_sigma2_
        ) : 
            y(y_),
            X(X_),
            K(K_),
            family(std::move(family_))
        {
            n = X.rows();
            p = X.cols();
            init_params(init_beta_, init_sigma2_);
        };

        void fit() {

            while (iter < max_iter) {
                update_Sigma();
                update_P();
                update_beta();
                update_eta();
                update_sigma2();
                update_mu();
                update_w();
                update_y_tilde();
                check_converge();
                if (is_converged) {
                    break;
                }
                iter += 1;
            }
        }
};

#endif
