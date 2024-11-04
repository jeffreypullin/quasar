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

double var(const Eigen::VectorXd& x) {
    double mean = x.mean();
    return (x.array() - mean).square().sum() / (x.size() - 1);
}

// Currently speclialised to Poisson GLMMs.
class GLMM {
	
	public:

        // Data.
        int n;
        int p;
        Eigen::MatrixXd y;
        Eigen::MatrixXd X;
        Eigen::MatrixXd lib_size;
        Eigen::MatrixXd Z;

        // Diagonal elements of D.
        Eigen::VectorXd d;
        Eigen::MatrixXd Sigma;
        Eigen::MatrixXd Sigma_inv;
        Eigen::MatrixXd Sigma_invX;
        Eigen::MatrixXd XtSigma_invX;
        Eigen::MatrixXd XtSigma_invX_inv;
        Eigen::MatrixXd P;
        Eigen::VectorXd y_tilde;
        Eigen::VectorXd eta;
        Eigen::VectorXd mu;
        // FIXME: Can we just use two individual matrices?
        std::vector<Eigen::MatrixXd> Ks;

        // Parameters
        Eigen::VectorXd beta;
        Eigen::VectorXd beta_prev;
        Eigen::VectorXd tau;
        Eigen::VectorXd tau_prev;

        // Control parameters.
        double tol = 1e-5;
        int max_iter = 10;
        bool is_converged;
        double step_size;
        int iter;

        void update_y_tilde() {
            y_tilde = eta - lib_size + (y - mu).cwiseQuotient(d);
        }

        void update_eta() {
            eta = y_tilde - (Sigma_inv * (y_tilde - X * beta)).cwiseQuotient(d);
            for (int i = 0; i < n; i++) {
                if (eta(i) < -30 || eta(i) > 30) {
                    eta(i) = std::numeric_limits<double>::epsilon();
                }
            }
        }

        void update_d() {
            mu = eta.array().exp();
            d = mu;
        }

        void update_Sigma() {
            Sigma = d.asDiagonal().inverse();
            for (int i = 0; i < tau.size(); i++) {
                Sigma += tau(i) * Ks[i];
            }
            Sigma_inv = Sigma.inverse();
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

        void update_tau() {
            tau_prev = tau;
            Eigen::VectorXd scores(tau.size());
            Eigen::MatrixXd AI(tau.size(), tau.size());

            Eigen::VectorXd Py_tilde = P * y_tilde;
            for(int i = 0; i < Ks.size(); i++) {

                Eigen::VectorXd KiPy_tilde = Ks[i] * Py_tilde;
                scores(i) = Py_tilde.dot(KiPy_tilde) - (P.cwiseProduct(Ks[i])).sum();
                for(int j = 0; j <= i; j++) {
                    Eigen::VectorXd KjPy_tilde = Ks[j] * Py_tilde;
                    // FIXME: Check this.
                    AI(i,j) = (KiPy_tilde.transpose() * KjPy_tilde)(0,0);
                    if(i != j) {
                        AI(j,i) = AI(i,j);
                    }
                }
            }
            for (int i = 0; i < tau.size(); i++) {
                tau(i) += step_size * (AI.ldlt().solve(scores))(i);
            }

            for (int i = 0; i < tau.size(); i++) {
                if (tau(i) < tol && tau_prev(i) < tol) {
                    tau(i) = 0.0;
                }
            }

            double ss = step_size;
            while((tau.array() < 0.0).any()) {
                ss *= 0.5;
                tau = ss * AI.ldlt().solve(scores) + tau_prev;
                
                for (int i = 0; i < tau.size(); i++) {
                    if (tau(i) < tol && tau_prev(i) < tol) {
                        tau(i) = 0.0;
                    }
                }
            }
            for (int i = 0; i < tau.size(); i++) {
                if (tau(i) < tol) {
                    tau(i) = 0.0;
                }
            }
        }

        void check_converge() {
            Eigen::VectorXd tol_vec_beta = Eigen::VectorXd::Constant(beta.size(), tol);
            Eigen::VectorXd tol_vec_tau = Eigen::VectorXd::Constant(tau.size(), tol);
            double diff1 = ((beta - beta_prev).cwiseAbs().cwiseQuotient(
                (beta).cwiseAbs() + (beta_prev).cwiseAbs() + tol_vec_beta)).maxCoeff();
            double diff2 = ((tau - tau_prev).cwiseAbs().cwiseQuotient(
                (tau).cwiseAbs() + (tau_prev).cwiseAbs() + tol_vec_tau)).maxCoeff();
            is_converged = (2 * std::max(diff1, diff2)) < tol;
        }

        void update_step_size() {

            if ((iter + 1) % 10 == 0) {
                step_size = 0.9 * step_size;
            }
        }

        void init_params(Eigen::Ref<Eigen::VectorXd> init_beta, Eigen::Ref<Eigen::VectorXd> init_tau) {
            beta = init_beta;
            tau = init_tau;
            eta = X * beta + lib_size;
            update_d();
            update_y_tilde();
            for (int i = 0; i < tau.size(); i++) {
                tau(i) = std::min(0.9, var(y_tilde) / tau.size());
            }
            update_Sigma();
            update_P();
            Eigen::VectorXd Py_tilde = P * y_tilde;
            for(int i = 0; i < tau.size(); i++) {
                // FIXME: Check this.
                tau(i) += (Py_tilde.dot(Ks[i] * Py_tilde) - P.cwiseProduct(Ks[i]).sum()) * tau(i) * tau(i);
                tau(i) = std::max(0.0, tau(i) / n);
            }
            iter = 0;
            step_size = 1;
            is_converged = false;
        }

        GLMM(
            const Eigen::Ref<Eigen::MatrixXd> X_, 
            const Eigen::Ref<Eigen::VectorXd> y_, 
            const std::vector<Eigen::MatrixXd> Ks_, 
            const Eigen::Ref<Eigen::VectorXd> lib_size_,
            const Eigen::Ref<Eigen::VectorXd> init_beta_,
            const Eigen::Ref<Eigen::VectorXd> init_tau_
        ) : 
            y(y_),
            X(X_),
            Ks(Ks_),
            lib_size(lib_size_)
        {
            n = X.rows();
            p = X.cols();
            init_params(init_beta_, init_tau_);
        };

        void fit() {
            int iter = 0;
            while (iter < max_iter) {
                update_Sigma();
                update_P();
                update_beta();
                update_eta();
                update_tau();
                update_d();
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
