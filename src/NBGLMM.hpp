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

#ifndef NBGLMM_H
#define NBGLMM_H

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "Family.hpp"
#include "Phi.hpp"
#include "GLMM.hpp"

class NBGLMM {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> offset;
        const Eigen::Ref<Eigen::MatrixXd> grm;
        double tol = 1e-5;
        int max_iter = 10;

    public:

        Eigen::VectorXd beta;
        Eigen::VectorXd beta_prev;
        Eigen::VectorXd mu;
        Eigen::VectorXd w;
        Eigen::MatrixXd P;

        double phi;
        double phi_prev;
        double sigma2;
        double sigma2_prev;
        bool is_converged;
        bool use_apl;
        bool phi_converged;

        NBGLMM(const Eigen::Ref<Eigen::MatrixXd> X_, 
               const Eigen::Ref<Eigen::VectorXd> y_,
               const Eigen::Ref<Eigen::VectorXd> offset_,
               const Eigen::Ref<Eigen::MatrixXd> grm_) : 
               X(X_),
               y(y_),
               offset(offset_),
               grm(grm_)
        {
            beta = Eigen::VectorXd::Zero(X.cols());
            mu = y.array() + 0.1;
        };

        void check_converge () {
            Eigen::VectorXd tol_vec_beta = Eigen::VectorXd::Constant(beta.size(), tol);
            double delta1 = ((beta - beta_prev).cwiseAbs().cwiseQuotient(
                (beta).cwiseAbs() + (beta_prev).cwiseAbs() + tol_vec_beta)).maxCoeff();
            double delta2 = std::abs(sigma2 - sigma2_prev) / (std::abs(sigma2) + std::abs(sigma2_prev) + tol);
            double delta3 = std::abs(phi - phi_prev) / (std::abs(phi) + std::abs(phi_prev) + tol);
            is_converged = (2 * std::max({delta1, delta2, delta3})) < tol;
        }

        void fit() {
            
            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLMM poisson_glmm(X, y, offset, std::move(poisson), grm);
            poisson_glmm.fit();
            mu = poisson_glmm.mu;
            beta_prev = beta;
            beta = poisson_glmm.beta;
            sigma2_prev = sigma2;
            sigma2 = poisson_glmm.sigma2;

            // TODO: How to best initialise this?
            phi_prev = 0.1;
            if (use_apl) {
               phi = estimate_phi_apl(y, mu, X);
            } else {
               phi = estimate_phi_ml(y, mu, phi_converged); 
            } 

            int iter = 0;
            while (iter < max_iter) {

                std::cout << "NBGLMM iteration: " << iter << std::endl;

                auto nb = std::unique_ptr<Family>(new NegativeBinomial(phi));
                GLMM nb_glmm(X, y, offset, std::move(nb), grm);
                nb_glmm.fit();

                sigma2_prev = sigma2;
                sigma2 = nb_glmm.sigma2;

                beta_prev = beta;
                beta = nb_glmm.beta;

                mu = nb_glmm.mu;
                w = nb_glmm.w; 
                P = nb_glmm.P;

                phi_prev = phi;
                if (use_apl) {
                    phi = estimate_phi_apl(y, mu, X);
                } else {
                    phi = estimate_phi_ml(y, mu, phi_converged); 
                } 

                check_converge();
                if (is_converged) {
                    break;
                }
                iter++;
            }
        }
};

#endif