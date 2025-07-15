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
        bool glmm_converged;
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
            glmm_converged = (2 * std::max({delta1, delta2, delta3})) < tol;
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
               phi = estimate_phi_apl(y, mu, X, phi_converged);
            } else {
               phi = estimate_phi_ml(y, mu, phi_converged); 
            } 

            int iter = 0;
            while (iter < max_iter) {

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
                    phi = estimate_phi_apl(y, mu, X, phi_converged);
                } else {
                    phi = estimate_phi_ml(y, mu, phi_converged); 
                } 

                check_converge();
                if (glmm_converged) {
                    break;
                }
                iter++;
            }
        }
};

#endif