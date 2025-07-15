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

#ifndef NBGLM_H
#define NBGLM_H

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "Family.hpp"
#include "Phi.hpp"
#include "GLM.hpp"

class NBGLM {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> offset;
        double tol = 1e-5;
        int max_iter = 25;

    public:

        Eigen::VectorXd w;
        Eigen::VectorXd beta;
        Eigen::VectorXd y_tilde;
        Eigen::VectorXd mu;
        Eigen::MatrixXd eta;
        double phi;
        bool use_apl;
        bool phi_converged;
        bool glm_converged;

        NBGLM(const Eigen::Ref<Eigen::MatrixXd> X_, const Eigen::Ref<Eigen::VectorXd> y_, 
              const Eigen::Ref<Eigen::MatrixXd> offset_, bool use_apl_) : 
              X(X_),
              y(y_),
              offset(offset_),
              use_apl(use_apl_)
        {
            beta = Eigen::VectorXd::Zero(X.cols());
            mu = y.array() + 0.1;
        };

        double ll() {
            double theta = 1 / phi;
            double ll = 0;

            for (int i = 0; i < y.size(); i++) {
                ll += std::lgamma(y(i) + theta) 
                    - std::lgamma(theta)
                    - std::lgamma(y(i) + 1.0)
                    + theta * std::log(theta)
                    + y(i) * std::log(mu(i))
                    - (theta + y(i)) * std::log(theta + mu(i));
            }
            return ll;
        }

        void fit() {
            
            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLM poisson_glm(X, y, offset, std::move(poisson));
            poisson_glm.fit();
            mu = poisson_glm.mu;

            if (use_apl) {
               phi = estimate_phi_apl(y, mu, X, phi_converged);
            } else {
               phi = estimate_phi_ml(y, mu, phi_converged); 
            }

            double theta_0;
            // FIXME: Use residuals to calculate d1.
            double d1 = std::sqrt(2);
            double theta_delta = 1;

            double ll_m = ll();
            double ll_0 = ll_m + 2 * d1;

            glm_converged = false;
            phi_converged = false;
            for (int i = 0; i < max_iter; ++i) {
                
                auto nb = std::unique_ptr<Family>(new NegativeBinomial(phi));
                GLM nb_glm(X, y, offset, std::move(nb));
                nb_glm.fit();
                mu = nb_glm.mu;
                beta = nb_glm.beta;

                theta_0 = 1 / phi;
                if (use_apl) {
                  phi = estimate_phi_apl(y, mu, X, phi_converged);
                } else {
                  phi = estimate_phi_ml(y, mu, phi_converged); 
                }
                theta_delta = theta_0 - 1 / phi;

                ll_0 = ll_m;
                ll_m = ll();

                if ((std::abs(ll_0 - ll_m) / d1 + std::abs(theta_delta)) < tol) {
                    glm_converged = true;
                    break;
                }
            }
        }
};

#endif
