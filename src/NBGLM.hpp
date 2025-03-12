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
        double tol = 1e-5;
        int max_iter = 10;

    public:

        Eigen::VectorXd w;
        Eigen::VectorXd beta;
        Eigen::VectorXd y_tilde;
        Eigen::VectorXd mu;
        Eigen::MatrixXd eta;
        double phi;

        NBGLM(const Eigen::Ref<Eigen::MatrixXd> X_, const Eigen::Ref<Eigen::VectorXd> y_) : 
              X(X_),
              y(y_)
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
            GLM poisson_glm(X, y, std::move(poisson));
            poisson_glm.fit();
            mu = poisson_glm.mu;
            
            phi = estimate_phi_ml(y, mu);

            double theta_0;
            // FIXME: Use residuals to calculate d1.
            double d1 = std::sqrt(2);
            double theta_delta = 1;

            double ll_m = ll();
            double ll_0 = ll_m + 2 * d1;

            int iter = 0;
            while ((std::abs(ll_0 - ll_m) / d1 + std::abs(theta_delta)) > tol) {
                
                std::cout << "NBGLM iteration: " << iter << std::endl;
                auto nb = std::unique_ptr<Family>(new NegativeBinomial(phi));
                GLM nb_glm(X, y, std::move(nb));
                nb_glm.fit();
                mu = nb_glm.mu;
                beta = nb_glm.beta;

                theta_0 = 1 / phi;
                phi = estimate_phi_ml(y, mu);
                theta_delta = theta_0 - 1 / phi;
                std::cout << "Theta: " << 1 / phi << std::endl;

                ll_0 = ll_m;
                ll_m = ll();

                iter++;

                if (iter > max_iter) {
                    break;
                }
            }
        }
};

#endif
