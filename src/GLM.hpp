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

#ifndef GLM_H
#define GLM_H

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "Family.hpp"

class GLM {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        const Eigen::Ref<Eigen::VectorXd> offset;
        std::unique_ptr<Family> family;
        int max_iter = 25;
        double tol = 1e-8;

    public:

        Eigen::VectorXd w;
        Eigen::VectorXd beta;
        Eigen::VectorXd y_tilde;
        Eigen::VectorXd mu;
        Eigen::MatrixXd eta;
        bool glm_converged;

        GLM(const Eigen::Ref<Eigen::MatrixXd> X_, 
            const Eigen::Ref<Eigen::VectorXd> y_,
            const Eigen::Ref<Eigen::VectorXd> offset_,
            std::unique_ptr<Family> family_) :
            X(X_),
            y(y_),
            offset(offset_),
            family(std::move(family_))
        {
            beta = Eigen::VectorXd::Zero(X.cols());
            mu = family->init(y);
            w = compute_w(mu);
        };

        Eigen::VectorXd compute_w(Eigen::VectorXd mu) {
            Eigen::VectorXd v_vec = family->var(mu);
            Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);

            return (v_vec.array() * mu_eta_vec.array().square()).inverse().matrix();
        }

        Eigen::VectorXd weighted_least_squares() {
            Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
            Eigen::MatrixXd XtWX_inv = XtWX.colPivHouseholderQr().inverse();
            Eigen::VectorXd XtWy = X.transpose() * w.asDiagonal() * y_tilde;
            return XtWX_inv * XtWy;
        }

        void fit() {

            double dev = 1; 
            double dev_old, delta;
            
            glm_converged = false;
            for (int i = 0; i < max_iter; ++i) {

                eta = family->link(mu);
                Eigen::VectorXd mu_eta_vec = family->mu_eta(mu);
                y_tilde = ((eta.array() - offset.array()) + mu_eta_vec.array() * (y - mu).array()).matrix();
                w = compute_w(mu);
    
                beta = weighted_least_squares();
                eta = X * beta;
                mu = family->invlink(eta + offset);

                dev_old = dev;
                dev = family->dev(y, mu);
                delta = std::abs(dev - dev_old) / (0.1 + std::abs(dev));

                if (delta < tol) {
                    glm_converged = true;
                    break;
                }
            }
        }
};

#endif
