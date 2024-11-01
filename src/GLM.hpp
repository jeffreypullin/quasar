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

// This is currently spcialised to Poisson GLMs.
class GLM {
	
    private:
		const Eigen::Ref<Eigen::MatrixXd> X;
		const Eigen::Ref<Eigen::VectorXd> y;
        int n;
        int p;
        double tol = 1e-5;

	public:
	
        Eigen::VectorXd w;
        Eigen::VectorXd beta;
        Eigen::VectorXd y_tilde;
        Eigen::VectorXd mu;
        Eigen::MatrixXd eta;

        GLM(const Eigen::Ref<Eigen::MatrixXd> X_, const Eigen::Ref<Eigen::VectorXd> y_) : 
			X(X_),
			y(y_)
		{
            n = y.size();
            p = X.cols();
            beta = Eigen::VectorXd::Zero(p);
            mu = y.array() + 0.1;
            w = mu.array().inverse();
        };

        Eigen::VectorXd weighted_least_squares() {
            Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
            Eigen::MatrixXd XtWX_inv = XtWX.colPivHouseholderQr().inverse();
            Eigen::VectorXd XtWz = X.transpose() * w.asDiagonal() * y_tilde;
            return XtWX_inv * XtWz;
        }

        double neg_ll() {
            Eigen::VectorXd ll = -mu.array() + y.array() * log(mu.array());
            return ll.sum();
        }

        // Fit the GLM using IRLS.
		void fit() {

            double delta = 1;
            int iter = 1;
            double ll = 0;
            
            while (delta > tol) {

                y_tilde = log(mu.array()) + (y.array() - mu.array()) / mu.array();
                w = mu.array().inverse();
    
                beta = weighted_least_squares();
                eta = X * beta;
                mu = eta.array().exp();
    
                double ll_old = ll;
                ll = neg_ll();
                delta = std::abs(ll - ll_old);

                iter++;
            }
		}
};

#endif
