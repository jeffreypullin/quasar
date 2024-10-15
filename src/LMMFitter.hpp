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

#ifndef LMMFITTER_H
#define LMMFITTER_H

#include <Eigen/Dense>
#include <brent_fmin.hpp>

class LMMFitter {
	
	private:
		const Eigen::Ref<Eigen::MatrixXd> X_tilde;
		const Eigen::Ref<Eigen::VectorXd> y_tilde;
		const Eigen::Ref<Eigen::VectorXd> lambda;
		double df_resid;
		
	public:
	
		Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_inv;
		double sigma2;
		double delta;

        LMMFitter(const Eigen::Ref<Eigen::MatrixXd> X_, const Eigen::Ref<Eigen::VectorXd> y_, const Eigen::Ref<Eigen::VectorXd> lambda_) : 
			X_tilde(X_),
			y_tilde(y_),
			lambda(lambda_)
		{
            df_resid = X_tilde.rows() - X_tilde.cols();
        };

		double neg_ll_reml(double delta_arg) {
			Eigen::VectorXd d_vals = lambda.array() + 1.00 * delta_arg;
			
			Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_inv = d_vals.asDiagonal().inverse();
			
			Eigen::MatrixXd XtDX = X_tilde.transpose() * D_inv * X_tilde;
			Eigen::VectorXd XtDy = X_tilde.transpose() * D_inv * y_tilde;
			Eigen::VectorXd beta = XtDX.colPivHouseholderQr().solve(XtDy);
			
			double sigma2_ = (y_tilde.dot(D_inv * y_tilde) - XtDy.dot(beta)) / df_resid;
			
			double ll = 0.5 * (
				df_resid * std::log(sigma2_) + 
				d_vals.array().log().sum() + 
                std::log(XtDX.determinant())
			);
			
			return ll;
		};
	
		void fit_reml() {
			
			std::function<double(double)> f = [this](double x) { 
                return neg_ll_reml(x);
            };
			delta = Brent_fmin(0.00, 10000, f, 2e-5);
			Eigen::VectorXd d_vals = lambda.array() + 1.00 * delta;
			
			D_inv = d_vals.asDiagonal().inverse();
			
			Eigen::MatrixXd XtDX = X_tilde.transpose() * D_inv * X_tilde;
			Eigen::VectorXd XtDy = X_tilde.transpose() * D_inv * y_tilde;
			Eigen::VectorXd beta = XtDX.colPivHouseholderQr().solve(XtDy);
			
			sigma2 = (y_tilde.dot(D_inv * y_tilde) - XtDy.dot(beta)) / df_resid;
			
			return;
		}
};

#endif
