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

#ifndef LMM_H
#define LMM_H

#include <Eigen/Dense>
#include <brent_fmin.hpp>

// Speclised to use the FastLMM algorithm.
class LMM {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X_tilde;
        const Eigen::Ref<Eigen::VectorXd> y_tilde;
        const Eigen::Ref<Eigen::VectorXd> lambda;
        double df_resid;

    public:
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_inv;
        double sigma2;
        double delta;
        Eigen::VectorXd beta;

        LMM(const Eigen::Ref<Eigen::MatrixXd> X_, const Eigen::Ref<Eigen::VectorXd> y_, const Eigen::Ref<Eigen::VectorXd> lambda_) : 
            X_tilde(X_),
            y_tilde(y_),
            lambda(lambda_)
        {
            df_resid = X_tilde.rows() - X_tilde.cols();
        };

        double neg_ll_reml(double delta_arg) {
            Eigen::VectorXd d_vals = delta_arg * lambda.array() + 1;
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

        void fit() {

            std::function<double(double)> f = [this](double x) { 
                return neg_ll_reml(x);
            };
            delta = Brent_fmin(0.00, 10000, f, 2e-5);
            Eigen::VectorXd d_vals = delta * lambda.array() + 1.00;
            
            D_inv = d_vals.asDiagonal().inverse();
            
            Eigen::MatrixXd XtDX = X_tilde.transpose() * D_inv * X_tilde;
            Eigen::VectorXd XtDy = X_tilde.transpose() * D_inv * y_tilde;
            beta = XtDX.colPivHouseholderQr().solve(XtDy);

            sigma2 = (y_tilde.dot(D_inv * y_tilde) - XtDy.dot(beta)) / df_resid;

            return;
        }
};

#endif
