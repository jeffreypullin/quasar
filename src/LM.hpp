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

#ifndef LM_H
#define LM_H

#include <Eigen/Dense>

class LM {

    private:
        const Eigen::Ref<Eigen::MatrixXd> X;
        const Eigen::Ref<Eigen::VectorXd> y;
        int n;
        int p;

    public:

        Eigen::VectorXd beta;
        double s;
        Eigen::MatrixXd V;

        LM(const Eigen::Ref<Eigen::MatrixXd> X_, const Eigen::Ref<Eigen::VectorXd> y_) : 
            X(X_),
            y(y_)
        {
            n = y.size();
            p = X.cols();
        };

        void fit() {

            Eigen::MatrixXd XtX = X.transpose()  * X;
            Eigen::MatrixXd XtX_inv = XtX.colPivHouseholderQr().inverse();
            beta = XtX_inv * X.transpose() * y;

            s = (y - X * beta).squaredNorm() / (n - p);
            V = s * XtX_inv;
        }
};

#endif