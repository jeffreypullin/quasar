/*
MIT License

Copyright (c) 2024 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "VarRatioApprox.hpp"

#include <random>
#include "Data.hpp"
#include <Eigen/Dense>
#include <iostream>

double estimate_r(Params& params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& GRM, const Eigen::MatrixXd& G, double sigma2_e, double delta){

    int n_rand_samples = params.main_n_rand_samples;

    int n_var = G.cols();
    int n_samples = GRM.rows();

    std::random_device rd;
    std::mt19937 gen(rd());

    double sigma2_g = delta * sigma2_e;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_samples, n_samples);
    Eigen::MatrixXd Omega = sigma2_g * GRM + sigma2_e * I;
    Eigen::MatrixXd Omega_inv = Omega.inverse();

    double var_num, var_denom;
    Eigen::VectorXd r_vec(n_rand_samples);

    Eigen::MatrixXd Omega_inv_X = Omega_inv * X;
    Eigen::MatrixXd Xt_Omega_inv_X = X.transpose() * Omega_inv_X;
    Eigen::MatrixXd Xt_Omega_inv_X_inv = Xt_Omega_inv_X.inverse();

    for (int i = 0; i < n_rand_samples; ++i) {	

        std::uniform_int_distribution<> dis(0, n_var - 1);
        int rand_ind = dis(gen);

        Eigen::VectorXd g_rand = G.col(rand_ind);
        //Eigen::VectorXd g_rand_tilde = g_rand - X * (X.transpose() * X).ldlt().solve(X.transpose() * g_rand);
        Eigen::VectorXd g_rand_tilde = g_rand;
        
        if (g_rand.cwiseAbs().sum() <= 20) {
            continue;
        }

        Eigen::VectorXd Omega_inv_g = Omega_inv * g_rand_tilde;
        double a = g_rand_tilde.dot(Omega_inv_g);
    
        Eigen::VectorXd Xt_Omega_inv_g = X.transpose() * Omega_inv_g;
        Eigen::VectorXd temp = Xt_Omega_inv_X_inv * Xt_Omega_inv_g;
        double b = g_rand_tilde.dot(Omega_inv_X * temp);
        var_num = a + b;
        var_denom = g_rand_tilde.squaredNorm();

        r_vec(i) = var_num / var_denom;
    }
    return r_vec.mean();
}

double estimate_r_glmm(const Eigen::MatrixXd& X, const Eigen::MatrixXd& GRM, const Eigen::MatrixXd& G, Eigen::MatrixXd& P, Eigen::VectorXd& d){

    int n_rand_samples = 30;

    int n_var = G.cols();
    int n_samples = GRM.rows();

    std::random_device rd;
    std::mt19937 gen(rd());

    double var_num, var_denom;
    Eigen::VectorXd r_vec(n_rand_samples);

    for (int i = 0; i < n_rand_samples; ++i) {	

        std::uniform_int_distribution<> dis(0, n_var - 1);
        int rand_ind = dis(gen);

        Eigen::VectorXd g_rand = G.col(rand_ind);

        // FIXME: Why are there negative values?
        if (g_rand.cwiseAbs().sum() <= 20) {
            continue;
        }

        Eigen::VectorXd Pg = P * g_rand;
        double var_num = g_rand.dot(Pg);

        Eigen::MatrixXd D = d.asDiagonal();
        double var_denom = g_rand.dot(D * g_rand);

        r_vec(i) = var_num / var_denom;
    }
    return r_vec.mean();
}
