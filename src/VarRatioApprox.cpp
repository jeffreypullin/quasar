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

double estimate_r(const Eigen::MatrixXd& X, const Eigen::MatrixXd& GRM, const Eigen::MatrixXd& G, double sigma2_g, double delta){

	int n_rand_samples = 30;
	
	int n_var = G.cols();
	int n_samples = GRM.rows();

	std::random_device rd;
	std::mt19937 gen(rd());

	double sigma2_e = delta * sigma2_g;
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
		Eigen::VectorXd g_rand_tilde = g_rand - X * (X.transpose() * X).ldlt().solve(X.transpose() * g_rand);

		// FIXME: Why are there negative values?
		if (g_rand_tilde.cwiseAbs().sum() <= 20) {
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

Eigen::MatrixXd estimate_R_int(
	const CovData& cov_data, 
	const Eigen::MatrixXd& GRM, 
	const Eigen::MatrixXd& G, 
	double sigma2_g, 
	double delta, 
	const std::vector<std::string>& int_cov_ids
){

	int n_rand_samples = 30;
	
	int n_var = G.cols();
	int n_samples = GRM.rows();

	std::random_device rd;
	std::mt19937 gen(rd());

	Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(n_samples, 1 + int_cov_ids.size());
	Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n_samples, cov_data.n_cov - int_cov_ids.size());
	// Reserve first column of Z for the genotype data.
	int j = 1, k = 0;
	for (int i = 0; i < cov_data.n_cov; ++i){
		if (std::find(int_cov_ids.begin(), int_cov_ids.end(), cov_data.cov_ids[i]) != int_cov_ids.end()) {
			Z.col(j) = cov_data.data.col(i);
			++j;
		} else {
			X.col(k) = cov_data.data.col(i);
			++k;
		}
	}

	double sigma2_e = delta * sigma2_g;
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_samples, n_samples);
	Eigen::MatrixXd Omega = sigma2_g * GRM + sigma2_e * I;
	Eigen::MatrixXd Omega_inv = Omega.inverse();

	Eigen::MatrixXd var_num(Z.cols(), Z.cols());
	Eigen::MatrixXd var_denom(Z.cols(), Z.cols());
	Eigen::VectorXd r_vec(n_rand_samples);

	Eigen::MatrixXd Omega_inv_X = Omega_inv * X;
	Eigen::MatrixXd Xt_Omega_inv_X = X.transpose() * Omega_inv_X;
	Eigen::MatrixXd Xt_Omega_inv_X_inv = Xt_Omega_inv_X.inverse();

	Eigen::MatrixXd R_mean(Z.cols(), Z.cols());
	for (int i = 0; i < n_rand_samples; ++i) {	

		std::uniform_int_distribution<> dis(0, n_var - 1);
		int rand_ind = dis(gen);

		Eigen::VectorXd g_rand = G.col(rand_ind);
		// FIXME: Why are there negative values?
		if (g_rand.cwiseAbs().sum() <= 20) {
			continue;
		}

		Z.col(0) = g_rand;
		Eigen::MatrixXd Z_tilde = Z - X * (X.transpose() * X).ldlt().solve(X.transpose() * Z);

		Eigen::MatrixXd Omega_inv_Z = Omega_inv * Z_tilde;
    	Eigen::MatrixXd a = Z_tilde.transpose() * Omega_inv_Z;

		Eigen::MatrixXd Xt_Omega_inv_Z = X.transpose() * Omega_inv_Z;
    	Eigen::MatrixXd temp = Xt_Omega_inv_X_inv * Xt_Omega_inv_Z;
    	Eigen::MatrixXd b = Z_tilde.transpose() * (Omega_inv_X * temp);

		var_num = a + b;
		var_denom = Z_tilde.transpose() * Z_tilde;

		Eigen::MatrixXd R = var_num * var_denom.inverse();
		R_mean += R;
	}
	R_mean /= n_rand_samples;
	return R_mean;
}

