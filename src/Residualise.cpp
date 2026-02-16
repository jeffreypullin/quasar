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

#include "Residualise.hpp"
#include "QTLMappingUtils.hpp"
#include "Data.hpp"
#include "ModelFit.hpp"
#include "GLM.hpp"
#include "LM.hpp"
#include "LMM.hpp"
#include "NBGLM.hpp"
#include "GLMM_GRM.hpp"
#include "GLMM_ID.hpp"
#include "GLMM_SC.hpp"
#include "NBGLMM.hpp"
#include "Phi.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>

void residualise(Params& params, ModelFit& model_fit, CovData& cov_data, PhenoData& pheno_data, GRM& grm) {

    Eigen::MatrixXd& Y = pheno_data.data;
    if (params.data_type == "single-cell") {
        Y = Eigen::MatrixXd::Zero(pheno_data.n_samples, pheno_data.n_pheno);
    }
    
    // We need to expand 'bulk' covariate data into single-cell covariate data.
    if (params.data_type == "single-cell" && cov_data.cov_data_type == "bulk") {
        cov_data.expand_cov_data(pheno_data.cell_counts);
    }

    Eigen::MatrixXd& X = (params.data_type == "single-cell") ? cov_data.sc_data : cov_data.data;

    int n_pheno = pheno_data.n_pheno;

    if (params.model == "lmm" || params.model == "lm") {
        std::cout << "\nPerforming rank normalization..." << std::endl;
        rank_normalize(Y); 
        std::cout << "Rank normalization finished." << std::endl;
    }

    // Compute offset for count-distribution models.
    Eigen::VectorXd offset;
    if (params.data_type == "single-cell") {
        offset = pheno_data.offset;
    } else {
        offset = Y.rowwise().sum().array().log();
    }
    
    Eigen::MatrixXd W(n_pheno, pheno_data.n_samples);

    std::vector<double> tr;
    std::vector<double> phi;

    std::vector<bool> glm_converged;
    std::vector<bool> phi_converged;
    std::vector<bool> glmm_converged;
    std::vector<Eigen::MatrixXd> XtWX_inv_vec;
    std::vector<Eigen::VectorXd> Xty_res_vec;
    std::vector<Eigen::MatrixXd> XtWZ_vec;

    Eigen::VectorXd ns(pheno_data.cell_counts.size());
    for (size_t idx = 0; idx < pheno_data.cell_counts.size(); ++idx) {
        ns(idx) = static_cast<double>(pheno_data.cell_counts[idx]);
    }

    if (params.model == "p_glmm_sc") {
        
        std::cout <<"\nFitting single-cell Poisson GLMMs..." << std::endl; 
        for (int i = 0; i < n_pheno; ++i) {

            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLMM_SC p_glmm_sc(X, pheno_data.sc_data.col(i), offset, std::move(poisson), ns);
            p_glmm_sc.fit();

            Y.col(i) = p_glmm_sc.y_out;
            W.row(i) = p_glmm_sc.mu_out;
            XtWX_inv_vec.push_back(p_glmm_sc.XtWX_inv);
            Xty_res_vec.push_back(p_glmm_sc.Xty_res);
            XtWZ_vec.push_back(p_glmm_sc.XtWZ);

            tr.push_back(p_glmm_sc.r_approx);
            glmm_converged.push_back(p_glmm_sc.glmm_converged);
        }
        std::cout << "Null single-cell Poisson GLMMs fitted." << std::endl;

    } else if (params.model == "lmm") {

        std::cout << "\nPerforming eigen decomposition of GRM..." << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(grm.mat);
        Eigen::MatrixXd Q = eig.eigenvectors();
        Eigen::VectorXd lambda = eig.eigenvalues();
        std::cout << "Eigen decomposition of GRM finished." << std::endl;

        if ((lambda.array() < 0).any()) {
            std::cerr << "\nError: GRM has negative eigenvalues. Please check that the GRM is positive semi-definite." << std::endl;
            exit(1);
        }

        Eigen::MatrixXd QtY, QtX;
        std::cout << "\nComputing rotated matrices..." << std::endl;
        QtY = (Q.transpose() * Y).eval();
        QtX = (Q.transpose() * X).eval();
        std::cout << "Rotated matrices computed." << std::endl;

        std::cout << "\nFitting null LMMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            LMM lmm(QtX, QtY.col(i), lambda);
            lmm.fit();

            Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_inv = lmm.D_inv;
            Y.col(i) = (Q * D_inv * (QtY.col(i) - QtX * lmm.beta)) / lmm.sigma2;
        
        }
        std::cout << "Null LMMs fitted." << std::endl;

    } else if (params.model == "lm") {

        std::cout << "\nFitting null LMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            LM lm(X, Y.col(i));
            lm.fit();

            Y.col(i) = (Y.col(i) - X * lm.beta) / std::sqrt(lm.s);
        }
        std::cout << "Null LMs fitted." << std::endl;

    } else if (params.model == "p_glmm" || params.model == "p_glmm_grm") {
        
        std::cout << "\nFitting null Poisson GLMMs with GRM..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLMM_GRM p_glmm_grm(X, Y.col(i), offset, std::move(poisson), grm.mat);
            p_glmm_grm.fit();
            
            Y.col(i) = (Y.col(i).array() - p_glmm_grm.mu.array()) / p_glmm_grm.mu.array();
            W.row(i) = p_glmm_grm.mu.array();

            Eigen::MatrixXd P = p_glmm_grm.P;
            Eigen::VectorXd w = p_glmm_grm.mu;
            tr.push_back(compute_r_approx(P, w, X));
            glmm_converged.push_back(p_glmm_grm.glmm_converged);
        }
        std::cout << "Null Poisson GLMMs fitted." << std::endl;

    } else if (params.model == "p_glmm_id") {
        
        std::cout << "\nFitting null Poisson identity GLMMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {
            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLMM_ID p_glmm_id(X, Y.col(i), offset, std::move(poisson));
            p_glmm_id.fit();
            
            Y.col(i) = (Y.col(i).array() - p_glmm_id.mu.array()) / p_glmm_id.mu.array();
            W.row(i) = p_glmm_id.mu.array();

            Eigen::MatrixXd P = p_glmm_id.P;
            Eigen::VectorXd w = p_glmm_id.mu;
            tr.push_back(compute_r_approx(P, w, X));
            glmm_converged.push_back(p_glmm_id.glmm_converged);
        }
        std::cout << "Null Poisson GLMMs fitted." << std::endl;

    } else if (params.model == "nb_glm") {

        std::cout << "\nFitting null NB-GLMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            bool use_apl = params.use_apl;
            NBGLM nb_glm(X, Y.col(i), offset, use_apl);
            nb_glm.fit();

            Y.col(i) = (Y.col(i).array() - (X * nb_glm.beta + offset).array().exp()) / nb_glm.mu.array();
            W.row(i) = nb_glm.mu.array() / (1 + nb_glm.phi * nb_glm.mu.array());

            phi.push_back(nb_glm.phi);
            phi_converged.push_back(nb_glm.phi_converged);
            glm_converged.push_back(nb_glm.glm_converged);
        }
        std::cout << "Null NB-GLMs fitted." << std::endl;

    } else if (params.model == "p_glm") {

        std::cout << "\nFitting null Poisson-GLMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLM p_glm(X, Y.col(i), offset, std::move(poisson));
            p_glm.fit();

            Y.col(i) = (Y.col(i).array() - (X * p_glm.beta + offset).array().exp()) / p_glm.mu.array();
            W.row(i) = p_glm.mu.array();
            glm_converged.push_back(p_glm.glm_converged);
        }
        std::cout << "Null Poisson GLMs fitted." << std::endl;

    } else if (params.model == "nb_glmm") {

        std::cout << "\nFitting null NB GLMMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            NBGLMM nb_glmm(X, Y.col(i), offset, grm.mat);
            nb_glmm.fit();

            Y.col(i) = (Y.col(i).array() - nb_glmm.mu.array()) / nb_glmm.mu.array();
            W.row(i) = nb_glmm.mu.array() / (1 + nb_glmm.phi * nb_glmm.mu.array());

            Eigen::MatrixXd P = nb_glmm.P;
            Eigen::VectorXd w = nb_glmm.mu;
            tr.push_back(compute_r_approx(P, w, X));
            phi.push_back(nb_glmm.phi);
            phi_converged.push_back(nb_glmm.phi_converged);
            glmm_converged.push_back(nb_glmm.glmm_converged);
        }
        std::cout << "Null NB GLMMs fitted." << std::endl;
    
    }

    model_fit.W = W;
    model_fit.phi = phi;
    model_fit.tr = tr;
    model_fit.XtWX_inv_vec = XtWX_inv_vec;
    model_fit.Xty_res_vec = Xty_res_vec;
    model_fit.XtWZ_vec = XtWZ_vec;

    model_fit.phi_converged = phi_converged;
    model_fit.glm_converged = glm_converged;
    model_fit.glmm_converged = glmm_converged;

    pheno_data.data = Y;
}

double compute_r_approx(
    const Eigen::MatrixXd& P,
    const Eigen::VectorXd& w,
    const Eigen::MatrixXd& X
) {
    Eigen::MatrixXd W = w.asDiagonal();
    Eigen::MatrixXd WX = W * X;

    double tr_P = P.trace();
    double tr_W = w.sum();

    Eigen::MatrixXd XtWX_inv = (X.transpose() * W * X).inverse();
    Eigen::MatrixXd XtW2X = (X.transpose() * (w.array() * w.array()).matrix().asDiagonal() * X);
    Eigen::MatrixXd XtW3X = (X.transpose() * (w.array() * w.array() * w.array()).matrix().asDiagonal() * X);

    double tr_WPw = tr_W - (XtW2X * XtWX_inv).trace();
    double a = tr_P / tr_WPw;
    
    double tr_PWPw = (P * W).trace() - (((P * WX) * XtWX_inv) * WX.transpose()).trace();
    double b = 2 * tr_PWPw / pow(tr_WPw, 2);

    double tmp = (XtW2X * XtWX_inv * XtW2X * XtWX_inv).trace();
    double tr_WPwWPw = (w.array() * w.array()).sum() - 2 * (XtW3X * XtWX_inv).trace() + tmp;
    double c = ((2 * tr_WPwWPw) * tr_P) / pow(tr_WPw, 3);

    return a - b + c;
}
