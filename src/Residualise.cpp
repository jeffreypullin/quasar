/*
  Copyright (C) 2024-26 Jeffrey Pullin

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
#include "LMM_SC.hpp"
#include "NBGLM.hpp"
#include "GLMM_GRM.hpp"
#include "GLMM_SC.hpp"
#include "GLMM_SC_INT.hpp"
#include "NBGLMM.hpp"
#include "Phi.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <limits>
#include <boost/math/special_functions/beta.hpp>

void residualise(Params& params, ModelFit& model_fit, CovData& cov_data, PhenoData& pheno_data, GRM& grm, CellGroups& cell_groups) {

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

    // Use precomputed / file-provided offset when available; otherwise log library size (bulk).
    Eigen::VectorXd offset;
    if (pheno_data.offset.size() > 0) {
        offset = pheno_data.offset;
    } else {
        offset = Y.rowwise().sum().array().log();
    }
    
    Eigen::MatrixXd W(n_pheno, pheno_data.n_samples);

    std::vector<double> tr;
    std::vector<double> tr_int;
    std::vector<double> phi;
    std::vector<double> sigma2;
    std::vector<double> tau0;
    std::vector<double> tau1;
    std::vector<double> tau2;

    std::vector<bool> glm_converged;
    std::vector<bool> phi_converged;
    std::vector<bool> glmm_converged;
    std::vector<bool> lmm_converged;
    std::vector<Eigen::MatrixXd> XtWX_inv_vec;
    std::vector<Eigen::VectorXd> Xty_res_vec;
    std::vector<Eigen::MatrixXd> XtWZ_vec;
    std::vector<Eigen::VectorXd> ZtSigma_invZ_diag_vec;
    std::vector<Eigen::MatrixXd> ZtSigma_invX_vec;
    std::vector<Eigen::MatrixXd> XtSigma_invX_inv_vec;
    std::vector<Eigen::VectorXd> ZtDSigma_invDZ_diag_vec;
    std::vector<Eigen::VectorXd> ZtDSigma_invZ_diag_vec;
    std::vector<Eigen::MatrixXd> ZtDSigma_invX_vec;
    std::vector<Eigen::VectorXd> ZtDy_res_vec;
    std::vector<Eigen::MatrixXd> XtWDZ_vec;
    std::vector<Eigen::VectorXd> Zty_res_vec;
    std::vector<Eigen::VectorXd> d_out_vec;
    std::vector<Eigen::VectorXd> dw_out_vec;
    std::vector<Eigen::VectorXd> dwd_out_vec;
    std::vector<std::vector<Eigen::MatrixXd>> XtWX_inv_g_vec;
    std::vector<std::vector<Eigen::VectorXd>> Xty_res_g_vec;
    std::vector<std::vector<Eigen::MatrixXd>> XtWZ_g_vec;
    std::vector<std::vector<Eigen::VectorXd>> y_out_g_vec;
    std::vector<std::vector<Eigen::VectorXd>> ZtSigma_invZ_diag_g_vec;
    std::vector<std::vector<Eigen::MatrixXd>> ZtSigma_invX_g_vec;
    std::vector<std::vector<Eigen::MatrixXd>> XtSigma_invX_inv_g_vec;
    std::vector<std::vector<double>> sigma2_g_vec;
    std::vector<std::vector<bool>> glmm_converged_g_vec;
    Eigen::VectorXd ns(pheno_data.cell_counts.size());
    for (size_t idx = 0; idx < pheno_data.cell_counts.size(); ++idx) {
        ns(idx) = static_cast<double>(pheno_data.cell_counts[idx]);
    }

    if ((params.model == "lmm_sc") & !params.do_interaction) {
        
        std::cout <<"\nFitting single-cell LMMs..." << std::endl;

        const int n_donors = ns.size();
        const int p = X.cols();
        Eigen::VectorXd cum_ns = Eigen::VectorXd::Zero(n_donors);
        for (int i = 1; i < n_donors; ++i) {
            cum_ns(i) = cum_ns(i - 1) + ns(i - 1);
        }
        Eigen::MatrixXd XtX = X.transpose() * X;
        Eigen::MatrixXd X_tilde = Eigen::MatrixXd::Zero(n_donors, p);
        Eigen::MatrixXd ZtX = Eigen::MatrixXd::Zero(n_donors, p);
        for (int j = 0; j < p; ++j) {
            for (int i = 0; i < n_donors; ++i) {
                const double sum = X.col(j).segment(cum_ns(i), ns(i)).sum();
                ZtX(i, j) = sum;
                X_tilde(i, j) = sum / std::sqrt(ns(i));
            }
        }

        for (int i = 0; i < n_pheno; ++i) {

            Eigen::VectorXd y = pheno_data.sc_data.col(i);
            rank_normalize_vec(y);
            LMM_SC lmm_sc(X, y, ns, XtX, X_tilde, ZtX);
            lmm_sc.fit();

            Y.col(i) = lmm_sc.y_out;
            W.row(i) = lmm_sc.mu_out;
            XtWX_inv_vec.push_back(lmm_sc.XtWX_inv);
            Xty_res_vec.push_back(lmm_sc.Xty_res);
            XtWZ_vec.push_back(lmm_sc.XtWZ);
            ZtSigma_invZ_diag_vec.push_back(lmm_sc.ZtSigma_invZ_diag);
            ZtSigma_invX_vec.push_back(lmm_sc.ZtSigma_invX);
            XtSigma_invX_inv_vec.push_back(lmm_sc.XtSigma_invX_inv);
            sigma2.push_back(lmm_sc.sigma2);
        }
        std::cout << "Null single-cell LMMs fitted." << std::endl;

    } else if ((params.model == "p_glmm_sc") & !params.do_interaction) {

        const size_t N_sc = static_cast<size_t>(X.rows());
        const size_t c_sc = static_cast<size_t>(X.cols());
        const size_t n_donors_sc = static_cast<size_t>(ns.size());
        std::vector<size_t> cum_ns_sc;
        if (cell_groups.n_groups > 0) {
            if (cell_groups.cell_to_group.size() != N_sc) {
                std::cerr << "Error: cell_to_group size does not match number of cells when fitting per-group GLMMs." << std::endl;
                std::exit(1);
            }
            cum_ns_sc.assign(n_donors_sc, 0);
            for (size_t d = 1; d < n_donors_sc; ++d) {
                cum_ns_sc[d] = cum_ns_sc[d - 1] + static_cast<size_t>(ns(static_cast<Eigen::Index>(d - 1)));
            }
        }

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
            ZtSigma_invZ_diag_vec.push_back(p_glmm_sc.ZtSigma_invZ_diag);
            ZtSigma_invX_vec.push_back(p_glmm_sc.ZtSigma_invX);
            XtSigma_invX_inv_vec.push_back(p_glmm_sc.XtSigma_invX_inv);

            if (cell_groups.n_groups > 0) {
                
                const double nan_val = std::numeric_limits<double>::quiet_NaN();
                std::vector<Eigen::MatrixXd> XtWX_inv_g(cell_groups.n_groups, Eigen::MatrixXd::Zero(c_sc, c_sc));
                std::vector<Eigen::VectorXd> Xty_res_g(cell_groups.n_groups, Eigen::VectorXd::Zero(c_sc));
                std::vector<Eigen::MatrixXd> XtWZ_g(cell_groups.n_groups, Eigen::MatrixXd::Zero(c_sc, n_donors_sc));
                std::vector<Eigen::VectorXd> y_out_g(cell_groups.n_groups, Eigen::VectorXd::Zero(n_donors_sc));
                std::vector<Eigen::VectorXd> ZtSigma_invZ_diag_g(cell_groups.n_groups, Eigen::VectorXd::Zero(n_donors_sc));
                std::vector<Eigen::MatrixXd> ZtSigma_invX_g(cell_groups.n_groups, Eigen::MatrixXd::Zero(n_donors_sc, c_sc));
                std::vector<Eigen::MatrixXd> XtSigma_invX_inv_g(cell_groups.n_groups, Eigen::MatrixXd::Zero(c_sc, c_sc));
                std::vector<double> sigma2_g(cell_groups.n_groups, nan_val);
                std::vector<bool> converged_g(cell_groups.n_groups, false);

                Eigen::VectorXd y_col = pheno_data.sc_data.col(i);

                for (size_t gi = 0; gi < cell_groups.n_groups; ++gi) {
                    std::vector<size_t> cell_inds_g;
                    cell_inds_g.reserve(N_sc);
                    std::vector<size_t> donor_keep_g;
                    std::vector<int> ns_g_kept;
                    donor_keep_g.reserve(n_donors_sc);
                    ns_g_kept.reserve(n_donors_sc);

                    for (size_t d = 0; d < n_donors_sc; ++d) {
                        const size_t start = cum_ns_sc[d];
                        const size_t end = start + static_cast<size_t>(ns(static_cast<Eigen::Index>(d)));
                        int n_cells_d_g = 0;
                        for (size_t k = start; k < end; ++k) {
                            if (cell_groups.cell_to_group[k] == static_cast<int>(gi)) {
                                cell_inds_g.push_back(k);
                                n_cells_d_g++;
                            }
                        }
                        if (n_cells_d_g > 0) {
                            donor_keep_g.push_back(d);
                            ns_g_kept.push_back(n_cells_d_g);
                        }
                    }

                    const size_t cells_g = cell_inds_g.size();
                    const size_t n_donors_g = donor_keep_g.size();

                    Eigen::MatrixXd X_g(cells_g, c_sc);
                    Eigen::VectorXd y_g(cells_g);
                    Eigen::VectorXd offset_g(cells_g);
                    for (size_t k = 0; k < cells_g; ++k) {
                        const Eigen::Index src = static_cast<Eigen::Index>(cell_inds_g[k]);
                        const Eigen::Index dst = static_cast<Eigen::Index>(k);
                        X_g.row(dst) = X.row(src);
                        y_g(dst) = y_col(src);
                        offset_g(dst) = offset(src);
                    }
                    Eigen::VectorXd ns_g(static_cast<Eigen::Index>(n_donors_g));
                    for (size_t d = 0; d < n_donors_g; ++d) {
                        ns_g(static_cast<Eigen::Index>(d)) = static_cast<double>(ns_g_kept[d]);
                    }

                    auto poisson_g = std::unique_ptr<Family>(new Poisson());
                    GLMM_SC fit_g(X_g, y_g, offset_g, std::move(poisson_g), ns_g);
                    fit_g.fit();

                    if (!fit_g.glmm_converged) {
                        XtWX_inv_g[gi].setConstant(nan_val);
                        Xty_res_g[gi].setConstant(nan_val);
                        XtWZ_g[gi].setConstant(nan_val);
                        y_out_g[gi].setConstant(nan_val);
                        ZtSigma_invZ_diag_g[gi].setConstant(nan_val);
                        ZtSigma_invX_g[gi].setConstant(nan_val);
                        XtSigma_invX_inv_g[gi].setConstant(nan_val);
                        converged_g[gi] = false;
                        continue;
                    }

                    XtWX_inv_g[gi] = fit_g.XtWX_inv;
                    Xty_res_g[gi] = fit_g.Xty_res;
                    XtSigma_invX_inv_g[gi] = fit_g.XtSigma_invX_inv;
                    for (size_t d = 0; d < n_donors_g; ++d) {
                        const Eigen::Index full_d = static_cast<Eigen::Index>(donor_keep_g[d]);
                        const Eigen::Index kept_d = static_cast<Eigen::Index>(d);
                        y_out_g[gi](full_d) = fit_g.y_out(kept_d);
                        ZtSigma_invZ_diag_g[gi](full_d) = fit_g.ZtSigma_invZ_diag(kept_d);
                        XtWZ_g[gi].col(full_d) = fit_g.XtWZ.col(kept_d);
                        ZtSigma_invX_g[gi].row(full_d) = fit_g.ZtSigma_invX.row(kept_d);
                    }
                    sigma2_g[gi] = fit_g.sigma2;
                    converged_g[gi] = fit_g.glmm_converged;
                }

                XtWX_inv_g_vec.push_back(XtWX_inv_g);
                Xty_res_g_vec.push_back(Xty_res_g);
                XtWZ_g_vec.push_back(XtWZ_g);
                y_out_g_vec.push_back(y_out_g);
                ZtSigma_invZ_diag_g_vec.push_back(ZtSigma_invZ_diag_g);
                ZtSigma_invX_g_vec.push_back(ZtSigma_invX_g);
                XtSigma_invX_inv_g_vec.push_back(XtSigma_invX_inv_g);
                sigma2_g_vec.push_back(sigma2_g);
                glmm_converged_g_vec.push_back(converged_g);
            }

            glmm_converged.push_back(p_glmm_sc.glmm_converged);
            sigma2.push_back(p_glmm_sc.sigma2);
        }
        std::cout << "Null single-cell Poisson GLMMs fitted." << std::endl;

    } else if ((params.model == "p_glmm_sc") & params.do_interaction) {

        std::cout <<"\nFitting null random slope single-cell Poisson GLMMs..." << std::endl; 
        for (int i = 0; i < n_pheno; ++i) {

            Eigen::VectorXd x = cov_data.sc_data.col(cov_data.interaction_ind);
            Eigen::VectorXd y = pheno_data.sc_data.col(i);
            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLMM_SC_INT p_glmm(X, y, x, offset, std::move(poisson), ns);
            p_glmm.fit();

            Y.col(i) = p_glmm.y_out;
            W.row(i) = p_glmm.mu_out;
            XtWX_inv_vec.push_back(p_glmm.XtWX_inv);
            Xty_res_vec.push_back(p_glmm.Xty_res);
            XtWZ_vec.push_back(p_glmm.XtWZ);
            ZtSigma_invZ_diag_vec.push_back(p_glmm.ZtSigma_invZ_diag);
            ZtSigma_invX_vec.push_back(p_glmm.ZtSigma_invX);
            XtSigma_invX_inv_vec.push_back(p_glmm.XtSigma_invX_inv);
            ZtDSigma_invDZ_diag_vec.push_back(p_glmm.ZtDSigma_invDZ_diag);
            ZtDSigma_invZ_diag_vec.push_back(p_glmm.ZtDSigma_invZ_diag);
            ZtDSigma_invX_vec.push_back(p_glmm.ZtDSigma_invX);
            ZtDy_res_vec.push_back(p_glmm.ZtDy_res);
            XtWDZ_vec.push_back(p_glmm.XtWDZ);
            Zty_res_vec.push_back(p_glmm.Zty_res);
            d_out_vec.push_back(p_glmm.d_out);
            dw_out_vec.push_back(p_glmm.dw_out);
            dwd_out_vec.push_back(p_glmm.dwd_out);
            glmm_converged.push_back(p_glmm.glmm_converged);
            tau0.push_back(p_glmm.tau(0));
            tau1.push_back(p_glmm.tau(1));
            tau2.push_back(p_glmm.tau(2));
        }
        std::cout << "Null random-slope single-cell Poisson GLMMs fitted." << std::endl;

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
            if (params.do_interaction) {
                Y.col(i) = Q * D_inv * (QtY.col(i) - QtX * lmm.beta) / std::sqrt(lmm.sigma2);
            } else {
                Y.col(i) = (Q * D_inv * (QtY.col(i) - QtX * lmm.beta)) / lmm.sigma2;
            }
        
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
            sigma2.push_back(p_glmm_grm.sigma2);
        }
        std::cout << "Null Poisson GLMMs fitted." << std::endl;

    } else if (params.model == "nb_glm") {

        std::cout << "\nFitting null NB-GLMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            bool use_apl = params.use_apl;
            NBGLM nb_glm(X, Y.col(i), offset, use_apl);
            nb_glm.fit();
            
            if (!params.use_quant_res) {
                Y.col(i) = (Y.col(i).array() - (X * nb_glm.beta + offset).array().exp()) / nb_glm.mu.array();
                W.row(i) = nb_glm.mu.array() / (1 + nb_glm.phi * nb_glm.mu.array());
            } else {
                // Mid-p quantile residual: Phi^{-1}((a+b)/2). Residual scale is
                // estimated in the score test via sigma2.
                const double size = 1.0 / nb_glm.phi;
                for (size_t j = 0; j < pheno_data.n_samples; ++j) {
                    const double y_raw = Y(j, i);
                    const double y = std::round(y_raw);
                    const double mu = nb_glm.mu(j);
                    const double p = size / (mu + size);
                    double a = 0.0;
                    if (y > 0.0) {
                        a = boost::math::ibeta(size, std::max(y, 1.0), p);
                    }
                    double b = boost::math::ibeta(size, y + 1.0, p);
                    a = std::max(0.0, std::min(a, 1.0));
                    b = std::max(a,   std::min(b, 1.0));
                    double u = 0.5 * (a + b);
                    u = std::max(1e-12, std::min(u, 1.0 - 1e-12));
                    Y(j, i) = qnorm(u, true);
                }
                W.row(i).setOnes();
            }

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

            bool use_apl = params.use_apl;
            NBGLMM nb_glmm(X, Y.col(i), offset, grm.mat, use_apl);
            nb_glmm.fit();

            Y.col(i) = (Y.col(i).array() - nb_glmm.mu.array()) / nb_glmm.mu.array();
            W.row(i) = nb_glmm.mu.array() / (1 + nb_glmm.phi * nb_glmm.mu.array());

            Eigen::MatrixXd P = nb_glmm.P;
            Eigen::VectorXd w = nb_glmm.mu;
            tr.push_back(compute_r_approx(P, w, X));
            phi.push_back(nb_glmm.phi);
            phi_converged.push_back(nb_glmm.phi_converged);
            glmm_converged.push_back(nb_glmm.glmm_converged);
            sigma2.push_back(nb_glmm.sigma2);
        }
        std::cout << "Null NB GLMMs fitted." << std::endl;
    
    }

    model_fit.W = W;
    model_fit.phi = phi;
    model_fit.tr = tr;
    model_fit.tr_int = tr_int;
    model_fit.sigma2 = sigma2;
    model_fit.tau0 = tau0;
    model_fit.tau1 = tau1;
    model_fit.tau2 = tau2;
    model_fit.XtWX_inv_vec = XtWX_inv_vec;
    model_fit.Xty_res_vec = Xty_res_vec;
    model_fit.XtWZ_vec = XtWZ_vec;
    model_fit.ZtSigma_invZ_diag_vec = ZtSigma_invZ_diag_vec;
    model_fit.ZtSigma_invX_vec = ZtSigma_invX_vec;
    model_fit.XtSigma_invX_inv_vec = XtSigma_invX_inv_vec;
    model_fit.ZtDSigma_invDZ_diag_vec = ZtDSigma_invDZ_diag_vec;
    model_fit.ZtDSigma_invZ_diag_vec = ZtDSigma_invZ_diag_vec;
    model_fit.ZtDSigma_invX_vec = ZtDSigma_invX_vec;
    model_fit.ZtDy_res_vec = ZtDy_res_vec;
    model_fit.XtWDZ_vec = XtWDZ_vec;
    model_fit.Zty_res_vec = Zty_res_vec;
    model_fit.d_out_vec = d_out_vec;
    model_fit.dw_out_vec = dw_out_vec;
    model_fit.dwd_out_vec = dwd_out_vec;

    model_fit.phi_converged = phi_converged;
    model_fit.glm_converged = glm_converged;
    model_fit.glmm_converged = glmm_converged;

    if (cell_groups.n_groups > 0) {
        model_fit.n_groups = cell_groups.n_groups;
        model_fit.group_ids = cell_groups.group_ids;
        model_fit.XtWX_inv_g_vec = XtWX_inv_g_vec;
        model_fit.Xty_res_g_vec = Xty_res_g_vec;
        model_fit.XtWZ_g_vec = XtWZ_g_vec;
        model_fit.y_out_g_vec = y_out_g_vec;
        model_fit.ZtSigma_invZ_diag_g_vec = ZtSigma_invZ_diag_g_vec;
        model_fit.ZtSigma_invX_g_vec = ZtSigma_invX_g_vec;
        model_fit.XtSigma_invX_inv_g_vec = XtSigma_invX_inv_g_vec;
        model_fit.sigma2_g_vec = sigma2_g_vec;
        model_fit.glmm_converged_g_vec = glmm_converged_g_vec;
    }

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
