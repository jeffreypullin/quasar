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

#include "ScoreTest.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <limits>
#include <Eigen/Eigenvalues>

void score_test(Params& params, ModelFit& model_fit, GenoData& geno_data, PhenoData& pheno_data, CovData& cov_data, CellGroups& cell_groups) {

    std::string model = params.model;
    std::string mode = params.mode;
    const size_t n_int = params.interaction_covs.size();
    const size_t m_int = n_int + 1;

    Eigen::MatrixXd& G = geno_data.genotype_matrix;
    // This X is only used in the case of bulk data.
    Eigen::MatrixXd& X = cov_data.data;
    Eigen::MatrixXd& Y = pheno_data.data;

    int n_samples = X.rows();
    int n_cov = X.cols();
    size_t n_snps = geno_data.n_snps;

    bool use_cell_groups = cell_groups.n_groups > 0;
    std::vector<double> group_linear_scores;
    if (use_cell_groups && cell_groups.has_values) {
        group_linear_scores = make_group_linear_scores_values(cell_groups.group_values);
    }

    std::ofstream variant_file(params.out + "-quasar-" + mode + "-variant.txt");
    std::string variant_header_line = make_variant_header_line(params, use_cell_groups ? cell_groups.group_ids : std::vector<std::string>{}, use_cell_groups && cell_groups.has_values);
    variant_file << variant_header_line;

    std::ofstream region_file;
    if (mode == "cis") {
        region_file.open(params.out + "-quasar-cis-region.txt");
        region_file << make_region_header_line(params, use_cell_groups ? cell_groups.group_ids : std::vector<std::string>{}, use_cell_groups && cell_groups.has_values);
    }

    bool is_glmm_model = 
        model == "p_glmm" || 
        model == "p_glmm_grm" || 
        model == "p_glmm_sc" || 
        model == "nb_glmm";

    // Iterate over features.
    for (size_t i = 0; i < pheno_data.n_pheno; ++i) {

        std::vector<double> main_pvals;
        std::vector<std::vector<double>> int_pvals(n_int);
        std::vector<std::vector<double>> group_pvals_cis(use_cell_groups ? cell_groups.n_groups : 0);
        std::vector<double> group_het_pvals_cis;
        std::vector<double> group_linear_pvals_cis;
        std::vector<double> group_combined_pvals_cis;

        bool mode_trans = mode == "trans";

        int cis_window_start = 0;
        int cis_window_end = 0;
        int chrom = 0;
        if (mode == "cis" || mode_trans) {
            cis_window_start = pheno_data.window_start[i];
            cis_window_end = pheno_data.window_end[i];
            chrom = pheno_data.chrom[i];
        }

        int window_start, window_end, window_n;
        if (mode == "cis") {
            window_start = cis_window_start;
            window_end = cis_window_end; 
            window_n = pheno_data.window_n[i];
        } else {
            window_start = 0;
            window_end = n_snps;
            window_n = n_snps;
        }
       
        const Eigen::MatrixXd& G_slice = G.middleCols(window_start, window_n);

        if (params.verbose) {
            if (mode == "cis") {
                std::cout << "Processing " << window_n << 
                    " SNPs for " << pheno_data.pheno_ids[i] <<
                    " in region: " << geno_data.chrom[window_start] <<
                    ":" << geno_data.pos[window_start] <<
                    "-" <<  geno_data.pos[window_end] << std::endl;
            } else {
                std::cout << "Processing SNPs for " << pheno_data.pheno_ids[i] << std::endl;
            }
        }

        const double y_col_sqnorm = Y.col(i).squaredNorm();
        double sigma2 = y_col_sqnorm / (n_samples - n_cov);
        Eigen::VectorXd w;
        if (model == "p_glm" || 
            model == "nb_glm" || 
            model == "p_glmm" || 
            model == "p_glmm_grm" || 
            model == "p_glmm_sc" || 
            model == "nb_glmm" ||
            model == "lmm_sc") {
            w = model_fit.W.row(i);
        } else {
            w = Eigen::VectorXd::Ones(n_samples);
        }

        Eigen::MatrixXd XtWX_inv, Xt, XtX_inv, XtWZ;
        Eigen::MatrixXd ZtSigma_invX, XtSigma_invX_inv;
        Eigen::VectorXd Xty_res, XtWX_inv_Xty_res, ZtSigma_invZ_diag;
        Eigen::MatrixXd ZtASigma_invAZ, ZtAy_res;
        std::vector<Eigen::MatrixXd> ZtAkSigma_invX, XtWAkZ;
        if (params.data_type != "single-cell") {
            XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
            XtX_inv = (X.transpose() * X).inverse();
            Xt = X.transpose();
        } else {
            XtWZ = model_fit.XtWZ_vec[i];
            XtWX_inv = model_fit.XtWX_inv_vec[i];
            Xty_res = model_fit.Xty_res_vec[i];
            XtWX_inv_Xty_res = XtWX_inv * Xty_res;
            if (model == "p_glmm_sc" || model == "lmm_sc") {
                ZtSigma_invZ_diag = model_fit.ZtSigma_invZ_diag_vec[i];
                ZtSigma_invX = model_fit.ZtSigma_invX_vec[i];
                XtSigma_invX_inv = model_fit.XtSigma_invX_inv_vec[i];
                if (params.do_interaction) {
                    ZtASigma_invAZ = model_fit.ZtASigma_invAZ_vec[i];
                    ZtAkSigma_invX = model_fit.ZtAkSigma_invX_vec[i];
                    ZtAy_res = model_fit.ZtAy_res_vec[i];
                    XtWAkZ = model_fit.XtWAkZ_vec[i];
                }
            }
        }

        Eigen::VectorXd g_s(n_samples);
        Eigen::VectorXd g(n_samples);
        Eigen::VectorXd ns;
        if (params.data_type == "single-cell") {
            ns.resize(static_cast<Eigen::Index>(pheno_data.cell_counts.size()));
            for (size_t j = 0; j < pheno_data.cell_counts.size(); ++j) {
                ns(static_cast<Eigen::Index>(j)) = static_cast<double>(pheno_data.cell_counts[j]);
            }
        }

        if (model == "p_glm" && params.do_interaction) {
            std::cerr << "Error: Interaction testing is not implemented for the Poisson GLM model." << std::endl;
            std::exit(1);
        }
        
        // Iterate over SNPs in the window.
        for (int k = window_start; k < window_end; ++k) {
            
            // Index into G_slice from 0 to window_n.
            int slice_ind = k - window_start; 

            std::stringstream variant_line;
            double u, v, gtg;
            double main_beta, main_se, main_zscore, main_pval_snp;
            const double nan_val = std::numeric_limits<double>::quiet_NaN();
            std::vector<double> int_beta(n_int, nan_val);
            std::vector<double> int_se(n_int, nan_val);
            std::vector<double> int_zscore(n_int, nan_val);
            std::vector<double> int_pval_snp(n_int, nan_val);

            std::vector<double> group_betas;
            std::vector<double> group_ses;
            std::vector<double> group_pvals;
            WeightedTrendResult group_linear_result{std::numeric_limits<double>::quiet_NaN(),
                                                    std::numeric_limits<double>::quiet_NaN(),
                                                    std::numeric_limits<double>::quiet_NaN()};
            double group_acat_pvalue = std::numeric_limits<double>::quiet_NaN();
            CochranQResult het_result{std::numeric_limits<double>::quiet_NaN(),
                                      std::numeric_limits<double>::quiet_NaN(),
                                      0};
            if (use_cell_groups) {
                group_betas.assign(cell_groups.n_groups, std::numeric_limits<double>::quiet_NaN());
                group_ses.assign(cell_groups.n_groups, std::numeric_limits<double>::quiet_NaN());
                group_pvals.assign(cell_groups.n_groups, std::numeric_limits<double>::quiet_NaN());
            }

            // Exclude variants in the cis window when in trans mode.
            if (mode_trans && chrom == geno_data.chrom[k] && k < cis_window_end && k > cis_window_start) {
                continue;
            }

            g = G_slice.col(slice_ind); 
            bool model_converged = false;
            if (is_glmm_model) {
                model_converged = model_fit.glmm_converged[i];
            } else if (model == "p_glm" || model == "nb_glm") {
                model_converged = model_fit.glm_converged[i];
            } else {
                model_converged = true;
            }

            if (geno_data.mac[k] < params.min_mac || !model_converged) {

                main_beta = main_se = main_zscore = main_pval_snp = std::numeric_limits<double>::quiet_NaN();
                std::fill(int_beta.begin(), int_beta.end(), nan_val);
                std::fill(int_se.begin(), int_se.end(), nan_val);
                std::fill(int_zscore.begin(), int_zscore.end(), nan_val);
                std::fill(int_pval_snp.begin(), int_pval_snp.end(), nan_val);
                if (use_cell_groups) {
                    std::fill(group_betas.begin(), group_betas.end(), std::numeric_limits<double>::quiet_NaN());
                    std::fill(group_ses.begin(), group_ses.end(), std::numeric_limits<double>::quiet_NaN());
                    std::fill(group_pvals.begin(), group_pvals.end(), std::numeric_limits<double>::quiet_NaN());
                }

            } else if (!params.do_interaction) {

                if (params.data_type == "bulk") {
                    g_s = g - X * (XtWX_inv * (Xt * g.cwiseProduct(w)));
                    u = g_s.cwiseProduct(w).dot(Y.col(i));
                    gtg = g_s.cwiseProduct(w).dot(g_s);
                    v = gtg;
                    if (model == "lmm" || (model == "nb_glm" && params.use_quant_res)) {
                        // Mid-p quantile residuals have Var < 1; scale by residual variance.
                        v *= sigma2;
                    } else if (is_glmm_model) {
                        v *= model_fit.tr[i];
                    }
                } else if (model == "p_glmm_sc" || model == "lmm_sc") {

                    Eigen::VectorXd t = XtWZ * g;
                    u = g.dot(Y.col(i)) - t.dot(XtWX_inv_Xty_res);
                    Eigen::VectorXd t_p = ZtSigma_invX.transpose() * g;
                    v = g.cwiseProduct(ZtSigma_invZ_diag).dot(g) - t_p.dot(XtSigma_invX_inv * t_p);

                } else {

                    Eigen::VectorXd t = XtWZ * g;
                    u = g.dot(Y.col(i)) - t.dot(XtWX_inv_Xty_res);
                    gtg = g.cwiseProduct(w).dot(g) - t.dot(XtWX_inv * t);
                    v = gtg;
                    if (is_glmm_model) {
                        v *= model_fit.tr[i];
                    }
                }

                main_beta = u / v;
                main_se = 1 / std::sqrt(v);
                main_zscore = main_beta / main_se;
                main_pval_snp = 2 * pnorm(std::abs(main_zscore), false);

                if (mode == "cis") {
                    main_pvals.push_back(main_pval_snp);
                }

                if (use_cell_groups) {
                    for (size_t gi = 0; gi < cell_groups.n_groups; ++gi) {
                        if (!model_fit.glmm_converged_g_vec[i][gi]) {
                            group_betas[gi] = group_ses[gi] = group_pvals[gi] = std::numeric_limits<double>::quiet_NaN();
                            continue;
                        }

                        const Eigen::MatrixXd& XtWZ_g = model_fit.XtWZ_g_vec[i][gi];
                        const Eigen::MatrixXd& XtWX_inv_g = model_fit.XtWX_inv_g_vec[i][gi];
                        const Eigen::VectorXd& Xty_res_g = model_fit.Xty_res_g_vec[i][gi];
                        const Eigen::VectorXd& y_g_donor = model_fit.y_out_g_vec[i][gi];
                        const Eigen::VectorXd& ZtSigma_invZ_diag_g = model_fit.ZtSigma_invZ_diag_g_vec[i][gi];
                        const Eigen::MatrixXd& ZtSigma_invX_g = model_fit.ZtSigma_invX_g_vec[i][gi];
                        const Eigen::MatrixXd& XtSigma_invX_inv_g = model_fit.XtSigma_invX_inv_g_vec[i][gi];

                        Eigen::VectorXd t_g = XtWZ_g * g;
                        double raw_u_g = g.dot(y_g_donor);
                        double correction_g = t_g.dot(XtWX_inv_g * Xty_res_g);
                        double u_g = raw_u_g - correction_g;
                        Eigen::VectorXd t_p_g = ZtSigma_invX_g.transpose() * g;
                        double v_g = g.cwiseProduct(ZtSigma_invZ_diag_g).dot(g) - t_p_g.dot(XtSigma_invX_inv_g * t_p_g);

                        if (v_g <= 0.0 || std::isnan(v_g)) {
                            group_betas[gi] = group_ses[gi] = group_pvals[gi] = std::numeric_limits<double>::quiet_NaN();
                        } else {
                            double beta_g = u_g / v_g;
                            double se_g = 1.0 / std::sqrt(v_g);
                            double z_g = beta_g / se_g;
                            group_betas[gi] = beta_g;
                            group_ses[gi] = se_g;
                            group_pvals[gi] = 2 * pnorm(std::abs(z_g), false);
                        }
                    }
                    het_result = compute_cochran_q(group_betas, group_ses);
                    if (cell_groups.has_values) {
                        group_linear_result = compute_weighted_trend(group_betas, group_ses, group_linear_scores);
                        group_acat_pvalue = ACAT({
                            group_linear_result.pvalue,
                            het_result.pvalue
                        });
                    }

                    if (mode == "cis") {
                        for (size_t gi = 0; gi < cell_groups.n_groups; ++gi) {
                            group_pvals_cis[gi].push_back(group_pvals[gi]);
                        }
                        group_het_pvals_cis.push_back(het_result.pvalue);
                        if (cell_groups.has_values) {
                            group_linear_pvals_cis.push_back(group_linear_result.pvalue);
                            group_combined_pvals_cis.push_back(group_acat_pvalue);
                        }
                    }
                }

            } else {
                if (params.data_type == "bulk") {
                    Eigen::MatrixXd Z(g.size(), static_cast<Eigen::Index>(m_int));
                    Z.col(0) = g - X * (XtX_inv * (Xt * g));
                    for (size_t ik = 0; ik < n_int; ++ik) {
                        Eigen::VectorXd g_int_raw = X.col(cov_data.interaction_inds[ik]).cwiseProduct(g);
                        Z.col(static_cast<Eigen::Index>(ik + 1)) =
                            g_int_raw - X * (XtX_inv * (Xt * g_int_raw));
                    }

                    Eigen::MatrixXd ZtZ = Z.transpose() * Z;
                    Eigen::VectorXd ZtY = Z.transpose() * Y.col(i);
                    Eigen::VectorXd beta = ZtZ.ldlt().solve(ZtY);
                    Eigen::MatrixXd cov_mat = ZtZ.inverse();

                    double full_rss = (Y.col(i) - Z * beta).squaredNorm();
                    double sigma_hat = full_rss / static_cast<double>(Y.rows() - static_cast<int>(m_int));
                    cov_mat *= sigma_hat;

                    main_beta = beta(0);
                    main_se = std::sqrt(cov_mat(0, 0));
                    main_zscore = main_beta / main_se;
                    bool bad_se = (main_se < 0) || std::isnan(main_zscore);
                    for (size_t ik = 0; ik < n_int; ++ik) {
                        int_beta[ik] = beta(static_cast<Eigen::Index>(ik + 1));
                        int_se[ik] = std::sqrt(cov_mat(static_cast<Eigen::Index>(ik + 1),
                                                       static_cast<Eigen::Index>(ik + 1)));
                        int_zscore[ik] = int_beta[ik] / int_se[ik];
                        if ((int_se[ik] < 0) || std::isnan(int_zscore[ik])) {
                            bad_se = true;
                        }
                    }

                    if (bad_se) {
                        main_beta = main_se = main_zscore = main_pval_snp = nan_val;
                        std::fill(int_beta.begin(), int_beta.end(), nan_val);
                        std::fill(int_se.begin(), int_se.end(), nan_val);
                        std::fill(int_zscore.begin(), int_zscore.end(), nan_val);
                        std::fill(int_pval_snp.begin(), int_pval_snp.end(), nan_val);
                    } else {
                        main_pval_snp = 2 * pnorm(std::abs(main_zscore), false);
                        for (size_t ik = 0; ik < n_int; ++ik) {
                            int_pval_snp[ik] = 2 * pnorm(std::abs(int_zscore[ik]), false);
                        }
                    }

                } else {

                    Eigen::VectorXd U = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(m_int));
                    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(m_int),
                                                              static_cast<Eigen::Index>(m_int));
                    std::vector<Eigen::VectorXd> t(m_int);
                    std::vector<Eigen::VectorXd> t_p(m_int);
                    for (size_t a = 0; a < m_int; ++a) {
                        t[a] = XtWAkZ[a] * g;
                        t_p[a] = ZtAkSigma_invX[a].transpose() * g;
                        U(static_cast<Eigen::Index>(a)) =
                            g.dot(ZtAy_res.col(static_cast<Eigen::Index>(a))) - t[a].dot(XtWX_inv_Xty_res);
                    }
                    for (size_t a = 0; a < m_int; ++a) {
                        for (size_t b = 0; b < m_int; ++b) {
                            Eigen::VectorXd diag_ab = ZtASigma_invAZ.col(static_cast<Eigen::Index>(a * m_int + b));
                            V(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) =
                                g.cwiseProduct(diag_ab).dot(g) - t_p[a].dot(XtSigma_invX_inv * t_p[b]);
                        }
                    }

                    double main_u = U(0);
                    double main_v = V(0, 0);
                    main_beta = main_u / main_v;
                    main_se = 1 / std::sqrt(main_v);
                    main_zscore = main_beta / main_se;
                    if ((main_se < 0) || std::isnan(main_zscore)) {
                        main_beta = main_se = main_zscore = main_pval_snp = nan_val;
                    } else {
                        main_pval_snp = 2 * pnorm(std::abs(main_zscore), false);
                    }

                    if (n_int > 0 && main_v > 0.0 && !std::isnan(main_v)) {
                        Eigen::VectorXd v_cross = V.col(0).tail(static_cast<Eigen::Index>(n_int));
                        Eigen::MatrixXd V_int = V.bottomRightCorner(static_cast<Eigen::Index>(n_int),
                                                                    static_cast<Eigen::Index>(n_int));
                        Eigen::VectorXd U_int = U.tail(static_cast<Eigen::Index>(n_int));
                        Eigen::MatrixXd V_cond = V_int - v_cross * v_cross.transpose() / main_v;
                        Eigen::VectorXd U_cond = U_int - v_cross * (main_u / main_v);
                        Eigen::MatrixXd V_cond_inv = V_cond.inverse();
                        Eigen::VectorXd beta_int = V_cond_inv * U_cond;
                        for (size_t ik = 0; ik < n_int; ++ik) {
                            int_beta[ik] = beta_int(static_cast<Eigen::Index>(ik));
                            int_se[ik] = std::sqrt(V_cond_inv(static_cast<Eigen::Index>(ik),
                                                              static_cast<Eigen::Index>(ik)));
                            int_zscore[ik] = int_beta[ik] / int_se[ik];
                            if ((int_se[ik] < 0) || std::isnan(int_zscore[ik])) {
                                int_beta[ik] = int_se[ik] = int_zscore[ik] = int_pval_snp[ik] = nan_val;
                            } else {
                                int_pval_snp[ik] = 2 * pnorm(std::abs(int_zscore[ik]), false);
                            }
                        }
                    }
                }

                if (mode == "cis") {
                    main_pvals.push_back(main_pval_snp);
                    for (size_t ik = 0; ik < n_int; ++ik) {
                        int_pvals[ik].push_back(int_pval_snp[ik]);
                    }
                }
            }

            variant_line << 
                pheno_data.pheno_ids[i] << "\t" <<
                geno_data.snp_id[k] << "\t" <<
                geno_data.chrom[k] << "\t" <<
                geno_data.pos[k] << "\t" <<
                geno_data.alt[k] << "\t" <<
                geno_data.ref[k] << "\t" <<
                geno_data.maf[k] << "\t" <<
                main_beta << "\t" << 
                main_se << "\t" <<
                main_pval_snp;

            if (model == "p_glm") {
                variant_line << "\t" << model_fit.glm_converged[i];
            } else if (model == "nb_glm") {
                variant_line << "\t" << model_fit.glm_converged[i] <<
                    "\t" << model_fit.phi[i] <<
                    "\t" << model_fit.phi_converged[i];
            } else if (model == "p_glmm" || 
                       model == "p_glmm_grm" || 
                       ((model == "p_glmm_sc") & !params.do_interaction)) {
                variant_line << "\t" << model_fit.glmm_converged[i] <<
                    "\t" << model_fit.sigma2[i];
            } else if ((model == "p_glmm_sc") & params.do_interaction) {
                variant_line << "\t" << model_fit.glmm_converged[i] <<
                    "\t" << model_fit.tau0[i] <<
                    "\t" << model_fit.tau1[i] <<
                    "\t" << model_fit.tau01[i];
            } else if (model == "nb_glmm") {
                variant_line << "\t" << model_fit.glmm_converged[i] << 
                    "\t" << model_fit.phi[i] <<
                    "\t" << model_fit.phi[i] <<
                    "\t" << model_fit.phi_converged[i];
            }

            if (params.do_interaction) {
                for (size_t ik = 0; ik < n_int; ++ik) {
                    variant_line << "\t" <<
                        int_beta[ik] << "\t" <<
                        int_se[ik] << "\t" <<
                        int_pval_snp[ik];
                }
            }

            if (use_cell_groups) {
                for (size_t gi = 0; gi < cell_groups.n_groups; ++gi) {
                    if (cell_groups.has_values) {
                        variant_line << "\t" << cell_groups.group_values[gi];
                    }
                    variant_line << "\t" << group_betas[gi]
                                 << "\t" << group_ses[gi]
                                 << "\t" << group_pvals[gi];
                }
                variant_line << "\t" << het_result.q
                             << "\t" << het_result.pvalue;
                if (cell_groups.has_values) {
                    variant_line << "\t" << group_linear_result.beta
                                 << "\t" << group_linear_result.se
                                 << "\t" << group_linear_result.pvalue
                                 << "\t" << group_acat_pvalue;
                }
            }

            variant_line << "\n";
            variant_file << variant_line.str();
        }

        if (mode == "cis") {
            std::stringstream region_line;

            region_line << 
                pheno_data.pheno_ids[i] << "\t" <<
                pheno_data.chrom[i] << "\t" <<
                pheno_data.start[i] << "\t" <<
                pheno_data.end[i] << "\t" <<
                ACAT(main_pvals);

            if (params.do_interaction) {
                for (size_t ik = 0; ik < n_int; ++ik) {
                    region_line << "\t" << ACAT(int_pvals[ik]);
                }
            } 

            if (use_cell_groups) {
                for (size_t gi = 0; gi < cell_groups.n_groups; ++gi) {
                    region_line << "\t" << ACAT(group_pvals_cis[gi]);
                }
                region_line << "\t" << ACAT(group_het_pvals_cis);
                if (cell_groups.has_values) {
                    region_line << "\t" << ACAT(group_linear_pvals_cis)
                                << "\t" << ACAT(group_combined_pvals_cis);
                }
            }

            region_line << "\n";
            region_file << region_line.str();
        }
    }

    if (mode == "cis") {
        region_file.close();
    }
    variant_file.close();
}
