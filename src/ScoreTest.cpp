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
#include <vector>
#include <string>
#include <numeric>
#include <Eigen/Eigenvalues>

void score_test(Params& params, ModelFit& model_fit, GenoData& geno_data, PhenoData& pheno_data, CovData& cov_data, CellGroups& cell_groups) {

    std::string model = params.model;
    std::string mode = params.mode;

    Eigen::MatrixXd& G = geno_data.genotype_matrix;
    // This X is only used in the case of bulk data.
    Eigen::MatrixXd& X = cov_data.data;
    Eigen::MatrixXd& Y = pheno_data.data;

    int n_samples = X.rows();
    int n_cov = X.cols();
    size_t n_snps = geno_data.n_snps;

    bool use_cell_groups = cell_groups.n_groups > 0;
    std::vector<double> group_linear_scores;
    std::vector<double> group_quadratic_scores;
    if (use_cell_groups && cell_groups.has_values) {
        group_linear_scores = make_group_linear_scores_values(cell_groups.group_values);
        group_quadratic_scores = make_group_quadratic_scores(group_linear_scores);
    }

    std::ofstream variant_file(params.out + "-quasar-" + mode + "-variant.txt");
    std::string variant_header_line = make_variant_header_line(params, use_cell_groups ? cell_groups.group_ids : std::vector<std::string>{}, use_cell_groups && cell_groups.has_values);
    variant_file << variant_header_line;

    std::ofstream region_file;
    if (mode == "cis") {
        region_file.open(params.out + "-quasar-cis-region.txt");
        std::string region_header_line = "feature_id\tchrom\tstart\tend\t";
        if (params.do_interaction) {
            region_header_line += "main_acat_pvalue\tint_acat_pvalue";
        } else {
            region_header_line += "pvalue";
        }
        region_header_line += "\n";
        region_file << region_header_line;
    }

    bool is_glmm_model = 
        model == "p_glmm" || 
        model == "p_glmm_grm" || 
        model == "p_glmm_sc" || 
        model == "nb_glmm";

    // Iterate over features.
    for (size_t i = 0; i < pheno_data.n_pheno; ++i) {

        std::vector<double> main_pvals;
        std::vector<double> int_pvals;

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
        Eigen::VectorXd Xty_res, XtWX_inv_Xty_res;
        if (params.data_type != "single-cell") {
            XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
            XtX_inv = (X.transpose() * X).inverse();
            Xt = X.transpose();
        } else {
            XtWZ = model_fit.XtWZ_vec[i];
            XtWX_inv = model_fit.XtWX_inv_vec[i];
            Xty_res = model_fit.Xty_res_vec[i];
            XtWX_inv_Xty_res = XtWX_inv * Xty_res;
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
            double int_beta, int_se, int_zscore, int_pval_snp;

            std::vector<double> group_betas;
            std::vector<double> group_ses;
            std::vector<double> group_pvals;
            WeightedTrendResult group_linear_result{std::numeric_limits<double>::quiet_NaN(),
                                                    std::numeric_limits<double>::quiet_NaN(),
                                                    std::numeric_limits<double>::quiet_NaN()};
            WeightedTrendResult group_quadratic_result{std::numeric_limits<double>::quiet_NaN(),
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
                int_beta = int_se = int_zscore = int_pval_snp = std::numeric_limits<double>::quiet_NaN();
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
                } else {
                    Eigen::VectorXd t = XtWZ * g;

                    u = g.dot(Y.col(i)) - t.dot(XtWX_inv_Xty_res);
                    gtg = g.cwiseProduct(w).dot(g) - t.dot(XtWX_inv * t);
                }
                v = gtg;

                if (model == "lmm") {
                    v *= sigma2;
                } else if (is_glmm_model || model == "lmm_sc") {
                    v *= model_fit.tr[i];
                }

                main_beta = u / v;
                main_se = 1 / std::sqrt(v);
                main_zscore = main_beta / main_se;
                main_pval_snp = 2 * pnorm(std::abs(main_zscore), true);

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
                        const Eigen::VectorXd& w_g_donor = model_fit.mu_out_g_vec[i][gi];

                        Eigen::VectorXd t_g = XtWZ_g * g;
                        double raw_u_g = g.dot(y_g_donor);
                        double correction_g = t_g.dot(XtWX_inv_g * Xty_res_g);
                        double u_g = raw_u_g - correction_g;
                        double gtg_g = g.cwiseProduct(w_g_donor).dot(g) - t_g.dot(XtWX_inv_g * t_g);
                        double v_g = model_fit.tr_g_vec[i][gi] * gtg_g;

                        if (v_g <= 0.0 || std::isnan(v_g)) {
                            group_betas[gi] = group_ses[gi] = group_pvals[gi] = std::numeric_limits<double>::quiet_NaN();
                        } else {
                            double beta_g = u_g / v_g;
                            double se_g = 1.0 / std::sqrt(v_g);
                            double z_g = beta_g / se_g;
                            group_betas[gi] = beta_g;
                            group_ses[gi] = se_g;
                            group_pvals[gi] = 2 * pnorm(std::abs(z_g), true);
                        }
                    }
                    het_result = compute_cochran_q(group_betas, group_ses);
                    if (cell_groups.has_values) {
                        group_linear_result = compute_weighted_trend(group_betas, group_ses, group_linear_scores);
                        group_quadratic_result = compute_weighted_trend(group_betas, group_ses, group_quadratic_scores);
                        group_acat_pvalue = ACAT({
                            group_linear_result.pvalue,
                            group_quadratic_result.pvalue,
                            het_result.pvalue
                        });
                    }
                }

            } else {
                Eigen::MatrixXd ZtZ = Eigen::MatrixXd::Zero(2, 2);
                Eigen::VectorXd ZtY = Eigen::VectorXd::Zero(2);
                Eigen::MatrixXd Z;

                if (params.data_type == "bulk") {
                    Eigen::VectorXd x_int = X.col(cov_data.interaction_ind);
                    Eigen::VectorXd g_main = g - X * (XtX_inv * (Xt * g));

                    Eigen::VectorXd g_int_raw = x_int.cwiseProduct(g);
                    Eigen::VectorXd g_int = g_int_raw - X * (XtX_inv * (Xt * g_int_raw));
                    Z.resize(g.size(), 2);
                    Z.col(0) = g_main;
                    Z.col(1) = g_int;
                   
                    ZtZ = Z.transpose() * Z;
                    ZtY = Z.transpose() * Y.col(i);
                    Eigen::VectorXd beta = ZtZ.ldlt().solve(ZtY);
                    Eigen::MatrixXd cov_mat = ZtZ.inverse();
         
                    double full_rss, sigma_hat;
                    full_rss = (Y.col(i) - Z * beta).squaredNorm();
                    sigma_hat = full_rss / static_cast<double>(Y.rows() - 2);
                    cov_mat *= sigma_hat;

                    main_beta = beta(0);
                    int_beta = beta(1); 
                    main_se = std::sqrt(cov_mat(0, 0));
                    int_se = std::sqrt(cov_mat(1, 1));
                    main_zscore = main_beta / main_se;
                    int_zscore = int_beta / int_se;

                    if (main_se < 0 || int_se < 0) {
                        main_beta = main_se = main_zscore = main_pval_snp = std::numeric_limits<double>::quiet_NaN();
                        int_beta = int_se = int_zscore = int_pval_snp = std::numeric_limits<double>::quiet_NaN();
                    } else {
                        main_pval_snp = 2 * pnorm(std::abs(main_zscore), true);
                        int_pval_snp = 2 * pnorm(std::abs(int_zscore), true);
                    }
                    
                } else {

                    Eigen::VectorXd ZtDy_res = model_fit.ZtDy_res_vec[i];
                    Eigen::VectorXd Zty_res = model_fit.Zty_res_vec[i];
                    Eigen::MatrixXd XtWDZ = model_fit.XtWDZ_vec[i];
                    Eigen::VectorXd d_out = model_fit.d_out_vec[i];
                    Eigen::VectorXd dw_out = model_fit.dw_out_vec[i];
                    Eigen::VectorXd dwd_out = model_fit.dwd_out_vec[i];

                    Eigen::VectorXd t = XtWZ * g;
                    Eigen::VectorXd Dt = XtWDZ * g;

                    // Compute main effect.
                    double main_u, main_v, gtg;
                    main_u = g.dot(Y.col(i)) - t.dot(XtWX_inv_Xty_res);
                    gtg = g.cwiseProduct(w).dot(g) - t.dot(XtWX_inv * t);
                    main_v = gtg * model_fit.tr[i];
                    main_beta = main_u / main_v;
                    main_se = 1 / std::sqrt(main_v);
                    main_zscore = main_beta / main_se;
                    if ((main_se < 0) | std::isnan(main_zscore)) {
                        main_beta = main_se = main_zscore = main_pval_snp = std::numeric_limits<double>::quiet_NaN();
                    } else {
                        main_pval_snp = 2 * pnorm(std::abs(main_zscore), true);
                    }

                    // Compute interaction effect.
                    double int_u, int_v, ztz;
                    int n = g.size();
                    int c = Xty_res.size();
                    Eigen::VectorXd Cty(c + 1);
                    Cty.head(c) = Xty_res;
                    Cty(c) = g.dot(Zty_res);
                    Eigen::MatrixXd CtWC(c + 1, c + 1);
                    CtWC.topLeftCorner(c, c) = XtWX_inv.inverse();
                    CtWC.topRightCorner(c, 1) = t;
                    CtWC.bottomLeftCorner(1, c) = t.transpose();
                    CtWC(c, c) = g.cwiseProduct(w).dot(g);
                    Eigen::MatrixXd CtWC_inv = CtWC.inverse();
                    Eigen::MatrixXd CtWDZ(c + 1, n);
                    CtWDZ.topRows(c) = XtWDZ;
                    CtWDZ.bottomRows(1) = g.cwiseProduct(dw_out).transpose();
                    Eigen::VectorXd CDt = CtWDZ * g;
                    
                    int_u = g.dot(ZtDy_res) - CDt.dot(CtWC_inv * Cty);
                    ztz = g.cwiseProduct(dwd_out).dot(g) - CDt.dot(CtWC_inv * CDt);
               
                    int_v = ztz * model_fit.tr_int[i];
                    int_beta = int_u / int_v;
                    int_se = 1 / std::sqrt(int_v);
                    int_zscore = int_beta / int_se;
                    if ((int_se < 0) | std::isnan(main_zscore)) {
                        int_beta = int_se = int_zscore = int_pval_snp = std::numeric_limits<double>::quiet_NaN();
                    } else {
                        int_pval_snp = 2 * pnorm(std::abs(int_zscore), true);
                    }

                    if (mode == "cis") {
                        main_pvals.push_back(main_pval_snp);
                        int_pvals.push_back(int_pval_snp);
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
                variant_line << "\t" << model_fit.glmm_converged[i];
            } else if (model == "nb_glmm") {
                variant_line << "\t" << model_fit.glmm_converged[i] << 
                    "\t" << model_fit.phi[i] <<
                    "\t" << model_fit.phi[i] <<
                    "\t" << model_fit.phi_converged[i];
            }

            if (params.do_interaction) {
                variant_line << "\t" <<
                    int_beta << "\t" << 
                    int_se << "\t" <<
                    int_pval_snp;
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
                                 << "\t" << group_quadratic_result.beta
                                 << "\t" << group_quadratic_result.se
                                 << "\t" << group_quadratic_result.pvalue
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
                region_line << "\t" << ACAT(int_pvals);
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
