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

void score_test(Params& params, ModelFit& model_fit, GenoData& geno_data, PhenoData& pheno_data, CovData& cov_data) {

    std::string model = params.model;
    std::string mode = params.mode;

    Eigen::MatrixXd& G = geno_data.genotype_matrix;
    // This X is only used in the case of bulk data.
    Eigen::MatrixXd& X = cov_data.data;
    Eigen::MatrixXd& Y = pheno_data.data;

    int n_samples = X.rows();
    int n_cov = X.cols();
    size_t n_snps = geno_data.n_snps;

    std::ofstream variant_file(params.out + "-quasar-" + mode + "-variant.txt");
    std::string variant_header_line = make_variant_header_line(params);
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

        int window_start, window_end, window_n;
        int cis_window_start = pheno_data.window_start[i];
        int cis_window_end = pheno_data.window_end[i];
        int cis_window_n = pheno_data.window_n[i];
        int chrom = pheno_data.chrom[i];

        bool mode_trans = mode == "trans";

        if (mode == "cis") {
            window_start = cis_window_start;
            window_end = cis_window_end; 
            window_n = cis_window_n;
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

        double sigma2 = Y.col(i).squaredNorm() / (n_samples - n_cov);
        Eigen::VectorXd w;
        if (model == "p_glm" || 
            model == "nb_glm" || 
            model == "p_glmm" || 
            model == "p_glmm_grm" || 
            model == "p_glmm_sc" || 
            model == "nb_glmm") {
            w = model_fit.W.row(i);
        } else {
            w = Eigen::VectorXd::Ones(n_samples);
        }

        Eigen::MatrixXd XtWX_inv, Xt, XtX_inv;
        if (params.data_type != "single-cell") {
            XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
            XtX_inv = (X.transpose() * X).inverse();
            Xt = X.transpose();
        }

        Eigen::VectorXd g_s(n_samples);
        Eigen::VectorXd g(n_samples);
        Eigen::VectorXd z(n_samples);
        Eigen::VectorXd z_s(n_samples);
        
        // Iterate over SNPs in the window.
        for (int k = window_start; k < window_end; ++k) {
            
            // Index into G_slice from 0 to window_n.
            int slice_ind = k - window_start; 

            std::stringstream variant_line;
            double u, v, gtg;
            double main_beta, main_se, main_zscore, main_pval_snp;
            double int_beta, int_se, int_zscore, int_pval_snp;

            // Exclude variants in the cis window when in trans mode.
            if (mode_trans && chrom == geno_data.chrom[k] && k < cis_window_end && k > cis_window_start) {
                continue;
            }

            bool model_converged = false;
            if (is_glmm_model) {
                model_converged = model_fit.glmm_converged[i];
            } else if (model == "p_glm" || model == "nb_glm") {
                model_converged = model_fit.glm_converged[i];
            } else {
                model_converged = true;
            }

            g = G_slice.col(slice_ind); 
            if (std::abs(geno_data.maf[k]) < 1e-8 || !model_converged) {

                main_beta = main_se = main_zscore = main_pval_snp = std::numeric_limits<double>::quiet_NaN();
                int_beta = int_se = int_zscore = int_pval_snp = std::numeric_limits<double>::quiet_NaN();

            } else if (!params.do_interaction) {

                if (params.data_type == "bulk") {
                    g_s = g - X * (XtWX_inv * (Xt * g.cwiseProduct(w)));
                    u = g_s.cwiseProduct(w).dot(Y.col(i));
                    gtg = g_s.cwiseProduct(w).dot(g_s);
                } else {
                    Eigen::MatrixXd XtWZ = model_fit.XtWZ_vec[i]; 
                    Eigen::MatrixXd XtWX_inv = model_fit.XtWX_inv_vec[i]; 
                    Eigen::VectorXd Xty_res = model_fit.Xty_res_vec[i]; 
                    Eigen::VectorXd t = XtWZ * g;

                    u = g.dot(Y.col(i)) - t.dot(XtWX_inv * Xty_res);
                    gtg = g.cwiseProduct(w).dot(g) - t.dot(XtWX_inv * t);
                }
                v = gtg;

                if (model == "lmm") {
                    v *= sigma2;
                } else if (is_glmm_model) {
                    v *= model_fit.tr[i];
                }

                main_beta = u / v;
                main_se = 1 / std::sqrt(v);
                main_zscore = main_beta / main_se;
                main_pval_snp = 2 * pnorm(std::abs(main_zscore), true);

                if (mode == "cis") {
                    main_pvals.push_back(main_pval_snp);
                }

            } else {
                if (params.data_type == "single-cell") {
                    std::cerr << "Error: Interaction testing is not implemented for single-cell mode." << std::endl;
                    exit(1);
                }
                if (params.model != "lm" && params.model != "nb_glm" && params.model != "lmm") {
                    std::cerr << "Error: Interaction testing is only implemented for the LM and NB-GLM models." << std::endl;
                    exit(1);
                }

                Eigen::VectorXd x_int = X.col(cov_data.interaction_ind);
                Eigen::VectorXd g_main = g - X * (XtX_inv * (Xt * g));

                Eigen::VectorXd g_int_raw = x_int.cwiseProduct(g);
                Eigen::VectorXd g_int = g_int_raw - X * (XtX_inv * (Xt * g_int_raw));

                Eigen::MatrixXd Z(g.size(), 2);
                Z.col(0) = g_main;
                Z.col(1) = g_int;
               
                Eigen::MatrixXd ZtZ = Z.transpose() * Z;
                Eigen::VectorXd ZtY = Z.transpose() * Y.col(i);
                Eigen::VectorXd beta = ZtZ.ldlt().solve(ZtY);
                Eigen::MatrixXd ZtZ_inv = ZtZ.ldlt().solve(
                    Eigen::MatrixXd::Identity(Z.cols(), Z.cols())
                );

                Eigen::MatrixXd cov_mat = ZtZ_inv;
                double full_rss = (Y.col(i) - Z * beta).squaredNorm();
                double sigma_hat = full_rss / static_cast<double>(Y.rows() - Z.cols());
                cov_mat *= sigma_hat;

                main_beta = beta(0);   
                int_beta = beta(1);   
                main_se = std::sqrt(cov_mat(0, 0));
                int_se = std::sqrt(cov_mat(1, 1));
                main_zscore = main_beta / main_se;
                int_zscore = int_beta / int_se;

                main_pval_snp = 2 * pnorm(std::abs(main_zscore), true);
                int_pval_snp = 2 * pnorm(std::abs(int_zscore), true);

                if (mode == "cis") {
                    main_pvals.push_back(main_pval_snp);
                    int_pvals.push_back(int_pval_snp);
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
                       model == "p_glmm_sc") {
                variant_line << "\t" << model_fit.glmm_converged[i] <<
                    "\t" << model_fit.sigma2[i];
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