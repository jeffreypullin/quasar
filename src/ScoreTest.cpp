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

#include "ScoreTest.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>

void score_test(Params& params, ModelFit& model_fit, GenoData& geno_data, PhenoData& pheno_data, CovData& cov_data) {

    std::string model = params.model;
    std::string mode = params.mode;

    Eigen::MatrixXd& G = geno_data.genotype_matrix;
    Eigen::MatrixXd& X = cov_data.data;
    Eigen::MatrixXd& Y = pheno_data.data;

    int n_samples = X.rows();
    int n_cov = X.cols();
    size_t n_snps = geno_data.n_snps;

    std::ofstream variant_file(params.out + "-quasar-" + mode + "-variant.txt");
    std::string variant_header_line = make_variant_header_line(model);
    variant_file << variant_header_line;

    std::ofstream region_file;
    if (mode == "cis") {
        region_file.open(params.out + "-quasar-cis-region.txt");
        std::string region_header_line = "feature_id\tchrom\tstart\tend\tpvalue\n";
        region_file << region_header_line;
    }

    // Iterate over features.
    for (size_t i = 0; i < pheno_data.n_pheno; ++i) {

        std::vector<double> pvals;

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
        if (model == "p_glm" || model == "nb_glm" || model == "p_glmm" || model == "nb_glmm") {
            w = model_fit.W.row(i);
        } else {
            w = Eigen::VectorXd::Ones(n_samples);
        }

        // Pre-calculate matrices used in covariate adjustment.
        Eigen::MatrixXd XtWX_inv, Xt;
        XtWX_inv = (X.transpose() * w.asDiagonal() * X).inverse();
        Xt = X.transpose();

        Eigen::VectorXd g_s(n_samples);
        Eigen::VectorXd g(n_samples);
        
        // Iterate over SNPs in the window.
        for (int k = window_start; k < window_end; ++k) {
            
            // Index into G_slice from 0 to window_n.
            int slice_ind = k - window_start; 

            std::stringstream variant_line;
            double beta, se, u, v, gtg, zscore, pval_esnp;

            // Exclude variants in the cis window when in trans mode.
            if (mode_trans && chrom == geno_data.chrom[k] && k < cis_window_end && k > cis_window_start) {
                continue;
            }

            g = G_slice.col(slice_ind); 
            g_s = g - X * (XtWX_inv * (Xt * g.cwiseProduct(w)));
            
            u = g_s.cwiseProduct(w).dot(Y.col(i));
            gtg = g_s.cwiseProduct(w).dot(g_s);
            v = gtg;

            if (model == "lmm") {
                v *= sigma2;
            }
            if (model == "p_glmm" || model == "nb_glmm") {
                v *= model_fit.tr[i];
            }
            beta = u / v;
            se = 1 / std::sqrt(v);
            if (v > 0) {
                zscore = beta / se;
                pval_esnp = 2 * pnorm(std::abs(zscore), true);
            } else {
                zscore = pval_esnp = std::numeric_limits<double>::quiet_NaN();
            }
            // MAF == 0 case.
            if (std::abs(geno_data.maf[k]) < 1e-8) {
                beta = se = zscore = pval_esnp = std::numeric_limits<double>::quiet_NaN();
            }
            if (mode == "cis") {
                pvals.push_back(pval_esnp);
            }

            variant_line << 
                pheno_data.pheno_ids[i] << "\t" <<
                geno_data.snp_id[k] << "\t" <<
                geno_data.chrom[k] << "\t" <<
                geno_data.pos[k] << "\t" <<
                geno_data.alt[k] << "\t" <<
                geno_data.ref[k] << "\t" <<
                geno_data.maf[k] << "\t" <<
                beta << "\t" << 
                se << "\t" <<
                pval_esnp;

            if (model == "p_glm") {
                variant_line << "\t" << model_fit.glm_converged[i];
            } else if (model == "nb_glm") {
                variant_line << "\t" << model_fit.glm_converged[i] <<
                    "\t" << model_fit.phi[i] <<
                    "\t" << model_fit.phi_converged[i];
            } else if (model == "p_glmm") {
                variant_line << "\t" << model_fit.glmm_converged[i];
            } else if (model == "nb_glmm") {
                variant_line << "\t" << model_fit.glmm_converged[i] << 
                    "\t" << model_fit.phi[i] <<
                    "\t" << model_fit.phi_converged[i];
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
                ACAT(pvals) << "\n";
            region_file << region_line.str();
        }
    }

    if (mode == "cis") {
        region_file.close();
    }
    variant_file.close();
}