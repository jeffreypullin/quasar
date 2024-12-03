/*
MIT License

Copyright (c) 2024 Jeffrey Pullin 

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

#include "IntQTLMapping.hpp"
#include "GLM.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <random>

void run_qtl_mapping_lmm_int(
    Params& params,
    GenoData& geno_data, 
    CovData& cov_data, 
    PhenoData& pheno_data, 
    GRM& grm, 
    Regions& regions
) {

    Eigen::MatrixXd& Y = pheno_data.data;
    Eigen::MatrixXd& X = cov_data.data;
    
    Y = Y.transpose().eval();
    int n_samples = Y.rows();

    Eigen::MatrixXd G = geno_data.genotype_matrix.transpose().eval();

    int n_pheno = Y.cols();
    int n_cov = X.cols();

    int int_cov_ind = -1;
    for (int i = 0; i < n_cov; ++i) {
        if (cov_data.cov_ids[i] == params.int_cov) {
            int_cov_ind = i;
            break;
        }
    }
    Eigen::VectorXd x = X.col(int_cov_ind);
    Eigen::MatrixXd QtY, QtX;

    std::cout << "\nPerforming rank normalization..." << std::endl;
    rank_normalize(Y); 
    std::cout << "Rank normalization finished." << std::endl;

    std::cout << "\nPerforming eigen decomposition of GRM..." << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(grm.mat);
    Eigen::MatrixXd Q = eig.eigenvectors();
    Eigen::VectorXd lambda = eig.eigenvalues();
    std::cout << "Eigen decomposition of GRM finished." << std::endl;

    if ((lambda.array() < 0).any()) {
        std::cerr << "\nError: GRM has negative eigenvalues. Please check that the GRM is positive semi-definite." << std::endl;
        exit(1);
    }

    std::cout << "\nComputing rotated matrices..." << std::endl;
    QtY = (Q.transpose() * Y).eval();
    QtX = (Q.transpose() * X).eval();
    std::cout << "Rotated matrices computed." << std::endl;
    
    std::vector<double> delta_vec(Y.cols());
    std::vector<double> sigma2_e_vec(Y.cols());
    std::vector<Eigen::VectorXd> alpha_vec(Y.cols());

    std::cout << "\nFitting null LMMs..." << std::endl;
    for (int i = 0; i < n_pheno; ++i) {

        // Fit the LMM.
        LMM lmm(QtX, QtY.col(i), lambda);
        lmm.fit();

        // Calculate Py.
        // NB: Here we perform the residualisation of covariates 
        // with parameters estiamted under the main-effect null
        // so y = W * alpha + u + e.
        double sigma2_e = lmm.sigma2;
        double sigma2_g = lmm.delta * sigma2_e;
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_samples, n_samples);
        Eigen::MatrixXd Omega = sigma2_g * grm.mat + sigma2_e * I;
        Eigen::MatrixXd Omega_inv = Omega.inverse();

        Eigen::MatrixXd Omega_inv_X = Omega_inv * X;
        Eigen::MatrixXd Xt_Omega_inv_X = X.transpose() * Omega_inv_X;
        Eigen::MatrixXd Xt_Omega_inv_X_inv = Xt_Omega_inv_X.inverse();

        Y.col(i) = Omega_inv * Y.col(i) - Omega_inv_X * Xt_Omega_inv_X_inv * Omega_inv_X.transpose() * Y.col(i);
        
        delta_vec[i] = lmm.delta;
        sigma2_e_vec[i] = lmm.sigma2;
    }
    std::cout << "Null LMMs fitted." << std::endl;

    std::cout << "\nEstimating variance ratio correction factors..." << std::endl;
    std::vector<double> r1_vec(n_pheno);
    std::vector<double> r0_vec(n_pheno);

    for (int i = 0; i < n_pheno; ++i){
        std::cout << pheno_data.pheno_ids[i] << std::endl;

        int n_rand_samples = params.int_n_rand_samples;

        int n_var = G.cols();
        int n_samples = grm.mat.rows();

        std::random_device rd;
        std::mt19937 gen(rd());
        
        int n_cis_window_samples = n_rand_samples * params.int_n_cis_window_prop;
        Eigen::VectorXd pheno_r1_vec(n_rand_samples);
        Eigen::VectorXd pheno_r0_vec(n_rand_samples);

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_samples, n_samples);
        Eigen::MatrixXd Omega0 = delta_vec[i] * sigma2_e_vec[i] * grm.mat + sigma2_e_vec[i] * I;
        Eigen::MatrixXd Omega0_inv = Omega0.inverse();
        Eigen::MatrixXd Omega0_inv_X = Omega0_inv * X; 
        Eigen::MatrixXd Xt_Omega0_inv_X = X.transpose() * Omega0_inv_X;
        Eigen::MatrixXd Xt_Omega0_inv_X_inv = Xt_Omega0_inv_X.inverse();
        for (int j = 0; j < n_rand_samples; ++j) {	

            // Sample a genetic variant.
            int rand_ind;
            if (j < n_cis_window_samples) {
                // Sample from the cis-window.
                std::uniform_int_distribution<> dis(
                    pheno_data.window_start[i],
                    pheno_data.window_end[i]
                );
                rand_ind = dis(gen);
            } else {
                // Sample uniformly from all variants.
                std::uniform_int_distribution<> dis(0, n_var - 1);
                rand_ind = dis(gen);
            }
            // Form Z = (X g).
            Eigen::VectorXd g_rand = G.col(rand_ind);
            Eigen::MatrixXd Z(n_samples, X.cols() + 1);
            Z.leftCols(X.cols()) = X;
            Z.rightCols(1) = g_rand; 
            Eigen::MatrixXd QtZ = (Q.transpose() * Z).eval();
            Eigen::VectorXd z = (g_rand.array() * x.array()).matrix();

            LMM lmm(QtZ, QtY.col(i), lambda);
            lmm.fit();

            double delta = lmm.delta;
            double sigma2_e = lmm.sigma2;
            double sigma2_g = delta * sigma2_e;
            
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_samples, n_samples);
            Eigen::VectorXd d_vals = delta * lambda.array() + 1;
            DiagonalXd D_inv = d_vals.asDiagonal().inverse();
            Eigen::MatrixXd Omega1_inv = (1 / sigma2_e) * Q * D_inv * Q.transpose();
            Eigen::MatrixXd Omega1_inv_Z = Omega1_inv * Z; 
            Eigen::MatrixXd Zt_Omega1_inv_Z = Z.transpose() * Omega1_inv_Z;
            Eigen::MatrixXd Zt_Omega1_inv_Z_inv = Zt_Omega1_inv_Z.inverse();
            
            double var_num1, var_num0, var_denom0, var_denom1;

            //Eigen::VectorXd g_rand_tilde = g_rand - X * (X.transpose() * X).ldlt().solve(X.transpose() * g_rand);
            if (g_rand.cwiseAbs().sum() <= 20) {
                continue;
            }

            Eigen::VectorXd Omega0_inv_g = Omega0_inv * g_rand;
            double a0 = g_rand.dot(Omega0_inv_g);

            Eigen::VectorXd Xt_Omega0_inv_g = X.transpose() * Omega0_inv_g;
            Eigen::VectorXd temp0 = Xt_Omega0_inv_X_inv * Xt_Omega0_inv_g;
            double b0 = g_rand.dot(Omega0_inv_X * temp0);
            
            Eigen::VectorXd Omega1_inv_z = Omega1_inv * z;
            double a1 = z.dot(Omega1_inv_z);

            Eigen::VectorXd Zt_Omega1_inv_z = Z.transpose() * Omega1_inv_z;
            Eigen::VectorXd temp1 = Zt_Omega1_inv_Z_inv * Zt_Omega1_inv_z;
            double b1 = z.dot(Omega1_inv_Z * temp1);

            var_num0 = a0 + b0;
            var_num1 = a1 + b1;

            var_denom1 = z.squaredNorm();
            var_denom0 = g_rand.squaredNorm();

            pheno_r1_vec(j) = var_num1 / var_denom1;
            pheno_r0_vec(j) = var_num0 / var_denom0;
        }
        r1_vec[i] = pheno_r1_vec.mean();
        r0_vec[i] = pheno_r0_vec.mean();
    }
    std::cout << "Correction factor estimation finished." << std::endl;

    std::cout << "\nCalculating variant significance..." << std::endl;

    std::string region_file_path = params.output_prefix + "_cis_region.txt";
    std::string variant_file_path = params.output_prefix + "_cis_variant.txt";

    std::ofstream region_file(region_file_path);
    std::ofstream variant_file(variant_file_path);

    std::stringstream variant_header_line;
    variant_header_line << 
        "int_cov" << "\t" <<
        "feature_id" << "\t" <<
        "chrom" << "\t" <<
        "pos" << "\t" <<
        "ref" << "\t" <<
        "alt" << "\t" <<
        "beta" << "\t" << 
        "se" << "\t" <<
        "pvalue" << "\n";
     variant_file << variant_header_line.str();

    // Iterate over features.
    for (int i = 0; i < pheno_data.n_pheno; ++i) {

        // Get the window parameters for the feature.
        int window_start = pheno_data.window_start[i];
        int window_end = pheno_data.window_end[i]; 
        int window_n = pheno_data.window_n[i];

        std::stringstream block_line;
        std::vector<double> pvals;
        const Eigen::MatrixXd& G_slice = G.middleCols(window_start, window_n);

        if (params.verbose) {
            std::cout << "Processing " << window_n << 
                " SNPs for " << pheno_data.pheno_ids[i] <<
                " in region: " << geno_data.chrom[window_start] <<
                ":" << geno_data.pos[window_start] <<
                "-" <<  geno_data.pos[window_end] << std::endl;
        }

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_samples, n_samples);
        Eigen::MatrixXd Omega = delta_vec[i] * sigma2_e_vec[i] * grm.mat + sigma2_e_vec[i] * I;
        Eigen::MatrixXd Omega_inv = Omega.inverse();

        // Iterate over SNPs in the window.
        for (int k = window_start; k < window_end; ++k) {
            
            // Index into G_slice from 0 to window_n.
            int slice_ind = k - window_start; 

            std::stringstream variant_line;
            double beta, se, u, v, gtg, zscore, pval_esnp;
            double main_effect_u, main_effect_v, main_effect_beta, ztz;

            if (geno_data.var[k] > 0.0) {
                
                Eigen::VectorXd g = G_slice.col(slice_ind);
                Eigen::VectorXd z = (g.array() * x.array()).matrix();
                //Eigen::VectorXd g_tilde = g - X * (X.transpose() * X).ldlt().solve(X.transpose() * g);
                gtg = g.squaredNorm();
                ztz = z.squaredNorm();
                main_effect_u = g.dot(Y.col(i));
                main_effect_v = r0_vec[i] * gtg;
                main_effect_beta = main_effect_u / main_effect_v;
                Eigen::VectorXd y = Omega_inv * (Y.col(i) - main_effect_beta * g);
                u = z.dot(y);
                v = r1_vec[i] * ztz;
                beta = u / v;
                se = 1 / std::sqrt(v);
                if (v > 0) {
                    zscore = u / std::sqrt(v);
                    pval_esnp = 2 * pnorm(std::abs(zscore), true);
                } else{
                    zscore = 0;
                    pval_esnp = -1;
                }
                pvals.push_back(pval_esnp);

                variant_line << 
                    params.int_cov << "\t" <<
                    pheno_data.pheno_ids[i] << "\t" <<
                    geno_data.chrom[k] << "\t" <<
                    geno_data.pos[k] << "\t" <<
                    geno_data.allele1[k] << "\t" <<
                    geno_data.allele2[k] << "\t" <<
                    beta << "\t" << 
                    se << "\t" <<
                    pval_esnp << "\n";
                variant_file << variant_line.str();
            }
        }
        if (pvals.size() > 0) {
            std::stringstream region_line;

            region_line << 
                pheno_data.pheno_ids[i] << "\t" <<
                pheno_data.chrom[i] << "\t" <<
                pheno_data.start[i] << "\t" <<
                ACAT(pvals) << "\n";
            region_file << region_line.str();
        }
    }
    region_file.close();
    variant_file.close();
} 

