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

#include "QTLMapping.hpp"
#include "GLM.hpp"
#include "VarRatioApprox.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>

void run_qtl_mapping_lmm(Params& params, GenoData& geno_data, CovData& cov_data, PhenoData& pheno_data, GRM& grm){

    Eigen::MatrixXd& Y = pheno_data.data;
    Eigen::MatrixXd& X = cov_data.data;

    Y = Y.transpose().eval();
    Eigen::MatrixXd G = geno_data.genotype_matrix.transpose().eval();
    
    int n_samples = X.rows();
    int n_pheno = Y.cols();
    int n_cov = X.cols();

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

    std::cout << "\nFitting null LMMs..." << std::endl;
    for (int i = 0; i < n_pheno; ++i) {

        // Fit the LMM.
        LMM lmm(QtX, QtY.col(i), lambda);
        lmm.fit();

        // Calculate Py.
        double sigma2_e = lmm.sigma2;
        double sigma2_g = lmm.delta * sigma2_e;
        
        Y.col(i) = QtX * lmm.beta;
        
        delta_vec[i] = lmm.delta;
        sigma2_e_vec[i] = lmm.sigma2;
    }
    std::cout << "Null LMMs fitted." << std::endl;

    /*
    std::cout << "\nEstimating variance ratio correction factors..." << std::endl;
    std::vector<double> r_vec(n_pheno);
    for (int i = 0; i < n_pheno; ++i) {
        if (params.verbose) {
            std::cout << "Estimating correction factor for feature " << i + 1 << " of " << n_pheno << std::endl;
        }
        r_vec[i] = estimate_r(params, X, grm.mat, G, sigma2_e_vec[i], delta_vec[i]);
    }
    std::cout << "Correction factor estimation finished." << std::endl;
    */

    std::cout << "\nCalculating variant significance..." << std::endl;

    std::string region_file_path = params.output_prefix + "-cis-region.txt";
    std::string variant_file_path = params.output_prefix + "-cis-variant.txt";

    std::ofstream region_file(region_file_path);
    std::ofstream variant_file(variant_file_path);

    std::stringstream variant_header_line;
    variant_header_line << 
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

        double sigma2 = Y.col(i).squaredNorm() / (n_samples - n_cov);
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

        // Iterate over SNPs in the window.
        for (int k = window_start; k < window_end; ++k) {
            
            // Index into G_slice from 0 to window_n.
            int slice_ind = k - window_start; 

            std::stringstream variant_line;
            double beta, se, u, v, gtg, zscore, pval_esnp;

            if (geno_data.var[k] > 0.0) {
                
                // Calculate the score statistic.
                Eigen::VectorXd g = G_slice.col(slice_ind);
                //Eigen::VectorXd g_tilde = g - X * (X.transpose() * X).ldlt().solve(X.transpose() * g);
                Eigen::VectorXd g_tilde = g;
                u = g_tilde.dot(Y.col(i));
                gtg = g_tilde.squaredNorm();
                v = sigma2 * gtg;
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
