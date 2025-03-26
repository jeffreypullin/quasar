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
#include "LM.hpp"
#include "LMM.hpp"
#include "NBGLM.hpp"
#include "NBGLMM.hpp"
#include "Phi.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>

void run_qtl_mapping(Params& params, GenoData& geno_data, CovData& cov_data, PhenoData& pheno_data, GRM& grm){

    Eigen::MatrixXd& Y = pheno_data.data;
    Eigen::MatrixXd& X = cov_data.data;
    Eigen::MatrixXd& G = geno_data.genotype_matrix;

    int n_samples = X.rows();
    int n_pheno = Y.cols();
    int n_cov = X.cols();

    if (params.model == "lmm" || params.model == "lm") {
        std::cout << "\nPerforming rank normalization..." << std::endl;
        rank_normalize(Y); 
        std::cout << "Rank normalization finished." << std::endl;
    }

    // Step 1 ------------------------------------------------------------
    std::cout << "\nResidualising features...." << std::endl;
    Eigen::MatrixXd W_mat(n_pheno, n_samples);
    if (params.model == "lmm") {

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
        for (int i = 0; i < Y.cols(); ++i) {

            LMM lmm(QtX, QtY.col(i), lambda);
            lmm.fit();

            double sigma2 = lmm.sigma2;
            Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_inv = lmm.D_inv;
            Y.col(i) = (Q.transpose() * D_inv * (QtY.col(i) - QtX * lmm.beta)) / std::sqrt(sigma2);
        
        }
        std::cout << "Null LMMs fitted." << std::endl;

    } else if (params.model == "lm") {

        std::cout << "\nFitting null LMs..." << std::endl;
        for (int i = 0; i < Y.cols(); ++i) {

            LM lm(X, Y.col(i));
            lm.fit();

            Y.col(i) = Y.col(i) - X * lm.beta;
        }
        std::cout << "Null LMs fitted." << std::endl;

    } else if (params.model == "glmm") {
        
        std::cout << "\nFitting null GLMMs..." << std::endl;
        for (int i = 0; i < Y.cols(); ++i) {

            std::cout << "Processing feature " << i + 1 << " of " << Y.cols() << std::endl;
            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLMM poisson_glmm(X, Y.col(i), std::move(poisson), grm.mat);
            poisson_glmm.fit();
            Y.col(i) = Y.col(i).array() - (X * poisson_glmm.beta).array().exp();
            W_mat.row(i) = poisson_glmm.w;
        }
        std::cout << "Null GLMMs fitted." << std::endl;

    } else if (params.model == "glm") {

        std::cout << "\nFitting null NB-GLMs..." << std::endl;
        for (int i = 0; i < Y.cols(); ++i) {

            bool use_apl = params.use_apl;
            NBGLM nb_glm(X, Y.col(i), use_apl);
            nb_glm.fit();

            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLM poisson_glm(X, Y.col(i), std::move(poisson)); 
            poisson_glm.fit();
        }
        std::cout << "Null NB-GLMs fitted." << std::endl;
    }

    // Step 2 ------------------------------------------------------------
    std::cout << "\nCalculating variant significance..." << std::endl;

    std::string region_file_path = params.output_prefix + "-cis-region.txt";
    std::string variant_file_path = params.output_prefix + "-cis-variant.txt";

    std::ofstream region_file(region_file_path);
    std::ofstream variant_file(variant_file_path);

    std::stringstream variant_header_line;
    variant_header_line << "feature_id\tchrom\tpos\tref\talt\tbeta\tse\tpvalue\n";
    variant_file << variant_header_line.str();

    std::stringstream region_header_line;
    region_header_line << "feature_id\tchrom\tpos\tpvalue\n";
    region_file << region_header_line.str();

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

        double sigma2 = Y.col(i).squaredNorm() / (n_samples - n_cov);
        Eigen::VectorXd w;
        if (params.model == "glm" || params.model == "glmm") {
            w = W_mat.row(i);
        } else {
            w = Eigen::VectorXd::Ones(n_samples);
        }
        // Pre-calculate matrices used in covariate adijustment.
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

            // Covariate-adjust and standardise g.
            g = G_slice.col(slice_ind); 
            g_s = g - X * (XtWX_inv * (Xt * g.cwiseProduct(w)));
            //standardise_vec(g_s);
            
            // Calculate the score statistic.
            u = g_s.cwiseProduct(w).dot(Y.col(i));
            gtg = g_s.cwiseProduct(w).dot(g_s);
            v = gtg;
            if (params.model == "lmm" || params.model == "glmm") {
                v *= sigma2;
            }
            beta = u / v;
            se = 1 / std::sqrt(v);
            if (v > 0) {
                zscore = u / std::sqrt(v);
                pval_esnp = 2 * pnorm(std::abs(zscore), true);
            } else {
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
        std::stringstream region_line;

        region_line << 
            pheno_data.pheno_ids[i] << "\t" <<
            pheno_data.chrom[i] << "\t" <<
            pheno_data.start[i] << "\t" <<
            ACAT(pvals) << "\n";
        region_file << region_line.str();
    }

    region_file.close();
    variant_file.close();
}
