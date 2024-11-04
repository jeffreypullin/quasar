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
#include "LMM.hpp"
#include "GLMM.hpp"
#include "QTLMappingUtils.hpp"
#include "VarRatioApprox.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>

void run_qtl_mapping_lmm(GenoData& geno_data, FeatData& feat_data, CovData& cov_data, PhenoData& pheno_data, GRM& grm, Regions& regions){

    Eigen::MatrixXd& Y = pheno_data.data;
	Eigen::MatrixXd& X = cov_data.data;

    Y = Y.transpose().eval();
    Eigen::MatrixXd G = geno_data.genotype_matrix.transpose().eval();
    
	int n_samples = X.rows();
    int n_pheno = Y.cols();
	int n_cov = X.cols();

	pheno_data.std_dev.resize(n_pheno);

    Eigen::MatrixXd QtY, QtX;

    std::cout << "\nPerforming rank normalization..." << std::endl;
    rank_normalize(Y); 
    std::cout << "Rank normalization finished." << std::endl;

    std::cout << "\nPerforming eigen decomposition of GRM..." << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(grm.mat);
    Eigen::MatrixXd Q = eig.eigenvectors();
    Eigen::VectorXd lambda = eig.eigenvalues();
    std::cout << "Eigen decomposition of GRM finished." << std::endl;

    std::cout << "\nComputing rotated matrices..." << std::endl;
    QtY = (Q.transpose() * Y).eval();
    QtX = (Q.transpose() * X).eval();
    std::cout << "Rotated matrices computed." << std::endl;
    
	std::vector<double> delta_vec(Y.cols());
	std::vector<double> sigma2_g_vec(Y.cols());

    std::cout << "\nFitting null LMMs..." << std::endl;
    for (int i = 0; i < n_pheno; ++i) {

        LMM lmm(QtX, QtY.col(i), lambda);
        lmm.fit();

		// Residualise phenotype.
		Eigen::VectorXd y_hat = X * lmm.beta;
		Y.col(i) = Y.col(i) - y_hat;

        delta_vec[i] = lmm.delta;
		sigma2_g_vec[i] = lmm.sigma2;
        pheno_data.std_dev[i] = std::sqrt(lmm.sigma2);
    }
    std::cout << "Null LMMs fitted." << std::endl;

	std::cout << "\nEstimating variance ratio correction factors..." << std::endl;
	std::vector<double> r_vec(n_pheno);
	for (int i = 0; i < n_pheno; ++i){
		r_vec[i] = estimate_r(X, grm.mat, G, sigma2_g_vec[i], delta_vec[i]);
	}
	std::cout << "Correction factor estimation finished." << std::endl;

    std::cout << "\nCalculating variant significance..." << std::endl;

    // Iterate over regions.
    for (int i = 0; i < regions.size(); ++i){

        std::cout << "Processing region " << (i + 1) << " (" << regions.regions_id[i] << ")..." << std::endl;

		int n_g = regions.geno_e[i] - regions.geno_s[i] + 1;
        int n_e = regions.feat_e[i] - regions.feat_s[i] + 1;

        int s_slice = 0;
		int n_slice = G.cols();

        // Iterate over phenotypes in a region.
        for (int jj = regions.feat_s[i], jm = 0; jj < regions.feat_s[i] + n_e; ++jj, ++jm) {

            std::cout << "Processing feature: " << feat_data.feat_id[jj] << std::endl;
            
			std::stringstream block_line;
				
			std::vector<double> pvals;
			std::vector<double> dist;
			double gene_tss = feat_data.start[jj];
				
			int idx_s = -1;
			int idx_n = 0;

			for (int ii = regions.geno_s[i], im = 0; ii < regions.geno_s[i] + n_g; ++ii, im++ ) {
				if (geno_data.chrom[ii] == feat_data.chrom[jj] && feat_data.start[jj] - geno_data.pos[ii] < 500000 && geno_data.pos[ii] - feat_data.end[jj] < 500000) {
					if (idx_s < 0) {
                        idx_s = im;
                    }
					idx_n++;
				} else {
					if (idx_s >= 0) {
                        break;
                    }
				}
			}
				
			if (idx_s < 0) {
                continue;
			} else if (idx_s >= G.cols()) {
				idx_s = G.cols() - 1;
			}
			if (idx_s + idx_n > G.cols()) {
				idx_n = G.cols() - idx_s;
			}
			if (n_slice <= 0 ){
				continue;
			}
				
			int v_s = idx_s + s_slice;
			int v_e = v_s + idx_n - 1;
			
			int pos_s = geno_data.pos[v_s];
			int pos_e = geno_data.pos[v_e];

            const Eigen::MatrixXd& G_slice = G.middleCols(idx_s, idx_n);

            std::cout << "Processing " << idx_n << " SNPs from " << pos_s << " to " << pos_e << std::endl;
            
            // Iterate over SNPs, calculating significance.
            for (int ii = v_s, im = idx_s, si = 0; im < idx_s + idx_n; ii++, im++, si++){
				
					std::stringstream long_line;

					double u, v, gtg, zscore, pval_esnp;

					if (geno_data.var[ii] > 0.0) {

						Eigen::VectorXd g = G_slice.col(si);
						// Project out covariates.
						Eigen::VectorXd g_tilde = g - X * (X.transpose() * X).ldlt().solve(X.transpose() * g);

						u = g_tilde.dot(Y.col(jj));
						gtg = g_tilde.squaredNorm();
						std::cout << "gtg: " << gtg << std::endl;
						std::cout << "r: " << r_vec[jj] << std::endl;
						v = r_vec[jj] * gtg;
						if (v > 0) {
							zscore = u / std::sqrt(v);
							pval_esnp = 2 * pnorm(std::abs(zscore), true);
						} else{
							zscore = 0;
							pval_esnp = -1;
						}
					} else {
						zscore = 0;
						pval_esnp = -1;
					}
					pvals.push_back(pval_esnp);

                    std::cout << "u: " << u << std::endl;
                    std::cout << "v: " << v << std::endl;
                    std::cout << "z-score: " << zscore << std::endl;
                    std::cout << "p-value: " << pval_esnp << std::endl;
            }
			
			// TODO Turn back on.
			if (false && pvals.size() > 0) {
						
				std::stringstream bed_block_line;

				std::cout << "Region results..." << std::endl;
				std::cout << feat_data.chrom[jj] << std::endl;
				std::cout << feat_data.start[jj] << std::endl;
				std::cout << feat_data.end[jj] << std::endl;
				std::cout << feat_data.feat_id[jj] << std::endl;
				std::cout << ACAT(pvals) << std::endl;
				std::cout << n_samples << std::endl;
				std::cout << n_cov << std::endl;
				std::cout << pheno_data.std_dev[jj] << std::endl;
				std::cout << idx_n << std::endl;
			}
		}
    }
}

void run_qtl_mapping_glmm(GenoData& geno_data, FeatData& feat_data, CovData& cov_data, PhenoData& pheno_data, GRM& grm, Regions& regions){

    Eigen::MatrixXd& Y = pheno_data.data;
	Eigen::MatrixXd& X = cov_data.data;

    Y = Y.transpose().eval();
    Eigen::MatrixXd G = geno_data.genotype_matrix.transpose().eval();
    
	int n_samples = X.rows();
    int n_pheno = Y.cols();
	int n_cov = X.cols();

	pheno_data.std_dev.resize(n_pheno);
	std::vector<double> sigma2_vec(Y.cols());

    std::cout << "\nFitting null GLMMs..." << std::endl;
	// Change back.
	std::vector<double> r_vec(n_pheno);
    for (int i = 0; i < n_pheno; ++i) {

        double tau1, tau2;

		Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_samples, n_samples);
		std::vector<Eigen::MatrixXd> Ks;
		Ks.push_back(grm.mat);
		Ks.push_back(I);

		Eigen::VectorXd lib_size = Eigen::VectorXd::Ones(n_samples);

		// Fit a GLM to get initial beta.
		Eigen::VectorXd y = Y.col(i).array() + 300;
		GLM glm(X, y);
		glm.fit();
		
		Eigen::VectorXd init_beta = glm.beta;	
		Eigen::VectorXd init_tau = Eigen::VectorXd::Ones(2);
		
		std::cout << "i: " << i << std::endl;
        GLMM glmm(X, y, Ks, lib_size, init_beta, init_tau);
        glmm.fit();
		tau1 = glmm.tau[0];
		std::cout << "tau1: " << tau1 << std::endl;
		tau2 = glmm.tau[1];
		std::cout << "tau2: " << tau2 << std::endl;
		
		Eigen::VectorXd y_hat = X * glmm.beta;
		Eigen::VectorXd y_res = Y.col(i) - y_hat;

		// FIXME: Check this.
		double sigma2 = glmm.tau.sum();
		if (sigma2 > 0) {
			Y.col(i) = (y_res / sigma2).eval();
		} else {
			Y.col(i) = y_res;
		}

		sigma2_vec[i] = sigma2;
        pheno_data.std_dev[i] = std::sqrt(sigma2);

		r_vec[i] = estimate_r_glmm(X, grm.mat, G, glmm.P, glmm.d);
    }
    std::cout << "Null GLMMs fitted." << std::endl;
	std::cout << "Correction factor estimation finished." << std::endl;

    std::cout << "\nCalculating variant significance..." << std::endl;

    // Iterate over regions.
    for (int i = 0; i < regions.size(); ++i){

        std::cout << "Processing region " << (i + 1) << " (" << regions.regions_id[i] << ")..." << std::endl;

		int n_g = regions.geno_e[i] - regions.geno_s[i] + 1;
        int n_e = regions.feat_e[i] - regions.feat_s[i] + 1;

        int s_slice = 0;
		int n_slice = G.cols();

        // Iterate over phenotypes in a region.
        for (int jj = regions.feat_s[i], jm = 0; jj < regions.feat_s[i] + n_e; ++jj, ++jm) {

            std::cout << "Processing feature: " << feat_data.feat_id[jj] << std::endl;
            
			std::stringstream block_line;
				
			std::vector<double> pvals;
			std::vector<double> dist;
			double gene_tss = feat_data.start[jj];
				
			int idx_s = -1;
			int idx_n = 0;

			for (int ii = regions.geno_s[i], im = 0; ii < regions.geno_s[i] + n_g; ++ii, im++ ) {
				if (geno_data.chrom[ii] == feat_data.chrom[jj] && feat_data.start[jj] - geno_data.pos[ii] < 500000 && geno_data.pos[ii] - feat_data.end[jj] < 500000) {
					if (idx_s < 0) {
                        idx_s = im;
                    }
					idx_n++;
				} else {
					if (idx_s >= 0) {
                        break;
                    }
				}
			}
				
			if (idx_s < 0) {
                continue;
			} else if (idx_s >= G.cols()) {
				idx_s = G.cols() - 1;
			}
			if (idx_s + idx_n > G.cols()) {
				idx_n = G.cols() - idx_s;
			}
			if (n_slice <= 0 ){
				continue;
			}
				
			int v_s = idx_s + s_slice;
			int v_e = v_s + idx_n - 1;
			
			int pos_s = geno_data.pos[v_s];
			int pos_e = geno_data.pos[v_e];

            const Eigen::MatrixXd& G_slice = G.middleCols(idx_s, idx_n);

            std::cout << "Processing " << idx_n << " SNPs from " << pos_s << " to " << pos_e << std::endl;

            Eigen::VectorXd U_vec = G_slice.transpose() * Y.col(jj);

            // Iterate over SNPs, computing statistics.
            for (int ii = v_s, im = idx_s, si = 0; im < idx_s + idx_n; ii++, im++, si++){
				
					std::stringstream long_line;
				
					double u = U_vec(si), v;
					double zscore, pval_esnp;
					
					if (geno_data.var[ii] > 0.0) {
						double gtg = G_slice.col(si).squaredNorm();
						std::cout << "gtg: " << gtg << std::endl;
						std::cout << "r: " << r_vec[jj] << std::endl;
						double v = r_vec[jj] * gtg;
						std::cout << "u: " << u << std::endl;
						std::cout << "v: " << v << std::endl;
						if (v > 0) {
							zscore = u / std::sqrt(v);
							pval_esnp = 2 * pnorm(std::abs(zscore), true);
						} else{
							zscore = 0;
							pval_esnp = -1;
						}
					} else {
						zscore = 0;
						pval_esnp = -1;
					}
					pvals.push_back(pval_esnp);

                    std::cout << "u: " << u << std::endl;
                    std::cout << "v: " << v << std::endl;
                    std::cout << "z-score: " << zscore << std::endl;
                    std::cout << "p-value: " << pval_esnp << std::endl;
            }
			
			// TODO Turn back on.
			if (false && pvals.size() > 0) {
						
				std::stringstream bed_block_line;

				std::cout << "Region results..." << std::endl;
				std::cout << feat_data.chrom[jj] << std::endl;
				std::cout << feat_data.start[jj] << std::endl;
				std::cout << feat_data.end[jj] << std::endl;
				std::cout << feat_data.feat_id[jj] << std::endl;
				std::cout << ACAT(pvals) << std::endl;
				std::cout << n_samples << std::endl;
				std::cout << n_cov << std::endl;
				std::cout << pheno_data.std_dev[jj] << std::endl;
				std::cout << idx_n << std::endl;
			}
		}
    }
}

void run_qtl_mapping_lmm_int(
	GenoData& geno_data, 
	FeatData& feat_data, 
	CovData& cov_data, 
	PhenoData& pheno_data, 
	GRM& grm, 
	Regions& regions, 
	std::vector<std::string> int_cov_ids
) {

 	Eigen::MatrixXd& Y = pheno_data.data;
    Y = Y.transpose().eval();
	int n_samples = Y.rows();

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

    Eigen::MatrixXd G = geno_data.genotype_matrix.transpose().eval();

    int n_pheno = Y.cols();
	int n_cov = X.cols();
	int n_int_covs = int_cov_ids.size();

	pheno_data.std_dev.resize(n_pheno);

    Eigen::MatrixXd QtY, QtX;

    std::cout << "\nPerforming rank normalization..." << std::endl;
    rank_normalize(Y); 
    std::cout << "Rank normalization finished." << std::endl;

    std::cout << "\nPerforming eigen decomposition of GRM..." << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(grm.mat);
    Eigen::MatrixXd Q = eig.eigenvectors();
    Eigen::VectorXd lambda = eig.eigenvalues();
    std::cout << "Eigen decomposition of GRM finished." << std::endl;

    std::cout << "\nComputing rotated matrices..." << std::endl;
    QtY = (Q.transpose() * Y).eval();
    QtX = (Q.transpose() * X).eval();
    std::cout << "Rotated matrices computed." << std::endl;
    
	std::vector<double> delta_vec(Y.cols());
	std::vector<double> sigma2_g_vec(Y.cols());

    std::cout << "\nFitting null LMMs..." << std::endl;
    for (int i = 0; i < n_pheno; ++i) {

        LMM lmm(QtX, QtY.col(i), lambda);
        lmm.fit();

		// Residualise phenotype.
		Eigen::VectorXd y_hat = X * lmm.beta;
		Y.col(i) = Y.col(i) - y_hat;

        delta_vec[i] = lmm.delta;
		sigma2_g_vec[i] = lmm.sigma2;
        pheno_data.std_dev[i] = std::sqrt(lmm.sigma2);
    }
    std::cout << "Null LMMs fitted." << std::endl;

	std::cout << "\nEstimating variance ratio correction factors..." << std::endl;
	std::vector<Eigen::MatrixXd> R_vec(n_pheno);
	for (int i = 0; i < n_pheno; ++i){
		R_vec[i] = estimate_R_int(cov_data, grm.mat, G, sigma2_g_vec[i], delta_vec[i], int_cov_ids);
	}
	std::cout << "Correction factor estimation finished." << std::endl;

    std::cout << "\nCalculating variant significance..." << std::endl;

    // Iterate over regions.
    for (int i = 0; i < regions.size(); ++i){

        std::cout << "Processing region " << (i + 1) << " (" << regions.regions_id[i] << ")..." << std::endl;

		int n_g = regions.geno_e[i] - regions.geno_s[i] + 1;
        int n_e = regions.feat_e[i] - regions.feat_s[i] + 1;

        int s_slice = 0;
		int n_slice = G.cols();

        // Iterate over phenotypes in a region.
        for (int jj = regions.feat_s[i], jm = 0; jj < regions.feat_s[i] + n_e; ++jj, ++jm) {

            std::cout << "Processing feature: " << feat_data.feat_id[jj] << std::endl;
            
			std::stringstream block_line;
				
			std::vector<double> pvals;
			std::vector<double> dist;
			double gene_tss = feat_data.start[jj];
				
			int idx_s = -1;
			int idx_n = 0;

			for (int ii = regions.geno_s[i], im = 0; ii < regions.geno_s[i] + n_g; ++ii, im++ ) {
				if (geno_data.chrom[ii] == feat_data.chrom[jj] && feat_data.start[jj] - geno_data.pos[ii] < 500000 && geno_data.pos[ii] - feat_data.end[jj] < 500000) {
					if (idx_s < 0) {
                        idx_s = im;
                    }
					idx_n++;
				} else {
					if (idx_s >= 0) {
                        break;
                    }
				}
			}
				
			if (idx_s < 0) {
                continue;
			} else if (idx_s >= G.cols()) {
				idx_s = G.cols() - 1;
			}
			if (idx_s + idx_n > G.cols()) {
				idx_n = G.cols() - idx_s;
			}
			if (n_slice <= 0 ){
				continue;
			}
				
			int v_s = idx_s + s_slice;
			int v_e = v_s + idx_n - 1;
			
			int pos_s = geno_data.pos[v_s];
			int pos_e = geno_data.pos[v_e];

            const Eigen::MatrixXd& G_slice = G.middleCols(idx_s, idx_n);

            std::cout << "Processing " << idx_n << " SNPs from " << pos_s << " to " << pos_e << std::endl;
            
           // Iterate over SNPs, calculating significance.
            for (int ii = v_s, im = idx_s, si = 0; im < idx_s + idx_n; ii++, im++, si++){
				
					std::stringstream long_line;

					Eigen::VectorXd U;
					Eigen::MatrixXd V;
					double chi, pval_esnp;

					if (geno_data.var[ii] > 0.0) {

						Z.col(0) = G_slice.col(si);
						// Project out covariates.
						Eigen::MatrixXd Z_tilde = Z - X * (X.transpose() * X).ldlt().solve(X.transpose() * Z);

						U = Z_tilde.transpose() * Y.col(jj);
						Eigen::MatrixXd ZtZ = Z_tilde.transpose() * Z_tilde;
						V = R_vec[jj] * ZtZ.trace();
						chi = U.transpose() * V * U;
						pval_esnp = pchisq(chi, Z_tilde.cols(), true);
					} else {
						chi = 0;
						pval_esnp = -1;
					}
					pvals.push_back(pval_esnp);

                    std::cout << "chi-square: " << chi << std::endl;
                    std::cout << "p-value: " << pval_esnp << std::endl;
            } 
			
			// TODO Turn back on.
			if (false && pvals.size() > 0) {
						
				std::stringstream bed_block_line;

				std::cout << "Region results..." << std::endl;
				std::cout << feat_data.chrom[jj] << std::endl;
				std::cout << feat_data.start[jj] << std::endl;
				std::cout << feat_data.end[jj] << std::endl;
				std::cout << feat_data.feat_id[jj] << std::endl;
				std::cout << ACAT(pvals) << std::endl;
				std::cout << n_samples << std::endl;
				std::cout << n_cov << std::endl;
				std::cout << pheno_data.std_dev[jj] << std::endl;
				std::cout << idx_n << std::endl;
			}
		}
    }
}