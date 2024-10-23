/*
 * MIT License
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "QTLmapping.hpp"
#include "LMMFitter.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <random>

void run_qtl_mapping_lmm(GenoData& geno_data, FeatData& feat_data, CovData& cov_data, PhenoData& pheno_data, GRM& grm, Regions& regions){

    Eigen::MatrixXd& Y = pheno_data.data;
	Eigen::MatrixXd& X = cov_data.data;

    Y = Y.transpose().eval();
    
	int n_samples = X.rows();
    int n_pheno = Y.cols();
	int n_cov = X.cols();

	pheno_data.std_dev.resize(n_pheno);

    Eigen::MatrixXd QY, QtY, QXtX, QtX, QUtG, QtG;

    std::cout << "\nPerforming rank normalization..." << std::endl;
    rank_normalize(Y); 
    std::cout << "Rank normalization finished." << std::endl;

    std::cout << "\nPerforming eigen decomposition of GRM..." << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(grm.mat);
    Eigen::MatrixXd Q = eig.eigenvectors();
    Eigen::VectorXd lambda = eig.eigenvalues();
    std::cout << "Eigen decomposition of GRM finished." << std::endl;

    std::cout << "\nComputing rotated matrices..." << std::endl;
    Eigen::MatrixXd G = geno_data.genotype_matrix.transpose().eval();

    QtY = (Q.transpose() * Y).eval();
    QXtX = (Q.transpose() * X).eval();
    QtX = (Q.transpose() * X).eval();
    QtG = (Q.transpose() * G).eval();
    QtG = (Q.transpose() * G).eval();

    std::cout << "Rotated matrices computed." << std::endl;
    
	std::vector<double> delta_vec(Y.cols());
	std::vector<double> hsq_vec(Y.cols());
	std::vector<double> sigma2_g_vec(Y.cols());
	std::vector<double> ssr_vec(Y.cols());

    std::cout << "\nEstimating variance components..." << std::endl;
    for (int i = 0; i < n_pheno; ++i) {

        DiagonalXd D_inv;
        double sigma2_g, sigma2_e, delta;

        LMMFitter fit(QtX, QtY.col(i), lambda);
        fit.fit_reml();
        D_inv = fit.D_inv;
		sigma2_g = fit.sigma2;
		sigma2_e = delta * sigma2_g;
		delta = fit.delta;

        double hsq = sigma2_g / (sigma2_e + sigma2_g);
		double scale = std::sqrt(sigma2_g + sigma2_e);

        // Residualise data?
        Eigen::MatrixXd XtDX = X.transpose() * D_inv * X;
		Eigen::VectorXd XtDy = X.transpose() * D_inv * Y.col(i);
		Eigen::VectorXd alpha = XtDX.colPivHouseholderQr().solve(XtDy);
		Eigen::VectorXd y_hat = X * alpha;
		Eigen::VectorXd y_res = Y.col(i) - y_hat;
			
		ssr_vec[i] = y_res.dot(D_inv * Y.col(i)) / sigma2_g;
		 
		Y.col(i) = (D_inv * y_res/std::sqrt(sigma2_g)).eval();

        delta_vec[i] = delta;
		hsq_vec[i] = hsq;
		sigma2_g_vec[i] = sigma2_g;
        pheno_data.std_dev[i] = std::sqrt(sigma2_g);
    }
    std::cout << "Variance components estimation finished." << std::endl;

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

double estimate_r(const Eigen::MatrixXd& X, const Eigen::MatrixXd& GRM, const Eigen::MatrixXd& G, double sigma2_g, double delta){

	int n_rand_samples = 30;
	
	int n_var = G.cols();
	int n_samples = GRM.rows();

	std::random_device rd;
	std::mt19937 gen(rd());

	double sigma2_e = delta * sigma2_g;
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_samples, n_samples);
	Eigen::MatrixXd Omega = sigma2_g * GRM + sigma2_e * I;
	Eigen::MatrixXd Omega_inv = Omega.inverse();

	double var_num, var_denom;
	Eigen::VectorXd r_vec(n_rand_samples);

	Eigen::MatrixXd Omega_inv_X = Omega_inv * X;
	Eigen::MatrixXd Xt_Omega_inv_X = X.transpose() * Omega_inv_X;
	Eigen::MatrixXd Xt_Omega_inv_X_inv = Xt_Omega_inv_X.inverse();

	for (int i = 0; i < n_rand_samples; ++i) {	

		std::uniform_int_distribution<> dis(0, n_var - 1);
		int rand_ind = dis(gen);

		Eigen::VectorXd g_rand = G.col(rand_ind);

		// FIXME: Why are there negative values?
		if (g_rand.cwiseAbs().sum() <= 20) {
			continue;
		}

		Eigen::VectorXd Omega_inv_g = Omega_inv * g_rand;
    	double a = g_rand.dot(Omega_inv_g);
    
		Eigen::VectorXd Xt_Omega_inv_g = X.transpose() * Omega_inv_g;
    	Eigen::VectorXd temp = Xt_Omega_inv_X_inv * Xt_Omega_inv_g;
    	double b = g_rand.dot(Omega_inv_X * temp);

		var_num = a + b;
		var_denom = g_rand.squaredNorm();

		r_vec(i) = var_num / var_denom;
	}
	return r_vec.mean();
}

double pnorm(double x, bool lower) {
	boost::math::normal N01(0.0, 1.0);
	if (lower) { 
        return boost::math::cdf(boost::math::complement(N01, x));
    } 
	return boost::math::cdf(N01, x);
}

double p_bd = 1e-300;
double q_bd = 3e+299;

double qnorm(double p, bool lower){
	boost::math::normal N01(0.0, 1.0);
	if (lower) { 
        return boost::math::quantile(boost::math::complement(N01, p));
    } 
	return boost::math::quantile(N01, p);
}

double qcauchy(double p, bool lower){
	p = p > p_bd ? p : p_bd;
	p = p < 1 - p_bd ? p : 1 - p_bd;
	
	boost::math::cauchy C01(0.0, 1.0);
	if (lower) {
		return boost::math::quantile(boost::math::complement(C01, p));
	}
	return boost::math::quantile(C01, p);
}

double pcauchy(double x, bool lower){
	x = x < q_bd ? x : q_bd;
	x = x > -q_bd ? x : -q_bd;
	
	boost::math::cauchy C01(0.0, 1.0);
	if (lower) {
		return boost::math::cdf(boost::math::complement(C01, x));
	}
	return boost::math::cdf(C01, x);
}

void rank_normalize(Eigen::MatrixXd& Y){
	double n = Y.rows();
	double p = Y.cols();
	
	std::vector<double> z((int) n);
	std::vector<double> rk((int) n);
	
	double mu = 0;
	double sd = 0;
	for (int i = 0; i < n; ++i){
		z[i] = qnorm(((double)i+1.0) / ((double)n+1.0), true);
		mu += z[i];
		sd += z[i] * z[i];
	}
	sd = std::sqrt(sd / (n - 1) - mu * mu / (n * (n - 1.0)));
	mu = mu / n;
	for (int i = 0; i < n; ++i){
		z[i] = (z[i] - mu) / sd;
	}
	
	for (int j = 0; j < p; ++j){

		std::vector<double> v(n);
		for (int i = 0; i < n; ++i){
			v[i] = Y(i, j);
		}
		
		std::vector<int> ranks = rank_vector(v);
		for (int i = 0; i < n; ++i){
			Y(i, j) = z[ranks[i] - 1];
		}
	}
}

std::vector<int> rank_vector(const std::vector<double>& v){
    
	std::vector<size_t> w(v.size());
    std::iota(w.begin(), w.end(), 0);
    std::sort(w.begin(), w.end(), [&v](size_t i, size_t j) { return v[i] < v[j]; });

    std::vector<int> r(w.size());
    for (size_t n, i = 0; i < w.size(); i += n) {
        n = 1;
        while (i + n < w.size() && v[w[i]] == v[w[i+n]]) ++n;
        for (size_t k = 0; k < n; ++k) {
            r[w[i+k]] = i + (n + 1) / 2.0;
        }
    }
    return r;
}

double ACAT(const std::vector<double>& pvals) {
	long double sum = 0.0;
	double n = pvals.size();
	for (const double& p: pvals) {
		if (p >= 1){
			sum += qcauchy(1 - 1 / n, true);
		} else if (p <= 0){
			std::cerr << "ACAT failed; input pval <= 0. \n";
			exit(1);
		} else {
			sum += qcauchy(p, true);
		}
	}
	return pcauchy(sum, true);
}
