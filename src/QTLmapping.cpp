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

void run_qtl_mapping(GenoData& geno_data, FeatData& feat_data, CovData& cov_data, PhenoData& pheno_data, GRM& grm, Regions& regions){

    Eigen::MatrixXd& Y = pheno_data.data;
	Eigen::MatrixXd& X = cov_data.data;

    Eigen::MatrixXd QY, QtY, QXtX, QtX, QUtG, QtG;

    std::cout << "\nScaling and centering phenotype data..." << std::endl;
    std::vector<double> y_scale;
    //scale_and_center(Y, y_scale);

    std::cout << "\nPerforming eigen decomposition of GRM..." << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(grm.mat);
    Eigen::MatrixXd Q = eig.eigenvectors();
    Eigen::VectorXd lambda = eig.eigenvalues();
    std::cout << "Eigen decomposition of GRM finished." << std::endl;

    std::cout << "\nComputing rotated matrices..." << std::endl;
    Y = Y.transpose().eval();
    Eigen::MatrixXd G = geno_data.genotype_matrix.transpose().eval();

    QtY = (Q.transpose() * Y).eval();
    QXtX = (Q.transpose() * X).eval();
    QtX = (Q.transpose() * X).eval();
    QtG = (Q.transpose() * G).eval();
    std::cout << "Rotated matrices computed." << std::endl;

    std::cout << "\nEstimating variance components..." << std::endl;
    for (int i = 0; i < Y.cols(); ++i) {
        LMMFitter lmm_fitter(QtX, QtY.col(i), lambda);
        lmm_fitter.fit_reml();
    }
    std::cout << "Variance components estimation finished." << std::endl;

    /*
    // Calculate P matrix naively.
    for (int i = 0; i < Y.cols(); ++i) {
        Eigen::MatrixXd P = QtX * (QtX.transpose() * QtX).ldlt().solve(QtX.transpose());
    }

    // Calculate variant significance. 
    for (int i = 0; i < regions.size(); ++i){

 		const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(s_slice, n_slice);

        for (int j = 0; j < Y.cols(); ++j) {
            Eigen::VectorXd U_vec = G_slice.transpose() * Y.col(jj);

        }
    }
    */

}

/*
void scale_and_center(Eigen::MatrixXd& Y, std::vector<double>& y_scale){

}
*/