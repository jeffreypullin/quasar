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

#include "Residualise.hpp"
#include "QTLMappingUtils.hpp"
#include "Data.hpp"
#include "ModelFit.hpp"
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

void residualise(Params& params, ModelFit& model_fit, CovData& cov_data, PhenoData& pheno_data, GRM& grm) {

    Eigen::MatrixXd& Y = pheno_data.data;
    Eigen::MatrixXd& X = cov_data.data;

    int n_pheno = pheno_data.n_pheno;

    if (params.model == "lmm" || params.model == "lm") {
        std::cout << "\nPerforming rank normalization..." << std::endl;
        rank_normalize(Y); 
        std::cout << "Rank normalization finished." << std::endl;
    }

    // Compute offset for count-distribution models.
    Eigen::VectorXd offset = Y.rowwise().sum().array().log();
    
    Eigen::MatrixXd W(n_pheno, pheno_data.n_samples);

    std::vector<double> tr;
    std::vector<double> phi;

    std::vector<bool> glm_converged;
    std::vector<bool> phi_converged;
    std::vector<bool> glmm_converged;

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
        for (int i = 0; i < n_pheno; ++i) {

            LMM lmm(QtX, QtY.col(i), lambda);
            lmm.fit();

            Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_inv = lmm.D_inv;
            Y.col(i) = (Q * D_inv * (QtY.col(i) - QtX * lmm.beta)) / std::sqrt(lmm.sigma2);
        
        }
        std::cout << "Null LMMs fitted." << std::endl;

    } else if (params.model == "lm") {

        std::cout << "\nFitting null LMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            LM lm(X, Y.col(i));
            lm.fit();

            Y.col(i) = (Y.col(i) - X * lm.beta) / std::sqrt(lm.s);
        }
        std::cout << "Null LMs fitted." << std::endl;

    } else if (params.model == "p_glmm") {
        
        std::cout << "\nFitting null Poisson GLMMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLMM p_glmm(X, Y.col(i), offset, std::move(poisson), grm.mat);
            p_glmm.fit();

            Y.col(i) = (Y.col(i).array() - (X * p_glmm.beta).array().exp()) / p_glmm.mu.array();

            W.row(i) = p_glmm.mu.array();
            tr.push_back(p_glmm.P.trace() / p_glmm.mu.sum());
            glmm_converged.push_back(p_glmm.glmm_converged);
        }
        std::cout << "Null Poisson GLMMs fitted." << std::endl;

    } else if (params.model == "nb_glm") {

        std::cout << "\nFitting null NB-GLMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            bool use_apl = params.use_apl;
            NBGLM nb_glm(X, Y.col(i), offset, use_apl);
            nb_glm.fit();

            Y.col(i) = (Y.col(i).array() - (X * nb_glm.beta).array().exp()) / nb_glm.mu.array();
            W.row(i) = nb_glm.mu.array() / (1 + nb_glm.phi * nb_glm.mu.array());

            phi.push_back(nb_glm.phi);
            phi_converged.push_back(nb_glm.phi_converged);
            glm_converged.push_back(nb_glm.glm_converged);
        }
        std::cout << "Null NB-GLMs fitted." << std::endl;

    } else if (params.model == "p_glm") {

        std::cout << "\nFitting null Poisson-GLMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            auto poisson = std::unique_ptr<Family>(new Poisson());
            GLM p_glm(X, Y.col(i), offset, std::move(poisson));
            p_glm.fit();

            Y.col(i) = (Y.col(i).array() - (X * p_glm.beta).array().exp()) / p_glm.mu.array();
            W.row(i) = p_glm.mu.array();
            glm_converged.push_back(p_glm.glm_converged);
        }
        std::cout << "Null Poisson GLMs fitted." << std::endl;

    } else if (params.model == "nb_glmm") {

        std::cout << "\nFitting null NB GLMMs..." << std::endl;
        for (int i = 0; i < n_pheno; ++i) {

            NBGLMM nb_glmm(X, Y.col(i), offset, grm.mat);
            nb_glmm.fit();

            Y.col(i) = (Y.col(i).array() - (X * nb_glmm.beta).array().exp()) / nb_glmm.mu.array();
            W.row(i) = nb_glmm.mu.array() / (1 + nb_glmm.phi * nb_glmm.mu.array());
            
            tr.push_back(nb_glmm.P.trace() / nb_glmm.mu.sum());
            phi.push_back(nb_glmm.phi);
            phi_converged.push_back(nb_glmm.phi_converged);
            glmm_converged.push_back(nb_glmm.glmm_converged);
        }
        std::cout << "Null NB GLMMs fitted." << std::endl;
    
    }

    model_fit.W = W;
    model_fit.phi = phi;
    model_fit.tr = tr;

    model_fit.phi_converged = phi_converged;
    model_fit.glm_converged = glm_converged;
    model_fit.glmm_converged = glmm_converged;

    pheno_data.data = Y;
}
