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

#include "Resid.hpp"
#include "GLM.hpp"
#include "LM.hpp"
#include "LMM.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>

void residualise(Params& params, Eigen::MatrixXd& Y, Eigen::MatrixXd& X, GRM& grm) {
    std::string model = params.model;
    if (model == "lmm") {
        residualise_lmm(Y, X, grm);
    } else if (model == "glmm") {
        residualise_glmm(Y, X, grm);
    } else if (model == "lm") {
        residualise_lm(Y, X);
    } else if (model == "glm") {
        residualise_glm(Y, X);
    }
}

void residualise_lmm(Eigen::MatrixXd& Y, Eigen::MatrixXd& X, GRM& grm) {

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

        // Fit the LMM.
        LMM lmm(QtX, QtY.col(i), lambda);
        lmm.fit();

        // Calculate Py.
        double sigma2 = lmm.sigma2;
        Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_inv = lmm.D_inv;
        Y.col(i) = (Q.transpose() * D_inv * (QtY.col(i) - QtX * lmm.beta)) / std::sqrt(sigma2);
        
    }
    std::cout << "Null LMMs fitted." << std::endl;

}

void residualise_glmm(Eigen::MatrixXd& Y, Eigen::MatrixXd& X, GRM& grm) {
    std::cerr << "\nError: Not yet implemented." << std::endl;
    exit(1);
}

void residualise_lm(Eigen::MatrixXd& Y, Eigen::MatrixXd& X) {
    std::cout << "\nFitting null LMs..." << std::endl;
    for (int i = 0; i < Y.cols(); ++i) {

        // Fit the LM.
        LM lm(X, Y.col(i));
        lm.fit();

        // Calculate Py.
        Y.col(i) = Y.col(i) - X * lm.beta;

    }
    std::cout << "Null LMs fitted." << std::endl;
}

void residualise_glm(Eigen::MatrixXd& Y, Eigen::MatrixXd& X) {
    std::cerr << "\nError: Not yet implemented." << std::endl;
    exit(1);
}




