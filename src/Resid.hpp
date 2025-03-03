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

void residualise(Params& params, Eigen::MatrixXd& Y, Eigen::MatrixXd& X, GRM& grm);
void residualise_lmm(Eigen::MatrixXd& Y, Eigen::MatrixXd& X, GRM& grm);
void residualise_glmm(Eigen::MatrixXd& Y, Eigen::MatrixXd& X, GRM& grm);
void residualise_lm(Eigen::MatrixXd& Y, Eigen::MatrixXd& X);
void residualise_glm(Eigen::MatrixXd& Y, Eigen::MatrixXd& X);
