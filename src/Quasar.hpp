/*
 * MIT License
 *
 * Copyright (c) 2024 Jeffrey Pullin
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

#ifndef QUASAR_H
#define QUASAR_H

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <string>

#include "Eigen/Dense"
#include "Eigen/StdVector"

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;
typedef Eigen::Array<bool,Eigen::Dynamic,1> ArrayXb;
typedef Eigen::Matrix<bool,Eigen::Dynamic,1> VectorXb;
typedef Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> MatrixXb;
typedef Eigen::Map<Eigen::ArrayXd > MapArXd;
typedef Eigen::Map<const Eigen::ArrayXd > MapcArXd;
typedef Eigen::Map<Eigen::ArrayXf > MapArXf;
typedef Eigen::Map<const Eigen::ArrayXf > MapcArXf;
typedef Eigen::Map<Eigen::MatrixXd > MapMatXd;
typedef Eigen::Map<const Eigen::MatrixXd > MapcMatXd;
typedef Eigen::Map<Eigen::MatrixXf > MapMatXf;
typedef Eigen::Map<const Eigen::MatrixXf > MapcMatXf;
typedef Eigen::Map<ArrayXb> MapArXb;
typedef Eigen::Map<const ArrayXb> MapcArXb;
typedef Eigen::Array<uint16_t,Eigen::Dynamic,1> ArrayXt;
typedef Eigen::Array<uint64,Eigen::Dynamic,1> ArrayXui;

struct Param {

    // Input files.
    std::string bed_prefix;
    std::string grm_file;
    std::string feat_anno_file;
    std::string cov_file;
    std::string pheno_file;

    // Parameters.
    int window_size = 500000;

    // Data.
    std::string output_prefix;
};

#endif