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

#include "Geno.hpp"
#include "Data.hpp"

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/StdVector"
#include "cxxopts.hpp"

struct Params {

    // Input files.
    std::string bed_prefix;
    std::string grm_file;
    std::string feat_file;
    std::string cov_file;
    std::string pheno_file;

    // Statistical model, can be either
    // "lmm": linear mixed model or
    // "glmm": generalised linear mixed model. 
    std::string model;
    // QTL mapping parameters.
    // By default we use a +/- 1Mb window.
    int window_size = 1000000;
    // The number of genotypes to sample
    // when estimating the r for
    // main effect QTL mapping.
    int main_n_rand_samples = 30;

    // Output parameters.
    std::string output_prefix;
    bool verbose = false;
};

#endif