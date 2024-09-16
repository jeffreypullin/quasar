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

#include "cxxopts.hpp"
#include <iostream>

int main(int argc, char* argv[]) {

    cxxopts::Options options("quasar", "QTL mapping software");

    // Define command-line arguments
    options.add_options()
        ("h,help", "Display help message")
        ("v,version", "Display version information")
        ("f,feat", "Feature annotation file", cxxopts::value<std::string>())
        ("c,cov", "Covariate file", cxxopts::value<std::string>())
        ("p,pheno", "Phenotype file", cxxopts::value<std::string>())
        ("b,bcf", "BCF file", cxxopts::value<std::string>())
        ("g,grm", "Genomic relatedness matrix", cxxopts::value<std::string>())
        ("o,output-prefix", "Output directory prefix", cxxopts::value<std::string>())
        ("w,window", "Cis window size in base pairs", cxxopts::value<int>());

    // Parse the arguments
    auto result = options.parse(argc, argv);

    // Handle arguments
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (result.count("version")) {
        std::cout << "quasar version 0.0.1" << std::endl;
        exit(0);
    }

    std::cout << "Execution finished." << std::endl;
    return 0;
}
