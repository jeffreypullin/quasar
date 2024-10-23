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


#ifndef GENO_H
#define GENO_H

#include <string>
#include <vector>
#include <cstdint>

#include "Eigen/Dense"

class GenoData {
  public:
    std::string bed_prefix;
    
    Eigen::MatrixXd genotype_matrix;
    
    int n_samples;
    std::vector<std::string> sample_ids;
    
    // SNP information.
    int n_snps;
    std::vector<int> chrom;
    std::vector<int> pos;
    std::vector<std::string> id;
    std::vector<std::string> allele1;
    std::vector<std::string> allele2;
    std::vector<size_t> index;
    std::vector<double> var;

    GenoData(std::string bed_prefix) {
      this->bed_prefix = bed_prefix;
    }

    void read_bim_file();
    void read_fam_file();
    void prepare_bed_file();
    void read_bed_file();
    void slice_samples(std::vector<std::string>& sample_ids);
    void compute_variant_var();
};

#endif
