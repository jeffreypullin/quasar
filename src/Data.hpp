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

#ifndef DATA_H
#define DATA_H

#include "Quasar.hpp"

class PhenoData {

    public:
      std::string pheno_file;
      size_t n_pheno;
      size_t n_samples;
      std::vector<std::string> pheno_ids;
      std::vector<std::string> sample_ids;
      Eigen::MatrixXd data;

      std::vector<int> chrom;
      std::vector<int> start;
      std::vector<int> end;
      std::vector<int> window_start;
      std::vector<int> window_end;
      std::vector<int> window_n;

      PhenoData(std::string pheno_file) {
        this->pheno_file = pheno_file;
      }
      void read_pheno_data();

      void construct_windows(GenoData& geno_data, int window_size, bool verbose); 

      void slice_chromosome(int chrom_id);
      void slice_samples(std::vector<std::string>& sample_ids);
};

class CovData {

    public:
      std::string cov_file;
      size_t n_cov;
      size_t n_samples;
      std::vector<std::string> cov_ids;
      std::vector<std::string> sample_ids;
      Eigen::MatrixXd data;
      CovData(std::string cov_file) {
        this->cov_file = cov_file;
      }
      void read_cov_data();
      void slice_samples(std::vector<std::string>& sample_ids);
};

class GRM {

  public:
    std::string grm_file;
    size_t n_samples;
    std::vector<std::string> sample_ids;
    Eigen::MatrixXd mat;
    GRM(std::string grm_file) {
      this->grm_file = grm_file;
    }
    void slice_samples(std::vector<std::string>& sample_ids);
    void read_grm();
};

#endif