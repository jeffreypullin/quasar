/*
  Copyright (C) 2024-25 Jeffrey Pullin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 or 3 of the License
  (at your option).

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is available at
   http://www.r-project.org/Licenses/
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
      void write_pheno_data(std::string out_file);

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