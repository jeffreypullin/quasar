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

#ifndef GENO_H
#define GENO_H

#include <string>
#include <vector>
#include <cstdint>
#include <map>

#include "Eigen/Dense"

class GenoData {
  public:
    std::string bed_prefix;
    
    Eigen::MatrixXd genotype_matrix;
    
    size_t n_samples;
    std::vector<std::string> sample_ids;
    
    // SNP information.
    size_t n_snps;
    std::vector<int> chrom;
    std::vector<int> pos;
    std::vector<std::string> id;
    std::vector<std::string> alt;
    std::vector<std::string> ref;
    std::vector<std::string> snp_id;
    std::vector<size_t> index;

    GenoData(std::string bed_prefix) {
      this->bed_prefix = bed_prefix;
    }

    void read_bim_file();
    void read_fam_file();
    void prepare_bed_file();
    void run_mean_imputation();
    void read_bed_file();
    void slice_samples(std::vector<std::string>& sample_ids);
};

#endif
