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

#ifndef MODELFIT_H
#define MODELFIT_H

#include "Data.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <map>

#include "Eigen/Dense"

class ModelFit {
  public:

    std::string model;
    std::string file;

    size_t n_pheno;
    size_t n_samples;
    std::vector<std::string> pheno_ids;
    std::vector<std::string> sample_ids;
    
    Eigen::MatrixXd W;
    
    std::vector<double> phi;
    std::vector<double> tr;
    
    std::vector<bool> phi_converged;
    std::vector<bool> glm_converged;
    std::vector<bool> glmm_converged;

    ModelFit(std::string model, std::string fit_file, PhenoData& pheno_data) {
      this->model = model;
      this->file = fit_file;
      this->pheno_ids = pheno_data.pheno_ids;
      this->sample_ids = pheno_data.sample_ids;
      this->n_pheno = pheno_data.n_pheno;
      this->n_samples = pheno_data.n_samples;
    }

    void write_model_fit(std::string out);
    void read_model_fit();
};

#endif