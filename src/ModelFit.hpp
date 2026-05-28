/*
  Copyright (C) 2024-26 Jeffrey Pullin

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
    std::vector<double> tr_int;
    std::vector<double> sigma2;
    std::vector<Eigen::MatrixXd> XtWX_inv_vec;
    std::vector<Eigen::VectorXd> Xty_res_vec;
    std::vector<Eigen::MatrixXd> XtWZ_vec;
    std::vector<Eigen::VectorXd> ZtDy_res_vec;
    std::vector<Eigen::VectorXd> Zty_res_vec;
    std::vector<Eigen::MatrixXd> XtWDZ_vec;
    std::vector<Eigen::VectorXd> d_out_vec;
    std::vector<Eigen::VectorXd> dw_out_vec;
    std::vector<Eigen::VectorXd> dwd_out_vec;
    
    std::vector<bool> phi_converged;
    std::vector<bool> glm_converged;
    std::vector<bool> glmm_converged;

    // Per-group quantities for single-cell --cell-groups score tests.
    size_t n_groups = 0;
    std::vector<std::string> group_ids;
    std::vector<std::vector<Eigen::MatrixXd>> XtWX_inv_g_vec;
    std::vector<std::vector<Eigen::VectorXd>> Xty_res_g_vec;
    std::vector<std::vector<Eigen::MatrixXd>> XtWZ_g_vec;
    std::vector<std::vector<Eigen::VectorXd>> y_out_g_vec;
    std::vector<std::vector<Eigen::VectorXd>> mu_out_g_vec;
    std::vector<std::vector<double>> tr_g_vec;
    std::vector<std::vector<double>> sigma2_g_vec;
    std::vector<std::vector<bool>> glmm_converged_g_vec;

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