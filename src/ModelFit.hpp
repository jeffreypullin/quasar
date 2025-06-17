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