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

#include "ModelFit.hpp"
#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void ModelFit::write_model_fit(std::string out) {

    std::ofstream lst_file(out + "-fit.lst");
    lst_file << model << "\n";
    if (model != "lm" && model != "lmm") {
      lst_file << out + "-fit-ww.tsv\n";
      lst_file << out + "-fit-info.tsv\n";
    }
    lst_file.close();

    if (model != "lm" && model != "lmm") {
        // Write working weight file.
        std::ofstream ww_file(out + "-fit-ww.tsv");
        ww_file << "feature_id";
        for (size_t i = 0; i < sample_ids.size(); ++i) {
            ww_file << "\t" << sample_ids[i];
        }
        ww_file << "\n";
        for (size_t i = 0; i < n_pheno; ++i) {
            ww_file << pheno_ids[i];
            for (size_t j = 0; j < n_samples; ++j) {
                ww_file << "\t" << W(i,j);
            }
            ww_file << "\n";
        }
        ww_file.close();

        // Write information file.
        std::ofstream info_file(out + "-fit-info.tsv");
        info_file << "feature_id";
        if (model == "p_glm") {
            info_file << "\tglm_converged";
        } else if (model == "nb_glm") {
            info_file << "\tglm_converged\tphi\tphi_converged";
        } else if (model == "p_glmm") {
            info_file << "\tglmm_converged\ttr";
        } else if (model == "nb_glmm") {
            info_file << "\tglmm_converged\ttr\tphi\tphi_converged";
        }
        info_file << "\n";

        for (size_t i = 0; i < n_pheno; ++i) {
            info_file << pheno_ids[i];
            if (model == "p_glm") {
                info_file << "\t" << glm_converged[i];
            } else if (model == "nb_glm") {
                info_file << "\t" << glm_converged[i] << "\t" << phi[i] << "\t" << phi_converged[i];
            } else if (model == "p_glmm") {
                info_file << "\t" << glmm_converged[i] << "\t" << tr[i];
            } else if (model == "nb_glmm") {
                info_file << "\t" << glmm_converged[i] << "\t" << tr[i] << "\t" << phi[i] << "\t" << phi_converged[i];
            }
            info_file << "\n";
        }
        info_file.close();
    }
} 

void ModelFit::read_model_fit() {

    std::size_t dir_pos = file.find_last_of("/");
    std::string dir = file.substr(0, dir_pos + 1);

    std::ifstream lst_file(file);
    std::string line;
    std::getline(lst_file, line);
    model = line;
    
    if (model != "lm" && model != "lmm") { 
        
        std::getline(lst_file, line);
        std::string ww_filename = dir + line;
        std::getline(lst_file, line); 
        std::string info_filename = dir + line;

        std::ifstream ww_file(ww_filename);

        if (!ww_file.is_open()) {
            std::cerr << "Error: Unable to open working weight file " << ww_filename << std::endl;
            exit(1);
        }    

        std::getline(ww_file, line);
        std::vector<std::string> header_fields = string_split(line, "\t");
        sample_ids = std::vector<std::string>(header_fields.begin() + 1, header_fields.end());
        
        W = Eigen::MatrixXd(n_pheno, n_samples);
        pheno_ids.resize(n_pheno);
        
        size_t row = 0;
        while (std::getline(ww_file, line)) {
            std::vector<std::string> fields = string_split(line, "\t");
            pheno_ids[row] = fields[0];
            
            for (size_t i = 1; i < fields.size(); i++) {
                W(row, i-1) = std::stod(fields[i]);
            }
            row++;
        }
        ww_file.close();

        // Read information file.
        std::ifstream info_file(info_filename);

        if (!info_file.is_open()) {
            std::cerr << "Error: Unable to open fit information file " << info_filename << std::endl;
            exit(1);
        }    
        
        std::getline(info_file, line);
    
        if (model == "p_glm") {
            glm_converged.resize(n_pheno);
        } else if (model == "nb_glm") {
            phi.resize(n_pheno);
            phi_converged.resize(n_pheno);
            glm_converged.resize(n_pheno);
        } else if (model == "p_glmm") {
            tr.resize(n_pheno);
            glmm_converged.resize(n_pheno);
        } else if (model == "nb_glmm") {
            phi.resize(n_pheno);
            phi_converged.resize(n_pheno);
            tr.resize(n_pheno);
            glmm_converged.resize(n_pheno);
        }

        row = 0;
        while (std::getline(info_file, line)) {
            std::vector<std::string> fields = string_split(line, "\t");
        
            if (model == "p_glm") {
                glm_converged[row] = fields[1] == "1";
            } else if (model == "nb_glm") {
                glm_converged[row] = fields[1] == "1";
                phi[row] = std::stod(fields[2]);
                phi_converged[row] = fields[3] == "1";
            } else if (model == "p_glmm") {
                glmm_converged[row] = fields[1] == "1";
                tr[row] = std::stod(fields[2]);
            } else if (model == "nb_glmm") {
                glmm_converged[row] = fields[1] == "1";
                tr[row] = std::stod(fields[2]);
                phi[row] = std::stod(fields[3]);
                phi_converged[row] = fields[4] == "1";
            }
            row++;
        }
        info_file.close();
    }
    lst_file.close();
}

