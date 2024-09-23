
#include "Data.hpp"
#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void read_pheno_data(Param* params, PhenoData* pheno_data) {

    std::ifstream file(params->pheno_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open phenotype file " << params->pheno_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, "\t");
        tokens = string_split(line, "\t");
        if (tokens.size() < 2 || tokens[0] != "feature_id") {
            std::cerr << "Error: Invalid header in covariate file. Expected 'feature_id' as the first column." << std::endl;
            return;
        }
        pheno_data->sample_ids = std::vector<std::string>(tokens.begin() + 1, tokens.end());
        pheno_data->n_samples = pheno_data->sample_ids.size();
    }

    int n_pheno = 0;
    while (std::getline(file, line)) {
        if (!line.empty()) {
          n_pheno++;
        }
    }
    pheno_data->n_pheno = n_pheno;

    file.clear();
    file.seekg(0);
    // Skip header.
    std::getline(file, line);
    
    pheno_data->data = Eigen::MatrixXd(pheno_data->n_pheno, pheno_data->n_samples);
    pheno_data->pheno_ids.reserve(pheno_data->n_pheno);

    int row = 0;
    while (std::getline(file, line) && row < pheno_data->n_pheno) {
        tokens = string_split(line, "\t");
        if (tokens.size() != pheno_data->n_samples + 1) {
            std::cerr << "Error: Inconsistent number of columns in phenotype file." << std::endl;
            return;
        }
        pheno_data->pheno_ids.push_back(tokens[0]);
        for (int col = 0; col < pheno_data->n_samples; ++col) {
            pheno_data->data(row, col) = std::stod(tokens[col + 1]);
        }
        row++;
    }

    file.close();

    std::cout << "Read " << pheno_data->n_pheno << " features and " << pheno_data->n_samples << " samples from phenotype file." << std::endl;
}

void read_cov_data(Param* params, CovData* cov_data) {
    std::ifstream file(params->cov_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open covariate file: " << params->cov_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, "\t");
        if (tokens.size() < 2 || tokens[0] != "sample_id") {
            std::cerr << "Error: Invalid header in covariate file. Expected 'sample_id' as the first column." << std::endl;
            return;
        }
        cov_data->cov_ids = std::vector<std::string>(tokens.begin() + 1, tokens.end());
        cov_data->n_cov = cov_data->cov_ids.size();
    }

    cov_data->n_samples = 0;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            cov_data->n_samples++;
        }
    }

    file.clear();
    file.seekg(0);
    // Skip header.
    std::getline(file, line);

    cov_data->data = Eigen::MatrixXd(cov_data->n_samples, cov_data->n_cov);
    cov_data->sample_ids.reserve(cov_data->n_samples);

    int row = 0;
    while (std::getline(file, line) && row < cov_data->n_samples) {
        tokens = string_split(line, "\t");
        if (tokens.size() != cov_data->n_cov + 1) {
            std::cerr << "Error: Inconsistent number of columns in covariate file at line " << row + 2 << std::endl;
            return;
        }
        cov_data->sample_ids.push_back(tokens[0]);
        for (int col = 0; col < cov_data->n_cov; ++col) {
            cov_data->data(row, col) = std::stod(tokens[col + 1]);
        }
        row++;
    }

    file.close();

    std::cout << "Read " << cov_data->n_samples << " samples and " << cov_data->n_cov << " covariates from covariate file." << std::endl;
}

void read_feat_data(Param* params, FeatData* feat_data) {

    std::ifstream file(params->feat_anno_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open feature annotation file: " << params->feat_anno_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, "\t");
        if (tokens.size() < 4 || tokens[0] != "feature_id" || tokens[1] != "chrom" || tokens[2] != "start" || tokens[3] != "end") {
            std::cerr << "Error: Invalid header in feature annotation file. Expected 'feature_id', 'chrom', 'start', 'end' as the first four columns." << std::endl;
            return;
        }
    }

    feat_data->n_feat = 0;
    while (std::getline(file, line)) {
        tokens = string_split(line, "\t");
        if (tokens.size() < 4) {
            std::cerr << "Error: Inconsistent number of columns in feature annotation file at line " << feat_data->n_feat + 2 << std::endl;
            return;
        }

        feat_data->feat_id.push_back(tokens[0]);
        feat_data->chrom.push_back(std::stoi(tokens[1]));
        feat_data->start.push_back(std::stoi(tokens[2]));
        feat_data->end.push_back(std::stoi(tokens[3]));

        // Check if gene_name column exists
        if (tokens.size() > 4) {
            feat_data->gene_name.push_back(tokens[4]);
        } else {
            feat_data->gene_name.push_back("");
        }

        feat_data->n_feat++;
    }

    file.close();

    std::cout << "Read " << feat_data->n_feat << " features from feature annotation file." << std::endl;
}

void read_grm(Param* params, GRM* grm) {
    std::ifstream file(params->grm_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open GRM file: " << params->grm_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, "\t");
        if (tokens.size() < 2 || tokens[0] != "sample_id") {
            std::cerr << "Error: Invalid header in GRM file. Expected 'sample_id' as the first column." << std::endl;
            return;
        }
        grm->n_samps = tokens.size() - 1;
        grm->samp_ids = std::vector<std::string>(tokens.begin() + 1, tokens.end());
    }

    grm->mat = Eigen::MatrixXd::Zero(grm->n_samps, grm->n_samps);

    int row = 0;
    while (std::getline(file, line)) {
        tokens = string_split(line, "\t");
        if (tokens.size() != grm->n_samps + 1) {
            std::cerr << "Error: Inconsistent number of columns in GRM file at line " << row + 2 << std::endl;
            return;
        }

        for (int col = 0; col < grm->n_samps; ++col) {
            grm->mat(row, col) = std::stod(tokens[col + 1]);
        }

        row++;
    }

    if (row != grm->n_samps) {
        std::cerr << "Error: Number of rows does not match number of samples in GRM file." << std::endl;
        return;
    }

    file.close();

    std::cout << "Read GRM with " << grm->n_samps << " samples." << std::endl;
}
